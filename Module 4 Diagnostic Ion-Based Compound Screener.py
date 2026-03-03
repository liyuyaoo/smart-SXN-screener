import pandas as pd
import numpy as np
import os


# ==================== 通用函数 ====================
def find_diagnostic_ions(msms_spectrum, diagnostic_ions, tolerance=0.01):
    """
    在MSMS spectrum中查找诊断离子（必须同时存在所有指定的离子）

    参数:
    msms_spectrum: 碎片离子字符串，如 "57.04081:26 77.03999:26 78.04160:26 215.01700:79"
    diagnostic_ions: 诊断离子列表，如 [234.0551, 170.0942]
    tolerance: 质量容差

    返回:
    bool: 是否找到所有诊断离子
    """
    if pd.isna(msms_spectrum):
        return False

    spectrum_str = str(msms_spectrum)
    if not spectrum_str or spectrum_str.strip() == "":
        return False

    fragments = []
    try:
        for frag in spectrum_str.strip().split():
            mz_part = frag.split(':')[0]
            try:
                fragments.append(float(mz_part))
            except ValueError:
                continue
    except Exception as e:
        print(f"解析碎片离子时出错: {e}")
        return False

    for diag_ion in diagnostic_ions:
        found = False
        for frag in fragments:
            if abs(frag - diag_ion) <= tolerance:
                found = True
                break
        if not found:
            return False
    return True


def calculate_ppm(measured_mz, theoretical_mz):
    """
    计算实测值与理论值的ppm偏差
    """
    try:
        measured = float(measured_mz)
        theoretical = float(theoretical_mz)
        return abs(measured - theoretical) / theoretical * 1e6
    except (ValueError, TypeError):
        return float('inf')


def read_csv_with_encoding(file_path):
    """
    尝试多种编码读取CSV文件，返回DataFrame
    """
    encodings = ['utf-8', 'gbk', 'gb2312', 'gb18030', 'latin1', 'iso-8859-1', 'cp1252']
    for encoding in encodings:
        try:
            print(f"尝试使用 {encoding} 编码读取文件...")
            df = pd.read_csv(file_path, encoding=encoding)
            print(f"成功使用 {encoding} 编码读取文件")
            return df
        except (UnicodeDecodeError, Exception):
            continue
    raise Exception(f"无法使用任何支持的编码读取文件 {file_path}")


# ==================== 主程序 ====================
def main():
    # 使用全局配置变量（在文件末尾定义）
    global PEAK_TABLE_PATH, COMPOUND_LIB_PATH, OUTPUT_RESULTS_PATH, OUTPUT_DIAGNOSTIC_PATH
    global DIAGNOSTIC_IONS, TOLERANCE, PPM_THRESHOLD, MARK_EACH_ION, RT_COLUMN_POSSIBLE_NAMES

    # ========== 1. 读取峰表文件 ==========
    print("正在读取峰表文件...")
    try:
        peak_df = read_csv_with_encoding(PEAK_TABLE_PATH)

        required_columns = ['Precursor m/z', 'MSMS spectrum']
        for col in required_columns:
            if col not in peak_df.columns:
                print(f"错误: 峰表文件中缺少 '{col}' 列")
                return

        # 检测保留时间列
        rt_column = None
        for col in RT_COLUMN_POSSIBLE_NAMES:
            if col in peak_df.columns:
                rt_column = col
                print(f"找到保留时间列: '{rt_column}'")
                break
        if not rt_column:
            print("警告: 未找到保留时间列(RT)，将不包含保留时间信息")

        # 转换数值类型
        print("正在转换峰表中的母离子列为数值类型...")
        peak_df['Precursor m/z'] = pd.to_numeric(peak_df['Precursor m/z'], errors='coerce')
        if rt_column:
            peak_df[rt_column] = pd.to_numeric(peak_df[rt_column], errors='coerce')

        nan_count = peak_df['Precursor m/z'].isna().sum()
        if nan_count > 0:
            print(f"警告: 峰表中有 {nan_count} 个母离子值无法转换为数值，将被忽略")
            peak_df = peak_df.dropna(subset=['Precursor m/z'])

    except FileNotFoundError:
        print(f"错误: 未找到文件 {PEAK_TABLE_PATH}")
        return
    except Exception as e:
        print(f"读取峰表文件时出错: {e}")
        return

    # ========== 2. 筛选含有所有诊断离子的母离子 ==========
    print(f"正在查找含有诊断离子 {DIAGNOSTIC_IONS} 的母离子...")
    diagnostic_mask = peak_df['MSMS spectrum'].apply(
        lambda x: find_diagnostic_ions(x, DIAGNOSTIC_IONS, TOLERANCE)
    )
    diagnostic_peaks = peak_df[diagnostic_mask].copy()
    diagnostic_count = len(diagnostic_peaks)
    print(f"同时含有所有指定诊断离子的化合物数量: {diagnostic_count}")

    # ========== 3. 读取化合物库 ==========
    print("正在读取化合物库文件...")
    try:
        compound_df = read_csv_with_encoding(COMPOUND_LIB_PATH)

        if 'M+H' not in compound_df.columns:
            print("错误: 化合物库文件中缺少 'M+H' 列")
            return

        print("正在转换化合物库中的M+H列为数值类型...")
        compound_df['M+H'] = pd.to_numeric(compound_df['M+H'], errors='coerce')
        nan_count_compound = compound_df['M+H'].isna().sum()
        if nan_count_compound > 0:
            print(f"警告: 化合物库中有 {nan_count_compound} 个M+H值无法转换为数值，将被忽略")
            compound_df = compound_df.dropna(subset=['M+H'])

    except FileNotFoundError:
        print(f"错误: 未找到文件 {COMPOUND_LIB_PATH}")
        return
    except Exception as e:
        print(f"读取化合物库文件时出错: {e}")
        return

    # ========== 4. 匹配母离子与化合物库，并（可选）标记单个诊断离子存在情况 ==========
    print("正在匹配母离子与化合物库...")
    matched_results = []
    match_count = 0
    error_count = 0

    # 如果需要标记每个离子，预先为每个诊断离子生成检查函数（避免重复解析谱图）
    if MARK_EACH_ION:
        # 对每个诊断离子分别检查的函数（重用谱图字符串）
        def check_each_ion(spectrum_str):
            # 解析一次碎片
            fragments = []
            for frag in spectrum_str.strip().split():
                mz_part = frag.split(':')[0]
                try:
                    fragments.append(float(mz_part))
                except ValueError:
                    continue
            # 检查每个离子
            ion_present = {}
            for ion in DIAGNOSTIC_IONS:
                found = any(abs(f - ion) <= TOLERANCE for f in fragments)
                ion_present[ion] = found
            return ion_present
    else:
        check_each_ion = None  # 不需要标记

    for idx, row in diagnostic_peaks.iterrows():
        precursor_mz = row['Precursor m/z']
        if pd.isna(precursor_mz):
            continue

        # 如果需要标记，获取每个离子的存在情况
        if MARK_EACH_ION:
            ion_present = check_each_ion(str(row['MSMS spectrum']))
        else:
            ion_present = {}

        for _, compound_row in compound_df.iterrows():
            theoretical_mz = compound_row['M+H']
            if pd.isna(theoretical_mz):
                continue

            ppm = calculate_ppm(precursor_mz, theoretical_mz)
            if ppm == float('inf'):
                error_count += 1
                continue

            if ppm < PPM_THRESHOLD:
                match_count += 1
                result_row = {
                    '母离子 (Precursor m/z)': precursor_mz,
                    '理论M+H': theoretical_mz,
                    'ppm': ppm,
                }

                # 添加每个诊断离子的标记列
                for ion, present in ion_present.items():
                    # 列名使用离子值，如 "Has_234.0551"
                    result_row[f'Has_{ion}'] = present

                # 添加保留时间
                if rt_column:
                    result_row['保留时间 (RT)'] = row[rt_column]
                else:
                    result_row['保留时间 (RT)'] = None

                # 添加化合物库中的其他列
                for col in compound_df.columns:
                    if col != 'M+H':
                        result_row[col] = compound_row[col]

                # 添加原始MSMS谱
                result_row['MSMS spectrum'] = row['MSMS spectrum']

                matched_results.append(result_row)

                if match_count % 10 == 0:
                    print(f"已匹配 {match_count} 个化合物...")

    print(f"匹配过程中遇到 {error_count} 个计算错误")

    # ========== 5. 创建结果DataFrame ==========
    matched_df = pd.DataFrame(matched_results)

    if not matched_df.empty:
        # 构建列顺序
        columns_order = ['母离子 (Precursor m/z)', '理论M+H', 'ppm']
        if MARK_EACH_ION:
            # 将标记列放在ppm之后
            for ion in DIAGNOSTIC_IONS:
                columns_order.append(f'Has_{ion}')
        if rt_column:
            columns_order.insert(2, '保留时间 (RT)')  # 放在母离子之后、理论M+H之前？按原习惯放在第二或第三位置，这里简单插入到第二
            # 注意：上面已经插入了，调整位置以符合原模块习惯：母离子、保留时间、理论M+H、ppm...
            # 简便起见，我们重新构建顺序
            columns_order = ['母离子 (Precursor m/z)', '保留时间 (RT)', '理论M+H', 'ppm']
            if MARK_EACH_ION:
                for ion in DIAGNOSTIC_IONS:
                    columns_order.append(f'Has_{ion}')

        # 添加化合物库的其他列
        for col in compound_df.columns:
            if col != 'M+H' and col not in columns_order:
                columns_order.append(col)

        # 添加MSMS spectrum列
        if 'MSMS spectrum' not in columns_order:
            columns_order.append('MSMS spectrum')

        # 确保所有存在的列都在columns_order中（防止遗漏）
        existing_cols = [c for c in columns_order if c in matched_df.columns]
        matched_df = matched_df[existing_cols]

    # ========== 6. 输出结果 ==========
    print("\n" + "=" * 50)
    print(f"同时含有诊断离子 {DIAGNOSTIC_IONS} 的化合物总数: {diagnostic_count}")
    print(f"与化合物库匹配的数量 (ppm < {PPM_THRESHOLD}): {len(matched_df)}")

    if not matched_df.empty:
        matched_df.to_csv(OUTPUT_RESULTS_PATH, index=False, encoding='utf-8-sig')
        print(f"匹配结果已保存到: {OUTPUT_RESULTS_PATH}")
        print("\n前5行匹配结果:")
        print(matched_df.head().to_string())
    else:
        print("未找到匹配的化合物")

    # 保存筛选出的诊断离子母离子列表
    if rt_column and rt_column not in diagnostic_peaks.columns:
        diagnostic_peaks[rt_column] = peak_df.loc[diagnostic_peaks.index, rt_column]
    diagnostic_peaks.to_csv(OUTPUT_DIAGNOSTIC_PATH, index=False, encoding='utf-8-sig')
    print(f"同时含有诊断离子的母离子列表已保存到: {OUTPUT_DIAGNOSTIC_PATH}")

    # ========== 7. 统计信息 ==========
    print("\n统计信息:")
    print(f"- 原始峰表中的总化合物数: {len(peak_df)}")
    print(f"- 同时含有诊断离子 {DIAGNOSTIC_IONS} 的化合物数: {diagnostic_count}")
    print(f"- 与化合物库匹配的化合物数: {len(matched_df)}")

    if not matched_df.empty:
        print(f"- 平均ppm: {matched_df['ppm'].mean():.2f}")
        print(f"- 最小ppm: {matched_df['ppm'].min():.2f}")
        print(f"- 最大ppm: {matched_df['ppm'].max():.2f}")

        print("\n最佳匹配 (ppm最小):")
        best_matches = matched_df.nsmallest(5, 'ppm')
        display_cols = ['母离子 (Precursor m/z)', '理论M+H', 'ppm']
        if MARK_EACH_ION:
            for ion in DIAGNOSTIC_IONS:
                display_cols.append(f'Has_{ion}')
        if rt_column:
            display_cols.insert(2, '保留时间 (RT)')
        # 只保留存在的列
        display_cols = [c for c in display_cols if c in best_matches.columns]
        print(best_matches[display_cols].to_string(index=False))

        if rt_column and matched_df['保留时间 (RT)'].notna().any():
            print(f"\n保留时间分布:")
            print(f"- 最小RT: {matched_df['保留时间 (RT)'].min():.2f}")
            print(f"- 最大RT: {matched_df['保留时间 (RT)'].max():.2f}")
            print(f"- 平均RT: {matched_df['保留时间 (RT)'].mean():.2f}")

    print("\n" + "=" * 50)
    print("分析完成！")


# ==================== 用户配置区域 ====================
# 请根据实际情况修改以下参数

# ----- 文件路径 -----
PEAK_TABLE_PATH = "F:/20251217-LYY-YSH/fengbiao-hrr-cl.csv"  # 峰表文件路径
COMPOUND_LIB_PATH = "G:/2026.1-各种AI预测/DNS-CL预测/氨基酸库-代码用.csv"  # 化合物库路径
OUTPUT_RESULTS_PATH = "C:/Users/liyuyao/Desktop/matching_results.csv"  # 匹配结果输出路径
OUTPUT_DIAGNOSTIC_PATH = "C:/Users/liyuyao/Desktop/diagnostic_peaks.csv"  # 诊断离子母离子列表输出路径

# ----- 诊断离子参数 -----
DIAGNOSTIC_IONS = [234.0551, 170.0942]  # 诊断离子列表（必须同时存在）
TOLERANCE = 0.01  # 质量容差（绝对偏差）
PPM_THRESHOLD = 10  # ppm匹配阈值

# ----- 输出选项 -----
MARK_EACH_ION = True  # 是否在结果中分别标记每个诊断离子的存在情况（True则添加Has_<ion>列）

# ----- 保留时间列可能名称（程序会自动检测第一个存在的列）-----
RT_COLUMN_POSSIBLE_NAMES = ['RT', 'rt', 'Retention Time', 'retention time', 'Retention_Time', 'retention_time']

# =====================================================

if __name__ == "__main__":
    main()