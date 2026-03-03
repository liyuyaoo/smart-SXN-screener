import re
import os
import csv
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict

# ================================
# 配置区域 - 请在这里填写所有路径和参数
# ================================

CONFIG = {
    # 基础数据文件路径
    "aglycones_file": "C:/Users/liyuyao/Desktop/aglycones.csv",
    "sugars_file": "C:/Users/liyuyao/Desktop/sugars.csv",
    "acyls_file": "C:/Users/liyuyao/Desktop/acyls.csv",

    # 实验数据文件路径
    "experimental_file": "D:/lilunku-xianji/MSDIAL-N.csv",

    # 输出文件路径
    "output_file": "C:/Users/liyuyao/Desktop/匹配结果-N.csv",

    # 离子类型选择 (填写 "M+H" 或 "M-H")
    "ion_type": "M-H",

    # PPM阈值
    "ppm_threshold": 10
}

# ================================
# 以下为程序代码，无需修改
# ================================

# 定义原子质量（使用最精确的数值）
ATOMIC_MASSES = {
    'H': 1.00782503223,
    'C': 12.0000000000,
    'O': 15.99491461957,
    'N': 14.00307400443,
    'S': 31.9720711744,
    'P': 30.97376199842,
}

# 质子和电子质量
PROTON_MASS = 1.007276466812
ELECTRON_MASS = 0.000548579909


# 从CSV文件读取数据的函数
def load_data_from_csv(file_path):
    data_dict = {}
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            # 读取表头
            header = next(reader, None)

            # 检查表头是否符合预期
            if header and len(header) >= 2:
                print(f"文件 {file_path} 的表头: {header}")

            for row in reader:
                if len(row) >= 2:
                    name = row[0].strip()
                    formula = row[1].strip()
                    if name and formula:
                        data_dict[name] = formula
    except FileNotFoundError:
        print(f"错误: 文件 {file_path} 未找到!")
        return None
    except Exception as e:
        print(f"读取文件 {file_path} 时出错: {e}")
        return None
    return data_dict


# 解析分子式并转换为原子计数字典的函数
def parse_formula(formula):
    elements = {}
    pattern = r'([A-Z][a-z]*)(\d*)'
    for element, count in re.findall(pattern, formula):
        count = int(count) if count else 1
        elements[element] = elements.get(element, 0) + count
    return elements


# 将原子计数字典转换回分子式字符串的函数
def format_formula(elements_dict):
    # 定义元素输出顺序
    element_order = ['C', 'H', 'O', 'N', 'S', 'P']
    formula_parts = []

    for element in element_order:
        if element in elements_dict and elements_dict[element] > 0:
            count = elements_dict[element]
            formula_parts.append(f"{element}{count if count > 1 else ''}")

    # 添加可能不在预设顺序中的元素
    for element, count in elements_dict.items():
        if element not in element_order and count > 0:
            formula_parts.append(f"{element}{count if count > 1 else ''}")

    return ''.join(formula_parts)


# 计算两个分子式相加的函数
def add_formulas(formula1, formula2):
    elem1 = parse_formula(formula1)
    elem2 = parse_formula(formula2)

    result = {}
    for element in set(elem1.keys()) | set(elem2.keys()):
        result[element] = elem1.get(element, 0) + elem2.get(element, 0)

    return format_formula(result)


# 计算化合物的精确分子量
def calculate_exact_mass(formula):
    elements = parse_formula(formula)
    mass = 0.0

    for element, count in elements.items():
        if element not in ATOMIC_MASSES:
            raise ValueError(f"未知元素: {element}")
        mass += ATOMIC_MASSES[element] * count

    return mass


# 计算化合物的[M+H]+准确分子量
def calculate_m_plus_h_mass(formula):
    exact_mass = calculate_exact_mass(formula)
    return exact_mass + PROTON_MASS


# 计算化合物的[M-H]-准确分子量
def calculate_m_minus_h_mass(formula):
    exact_mass = calculate_exact_mass(formula)
    return exact_mass - PROTON_MASS


# 生成理论库（糖苷和酰化糖苷合并）
def generate_theoretical_library(aglycones, sugars, acyls):
    results = []

    print("正在生成理论库...")

    # 生成糖苷（苷元+糖）
    print("生成糖苷...")
    for aglycone_name, aglycone_formula in aglycones.items():
        for sugar_name, sugar_formula in sugars.items():
            combined_formula = add_formulas(aglycone_formula, sugar_formula)
            compound_name = f"{aglycone_name}-{sugar_name}"
            results.append((compound_name, combined_formula))

    # 生成酰化糖苷（苷元+糖+酰基）
    print("生成酰化糖苷...")
    for aglycone_name, aglycone_formula in aglycones.items():
        for sugar_name, sugar_formula in sugars.items():
            for acyl_name, acyl_formula in acyls.items():
                combined_formula = add_formulas(aglycone_formula, sugar_formula)
                combined_formula = add_formulas(combined_formula, acyl_formula)
                compound_name = f"{aglycone_name}-{sugar_name}-{acyl_name}"
                results.append((compound_name, combined_formula))

    print(f"理论库生成完成，共 {len(results)} 个化合物")
    return results


# 计算理论分子量
def calculate_theoretical_masses(compounds, ion_type="M+H"):
    results = []

    print(f"正在计算理论分子量 ({ion_type})...")

    for compound_name, formula in compounds:
        try:
            if ion_type == "M+H":
                mass = calculate_m_plus_h_mass(formula)
            elif ion_type == "M-H":
                mass = calculate_m_minus_h_mass(formula)
            else:
                mass = calculate_exact_mass(formula)

            results.append((compound_name, formula, mass))
        except Exception as e:
            print(f"计算化合物 {compound_name} ({formula}) 时出错: {e}")

    return results


# 计算PPM差异
def calculate_ppm_difference(theoretical_mass, experimental_mass):
    return abs(theoretical_mass - experimental_mass) / theoretical_mass * 1e6


# 匹配理论分子量和实验分子量
def match_masses(theoretical_data, experimental_data, ppm_threshold=10):
    matches = []

    print("正在匹配理论分子量和实验分子量...")

    for exp_no, exp_mass in experimental_data:
        for theo_name, theo_formula, theo_mass in theoretical_data:
            ppm_diff = calculate_ppm_difference(theo_mass, exp_mass)

            if ppm_diff <= ppm_threshold:
                matches.append({
                    'NO.': exp_no,
                    'Experimental_mz': exp_mass,
                    'Compound_Name': theo_name,
                    'Molecular_Formula': theo_formula,
                    'Theoretical_Mass': theo_mass,
                    'PPM_Difference': ppm_diff
                })

    return matches


# 去重并整合化合物名称
def deduplicate_and_merge(matches):
    # 按分子式分组
    formula_groups = {}

    for match in matches:
        formula = match['Molecular_Formula']
        if formula not in formula_groups:
            formula_groups[formula] = []
        formula_groups[formula].append(match)

    # 对每个分子式组，选择PPM差异最小的匹配，并整合化合物名称
    results = []

    for formula, group_matches in formula_groups.items():
        # 找到PPM差异最小的匹配
        best_match = min(group_matches, key=lambda x: x['PPM_Difference'])

        # 整合所有化合物名称
        all_names = list(set([match['Compound_Name'] for match in group_matches]))
        merged_names = '; '.join(all_names)

        # 创建新的结果条目
        result_entry = best_match.copy()
        result_entry['Compound_Name'] = merged_names
        result_entry['All_Compound_Count'] = len(all_names)

        results.append(result_entry)

    return results


# 读取实验数据
def load_experimental_data(file_path):
    experimental_data = []

    try:
        df = pd.read_csv(file_path)

        # 检查必要的列
        if 'NO.' not in df.columns or 'm/z' not in df.columns:
            print("错误: 实验数据文件必须包含 'NO.' 和 'm/z' 列")
            print(f"当前文件的列名: {list(df.columns)}")
            return None

        for _, row in df.iterrows():
            no = row['NO.']
            mz = float(row['m/z'])
            experimental_data.append((no, mz))

        print(f"成功加载 {len(experimental_data)} 个实验数据点")
        return experimental_data

    except Exception as e:
        print(f"读取实验数据文件时出错: {e}")
        return None


# 主程序
def main():
    # 使用CONFIG中的配置
    config = CONFIG

    print(f"配置信息:")
    print(f"苷元文件: {config['aglycones_file']}")
    print(f"糖文件: {config['sugars_file']}")
    print(f"酰基文件: {config['acyls_file']}")
    print(f"实验数据文件: {config['experimental_file']}")
    print(f"输出文件: {config['output_file']}")
    print(f"离子类型: {config['ion_type']}")
    print(f"PPM阈值: {config['ppm_threshold']}")

    # 从CSV文件加载数据
    print("\n正在加载基础数据...")
    aglycones = load_data_from_csv(config['aglycones_file'])
    sugars = load_data_from_csv(config['sugars_file'])
    acyls = load_data_from_csv(config['acyls_file'])

    # 检查数据是否成功加载
    if aglycones is None or sugars is None or acyls is None:
        print("错误: 无法加载所有必需的数据文件，程序终止。")
        return

    print(f"加载了 {len(aglycones)} 个苷元, {len(sugars)} 个糖, {len(acyls)} 个酰基")

    # 检查是否有数据
    if not aglycones or not sugars or not acyls:
        print("错误: 基础数据为空!")
        return

    # 加载实验数据
    experimental_data = load_experimental_data(config['experimental_file'])
    if experimental_data is None:
        return

    # 生成理论库
    theoretical_compounds = generate_theoretical_library(aglycones, sugars, acyls)

    # 计算理论分子量
    theoretical_data = calculate_theoretical_masses(theoretical_compounds, config['ion_type'])

    # 匹配理论分子量和实验分子量
    matches = match_masses(theoretical_data, experimental_data, config['ppm_threshold'])

    if not matches:
        print("没有找到匹配的化合物")
        return

    print(f"找到 {len(matches)} 个初步匹配项")

    # 去重并整合化合物名称
    final_results = deduplicate_and_merge(matches)

    print(f"去重后保留 {len(final_results)} 个唯一分子式")

    # 保存结果
    try:
        # 创建输出目录（如果不存在）
        os.makedirs(os.path.dirname(config['output_file']), exist_ok=True)

        # 转换为DataFrame并保存
        result_df = pd.DataFrame(final_results)

        # 重新排列列的顺序
        columns_order = ['NO.', 'Experimental_mz', 'Molecular_Formula', 'Compound_Name',
                         'All_Compound_Count', 'Theoretical_Mass', 'PPM_Difference']

        # 只保留存在的列
        columns_order = [col for col in columns_order if col in result_df.columns]
        result_df = result_df[columns_order]

        # 按PPM差异排序
        result_df = result_df.sort_values('PPM_Difference')

        result_df.to_csv(config['output_file'], index=False)
        print(f"\n结果已保存到: {config['output_file']}")

        # 显示统计信息
        avg_ppm = np.mean(result_df['PPM_Difference'])
        max_ppm = np.max(result_df['PPM_Difference'])
        min_ppm = np.min(result_df['PPM_Difference'])

        print(f"\n匹配统计:")
        print(f"平均PPM差异: {avg_ppm:.4f}")
        print(f"最大PPM差异: {max_ppm:.4f}")
        print(f"最小PPM差异: {min_ppm:.4f}")

        # 显示前几个匹配结果
        print(f"\n前5个匹配结果:")
        for i, (_, row) in enumerate(result_df.head(5).iterrows()):
            print(f"{i + 1}. 分子式: {row['Molecular_Formula']}")
            print(f"   化合物名称: {row['Compound_Name']}")
            print(
                f"   实验m/z: {row['Experimental_mz']:.6f}, 理论质量: {row['Theoretical_Mass']:.6f}, PPM: {row['PPM_Difference']:.4f}")

    except Exception as e:
        print(f"保存结果时出错: {e}")


if __name__ == "__main__":
    main()