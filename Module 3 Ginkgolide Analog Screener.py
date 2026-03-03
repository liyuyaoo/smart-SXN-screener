import pandas as pd
import re

# 常量定义 - 精确质量值
CO_LOSS = 27.99546
TWO_CO_LOSS = 2 * CO_LOSS  # 55.99092
CO2_LOSS = 43.99038
H2O_LOSS = 18.01111
CO2_H2O_LOSS = 62.00094  # CO₂+H₂O
FEATURE_LOSS = 72.99312  # 特征中性丢失
MZ_125_02442 = 125.02442
MZ_141_01933 = 141.01933
TOLERANCE = 0.01  # 容差


def parse_msms(msms_str):
    """
    解析MS/MS字符串，返回碎片离子的m/z值列表
    格式示例: "57.04081:26 77.03999:26 78.04160:26 215.01700:79"
    冒号后面是响应强度，我们只需要m/z值
    """
    if pd.isna(msms_str) or not isinstance(msms_str, str) or msms_str.strip() == "":
        return []

    fragments = []
    # 使用正则表达式匹配所有"数字:数字"格式的片段
    pattern = r'([\d\.]+):\d+'
    matches = re.findall(pattern, msms_str)

    for match in matches:
        try:
            fragment_mz = float(match)
            fragments.append(fragment_mz)
        except ValueError:
            continue

    return fragments


def check_first_level(parent_mz, fragments):
    """
    第一级筛选：检查是否存在特征中性丢失72.99312以及母离子是否丢失CO或2CO
    """
    # 检查特征中性丢失72.99312
    has_feature_loss = any(
        abs(parent_mz - fragment - FEATURE_LOSS) <= TOLERANCE
        for fragment in fragments
    )

    # 检查母离子是否丢失CO (27.99546) 或 2CO (55.99092)
    has_co_loss = any(
        abs(parent_mz - fragment - CO_LOSS) <= TOLERANCE or
        abs(parent_mz - fragment - TWO_CO_LOSS) <= TOLERANCE
        for fragment in fragments
    )

    return has_feature_loss and has_co_loss


def check_neutral_loss(parent_mz, fragments):
    """
    检查是否发生进一步的中性丢失，如丢失CO₂、H₂O或两者
    """
    losses_to_check = [CO2_LOSS, H2O_LOSS, CO2_H2O_LOSS]

    for loss in losses_to_check:
        if any(abs(parent_mz - fragment - loss) <= TOLERANCE for fragment in fragments):
            return True

    return False


def ginkgolide_screening(parent_mz, fragments):
    """
    完整的银杏内酯筛选流程
    返回: (是否为银杏内酯, 分类结果)
    """
    # 第一级筛选
    if not check_first_level(parent_mz, fragments):
        return (False, "非银杏内酯类化合物")

    # 检查是否存在碎片离子 m/z 125.02442
    has_125 = any(abs(fragment - MZ_125_02442) <= TOLERANCE for fragment in fragments)

    if not has_125:
        # 检查是否发生进一步的中性丢失
        if check_neutral_loss(parent_mz, fragments):
            return (True, "银杏内酯类似物")
        else:
            return (True, "未分类银杏内酯（有特征丢失但无125.02442和无中性丢失）")
    else:
        # 检查是否同时存在 m/z 141.01933
        has_141 = any(abs(fragment - MZ_141_01933) <= TOLERANCE for fragment in fragments)

        if has_141:
            return (True, "R₁=OH和R₃=OH型银杏内酯")
        else:
            return (True, "可能为R1=OH,C-C16为双键的银杏内酯类似物")


def main():
    # 文件路径
    input_file = r'C:\Users\liyuyao\Desktop\SXN-N-30-fengbiao-hrr.csv'
    output_file = r'C:\Users\liyuyao\Desktop\N-内酯-30CE\ginkgolide_3.csv'

    print("开始筛选银杏内酯类化合物...")
    print(f"使用以下精确质量值:")
    print(f"  CO丢失: {CO_LOSS}")
    print(f"  2CO丢失: {TWO_CO_LOSS}")
    print(f"  CO₂丢失: {CO2_LOSS}")
    print(f"  H₂O丢失: {H2O_LOSS}")
    print(f"  CO₂+H₂O丢失: {CO2_H2O_LOSS}")
    print(f"  特征丢失: {FEATURE_LOSS}")
    print(f"  碎片m/z 125.02442")
    print(f"  碎片m/z 141.01933")
    print(f"  容差: ±{TOLERANCE}")

    # 读取CSV文件
    try:
        df = pd.read_csv(input_file, encoding='utf-8')
        print(f"\n成功读取文件: {input_file}")
        print(f"数据形状: {df.shape}")
        print(f"列名: {list(df.columns)}")
    except Exception as e:
        print(f"读取文件失败: {e}")
        # 尝试其他编码
        try:
            df = pd.read_csv(input_file, encoding='gbk')
            print(f"使用GBK编码成功读取文件: {input_file}")
        except:
            print("请检查文件路径和编码格式")
            return

    # 检查必要的列是否存在
    required_columns = ['m/z', 'MS/MS', 'RT']
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        print(f"缺少必要的列: {missing_columns}")
        print(f"现有列: {list(df.columns)}")
        return

    # 确定MS/MS列的列名（根据您的描述，应该是'MS/MS'）
    msms_column = 'MS/MS'
    if msms_column not in df.columns:
        # 尝试其他可能的列名
        possible_msms_names = ['MS/MS', 'MSMS', 'MS2', 'msms', 'ms/ms', 'Fragments', 'fragments', '二级碎片']
        for col in df.columns:
            if col in possible_msms_names:
                msms_column = col
                break

        if msms_column not in df.columns:
            print("未找到MS/MS数据列，请检查数据文件")
            print(f"现有列: {list(df.columns)}")
            return

    print(f"使用 '{msms_column}' 作为MS/MS数据列")

    # 初始化结果列
    df['是否为银杏内酯'] = '否'
    df['分类结果'] = ''
    df['第一级筛选'] = '不满足'

    # 添加详细筛选信息列
    df['是否有特征丢失和CO/2CO丢失'] = '否'
    df['是否有125.02442'] = '否'
    df['是否有141.01933'] = '否'
    df['是否有中性丢失(CO2/H2O/两者)'] = '否'

    # 进度跟踪
    total_rows = len(df)
    print(f"\n开始处理 {total_rows} 个化合物...")

    # 遍历每一行数据
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"处理进度: {idx + 1}/{total_rows}")

        try:
            # 获取母离子m/z
            parent_mz = float(row['m/z'])

            # 解析MS/MS数据
            fragments = parse_msms(row[msms_column])

            if not fragments:
                continue  # 如果没有碎片数据，跳过

            # 第一级筛选
            first_level_pass = check_first_level(parent_mz, fragments)
            df.at[idx, '是否有特征丢失和CO/2CO丢失'] = '是' if first_level_pass else '否'
            df.at[idx, '第一级筛选'] = '满足' if first_level_pass else '不满足'

            # 如果不是银杏内酯，跳过后续检查
            if not first_level_pass:
                df.at[idx, '分类结果'] = "非银杏内酯类化合物"
                continue

            # 检查是否有125.02442
            has_125 = any(abs(fragment - MZ_125_02442) <= TOLERANCE for fragment in fragments)
            df.at[idx, '是否有125.02442'] = '是' if has_125 else '否'

            if not has_125:
                # 检查中性丢失
                has_neutral_loss = check_neutral_loss(parent_mz, fragments)
                df.at[idx, '是否有中性丢失(CO2/H2O/两者)'] = '是' if has_neutral_loss else '否'

                if has_neutral_loss:
                    df.at[idx, '是否为银杏内酯'] = '是'
                    df.at[idx, '分类结果'] = "银杏内酯类似物"
                else:
                    df.at[idx, '是否为银杏内酯'] = '是'  # 有特征丢失但无125.02442和无中性丢失
                    df.at[idx, '分类结果'] = "未分类银杏内酯（有特征丢失但无125.02442和无中性丢失）"
            else:
                # 检查是否有141.01933
                has_141 = any(abs(fragment - MZ_141_01933) <= TOLERANCE for fragment in fragments)
                df.at[idx, '是否有141.01933'] = '是' if has_141 else '否'

                if has_141:
                    df.at[idx, '是否为银杏内酯'] = '是'
                    df.at[idx, '分类结果'] = "R₁=OH和R₃=OH型银杏内酯"
                else:
                    df.at[idx, '是否为银杏内酯'] = '是'
                    df.at[idx, '分类结果'] = "可能为R1=OH,C-C16为双键的银杏内酯类似物"

        except Exception as e:
            print(f"处理第{idx + 1}行时出错: {e}")
            continue

    print(f"处理完成，已处理 {total_rows} 个化合物")

    # 保存结果
    try:
        # 重新排列列的顺序，使结果列更容易阅读
        columns_order = ['m/z', 'RT', 'MS/MS', '是否为银杏内酯', '分类结果', '第一级筛选',
                         '是否有特征丢失和CO/2CO丢失', '是否有125.02442', '是否有141.01933',
                         '是否有中性丢失(CO2/H2O/两者)']

        # 确保所有列都存在
        existing_columns = [col for col in columns_order if col in df.columns]
        other_columns = [col for col in df.columns if col not in existing_columns]
        final_columns_order = existing_columns + other_columns

        df = df[final_columns_order]
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\n结果已保存到: {output_file}")

        # 统计结果
        total_compounds = len(df)
        ginkgolide_count = len(df[df['是否为银杏内酯'] == '是'])

        print(f"\n筛选结果统计:")
        print(f"总化合物数量: {total_compounds}")
        print(f"银杏内酯类化合物数量: {ginkgolide_count}")

        if total_compounds > 0:
            print(f"银杏内酯占比: {ginkgolide_count / total_compounds * 100:.2f}%")

        # 分类统计
        if ginkgolide_count > 0:
            print(f"\n分类结果统计:")
            classification_counts = df['分类结果'].value_counts()
            for category, count in classification_counts.items():
                if category and category != '':
                    print(f"  {category}: {count}")

            print(f"\n第一级筛选结果:")
            first_level_counts = df['第一级筛选'].value_counts()
            for result, count in first_level_counts.items():
                if result and result != '':
                    print(f"  {result}: {count}")

        # 输出前10个银杏内酯类化合物的信息
        ginkgolides = df[df['是否为银杏内酯'] == '是'].head(10)
        if not ginkgolides.empty:
            print(f"\n前10个银杏内酯类化合物:")
            for idx, row in ginkgolides.iterrows():
                print(f"  m/z: {row['m/z']}, RT: {row['RT']}, 分类: {row['分类结果']}")
        else:
            print("\n未筛选到银杏内酯类化合物")

    except Exception as e:
        print(f"保存结果失败: {e}")


if __name__ == "__main__":
    main()