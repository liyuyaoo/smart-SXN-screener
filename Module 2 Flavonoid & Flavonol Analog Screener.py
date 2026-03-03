import pandas as pd
import numpy as np
import os
import re
from datetime import datetime


class FlavonoidScreener:
    def __init__(self, mass_tolerance=0.01):
        """
        初始化黄酮筛选器

        Parameters:
        -----------
        mass_tolerance : float
            质量容差，默认±0.01
        """
        self.mass_tolerance = mass_tolerance

        # 文件路径
        self.aglycone_file = r"G:/2026.1-各种AI预测/黄酮预测/66苷元信息代码.csv"
        self.peak_file = r"G:/2026.1-各种AI预测/黄酮预测/SXN-P-45-fengbiao-hrr.csv"
        self.output_dir = r"G:/2026.1-各种AI预测/黄酮预测/results-45/5"

        # 创建输出文件夹
        os.makedirs(self.output_dir, exist_ok=True)

        # 黄酮特征碎片 (m/z值)
        self.flavonoid_fragments = {
            1: 137.02332,  # 1个羟基
            2: 153.01824,  # 2个羟基
            3: 169.01315,  # 3个羟基
            4: 185.00806  # 4个羟基
        }

        # 黄酮醇特征碎片 (m/z值)
        self.flavonol_fragments = {
            1: 149.02332,  # 1个羟基
            2: 165.01824,  # 2个羟基
            3: 181.01315,  # 3个羟基
            4: 197.00806  # 4个羟基
        }

    def load_data(self):
        """加载数据文件"""
        print("=" * 70)
        print("加载数据文件...")
        print("=" * 70)

        # 1. 加载苷元数据
        print(f"1. 加载苷元数据: {self.aglycone_file}")
        try:
            # 尝试不同的编码
            encodings = ['utf-8', 'gbk', 'gb2312', 'latin1', 'cp1252']
            for encoding in encodings:
                try:
                    self.aglycone_df = pd.read_csv(self.aglycone_file, encoding=encoding)
                    print(f"   使用 {encoding} 编码成功加载")
                    break
                except:
                    continue
            else:
                print("   错误: 无法读取苷元文件")
                return False

            # 检查列名
            self.aglycone_df.columns = self.aglycone_df.columns.str.strip()
            print(f"   列名: {list(self.aglycone_df.columns)}")
            print(f"   数据形状: {self.aglycone_df.shape}")

            # 显示Aglycone6的信息
            if 'Aglycone6' in self.aglycone_df['Name'].values:
                aglycone6 = self.aglycone_df[self.aglycone_df['Name'] == 'Aglycone6'].iloc[0]
                print(f"\n   苷元Aglycone6信息:")
                print(f"     名称: {aglycone6['Name']}")
                print(f"     分子式: {aglycone6['Formula']}")
                print(f"     m/z: {aglycone6['aglycone m/z']}")

        except Exception as e:
            print(f"   加载苷元数据时出错: {e}")
            return False

        # 2. 加载峰表数据
        print(f"\n2. 加载峰表数据: {self.peak_file}")
        try:
            for encoding in encodings:
                try:
                    self.peak_df = pd.read_csv(self.peak_file, encoding=encoding)
                    print(f"   使用 {encoding} 编码成功加载")
                    break
                except:
                    continue
            else:
                print("   错误: 无法读取峰表文件")
                return False

            # 检查列名
            self.peak_df.columns = self.peak_df.columns.str.strip()
            print(f"   列名: {list(self.peak_df.columns)}")
            print(f"   数据形状: {self.peak_df.shape}")

        except Exception as e:
            print(f"   加载峰表数据时出错: {e}")
            return False

        print("\n数据加载完成!")
        return True

    def parse_msms(self, msms_str):
        """
        解析MS/MS字符串，提取碎片m/z和强度
        支持多种分隔符格式：空格、逗号、分号

        Parameters:
        -----------
        msms_str : str
            MS/MS字符串，格式如 "57.04081:26 77.03999:26" 或 "57.04081:26, 77.03999:26"

        Returns:
        --------
        dict
            碎片m/z和强度的字典
        """
        fragments = {}

        if pd.isna(msms_str) or not isinstance(msms_str, str):
            return fragments

        try:
            # 清理字符串
            msms_str = msms_str.strip()
            if not msms_str:
                return fragments

            # 支持多种分隔符：先统一替换为空格
            # 替换逗号、分号、制表符为空格
            msms_str = re.sub(r'[,\t;]+', ' ', msms_str)

            # 分割字符串
            items = msms_str.split()

            for item in items:
                if ':' in item:
                    try:
                        parts = item.split(':')
                        if len(parts) >= 2:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            fragments[mz] = intensity
                    except:
                        continue

        except Exception as e:
            print(f"解析MS/MS时出错: {e}")

        return fragments

    def is_match(self, value1, value2):
        """检查两个值是否在容差范围内匹配"""
        return abs(value1 - value2) <= self.mass_tolerance

    def find_aglycone_match(self, peak_mz):
        """查找匹配的苷元"""
        matches = []

        for _, row in self.aglycone_df.iterrows():
            try:
                aglycone_mz = float(row['aglycone m/z'])
                if self.is_match(peak_mz, aglycone_mz):
                    matches.append({
                        'Name': row['Name'],
                        'Formula': row['Formula'],
                        'aglycone_m/z': aglycone_mz,
                        'Mass_Difference': abs(peak_mz - aglycone_mz),
                        'Matched_Fragment_m/z': peak_mz,  # 添加匹配的碎片信息
                        'Match_Type': '母离子匹配'  # 标记匹配类型
                    })
            except:
                continue

        return matches

    def find_aglycone_in_fragments(self, fragments):
        """在MS/MS碎片中查找匹配的苷元"""
        matches = []

        for _, row in self.aglycone_df.iterrows():
            try:
                aglycone_mz = float(row['aglycone m/z'])
                # 检查每个碎片是否与苷元匹配
                for frag_mz in fragments.keys():
                    if self.is_match(frag_mz, aglycone_mz):
                        matches.append({
                            'Name': row['Name'],
                            'Formula': row['Formula'],
                            'aglycone_m/z': aglycone_mz,
                            'Fragment_m/z': frag_mz,
                            'Mass_Difference': abs(frag_mz - aglycone_mz),
                            'Matched_Fragment_m/z': frag_mz,
                            'Match_Type': '碎片匹配'  # 标记匹配类型
                        })
                        break  # 找到一个匹配就停止，避免重复
            except:
                continue

        return matches

    def find_all_matching_fragments(self, fragments, target_fragments):
        """
        在碎片中查找所有匹配的特征碎片

        Parameters:
        -----------
        fragments : dict
            碎片字典 {m/z: 强度}
        target_fragments : dict
            目标碎片字典 {羟基个数: m/z值}

        Returns:
        --------
        list
            匹配的碎片列表，每个元素为(羟基个数, 碎片m/z, 强度)
        """
        matches = []

        for oh_count, target_mz in target_fragments.items():
            for frag_mz, intensity in fragments.items():
                if self.is_match(frag_mz, target_mz):
                    matches.append((oh_count, frag_mz, intensity))

        return matches

    def determine_flavonoid_type(self, fragments):
        """
        根据碎片确定黄酮类型

        Parameters:
        -----------
        fragments : dict
            碎片字典 {m/z: 强度}

        Returns:
        --------
        tuple
            (黄酮类型, 羟基个数, 特征碎片m/z, 碎片强度, 置信度, 匹配碎片信息)
        """
        # 查找黄酮特征碎片
        flavonoid_matches = self.find_all_matching_fragments(fragments, self.flavonoid_fragments)
        # 查找黄酮醇特征碎片
        flavonol_matches = self.find_all_matching_fragments(fragments, self.flavonol_fragments)

        # 统计匹配数量
        flavonoid_count = len(flavonoid_matches)
        flavonol_count = len(flavonol_matches)

        # 准备匹配碎片信息
        all_matches_info = []

        # 黄酮匹配信息
        for oh_count, mz, intensity in flavonoid_matches:
            all_matches_info.append(f"黄酮碎片{oh_count}OH:{mz:.5f}({intensity})")

        # 黄酮醇匹配信息
        for oh_count, mz, intensity in flavonol_matches:
            all_matches_info.append(f"黄酮醇碎片{oh_count}OH:{mz:.5f}({intensity})")

        matches_info_str = "; ".join(all_matches_info)

        # 判断逻辑
        if flavonol_count >= 1:
            # 有黄酮醇专属碎片
            if flavonoid_count >= 1:
                # 同时有黄酮通用碎片和黄酮醇专属碎片 -> 高度可信的黄酮醇
                # 选择响应最高的碎片作为特征碎片
                all_matches = flavonol_matches + flavonoid_matches
                best_match = max(all_matches, key=lambda x: x[2])
                oh_count, frag_mz, frag_intensity = best_match
                return "黄酮醇", oh_count, frag_mz, frag_intensity, "高度可信", matches_info_str
            else:
                # 只有黄酮醇专属碎片 -> 可信的黄酮醇
                best_match = max(flavonol_matches, key=lambda x: x[2])
                oh_count, frag_mz, frag_intensity = best_match
                return "黄酮醇", oh_count, frag_mz, frag_intensity, "可信", matches_info_str
        elif flavonoid_count >= 1:
            # 只有黄酮通用碎片 -> 可能为黄酮
            best_match = max(flavonoid_matches, key=lambda x: x[2])
            oh_count, frag_mz, frag_intensity = best_match
            return "黄酮", oh_count, frag_mz, frag_intensity, "可能", matches_info_str
        else:
            # 没有匹配的特征碎片
            return "未知", 0, 0, 0, "未知", ""

    def find_best_fragment(self, fragments, target_fragments):
        """
        在碎片中查找最佳匹配的特征碎片

        Parameters:
        -----------
        fragments : dict
            碎片字典 {m/z: 强度}
        target_fragments : dict
            目标碎片字典 {羟基个数: m/z值}

        Returns:
        --------
        tuple or None
            (羟基个数, 碎片m/z, 强度) 或 None
        """
        best_match = None
        best_intensity = 0

        for oh_count, target_mz in target_fragments.items():
            for frag_mz, intensity in fragments.items():
                if self.is_match(frag_mz, target_mz) and intensity > best_intensity:
                    best_match = (oh_count, frag_mz, intensity)
                    best_intensity = intensity

        return best_match

    def screen_flavonoids(self):
        """筛选黄酮化合物"""
        print("=" * 70)
        print("开始筛选黄酮化合物...")
        print("=" * 70)

        results = []

        # 处理每个峰
        for idx, row in self.peak_df.iterrows():
            # 显示进度
            if (idx + 1) % 100 == 0:
                print(f"  正在处理第 {idx + 1}/{len(self.peak_df)} 个峰...")

            try:
                # 获取峰信息
                peak_mz = float(row['m/z'])
                rt = row['RT']
                msms_str = row['MS/MS'] if 'MS/MS' in row else ''

                # 1. 解析MS/MS
                fragments = self.parse_msms(msms_str)

                # 2. 首先用母离子m/z匹配苷元
                aglycone_matches_from_precursor = self.find_aglycone_match(peak_mz)

                # 3. 在MS/MS碎片中查找匹配的苷元
                aglycone_matches_from_fragments = self.find_aglycone_in_fragments(fragments)

                # 4. 合并所有匹配结果
                all_aglycone_matches = aglycone_matches_from_precursor + aglycone_matches_from_fragments

                # 5. 如果没有匹配到苷元，跳过该峰
                if not all_aglycone_matches:
                    continue

                # 6. 根据新的逻辑确定黄酮类型
                flav_type, oh_count, frag_mz, frag_intensity, confidence, matches_info = self.determine_flavonoid_type(
                    fragments)

                # 7. 为每个匹配的苷元创建结果记录
                for match in all_aglycone_matches:
                    # 检查是否是Aglycone6匹配
                    is_aglycone6 = match['Name'] == 'Aglycone6'

                    if is_aglycone6:
                        print(f"\n发现Aglycone6匹配!")
                        print(f"  峰m/z: {peak_mz}")
                        print(f"  苷元m/z: {match['aglycone_m/z']}")
                        print(f"  匹配类型: {match.get('Match_Type', '未知')}")
                        print(f"  匹配碎片m/z: {match.get('Matched_Fragment_m/z', 'N/A')}")
                        print(f"  质量差异: {match['Mass_Difference']}")

                    result = {
                        'Peak_ID': idx + 1,
                        'Peak_m/z': peak_mz,
                        'RT': rt,
                        'Original_MS/MS': msms_str,
                        'Aglycone_Name': match['Name'],
                        'Aglycone_Formula': match['Formula'],
                        'Aglycone_m/z': match['aglycone_m/z'],
                        'Mass_Difference': match['Mass_Difference'],
                        'Match_Type': match.get('Match_Type', '未知'),  # 添加匹配类型
                        'Flavonoid_Type': flav_type,
                        'Confidence': confidence,  # 添加置信度
                        'OH_Count': oh_count,
                        'Characteristic_Fragment_m/z': frag_mz,
                        'Fragment_Intensity': frag_intensity,
                        'MS/MS_Fragment_Count': len(fragments),
                        'Matched_Fragment_m/z': match.get('Matched_Fragment_m/z', ''),  # 添加匹配的碎片m/z
                        'Matched_Fragments_Info': matches_info  # 添加所有匹配碎片的详细信息
                    }

                    # 如果是Aglycone6，显示详细信息
                    if is_aglycone6:
                        print(f"  已添加到结果: {result['Aglycone_Name']}")
                        print(f"  黄酮类型: {flav_type} ({confidence})")
                        if matches_info:
                            print(f"  匹配碎片: {matches_info}")

                    results.append(result)

            except Exception as e:
                print(f"处理峰 {idx + 1} 时出错: {e}")
                continue

        # 转换为DataFrame
        if results:
            self.results_df = pd.DataFrame(results)
            print(f"\n筛选完成! 共找到 {len(results)} 个匹配结果")

            # 统计信息
            flavonoid_count = len(self.results_df[self.results_df['Flavonoid_Type'] == '黄酮'])
            flavonol_count = len(self.results_df[self.results_df['Flavonoid_Type'] == '黄酮醇'])
            unknown_count = len(self.results_df[self.results_df['Flavonoid_Type'] == '未知'])

            # 置信度统计
            high_confidence = len(self.results_df[self.results_df['Confidence'] == '高度可信'])
            medium_confidence = len(self.results_df[self.results_df['Confidence'] == '可信'])
            low_confidence = len(self.results_df[self.results_df['Confidence'] == '可能'])
            unknown_confidence = len(self.results_df[self.results_df['Confidence'] == '未知'])

            # 匹配类型统计
            precursor_match_count = len(self.results_df[self.results_df['Match_Type'] == '母离子匹配'])
            fragment_match_count = len(self.results_df[self.results_df['Match_Type'] == '碎片匹配'])

            print(f"\n黄酮类型分布:")
            print(f"  黄酮: {flavonoid_count} 个")
            print(f"  黄酮醇: {flavonol_count} 个")
            print(f"  未知: {unknown_count} 个")

            print(f"\n置信度分布:")
            print(f"  高度可信: {high_confidence} 个")
            print(f"  可信: {medium_confidence} 个")
            print(f"  可能: {low_confidence} 个")
            print(f"  未知: {unknown_confidence} 个")

            print(f"\n匹配方式分布:")
            print(f"  母离子匹配: {precursor_match_count} 个")
            print(f"  碎片匹配: {fragment_match_count} 个")

            # 检查Aglycone6是否在结果中
            if 'Aglycone6' in self.results_df['Aglycone_Name'].values:
                aglycone6_results = self.results_df[self.results_df['Aglycone_Name'] == 'Aglycone6']
                print(f"\nAglycone6匹配结果: {len(aglycone6_results)} 个")
                for i, result in aglycone6_results.iterrows():
                    print(f"  结果 {i + 1}: 峰m/z={result['Peak_m/z']}, 匹配类型={result['Match_Type']}, "
                          f"黄酮类型={result['Flavonoid_Type']}({result['Confidence']}), "
                          f"匹配碎片={result['Matched_Fragment_m/z']}")
        else:
            print("\n未找到任何匹配的化合物")
            self.results_df = pd.DataFrame()

        return len(results) > 0

    def save_results(self):
        """保存结果到文件（仅保存完整结果）"""
        print("\n" + "=" * 70)
        print("保存结果文件...")
        print("=" * 70)

        if not hasattr(self, 'results_df') or self.results_df.empty:
            print("没有结果数据可保存")
            return False

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # 仅保存完整结果
        full_results_file = os.path.join(self.output_dir, f"flavonoid_results_full_{timestamp}.csv")
        try:
            # 重新排列列的顺序
            columns_order = [
                'Peak_ID', 'Peak_m/z', 'RT', 'Original_MS/MS',
                'Aglycone_Name', 'Aglycone_Formula', 'Aglycone_m/z',
                'Match_Type', 'Matched_Fragment_m/z', 'Mass_Difference',
                'Flavonoid_Type', 'Confidence', 'OH_Count',
                'Characteristic_Fragment_m/z', 'Fragment_Intensity',
                'Matched_Fragments_Info', 'MS/MS_Fragment_Count'
            ]

            # 确保所有列都存在
            available_columns = [col for col in columns_order if col in self.results_df.columns]
            remaining_columns = [col for col in self.results_df.columns if col not in available_columns]

            # 重新排列DataFrame
            sorted_df = self.results_df[available_columns + remaining_columns]

            sorted_df.to_csv(full_results_file, index=False, encoding='utf-8-sig')
            print(f"1. 完整结果已保存: {full_results_file}")
            print(f"   包含 {len(sorted_df)} 条记录")
        except Exception as e:
            print(f"保存完整结果时出错: {e}")
            return False

        print(f"\n结果文件已保存到: {self.output_dir}")
        return True

    def run(self):
        """运行完整的筛选流程"""
        print("=" * 70)
        print("黄酮化合物筛选程序")
        print("=" * 70)
        print(f"质量容差: ±{self.mass_tolerance}")
        print(f"输出文件夹: {self.output_dir}")
        print("=" * 70 + "\n")

        # 1. 加载数据
        if not self.load_data():
            print("数据加载失败，程序终止")
            return False

        # 2. 筛选黄酮
        if not self.screen_flavonoids():
            print("筛选过程失败或未找到匹配")
            return False

        # 3. 保存结果
        if not self.save_results():
            print("结果保存失败")
            return False

        return True


# 主程序
if __name__ == "__main__":
    print("=" * 70)
    print("黄酮化合物筛选程序 v3.0 (增强黄酮醇识别)")
    print("=" * 70)

    try:
        # 创建筛选器实例
        screener = FlavonoidScreener(mass_tolerance=0.01)

        # 运行筛选
        success = screener.run()

        if success:
            print("\n程序执行成功!")
        else:
            print("\n程序执行失败!")

    except Exception as e:
        print(f"\n程序执行时发生错误: {e}")
        import traceback
        traceback.print_exc()

    # 保持窗口打开
    input("\n按Enter键退出...")