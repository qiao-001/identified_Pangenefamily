import pandas as pd
import os
input_dir = "./5stringtie"  # StringTie输出文件的目录
#input_dir = "/mnt/e/24.4.-转录组数据备份/0-Ara测试/7ncRNAtab"  # StringTie输出文件的目录#
sample_files = [f for f in os.listdir(input_dir) if f.endswith(".tab")]  # 获取所有TSV文件
def process_sample(file_path, sample_name):
    # 读取TSV文件，假设列名包含 'Gene_ID' 和 'TPM'
    df = pd.read_csv(file_path, sep='\t', comment='#')  # 跳过注释行（以#开头）
    
    # 处理重复基因：按Gene_ID分组，对TPM求和（可选均值、最大值等）
    df_agg = df.groupby('Gene ID', as_index=False)['TPM'].sum()
    
    # 重命名TPM列为样本名
    df_agg.rename(columns={'TPM': sample_name}, inplace=True)
    return df_agg

# 初始化表达矩阵（以第一个样本为基准）
merged_df = None

for file in sample_files:
    file_path = os.path.join(input_dir, file)
    sample_name = os.path.splitext(file)[0]  # 从文件名提取样本名（例如 "sample1.tsv" → "sample1"）
    
    # 处理当前样本
    df_sample = process_sample(file_path, sample_name)
    
    # 逐步合并
    if merged_df is None:
        merged_df = df_sample
    else:
        merged_df = pd.merge(merged_df, df_sample, on='Gene ID', how='outer')

# 填充缺失值为0
merged_df.fillna(0, inplace=True)
output_file = "TPM_expression_matrix.csv"
merged_df.to_csv(output_file, index=False)
print(f"表达矩阵已保存至 {output_file}")