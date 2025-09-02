import os
import re
import pandas as pd

# 定义函数：从单个GTF文件提取转录本表达量（FPKM或FPKM）
def extract_expression(gtf_file, metric="FPKM"):
    expr_dict = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # 跳过注释行
                continue
            cols = line.strip().split('\t')
            if cols[2] == 'transcript':  # 只处理转录本行
                attr = cols[8]
                # 提取transcript_id
                tid_match = re.search(r'transcript_id "([^"]+)"', attr)
                # 提取指定的表达量指标（FPKM或FPKM）
                metric_match = re.search(f'{metric} "([^"]+)"', attr)
                if tid_match and metric_match:
                    tid = tid_match.group(1)
                    value = float(metric_match.group(1))
                    expr_dict[tid] = value
    return expr_dict

# 批量处理所有GTF文件
gtf_dir = "./5B-stringtie/"  # GTF文件所在目录
metric = "FPKM"           # 选择提取FPKM或FPKM
gtf_files = [f for f in os.listdir(gtf_dir) if f.endswith('.gtf')]

# 初始化表达量矩阵（字典的字典）
expr_matrix = {}
for gtf in gtf_files:
    sample_name = gtf.replace('.gtf', '')  # 样本名（从文件名提取）
    expr_dict = extract_expression(os.path.join(gtf_dir, gtf), metric)
    # 写入矩阵
    for tid, value in expr_dict.items():
        if tid not in expr_matrix:
            expr_matrix[tid] = {}
        expr_matrix[tid][sample_name] = value

# 转换为DataFrame并保存为CSV
df = pd.DataFrame.from_dict(expr_matrix, orient='index').fillna(0)  # 缺失值填0
df.to_csv(f'transcript_{metric}_matrix.csv')
print(f"已生成转录本{metric}表达矩阵：transcript_{metric}_matrix.csv")