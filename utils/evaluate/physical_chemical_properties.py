# import json
# import re
# from Bio.SeqUtils.ProtParam import ProteinAnalysis
# from collections import Counter
# '''计算多条物理化学性质数值'''
# # 读取 JSON 文件
# with open('/root/autodl-tmp/PP_generate_v1/data/output/RFdiffusion_sequence.json', 'r') as file:
#     data = json.load(file)

# # 用于存储计算结果的字典
# results = {}

# # 遍历 JSON 数据中的每个条目
# for key, value in data.items():
#     clean_peptide = value['sequence']
#     # 移除序列中的数字
#     # clean_peptide = re.sub(r'\d', '', peptides[0])  # 使用正则表达式删除数字
#     print('多肽序列：', clean_peptide)
#     peptide_properties = []

#     # 使用 ProteinAnalysis 对象分析序列
#     prot_analysis = ProteinAnalysis(clean_peptide)

#     # 计算艾森伯格疏水性（使用平均疏水性代替）
#     hydrophobicity = prot_analysis.gravy()
#     # 计算艾森伯格疏水矩（假设用氨基酸的疏水性偏差表示）
#     eisenberg_scale = {
#         'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
#         'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
#         'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
#         'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
#     }
#     window_size = 7
#     # hydrophobic_moment = 

#     hydrophobic_moment = sum([abs(h - hydrophobicity) for h in prot_analysis.protein_scale(eisenberg_scale, window_size)]) / len(clean_peptide)

#     # 计算氨基酸频率
#     amino_acid_frequency = prot_analysis.get_amino_acids_percent()

#     # 计算总电荷
#     total_charge = prot_analysis.charge_at_pH(7.4)

#     # 计算序列长度
#     sequence_length = len(clean_peptide)

#     # 存储每个多肽的计算结果
#     peptide_properties.append({
#         'Hydrophobicity (Eisenberg)': hydrophobicity,
#         'Hydrophobic moment (Eisenberg)': hydrophobic_moment,
#         'Amino acid frequency': amino_acid_frequency,
#         'Total charge at pH 7.4': total_charge,
#         'Sequence length': sequence_length
#     })

#     # 将每个条目的结果添加到结果字典中
#     results[key] = {
#         # 'protein_sequence': value['protein_sequence'],
#         'generated_peptide': clean_peptide,
#         'properties': peptide_properties
#     }

# # 将结果字典保存到新的 JSON 文件
# with open('/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/RFdiffusion_physical.json', 'w') as outfile:
#     json.dump(results, outfile, indent=4)



import json
import numpy as np
import matplotlib.pyplot as plt

# ——— 1. 文件与标签 —————————————————————————————
json_files = [
    '/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/test_physical.json',
    '/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/peptides_physical.json',
    '/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/ppflow_physical.json',
    '/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/RFdiffusion_physical.json',
    # '/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/quvina_physical.json'
]
labels = ['test', 'peptide', 'ppflow','RFdiffusion']
# labels = ['test', 'peptide']

# ——— 2. 参考示例图的配色（Uniprot 粉／Real-AMPs 橙／GPT 黄）—————
# Uniprot（粉色）：#e8767e
# Real-AMPs（橙色）：#f1a21d
# GPT（黄色）：#ffd324
# Prompt（蓝色）：#4da6ff
# Prompt-TopK（青绿色）：#17becf
# Distillation（紫色）：#9467bd
colors = {
    'test':         '#fac55d',  # 对应 Uniprot 粉
    'peptide':      '#745cac',  # 对应 Real-AMPs 橙
    'ppflow':       '#366b82',  # 对应 GPT 黄
    'RFdiffusion':  '#3570b3',  # 对应 GPT 黄
}


# ——— 3. 读取并整理数据 —————————————————————————————
aa_list = 'ACDEFGHIKLMNPQRSTVWY'
aa_freqs = {lbl: [] for lbl in labels}
charges = {lbl: [] for lbl in labels}
hydrophobicity = {lbl: [] for lbl in labels}
hydrophobic_moment = {lbl: [] for lbl in labels}
seq_lengths = {lbl: [] for lbl in labels}

for fp, lbl in zip(json_files, labels):
    with open(fp, 'r') as f:
        data = json.load(f)
    for pid, rec in data.items():
        p = rec['properties'][0]
        # 氨基酸频率
        aa_freqs[lbl].append([p['Amino acid frequency'][aa] for aa in aa_list])
        # 总电荷
        charges[lbl].append(p['Total charge at pH 7.4'])
        # 疏水性与疏水矩
        hydrophobicity[lbl].append(p['Hydrophobicity (Eisenberg)'])
        hydrophobic_moment[lbl].append(p['Hydrophobic moment (Eisenberg)'])
        # 序列长度
        seq_lengths[lbl].append(p['Sequence length'])

# ——— 4. 绘图 ————————————————————————————————————————

# A. 氨基酸频率条形图（采用平均值）
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(aa_list))
width = 0.25
for i, lbl in enumerate(labels):
    ax.bar(x + i*width,
           np.mean(aa_freqs[lbl], axis=0),
           width=width,
           color=colors[lbl],
           label=lbl)
ax.set_xticks(x + width*(len(labels)-1)/2)
ax.set_xticklabels(list(aa_list))
ax.set_xlabel('Amino acids')
ax.set_ylabel('Average frequency')
ax.set_title('Amino Acid Frequencies')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig('/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/Frequencies.png', dpi=300)
plt.close(fig)


# B. 总电荷分布直方图
fig, ax = plt.subplots(figsize=(10, 6))
bins = np.arange(-5, 21, 1)
for lbl in labels:
    ax.hist(charges[lbl],
            bins=bins,
            density=True,
            alpha=0.6,
            color=colors[lbl],
            label=lbl)
ax.set_xlabel('Global charge (pH 7.4)')
ax.set_ylabel('Density')
ax.set_title('Global Charge Distribution')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig('/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/Charge.png', dpi=300)
plt.close(fig)


# C. 序列长度柱状图图（whis=(0,100) 强制展示最小/最大）
means = [np.mean(seq_lengths[lbl]) for lbl in labels]

plt.figure(figsize=(8, 6))
bars = plt.bar(labels, means, color=[colors[lbl] for lbl in labels], width=0.2)
# 如果某组全为同一个值，可以加文本说明
for i, lbl in enumerate(labels):
    vals = seq_lengths[lbl]
    if len(set(vals)) == 1:
        plt.text(i, means[i], f"All={vals[0]}", ha='center', color='red')

plt.ylabel('Sequence length')
plt.title('Sequence Length Mean Value')
plt.tight_layout()
plt.show()
plt.savefig('/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/Sequence_Length.png', dpi=300)

# D. Eisenberg 疏水性小提琴图（含中位线/极值线）
fig, ax = plt.subplots(figsize=(8, 6))
data_h = [hydrophobicity[lbl] for lbl in labels]
vp = ax.violinplot(data_h,
                   showmedians=True,
                   showextrema=True)
for i, lbl in enumerate(labels):
    vp['bodies'][i].set_facecolor(colors[lbl])
    vp['bodies'][i].set_edgecolor('black')
    vp['bodies'][i].set_alpha(0.7)
# 设置中位线和极值线颜色
for key in ('cmedians','cmins','cmaxes'):
    plt.setp(vp[key], color='black')
ax.set_xticks(np.arange(1, len(labels)+1))
ax.set_xticklabels(labels)
ax.set_ylabel('Eisenberg hydrophobicity')
ax.set_title('Eisenberg Hydrophobicity')
fig.tight_layout()
fig.savefig('/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/Hydrophobicity_Violin.png', dpi=300)
plt.close(fig)


# E. Eisenberg 疏水矩小提琴图（含中位线/极值线）
fig, ax = plt.subplots(figsize=(8, 6))
data_m = [hydrophobic_moment[lbl] for lbl in labels]
vp2 = ax.violinplot(data_m,
                    showmedians=True,
                    showextrema=True)
for i, lbl in enumerate(labels):
    vp2['bodies'][i].set_facecolor(colors[lbl])
    vp2['bodies'][i].set_edgecolor('black')
    vp2['bodies'][i].set_alpha(0.7)
for key in ('cmedians','cmins','cmaxes'):
    plt.setp(vp2[key], color='black')
ax.set_xticks(np.arange(1, len(labels)+1))
ax.set_xticklabels(labels)
ax.set_ylabel('Eisenberg hydrophobic moment')
ax.set_title('Eisenberg Hydrophobic Moment')
fig.tight_layout()
fig.savefig('/root/autodl-tmp/PP_generate_v1/data/physical_chemical_properties/Hydrophobicity_Moment_Violin.png', dpi=300)
plt.close(fig)

print("所有图表已生成并保存到当前目录。")
