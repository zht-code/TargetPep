# import json
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd

# # 读取数据
# file_path = '/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json'
# file_path1 = '/root/autodl-tmp/PP_generate_v1/data/output/quvina_final_hdock.json'

# with open(file_path, 'r') as f:
#     data = json.load(f)
# with open(file_path1, 'r') as f:
#     data1 = json.load(f)

# # 提取亲和力数据
# affinities = {
#     #"natural peptide": [],
#     "generated": [],
#     "quvina": []
# }

# for key, value in data.items():
#     for prop in value['properties']:
#         #affinities["natural peptide"].append(prop["test Affinity (kcal/mol)"])
#         affinities["generated"].append(prop["generate Affinity (kcal/mol)"])
# for key, value in data1.items():
#     for prop in value['properties']:
#         affinities["quvina"].append(prop["generate Affinity (kcal/mol)"])

# df = pd.DataFrame(affinities)

# # 画小提琴图
# plt.figure(figsize=(10, 6))
# #order = ["natural peptide", "generated", "quvina"]
# order = [ "generated", "quvina"]
# ax = sns.violinplot(data=df, order=order)

# # ✨计算均值和中位数并标注
# for i, model in enumerate(order):
#     mean_value = df[model].mean()
#     median_value = df[model].median()
#     # 标注均值
#     ax.text(i, mean_value + 0.5, f"Mean: {mean_value:.2f}", 
#             horizontalalignment='center', size='small', color='black', weight='semibold')
#     # # 标注中位数
#     # ax.text(i, median_value - 1.0, f"Median: {median_value:.2f}", 
#     #         horizontalalignment='center', size='small', color='blue', weight='semibold')

# # 设置标题和标签
# plt.title("Affinity Comparison Across Different Models")
# plt.ylabel("Affinity (kcal/mol)")
# plt.xlabel("Models")

# # 保存并显示
# plt.tight_layout()
# plt.savefig("/root/autodl-tmp/PP_generate_v1/utils/evaluate/quvina_violin3.png", dpi=300)
# plt.show()

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 替换为你的文件路径
peptides_file_path = '/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json'
quvina_file_path1 = '/root/autodl-tmp/PP_generate_v1/data/output/quvina_final_hdock.json'

# 读取亲和力数据
def load_affinity_from_json(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    affinities = []
    for entry in data.values():
        try:
            affinity = entry["properties"][0]["generate Affinity (kcal/mol)"]
            affinities.append(affinity)  # 取绝对值用于可视化
        except (KeyError, IndexError):
            continue
    return affinities

# 每5个取平均
def group_average(data, group_size=10):
    return [np.mean(data[i:i + group_size]) for i in range(0, len(data), group_size)]

# 加载数据
peptide_affinity = load_affinity_from_json(peptides_file_path)
quvina_affinity = load_affinity_from_json(quvina_file_path1)

# 每5个取平均
peptide_avg = group_average(peptide_affinity)
quvina_avg = group_average(quvina_affinity)

# 绘图
plt.figure(figsize=(6, 4))
sns.kdeplot(quvina_avg, fill=True, label="quvina", color='green', alpha=0.3)
sns.kdeplot(peptide_avg, fill=True, label="TargetPep", color='orange', alpha=0.6)

plt.legend()
# plt.xlim(6, 31)
# plt.ylim(0, 35)
plt.tight_layout()
plt.savefig("/root/autodl-tmp/PP_generate_v1/data/Ablation/ablation.png", dpi=300)
plt.show()
