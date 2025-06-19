
# import os
# import glob
# import json
# import freesasa
# from Bio.PDB import PDBParser
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd
# import traceback
# '''绘制小提琴图'''
# # 提取PDB中氨基酸
# def get_amino_acids_from_pdb(pdb_file):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("peptide", pdb_file)
#     amino_acids = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.has_id("CA"):
#                     amino_acids.append(residue.get_resname())
#     return amino_acids

# # 计算SASA
# def calculate_sasa_with_freesasa(pdb_file):
#     structure = freesasa.Structure(pdb_file)
#     result = freesasa.calc(structure)
#     total_sasa = result.totalArea()
#     return total_sasa

# # 计算平均疏水性
# def calculate_hydrophobicity(amino_acids):
#     hydrophobicity = {
#         'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
#         'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
#         'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
#         'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
#     }
#     total_hydro = sum(hydrophobicity.get(aa, 0) for aa in amino_acids)
#     return total_hydro / len(amino_acids) if amino_acids else 0

# # 溶解性预测
# def predict_solubility(hydrophobicity, sasa):
#     return sasa - (hydrophobicity * 100)

# # 处理一个文件夹
# def process_folder(folder_path, ids, label_name):
#     results = []
#     for protein_id in ids:
#         pdb_filename = f"{protein_id}_unrelaxed_rank_001*.pdb"
#         pdb_path = os.path.join(folder_path, pdb_filename)
#         ligand_files = glob.glob(pdb_path)
#         ligand_pdb_file = ligand_files[0]
#         if not os.path.exists(ligand_pdb_file):
#             print(f"Missing PDB file for {protein_id} in {label_name}")
#             continue
#         try:
#             amino_acids = get_amino_acids_from_pdb(ligand_pdb_file)
#             sasa = calculate_sasa_with_freesasa(ligand_pdb_file)
#             hydro = calculate_hydrophobicity(amino_acids)
#             solubility = predict_solubility(hydro, sasa)
#             results.append({
#                 "protein_id": protein_id,
#                 "solubility_score": solubility,
#                 "source": label_name
#             })
#         except Exception as e:
#             print(f"Error processing {ligand_pdb_file}: {e}\n{traceback.format_exc()}")
#     return results

# # 画小提琴图
# def plot_violin(data, output_path="/root/autodl-tmp/PP_generate_v1/utils/evaluate/solubility_violin_plot.png"):
#     plt.figure(figsize=(8, 6))
#     sns.violinplot(x="source", y="solubility_score", data=data, inner="box")
    
#     # 计算均值并标注
#     means = data.groupby('source')['solubility_score'].mean()
#     ax = plt.gca()
#     for i, (label, mean_value) in enumerate(means.items()):
#         ax.text(i, mean_value + 50, f"{mean_value:.1f}", 
#                 horizontalalignment='center', size='medium', color='black', weight='semibold')
    
#     plt.xlabel('Source')
#     plt.ylabel('Solubility Score')
#     plt.title('Comparison of Solubility Scores: cross vs qusolubility')
#     plt.tight_layout()

#     # 自动创建输出目录
#     os.makedirs(os.path.dirname(output_path), exist_ok=True)
#     plt.savefig(output_path, dpi=300)
#     plt.show()

# if __name__ == "__main__":
#     # 配置你的路径
#     json_file = "/root/autodl-tmp/PP_generate_v1/data/output/Hdock_generate_vina_cross1.json"  # 改成你的id json路径
#     cross_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/cross"
#     qusolubility_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/qusolubility"

#     # 读取ids
#     with open(json_file, "r") as f:
#         id_data = json.load(f)
#     ids = list(id_data.keys())

#     # 分别处理两个文件夹
#     cross_data = process_folder(cross_folder, ids, "generate")
#     qusolubility_data = process_folder(qusolubility_folder, ids, "qusolubility")

#     all_data = pd.DataFrame(cross_data + qusolubility_data)
#     print(f"Collected {len(all_data)} solubility scores.")

#     plot_violin(all_data)




# import os
# import glob
# import json
# import freesasa
# from Bio.PDB import PDBParser
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd
# import traceback
# '''绘制散点图'''
# # 提取PDB中氨基酸
# def get_amino_acids_from_pdb(pdb_file):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("peptide", pdb_file)
#     amino_acids = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.has_id("CA"):
#                     amino_acids.append(residue.get_resname())
#     return amino_acids

# # 计算SASA
# def calculate_sasa_with_freesasa(pdb_file):
#     structure = freesasa.Structure(pdb_file)
#     result = freesasa.calc(structure)
#     total_sasa = result.totalArea()
#     return total_sasa

# # 计算平均疏水性
# def calculate_hydrophobicity(amino_acids):
#     hydrophobicity = {
#         'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
#         'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
#         'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
#         'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
#     }
#     total_hydro = sum(hydrophobicity.get(aa, 0) for aa in amino_acids)
#     return total_hydro / len(amino_acids) if amino_acids else 0

# # 溶解性预测
# def predict_solubility(hydrophobicity, sasa):
#     return sasa - (hydrophobicity * 100)

# # 处理一个文件夹
# def process_folder(folder_path, ids, label_name):
#     results = []
#     for protein_id in ids:
#         pdb_filename = f"{protein_id}_unrelaxed_rank_001*.pdb"
#         pdb_path = os.path.join(folder_path, pdb_filename)
#         ligand_files = glob.glob(pdb_path)
#         if not ligand_files:
#             print(f"Missing PDB file for {protein_id} in {label_name}")
#             continue
#         ligand_pdb_file = ligand_files[0]
#         try:
#             amino_acids = get_amino_acids_from_pdb(ligand_pdb_file)
#             sasa = calculate_sasa_with_freesasa(ligand_pdb_file)
#             hydro = calculate_hydrophobicity(amino_acids)
#             solubility = predict_solubility(hydro, sasa)
#             results.append({
#                 "protein_id": protein_id,
#                 "solubility_score": solubility,
#                 "hydrophobicity": hydro,
#                 "source": label_name
#             })
#         except Exception as e:
#             print(f"Error processing {ligand_pdb_file}: {e}\n{traceback.format_exc()}")
#     return results

# # 对solubility_score进行归一化 (新增函数)
# def normalize_solubility(df):
#     min_score = df['solubility_score'].min()
#     max_score = df['solubility_score'].max()
#     if max_score - min_score == 0:
#         df['normalized_solubility_score'] = 0.5  # 防止除0，全部设为中间值
#     else:
#         df['normalized_solubility_score'] = (df['solubility_score'] - min_score) / (max_score - min_score)
#     return df

# # 画小提琴图（修改：画归一化后的）
# # def plot_violin(data, output_path="/root/autodl-tmp/PP_generate_v1/utils/evaluate/solubility_violin_plot1.png"):
# #     plt.figure(figsize=(8, 6))
# #     sns.violinplot(x="source", y="normalized_solubility_score", data=data, inner="box", palette={"generate": "orange", "qusolubility": "green"})
# #     # sns.violinplot(
# #     #     x="source",
# #     #     y="normalized_solubility_score",
# #     #     data=data,
# #     #     inner="box",
# #     #     }  # 添加颜色映射
# #     # )

# #     # 计算均值并标注
# #     means = data.groupby('source')['normalized_solubility_score'].mean()
# #     ax = plt.gca()
# #     for i, (label, mean_value) in enumerate(means.items()):
# #         ax.text(i, mean_value + 0.02, f"{mean_value:.2f}", 
# #                 horizontalalignment='left', size='medium', color='black', weight='semibold')
    
# #     plt.xlabel('Source')
# #     plt.ylabel('Normalized Solubility Score (0-1)')
# #     plt.title('Comparison of Normalized Solubility Scores: cross vs qusolubility')
# #     plt.tight_layout()

# #     # 自动创建输出目录
# #     os.makedirs(os.path.dirname(output_path), exist_ok=True)
# #     plt.savefig(output_path, dpi=300)
# #     plt.show()
# # 绘制散点图
# def plot_scatter(data, output_path="/root/autodl-tmp/PP_generate_v1/utils/evaluate/solubility_scatter_plot.png"):
#     plt.figure(figsize=(8, 6))
#     sns.scatterplot(
#         data=data,
#         x="hydrophobicity",
#         y="normalized_solubility_score",
#         hue="source",
#         palette={"generate": "orange", "qusolubility": "green"},
#         alpha=0.7,
#         s=50  # 点的大小
#     )
#     plt.xlabel("Average Hydrophobicity")
#     plt.ylabel("Solubility Score")
#     plt.title("Hydrophobicity vs Solubility Score")
#     plt.tight_layout()
#     os.makedirs(os.path.dirname(output_path), exist_ok=True)
#     plt.savefig(output_path, dpi=300)
#     plt.show()



# if __name__ == "__main__":
#     # 配置你的路径
#     json_file = "/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json"  # 改成你的id json路径
#     cross_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/peptides_final1"
#     qusolubility_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/qusolubility_final"

#     # 读取ids
#     with open(json_file, "r") as f:
#         id_data = json.load(f)
#     ids = list(id_data.keys())

#     # 分别处理两个文件夹
#     cross_data = process_folder(cross_folder, ids, "generate")
#     qusolubility_data = process_folder(qusolubility_folder, ids, "qusolubility")

#     # 合并数据
#     all_data = pd.DataFrame(cross_data + qusolubility_data)
#     print(f"Collected {len(all_data)} solubility scores.")

#     # 归一化处理
#     all_data = normalize_solubility(all_data)

#     # 绘图
#     plot_scatter(all_data)






import os
import glob
import json
import freesasa
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import traceback
'''绘制散点图'''
# 提取PDB中氨基酸
def get_amino_acids_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("peptide", pdb_file)
    amino_acids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    amino_acids.append(residue.get_resname())
    return amino_acids

# 计算SASA
def calculate_sasa_with_freesasa(pdb_file):
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    total_sasa = result.totalArea()
    return total_sasa

# 计算平均疏水性
def calculate_hydrophobicity(amino_acids):
    hydrophobicity = {
        'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
        'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
        'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
        'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
    }
    total_hydro = sum(hydrophobicity.get(aa, 0) for aa in amino_acids)
    return total_hydro / len(amino_acids) if amino_acids else 0

# 溶解性预测
def predict_solubility(hydrophobicity, sasa):
    return sasa - (hydrophobicity * 100)

# 处理一个文件夹
def process_folder(folder_path, ids, label_name):
    results = []
    for protein_id in ids:
        pdb_filename = f"{protein_id}_unrelaxed_rank_001*.pdb"
        pdb_path = os.path.join(folder_path, pdb_filename)
        ligand_files = glob.glob(pdb_path)
        if not ligand_files:
            print(f"Missing PDB file for {protein_id} in {label_name}")
            continue
        ligand_pdb_file = ligand_files[0]
        try:
            amino_acids = get_amino_acids_from_pdb(ligand_pdb_file)
            sasa = calculate_sasa_with_freesasa(ligand_pdb_file)
            hydro = calculate_hydrophobicity(amino_acids)
            solubility = predict_solubility(hydro, sasa)
            results.append({
                "protein_id": protein_id,
                "solubility_score": solubility,
                "hydrophobicity": hydro,
                "source": label_name
            })
        except Exception as e:
            print(f"Error processing {ligand_pdb_file}: {e}\n{traceback.format_exc()}")
    return results

# 对solubility_score进行归一化 (新增函数)
def normalize_solubility(df):
    min_score = df['solubility_score'].min()
    max_score = df['solubility_score'].max()
    if max_score - min_score == 0:
        df['normalized_solubility_score'] = 0.5  # 防止除0，全部设为中间值
    else:
        df['normalized_solubility_score'] = (df['solubility_score'] - min_score) / (max_score - min_score)
    return df

# 画对比散点图（修改：画归一化后的）
def plot_method_comparison_scatter(df1, df2, output_path="/root/autodl-tmp/PP_generate_v1/data/Ablation/compare_methods_scatter.png"):
    # 设置索引
    df1 = pd.DataFrame(df1).set_index("protein_id")
    df2 = pd.DataFrame(df2).set_index("protein_id")

    # 重命名列
    df1 = df1.rename(columns={"normalized_solubility_score": "score_generate"})
    df2 = df2.rename(columns={"normalized_solubility_score": "score_qusolubility"})

    # 合并
    df_merged = pd.merge(df1, df2, left_index=True, right_index=True)

    # 分类颜色标签
    def compare(row):
        if row["score_generate"] > row["score_qusolubility"]:
            return "generate_better"
        elif row["score_generate"] < row["score_qusolubility"]:
            return "qusolubility_better"
        else:
            return "equal"

    df_merged["comparison"] = df_merged.apply(compare, axis=1)

    # 绘图
    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=df_merged,
        x="score_generate",
        y="score_qusolubility",
        hue="comparison",
        palette={
            "generate_better": "orange",
            "qusolubility_better": "green",
            "equal": "gray"
        },
        alpha=0.7,
        s=50
    )
    # 统计每类点的数量
    gen_count = (df_merged["comparison"] == "generate_better").sum()
    quso_count = (df_merged["comparison"] == "qusolubility_better").sum()

    # 添加文字到图中显示（右下角，避免和 legend 冲突）
    plt.text(
        0.99, 0.01,
        f"Generated better: {gen_count}\nQusolubility better: {quso_count}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
    )

    plt.plot([0, 1], [0, 1], 'r--', linewidth=1, label="y = x")
    plt.xlabel("Generate Method (Normalized Solubility)")
    plt.ylabel("Qusolubility Method (Normalized Solubility)")
    plt.title("Method Comparison: Generate vs Qusolubility")
    plt.legend(title="Which is better")
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.show()



if __name__ == "__main__":
    # 配置你的路径
    json_file = "/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json"  # 改成你的id json路径
    cross_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/peptides_final1"
    qusolubility_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/qusolubility_final"

    # 读取ids
    with open(json_file, "r") as f:
        id_data = json.load(f)
    ids = list(id_data.keys())

    # 分别处理两个文件夹
    cross_data = process_folder(cross_folder, ids, "generate")
    qusolubility_data = process_folder(qusolubility_folder, ids, "qusolubility")

    # 合并数据
    cross_df = pd.DataFrame(cross_data)
    qusol_df = pd.DataFrame(qusolubility_data)

    # 归一化处理
    cross_df = normalize_solubility(cross_df)
    qusol_df = normalize_solubility(qusol_df)

    # 对比散点图
    plot_method_comparison_scatter(cross_df.to_dict(orient="records"), qusol_df.to_dict(orient="records"))
