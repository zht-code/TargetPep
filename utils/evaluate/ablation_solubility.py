import os
import glob
import json
import freesasa
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import traceback

'''绘制散点图 + 使用 GRAVY 作为溶解性指标'''

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

# 计算GRAVY指数（平均疏水性）
def calculate_gravy(amino_acids):
    hydrophobicity = {
        'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
        'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
        'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
        'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
    }
    total_hydro = sum(hydrophobicity.get(aa, 0) for aa in amino_acids)
    return total_hydro / len(amino_acids) if amino_acids else 0

# 使用 GRAVY 指数直接作为溶解性得分（负相关）
def predict_solubility_by_gravy(gravy_score):
    return -gravy_score  # 越小越亲水，越可溶；取负变成正向评分

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
            gravy = calculate_gravy(amino_acids)
            solubility = predict_solubility_by_gravy(gravy)
            results.append({
                "protein_id": protein_id,
                "solubility_score": solubility,
                "gravy_score": gravy,
                "source": label_name
            })
        except Exception as e:
            print(f"Error processing {ligand_pdb_file}: {e}\n{traceback.format_exc()}")
    return results

# 对 solubility_score 归一化
def normalize_solubility(df):
    min_score = df['solubility_score'].min()
    max_score = df['solubility_score'].max()
    if max_score - min_score == 0:
        df['normalized_solubility_score'] = 0.5
    else:
        df['normalized_solubility_score'] = (df['solubility_score'] - min_score) / (max_score - min_score)
    return df

# 绘图
def plot_method_comparison_scatter(df1, df2, output_path="/root/autodl-tmp/PP_generate_v1/data/Ablation/compare_methods_scatter1.png"):
    df1 = pd.DataFrame(df1).set_index("protein_id")
    df2 = pd.DataFrame(df2).set_index("protein_id")

    df1 = df1.rename(columns={"normalized_solubility_score": "score_TargetPep"})
    df2 = df2.rename(columns={"normalized_solubility_score": "score_qusolubility"})

    df_merged = pd.merge(df1, df2, left_index=True, right_index=True)

    def compare(row):
        if row["score_TargetPep"] > row["score_qusolubility"]:
            return "TargetPep_better"
        elif row["score_TargetPep"] < row["score_qusolubility"]:
            return "qusolubility_better"
        else:
            return "equal"

    df_merged["comparison"] = df_merged.apply(compare, axis=1)

    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=df_merged,
        x="score_TargetPep",
        y="score_qusolubility",
        hue="comparison",
        palette={
            "TargetPep_better": "orange",
            "qusolubility_better": "green",
            "equal": "gray"
        },
        alpha=0.7,
        s=50
    )
    gen_count = (df_merged["comparison"] == "TargetPep_better").sum()
    quso_count = (df_merged["comparison"] == "qusolubility_better").sum()

    plt.text(
        0.99, 0.01,
        f"TargetPep better: {gen_count}\nQusolubility better: {quso_count}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
    )

    plt.plot([0, 1], [0, 1], 'r--', linewidth=1, label="y = x")
    plt.xlabel("TargetPep Method (Normalized Solubility)")
    plt.ylabel("Qusolubility Method (Normalized Solubility)")
    plt.title("Method Comparison (by GRAVY-derived Solubility)")
    plt.legend(title="Which is better")
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.show()

# 主执行逻辑
if __name__ == "__main__":
    json_file = "/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json"
    cross_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/peptides_final1"
    qusolubility_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/qusolubility_final"

    with open(json_file, "r") as f:
        id_data = json.load(f)
    ids = list(id_data.keys())

    cross_data = process_folder(cross_folder, ids, "generate")
    qusolubility_data = process_folder(qusolubility_folder, ids, "qusolubility")

    cross_df = pd.DataFrame(cross_data)
    qusol_df = pd.DataFrame(qusolubility_data)

    cross_df = normalize_solubility(cross_df)
    qusol_df = normalize_solubility(qusol_df)

    plot_method_comparison_scatter(cross_df.to_dict(orient="records"), qusol_df.to_dict(orient="records"))
