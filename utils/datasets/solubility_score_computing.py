import freesasa
from Bio.PDB import PDBParser
import os
import json
import traceback

# 从 PDB 文件中提取氨基酸序列
def get_amino_acids_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("peptide", pdb_file)
    amino_acids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):  # 确定是氨基酸
                    amino_acids.append(residue.get_resname())  # 获取氨基酸名称
    return amino_acids

# 使用 FreeSASA 计算 SASA
def calculate_sasa_with_freesasa(pdb_file):
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    total_sasa = result.totalArea()  # 总的溶剂可及表面积
    return total_sasa

# 结合使用 FreeSASA 计算的SASA和疏水性进行溶解性预测
def calculate_hydrophobicity(amino_acids):
    hydrophobicity = {
        'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
        'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
        'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
        'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
    }
    total_hydrophobicity = sum(hydrophobicity.get(aa, 0) for aa in amino_acids)
    return total_hydrophobicity / len(amino_acids) if amino_acids else 0

def predict_solubility(hydrophobicity, sasa):
    #结合sasa和疏水性预测溶解性
    solubility_score = sasa - (hydrophobicity * 100)
    return solubility_score

def process_and_save(folder_path, batch_size=10, output_file='./solubility_results.json'):
    processed_count = 0
    batch_number = 1
    results = load_existing_results(output_file)
    skipped_folders = []

    for dataSetFolder in os.listdir(folder_path):
        path_to_dataSet_folder = os.path.join(folder_path, dataSetFolder)
        for subFolder in os.listdir(path_to_dataSet_folder):
            if subFolder in ["ligands133", "receptors133", "results"]:
                skipped_folders.append(f"Skipped {subFolder}: Excluded folder")
                continue
            if subFolder in results:
                skipped_folders.append(f"Skipped {subFolder}: Already processed")
                continue
            
            path_to_single_folder = os.path.join(path_to_dataSet_folder, subFolder)
            if not os.path.isdir(path_to_single_folder):
                skipped_folders.append(f"Skipped {subFolder}: Not a directory")
                continue
            
            pdb_file = os.path.join(path_to_single_folder, "peptide.pdb")
            if not os.path.exists(pdb_file):
                skipped_folders.append(f"Skipped {subFolder}: No peptide.pdb file")
                continue

            try:
                # 处理数据
                amino_acids = get_amino_acids_from_pdb(pdb_file)
                sasa_value = calculate_sasa_with_freesasa(pdb_file)
                avg_hydrophobicity = calculate_hydrophobicity(amino_acids)
                solubility_prediction = predict_solubility(avg_hydrophobicity, sasa_value)

                results[subFolder] = {
                    "sasa": sasa_value,
                    "hydrophobicity": avg_hydrophobicity,
                    "solubility_score": solubility_prediction
                }

                processed_count += 1
                print(f"Processed {subFolder}")

                # 当处理的项目数达到batch_size时，进行持久化
                if processed_count % batch_size == 0:
                    save_to_file(results, output_file)
                    print(f"Completed and saved batch {batch_number}")
                    batch_number += 1
                    # 重新加载结果，以获取可能的外部更改
                    results = load_existing_results(output_file)
            except Exception as e:
                error_message = f"Error processing {subFolder}: {str(e)}\n{traceback.format_exc()}"
                print(error_message)
                skipped_folders.append(error_message)

    # 保存最后一批数据（如果有的话）
    if processed_count % batch_size != 0:
        save_to_file(results, output_file)
        print(f"Saved final batch")

    print("All data processed and saved.")
    print("\nSkipped folders:")
    for folder in skipped_folders:
        print(folder)

def load_existing_results(filename):
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            return json.load(f)
    return {}

def save_to_file(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"Saved data to {filename}")

# 主程序
if __name__ == "__main__":
    #数据集目录
    dataset_path = '../../dataset'
    process_and_save(dataset_path, batch_size=10)

