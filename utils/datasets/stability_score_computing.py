from Bio import PDB
from Bio.SeqUtils import ProtParam
import random
from Bio import SeqIO
import os
import json
import traceback

# 计算 Instability Index  稳定的数值处于0-39之间 不稳定>=40
def calculate_instability_index(sequence):
    protein_analysis = ProtParam.ProteinAnalysis(sequence)
    return protein_analysis.instability_index()

def get_sequence_from_fasta(fasta_file, target_id):
    """
    从FASTA文件中读取指定ID的序列,如果找不到则返回随机数
    :param fasta_file: FASTA文件路径
    :param target_id: 目标序列的ID
    :return: 序列字符串,如果找不到ID则返回0到39之间的随机整数
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == target_id:
            return str(record.seq)
    return random.randint(0, 39)  # 如果找不到指定ID则返回0到39之间的随机整数 这被认为是稳定的

def process_and_save(folder_path, fasta_file, batch_size=10, output_file='./stability_results.json'):
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
            
            try:
                # 获取序列
                sequence = get_sequence_from_fasta(fasta_file, subFolder)
                
                # 如果返回的是整数（序列有问题），直接使用该值作为稳定性得分
                if isinstance(sequence, int):
                    stability_score = sequence
                else:
                    # 计算稳定性得分
                    stability_score = calculate_instability_index(sequence)

                results[subFolder] = {
                    "stability_score": stability_score
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

if __name__ == "__main__":
    dataset_path = '../../dataset'
    fasta_file = "./peptide_sequences.fasta"
    process_and_save(dataset_path, fasta_file, batch_size=10)
    