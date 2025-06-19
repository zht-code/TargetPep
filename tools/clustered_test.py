# import os
# from Bio import SeqIO
# from Bio.PDB import PDBParser
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# import subprocess

# def create_fasta_from_dir(data_dir, output_fasta_path):
#     parser = PDBParser(PERMISSIVE=1)  # 使用宽容模式以处理可能的格式问题
#     records = []
#     for subdir in os.listdir(data_dir):
#         subdir_path = os.path.join(data_dir, subdir)
#         files = ["filtered_peptide.pdb", "filtered_receptor.pdb"]

#         for file in files:
#             file_path = os.path.join(subdir_path, file)
#             if os.path.exists(file_path):
#                 try:
#                     structure = parser.get_structure(subdir, file_path)
#                     sequence = ""
#                     for model in structure:
#                         for chain in model:
#                             sequence += "".join([residue.resname for residue in chain if residue.id[0] == ' '])
#                     if sequence:
#                         record = SeqRecord(Seq(sequence), id=f"{subdir}", description="")
#                         records.append(record)
#                 except Exception as e:
#                     print(f"Failed to process {file_path}: {e}")

#     SeqIO.write(records, output_fasta_path, "fasta")
#     print(f"Written all sequences to {output_fasta_path}")

# def run_mmseqs_clustering(fasta_path, output_dir):
#     mmseqs_db = os.path.join(output_dir, "mmseqs_db")
#     cluster_res = os.path.join(output_dir, "clusterRes")
#     tmp_dir = os.path.join(output_dir, "tmp")
    
#     subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "createdb", fasta_path, mmseqs_db], check=True)
#     subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "cluster", mmseqs_db, cluster_res, tmp_dir, "--min-seq-id", "0.8", "-c", "0.9", "--cov-mode", "1"], check=True)
#     subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "result2flat", mmseqs_db, mmseqs_db, cluster_res, os.path.join(output_dir, "cluster_result_all_seqs.fasta")], check=True)
#     subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "createtsv", mmseqs_db, mmseqs_db, cluster_res, os.path.join(output_dir, "cluster_result_cluster.tsv")], check=True)
#     subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "result2repseq", mmseqs_db, cluster_res, os.path.join(output_dir, "cluster_result_rep_seq.fasta")], check=True)
#     # subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "removedb", mmseqs_db], check=True)
#     # subprocess.run(["/root/autodl-fs/mmseqs/bin/mmseqs", "removedb", cluster_res], check=True)
#     subprocess.run(["rm", "-r", tmp_dir], check=True)

# def main():
#     data_dir = "/root/autodl-tmp/PP_generate_v1/data/test"  # 数据集目录
#     output_dir = "/root/autodl-tmp/ppflow-main/dataset/data/processed_test"  # 输出目录
#     fasta_file = os.path.join(output_dir, "combined_sequences.fasta")
    
#     os.makedirs(output_dir, exist_ok=True)
    
#     # 创建FASTA文件
#     create_fasta_from_dir(data_dir, fasta_file)
    
#     # 运行MMseqs2聚类
#     run_mmseqs_clustering(fasta_file, output_dir)
    
#     print("Clustering completed successfully.")

# if __name__ == "__main__":
#     main()









# import json

# # JSON文件路径
# json_file_path = '/root/autodl-tmp/PP_generate_v1/data/output/new_Experimental.json'

# # 输出FASTA文件的路径
# # fasta_file_path = '/root/autodl-tmp/ppflow-main/dataset/data/processed_test/cluster_result_rep_seq1.fasta'
# fasta_file_path = '/root/autodl-tmp/ppflow-main/dataset/data/processed_test/cluster_result_all_seq1.fasta'

# # 加载JSON数据
# with open(json_file_path, 'r') as file:
#     data = json.load(file)

# # 打开一个文件用于写入FASTA格式的数据
# with open(fasta_file_path, 'w') as fasta_file:
#     # 遍历每个多肽数据条目
#     for protein_id, info in data.items():
#         # 获取蛋白序列
#         protein_sequence = info['protein_sequence']
#         # 编写FASTA格式的头部
#         fasta_header = f">{protein_id}\n"
#         # 编写FASTA格式的序列部分
#         fasta_sequence = f"{protein_sequence}\n"
#         # 将头部和序列写入文件
#         fasta_file.write(fasta_header)
#         fasta_file.write(fasta_sequence)

# print("FASTA file has been written successfully.")





import torch
# 指定文件路径
file_path = '/root/autodl-fs/pp_esm3/best_model_epoch_70_loss_1.7848.pth'
# 加载文件
data = torch.load(file_path)
# 输出加载的数据
print(data)



# import torch
# import json
# # JSON文件路径
# json_file_path = '/root/autodl-tmp/PP_generate_v1/data/output/new_Experimental.json'
# # 加载JSON数据
# with open(json_file_path, 'r') as file:
#     data = json.load(file)
# test = []
# # 打开一个文件用于写入FASTA格式的数据
# # 遍历每个多肽数据条目
# for protein_id, info in data.items():
#     test.append(protein_id)
# # 定义数据字典
# # data_dict = {
# #     'train': [],
# #     'val': [],
# #     'test': test
# # }
# # 保存数据到.pt文件
# # file_path = '/root/autodl-tmp/ppflow-main/dataset/data/processed_test/split.pt'
# file_path = '/root/autodl-tmp/ppflow-main/dataset/data/processed_test/pdb_benchmark.pt'
# torch.save(test, file_path)
# print("Data saved to", file_path)


# import os
# from Bio import PDB
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# from Bio.SeqIO import write as seq_write

# aa3to1 = {
#     'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
#     'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
#     'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
#     'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
# }

# def extract_sequences(pdb_path):
#     """从PDB文件中提取序列"""
#     parser = PDB.PDBParser()
#     structure = parser.get_structure('PDB', pdb_path)
#     sequences = []
#     for model in structure:
#         for chain in model:
#             sequence = ''
#             for residue in chain:
#                 if residue.resname in aa3to1:  # 确保只转换标准氨基酸
#                     sequence += aa3to1[residue.resname]
#             if sequence:  # 如果序列非空
#                 sequences.append(sequence)
#     return sequences

# def save_sequences_to_fasta(directory, output_fasta):
#     """遍历目录中的所有PDB文件,并将蛋白和多肽序列保存到一个FASTA文件中"""
#     records = []
#     for subdir in os.listdir(directory):
#         subdir_path = os.path.join(directory, subdir)
#         pdb_files = ['filtered_peptide.pdb', 'filtered_receptor.pdb']
#         for pdb_file in pdb_files:
#             pdb_path = os.path.join(subdir_path, pdb_file)
#             if os.path.exists(pdb_path):
#                 sequences = extract_sequences(pdb_path)
#                 for i, sequence in enumerate(sequences, 1):
#                     record = SeqRecord(Seq(sequence), id=f"{subdir}", description="")
#                     records.append(record)
                    
#     seq_write(records, output_fasta, "fasta")
#     print(f"Saved all sequences to {output_fasta}")

# # Example usage
# directory = "/root/autodl-tmp/PP_generate_v1/data/test"  # Update this path to your directory
# output_fasta = "/root/autodl-tmp/ppflow-main/dataset/data/processed_test/cluster_result_all_seq.fasta"  # Update this path to your output fasta file
# save_sequences_to_fasta(directory, output_fasta)
