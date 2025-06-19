# import os
# import shutil
# '''批量修改蛋白、多肽文件名'''
# # 路径设置为你的文件夹所在的上级目录
# base_path = '/root/autodl-tmp/PP_generate_v1/data/test'
# # 遍历base_path下的所有文件夹
# for folder_name in os.listdir(base_path):
#     folder_path = os.path.join(base_path, folder_name)
    
#     # 检查是否是文件夹
#     if os.path.isdir(folder_path):
#         # 定义旧文件和新文件的名字
#         old_protein_name = os.path.join(folder_path, f'{folder_name}_protein.pdb')
#         new_protein_name = os.path.join(folder_path, 'receptor.pdb')
#         old_pocket_name = os.path.join(folder_path, f'{folder_name}_pocket.pdb')
#         new_pocket_name = os.path.join(folder_path, 'peptide.pdb')

#         # 重命名文件
#         if os.path.exists(old_protein_name):
#             shutil.move(old_protein_name, new_protein_name)
#         if os.path.exists(old_pocket_name):
#             shutil.move(old_pocket_name, new_pocket_name)





# import os
# import shutil
# '''将test文件中与ppbench2024 id相同的删除掉'''
# # 设置两个主文件夹的路径
# ppbench2024_path = '/root/autodl-tmp/PP_generate_v1/data/ppbench2024'
# test_path = '/root/autodl-tmp/PP_generate_v1/data/test'
# # 获取ppbench2024中所有子文件夹的名称，移除后面的"_A"后缀
# ppbench2024_folders = set(folder.split('_')[0] for folder in os.listdir(ppbench2024_path))

# # 遍历test文件夹中的每个子文件夹
# for test_folder in os.listdir(test_path):
#     # 检查test_folder是否在ppbench2024_folders集合中
#     if test_folder in ppbench2024_folders:
#         folder_to_delete = os.path.join(test_path, test_folder)
        
#         # 确保路径是一个文件夹
#         if os.path.isdir(folder_to_delete):
#             # 使用shutil.rmtree来删除文件夹及其所有内容
#             shutil.rmtree(folder_to_delete)
#             print(f"Deleted folder: {folder_to_delete}")
#         else:
#             print(f"No directory found at {folder_to_delete}")
#     else:
#         print(f"Folder {test_folder} is not matched in ppbench2024 and is not deleted.")




# import os
# '''将结构文件处理成fasta'''
# def extract_sequence_from_pdb(pdb_file):
#     """从PDB文件中提取蛋白质序列（只考虑CA原子）"""
#     sequence = ''
#     with open(pdb_file, 'r') as file:
#         for line in file:
#             if line.startswith('ATOM') and line[12:16].strip() == 'CA':
#                 amino_acid = three_to_one(line[17:20].strip())
#                 sequence += amino_acid
#     return sequence

# def three_to_one(amino_acid):
#     """将三字母代码的氨基酸转换为一字母代码"""
#     conversion = {
#         'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
#         'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
#         'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
#         'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
#     }
#     return conversion.get(amino_acid, '?')  # 使用问号作为未知氨基酸的占位符

# def write_fasta(output_file, folder_path):
#     """遍历文件夹，并将所有PDB中的序列写入一个FASTA文件"""
#     with open(output_file, 'w') as fasta_file:
#         for subdir, dirs, files in os.walk(folder_path):
#             for file in files:
#                 if file.endswith('peptide.pdb'):
#                     pdb_path = os.path.join(subdir, file)
#                     sequence = extract_sequence_from_pdb(pdb_path)
#                     fasta_header = f'>{os.path.basename(subdir)}'
#                     fasta_file.write(f'{fasta_header}\n{sequence}\n')

# if __name__ == '__main__':
#     folder_path = '/root/autodl-tmp/PP_generate_v1/data/test'  # 修改为你的文件夹路径
#     output_file = '/root/autodl-tmp/PP_generate_v1/data/test_peptide.fasta'  # 设置输出FASTA文件的名字
#     write_fasta(output_file, folder_path)
#     print(f'FASTA file has been created at {output_file}')




# import os
# import subprocess
# "删除多肽序列之间相似度高于30的文件夹"
# def extract_sequence_from_pdb(pdb_file):
#     sequence = ''
#     with open(pdb_file, 'r') as file:
#         for line in file:
#             if line.startswith('ATOM') and line[12:16].strip() == 'CA':
#                 amino_acid = line[17:20].strip()
#                 sequence += three_to_one(amino_acid)
#     return sequence

# def three_to_one(amino_acid):
#     conversion = {
#         'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
#         'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
#         'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
#         'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
#     }
#     return conversion.get(amino_acid, '')

# def create_fasta_file(folder_path, output_fasta):
#     with open(output_fasta, 'w') as fasta_file:
#         for subdir in os.listdir(folder_path):
#             pdb_path = os.path.join(folder_path, subdir, 'peptide.pdb')
#             if os.path.exists(pdb_path):
#                 sequence = extract_sequence_from_pdb(pdb_path)
#                 fasta_file.write(f'>{subdir}\n{sequence}\n')

# def make_blast_db(fasta_file, db_name):
#     cmd = ['/root/autodl-tmp/PP_generate_v1/methods/ncbi-blast-2.16.0+/bin/makeblastdb', '-in', fasta_file, '-dbtype', 'prot', '-out', db_name]
#     subprocess.run(cmd, check=True)

# def run_blastp(db_name, fasta_file, output_file):
#     cmd = ['/root/autodl-tmp/PP_generate_v1/methods/ncbi-blast-2.16.0+/bin/blastp', '-db', db_name, '-query', fasta_file, '-outfmt', '6 qacc sacc pident', '-out', output_file]
#     subprocess.run(cmd, check=True)

# def parse_blast_results(output_file):
#     to_remove = set()
#     with open(output_file, 'r') as results:
#         for line in results:
#             parts = line.strip().split()
#             if float(parts[2]) > 30:  # parts[2] is the percentage identity
#                 to_remove.add(parts[0])  # parts[0] is the query accession
#                 to_remove.add(parts[1])  # parts[1] is the subject accession
#     return to_remove

# def remove_folders(folder_path, to_remove):
#     for subdir in os.listdir(folder_path):
#         if subdir in to_remove:
#             full_path = os.path.join(folder_path, subdir)
#             subprocess.run(['rm', '-rf', full_path])
#             print(f'Removed folder: {full_path}')

# def main(folder_path):
#     fasta_file = '/root/autodl-tmp/PP_generate_v1/data/test_peptide.fasta'
#     db_name = '/root/autodl-tmp/PP_generate_v1/data/peptides_db'
#     output_file = '/root/autodl-tmp/PP_generate_v1/data/blast_results.txt'

#     # create_fasta_file(folder_path, fasta_file)
#     # make_blast_db(fasta_file, db_name)
#     # run_blastp(db_name, fasta_file, output_file)
#     to_remove = parse_blast_results(output_file)
#     remove_folders(folder_path, to_remove)
#     print(f'BLASTP results processed and folders removed if needed.')

# if __name__ == '__main__':
#     folder_path = '/root/autodl-tmp/PP_generate_v1/data/test'  # 修改为你的文件夹路径
#     main(folder_path)



import os
import glob
from Bio import PDB
'''处理多肽肽链大于1的多肽文件'''
def filter_single_chain_peptide(input_pdb_path, output_pdb_path):
    """
    函数将从输入的PDB文件中只保留第一条链，并将结果保存到输出路径。
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('Peptide', input_pdb_path)
    model = structure[0]  # 默认处理第一个模型
    chain_ids = [chain.id for chain in model.get_chains()]

    if len(chain_ids) > 1:
        # 选取第一条链
        chain_to_keep = model[chain_ids[0]]
        # 创建新的结构来保存过滤后的链
        io = PDB.PDBIO()
        io.set_structure(chain_to_keep)
        io.save(output_pdb_path)
        print(f"Filtered {input_pdb_path}, retained chain {chain_ids[0]}, saved to {output_pdb_path}")
    else:
        # 如果只有一条链，跳过这个文件
        print(f"Only one chain in {input_pdb_path}, skipping.")
        return None  # 这将表示此文件不需要进一步操作

def process_folder(folder_path, RFdiffusion_path):
    """
    遍历给定文件夹内的所有子文件夹，处理其中的filtered_peptide.pdb文件。
    """
    # for subdir in os.listdir(folder_path):
    #     subdir_path = os.path.join(folder_path, subdir)
    #     peptide_path = os.path.join(subdir_path, 'filtered_peptide.pdb')
    #     output_path = os.path.join(subdir_path, 'filtered_peptide.pdb')

    #     if os.path.exists(peptide_path):
    #         result = filter_single_chain_peptide(peptide_path, output_path)
    #         if result is None:
    #             continue  # 如果函数返回None，跳过当前文件
    #     else:
    #         print(f"No peptide file found in {subdir_path}")
    # a = []
    for subdir in os.listdir(folder_path):
        peptide_path = os.path.join(RFdiffusion_path, f'design_{subdir}_flexible_peptide_0.pdb')
        ligand_files = glob.glob(peptide_path)
        # if len(ligand_files)==0 :
        #     a.append(subdir)
        #     continue
        ligand_pdb_file = ligand_files[0]
        output_path = os.path.join(RFdiffusion_path, f'new_design_{subdir}_flexible_peptide_0.pdb')

        if os.path.exists(peptide_path):
            result = filter_single_chain_peptide(ligand_pdb_file, output_path)
            if result is None:
                continue  # 如果函数返回None，跳过当前文件
        # else:
        #     print(f"No peptide file found in {subdir_path}")

if __name__ == "__main__":
    folder_path = '/root/autodl-tmp/PP_generate_v1/data/test'  # 修改为你的文件夹路径
    RFdiffusion_path = '/root/autodl-tmp/PP_generate_v1/data/downloads/RFdiffusion_outputs1'
    process_folder(folder_path, RFdiffusion_path)




# import os
# import tqdm
# import glob
# from Bio import PDB
# '''处理多肽肽链大于1的多肽文件'''
# def filter_single_chain_peptide(input_pdb_path, output_pdb_path):
#     """
#     函数将从输入的PDB文件中只保留第一条链，并将结果保存到输出路径。
#     """
#     parser = PDB.PDBParser()
#     structure = parser.get_structure('Peptide', input_pdb_path)
#     model = structure[0]  # 默认处理第一个模型
#     chain_ids = [chain.id for chain in model.get_chains()]

#     if len(chain_ids) > 1:
#         # 选取第一条链
#         chain_to_keep = model[chain_ids[0]]
#         # 创建新的结构来保存过滤后的链
#         io = PDB.PDBIO()
#         io.set_structure(chain_to_keep)
#         io.save(output_pdb_path)
#         print(f"Filtered {input_pdb_path}, retained chain {chain_ids[0]}, saved to {output_pdb_path}")
#     else:
#         # 如果只有一条链，跳过这个文件
#         print(f"Only one chain in {input_pdb_path}, skipping.")
#         return None  # 这将表示此文件不需要进一步操作

# def process_folder(folder_path):
#     """
#     遍历给定文件夹内的所有子文件夹，处理其中的filtered_peptide.pdb文件。
#     """
#     for i in range(10):
#         peptide_path = os.path.join(folder_path,  f"design_1t1f{i+1}_*.pdb")
#         ligand_files = glob.glob(peptide_path)
#         ligand_pdb_file = ligand_files[0]
#         output_path = os.path.join(folder_path, f'new_design_1t1f{i+1}_peptide.pdb')


#         if os.path.exists(ligand_pdb_file):
#             result = filter_single_chain_peptide(ligand_pdb_file, output_path)
#             if result is None:
#                 continue  # 如果函数返回None，跳过当前文件


# if __name__ == "__main__":
#     folder_path = '/root/autodl-tmp/PP_generate_v1/data/wet_test/1t1f_gene'  # 修改为你的文件夹路径
#     process_folder(folder_path)





# import os
# import numpy as np
# from Bio.PDB import PDBParser, PDBIO, Select,NeighborSearch
# from Bio.PDB.Polypeptide import is_aa
# '''将蛋白和多肽pdb文件按照ppflow的方法进行处理'''
# def read_pdb(file_path):
#     parser = PDBParser()
#     structure = parser.get_structure('PDB', file_path)
#     return structure

# def filter_peptide_length(peptide_structure, max_length=30):
#     residues = [res for res in peptide_structure.get_residues() if is_aa(res, standard=True)]
#     return len(residues) <= max_length

# def filter_water_and_heteroatoms(structure):
#     # Create a new empty structure to store filtered model
#     filtered_structure = structure.copy()
#     for model in filtered_structure:
#         for chain in list(model):
#             for res in list(chain):
#                 # Remove non-standard amino acids and water
#                 if res.id[0] != ' ' or not is_aa(res, standard=True):
#                     chain.detach_child(res.id)
#     return filtered_structure

# def write_filtered_pdb(structure, output_file):
#     io = PDBIO()
#     io.set_structure(structure)
#     io.save(output_file)

# def minimum_distance(receptor_structure, peptide_structure):
#     ns = NeighborSearch(list(receptor_structure.get_atoms()))
#     min_distance = float('inf')
#     for atom in peptide_structure.get_atoms():
#         close_atoms = ns.search(atom.coord, 5.0, 'A')
#         for close_atom in close_atoms:
#             distance = np.linalg.norm(atom.coord - close_atom.coord)
#             if distance < min_distance:
#                 min_distance = distance
#     return min_distance

# def process_folder(folder_path):
#     valid_peptides = []
#     for subdir in os.listdir(folder_path):
#         subdir_path = os.path.join(folder_path, subdir)
#         receptor_path = os.path.join(subdir_path, 'receptor.pdb')
#         peptide_path = os.path.join(subdir_path, 'peptide.pdb')

#         if os.path.exists(receptor_path) and os.path.exists(peptide_path):
#             receptor_structure = read_pdb(receptor_path)
#             peptide_structure = read_pdb(peptide_path)

#             filtered_receptor = filter_water_and_heteroatoms(receptor_structure)
#             filtered_peptide = filter_water_and_heteroatoms(peptide_structure)

#             receptor_output_path = os.path.join(subdir_path, 'filtered_receptor.pdb')
#             peptide_output_path = os.path.join(subdir_path, 'filtered_peptide.pdb')
#             write_filtered_pdb(filtered_receptor, receptor_output_path)
#             write_filtered_pdb(filtered_peptide, peptide_output_path)

#             if filter_peptide_length(filtered_peptide) and minimum_distance(filtered_receptor, filtered_peptide) <= 5.0:
#                 valid_peptides.append(subdir)

#     return valid_peptides

# def main():
#     folder_path = '/root/autodl-tmp/PP_generate_v1/data/test'
#     valid_peptides = process_folder(folder_path)
#     print("Valid peptides folders:", valid_peptides)

# if __name__ == '__main__':
#     main()
