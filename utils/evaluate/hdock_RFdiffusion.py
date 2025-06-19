# import subprocess
# import os
# from Bio import PDB
# import numpy as np
# import json
# from tqdm import tqdm  # 用于显示进度条，可选
# import glob
# '''
# 批量计算蛋白和多肽的亲和力,需使用py文件运行
# '''
# # 自定义选择器 用于过滤
# class MySelect(PDB.Select):
#     '''预处理蛋白质---选择类'''
#     def accept_residue(self, residue):
#         '''过滤 保留主要的氨基酸残基'''
#         return residue.get_resname() not in ['HOH', 'Na', 'Cl', 'K', 'Mg', 'Ca']

# # 过滤函数    
# def filter_structure_without_h2oNion(protein_structure, protein_name):
#     '''预处理蛋白质---过滤水和离子'''
#     filtered_protein_structure = PDB.Structure.Structure("filtered_" + protein_name)
#     for model in protein_structure:
#         new_model = PDB.Model.Model(model.id)
#         for chain in model:
#             new_chain = PDB.Chain.Chain(chain.id)
#             for residue in chain:
#                 if MySelect().accept_residue(residue):
#                     new_chain.add(residue.copy())  # 添加过滤后的残基
#             if len(new_chain) > 0:
#                 new_model.add(new_chain)
#         filtered_protein_structure.add(new_model)
#     return filtered_protein_structure

# # 转换PDB为PDBQT格式
# def pdb_to_pdbqt(pdb_file, pdbqt_file, is_receptor=True):
#     """将 PDB 文件转换为 PDBQT 文件"""
#     script = 'prepare_receptor4.py' if is_receptor else 'prepare_ligand4.py'
#     command = [
#         '/root/mgltools_x86_64Linux2_1.5.7/bin/pythonsh',
#         f'/root/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/{script}',
#         '-r' if is_receptor else '-l', pdb_file,
#         '-o', pdbqt_file
#     ]
#     print(f"Running command: {' '.join(command)}")

#     # 运行命令并捕获输出和错误信息
#     # subprocess.run(command, capture_output=True, text=True)
#     # 运行命令并捕获输出和错误信息
#     env = os.environ.copy()
#     env.pop('PYTHONHOME', None)
#     env.pop('PYTHONPATH', None)
#     subprocess.run(command, capture_output=True, text=True, env=env)
#     # subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#     # result = subprocess.run(command, shell=True, check=True)
#     # print(f"Command error: {result.stderr}")
    
#     # print(f'Wrote converted file to {pdbqt_file}')
#     return pdbqt_file

# # 运行 QVina 并解析日志文件，获取最好的 affinity
# def run_qvina_and_parse(path_to_single_folder, receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z):
#     """运行 QVina 并解析日志文件，获取最好的 affinity"""
#     output_pdbqt = os.path.join(path_to_single_folder, "RFdiffusion_vina.pdbqt")
    
#     # 构建命令
#     command = [
#         '/root/qvina/qvina-master/bin/qvina02',
#         '--receptor', receptor_pdbqt,
#         '--ligand', ligand_pdbqt,
#         '--center_x', str(center_x),
#         '--center_y', str(center_y),
#         '--center_z', str(center_z),
#         '--size_x', '20',  # 可以根据需要调整
#         '--size_y', '20',  # 可以根据需要调整
#         '--size_z', '20',  # 可以根据需要调整
#         '--out', output_pdbqt,
#         '--exhaustiveness', '8'
#     ]
    
#     # print(f"Running command: {' '.join(command)}")
    

#     # 打印详细的输出信息
#     # if result.returncode != 0:
#     #     print(f"Error running QVina: {result.stderr}")
#     #     raise RuntimeError(f"QVina failed for receptor {receptor_pdbqt} and ligand {ligand_pdbqt}")
    
#     print("QVina finished. Parsing the results...")

#     subprocess.run(command, capture_output=True, text=True)
#     with open(output_pdbqt, 'r') as f:
#         lines = f.readlines()
    
#     # 亲和力值通常在Vina的输出文件中以以下格式出现
#     for line in lines:
#         if line.startswith("REMARK VINA RESULT"):
#             # 亲和力值在"REMARK VINA RESULT"行之后，通常是以负值表示亲和力
#             affinity = float(line.split()[3])  # 亲和力值位于该行的第四个字段
#             return affinity

# # 保存过滤后的结构到 PDB 文件
# def save_filtered_structure(filtered_structure, output_filename):
#     io = PDB.PDBIO()
#     io.set_structure(filtered_structure)
#     io.save(output_filename)

# def main():
#     # import argparse

#     # parser = argparse.ArgumentParser(description="批量计算蛋白和多肽的亲和力并更新 JSON 文件")
#     # parser.add_argument('--json_file', type=str, default='/root/autodl-tmp/PP_generate_v1/data/result/peptides_physicochemical.json', help='已有的 JSON 文件路径')
#     # parser.add_argument('--output_json', type=str, default='/root/autodl-tmp/PP_generate_v1/data/result/peptides_vina.json', help='输出的更新后的 JSON 文件路径')
#     # parser.add_argument('--receptor_dir', type=str, default='/root/autodl-tmp/PP_generate_v1/data/PPDbench', help='蛋白 PDB 文件所在的根目录')
#     # parser.add_argument('--ligand_dir', type=str, default='/root/autodl-tmp/PP_generate_v1/data/downloads/result', help='多肽 PDB 文件所在的根目录')
#     # parser.add_argument('--center', type=float, nargs=3, default=[0, 0, 0], help='Qvina 的中心坐标 (x y z)')
    
#     # args = parser.parse_args()

#     # json_file = '/root/autodl-tmp/PP_generate_v1/data/result/peptides_vina3.json'
#     # # json_file = '/root/autodl-tmp/PP_generate_v1/data/result/peptides_vina.json'
#     # output_json = '/root/autodl-tmp/PP_generate_v1/data/result/peptides_vina4.json'
#     # receptor_dir = '/root/autodl-tmp/PP_generate_v1/data/PPDbench'
#     # ligand_dir = '/root/autodl-tmp/PP_generate_v1/data/downloads/pep_outputs'


#     json_file = '/root/autodl-tmp/PP_generate_v1/data/output/Experimental_vina3.json'
#     output_json = '/root/autodl-tmp/PP_generate_v1/data/output/Experimental_vina4.json'
#     receptor_dir = '/root/autodl-tmp/PP_generate_v1/data/Peptide_test'
#     ligand_dir =  '/root/autodl-tmp/PP_generate_v1/data/downloads/pep_outputs'
#     # 新的输出目录，用于存放生成的 PDBQT 和 Vina 输出文件
#     output_dir = '/root/autodl-tmp/PP_generate_v1/data/geneta_pdbqt'
#     center_x, center_y, center_z = [0, 0, 0]

#     # 读取已有的 JSON 文件
#     with open(json_file, 'r') as f:
#         data = json.load(f)
    
#     # 用于加载结构的PDBParser
#     pdb_parser = PDB.PDBParser(QUIET=True)
    
#     # 遍历每个蛋白 ID
#     for protein_id in tqdm(data.keys(), desc="Processing proteins"):
#         try:
#             print(f"\nProcessing protein ID: {protein_id}")
            
#             # 构建蛋白和多肽的 PDB 文件路径
#             receptor_pdb_file = os.path.join(receptor_dir, f"{protein_id}/receptor.pdb")
#             # 假设多肽文件名遵循一定的模式，如 "{protein_id}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
#             ligand_pdb_file = os.path.join(ligand_dir, f"*_{protein_id}_*_0.pdb")
#             ligand_files = glob.glob(ligand_pdb_file)
#             if not os.path.exists(receptor_pdb_file):
#                 print(f"Receptor PDB file does not exist: {receptor_pdb_file}. Skipping.")
#                 continue
#             if not ligand_files:
#                 print(f"No ligand PDB files match the pattern {ligand_pdb_file}. Skipping.")
#                 continue
#             elif len(ligand_files) > 1:
#                 print(f"Multiple ligand PDB files found for {protein_id}. Using the first one: {ligand_files[0]}")
#             ligand_pdb_file = ligand_files[0]  # 选择第一个匹配的文件
#             # path_to_single_folder = os.path.dirname(receptor_pdb_file)
# #             # 创建目标输出文件夹（以 protein_id 命名）
#             target_folder = os.path.join(output_dir, protein_id)
#             os.makedirs(target_folder, exist_ok=True)
#             # 处理和过滤结构
#             receptor_structure = pdb_parser.get_structure("receptor", receptor_pdb_file)
#             ligand_structure = pdb_parser.get_structure("peptide", ligand_pdb_file)
            
#             filtered_receptor_structure = filter_structure_without_h2oNion(receptor_structure, "receptor")
#             filtered_ligand_structure = filter_structure_without_h2oNion(ligand_structure, "peptide")
            
#             # 保存过滤后的结构
#             filtered_receptor_file = os.path.join(target_folder, "RFdiffusion_receptor_filtered.pdb")
#             save_filtered_structure(filtered_receptor_structure, filtered_receptor_file)
#             filtered_peptide_file = os.path.join(target_folder, "RFdiffusion_peptide_filtered.pdb")
#             save_filtered_structure(filtered_ligand_structure, filtered_peptide_file)
            
#             # 将PDB转换为PDBQT
#             receptor_pdbqt = pdb_to_pdbqt(filtered_receptor_file, os.path.join(target_folder, "RFdiffusion_receptor.pdbqt"), is_receptor=True)
#             ligand_pdbqt = pdb_to_pdbqt(filtered_peptide_file, os.path.join(target_folder, "RFdiffusion_peptide.pdbqt"), is_receptor=False)
            
#             # 运行Qvina并获取亲和力 
#             affinity = run_qvina_and_parse(target_folder, receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z)
#             print(f"Best affinity for {protein_id}: {affinity} kcal/mol")
            
#             # 将亲和力添加到 JSON 数据中
#             if "properties" in data[protein_id]:
#                 data[protein_id]["properties"][0]["RFdiffusion Affinity (kcal/mol)"] = affinity
        
#         except Exception as e:
#             print(f"Error processing protein ID {protein_id}: {e}")
#             # 可以选择继续处理下一个蛋白 ID，或者终止程序
#             continue
    
#     # 保存更新后的 JSON 文件
#     with open(output_json, 'w') as f:
#         json.dump(data, f, indent=4)
#     print(f"\nUpdated JSON saved to {output_json}")

# if __name__ == "__main__":
#     main()









# import subprocess
# import os
# from Bio import PDB
# import numpy as np
# import json
# from tqdm import tqdm  # 用于显示进度条，可选
# import glob
# from Bio.PDB.NeighborSearch import NeighborSearch

# '''
# 批量计算蛋白和多肽的亲和力, 需使用 py 文件运行
# '''

# # 自定义选择器 用于过滤
# class MySelect(PDB.Select):
#     '''预处理蛋白质---选择类'''
#     def accept_residue(self, residue):
#         '''过滤 保留主要的氨基酸残基'''
#         return residue.get_resname() not in ['HOH', 'Na', 'Cl', 'K', 'Mg', 'Ca']

# # 过滤函数    
# def filter_structure_without_h2oNion(protein_structure, protein_name):
#     '''预处理蛋白质---过滤水和离子'''
#     filtered_protein_structure = PDB.Structure.Structure("filtered_" + protein_name)
#     for model in protein_structure:
#         new_model = PDB.Model.Model(model.id)
#         for chain in model:
#             new_chain = PDB.Chain.Chain(chain.id)
#             for residue in chain:
#                 if MySelect().accept_residue(residue):
#                     new_chain.add(residue.copy())  # 添加过滤后的残基
#             if len(new_chain) > 0:
#                 new_model.add(new_chain)
#         filtered_protein_structure.add(new_model)
#     return filtered_protein_structure

# # 转换PDB为PDBQT格式
# def pdb_to_pdbqt(pdb_file, pdbqt_file, is_receptor=True):
#     """将 PDB 文件转换为 PDBQT 文件"""
#     script = 'prepare_receptor4.py' if is_receptor else 'prepare_ligand4.py'
#     command = [
#         '/root/mgltools_x86_64Linux2_1.5.7/bin/pythonsh',
#         f'/root/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/{script}',
#         '-r' if is_receptor else '-l', pdb_file,
#         '-o', pdbqt_file
#     ]
#     print(f"Running command: {' '.join(command)}")
#     env = os.environ.copy()
#     env.pop('PYTHONHOME', None)
#     env.pop('PYTHONPATH', None)
#     subprocess.run(command, capture_output=True, text=True, env=env)
#     return pdbqt_file

# # 运行 QVina 并解析日志文件，获取最好的 affinity
# def run_qvina_and_parse(path_to_single_folder, receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z):
#     """运行 QVina 并解析日志文件，获取最好的 affinity"""
#     output_pdbqt = os.path.join(path_to_single_folder, "RFdiffusion_vina.pdbqt")
#     command = [
#         '/root/qvina/qvina-master/bin/qvina02',
#         '--receptor', receptor_pdbqt,
#         '--ligand', ligand_pdbqt,
#         '--center_x', str(center_x),
#         '--center_y', str(center_y),
#         '--center_z', str(center_z),
#         '--size_x', '20',
#         '--size_y', '20',
#         '--size_z', '20',
#         '--out', output_pdbqt,
#         '--exhaustiveness', '8'
#     ]
#     print("Running QVina...")
#     subprocess.run(command, capture_output=True, text=True)
#     print("QVina finished. Parsing the results...")
#     with open(output_pdbqt, 'r') as f:
#         lines = f.readlines()
#     for line in lines:
#         if line.startswith("REMARK VINA RESULT"):
#             affinity = float(line.split()[3])
#             return affinity

# # 保存过滤后的结构到 PDB 文件
# def save_filtered_structure(filtered_structure, output_filename):
#     io = PDB.PDBIO()
#     io.set_structure(filtered_structure)
#     io.save(output_filename)

# # 自动检测结合位点，计算几何中心作为对接盒中心
# def get_binding_site_center(filtered_receptor_structure, filtered_peptide_structure, threshold=5.0):
#     """
#     利用 NeighborSearch 找到蛋白和多肽之间距离小于 threshold (Å) 的原子，
#     并计算这些原子的几何中心作为对接盒中心坐标。
#     """
#     # 使用多肽中的所有原子建立搜索树
#     peptide_atoms = list(filtered_peptide_structure.get_atoms())
#     ns = NeighborSearch(peptide_atoms)
#     coords = []
#     # 遍历蛋白中的每个原子，找出与多肽原子距离小于 threshold 的邻居
#     for atom in filtered_receptor_structure.get_atoms():
#         neighbors = ns.search(atom.get_coord(), threshold)
#         if neighbors:
#             coords.append(atom.get_coord())
#             for neighbor in neighbors:
#                 coords.append(neighbor.get_coord())
#     if coords:
#         center = np.mean(np.array(coords), axis=0)
#         return tuple(center)
#     else:
#         print("Warning: 未检测到结合位点，使用默认中心 (0,0,0)")
#         return (0.0, 0.0, 0.0)

# def main():
#     json_file = '/root/autodl-tmp/PP_generate_v1/data/output/new_Experimental_vina_proteingenerator.json'
#     output_json = '/root/autodl-tmp/PP_generate_v1/data/output/new_Experimental_vina_RFdiffusion.json'
#     receptor_dir = '/root/autodl-tmp/PP_generate_v1/data/test'
#     ligand_dir =  '/root/autodl-tmp/PP_generate_v1/data/downloads/RFdiffusion_outputs1'
#     # 新的输出目录，用于存放生成的 PDBQT 和 Vina 输出文件
#     output_dir = '/root/autodl-tmp/PP_generate_v1/data/new_geneta_pdbqt'
    
#     with open(json_file, 'r') as f:
#         data = json.load(f)
    
#     pdb_parser = PDB.PDBParser(QUIET=True)
    
#     for protein_id in tqdm(data.keys(), desc="Processing proteins"):
#         try:
#             print(f"\nProcessing protein ID: {protein_id}")
            
#             # 构建蛋白和多肽的 PDB 文件路径
#             receptor_pdb_file = os.path.join(receptor_dir, f"{protein_id}/filtered_receptor.pdb")
#             ligand_pdb_file = os.path.join(ligand_dir, f"new_*_{protein_id}_*_0.pdb")
#             ligand_files = glob.glob(ligand_pdb_file)
#             if not os.path.exists(receptor_pdb_file):
#                 print(f"Receptor PDB file does not exist: {receptor_pdb_file}. Skipping.")
#                 continue
#             if not ligand_files:
#                 print(f"No ligand PDB files match the pattern {ligand_pdb_file}. Skipping.")
#                 continue
#             elif len(ligand_files) > 1:
#                 print(f"Multiple ligand PDB files found for {protein_id}. Using the first one: {ligand_files[0]}")
#             ligand_pdb_file = ligand_files[0]
#             # 创建目标输出文件夹（以 protein_id 命名）
#             target_folder = os.path.join(output_dir, protein_id)
#             os.makedirs(target_folder, exist_ok=True)
            
#             # 加载蛋白和多肽结构
#             receptor_structure = pdb_parser.get_structure("receptor", receptor_pdb_file)
#             ligand_structure = pdb_parser.get_structure("peptide", ligand_pdb_file)
            
#             # 过滤结构
#             filtered_receptor_structure = filter_structure_without_h2oNion(receptor_structure, "receptor")
#             filtered_ligand_structure = filter_structure_without_h2oNion(ligand_structure, "peptide")
            
#             # 自动检测结合位点，计算对接盒中心坐标
#             center_x, center_y, center_z = get_binding_site_center(filtered_receptor_structure, filtered_ligand_structure, threshold=5.0)
#             print(f"自动检测的对接中心: x={center_x:.2f}, y={center_y:.2f}, z={center_z:.2f}")
            
#             # 保存过滤后的结构
#             filtered_receptor_file = os.path.join(target_folder, "RFdiffusion_receptor_filtered.pdb")
#             save_filtered_structure(filtered_receptor_structure, filtered_receptor_file)
#             filtered_peptide_file = os.path.join(target_folder, "RFdiffusion_peptide_filtered.pdb")
#             save_filtered_structure(filtered_ligand_structure, filtered_peptide_file)
            
#             # 将PDB转换为PDBQT文件
#             receptor_pdbqt = pdb_to_pdbqt(filtered_receptor_file, os.path.join(target_folder, "RFdiffusion_receptor.pdbqt"), is_receptor=True)
#             ligand_pdbqt = pdb_to_pdbqt(filtered_peptide_file, os.path.join(target_folder, "RFdiffusion_peptide.pdbqt"), is_receptor=False)
            
#             # 使用自动计算的对接中心运行 QVina 并获取亲和力 
#             affinity = run_qvina_and_parse(target_folder, receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z)
#             print(f"Best affinity for {protein_id}: {affinity} kcal/mol")
            
#             # 将亲和力写入 JSON 数据中
#             if "properties" in data[protein_id]:
#                 data[protein_id]["properties"][0]["RFdiffusion Affinity (kcal/mol)"] = affinity
        
#         except Exception as e:
#             print(f"Error processing protein ID {protein_id}: {e}")
#             continue
    
#     with open(output_json, 'w') as f:
#         json.dump(data, f, indent=4)
#     print(f"\nUpdated JSON saved to {output_json}")

# if __name__ == "__main__":
#     main()




import subprocess
import os
import shutil
from Bio import PDB
import numpy as np
import json
from tqdm import tqdm  # 用于显示进度条
import glob
from Bio.PDB.NeighborSearch import NeighborSearch

class MySelect(PDB.Select):
    '''过滤保留主要的氨基酸残基'''
    def accept_residue(self, residue):
        return residue.get_resname() not in ['HOH', 'Na', 'Cl', 'K', 'Mg', 'Ca']

def filter_structure_without_h2oNion(protein_structure, protein_name):
    '''过滤水和离子'''
    filtered_protein_structure = PDB.Structure.Structure("filtered_" + protein_name)
    for model in protein_structure:
        new_model = PDB.Model.Model(model.id)
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            for residue in chain:
                if MySelect().accept_residue(residue):
                    new_chain.add(residue.copy())
            if len(new_chain) > 0:
                new_model.add(new_chain)
        filtered_protein_structure.add(new_model)
    return filtered_protein_structure

def save_filtered_structure(filtered_structure, output_filename):
    io = PDB.PDBIO()
    io.set_structure(filtered_structure)
    io.save(output_filename)

def run_hdock_and_parse(path_to_single_folder, receptor_file, ligand_file):
    """运行 HDOCK 并解析日志文件，获取最好的模型"""
    output_dir = os.path.join(path_to_single_folder, "RFdiffusion5")
    os.makedirs(output_dir, exist_ok=True)
    command = [
        '/root/autodl-fs/HDOCKlite/hdock',
        ligand_file,
        receptor_file
        # '-out',f'{output_dir}/hdock.out'
    ]
    command1 = [
        '/root/autodl-fs/HDOCKlite/createpl',
        '/root/autodl-tmp/hdock/Hdock.out',
        # '/root/Hdock.out',
        'top10.pdb',
        '-nmax', '10',
        '-complex', '-models'
    ]
    print("Running HDOCK...")
    subprocess.run(command, capture_output=True, text=True)
    subprocess.run(command1, capture_output=True, text=True)
    # result_file = os.path.join(output_dir, "model_1.pdb")  # 假设输出结果是最佳模型
     # 例如移动生成的文件到新位置
    shutil.move('/root/autodl-tmp/hdock/Hdock.out', output_dir)  # 修改路径为实际需要的路径
    shutil.move('/root/autodl-tmp/hdock/model_1.pdb', output_dir)  # 修改路径为实际需要的路径
    # with open(result_file, 'r') as f:
    with open(f'{output_dir}/model_1.pdb', 'r') as f:
    # with open(f'/root/model_1.pdb', 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("REMARK Score"):
            score = float(line.split()[2])  # 提取评分
            return score

def main():
    json_file = '/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json'
    # json_file = '/root/autodl-tmp/PP_generate_v1/data/output/new_Experimental.json'
    # output_json = '/root/autodl-tmp/PP_generate_v1/data/output/Hdock_RFdiffusion_vina1.json'
    output_json = '/root/autodl-tmp/PP_generate_v1/data/SOTA/Rdiffusion5_hdock.json'
    receptor_dir = '/root/autodl-fs/PPDbench'
    ligand_dir = '/root/autodl-tmp/RFdiffusion_top5'
    output_dir = '/root/autodl-tmp/PP_generate_v1/data/geneta_pdbqt'
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    pdb_parser = PDB.PDBParser(QUIET=True)
    
    for protein_id in tqdm(data.keys(), desc="Processing proteins"):
        try:
            receptor_pdb_file = os.path.join(receptor_dir, f"{protein_id}/receptor.pdb")
            ligand_pdb_file = os.path.join(ligand_dir, f"{protein_id}_5.pdb")
            ligand_files = glob.glob(ligand_pdb_file)
            if not os.path.exists(receptor_pdb_file):
                print(f"Receptor PDB file does not exist: {receptor_pdb_file}. Skipping.")
                continue
            if not ligand_files:
                print(f"No ligand PDB files match the pattern {ligand_pdb_file}. Skipping.")
                continue
            elif len(ligand_files) > 1:
                print(f"Multiple ligand PDB files found for {protein_id}. Using the first one: {ligand_files[0]}")
            ligand_pdb_file = ligand_files[0]
            target_folder = os.path.join(output_dir, protein_id)
            score = run_hdock_and_parse(target_folder, receptor_pdb_file, ligand_pdb_file)
            print(f"Best score for {protein_id}: {score}")
            
            if "properties" in data[protein_id]:
                data[protein_id]["properties"][0]["RFdiffusion5 Affinity (kcal/mol)"] = score
        
        except Exception as e:
            print(f"Error processing protein ID {protein_id}: {e}")
    
    with open(output_json, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"\nUpdated JSON saved to {output_json}")

if __name__ == "__main__":
    main()
