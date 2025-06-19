# import os
# '''修改ppflow文件名'''
# def rename_folders(parent_directory):
#     # 获取指定目录下的所有条目
#     for item in os.listdir(parent_directory):
#         # 构造完整的路径
#         full_path = os.path.join(parent_directory, item)
#         # 检查这个条目是否是文件夹
#         if os.path.isdir(full_path):
#             # 分割文件夹名称并找到需要的部分
#             parts = item.split('_')
#             # 假设 '3rul' 总是在第二部分
#             new_name = parts[1]  # 修改这个索引，如果 '3rul' 的位置不同
#             new_full_path = os.path.join(parent_directory, new_name)
#             # 重命名文件夹
#             os.rename(full_path, new_full_path)
#             print(f"Renamed '{full_path}' to '{new_full_path}'")

# # 调用函数
# parent_directory = '/root/autodl-tmp/ppflow4'  # 更改为你的文件夹路径
# rename_folders(parent_directory)



# import os
# import json
# import glob
# import tqdm
# def rename_files(directory):
#     # 调用函数
#     json_file = '/root/autodl-tmp/PP_generate_v1/data/output/new_Experimental_vina.json'

#     with open(json_file, 'r') as f:
#         data = json.load(f)
#     # 遍历指定目录下的所有文件
#     for protein_id in data.keys():
#         ligand_pdb_file = os.path.join(directory, f"{protein_id}_unrelaxed_rank_001*.pdb")
#         ligand_files = glob.glob(ligand_pdb_file)
#         ligand_pdb_file = ligand_files[0]  # 选择第一个匹配的文件
#         # 构建完整的文件路径
#         full_path = os.path.join(directory, ligand_pdb_file)
#         # 确保它是一个文件而不是目录
#         if os.path.isfile(full_path):
#             # 分割文件名以获取所需的部分
#             parts = ligand_pdb_file.split('/')[-1].split("_")
#             # 查找具有特定前缀的部分并组合
#             new_filename = '_'.join([parts[0], parts[3]]) + os.path.splitext(full_path)[1]  # 添加文件扩展名
#             # 构建新文件的完整路径
#             new_full_path = os.path.join(directory, new_filename)
#             # 重命名文件
#             os.rename(full_path, new_full_path)
#             print(f"Renamed '{protein_id}' to '{new_filename}'")


# directory = '/root/autodl-fs/new_generated_peptide'  # 替换为你的文件夹路径
# rename_files(directory)





# import os
# import json
# import glob
# import tqdm
# '''修改生成多肽的文件名字'''
# def rename_files(directory):
#     # 调用函数
#     json_file = '/root/autodl-tmp/PP_generate_v1/data/output/ppbench.json'

#     with open(json_file, 'r') as f:
#         data = json.load(f)
#     # 遍历指定目录下的所有文件
#     for protein_id in data.keys():
#         ligand_pdb_file = os.path.join(directory, f"{protein_id}_unrelaxed_*_001_*.pdb")
#         ligand_files = glob.glob(ligand_pdb_file)
#         ligand_pdb_file = ligand_files[0]  # 选择第一个匹配的文件
#         # 构建完整的文件路径
#         full_path = os.path.join(directory, ligand_pdb_file)
#         # 确保它是一个文件而不是目录
#         if os.path.isfile(full_path):
#             # 分割文件名以获取所需的部分
#             parts = ligand_pdb_file.split('/')[-1].split("_")
#             # 查找具有特定前缀的部分并组合
#             new_filename = parts[0] + os.path.splitext(full_path)[1]  # 添加文件扩展名
#             # 构建新文件的完整路径
#             new_full_path = os.path.join(directory, new_filename)
#             # 重命名文件
#             os.rename(full_path, new_full_path)
#             print(f"Renamed '{protein_id}' to '{new_filename}'")


# directory = '/root/autodl-tmp/peptides5'  # 替换为你的文件夹路径
# rename_files(directory)




import os
import json
import glob
from tqdm import tqdm

'''修改RFdiffusion文件名字'''

def rename_files(directory, json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    # 遍历JSON文件中的所有protein_id
    for protein_id in tqdm(data.keys(), desc="Renaming files"):
        # 搜索与当前protein_id匹配的文件
        pattern = os.path.join(directory, f"design_{protein_id}_flexible_peptide_*.pdb")
        matching_files = glob.glob(pattern)

        for i, old_file_path in enumerate(matching_files):
            # 定义新文件名（若有多个文件则分别标记）
            new_file_name = f"{protein_id}_{i+1}.pdb" if len(matching_files) > 1 else f"{protein_id}.pdb"
            new_file_path = os.path.join(directory, new_file_name)

            # 执行文件重命名
            os.rename(old_file_path, new_file_path)
            print(f"Renamed '{old_file_path}' to '{new_file_name}'")


# 替换成你的实际路径
directory = '/root/autodl-tmp/RFdiffusion_top7'
json_file = '/root/autodl-tmp/PP_generate_v1/data/output/ppbench.json'
rename_files(directory, json_file)