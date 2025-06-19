import os
import shutil

def delete_folders_without_atom(base_path):
    # 遍历base_path下的所有子文件夹
    for folder in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder)
        if os.path.isdir(folder_path):
            pdb_file_path = os.path.join(folder_path, 'peptide.pdb')
            try:
                # 打开并读取peptide.pdb文件
                with open(pdb_file_path, 'r') as file:
                    contents = file.readlines()
                # 检查文件中是否有包含"ATOM"的行
                has_atom = any("ATOM" in line for line in contents)
                # 如果没有"ATOM"，删除整个文件夹
                if not has_atom:
                    shutil.rmtree(folder_path)
                    print(f"Deleted folder: {folder_path}")
            except FileNotFoundError:
                print(f"No peptide.pdb file found in {folder_path}")
            except Exception as e:
                print(f"Error processing {folder_path}: {str(e)}")

# 调用函数，替换下面的路径为你的实际路径
delete_folders_without_atom('/root/autodl-tmp/PP_generate_v1/data/ppbench2024')
