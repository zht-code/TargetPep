# from Bio.PDB import MMCIFParser, PDBIO
# import sys

# def convert_cif_to_pdb(input_cif, output_pdb):
#     """
#     Convert a .cif file to .pdb format using Biopython.
    
#     :param input_cif: Path to the input CIF file.
#     :param output_pdb: Path to the output PDB file.
#     """
#     parser = MMCIFParser(QUIET=True)  # QUIET=True suppresses warnings
#     structure = parser.get_structure("protein", input_cif)
    
#     io = PDBIO()
#     io.set_structure(structure)
#     io.save(output_pdb)

# if __name__ == "__main__":
#     # Example usage
#     input_cif = "/root/autodl-tmp/PP_generate_v1/data/Ablation/quvina/1cjr/seed_101/predictions/1cjr_seed_101_sample_0.cif"  # Replace with your CIF file
#     output_pdb = "/root/autodl-tmp/PP_generate_v1/data/Ablation/quvina/1cjr/1cjr.pdb"  # Replace with desired PDB file name
#     convert_cif_to_pdb(input_cif, output_pdb)
#     print(f"Conversion completed: {output_pdb}")



import os
from Bio.PDB import MMCIFParser, PDBIO
def convert_cif_to_pdb(input_cif, output_pdb):
    """
    Convert a .cif file to .pdb format using Biopython.
    
    :param input_cif: Path to the input CIF file.
    :param output_pdb: Path to the output PDB file.
    """
    parser = MMCIFParser(QUIET=True)  # QUIET=True suppresses warnings
    structure = parser.get_structure("protein", input_cif)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
def batch_convert_cif_to_pdb(root_folder):
    """
    Batch convert all .cif files under subfolders of the root_folder to .pdb files.
    
    :param root_folder: The root directory containing subfolders with .cif files.
    """
    for subdir, dirs, files in os.walk(root_folder):
        for dir in dirs:
            input_cif = os.path.join(subdir, f'{dir}/seed_101/predictions/{dir}_seed_101_sample_0.cif')
            # 生成对应的pdb文件路径，例如：1cka/fold_xxx_model_0.cif -> 1cka/1cka_af3.pdb
            output_pdb = os.path.join(subdir, f"{dir}/{dir}.pdb")
            convert_cif_to_pdb(input_cif, output_pdb)
            print(f"Conversion completed: {output_pdb}")


if __name__ == "__main__":
    root_folder = "/root/autodl-tmp/PP_generate_v1/data/Ablation/quvina"  # 改成你的总文件夹路径
    batch_convert_cif_to_pdb(root_folder)
