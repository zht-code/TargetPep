import os,sys
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# 确保正确加载模块路径
sys.path.append("/home/zht/programs/PP_generate_v1")
from models.esm.utils.structure.protein_chain import ProteinChain
def extract_sequence_from_pdb(pdb_path):
    protein_chain = ProteinChain.from_pdb(pdb_path)
    return protein_chain.sequence
def extract_sequences_from_pdb(pdb_folder, output_fasta):
    # 创建PDB解析器
    parser = PDB.PDBParser(QUIET=True)
    # 创建一个空列表来存储提取的序列
    sequence_records = []

    # 遍历存放PDB文件的文件夹
    for folder_name in os.listdir(pdb_folder):
        folder_path = os.path.join(pdb_folder, folder_name)

        if os.path.isdir(folder_path):
            pdb_file = os.path.join(folder_path, f"{folder_name}.pdb")

            if os.path.exists(pdb_file):
                # 解析PDB文件
                structure = parser.get_structure(folder_name, pdb_file)

                for model in structure:
                    for chain in model:
                        # 提取链的氨基酸序列
                        seq = []
                        for residue in chain:
                            if PDB.is_aa(residue, standard=True):
                                seq.append(PDB.Polypeptide.three_to_one(residue.resname))
                        sequence_str = ''.join(seq)
                        
                        # 创建FASTA序列记录
                        sequence_record = SeqRecord(Seq(sequence_str), id=folder_name, description=f"Chain {chain.id}")
                        sequence_records.append(sequence_record)

    # 将所有序列写入FASTA文件
    with open(output_fasta, "w") as fasta_output:
        SeqIO.write(sequence_records, fasta_output, "fasta")

    print(f"FASTA sequences saved to {output_fasta}")
def batch_extract_sequences(pdb_folder, output_fasta):
    """
    批量提取PDB文件中的多肽序列，并保存到FASTA文件中
    """
    sequence_records = []

    # 遍历PDB文件所在的文件夹
    for folder_name in os.listdir(pdb_folder):
        folder_path = os.path.join(pdb_folder, folder_name)

        if os.path.isdir(folder_path):
            pdb_file = os.path.join(folder_path, f"peptide.pdb")

            if os.path.exists(pdb_file):
                try:
                    # 使用自定义的提取函数提取序列
                    sequence_str = extract_sequence_from_pdb(pdb_file)
                    
                    # 创建FASTA格式的序列记录
                    sequence_record = SeqRecord(Seq(sequence_str), id=folder_name, description=f"Extracted from {folder_name}")
                    sequence_records.append(sequence_record)
                except Exception as e:
                    print(f"Error processing {pdb_file}: {e}")

    # 将所有序列保存为FASTA文件
    with open(output_fasta, "w") as fasta_output:
        SeqIO.write(sequence_records, fasta_output, "fasta")

    print(f"FASTA sequences saved to {output_fasta}")
# 设定PDB文件所在的主文件夹路径和输出的FASTA文件路径
pdb_folder = '/home/zht/programs/PP_generate_v1/data/ppbench2024'  # 修改为实际PDB文件的主文件夹路径
output_fasta = '/home/zht/programs/PP_generate_v1/data/peptide_sequences.fasta'


# 调用函数进行批量提取并保存序列
batch_extract_sequences(pdb_folder, output_fasta)
