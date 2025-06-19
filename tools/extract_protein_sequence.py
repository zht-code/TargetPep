# import sys
# sys.path.append("/root/autodl-tmp/PP_generate_v1/")
# from models.esm.utils.structure.protein_chain import ProteinChain

# def extract_sequence_from_pdb(pdb_path):
#     try:
#         protein_chain = ProteinChain.from_pdb(pdb_path)
#         return protein_chain.sequence if protein_chain else None
#     except Exception as e:
#         print(f"Error reading {pdb_path}: {e}")
#         return None


# if __name__ == "__main__":
#     protein_path = "/root/autodl-tmp/PP_generate_v1/tools/4hfz.pdb"
#     sequence = extract_sequence_from_pdb(protein_path)
#     print('protein sequence:',sequence)

# -----------------------------------------------------------------------------------------------------------

# import os
# import json
# from Bio.PDB import PDBParser, PPBuilder
# '''
# 提取peptide2结构文件序列到新的json文件中,
# peptide字段是peptides_sota的generated_peptide序列,
# sequence字段是peptide2结构文件序列
# '''
# # 文件路径
# json_path = '/root/autodl-tmp/PP_generate_v1/data/SOTA/peptides_sota.json'
# pdb_folder = '/root/autodl-fs/RFdiffusion1'
# output_json = '/root/autodl-tmp/PP_generate_v1/data/output/RFdiffusion_sequence.json'

# # 加载原始json
# with open(json_path, 'r') as f:
#     data = json.load(f)

# output = {}

# for pdb_id in data:
#     pdb_file = os.path.join(pdb_folder, f"{pdb_id}.pdb")
#     entry = {}

#     # 复制原json里peptide/generate_peptide字段（防止名字不同）
#     if "peptide" in data[pdb_id]:
#         entry["peptide"] = data[pdb_id]["peptide"]
#     elif "generated_peptide" in data[pdb_id]:
#         entry["peptide"] = data[pdb_id]["generated_peptide"]
#     else:
#         entry["peptide"] = []

#     # 如果pdb文件存在，提取序列
#     if os.path.exists(pdb_file):
#         try:
#             parser = PDBParser(QUIET=True)
#             structure = parser.get_structure(pdb_id, pdb_file)
#             ppb = PPBuilder()
#             sequences = []
#             for model in structure:
#                 for chain in model:
#                     for pp in ppb.build_peptides(chain):
#                         seq = str(pp.get_sequence())
#                         sequences.append(seq)
#                     break  # 只取第一个chain
#                 break  # 只取第一个model
#             # 拼接所有片段
#             entry["sequence"] = ''.join(sequences)
#         except Exception as e:
#             print(f"Error parsing {pdb_file}: {e}")
#             entry["sequence"] = ""
#     else:
#         entry["sequence"] = ""
#         print(f"Warning: {pdb_file} not found!")

#     output[pdb_id] = entry

# # 保存到新json
# with open(output_json, 'w') as f:
#     json.dump(output, f, indent=4)

# -----------------------------------------------------------------------------------------------------------


#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量把 PDB 提取成 FASTA

用法（示例）：
    python pdb2fasta.py \
        --json  /root/autodl-tmp/PP_generate_v1/data/SOTA/peptides_sota.json \
        --pdbdir /root/autodl-tmp/peptides2 \
        --out    /root/autodl-tmp/peptides2/all_peptides.fa
"""

import argparse, json, os, sys
from pathlib import Path
from typing import Dict

from Bio.PDB import PDBParser

# ----------- 三字母→一字母 ----------------------------------------------------
AA3_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O","MSE":"M"
}

def pdb_to_seq(pdb_path: Path) -> str:
    """返回该 PDB 中最长链的一级序列（若无标准残基则抛错）"""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("S", pdb_path)

    chains = {}
    for ch in next(struct.get_models()):        # 只取 model 0
        aa = [AA3_1.get(res.get_resname().upper(), "")
               for res in ch if res.id[0] == " "]
        if aa:
            chains[ch.id] = "".join(aa)
    if not chains:
        raise ValueError(f"No standard residues in {pdb_path}")
    return max(chains.values(), key=len)


# ------------------------- CLI & 主流程 ---------------------------------------
def main(args):
    # 1. 读 JSON，拿到所有 ID
    with open(args.json) as f:
        meta: Dict[str, Dict] = json.load(f)
    pdb_ids = list(meta.keys())
    print(f"JSON 中共有 {len(pdb_ids)} 个 PDB ID")

    # 2. 遍历，每个提序列
    records = {}
    missing = []

    for pid in pdb_ids:
        pdb_file = Path(args.pdbdir) / f"{pid}/0004.pdb"
        if not pdb_file.is_file():
            missing.append(pid)
            continue
        try:
            records[pid] = pdb_to_seq(pdb_file)
        except Exception as e:
            print(f"[WARN] {pid} 解析失败: {e}", file=sys.stderr)

    print(f"成功提取 {len(records)} 条序列；缺失或出错 {len(missing)} 条")

    # 3-A. 写到一个大 FASTA
    if args.out:
        with open(args.out, "w") as fout:
            for pid, seq in records.items():
                fout.write(f">{pid}\n{seq}\n")
        print("已写入:", args.out)

    # 3-B. (可选) 每条单独写 FASTA
    if args.split_dir:
        Path(args.split_dir).mkdir(parents=True, exist_ok=True)
        for pid, seq in records.items():
            with open(Path(args.split_dir)/f"{pid}.fa", "w") as f:
                f.write(f">{pid}\n{seq}\n")
        print("单独 FASTA 写入目录:", args.split_dir)

    if missing:
        print("以下 ID 未找到或解析失败：", ", ".join(missing))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch PDB → FASTA")
    parser.add_argument("--json", default='/root/autodl-tmp/PP_generate_v1/data/SOTA/peptides_sota.json', help="包含 PDB ID 的 JSON 文件")
    parser.add_argument("--pdbdir", default='/root/autodl-tmp/ppflow', help="存放 *.pdb 的目录")
    parser.add_argument("--out", default="/root/autodl-tmp/PP_generate_v1/data/Diversity/ppflow_sequence5.fasta",
                        help="合并输出 FASTA 路径；设为空则跳过")
    parser.add_argument("--split_dir", default="",
                        help="若指定则把每条序列单独写到该目录")
    main(parser.parse_args())
