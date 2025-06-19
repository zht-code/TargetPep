#!/usr/bin/env python3
import os
import glob
import argparse
'''将proteinMPNN恢复RFdiffusion的序列整合到一个fasta文件中'''
def extract_and_merge_seqs(input_dir: str, output_fasta: str):
    """
    从 input_dir 中的所有 .fa 文件提取第4行序列，
    并按 FASTA 格式写入 output_fasta。
    """
    fa_paths = sorted(glob.glob(os.path.join(input_dir, "*_5.fa")))
    if not fa_paths:
        print(f"[!] 在 {input_dir} 中未找到任何 .fa 文件。")
        return

    with open(output_fasta, "w") as out_f:
        for fa in fa_paths:
            name = os.path.splitext(os.path.basename(fa))[0]  # 去掉 .fa
            with open(fa) as f:
                lines = f.readlines()
            if len(lines) < 4:
                print(f"[!] 文件 {fa} 少于 4 行，跳过。")
                continue
            seq = lines[3].strip()
            if not seq:
                print(f"[!] 文件 {fa} 第4行为空,跳过。")
                continue

            # 写入合并文件
            out_f.write(f">{name}\n")
            out_f.write(f"{seq}\n")

    print(f"[✔] 已将 {len(fa_paths)} 个文件的序列合并到 {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="从多个 .fa 文件提取第4行并合并为一个 FASTA 文件"
    )
    parser.add_argument(
        "-i", "--input_dir", default='/root/autodl-tmp/ProteinMPNN/outputs/RFdiffusion_top7/seqs',
        help="包含 .fa 文件的文件夹路径"
    )
    parser.add_argument(
        "-o", "--output_fasta",default='/root/autodl-tmp/PP_generate_v1/data/Diversity/RFdiffusion_sequence555.fasta',
        help="输出的合并 fasta 文件路径"
    )

    args = parser.parse_args()
    extract_and_merge_seqs(args.input_dir, args.output_fasta)
