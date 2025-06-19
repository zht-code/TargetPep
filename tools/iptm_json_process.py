# import json
# import os

# # 假设你原始文件名为 input.json
# with open("/root/autodl-tmp/quvina.json", "r") as f:
#     data = json.load(f)

# # 创建输出目录（可选）
# os.makedirs("/root/autodl-tmp/quvina", exist_ok=True)

# for entry_id, content in data.items():
#     protein_seq = content["protein_sequence"]
#     peptide_seq = content["generated_peptide"][0]  # 只取第一个生成的多肽

#     output_data = [{
#         "sequences": [
#             {
#                 "proteinChain": {
#                     "sequence": protein_seq,
#                     "count": 1
#                 }
#             },
#             {
#                 "proteinChain": {
#                     "sequence": peptide_seq,
#                     "count": 1
#                 }
#             }
#         ],
#         "name": entry_id
#     }]

#     # 保存到以 ID 命名的 JSON 文件
#     output_path = os.path.join("/root/autodl-tmp/quvina", f"{entry_id}.json")
#     with open(output_path, "w") as out_f:
#         json.dump(output_data, out_f, indent=4)

# print("✅ 所有 JSON 文件已生成完成。")



import json
import os
'''根据json和fasra文件'''
# 载入 JSON 数据
with open("/root/autodl-tmp/quvina.json", "r") as f:
    data = json.load(f)

# 载入参考 FASTA 文件
fasta_path = "/root/autodl-tmp/reference.fasta"
fasta_sequences = {}
with open(fasta_path, "r") as fasta_file:
    current_id = None
    current_seq = []
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                fasta_sequences[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0].lower()  # 统一转小写与 entry_id 匹配
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        fasta_sequences[current_id] = "".join(current_seq)

# 创建输出目录
os.makedirs("/root/autodl-tmp/quvina", exist_ok=True)

# 构建新的 JSON 文件
for entry_id, content in data.items():
    # 使用 reference.fasta 中的蛋白序列（如果存在）
    if entry_id.lower() in fasta_sequences:
        protein_seq = fasta_sequences[entry_id.lower()]
    else:
        print(f"⚠️ 未找到 entry_id={entry_id} 的参考序列，使用原始序列。")
        protein_seq = content["protein_sequence"]

    peptide_seq = content["generated_peptide"][0]  # 只取第一个生成的多肽

    output_data = [{
        "sequences": [
            {
                "proteinChain": {
                    "sequence": protein_seq,
                    "count": 1
                }
            },
            {
                "proteinChain": {
                    "sequence": peptide_seq,
                    "count": 1
                }
            }
        ],
        "name": entry_id
    }]

    output_path = os.path.join("/root/autodl-tmp/quvina", f"{entry_id}.json")
    with open(output_path, "w") as out_f:
        json.dump(output_data, out_f, indent=4)

print("✅ 所有 JSON 文件已生成完成。")
