import json

# 加载JSON文件
with open('/root/autodl-tmp/PP_generate_v1/data/output/Hdock_generate_vina_cross1.json', 'r') as f:
    data = json.load(f)

# 收集所有generate Affinity值
affinities = []

for entry in data.values():
    properties = entry.get("properties", [])
    for prop in properties:
        affinity = prop.get("generate Affinity (kcal/mol)")
        # affinity = prop.get("protein_generator Affinity (kcal/mol)")
        if affinity is not None:
            affinities.append(affinity)

# 计算平均值
if affinities:
    average_affinity = sum(affinities) / len(affinities)
    print(f"平均 generate Affinity (kcal/mol): {average_affinity:.2f}")
else:
    print("没有找到任何 'generate Affinity (kcal/mol)' 字段。")
