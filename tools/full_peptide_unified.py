import json

# 输入文件路径
base_file = "/root/autodl-tmp/PP_generate_v1/data/output/peptides_final1_hdock.json"
ppflow_file = "/root/autodl-tmp/PP_generate_v1/data/output/Hdock_generate_vina_cross1.json"
rf_file = "/root/autodl-tmp/PP_generate_v1/data/output/Hdock_RFdiffusion_vina1.json"
output_file = "/root/autodl-tmp/PP_generate_v1/data/SOTA/peptides_sota.json"

# 加载所有JSON文件
with open(base_file, 'r') as f:
    base_data = json.load(f)
with open(ppflow_file, 'r') as f:
    ppflow_data = json.load(f)
with open(rf_file, 'r') as f:
    rf_data = json.load(f)

# 合并字段
for pid in base_data:
    if "properties" in base_data[pid] and isinstance(base_data[pid]["properties"], list):
        props_list = base_data[pid]["properties"]
        if props_list:
            props = props_list[0]

            # 提取 Affinity 值
            ppflow_affinity = ppflow_data.get(pid, {}).get("properties", [{}])[0].get("PPFlow Affinity (kcal/mol)", None)
            rf_affinity = rf_data.get(pid, {}).get("properties", [{}])[0].get("RFdiffusion Affinity (kcal/mol)", None)

            # 插入字段
            new_props = {}
            for key, val in props.items():
                new_props[key] = val
                if key == "generate Affinity (kcal/mol)":
                    if ppflow_affinity is not None:
                        new_props["PPFlow Affinity (kcal/mol)"] = ppflow_affinity
                    else:
                        print(f"[缺失] PPFlow affinity not found for ID: {pid}")
                    if rf_affinity is not None:
                        new_props["RFdiffusion Affinity (kcal/mol)"] = rf_affinity
                    else:
                        print(f"[缺失] RFdiffusion affinity not found for ID: {pid}")

            # 更新属性
            base_data[pid]["properties"][0] = new_props
        else:
            print(f"[跳过] Empty properties list for ID: {pid}")
    else:
        print(f"[跳过] No properties found or not a list for ID: {pid}")

# 写入新文件
with open(output_file, 'w') as f:
    json.dump(base_data, f, indent=4)

print(f"合并完成，输出保存至: {output_file}")
