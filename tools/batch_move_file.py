import shutil
import json
import os
from tqdm import tqdm
json_file = '/root/autodl-tmp/PP_generate_v1/data/output/final_sequence.json'
dir_path = '/root/quvina'
with open(json_file, 'r') as f:
    data = json.load(f)

for protein_id in tqdm(data.keys(), desc="Processing proteins"):
    pdb_file = os.path.join(dir_path, f"{protein_id}/{protein_id}.pdb")
    shutil.move(pdb_file, '/root/quvina')