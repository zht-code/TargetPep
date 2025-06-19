import os
import torch
import json
from torch.utils.data import Dataset, DataLoader
from tokenizers import Tokenizer
from models.esm.utils.structure.protein_chain import ProteinChain
from models.esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
from models.esm.utils.encoding import tokenize_sequence
# 自定义Dataset类
class ProteinPeptideDataset(Dataset):
    def __init__(self, id_list, pdb_dir, json_files):
        """
        Args:
        - id_list: 一个ID列表，例如 ['1a0m_A', '2b0n_B']，每个ID对应一对多肽和蛋白
        - peptide_pdb_dir: 多肽PDB文件所在的目录
        - protein_pdb_dir: 蛋白PDB文件所在的目录
        """
        self.id_list = id_list
        self.pdb_dir = pdb_dir
        self.decoder_tokenizer = Tokenizer.from_file("/root/autodl-tmp/PP_generate_v1/logs/progen2/tokenizer_pep.json")
        # with open("/home/zht/programs/PP_generate_v1/data/stability_results.json", 'r') as f:
        #     self.peptide_attributes = json.load(f)

        # 初始化一个字典来存储多个JSON文件的多肽属性
        self.peptide_attributes = {}
        
        # 依次加载多个JSON文件，将内容合并到self.peptide_attributes
        for json_file in json_files:
            with open(json_file, 'r') as f:
                json_data = json.load(f)
                # 合并json数据到self.peptide_attributes
                for key, value in json_data.items():
                    if key not in self.peptide_attributes:
                        self.peptide_attributes[key] = {}
                    self.peptide_attributes[key].update(value)
        # Preprocess and filter out empty sequences
        # self.valid_samples = []
        # for sample_id in self.id_list:
        #     peptide_pdb_file = os.path.join(self.pdb_dir, f"{sample_id}/peptide.pdb")
        #     protein_pdb_file = os.path.join(self.pdb_dir, f"{sample_id}/receptor.pdb")
        #     peptide_seq = extract_sequence_from_pdb(peptide_pdb_file)
        #     receptor_seq = extract_sequence_from_pdb(protein_pdb_file)
        #     if peptide_seq and receptor_seq:
        #         self.valid_samples.append(sample_id)

    def __len__(self):
        # return len(self.valid_samples)
        return len(self.id_list)

    def __getitem__(self, idx):
        # 获取ID
        # sample_id = self.valid_samples[idx]
        sample_id = self.id_list[idx]

        # 构建PDB文件路径
        peptide_pdb_file = os.path.join(self.pdb_dir, f"{sample_id}/peptide.pdb")
        protein_pdb_file = os.path.join(self.pdb_dir, f"{sample_id}/receptor.pdb")

        # 提取多肽和蛋白的序列
        receptor_seq = extract_sequence_from_pdb(protein_pdb_file)  # 形状 [num_residues_peptide, feature_dim]
        peptide_seq = extract_sequence_from_pdb(peptide_pdb_file)  # 形状 [num_residues_protein, feature_dim]
        if peptide_seq is None:
            raise ValueError("Missing peptide sequence for index {}".format(sample_id))
        peptide_seq = "<|bos|>" + peptide_seq + "<|eos|>" 
        # 创建序列tokenizer
        sequence_tokenizer = EsmSequenceTokenizer()
        # 蛋白序列张量token
        receptor_seq_tensor = tokenize_sequence(receptor_seq, sequence_tokenizer)
        # 多肽序列张量token
        peptide_seq_tensor = torch.tensor(self.decoder_tokenizer.encode(peptide_seq).ids)
        # 获取多肽的属性特征
        peptide_attr = self.peptide_attributes.get(sample_id, {})
        binding_pocket_tensor = torch.tensor(peptide_attr.get('interaction_windows', [[0.0, 0.0, 0.0]]))
        # MAX_LENGTH = 28
        # # 确保长度一致性
        # if len(binding_pocket_data) < MAX_LENGTH:
        #     # 用0填充到最大长度
        #     binding_pocket_data.extend([[0.0, 0.0, 0.0]] * (MAX_LENGTH - len(binding_pocket_data)))
        # elif len(binding_pocket_data) > MAX_LENGTH:
        #     # 如果数据超过最大长度，进行截断
        #     binding_pocket_data = binding_pocket_data[:MAX_LENGTH]
        # binding_pocket_tensor = torch.tensor(binding_pocket_data)
        vina_affinity_tensor = torch.tensor(peptide_attr.get('score', 0.0)).unsqueeze(0)
        stability_tensor = torch.tensor(peptide_attr.get('stability_score', 0.0)).unsqueeze(0)
        solubility_tensor = torch.tensor(peptide_attr.get('solubility_score', 0.0)).unsqueeze(0)
        # if binding_pocket_tensor.size(0) == 0:
        #     print(f"Warning: No binding pocket data available for sample {sample_id}.")
        # if stability_tensor.nelement() == 0 or solubility_tensor.nelement() == 0:
        #     print(f"Warning: Stability or solubility score missing for sample {sample_id}.")

        return {
            "id": sample_id,
            "peptide_seq": peptide_seq,
            "receptor_seq": receptor_seq,
            "peptide_seq_tensor": peptide_seq_tensor,
            "receptor_seq_tensor": receptor_seq_tensor,
            "binding_pocket": binding_pocket_tensor.to(torch.float),
            "vina_affinity": vina_affinity_tensor.to(torch.float),
            "stability": stability_tensor.to(torch.float),
            "solubility": solubility_tensor.to(torch.float)
        }

# 加载PDB文件并提取序列
# def extract_sequence_from_pdb(pdb_path):
#     protein_chain = ProteinChain.from_pdb(pdb_path)
#     return protein_chain.sequence
def extract_sequence_from_pdb(pdb_path):
    try:
        protein_chain = ProteinChain.from_pdb(pdb_path)
        if protein_chain is None or not protein_chain.sequence:
            return None
        return protein_chain.sequence
    except Exception as e:
        print(f"Error reading {pdb_path}: {e}")
        return None