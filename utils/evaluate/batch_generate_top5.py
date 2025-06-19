import torch
import argparse
import sys
import torch.nn.functional as F
import os
import json
sys.path.append("/root/autodl-tmp/PP_generate_v1/")
from models.peptide_generation_eval_v9 import ProteinPeptideModel
from models.esm.utils.structure.protein_chain import ProteinChain
from models.esm.utils.encoding import tokenize_sequence
from torch.nn.functional import softmax
from models.esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
from tokenizers import Tokenizer
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"
os.environ["TORCH_USE_CUDA_DSA"] = "1"

def load_model(model_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = ProteinPeptideModel(device)
    model.load_state_dict(torch.load(model_path), strict=False)
    model.to(device)
    model.eval()
    return model

def extract_sequence_from_pdb(pdb_path):
    try:
        protein_chain = ProteinChain.from_pdb(pdb_path)
        return protein_chain.sequence if protein_chain else None
    except Exception as e:
        print(f"Error reading {pdb_path}: {e}")
        return None

def prepare_sequence(sequence, tokenizer, device):
    if sequence:
        sequence_tensor = tokenize_sequence(sequence, tokenizer).to(device)
        return sequence_tensor.unsqueeze(0), sequence  # Add batch dimension
    return None, None
def top_k_top_p_filtering(logits, top_k=0, top_p=0.0, filter_value=-float('Inf')):
    assert logits.dim() == 1  # logits are one-dimensional
    sorted_logits, sorted_indices = torch.sort(logits, descending=True)
    cumulative_probs = torch.cumsum(F.softmax(sorted_logits, dim=-1), dim=0)

    # Remove tokens with cumulative probability above the threshold (nucleus sampling)
    sorted_indices_to_remove = cumulative_probs > top_p
    sorted_indices_to_remove[1:] = sorted_indices_to_remove[:-1].clone()
    sorted_indices_to_remove[0] = False
    
    # Remove tokens by setting their logits to a large negative value
    if top_k > 0:
        sorted_indices_to_remove[top_k:] = True  # Top-k filtering
    indices_to_remove = sorted_indices[sorted_indices_to_remove]
    logits[indices_to_remove] = filter_value
    return logits
def generate_peptide(model, receptor_seq_tensor, tokenizer, device, max_length=30, top_k=10, top_p=0.95):
    batch = {}
    batch['receptor_seq_tensor'] = receptor_seq_tensor
    batch['input_ids'] = torch.tensor([tokenizer.encode("1").ids]).squeeze(0).to(device)
    
    with torch.no_grad():
        for step in range(max_length):
            # input_ids = torch.tensor(generated_ids).to(device)
            peptide_logits = model(batch, mode = 'test')
            next_token_logits = peptide_logits[-1, :]
            filtered_logits = top_k_top_p_filtering(next_token_logits, top_k=top_k, top_p=top_p)
            next_token_probs = F.softmax(filtered_logits, dim=-1)
            next_token_id = torch.multinomial(next_token_probs, num_samples=1)
            batch['input_ids'] = torch.cat([batch['input_ids'], next_token_id], dim=-1)
            if next_token_id == tokenizer.token_to_id('<|eos|>'):
                break
    generated_sequence = [tokenizer.decode(row.detach().cpu().numpy().tolist()) for row in batch['input_ids'].unsqueeze(0)]
    return generated_sequence

def batch_generate_peptides(model, pdb_folder, output_file, replicate_id):
    sequence_tokenizer = EsmSequenceTokenizer()
    decoder_tokenizer = Tokenizer.from_file("/root/autodl-tmp/PP_generate_v1/logs/progen2/tokenizer_pep.json")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    peptides_data = {}
    pdb_ids = [name for name in os.listdir(pdb_folder) if os.path.isdir(os.path.join(pdb_folder, name))]
    
    for pdb_id in pdb_ids:
        receptor_pdb_path = os.path.join(pdb_folder, f"{pdb_id}/receptor.pdb")
        peptide_pdb_path = os.path.join(pdb_folder, f"{pdb_id}/peptide.pdb")

        sequence = extract_sequence_from_pdb(receptor_pdb_path)
        test_sequence = extract_sequence_from_pdb(peptide_pdb_path)
        receptor_seq_tensor, _ = prepare_sequence(sequence, sequence_tokenizer, device)

        if receptor_seq_tensor is not None:
            peptide = generate_peptide(model, receptor_seq_tensor, decoder_tokenizer, device)
            peptides_data[pdb_id] = {
                "protein_sequence": sequence,
                "generated_peptide": peptide,
                "test_peptide": test_sequence,
            }

    # 输出文件名加入 replicate_id
    output_path = output_file.replace(".json", f"_rep{replicate_id}.json")
    with open(output_path, 'w') as f:
        json.dump(peptides_data, f, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_folder", type=str, default="/root/autodl-tmp/PP_generate_v1/data/PPDbench", help="Protein dir")
    parser.add_argument("--output_file", type=str, default="/root/autodl-tmp/PP_generate_v1/data/top_k/peptides_final.json", help="Base output JSON file name")
    parser.add_argument("--model_path", type=str, default="/root/autodl-fs/pp_esm3/best_model_epoch_25_loss_0.4112.pth", help="Model path")
    args = parser.parse_args()

    model = load_model(args.model_path)
    model.eval()

    # 生成5次，每次保存一个JSON
    for i in range(1, 6):
        batch_generate_peptides(model, args.pdb_folder, args.output_file, replicate_id=i)

