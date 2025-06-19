import torch.nn as nn
import torch
import os
import warnings
import logging
import json
import psutil
MEMORY_THRESHOLD = 20 * 1024 ** 3 
LOSS_RECORD_FILE = "/root/autodl-tmp/PP_generate_v1/model_loss_record.json"
from ranger import Ranger
from models.esm.utils.structure.protein_chain import ProteinChain
import torch.nn.functional as F
from utils.datasets.dataset_process import DatasetProcess
from models.peptide_generation_eval import ProteinPeptideModel
from modules.Reinforce_learning import *
from transformers import get_cosine_schedule_with_warmup  
from modules.paddingcollate import *
from torch.optim.lr_scheduler import CosineAnnealingLR, ReduceLROnPlateau
from torch.optim.lr_scheduler import StepLR
def get_disk_available_gb(path="/root/autodl-tmp"):
    usage = psutil.disk_usage(path)
    available_gb = usage.free / (1024 ** 3)
    return available_gb


def update_model_loss_record(model_path, val_loss):
    record = {}
    if os.path.exists(LOSS_RECORD_FILE):
        with open(LOSS_RECORD_FILE, 'r') as f:
            record = json.load(f)
    record[model_path] = val_loss
    with open(LOSS_RECORD_FILE, 'w') as f:
        json.dump(record, f, indent=2)

def delete_worst_model():
    if not os.path.exists(LOSS_RECORD_FILE):
        return
    with open(LOSS_RECORD_FILE, 'r') as f:
        record = json.load(f)
    if not record:
        return
    worst_model = max(record.items(), key=lambda x: x[1])[0]
    try:
        os.remove(worst_model)
        logging.info(f"Deleted worst model due to low memory: {worst_model}")
        del record[worst_model]
        with open(LOSS_RECORD_FILE, 'w') as f:
            json.dump(record, f, indent=2)
    except Exception as e:
        logging.warning(f"Failed to delete model: {worst_model}, Error: {e}")

warnings.filterwarnings("ignore")


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

logging.basicConfig(filename='/root/autodl-tmp/PP_generate_v1/training.log', level=logging.INFO,
                    format='%(asctime)s:%(levelname)s:%(message)s')
def extract_sequence_from_pdb(pdb_path):
    protein_chain = ProteinChain.from_pdb(pdb_path)
    return protein_chain.sequence

def adaptive_gradient_clipping(parameters, clip_factor=0.01, eps=1e-3):
    for p in parameters:
        if p.grad is None:
            continue
        param_norm = torch.max(torch.norm(p, dim=-1, keepdim=True), torch.tensor(eps).to(p.device))
        grad_norm = torch.norm(p.grad, dim=-1, keepdim=True)
        max_norm = param_norm * clip_factor
        clip_coef = torch.min(max_norm / (grad_norm + eps), torch.tensor(1.).to(p.device))
        p.grad.data.mul_(clip_coef)


def train_model(data_loader, model_save_path):
    model = ProteinPeptideModel(device).to(device)
    optimizer = Ranger(model.parameters(), lr=1e-5)
    scheduler = StepLR(optimizer, step_size=3, gamma=0.5)
    
    min_loss = float('inf')
    total_epochs = 1000
    for epoch in range(total_epochs):
        model.train()
        total_loss = 0
        count = 0
        for i, batch in enumerate(data_loader):
            optimizer.zero_grad(set_to_none=True)
            batch['peptide_seq_tensor'] = batch['peptide_seq_tensor'].to(device)
            batch['binding_pocket'] = batch['binding_pocket'].to(device)
            batch['stability'] = batch['stability'].to(device)
            batch['solubility'] = batch['solubility'].to(device)
            batch['vina_affinity'] = batch['vina_affinity'].to(device)
            batch['receptor_seq_tensor'] = batch['receptor_seq_tensor'].to(device)
            loss = model(batch)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            
            total_loss += loss.item()
            count += 1
        avg_loss = total_loss / count
        logging.info(f"Epoch {epoch+1}, Average Loss: {avg_loss}")
        scheduler.step()
        logging.info(f"Learning rate after epoch {epoch+1}: {optimizer.param_groups[0]['lr']}")
        if avg_loss < min_loss:
            min_loss = avg_loss
            if get_disk_available_gb() < MEMORY_THRESHOLD:
                logging.warning("Memory low. Attempting to delete worst model before saving new one.")
                delete_worst_model()
            
            model_path = os.path.join(model_save_path, f"best_model_epoch_{epoch+1}_loss_{min_loss:.4f}.pth")
            torch.save(model.state_dict(), model_path)
            update_model_loss_record(model_path, min_loss)
            logging.info(f"Saved new best model with loss {min_loss:.4f}")

if __name__ == "__main__":
    pdb_dir = '/root/autodl-tmp/PP_generate_v1/data/ppbench2024'
    model_save_path = '/root/autodl-tmp/PP_generate_v1/logs/pp_esm3/qu_vina'
    json_files = [
        '/root/autodl-tmp/PP_generate_v1/data/final_vina_results.json',
        '/root/autodl-tmp/PP_generate_v1/data/solubility_results.json',
        '/root/autodl-tmp/PP_generate_v1/data/interaction_windows_results.json'
    ]
    dataset_processor = DatasetProcess(
        pdb_dir, 
        json_files, 
        batch_size=1,  
        collate_fn=PaddingCollate(pad_keys=['peptide_seq_tensor', 'receptor_seq_tensor', 'binding_pocket', 'vina_affinity', 'stability', 'solubility']),
        shuffle=False, 
        num_workers=4
    )
    data_loader = dataset_processor.create_data_loader()
    train_model(data_loader, model_save_path)







    