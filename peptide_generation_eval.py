import torch
import torch.nn as nn
import torch
import sys
# 确保正确加载模块路径
sys.path.append("/root/autodl-tmp/PP_generate_v1/logs/esm_mol/weights/")
from tokenizers import Tokenizer
from models.esm.pretrained import load_local_model
from models.esm.utils.encoding import tokenize_sequence
from models.esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
import torch.nn.functional as F
from models.progen2.ProGen2.src.models.progen.attr_cross_modeling_progen import ProGenForCausalLM
import json
# 定义一个类来转换字典为对象属性
class Config:
    def __init__(self, data):
        for key, value in data.items():
            setattr(self, key, value)

def load_config(config_path):
    """ 加载JSON配置文件 """
    with open(config_path, 'r', encoding='utf-8') as file:
        config = json.load(file)
    return Config(config)
# 创建序列掩码
def create_sequence_mask(seq):
    seq_len = seq.size(0)
    mask = torch.triu(torch.ones((seq_len, seq_len)), diagonal=1)
    return mask # (seq_len, seq_len)

# 定义特征融合函数
# class FeatureFusionLayer(nn.Module):
#     def __init__(self, fused_dim):
#         super(FeatureFusionLayer, self).__init__()
#         self.fc = nn.Linear(1664, fused_dim)
    
#     def forward(self, esm_features, struct_features):
#         # 将ESM特征和3D结构特征在最后一维拼接
#         fused_features = torch.cat((esm_features, struct_features), dim=-1)
#         # 使用线性层融合
#         fused_features = self.fc(fused_features)
#         return fused_features
class MLP(nn.Module):
    def __init__(
        self, intermediate_size, config
    ):  # in MLP: intermediate_size= 4 * embed_dim
        super().__init__()
        embed_dim = config.embed_dim
        self.mlp = nn.Sequential(
            # nn.Linear(embed_dim, intermediate_size),
            nn.LeakyReLU(0.01) ,
            nn.Dropout(config.resid_pdrop),
            nn.Linear(intermediate_size, embed_dim),
            nn.LeakyReLU(0.01) ,
            nn.Dropout(config.resid_pdrop),
            nn.LayerNorm(embed_dim, config.layer_norm_epsilon)
        )
    def reparameterize(self, mu, log_sigma):
        """Reparameterization trick to sample from N(mu, sigma^2)."""
        std = torch.exp(0.5 * log_sigma)  # sigma = exp(0.5 * log_sigma)
        eps = torch.randn_like(std)
        return mu + eps * std  # z = mu + sigma * epsilon
    
    def forward(self, hidden_states):
        '''执行编码器前向传递，返回重新参数化的潜变量 z。'''
        mu = self.mlp(hidden_states)
        log_sigma = self.mlp(hidden_states)
        z = self.reparameterize(mu, log_sigma).mean(dim=1, keepdim=True)  # 取平均，维度变为 [1, 1, 1536]
        return z, mu, log_sigma
# 定义整体模型
class ProteinPeptideModel(nn.Module):
    def __init__(self, device):
        super(ProteinPeptideModel, self).__init__()
        # 加载encoder:esm3模型以及tokenizer
        esm3_model = load_local_model(model_name="esm3_sm_open_v1")
        self.esm_tokenizer = EsmSequenceTokenizer()
        self.device = device
        self.esm3_model = esm3_model
        # 加载decoder:progen2-medium以及tokenizer
        self.decoder = ProGenForCausalLM.from_pretrained("/root/autodl-tmp/PP_generate_v1/logs/progen2").to(device)
        self.decoder_tokenizer = Tokenizer.from_file("/root/autodl-tmp/PP_generate_v1/logs/progen2/tokenizer_pep.json")
        vocab_size = self.decoder_tokenizer.get_vocab_size()
        config = load_config('/root/autodl-tmp/PP_generate_v1/logs/progen2/config.json')
        self.config = config
        self.attr = MLP(config.embed_dim, config)
        self.embedding_layer = nn.Embedding(num_embeddings=vocab_size, embedding_dim=config.embed_dim)
        self.pocket_embed = nn.Linear(3, config.embed_dim)
        self.vina_embed = nn.Linear(1, config.embed_dim)
        self.solubility_embed = nn.Linear(1, config.embed_dim)
        self.stability_embed = nn.Linear(1, config.embed_dim)
        # self.fusion_layer = FeatureFusionLayer(fused_dim=fusion_dim)
        # self.cross_attention = CrossAttentionLayer(config)
        # self.affinity_predictor = nn.Sequential(
        #     nn.Linear(config.embed_dim, 1),
        #     nn.Sigmoid()  # 将输出限制在 [0, 1] 范围内
        # )
        # 冻结esm3的全部参数，只解冻最后两层
        for param in self.esm3_model.parameters():
            param.requires_grad = False
        for param in self.decoder.parameters():
            param.requires_grad = False
        # for param in self.esm3_model.transformer.norm.parameters():
        #     param.requires_grad = True
        for param in self.esm3_model.transformer.blocks[-3:].parameters():
            param.requires_grad = True
        # 冻结decoder的全部参数，只解冻最后两层,解冻cross_attention层
        # for param in self.decoder.parameters():
        #     param.requires_grad = False
        for param in self.decoder.transformer.h[-2:].parameters():
            param.requires_grad = True
        # for block in self.decoder.transformer.h:
        #     for param in block.cross_attn.parameters():
        #         param.requires_grad = True
        # 解冻decoder中cross_attn层的参数
        for param in self.decoder.transformer.ln_f.parameters():
            param.requires_grad = True
        for block in self.decoder.transformer.h:
            for param in block.attr_encode.parameters():
                param.requires_grad = True
        for block in self.decoder.transformer.h:
            for param in block.ln_encode.parameters():
                param.requires_grad = True
        for block in self.decoder.transformer.h:
            for param in block.attr_cross_attn.parameters():
                param.requires_grad = True
        for block in self.decoder.transformer.h:
            for param in block.cross_attn.parameters():
                param.requires_grad = True
        # for param in self.decoder.transformer.ln_2.parameters():
        #     param.requires_grad = True
        # for param in self.decoder.transformer.ln_3.parameters():
        #     param.requires_grad = True
        for param in self.decoder.lm_head.parameters():
            param.requires_grad = True
        # 亲和力预测器层
        # for param in self.decoder.ln_encode.parameters():
        #     param.requires_grad = True
        # for param in self.decoder.cross_attn.parameters():
        #     param.requires_grad = True
        # for param in self.decoder.ln_aff.parameters():
        #     param.requires_grad = True
        # for param in self.decoder.affinity_predictor.parameters():
        #     param.requires_grad = True

    
    def forward(self, batch ,mode='train'):
        # 1. 使用ESM3编码蛋白质序列
        esm_output = self.esm3_model(sequence_tokens=batch['receptor_seq_tensor'])
        esm_features = esm_output.embeddings
        # esm_features = esm_output.sequence_logits
        # 调整ESM特征以适应decoder的输入尺寸
        # esm_features_resized = F.interpolate(esm_features.permute(0, 2, 1), size=self.config.embed_dim, mode='linear').permute(0, 2, 1)
        # past_key_values = None
        if mode == 'train':
            # 训练模式，处理多肽的属性和序列
            pep_ids = self.embedding_layer(batch['peptide_seq_tensor'])
            pocket_coords = self.pocket_embed(batch['binding_pocket'])
            affinity = self.vina_embed(batch['vina_affinity'].unsqueeze(-1))
            # stability = self.stability_embed(batch['stability'].unsqueeze(-1))
            solubility = self.solubility_embed(batch['solubility'].unsqueeze(-1))
            # 所有属性特征融合
            # attributes = torch.cat([pocket_coords, affinity, solubility], dim=1)
            # 去掉亲和力属性
            attributes = torch.cat([pocket_coords, solubility], dim=1)
            # 去掉稳定性属性
            # attributes = torch.cat([pocket_coords, affinity, solubility], dim=1)
            # 去掉溶解性属性
            # attributes = torch.cat([pocket_coords, affinity], dim=1)
            
            # 以平均值扩展属性
            mean_attributes = attributes.mean(dim=1, keepdim=True)  # 取平均，维度变为 [1, 1, 1536]
            mean_attributes = mean_attributes.expand_as(pep_ids[:,1:-1,:])
            '''
            计算 KL 散度,用于衡量编码后的潜在变量 z 的分布与其先验分布之间的差异,这个KL散度损失函数是VAE模型中非常关键的正则化项
            用于确保学习到的潜在空间具有良好的结构性，使得模型可以通过采样潜在变量来生成具有多样性的新数据。这个损失推动潜在空间的分布向
            先验分布（通常是高斯分布）靠拢，有助于生成过程的稳定性和多样性。
            '''
            # z, mu, log_sigma = self.attr(attributes)  # 取平均，维度变为 [1, 1, 1536]
            # expanded_attributes = z.expand_as(pep_ids[:,1:-1,:])
            # kld_loss = -0.5 * torch.sum(1 + log_sigma - mu.pow(2) - log_sigma.exp())
            # kld_loss = kld_loss / pep_ids.size(0)  # 按批次平均
            # 计算加权和
            # alpha = 0.5  # 权重因子可以根据需要调整
            # pep_input_id = pep_ids.clone()
            # pep_input_id[:,1:,:] = alpha * pep_ids[:,1:,:] + (1 - alpha) * expanded_attributes
            # 创建attention_mask
            # self_attention_mask_length = pep_input_id.size(1)  # 总长度
            # cross_attention_mask_length = pep_input_id.size(1)  # 总长度
            # attention_mask = torch.ones((pep_input_id.size(1), self_attention_mask_length), dtype=torch.float, device=self.device)
            # attention_mask = create_sequence_mask(batch['peptide_seq_tensor']).to(self.device)
            # cross_attention_mask = torch.zeros((1, cross_attention_mask_length, esm_features.size(1)), dtype=torch.float, device=self.device)
            attr_cross_attention_mask = torch.zeros((batch['peptide_seq_tensor'].size(1), mean_attributes.size(1)), dtype=torch.float, device=self.device)
            cross_attention_mask = torch.zeros((batch['peptide_seq_tensor'].size(1), esm_features.size(1)), dtype=torch.float, device=self.device)
            # 计算需要被忽略的属性嵌入的长度
            # ignore_length = pocket_coords.size(1) + affinity.size(1) + stability.size(1) + solubility.size(1)
            # 将属性部分的mask设置为0,将多肽序列第一个以及最后一个掩盖掉
            attr_cross_attention_mask[0, :] = -1e9
            attr_cross_attention_mask[-1, :] = -1e9
            cross_attention_mask[0, :] = -1e9
            cross_attention_mask[-1, :] = -1e9
            labels = batch
            generated_peptide = self.decoder(input_ids=batch['peptide_seq_tensor'],
                                             labels = labels,
                                             attributes = mean_attributes,
                                            #  ignore_length = ignore_length,
                                              encoder_embeddings = esm_features,
                                            #   past_key_values = past_key_values,
                                            #   attention_mask = attention_mask,
                                              attr_cross_attention_mask = attr_cross_attention_mask,
                                              cross_attention_mask = cross_attention_mask,
                                            #   output_hidden_states=True
                                            )
            loss = generated_peptide.loss
            # total_loss = loss + 0.1 * kld_loss
            # past_key_values = generated_peptide.past_key_values
            # 检查是否有 NaN
            # if torch.isnan(generated_peptide_logit).any(): 
            #     print("NaN detected in predicted_peptide_logits")
            # peptide_pro = generated_peptide.hidden_states
            # cross_attn_output = self.cross_attention(
            #     query=peptide_pro,
            #     key=esm_features,
            #     value=esm_features
            # )
            # # 使用cross-attention的输出进行亲和力预测
            # affinity_score = self.affinity_predictor(cross_attn_output[:, -1, :]) 
            # return generated_peptide_logit, affinity_score
            return loss

        elif mode == 'test':
            # peptide_seq = "<|bos|>"
            # input_ids = torch.tensor(self.decoder_tokenizer.encode(peptide_seq).ids).unsqueeze(0).to(self.device)
            # input_embeds = self.embedding_layer(batch['input_ids'] )
            attention_mask = create_sequence_mask(batch['input_ids']).to(self.device)
            cross_attention_mask = torch.zeros((31, esm_features.size(1)), dtype=torch.float, device=self.device)
            # cross_attention_mask = torch.zeros((batch['input_ids'].size(0), esm_features.size(1)), dtype=torch.float, device=self.device)
            # 测试模式，生成多肽序列
            cross_attention_mask[:1, :] = -1e9
            generated_peptide_ids = self.decoder(
                input_ids = batch['input_ids'].unsqueeze(0),  # 在生成时无需提供input_ids
                encoder_embeddings=esm_features,
                # attention_mask = attention_mask,
                cross_attention_mask = cross_attention_mask,
                # past_key_values = past_key_values,
                # use_cache=True,
                # max_new_tokens=50,  # 假设多肽的最大长度
                # max_length = 50,
                # # num_beams=5,    # 使用beam search搜索的宽度
                # # num_return_sequences=3,   # 生成序列的数量
                # no_repeat_ngram_size=2     # 防止重复的n-gram
                
            )
            # past_key_values = generated_peptide_ids.past_key_values
            generated_peptide_ligit = generated_peptide_ids.logits
            # 可以选择解码生成的token以获取实际的氨基酸序列
            # generated_peptides = [self.tokenizer.decode(g, skip_special_tokens=True) for g in generated_peptide_ids]
            
            return generated_peptide_ligit.squeeze(0)

class CrossAttentionLayer(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.attention = nn.MultiheadAttention(config.embed_dim, 8)
        self.linear = nn.Linear(config.embed_dim, config.embed_dim)
        self.dropout = nn.Dropout(config.resid_pdrop)
        self.norm = nn.LayerNorm(config.embed_dim)

    def forward(self, query, key, value, key_padding_mask=None):
        # MultiheadAttention expects [seq_len, batch_size, embedding_dim]
        query = query.transpose(0, 1)  
        key = key.transpose(0, 1)  
        value = value.transpose(0, 1)

        attn_output, _ = self.attention(query, key, value, key_padding_mask=key_padding_mask)
        attn_output = attn_output.transpose(0, 1)  # back to [batch_size, seq_len, embedding_dim]
        query = query.transpose(0, 1)
        
        # Apply linear layer and normalization
        output = self.dropout(self.linear(attn_output))
        output = self.norm(output + query)  # Add & Norm step

        return output