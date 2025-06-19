# coding=utf-8
# Copyright 2021 The EleutherAI and HuggingFace Teams. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License atí
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Modified forward-pass implementation based on https://github.com/huggingface/transformers/blob/main/src/transformers/models/gptj/modeling_gptj.py

from typing import Tuple

import numpy as np

import torch
import torch.utils.checkpoint
from torch import nn
from torch.nn import CrossEntropyLoss
import torch.nn.functional as F
import math
from transformers.activations import ACT2FN
from transformers.modeling_outputs import (
    BaseModelOutputWithPast,
    CausalLMOutputWithPast,
)
from transformers.modeling_utils import PreTrainedModel
from transformers.utils import logging
from transformers.utils.model_parallel_utils import assert_device_map, get_device_map
from .configuration_progen import ProGenConfig


logger = logging.get_logger(__name__)


def fixed_pos_embedding(x, seq_dim=1, seq_len=None):
    dim = x.shape[-1]
    if seq_len is None:
        seq_len = x.shape[seq_dim]
    inv_freq = 1.0 / (10000 ** (torch.arange(0, dim, 2) / dim))
    sinusoid_inp = (
        torch.einsum("i , j -> i j", torch.arange(seq_len), inv_freq)
        .to(x.device)
        .float()
    )
    return torch.sin(sinusoid_inp), torch.cos(sinusoid_inp)


def rotate_every_two(x: torch.Tensor):
    x1 = x[:, :, :, ::2]
    x2 = x[:, :, :, 1::2]
    x = torch.stack((-x2, x1), axis=-1)
    return x.flatten(-2)

def apply_rotary_pos_emb(x, sincos, offset=0):
    sin, cos = map(
        lambda t: t[None, offset : x.shape[1] + offset, None, :].repeat_interleave(
            2, 3
        ),
        sincos,
    )
    # einsum notation for lambda t: repeat(t[offset:x.shape[1]+offset,:], "n d -> () n () (d j)", j=2)
    return (x * cos) + (rotate_every_two(x) * sin)


class ProGenAttention(nn.Module):
    def __init__(self, config):
        super().__init__()

        max_positions = config.n_positions
        self.register_buffer(
            "bias",
            torch.tril(
                torch.ones((max_positions, max_positions), dtype=torch.bool)
            ).view(1, 1, max_positions, max_positions),
            persistent=False
        )
        self.register_buffer("masked_bias", torch.tensor(-1e9), persistent=False)     # approx. -inf

        self.attn_dropout = nn.Dropout(config.attn_pdrop)
        self.resid_dropout = nn.Dropout(config.resid_pdrop)

        self.embed_dim = config.embed_dim
        self.num_attention_heads = config.n_head
        self.head_dim = self.embed_dim // self.num_attention_heads
        if self.head_dim * self.num_attention_heads != self.embed_dim:
            raise ValueError(
                f"embed_dim must be divisible by num_attention_heads (got `embed_dim`: {self.embed_dim} and `num_attention_heads`: {self.num_attention_heads})."
            )
        self.scale_attn = torch.sqrt(
            torch.tensor(self.head_dim, dtype=torch.float32)
        ).to(torch.get_default_dtype())
        self.qkv_proj = nn.Linear(self.embed_dim, self.embed_dim * 3, bias=False)

        self.out_proj = nn.Linear(self.embed_dim, self.embed_dim, bias=False)
        self.rotary_dim = None
        if config.rotary_dim is not None:
            self.rotary_dim = config.rotary_dim

    def _split_heads(self, x: torch.Tensor, n_head, dim_head) -> torch.Tensor:
        x = x.reshape(x.shape[:-2] + (-1,))                             # (B, T, 8 * E // 8)
        x = x.reshape(x.shape[:-1] + (n_head, dim_head))             # (B, T, n_heads, dim_head)
        return x

    def _merge_heads(self, tensor, num_attention_heads, attn_head_size) -> torch.Tensor:
        """
        Merges attn_head_size dim and num_attn_heads dim into n_positions
        """
        if len(tensor.shape) == 5:
            tensor = tensor.permute(0, 1, 3, 2, 4).contiguous()
        elif len(tensor.shape) == 4:
            tensor = tensor.permute(0, 2, 1, 3).contiguous()
        else:
            raise ValueError(
                f"Input tensor rank should be one of [4, 5], but is: {len(tensor.shape)}"
            )
        new_shape = tensor.size()[:-2] + (num_attention_heads * attn_head_size,)
        return tensor.view(new_shape)

    # def _attn(
    #     self,
    #     query,
    #     key,
    #     value,
    #     attention_mask=None,
    #     head_mask=None,
    # ):
    #     # compute causal mask from causal mask buffer
    #     query_length, key_length = query.size(-2), key.size(-2)
    #     causal_mask = self.bias[
    #         :, :, key_length - query_length : key_length, :key_length
    #     ]

    #     # Keep the attention weights computation in fp32 to avoid overflow issues
    #     query = query.to(torch.float32)
    #     key = key.to(torch.float32)

    #     attn_weights = query @ key.transpose(-1, -2)                 # (B, n_heads, T, T)

    #     attn_weights = attn_weights / self.scale_attn

    #     # attend only to previous positions
    #     attn_weights = torch.where(
    #         causal_mask, attn_weights, self.masked_bias.to(attn_weights.dtype)
    #     )

    #     if attention_mask is not None:
    #         attn_weights = attn_weights + attention_mask

    #     attn_weights = F.softmax(attn_weights, dim=-1)
    #     attn_weights = attn_weights.to(value.dtype)
    #     attn_weights = self.attn_dropout(attn_weights)

    #     if head_mask is not None:
    #         attn_weights = attn_weights * head_mask

    #     attn_output = attn_weights @ value                # (B, n_heads, T, dim_head)

    #     return attn_output, attn_weights
    def _attn(
        self,
        query,
        key,
        value,
        attention_mask=None,
        head_mask=None,
    ):
        # 获取维度信息
        batch_size, num_heads, query_length, head_dim = query.size()
        key_length = key.size(-2)

        # 生成正确的因果遮掩
        causal_mask = self.bias[:, :, :query_length, :key_length].to(query.device)

        # 计算注意力得分
        attn_scores = torch.matmul(query, key.transpose(-1, -2))  # (B, n_heads, Q_len, K_len)
        attn_scores = attn_scores / self.scale_attn

        # 应用因果遮掩
        attn_scores = attn_scores.masked_fill(causal_mask == 0, self.masked_bias)

        # 如果提供了 attention_mask，添加到得分上
        if attention_mask is not None:
            attn_scores = attn_scores + attention_mask

        # 计算注意力权重
        attn_weights = F.softmax(attn_scores, dim=-1)
        attn_weights = self.attn_dropout(attn_weights)

        # 如果提供了 head_mask，应用到注意力权重上
        if head_mask is not None:
            attn_weights = attn_weights * head_mask

        # 计算注意力输出
        attn_output = torch.matmul(attn_weights, value)  # (B, n_heads, Q_len, head_dim)

        return attn_output, attn_weights

    def forward(
        self,
        hidden_states,
        attention_mask=None,
        layer_past=None,
        head_mask=None,
        use_cache=False,
        output_attentions=False,
    ):
        qkv = self.qkv_proj(hidden_states)                                         # (B, T, 3 * E)

        mp_num = 8
        qkv_split = qkv.reshape(qkv.shape[:-1] + (mp_num, -1))                     # (B, T, 8, 3 * E // 8)

        query, key, value = torch.split(qkv_split, self.embed_dim // mp_num, dim=-1) # 3 * (B, T, 8, E // 8)

        query = self._split_heads(query, self.num_attention_heads, self.head_dim)  # (B, T, n_heads, dim_head)
        key = self._split_heads(key, self.num_attention_heads, self.head_dim)      # (B, T, n_heads, dim_head)
        value = self._split_heads(value, self.num_attention_heads, self.head_dim)  # (B, T, n_heads, dim_head)
        value = value.permute(0, 2, 1, 3)

        seq_len = key.shape[1]
        offset = 0

        if layer_past is not None:
            offset = layer_past[0].shape[-2]
            seq_len = seq_len + offset

        if self.rotary_dim is not None:
            k_rot = key[:, :, :, : self.rotary_dim]
            k_pass = key[:, :, :, self.rotary_dim :]

            q_rot = query[:, :, :, : self.rotary_dim]
            q_pass = query[:, :, :, self.rotary_dim :]

            sincos = fixed_pos_embedding(k_rot, 1, seq_len=seq_len)
            k_rot = apply_rotary_pos_emb(k_rot, sincos, offset=offset)
            q_rot = apply_rotary_pos_emb(q_rot, sincos, offset=offset)

            key = torch.cat([k_rot, k_pass], dim=-1)
            query = torch.cat([q_rot, q_pass], dim=-1)
        else:
            sincos = fixed_pos_embedding(key, 1, seq_len=seq_len)
            key = apply_rotary_pos_emb(key, sincos, offset=offset)
            query = apply_rotary_pos_emb(query, sincos, offset=offset)

        key = key.permute(0, 2, 1, 3)
        query = query.permute(0, 2, 1, 3)

        if layer_past is not None:
            past_key = layer_past[0]
            past_value = layer_past[1]
            key = torch.cat((past_key, key), dim=-2)
            value = torch.cat((past_value, value), dim=-2)

        if use_cache is True:
            present = (key, value)
        else:
            present = None

        # compute self-attention: softmax((Q @ K.T) / sqrt(dim_head)) @ V
        attn_output, attn_weights = self._attn(
            query, key, value, attention_mask, head_mask
        )

        attn_output = self._merge_heads(                               # (B, T, E)  
            attn_output, self.num_attention_heads, self.head_dim
        )

        attn_output = self.out_proj(attn_output)
        attn_output = self.resid_dropout(attn_output)

        outputs = (attn_output, present)
        if output_attentions:
            outputs = outputs + (attn_weights,)

        return outputs  # a, present, (attentions)

def softmax(x, axis=-1):
    """Softmax函数，用于计算注意力权重"""
    e_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
    return e_x / e_x.sum(axis=axis, keepdims=True)
 
def scaled_dot_product_attention(q, k, v, mask=None):
    """缩放点积注意力机制，用于计算输出和注意力权重"""
    # print(q.shape)
    # print(k.transpose().shape)
    matmul_qk = torch.matmul(q, k.transpose(3,2))  # 计算查询和键的矩阵乘积
    d_k = k.shape[-1]  # 键的维度
    scaled_attention_logits = matmul_qk / math.sqrt(d_k)  # 缩放注意力分数
 
    if mask is not None:  # 如果有注意力掩码，将其添加到分数上
        scaled_attention_logits += (mask * -1e9)
 
    attention_weights = F.softmax(scaled_attention_logits, dim=-1)  # 计算注意力权重
    output = torch.matmul(attention_weights, v)  # 计算输出
    return output, attention_weights

class ProGenCrossAttention(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.num_attention_heads = config.n_head
        self.head_dim = config.embed_dim // self.num_attention_heads
        # self.multihead_attn = nn.MultiheadAttention(embed_dim=config.embed_dim, num_heads=self.num_attention_heads)

        self.query_proj = nn.Linear(config.embed_dim, config.embed_dim)
        self.key_proj = nn.Linear(config.embed_dim, config.embed_dim)
        self.value_proj = nn.Linear(config.embed_dim, config.embed_dim)
        # self.out_proj = nn.Linear(config.embed_dim, config.embed_dim)

        # self.dropout = nn.Dropout(config.attn_pdrop)

    def forward(self, query, key, value, attention_mask=None):
        # Define batch size and sequence lengths
        batch_size = query.size(0)
        seq_len_query = query.size(1)
        seq_len_key = key.size(1)
        # Split the embedding dimension into multiple heads
        query = self.query_proj(query).view(-1, query.size(-2), self.num_attention_heads, self.head_dim).transpose(1, 2)
        key = self.key_proj(key).view(-1, key.size(-2), self.num_attention_heads, self.head_dim).transpose(1, 2)
        value = self.value_proj(value).view(-1, value.size(-2), self.num_attention_heads, self.head_dim).transpose(1, 2)

        # # Calculate attention scores
        # attention_scores = torch.matmul(query_layer, key_layer.transpose(-1, -2))
        # attention_scores = attention_scores / math.sqrt(self.head_dim)
        # if attention_mask is not None:
        #     attention_scores = attention_scores + attention_mask

        # # Apply softmax to get probabilities
        # attention_probs = F.softmax(attention_scores, dim=-1)
        # attention_probs = self.dropout(attention_probs)

        # # Multiply weights by values
        # context_layer = torch.matmul(attention_probs, value_layer)

        # # Concatenate heads and project back to the embedding dimension
        # context_layer = context_layer.transpose(1, 2).contiguous().view(-1, context_layer.size(-2), self.num_attention_heads * self.head_dim)
        # attention_output = self.out_proj(context_layer)


        # query: [seq_len, batch_size, embed_size] from decoder
        # key, value: [seq_len, batch_size, embed_size] from encoder
        # key_padding_mask: [batch_size, seq_len] if you have padding tokens in the encoder outputs
        
        # Apply cross attention
        # Apply cross attention
        # 调整维度顺序
        # query = query.transpose(0, 1)  # [seq_len_query, batch_size, embed_dim]
        # key = key.transpose(0, 1)      # [seq_len_key_value, batch_size, embed_dim]
        # value = value.transpose(0, 1)  # [seq_len_key_value, batch_size, embed_dim]

        
        # attn_output, attn_output_weights = self.multihead_attn(
        #     query, key, value, 
        #     # key_padding_mask=key_padding_mask,
        #     attn_mask=attention_mask  # Apply additional attention mask if provided
        # )
        attention_output, attention_weights = scaled_dot_product_attention(query, key, value, attention_mask)
        attn_output = attention_output.transpose(1, 2).contiguous().view(batch_size, seq_len_query, -1)
        # attn_output = attn_output.squeeze(0)
        return attn_output          # [batch_size，seq_len_query, embed_dim]


class ProGenMLP(nn.Module):
    def __init__(
        self, intermediate_size, config
    ):  # in MLP: intermediate_size= 4 * embed_dim
        super().__init__()
        embed_dim = config.embed_dim

        self.fc_in = nn.Linear(embed_dim, intermediate_size)
        self.fc_out = nn.Linear(intermediate_size, embed_dim)

        self.act = ACT2FN[config.activation_function]
        # self.act = nn.GELU()
        self.dropout = nn.Dropout(config.resid_pdrop)
        

    def forward(self, hidden_states):
        hidden_states = self.fc_in(hidden_states)
        hidden_states = self.act(hidden_states)
        hidden_states = self.fc_out(hidden_states)
        hidden_states = self.dropout(hidden_states)
        return hidden_states


class ProGenBlock(nn.Module):
    def __init__(self, config):
        super().__init__()
        inner_dim = config.n_inner if config.n_inner is not None else 4 * config.embed_dim
        # self.ln_1 = nn.LayerNorm(config.embed_dim, eps=config.layer_norm_epsilon)
        self.attn = ProGenAttention(config)
        # self.cross_attn = ProGenCrossAttention(config)
        self.mlp = ProGenMLP(inner_dim, config)
        # self.ln_2 = nn.LayerNorm(config.embed_dim, eps=config.layer_norm_epsilon)
        # self.ln_3 = nn.LayerNorm(config.embed_dim, eps=config.layer_norm_epsilon)
        # self.ln_4 = nn.LayerNorm(config.embed_dim, eps=config.layer_norm_epsilon)
        # self.cross_attn_weight = nn.Parameter(torch.tensor(0.5))  # 初始权重为0.5
        # self.self_attn_weight = nn.Parameter(torch.tensor(1.0))

    def forward(
        self,
        hidden_states,
        ignore_length,
        layer_past=None,
        encoder_hidden_states=None,  # Add encoder_hidden_states argument
        attention_mask=None,
        cross_attention_mask=None,
        head_mask=None,
        use_cache=False,
        output_attentions=False,
    ):
        residual = hidden_states
        # hidden_states = self.ln_1(hidden_states)
        attn_outputs = self.attn(
            hidden_states,
            layer_past = layer_past,
            attention_mask=attention_mask,
            head_mask=head_mask,
            use_cache=use_cache,
            output_attentions=output_attentions,
        )
        attn_output = attn_outputs[0]
        outputs = attn_outputs[1:]
        # hidden_states = attn_output + residual
        # hidden_states = residual + self.self_attn_weight * attn_output
        # hidden_states = hidden_states[:,ignore_length:,:]
        # 归一化蛋白序列的embedding
        # encoder_embeddings_normalized = encoder_hidden_states / encoder_hidden_states.norm(dim=1, keepdim=True)


        # hidden_states = self.ln_3(hidden_states)
        # encoder_hidden_states = self.ln_2(encoder_hidden_states)
        # if encoder_hidden_states is not None:
        #     cross_attn_output = self.cross_attn(hidden_states, encoder_hidden_states, encoder_hidden_states, cross_attention_mask)
        #     hidden_states = hidden_states + cross_attn_output  # Combine outputs from self-attn and cross-attn
        #     # hidden_states = hidden_states + self.cross_attn_weight * cross_attn_output
        # hidden_states = self.ln_4(hidden_states)
        feed_forward_hidden_states = self.mlp(hidden_states)
        hidden_states = attn_output + feed_forward_hidden_states + residual

        if use_cache:
            outputs = (hidden_states,) + attn_outputs
        else:
            outputs = (hidden_states,) + attn_outputs[1:]
        return outputs


class ProGenPreTrainedModel(PreTrainedModel):
    """
    An abstract class to handle weights initialization and a simple interface for downloading and loading pretrained
    models.
    """

    config_class = ProGenConfig
    base_model_prefix = "transformer"
    is_parallelizable = True

    def __init__(self, *inputs, **kwargs):
        super().__init__(*inputs, **kwargs)

    def _init_weights(self, module):
        """Initialize the weights."""
        if isinstance(module, (nn.Linear,)):
            # Slightly different from Mesh Transformer JAX which uses truncated_normal for initialization
            # cf https://github.com/pytorch/pytorch/pull/5617
            # module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            torch.nn.init.xavier_uniform_(module.weight)
            if module.bias is not None:
                module.bias.data.zero_()
        elif isinstance(module, nn.Embedding):
            # module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            torch.nn.init.xavier_uniform_(module.weight)
            if module.padding_idx is not None:
                module.weight.data[module.padding_idx].zero_()
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)


class ProGenModel(ProGenPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)
        self.vocab_size_emb = config.vocab_size_emb
        self.embed_dim = config.embed_dim
        self.max_position_embeddings = config.max_position_embeddings
        self.n_head = config.n_head
        self.wte = nn.Embedding(config.vocab_size_emb, self.embed_dim)
        self.wpe = nn.Embedding(config.n_positions, self.embed_dim)
        self.drop = nn.Dropout(config.embd_pdrop)
        self.h = nn.ModuleList([ProGenBlock(config) for _ in range(config.n_layer)])
        self.cross_attn = ProGenCrossAttention(config)
        self.ln_f = nn.LayerNorm(self.embed_dim, eps=config.layer_norm_epsilon)
        self.ln_2 = nn.LayerNorm(self.embed_dim, eps=config.layer_norm_epsilon)
        self.ln_3 = nn.LayerNorm(self.embed_dim, eps=config.layer_norm_epsilon)
        # self.ln_attr = nn.LayerNorm(self.embed_dim, eps=config.layer_norm_epsilon)
        self.input_weight = nn.Parameter(torch.tensor(1.0))
        self.rotary_dim = min(
            config.rotary_dim, config.n_positions // config.n_head
        )
        self.init_weights()

        # Model parallel
        self.model_parallel = False
        self.device_map = None

    def parallelize(self, device_map=None):
        # Check validity of device_map
        self.device_map = (
            get_device_map(len(self.h), range(torch.cuda.device_count())) if device_map is None else device_map
        )
        assert_device_map(self.device_map, len(self.h))
        self.model_parallel = True
        self.first_device = "cpu" if "cpu" in self.device_map.keys() else "cuda:" + str(min(self.device_map.keys()))
        self.last_device = "cuda:" + str(max(self.device_map.keys()))
        self.wte = self.wte.to(self.first_device)
        # Load onto devices
        for k, v in self.device_map.items():
            for block in v:
                cuda_device = "cuda:" + str(k)
                self.h[block] = self.h[block].to(cuda_device)
        # ln_f to last
        self.ln_f = self.ln_f.to(self.last_device)

    def deparallelize(self):
        self.model_parallel = False
        self.device_map = None
        self.first_device = "cpu"
        self.last_device = "cpu"
        self.wte = self.wte.to("cpu")
        for index in range(len(self.h)):
            self.h[index] = self.h[index].to("cpu")
        self.ln_f = self.ln_f.to("cpu")
        torch.cuda.empty_cache()

    def forward(
        self,
        input_ids=None,
        attributes = None,
        ignore_length = None,
        encoder_hidden_states = None,
        cross_attention_mask = None,
        past_key_values=None,
        attention_mask=None,
        token_type_ids=None,
        position_ids=None,
        head_mask=None,
        inputs_embeds=None,
        use_cache=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
    ):
        output_attentions = (
            output_attentions
            if output_attentions is not None
            else self.config.output_attentions
        )
        output_hidden_states = (
            output_hidden_states
            if output_hidden_states is not None
            else self.config.output_hidden_states
        )
        use_cache = use_cache if use_cache is not None else self.config.use_cache
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        if input_ids is not None and inputs_embeds is not None:
            raise ValueError(
                "You cannot specify both input_ids and inputs_embeds at the same time"
            )
        elif input_ids is not None:
            input_shape = input_ids.size()
            input_ids = input_ids.view(-1, input_shape[-1])
            batch_size = input_ids.shape[0]
        elif inputs_embeds is not None:
            input_shape = inputs_embeds.size()[:-1]
            batch_size = inputs_embeds.shape[0]
        else:
            raise ValueError("You have to specify either input_ids or inputs_embeds")

        device = input_ids.device if input_ids is not None else inputs_embeds.device

        if token_type_ids is not None:
            token_type_ids = token_type_ids.view(-1, input_shape[-1])

        if position_ids is not None:
            position_ids = position_ids.view(-1, input_shape[-1])

        if past_key_values is None:
            past_length = 0
            past_key_values = tuple([None] * len(self.h))
        else:
            past_length = past_key_values[0][0].size(-2)

        if position_ids is None:
            position_ids = torch.arange(
                past_length,
                inputs_embeds.size(1),
                dtype=torch.long,
                device=device,
            )
            position_ids = position_ids.unsqueeze(0).expand(inputs_embeds.size(0), -1)

        # Attention mask.
        if attention_mask is not None:
            assert batch_size > 0, "batch_size has to be defined and > 0"
            # attention_mask = attention_mask.view(batch_size, -1)
            # We create a 3D attention mask from a 2D tensor mask.
            # Sizes are [batch_size, 1, 1, to_seq_length]
            # So we can broadcast to [batch_size, num_heads, from_seq_length, to_seq_length]
            # this attention mask is more simple than the triangular masking of causal attention
            # used in OpenAI GPT, we just need to prepare the broadcast dimension here.
            attention_mask = attention_mask[None, None,:,  :]

            # Since attention_mask is 1.0 for positions we want to attend and 0.0 for
            # masked positions, this operation will create a tensor which is 0.0 for
            # positions we want to attend and -10000.0 for masked positions.
            # Since we are adding it to the raw scores before the softmax, this is
            # effectively the same as removing these entirely.
            attention_mask = attention_mask.to(dtype=self.dtype)  # fp16 compatibility
            attention_mask = (1.0 - attention_mask) * -10000.0
        if cross_attention_mask is not None:
                # 扩展掩码维度为 [num_heads, seq_len_query, seq_len_key]
            cross_attention_mask = cross_attention_mask.unsqueeze(0).unsqueeze(0)  # shape: [1, 1, seq_len_query, seq_len_key]
            cross_attention_mask = cross_attention_mask.expand(self.n_head, -1, cross_attention_mask.size(2), cross_attention_mask.size(3))  # shape: [num_heads, 1, seq_len_query, seq_len_key]
            cross_attention_mask = cross_attention_mask.squeeze(1)  # 去掉 batch 维度，shape: [num_heads, seq_len_query, seq_len_key]
        #     # 扩展维度以适配多头注意力的形状要求
        #     cross_attention_mask = cross_attention_mask.unsqueeze(0) # 扩展维度
        #     # 转换值，0变为不关注（-10000），1保持关注（0）
        #     cross_attention_mask = (1.0 - cross_attention_mask) * -10000.0
        #     # 转换数据类型（如果使用FP16或其他特殊类型）
        #     cross_attention_mask = cross_attention_mask.to(dtype=self.dtype)
        # Prepare head mask if needed
        # 1.0 in head_mask indicate we keep the head
        # attention_probs has shape bsz x num_attention_heads x N x N
        # head_mask has shape n_layer x batch x num_attention_heads x N x N
        head_mask = self.get_head_mask(head_mask, self.config.n_layer)

        if inputs_embeds is None:
            inputs_embeds = self.wte(input_ids)
            if position_ids is not None:
                position_embeds = self.wpe(position_ids)
            if attributes is not None :
                self.wpe = nn.Linear(self.max_position_embeddings, self.embed_dim).to(attributes.device)
                attributes = self.wpe(attributes)
                #    对属性嵌入在序列维度上进行平均池化
                pooled_attributes = attributes.mean(dim=1, keepdim=True)  # 形状：[1, 1, 1536]
                # 扩展维度以匹配多肽序列长度
                expanded_attributes = pooled_attributes.expand(-1, inputs_embeds.size(1), -1)
                # inputs_embeds = torch.cat([attributes, inputs_embeds], dim=1)
                # 为 inputs_embeds 添加权重
                # inputs_embeds = self.input_weight * inputs_embeds
            input_shape = inputs_embeds.size()[:-1]
        if position_ids is not None:
            position_embeds = self.wpe(position_ids)
        # if attributes is not None :
        #     self.wpe = nn.Linear(self.embed_dim, self.embed_dim).to(attributes.device)
        #     attributes = self.wpe(attributes)
        #     #    对属性嵌入在序列维度上进行平均池化
        #     pooled_attributes = attributes.mean(dim=1, keepdim=True)  # 形状：[1, 1, 1536]
        #     # 扩展维度以匹配多肽序列长度
        #     expanded_attributes = pooled_attributes.expand(-1, inputs_embeds.size(1), -1)

        hidden_states = inputs_embeds + position_embeds


        if token_type_ids is not None:
            token_type_embeds = self.wte(token_type_ids)
            hidden_states = hidden_states + token_type_embeds
            

        hidden_states = self.drop(hidden_states)

        output_shape = input_shape + (hidden_states.size(-1),)

        presents = () if use_cache else None
        all_self_attentions = () if output_attentions else None
        all_hidden_states = () if output_hidden_states else None
        for i, (block, layer_past) in enumerate(zip(self.h, past_key_values)):
            # Model parallel
            if self.model_parallel:
                torch.cuda.set_device(hidden_states.device)
                # Ensure layer_past is on same device as hidden_states (might not be correct)
                if layer_past is not None:
                    layer_past = tuple(past_state.to(hidden_states.device) for past_state in layer_past)
                # Ensure that attention_mask is always on the same device as hidden_states
                if attention_mask is not None:
                    attention_mask = attention_mask.to(hidden_states.device)
                if isinstance(head_mask, torch.Tensor):
                    head_mask = head_mask.to(hidden_states.device)

            if output_hidden_states:
                all_hidden_states = all_hidden_states + (hidden_states,)

            if getattr(self.config, "gradient_checkpointing", False) and self.training:
                if use_cache:
                    logger.warning(
                        "`use_cache=True` is incompatible with `config.gradient_checkpointing=True`. Setting "
                        "`use_cache=False`..."
                    )
                    use_cache = False

                def create_custom_forward(module):
                    def custom_forward(*inputs):
                        # None for past_key_value
                        return module(*inputs, use_cache, output_attentions)

                    return custom_forward

                outputs = torch.utils.checkpoint.checkpoint(
                    create_custom_forward(block),
                    hidden_states,
                    None,
                    attention_mask,
                    head_mask[i],
                )
            else:
                outputs = block(
                    hidden_states,
                    layer_past=layer_past,
                    ignore_length = ignore_length,
                    encoder_hidden_states = encoder_hidden_states,
                    cross_attention_mask = cross_attention_mask,
                    attention_mask=attention_mask,
                    head_mask=head_mask[i],
                    use_cache=use_cache,
                    output_attentions=output_attentions,
                )

            hidden_states = outputs[1]
            if use_cache is True:
                presents = presents + (outputs[1],)

            if output_attentions:
                all_self_attentions = all_self_attentions + (
                    outputs[2 if use_cache else 1],
                )

            # Model Parallel: If it's the last layer for that device, put things on the next device
            if self.model_parallel:
                for k, v in self.device_map.items():
                    if i == v[-1] and "cuda:" + str(k) != self.last_device:
                        hidden_states = hidden_states.to("cuda:" + str(k + 1))

        hidden_states = self.ln_f(hidden_states)
        encoder_hidden_states = self.ln_2(encoder_hidden_states)
        if encoder_hidden_states is not None:
            cross_attn_output = self.cross_attn(hidden_states, encoder_hidden_states, encoder_hidden_states, cross_attention_mask)
            hidden_states = hidden_states + cross_attn_output  # Combine outputs from self-attn and cross-attn
            # hidden_states = hidden_states + self.cross_attn_weight * cross_attn_output
        hidden_states = self.ln_3(hidden_states)
        hidden_states = hidden_states.view(*output_shape)
        # Add last hidden state
        if output_hidden_states:
            all_hidden_states = all_hidden_states + (hidden_states,)

        if not return_dict:
            return tuple(
                v
                for v in [
                    hidden_states,
                    presents,
                    all_hidden_states,
                    all_self_attentions,
                ]
                if v is not None
            )

        return BaseModelOutputWithPast(
            last_hidden_state=hidden_states,
            past_key_values=presents,
            hidden_states=all_hidden_states,
            attentions=all_self_attentions,
        )


class ProGenForCausalLM(ProGenPreTrainedModel):
    _keys_to_ignore_on_load_missing = [
        r"h\.\d+\.attn\.masked_bias",
        r"h\.\d+\.attn\.bias",
        r"lm_head\.weight",
    ]

    def __init__(self, config):
        super().__init__(config)
        self.transformer = ProGenModel(config)
        self.lm_head = nn.Linear(config.embed_dim, config.vocab_size_lm_head)
        self.init_weights()
        # Model parallel
        self.model_parallel = False
        self.device_map = None

    def parallelize(self, device_map=None):
        self.device_map = (
            get_device_map(len(self.transformer.h), range(torch.cuda.device_count()))
            if device_map is None
            else device_map
        )
        assert_device_map(self.device_map, len(self.transformer.h))
        self.transformer.parallelize(self.device_map)
        self.lm_head = self.lm_head.to(self.transformer.first_device)
        self.model_parallel = True

    def deparallelize(self):
        self.transformer.deparallelize()
        self.transformer = self.transformer.to("cpu")
        self.lm_head = self.lm_head.to("cpu")
        self.model_parallel = False
        torch.cuda.empty_cache()
    
    def prepare_inputs_for_generation(self, input_ids=None, input_embeds =None, encoder_hidden_states=None, attention_mask=None, cross_attention_mask=None, past_key_values=None, **kwargs):
        # 获取来自上一步的状态
        past_key_values = past_key_values[0] if past_key_values is not None else None

        esm_features = encoder_hidden_states

        # 获得注意力掩码
        attention_mask = attention_mask

        # 获得交叉注意力掩码，初始化或来自上一步
        cross_attention_mask = cross_attention_mask
        if cross_attention_mask is None:
            # 假设交叉注意力掩码全1，实际使用时根据需要调整
            cross_attention_mask = torch.ones((input_ids.size(0), esm_features.size(1)), dtype=torch.float, device=self.device)

        inputs = {
            "input_ids": None,
            "inputs_embeds":input_embeds,
            "past_key_values": past_key_values,
            "attention_mask": attention_mask,
            "encoder_hidden_states": esm_features,
            "cross_attention_mask": cross_attention_mask
        }
        return inputs

    def forward(
        self,
        input_ids=None,
        attributes = None,
        ignore_length = None,
        encoder_hidden_states = None,
        past_key_values=None,
        attention_mask=None,
        cross_attention_mask = None,
        token_type_ids=None,
        position_ids=None,
        head_mask=None,
        inputs_embeds=None,
        labels=None,
        use_cache=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
    ):
        r"""
        labels (:obj:`torch.LongTensor` of shape :obj:`(batch_size, sequence_length)`, `optional`):
            Labels for language modeling. Note that the labels **are shifted** inside the model, i.e. you can set
            ``labels = input_ids`` Indices are selected in ``[-100, 0, ..., config.vocab_size]`` All labels set to
            ``-100`` are ignored (masked), the loss is only computed for labels in ``[0, ..., config.vocab_size]``
        """
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        transformer_outputs = self.transformer(
            input_ids,
            attributes = attributes,
            ignore_length = ignore_length,
            encoder_hidden_states = encoder_hidden_states,
            cross_attention_mask = cross_attention_mask,
            past_key_values=past_key_values,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeds=inputs_embeds,
            use_cache=use_cache,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )
        hidden_states = transformer_outputs["last_hidden_state"]
        past_key_values = transformer_outputs["past_key_values"]

        # make sure sampling in fp16 works correctly and
        # compute loss in fp32 to match with mesh-tf version
        # https://github.com/EleutherAI/gpt-neo/blob/89ce74164da2fb16179106f54e2269b5da8db333/models/gpt2/gpt2.py#L179
        lm_logits = self.lm_head(hidden_states).to(torch.float32)

        loss = None
        if labels is not None:
            # Shift so that tokens < n predict n
            shift_logits = lm_logits[..., :-1, :].contiguous()
            shift_labels = labels[..., 1:].contiguous()
            loss_fct = CrossEntropyLoss()
            loss = loss_fct(
                shift_logits.view(-1, shift_logits.size(-1)), shift_labels.view(-1)
            )
            loss = loss.to(hidden_states.dtype)

        if not return_dict:
            output = (lm_logits,) + transformer_outputs[1:]
            return ((loss,) + output) if loss is not None else output

        return CausalLMOutputWithPast(
            loss=loss,
            logits=lm_logits,
            past_key_values=transformer_outputs["past_key_values"],
            hidden_states=transformer_outputs["last_hidden_state"],
            attentions=transformer_outputs.attentions,
        )

    @staticmethod
    def _reorder_cache(
        past: Tuple[Tuple[torch.Tensor]], beam_idx: torch.Tensor
    ) -> Tuple[Tuple[torch.Tensor]]:
        """
        This function is used to re-order the :obj:`past_key_values` cache if
        :meth:`~transformers.PretrainedModel.beam_search` or :meth:`~transformers.PretrainedModel.beam_sample` is
        called. This is required to match :obj:`past_key_values` with the correct beam_idx at every generation step.
        """
        # 验证 beam_idx 是否在正确的范围内
        # 假设你已经知道 beams 的数量，例如，通过配置或直接传参
        # num_beams = past[0][0].size(0)  # 假设 past_key_values 中的每个 tensor 第一维大小等于 beam 的数量
        # if torch.any(beam_idx >= num_beams) or torch.any(beam_idx < 0):
        #     raise ValueError(f"Beam index out of expected range [0, {num_beams-1}]")
        # if past[0][0].size(0) == 1 and num_beams > 1:
        #     # 扩展 past_key_values 以匹配每个 beam
        #     past = tuple(
        #         tuple(past_state.expand(num_beams, *past_state.shape[1:]) for past_state in layer_past)
        #         for layer_past in past
        #     )
        return tuple(
            tuple(
                past_state.index_select(0, beam_idx.to(past_state.device))
                for past_state in layer_past
            )
            for layer_past in past
        )
