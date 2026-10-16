"""Pinned last-four q/v FP32 master weights with the original BF16 projection."""
import torch
from torch import nn
from torch.nn import functional as F
from src.m28.joint_encoder_lora import encode_joint_window

LAYERS = (22, 23, 24, 25)
PROJECTIONS = ('q_proj', 'v_proj')
PARAMETERS = 18874368
NATIVE_LR = 1e-5


class FullRankLinear(nn.Module):
    def __init__(self, base):
        super().__init__()
        base.requires_grad_(False)
        self.weight = nn.Parameter(base.weight.detach().float().clone())

    def forward(self, x):
        return F.linear(x, self.weight.to(x.dtype))


def attach_dense_qv(backbone):
    config = backbone.config
    if (len(backbone.layers) != 26 or config.hidden_size != 2048
            or config.num_attention_heads != 32 or config.num_key_value_heads != 4
            or config.pretraining_tp != 1 or config._attn_implementation != 'sdpa'):
        raise ValueError('Pinned26-layer GENERanno SDPA interface required')
    bases = {}
    for index in LAYERS:
        for projection in PROJECTIONS:
            module = getattr(backbone.layers[index].self_attn, projection)
            output = 2048 if projection == 'q_proj' else 256
            if (not isinstance(module, nn.Linear) or module.in_features != 2048
                    or module.out_features != output or module.bias is not None
                    or module.weight.dtype != torch.bfloat16):
                raise ValueError('Original BF16 bias-free q/v geometry required')
            bases[(index, projection)] = module
    backbone.eval().requires_grad_(False)
    projections = {}
    for (index, projection), base in bases.items():
        replacement = FullRankLinear(base)
        setattr(backbone.layers[index].self_attn, projection, replacement)
        projections[f'layers.{index}.self_attn.{projection}'] = replacement
    if sum(p.numel() for p in dense_parameters(projections)) != PARAMETERS:
        raise ValueError('Exactly18874368 selected q/v weight parameters required')
    return projections


def dense_parameters(projections):
    return [module.weight for module in projections.values()]


def dense_state(projections):
    return {name+'.weight': module.weight.detach().cpu().clone()
            for name, module in projections.items()}


@torch.no_grad()
def restore_dense(projections, state):
    expected = {name+'.weight': module.weight for name, module in projections.items()}
    if set(state) != set(expected):
        raise ValueError('All eight named full-rank q/v weights required')
    for name, parameter in expected.items():
        value = state[name]
        if (value.shape != parameter.shape or value.dtype != torch.float32
                or not bool(torch.isfinite(value).all())):
            raise ValueError('Matching finite FP32 q/v state required: '+name)
    for name, parameter in expected.items():
        parameter.copy_(state[name].to(parameter.device))
