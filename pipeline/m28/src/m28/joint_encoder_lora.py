"""Fixed last-four q/v low-rank residuals for the pinned GENERanno encoder."""
import math
import torch
from torch import nn
from torch.nn import functional as F
from torch.nn.attention import SDPBackend, sdpa_kernel

LAYERS = (22, 23, 24, 25)
PROJECTIONS = ('q_proj', 'v_proj')
RANK = 4
SCALE = 1.0
SEED = 0
PARAMETERS = 102400


class LowRankLinear(nn.Module):
    """BF16 base output plus an FP32 low-rank residual rounded to BF16."""
    def __init__(self, base, generator):
        super().__init__()
        self.base = base.requires_grad_(False)
        a = torch.empty(RANK, base.in_features, dtype=torch.float32)
        bound = 1 / math.sqrt(base.in_features)
        a.uniform_(-bound, bound, generator=generator)
        self.A = nn.Parameter(a.to(base.weight.device))
        self.B = nn.Parameter(torch.zeros(base.out_features, RANK,
            dtype=torch.float32, device=base.weight.device))

    def forward(self, x):
        value = self.base(x)
        residual = F.linear(F.linear(x.float(), self.A), self.B) * SCALE
        return value + residual.to(value.dtype)


def attach_adapters(backbone):
    config = backbone.config
    if (len(backbone.layers) != 26 or config.hidden_size != 2048
            or config.num_attention_heads != 32 or config.num_key_value_heads != 4
            or config.pretraining_tp != 1 or config._attn_implementation != 'sdpa'):
        raise ValueError('Pinned 26-layer GENERanno SDPA interface required')
    bases = {}
    for index in LAYERS:
        for projection in PROJECTIONS:
            module = getattr(backbone.layers[index].self_attn, projection)
            out_features = 2048 if projection == 'q_proj' else 256
            if (not isinstance(module, nn.Linear) or module.in_features != 2048
                    or module.out_features != out_features or module.bias is not None
                    or module.weight.dtype != torch.bfloat16):
                raise ValueError('Original BF16 q/v projection geometry required')
            bases[(index, projection)] = module
    backbone.eval().requires_grad_(False)
    generator = torch.Generator(device='cpu').manual_seed(SEED)
    adapters = {}
    for (index, projection), module in bases.items():
        replacement = LowRankLinear(module, generator)
        setattr(backbone.layers[index].self_attn, projection, replacement)
        adapters[f'layers.{index}.self_attn.{projection}'] = replacement
    return adapters


def adapter_parameters(adapters):
    return [p for module in adapters.values() for p in (module.A, module.B)]


def adapter_state(adapters):
    return {name+'.'+key: getattr(module, key).detach().cpu().clone()
        for name, module in adapters.items() for key in ('A', 'B')}


@torch.no_grad()
def restore_adapters(adapters, state):
    expected = {name+'.'+key: getattr(module, key)
        for name, module in adapters.items() for key in ('A', 'B')}
    if set(state) != set(expected):
        raise ValueError('All sixteen named adapter tensors required')
    for name, parameter in expected.items():
        value = state[name]
        if (value.shape != parameter.shape or value.dtype != torch.float32
                or not bool(torch.isfinite(value).all())):
            raise ValueError('Matching finite FP32 adapter state required: '+name)
    for name, parameter in expected.items():
        parameter.copy_(state[name].to(parameter.device))


def encode_joint_window(backbone, tokenizer, sequence, device):
    """Same joint24576 literal encoding, without the frozen inference decorator."""
    if not 0 < len(sequence) <= 24576:
        raise ValueError('Nonempty original-width oriented window required')
    padded = sequence+'N'*(24576-len(sequence))
    ids = tokenizer(padded, add_special_tokens=False, return_tensors='pt',
        truncation=False)['input_ids']
    if tuple(ids.shape) != (1, 4096):
        raise ValueError('All4096 six-mer tokens without specials or truncation required')
    ids = ids.to(device)
    with sdpa_kernel(SDPBackend.FLASH_ATTENTION):
        hidden = backbone(input_ids=ids, attention_mask=torch.ones_like(ids)).last_hidden_state
    if tuple(hidden.shape) != (1, 4096, 2048) or not bool(torch.isfinite(hidden).all()):
        raise ValueError('Finite full-window native features required')
    return hidden.to(torch.bfloat16).float()
