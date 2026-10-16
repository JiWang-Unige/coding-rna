"""Fixed q/v LoRA plus last-four attention-output residuals; same base encoder."""
import torch
from torch import nn
from src.m28 import joint_encoder_lora as QV

LAYERS = QV.LAYERS
PROJECTIONS = ('q_proj', 'v_proj', 'o_proj')
RANK, SCALE, SEED = QV.RANK, QV.SCALE, QV.SEED
PARAMETERS = 167936
LowRankLinear = QV.LowRankLinear
adapter_parameters = QV.adapter_parameters
adapter_state = QV.adapter_state
encode_joint_window = QV.encode_joint_window


def attach_adapters(backbone):
    # Inspect the extra target before the established q/v injector mutates anything.
    outputs = {}
    for index in LAYERS:
        module = backbone.layers[index].self_attn.o_proj
        if (not isinstance(module, nn.Linear) or module.in_features != 2048
                or module.out_features != 2048 or module.bias is not None
                or module.weight.dtype != torch.bfloat16):
            raise ValueError('Original bias-free BF162048x2048 output projections required')
        outputs[index] = module
    adapters = QV.attach_adapters(backbone)
    # Keep every existing q/v initial A/B exactly as before; o has its own seed0 stream.
    generator = torch.Generator(device='cpu').manual_seed(SEED)
    for index, module in outputs.items():
        replacement = LowRankLinear(module, generator)
        backbone.layers[index].self_attn.o_proj = replacement
        adapters[f'layers.{index}.self_attn.o_proj'] = replacement
    return {f'layers.{index}.self_attn.{projection}':
        adapters[f'layers.{index}.self_attn.{projection}']
        for index in LAYERS for projection in PROJECTIONS}


@torch.no_grad()
def restore_adapters(adapters, state):
    expected = {name+'.'+key: getattr(module, key)
        for name, module in adapters.items() for key in ('A', 'B')}
    if set(state) != set(expected):
        raise ValueError('All twenty-four named q/v/o adapter tensors required')
    for name, parameter in expected.items():
        value = state[name]
        if (value.shape != parameter.shape or value.dtype != torch.float32
                or not bool(torch.isfinite(value).all())):
            raise ValueError('Matching finite FP32 q/v/o adapter state required: '+name)
    for name, parameter in expected.items():
        parameter.copy_(state[name].to(parameter.device))
