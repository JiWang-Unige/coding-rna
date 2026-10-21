"""A2 residual head: C0 with zero-initialised, BIAS-FREE role layer (366,176 parameters)."""
import torch
import a2_common as C2     # noqa: F401
import a1_model as MD


def new_model(seed=0):
    torch.manual_seed(seed)
    m = MD.G.M28Core('C0', 2048, 128)
    torch.nn.init.zeros_(m.segmentation.weight)
    m.segmentation.bias = None
    assert sum(p.numel() for p in m.parameters()) == 366191 - 15
    return m
