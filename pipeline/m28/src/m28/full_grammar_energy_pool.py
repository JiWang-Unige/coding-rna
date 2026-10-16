"""Fixed equal logarithmic pool of aligned original fifteen-role energies."""
import torch


@torch.no_grad()
def equal_role_energy_pool(left, right):
    """Average two CPU FP32 role arrays in FP64, without modifying either.

    For a shared grammar this averages complete-set energies, including each
    model's background subtraction. It does not average decoded chains or
    marginal probabilities. DNA and orientation alignment belong to the caller.
    """
    if (left.device.type != 'cpu' or right.device.type != 'cpu'
            or left.dtype != torch.float32 or right.dtype != torch.float32
            or left.ndim != 2 or left.shape[1] != 15 or left.shape != right.shape):
        raise ValueError('Two aligned CPU FP32 N by fifteen role arrays required')
    result = left.double()  # Fresh copy; promote before adding FP32 operands.
    return result.add_(right).mul_(0.5)


@torch.no_grad()
def minimum_centered_role_energy_pool(left, right):
    """Fixed per-base agreement penalty, not a minimum of whole-gene scores.

    Subtract each branch's own background in FP64 before taking the minimum.
    Thus every local score is the mean centered energy minus half the absolute
    branch disagreement. Background remains exactly zero; the original arrays
    remain unchanged. No fitted coefficient or posterior threshold is used.
    """
    if (left.device.type != 'cpu' or right.device.type != 'cpu'
            or left.dtype != torch.float32 or right.dtype != torch.float32
            or left.ndim != 2 or left.shape[1] != 15 or left.shape != right.shape):
        raise ValueError('Two aligned CPU FP32 N by fifteen role arrays required')
    a, b = left.double(), right.double()
    # Clone the background column because its view overlaps the in-place target.
    a.sub_(a[:, :1].clone()); b.sub_(b[:, :1].clone())
    return torch.minimum(a, b)
