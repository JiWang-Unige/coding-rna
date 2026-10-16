"""Zero-start affine role margins on the frozen A3 context and A31 potential."""
import torch
from torch import nn
from .local_path_score import local_path_scores
from .structured_emission import trainable_path_scores


class EmissionResidual(nn.Module):
    def __init__(self, hidden=128):
        super().__init__()
        # Intergenic correction is fixed zero; learn only the 14 role margins.
        self.margin = nn.Linear(hidden, 14)
        nn.init.zeros_(self.margin.weight)
        nn.init.zeros_(self.margin.bias)

    def forward(self, frozen_hidden, frozen_logits, chains, valid_bases=None):
        length = len(frozen_hidden) if valid_bases is None else int(valid_bases)
        if (frozen_hidden.ndim != 2 or frozen_logits.ndim != 2 or
                frozen_hidden.dtype != torch.float32 or frozen_logits.dtype != torch.float32 or
                frozen_hidden.requires_grad or frozen_logits.requires_grad or
                not 0 < length <= min(len(frozen_hidden), len(frozen_logits))):
            raise ValueError('Detached original oriented FP32 context and local logits required')
        margins = self.margin(frozen_hidden[:length])
        correction = torch.cat([margins.new_zeros((length, 1)), margins], dim=1)
        base = local_path_scores(frozen_logits, chains, length)
        residual = trainable_path_scores(correction, chains, length)
        # Add the separately accumulated potentials, not rounded FP32 logits.
        return dict(logits=base+residual, null=margins.sum()*0,
                    frozen_path_gain=base, residual_path_gain=residual,
                    residual_local_logits=correction)
