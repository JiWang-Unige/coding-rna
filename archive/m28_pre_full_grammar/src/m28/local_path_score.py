"""Fixed complete-CDS role emission gain; no head, labels or calibration."""
import torch
from .phase_evidence_chain import candidate_roles


@torch.no_grad()
def local_path_scores(local_logits, chains, valid_bases=None):
    """Sum log P(candidate role) - log P(intergenic) over each gene span.

    Softmax normalizers cancel exactly. Roles use the original ascending phase
    convention and start/stop precedence. This is a fixed emission potential,
    not a calibrated gene probability or a full HMM likelihood. Float64 only
    accumulates the frozen local logits; the model still runs in FP32.
    """
    length = len(local_logits) if valid_bases is None else int(valid_bases)
    if local_logits.ndim != 2 or local_logits.shape[1] != 15 or not 0 < length <= len(local_logits):
        raise ValueError('Expected valid oriented 15-class local logits')
    logits = local_logits[:length].double()
    if not bool(torch.isfinite(logits).all()):
        raise ValueError('Nonfinite valid local emission')
    rows = []
    for chain in chains:
        roles = candidate_roles(chain, length, logits.device)
        a, b = chain[0][0], chain[-1][1]
        positions = torch.arange(a, b, device=logits.device)
        rows.append((logits[positions, roles[a:b]]-logits[positions, 0]).sum())
    return torch.stack(rows) if rows else logits.new_empty((0,))
