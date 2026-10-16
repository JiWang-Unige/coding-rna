"""Trainable version of the fixed complete-CDS emission potential."""
import torch
from .chain_set import chain_set_loss
from .phase_evidence_chain import candidate_roles


def trainable_path_scores(local_logits, chains, valid_bases=None):
    """Identical A31 unit role-vs-intergenic sum, retaining its gradients.

    The network remains FP32; FP64 accumulation is not a new role weight,
    length normalization, transition model or calibrated gene probability.
    """
    length = len(local_logits) if valid_bases is None else int(valid_bases)
    if local_logits.ndim != 2 or local_logits.shape[1] != 15 or not 0 < length <= len(local_logits):
        raise ValueError('Expected valid oriented 15-class local logits')
    logits = local_logits[:length].double()
    if not bool(torch.isfinite(logits).all()):
        raise ValueError('Nonfinite valid local emission')
    values = []
    for chain in chains:
        roles = candidate_roles(chain, length, logits.device)
        a, b = chain[0][0], chain[-1][1]
        positions = torch.arange(a, b, device=logits.device)
        values.append((logits[positions, roles[a:b]]-logits[positions, 0]).sum())
    return torch.stack(values) if values else logits.new_empty((0,))


def structured_emission_loss(local_logits, chains, targets, valid_bases=None):
    """Existing compatible-known-set likelihood, now through local emissions.

    Unknown candidates are excluded from supervision, not from later natural
    prediction. All-unknown/empty pools retain a differentiable zero.
    """
    gain = trainable_path_scores(local_logits, chains, valid_bases)
    return chain_set_loss({'logits': gain, 'null': local_logits.sum()*0}, chains, targets)
