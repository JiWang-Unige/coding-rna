"""Independent per-base likelihood with the existing fine/coarse channel."""
import torch
from .full_grammar_observed_partition import allowed_roles_from_local_targets


def local_coarse_channel_nll(logits, targets, epsilon=0.1):
    n = len(logits)
    allowed = allowed_roles_from_local_targets(targets, n)
    if logits.shape != (n, 15):
        raise ValueError('One fifteen-role row per original core base required')
    epsilon = torch.as_tensor(epsilon, dtype=torch.float64, device=logits.device)
    if epsilon.ndim or not bool(torch.isfinite(epsilon)) or not 0 < float(epsilon.detach()) < 1:
        raise ValueError('One finite observation error probability between zero and one required')
    fine = torch.as_tensor(targets['seg'], device=logits.device)[:n]
    region = torch.as_tensor(targets['region'], device=logits.device)[:n]
    exact = fine >= 0
    coarse = (fine < 0) & (region >= 0)
    logp = logits.double().log_softmax(dim=1)
    # The zero term retains a zero VJP when the complete core is unknown.
    loss = logp.sum() * 0
    loss = loss - logp[torch.where(exact)[0], fine[exact].long()].sum()
    logq = torch.where(allowed.to(logits.device)[coarse], torch.log1p(-epsilon), torch.log(epsilon/2))
    loss = loss - torch.logsumexp(logp[coarse] + logq, dim=1).sum()
    return dict(loss=loss, fine_positions=int(exact.sum()), coarse_positions=int(coarse.sum()),
                fully_unknown_positions=n-int(exact.sum())-int(coarse.sum()))
