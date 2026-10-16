"""A20 preparation only: marginal likelihood of partially observed chain sets.

Unknown candidates are latent, not assigned label zero. This is NOT the old
zero-gradient unknown policy: span conflicts can give unknown gains gradients.
No experiment runner selects this loss for training.
"""
from src.m28.chain_set import chain_set_loss, interval_log_partition


def span_overlap(a, b):
    return a[0][0] < b[-1][1] and b[0][0] < a[-1][1]


def partial_chain_set_loss(scores, chains, targets):
    prior = chain_set_loss(scores, chains, targets)
    gain = scores['logits']
    positive = [i for i, y in enumerate(targets) if y == 1]
    unknown = [i for i, y in enumerate(targets) if y == -1]
    allowed = [i for i in unknown if not any(span_overlap(chains[i], chains[j]) for j in positive)]
    full = interval_log_partition(gain, chains)
    latent = interval_log_partition(gain[allowed], [chains[i] for i in allowed])
    observed = prior['gold_gain'] + latent
    known = prior['known_candidates']
    # When every variable is unobserved, the observed event is the full space.
    # Reuse the same partition and the original differentiable zero explicitly.
    if not known:
        observed, loss = full, prior['loss']
    else:
        loss = (full-observed)/known
    return dict(loss=loss, log_partition=full, observed_log_partition=observed,
                gold_gain=prior['gold_gain'], known_candidates=known,
                known_positives=len(positive), known_nonmatches=targets.count(0),
                latent_unknown=len(unknown), latent_unknown_compatible_with_gold=len(allowed),
                unknown_conflicting_with_gold=len(unknown)-len(allowed))
