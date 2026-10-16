"""Exact candidate-set likelihood with the existing half-open span constraint."""
from bisect import bisect_right
import torch


def interval_log_partition(gain, chains):
    """Log sum exp of summed gains over all compatible sets, including empty.

    This is the sum-product counterpart of core.select_nonoverlapping, not a
    softmax choosing one gene per window. Candidate geometry is already closed.
    """
    if len(chains) != len(gain):
        raise ValueError('One gain per candidate required')
    order = sorted(range(len(chains)), key=lambda i: (chains[i][-1][1], chains[i][0][0], i))
    ends = [chains[i][-1][1] for i in order]
    states = [gain.sum()*0]
    for k, i in enumerate(order):
        previous = bisect_right(ends, chains[i][0][0], 0, k)
        states.append(torch.logaddexp(states[-1], states[previous]+gain[i]))
    return states[-1]


def chain_set_loss(scores, chains, targets):
    """Known candidates only; every known positive must belong to the gold set.

    Division by the known population matches A6's per-candidate BCE scale.
    With no span conflicts this is exactly that BCE, including multiple genes.
    Unknown/protected candidates never enter the partition or gold constraints.
    """
    gain = scores['logits']
    if len(chains) != len(targets) or len(chains) != len(gain):
        raise ValueError('One original label and gain per candidate required')
    if any(y not in (-1, 0, 1) for y in targets):
        raise ValueError('Expected original 1/0/-1 labels')
    known = [i for i, y in enumerate(targets) if y >= 0]
    positive = [i for i, y in enumerate(targets) if y == 1]
    keys = [tuple(map(tuple, chains[i])) for i in known]
    if len(set(keys)) != len(keys):
        raise ValueError('Known candidate geometry must be unique')
    ordered = sorted(positive, key=lambda i: chains[i][0][0])
    if any(chains[a][-1][1] > chains[b][0][0] for a, b in zip(ordered, ordered[1:])):
        raise ValueError('Gold positives conflict under the actual span selector')
    zero = gain.sum()*0+scores['null']*0
    logz = interval_log_partition(gain[known], [chains[i] for i in known]) if known else zero
    gold = gain[positive].sum() if positive else zero
    return {'loss': (logz-gold)/len(known) if known else zero,
            'log_partition': logz, 'gold_gain': gold,
            'known_candidates': len(known), 'known_positives': len(positive),
            'known_nonmatches': targets.count(0), 'unknown_excluded': targets.count(-1)}
