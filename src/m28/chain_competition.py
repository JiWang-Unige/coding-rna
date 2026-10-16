"""TRAIN-only local chain-vs-null/conflict objective; not a window-wide softmax."""
import torch
from torch.nn import functional as F


def chain_competition_loss(scores, chains, targets):
    gain = scores['logits']
    if len(chains) != len(targets) or len(chains) != len(gain):
        raise ValueError('One label and score per original candidate required')
    positives = [i for i, y in enumerate(targets) if y == 1]
    negatives = [i for i, y in enumerate(targets) if y == 0]
    terms, positive_losses = [], []
    conflict_pairs = 0
    for p in positives:
        a, b = chains[p][0][0], chains[p][-1][1]
        rivals = [n for n in negatives if chains[n][0][0] < b and chains[n][-1][1] > a]
        # Each positive has its own alternatives. Other correct genes are never
        # mutually exclusive classes; protected/unknown chains are never rivals.
        alternatives = torch.cat([gain.new_zeros(1), gain[p:p+1], gain[rivals]])
        positive_losses.append(torch.logsumexp(alternatives, dim=0)-gain[p])
        conflict_pairs += len(rivals)
    positive = torch.stack(positive_losses).mean() if positive_losses else None
    negative = F.softplus(gain[negatives]).mean() if negatives else None
    if positive is not None:
        terms.append(positive)
    if negative is not None:
        terms.append(negative)
    zero = gain.sum()*0+scores['null']*0
    return {'loss': torch.stack(terms).mean() if terms else zero,
            'positive_competition': positive if positive is not None else zero,
            'negative_null_anchor': negative if negative is not None else zero,
            'known_positives': len(positives), 'known_nonmatches': len(negatives),
            'conflict_pairs': conflict_pairs, 'unknown_excluded': targets.count(-1)}
