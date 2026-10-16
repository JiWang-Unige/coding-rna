"""A39 fixed equal composite of absolute null and compatible-set constraints."""
from .core import chain_loss
from .chain_set import chain_set_loss


def null_anchored_chain_set_loss(scores, chains, targets):
    structured = chain_set_loss(scores, chains, targets)
    pointwise = chain_loss(scores, targets)
    return {**structured, 'loss': (pointwise+structured['loss'])/2,
            'mean_known_BCE': pointwise, 'mean_known_set_NLL': structured['loss']}
