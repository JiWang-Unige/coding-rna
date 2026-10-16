"""Training-only local and finite-chain supervision; no inference routing."""
import math
from collections import Counter
import torch

from .labels import local_loss, oriented_parts
from .structured_emission import structured_emission_loss


def local_structured_emission_loss(local_logits, chains, chain_targets, local_targets,
                                   valid_bases=None, *, local_coefficient,
                                   structured_coefficient):
    """Keep each existing loss's own mask and denominator on the same logits.

    Local targets must come from the original R4 build_targets, not from chain
    positives or an old cache window. Unknown chains are excluded only from
    the finite partition; they do not manufacture per-base negative labels.
    Coefficients have no defaults: a future fit must declare them explicitly.
    The caller applies the original window importance once to the total.
    """
    weights = (float(local_coefficient), float(structured_coefficient))
    if any(not math.isfinite(w) or w < 0 for w in weights) or not any(weights):
        raise ValueError('Explicit finite nonnegative loss coefficients required')
    structured = structured_emission_loss(local_logits, chains, chain_targets, valid_bases)
    length = len(local_logits) if valid_bases is None else int(valid_bases)
    fine = torch.as_tensor(local_targets['seg'], device=local_logits.device)
    region = torch.as_tensor(local_targets['region'], device=local_logits.device)
    if fine.shape != (len(local_logits),) or region.shape != fine.shape:
        raise ValueError('Local targets must match the full oriented logit window')
    if bool((fine[length:] >= 0).any()) or bool((region[length:] >= 0).any()):
        raise ValueError('Original R4 local targets must ignore padded positions')
    known = (fine >= 0) | ((fine < 0) & (region >= 0) & (region < 3))
    local = local_loss(local_logits, local_targets)
    return dict(structured, loss=weights[0]*local+weights[1]*structured['loss'],
                local_loss=local, structured_loss=structured['loss'],
                local_supervised_positions=int(known.sum()), loss_coefficients=list(weights))


def reconcile_local_supervision(record, window, catalog, support_artifact, local_artifact):
    """Join already validated original R4 inputs to the closed finite support.

    This does not generate targets or classify unknown chains. The original
    loader retains responsibility for TRAIN filtering and transcript identity.
    Geometry, the actual cache source and complete gold must agree before a
    caller combines the two supervision granularities.
    """
    fields = ('id', 'species', 'split', 'seqid', 'strand', 'start', 'end',
              'width', 'valid_bases')
    if (window['split'] != 'train' or window['label_policy'] != 'CDS_assessable_primary_R4'
            or any(record['window'][k] != window[k] for k in fields)):
        raise ValueError('Original R4 local and finite-support window geometry differs')
    if support_artifact != local_artifact:
        raise ValueError('Original R4 local and finite-support cache source differs')
    original = Counter(tuple(oriented_parts(catalog[ident], window))
                       for ident in window['positive_ids'])
    finite = Counter(tuple(tuple(exon) for exon in row['chain'])
                     for row in record['candidates'] if row['label'] == 1)
    if original != finite:
        raise ValueError('Complete original R4 and finite-support gold differs')


def position_normalized_local_structured_loss(local_logits, chains, chain_targets,
        local_targets, valid_bases, *, local_coefficient, structured_coefficient):
    """A44 fit units: original local mean plus finite-set NLL per valid base.

    This is a declared new fit objective, not a change to the closed A43 loss
    or A44's original known-candidate-mean probe. Unknown chains remain absent
    from logZ and gold; local unknown positions keep their original masks.
    """
    result = local_structured_emission_loss(local_logits, chains, chain_targets,
        local_targets, valid_bases, local_coefficient=local_coefficient,
        structured_coefficient=structured_coefficient)
    structured = (result['log_partition']-result['gold_gain'])/int(valid_bases)
    return dict(result, loss=float(local_coefficient)*result['local_loss']
                +float(structured_coefficient)*structured,
                structured_candidate_mean_loss=result['structured_loss'],
                structured_loss=structured, structured_normalizer_valid_bases=int(valid_bases))
