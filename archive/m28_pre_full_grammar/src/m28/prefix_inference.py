"""B2 paired final readout and frozen additive-control reconciliation."""
import numpy as np
import torch

def chain_key(chain):
    return tuple(tuple(exon) for exon in chain)

@torch.no_grad()
def score_paired_pools(model,features,early,late):
    """One fresh sorted union readout; common identities receive identical values."""
    union=sorted({chain_key(c) for pool in (early,late) for c in pool})
    scorer=model.for_window(features)
    logits=scorer.score('final',None,[(0.,c,'') for c in union])
    if not bool(torch.isfinite(logits).all()):
        raise ValueError('Nonfinite paired final gain')
    values=dict(zip(union,logits.detach().cpu().tolist()))
    return {name:[values[chain_key(c)] for c in pool]
            for name,pool in (('early',early),('late',late))}

def validate_additive_replay(actual,saved,tolerance=1e-5):
    if actual['reference_used'] is not False or saved['reference_used'] is not False:
        raise ValueError('Reference entered additive inference control')
    if ([chain_key(c) for c in actual['chains']]!=[chain_key(c) for c in saved['chains']] or
        [tuple(p) for p in actual['links']]!=[tuple(p) for p in saved['links']] or
        actual['counts']!=saved['counts'] or actual['budget']!=saved['budget']):
        raise ValueError('R5 additive candidate/control replay mismatch')
    a=np.asarray(actual['proposal_scores'],dtype=float)
    b=np.asarray(saved['proposal_scores'],dtype=float)
    if a.shape!=b.shape or a.shape!=(len(actual['chains']),) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('R5 proposal alignment or finite-score mismatch')
    delta=float(np.abs(a-b).max()) if len(a) else 0.
    if delta>tolerance:
        raise ValueError('R5 additive proposal error exceeds frozen tolerance: '+str(delta))
    return delta
