"""Candidate coarse-observation channel; preparation, not an approved fit.

Fine labels remain exact. Every observed coarse class contributes a normalized
three-class likelihood: correct 1-epsilon, each other class epsilon/2. Unknown
rows contribute one. Original targets and the complete-gene universe do not
change. Epsilon is supplied by the caller; no fitted value is selected here.
"""
import torch
from .full_grammar_observed_partition import allowed_roles_from_local_targets
from .full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition


def coarse_channel_nll(sequence,logits,max_exons,core_start,core_end,targets,epsilon,
                       *,cuts=(),checkpoint_bp=256):
    n=core_end-core_start; length=len(sequence)
    if not 0<=core_start<core_end<=length:
        raise ValueError('Original core must lie inside the supplied fixed DNA view')
    original=allowed_roles_from_local_targets(targets,n)
    fine=torch.as_tensor(targets['seg'],device='cpu')[:n]
    region=torch.as_tensor(targets['region'],device='cpu')[:n]
    epsilon=torch.as_tensor(epsilon,dtype=torch.float64,device='cpu')
    if epsilon.ndim or not bool(torch.isfinite(epsilon)) or not 0<float(epsilon.detach())<1:
        raise ValueError('One finite scalar observation error probability strictly between zero and one required')
    exact=fine>=0; coarse=(fine<0)&(region>=0)
    allowed=torch.ones((length,15),dtype=torch.bool)
    allowed[core_start:core_end][exact]=original[exact]
    logq=torch.zeros((length,15),dtype=torch.float64)
    positions=torch.where(coarse)[0]+core_start
    logq[positions]=torch.where(original[coarse],torch.log1p(-epsilon),torch.log(epsilon/2))
    denominator=compiled_boundary_full_grammar_log_partition(sequence,logits,max_exons,
        cuts=cuts,checkpoint_bp=checkpoint_bp)
    weighted=compiled_boundary_full_grammar_log_partition(sequence,logits.double()+logq,max_exons,
        allowed=allowed,cuts=cuts,checkpoint_bp=checkpoint_bp)
    # The existing energy is logit(role)-logit(background). Adding logq also
    # subtracts logq(background), so restore that absolute observation factor.
    background_correction=logq[:,0].sum()
    observed=weighted['log_partition']+background_correction
    return dict(loss=denominator['log_partition']-observed,
        log_partition_all=denominator['log_partition'],log_partition_observed=observed,
        weighted_relative_log_partition=weighted['log_partition'],background_correction=background_correction,
        fine_positions=int(exact.sum()),coarse_positions=int(coarse.sum()),
        fully_unknown_positions=n-int(exact.sum())-int(coarse.sum()),
        epsilon=epsilon,statistics_all=denominator['statistics'],statistics_observed=weighted['statistics'],
        original_labels_changed=False,coarse_observation_rows_dropped=False,
        inference_grammar_changed=False,fit_error_probability_selected=False)
