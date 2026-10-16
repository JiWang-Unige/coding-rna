#!/usr/bin/env python3
"""A25 fixed zero-update real TRAIN preflight; no candidate-quality claim."""
import math
import time
import torch
import run_query_known_mean as A
import run_query_local_context as L
import run_chain_quality as C
from src.m28.labels import build_targets
from src.m28.query_start_weighted import StartWeightedGeneQueryModel, joint_start_weighted_losses

ROOT, Q = A.ROOT, A.Q
OUT = ROOT/'outputs/M28-A25-START-WEIGHTED-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a25_start_weighted_protocol.md'
CONFIG = {**A.CONFIG, 'initial_context': 'A17_local_attention_plus_shared_initiation_logits',
    'pointer_domain': 'unchanged_global_geometry', 'initiation_head': 'linear_hidden_to_one_no_bias',
    'initiation_weight_initialization': 'zero', 'initiation_ce_coefficient': 1.,
    'initiation_ce_support': 'known_fine_non13_plus_eligible_positive_starts', 'new_parameters': 128}


def new_model(device='cuda'):
    torch.manual_seed(CONFIG['seed'])
    return StartWeightedGeneQueryModel(CONFIG['feature_dim'], CONFIG['num_queries'], CONFIG['hidden']).to(device)


def supervised_readout(model, outputs, window, catalog, sequence, weight):
    targets = build_targets(window, catalog, sequence)
    result = joint_start_weighted_losses(model, outputs, targets, window['valid_bases'], weight,
                                        CONFIG['loss_coefficients'])
    if not bool(torch.isfinite(result['total'])):
        raise ValueError('Nonfinite A25 training loss')
    return result


def main():
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population required')
    selected = Q.select_timing_draws(draws, windows, catalogs)
    model, original = new_model().train(), L.new_model().train()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if sum(p.numel() for p in model.parameters()) != 649585 or any(
        not torch.equal(v, before[k]) for k, v in original.state_dict().items()):
        raise ValueError('Fresh complete A17 common parameters plus128 required')
    records = []
    began = time.perf_counter()
    for witness in selected:
        draw = witness['draw']; w = windows[draw['window_id']]
        cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        features, dna = Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Cached backbone and DNA must remain fixed')
        model.zero_grad(set_to_none=True)
        started = time.perf_counter()
        outputs = model(features, dna, [w['valid_bases']])
        weight = draw['importance_weight_uniform_within_stratum']
        result = supervised_readout(model, outputs, w, catalogs[w['species']], cached['sequence'], weight)
        with torch.no_grad():
            prior = A.supervised_readout(original, original(features, dna, [w['valid_bases']]), w,
                                        catalogs[w['species']], cached['sequence'], weight)
        if result['query']['gene_to_query'] != prior['query']['gene_to_query']:
            raise ValueError('Zero-prior fresh assignment differs from A17')
        for key, value in prior['losses'].items():
            torch.testing.assert_close(result['losses'][key], value, rtol=1e-6, atol=1e-6)
        result['total'].backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(), float('inf'), error_if_nonfinite=True)
        gradients = {}
        for name in ('shared', 'segmentation', 'query'):
            values = [p.grad for p in getattr(model, name).parameters() if p.grad is not None]
            if any(not bool(torch.isfinite(v).all()) for v in values):
                raise ValueError('Nonfinite '+name+' gradient')
            gradients[name] = math.sqrt(sum(float(v.double().square().sum()) for v in values))
        gradient = model.query.initiation.weight.grad
        initiation_norm = 0. if gradient is None else float(gradient.double().norm())
        positive = bool(result['initiation_starts'])
        if positive and (initiation_norm <= 0 or not math.isfinite(float(result['losses']['initiation_ce']))):
            raise ValueError('Actual positive initiation supervision has no finite active gradient')
        if not positive and float(result['losses']['initiation_ce']) != 0:
            raise ValueError('No-positive window received auxiliary absence supervision')
        torch.cuda.synchronize()
        labels = result['query']['query_labels']
        records.append(dict(window_id=w['id'], species=w['species'], category=witness['category'],
            draw_ordinal=witness['ordinal'], valid_bases=w['valid_bases'],
            weighted_total=float(result['total'].detach()), gradient_norm=float(norm), gradient_norms=gradients,
            initiation_starts=result['initiation_starts'], initiation_support_bases=result['initiation_support_bases'],
            initiation_ce=float(result['losses']['initiation_ce'].detach()), initiation_gradient_norm=initiation_norm,
            query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
            protection_reasons=result['protection_reasons'], zero_prior_A17_losses_and_assignment_match=True,
            forward_two_models_supervision_backward_seconds=time.perf_counter()-started))
        del cached, features, dna, outputs, result, prior
    if any(not torch.equal(before[k], v) for k, v in model.state_dict().items()):
        raise ValueError('Preflight changed model weights')
    if any(not any(r['gradient_norms'][k] > 0 for r in records) for k in ('shared', 'segmentation', 'query')):
        raise ValueError('Required common branch has no actual gradient')
    if any(not any(r['species'] == s and r['initiation_starts'] and r['initiation_gradient_norm'] > 0
                   for r in records) for s in A.ALLOW):
        raise ValueError('Both species need actual eligible initiation supervision')
    C.R.write_summary(OUT, dict(completed=True, arm='A25', configuration=CONFIG, protocol=PROTOCOL,
        original_draws=1536, unique_TRAIN_windows=1491, cases=len(records), records=records,
        trainable_parameters=649585, exact_A17_common_fresh_parameters=True,
        zero_prior_actual_old_losses_and_assignment_verified=True, model_weights_unchanged=True,
        optimizer_steps_this_run=0, supervised_checkpoint_inheritance=False,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began,
        limitation='Zero-update engineering readiness only; no candidate-quality or publication conclusion.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
