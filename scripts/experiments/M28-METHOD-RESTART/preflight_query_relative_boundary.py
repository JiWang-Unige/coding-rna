#!/usr/bin/env python3
"""A23 fixed zero-update real TRAIN preflight, not a candidate-quality result."""
import math
import time
import torch
import run_query_known_mean as A
import run_chain_quality as C
from src.m28.query_relative_boundary import RelativeBoundaryGeneQueryModel

ROOT, Q = A.ROOT, A.Q
OUT = ROOT/'outputs/M28-A23-RELATIVE-BOUNDARY-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a23_relative_boundary_protocol.md'
CONFIG = {**A.CONFIG, 'relative_boundary_prior': 'distance_over_valid_length_and_log_distance_2_32_2',
    'relative_boundary_events': ['end', 'next_start'], 'relative_boundary_output_bias': False,
    'relative_boundary_last_weight_initialization': 'zero', 'new_parameters': 160,
    'initial_context': 'unchanged_A3_global', 'pointer_domain': 'unchanged_global_geometry'}


def new_model(device='cuda'):
    torch.manual_seed(CONFIG['seed'])
    return RelativeBoundaryGeneQueryModel(CONFIG['feature_dim'], CONFIG['num_queries'], CONFIG['hidden']).to(device)


def main():
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population required')
    selected = Q.select_timing_draws(draws, windows, catalogs)
    model, original = new_model().train(), A.new_model().train()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if sum(p.numel() for p in model.parameters()) != 649617 or any(
        not torch.equal(v, before[k]) for k, v in original.state_dict().items()):
        raise ValueError('Same complete fresh A3 common parameters plus160 new parameters required')
    records = []
    began = time.perf_counter()
    for witness in selected:
        draw = witness['draw']
        w = windows[draw['window_id']]
        cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        features, dna = Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Cached backbone and DNA must remain fixed')
        model.zero_grad(set_to_none=True)
        started = time.perf_counter()
        outputs = model(features, dna, [w['valid_bases']])
        result = A.supervised_readout(model, outputs, w, catalogs[w['species']], cached['sequence'],
                                      draw['importance_weight_uniform_within_stratum'])
        with torch.no_grad():
            prior = A.supervised_readout(original, original(features, dna, [w['valid_bases']]), w,
                catalogs[w['species']], cached['sequence'], draw['importance_weight_uniform_within_stratum'])
        if result['query']['gene_to_query'] != prior['query']['gene_to_query']:
            raise ValueError('Zero-delta fresh assignment differs from A3')
        for key, value in result['losses'].items():
            torch.testing.assert_close(value, prior['losses'][key], rtol=1e-6, atol=1e-6)
        result['total'].backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(), float('inf'), error_if_nonfinite=True)
        gradients = {}
        for name in ('shared', 'segmentation', 'query'):
            values = [p.grad for p in getattr(model, name).parameters() if p.grad is not None]
            if any(not bool(torch.isfinite(v).all()) for v in values):
                raise ValueError('Nonfinite '+name+' gradient')
            gradients[name] = math.sqrt(sum(float(v.double().square().sum()) for v in values))
        last = model.query.relative_boundary[-1].weight.grad
        event_gradients = [0., 0.] if last is None else [float(v.double().norm()) for v in last]
        assigned = result['query']['positive_chains']
        applicable = [bool(assigned), any(len(c) > 1 for c in assigned.values())]
        if any((g <= 0 if use else g != 0) for g, use in zip(event_gradients, applicable)):
            raise ValueError('Relative-event output gradient differs from actual supervision support')
        first = [p.grad for p in model.query.relative_boundary[0].parameters()]
        if any(g is not None and bool((g != 0).any()) for g in first):
            raise ValueError('Zero last weight should give zero first-layer gradient at initialization')
        torch.cuda.synchronize()
        labels = result['query']['query_labels']
        records.append(dict(window_id=w['id'], species=w['species'], category=witness['category'],
            draw_ordinal=witness['ordinal'], valid_bases=w['valid_bases'],
            weighted_total=float(result['total'].detach()), gradient_norm=float(norm), gradient_norms=gradients,
            relative_event_applicable=applicable, relative_event_output_gradient_norms=event_gradients,
            query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
            protection_reasons=result['protection_reasons'], zero_delta_A3_losses_and_assignment_match=True,
            forward_two_models_supervision_backward_seconds=time.perf_counter()-started))
        del cached, features, dna, outputs, result, prior
    if any(not torch.equal(before[k], v) for k, v in model.state_dict().items()):
        raise ValueError('Preflight changed model weights')
    if any(not any(r['gradient_norms'][k] > 0 for r in records) for k in ('shared', 'segmentation', 'query')):
        raise ValueError('A required common branch has no actual gradient')
    if any(not any(r['species'] == s and r['relative_event_output_gradient_norms'][i] > 0 for r in records)
           for s in A.ALLOW for i in (0, 1)):
        raise ValueError('Both species must cover both relative events in the fixed metadata cases')
    C.R.write_summary(OUT, dict(completed=True, arm='A23', configuration=CONFIG, protocol=PROTOCOL,
        original_draws=1536, unique_TRAIN_windows=1491, cases=len(records), records=records,
        trainable_parameters=649617, exact_A3_common_fresh_parameters=True,
        zero_delta_actual_losses_and_assignment_verified=True, model_weights_unchanged=True,
        optimizer_steps_this_run=0, supervised_checkpoint_inheritance=False,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began, limitation='Zero-update engineering readiness only; no candidate-quality or publication conclusion.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
