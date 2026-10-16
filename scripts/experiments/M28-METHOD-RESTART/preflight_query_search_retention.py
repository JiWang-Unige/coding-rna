#!/usr/bin/env python3
"""A30 fixed eight-case zero-update actual TRAIN search/gradient preflight."""
import json
import math
import time
import torch
import run_query_known_mean as A
import run_chain_quality as C
from query_search_retention_training import search_retention_loss
from trace_frozen_query_beam import trace_query, original_query, reconcile

ROOT, Q = A.ROOT, A.Q
OUT = ROOT/'outputs/M28-A30-SEARCH-RETENTION-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a30_search_retention_protocol.md'
CONFIG = {**A.CONFIG, 'search_retention': {
    'coefficient': 1., 'margin_log_score_units': 1.,
    'collection': 'current_original_beam_assigned_positive_query_first_loss',
    'negative_qualification': 'complete_R4_known_nonmatch_only',
    'reduction': 'mean_over_all_assigned_positive_queries_including_zero',
    'raw_unsupported': 'original_teacher_supervision_only',
    'beam_width': 4, 'completed_limit_per_query': 4, 'max_exons': 128,
    'minimum_raw_intron_bp': 5, 'reference_injected_candidates': 0},
    'new_parameters': 0, 'initial_context': 'unchanged_A3_global'}


def gradient_norm(values):
    values = [v for v in values if v is not None]
    if any(not bool(torch.isfinite(v).all()) for v in values):
        raise ValueError('Nonfinite actual gradient')
    return math.sqrt(sum(float(v.double().square().sum()) for v in values))


def main():
    for name, job in (('M28-SEARCH-SCORE-REPLAY-UNITCHECK', 13393178),
                      ('M28-SEARCH-RETENTION-UNITCHECK', 13393303)):
        source = ROOT/'outputs'/name
        if ((source/'STATUS').read_text().strip() != 'COMPLETED' or
                int((source/'JOBID').read_text().strip()) != job):
            raise ValueError('Actual completed score/collection tests required')
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population required')
    selected = Q.select_timing_draws(draws, windows, catalogs)
    if len(selected) != 8:
        raise ValueError('Original eight metadata cases required')
    model, original = A.new_model().train(), A.new_model().train()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if sum(p.numel() for p in model.parameters()) != 649457 or any(
            not torch.equal(v, before[k]) for k, v in original.state_dict().items()):
        raise ValueError('Complete same fresh A3 parameters required')
    group_parameters = {name: tuple(getattr(model, name).parameters()) for name in ('shared', 'query')}
    flat = group_parameters['shared']+group_parameters['query']
    records = []
    began = time.perf_counter()
    torch.cuda.reset_peak_memory_stats()
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
        weight = draw['importance_weight_uniform_within_stratum']
        base = A.supervised_readout(model, outputs, w, catalogs[w['species']], cached['sequence'], weight)
        with torch.no_grad():
            prior = A.supervised_readout(original, original(features, dna, [w['valid_bases']]),
                w, catalogs[w['species']], cached['sequence'], weight)
        for key in ('gene_to_query', 'query_labels', 'positive_chains'):
            if base['query'][key] != prior['query'][key]:
                raise ValueError('Fresh A3 original assignment/labels differ')
        if base['protection_reasons'] != prior['protection_reasons']:
            raise ValueError('Original protection differs')
        for key, value in base['losses'].items():
            torch.testing.assert_close(value, prior['losses'][key], rtol=1e-6, atol=1e-6)
        torch.testing.assert_close(base['total'], prior['total'], rtol=1e-6, atol=1e-6)
        targets = Q.build_targets(w, catalogs[w['species']], cached['sequence'])
        window = model.query_window(outputs['features'][0], w['valid_bases'])
        t = time.perf_counter()
        extra = search_retention_loss(window, cached['sequence'], base['query']['positive_chains'],
            targets, CONFIG['search_retention']['margin_log_score_units'])
        torch.cuda.synchronize()
        collection_seconds = time.perf_counter()-t
        if not bool(torch.isfinite(extra['loss'])):
            raise ValueError('Nonfinite current-search margin')
        reconciled = 0
        for row in extra['records']:
            if not row['raw_grammar_supported']:
                continue
            q, chain = row['query_index'], row['chain']
            natural, observed = trace_query(window, cached['sequence'], q, {0: chain})
            reconcile(original_query(window, cached['sequence'], q), natural)
            if observed[0]['returned'] != row['returned'] or observed[0]['first_loss'] != row['first_loss']:
                raise ValueError('Collection changed current original first loss')
            reconciled += 1
        values = torch.autograd.grad(extra['loss'], flat, retain_graph=True, allow_unused=True)
        split = len(group_parameters['shared'])
        extra_gradients = dict(shared=gradient_norm(values[:split]), query=gradient_norm(values[split:]))
        total = base['total']+weight*CONFIG['search_retention']['coefficient']*extra['loss']
        if not bool(torch.isfinite(total)):
            raise ValueError('Nonfinite combined loss')
        total.backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(), float('inf'), error_if_nonfinite=True)
        gradients = {name: gradient_norm(p.grad for p in getattr(model, name).parameters())
                     for name in ('shared', 'segmentation', 'query')}
        torch.cuda.synchronize()
        labels = base['query']['query_labels']
        record = dict(window_id=w['id'], species=w['species'], category=witness['category'],
            draw_ordinal=witness['ordinal'], valid_bases=w['valid_bases'], importance_weight=weight,
            original_four_losses={k: float(v.detach()) for k, v in base['losses'].items()},
            weighted_original_total=float(base['total'].detach()), extra_margin=float(extra['loss'].detach()),
            weighted_combined_total=float(total.detach()), gradient_norm=float(norm),
            combined_gradient_norms=gradients, extra_gradient_norms=extra_gradients,
            query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
            protection_reasons=base['protection_reasons'], positive_queries=extra['positive_queries'],
            collection_records=extra['records'], original_queries_reconciled=reconciled,
            original_four_losses_assignment_and_protection_match=True,
            current_search_collection_and_live_score_seconds=collection_seconds,
            full_two_model_comparison_collection_reconciliation_backward_seconds=time.perf_counter()-started)
        records.append(record)
        print(json.dumps({k: v for k, v in record.items() if k not in ('collection_records', 'protection_reasons')},
                         allow_nan=False), flush=True)
        del cached, features, dna, outputs, base, prior, targets, window, extra, values, total
    if any(not torch.equal(before[k], v) for k, v in model.state_dict().items()):
        raise ValueError('Zero-update preflight changed weights')
    if any(not any(r['combined_gradient_norms'][name] > 0 for r in records)
           for name in ('shared', 'segmentation', 'query')):
        raise ValueError('Original common branch gradient missing')
    if any(not any(r['species'] == species and r['extra_gradient_norms'][name] > 0 for r in records)
           for species in A.ALLOW for name in ('shared', 'query')):
        raise ValueError('Fixed actual cases do not establish both-species shared/query margin gradients')
    C.R.write_summary(OUT, dict(completed=True, arm='A30', configuration=CONFIG, protocol=PROTOCOL,
        original_draws=1536, unique_TRAIN_windows=1491, cases=len(records), records=records,
        trainable_parameters=649457, exact_complete_A3_fresh_parameters=True,
        original_four_losses_assignment_and_protection_verified=True,
        original_search_geometry_order_scores_and_first_loss_verified=True,
        finite_actual_combined_and_extra_gradients=True,
        both_species_actual_shared_and_query_extra_gradients_verified=True,
        model_weights_unchanged=True, optimizer_steps_this_run=0,
        supervised_checkpoint_inheritance=False, reference_injected_candidates=0,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        seconds=time.perf_counter()-began,
        limitation='Zero-update engineering readiness only; no candidate-quality or publication conclusion.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
