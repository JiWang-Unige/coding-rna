#!/usr/bin/env python3
"""A43 closed finite support through fresh emissions, no optimizer or checkpoint."""
from collections import Counter
import json
import math
import time
import torch

import build_natural_competitor_support as B
import preflight_structured_emission as A36
from src.m28.chain_set import chain_set_loss
from src.m28.structured_emission import structured_emission_loss, trainable_path_scores

ROOT, PROTOCOL, C = B.ROOT, B.PROTOCOL, A36.C
OUT = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-PREFLIGHT'
CONFIG = {**A36.CONFIG,
    'candidate_pool': 'closed_A43_19189_original_supervision_plus_all_A42_natural_geometry',
    'positive_targets': 'unchanged_original_A12_2916_positive_window_chain_pairs',
    'unknown_candidates': 1844, 'known_nonmatch_candidates': 14429,
    'source_support_configuration': B.CONFIG,
    'comparison': 'frozen_A31_A3_original1491_TRAIN_full_grammar_natural_readout_not_finite_pool_oracle',
    'support_scope': 'offline_frozen_A3_natural_competitors_not_online_mining_or_full_grammar_partition'}


def population(support_job):
    B.D.B.F.A.require_jobs([support_job])
    summary = B.D.R.completed(B.OUT, support_job)
    proofs = ('actual_complete_original_supervision_and_natural_additions_reconciled',
              'actual_all_A42_correct_geometry_in_original_gold_and_no_positive_added',
              'actual_original_gold_compatibility_and_reference_denominators_preserved',
              'actual_original_draw_order_species_and_importance_weights_reconciled',
              'all_support_files_closed_before_summary')
    if (summary['protocol'] != PROTOCOL or summary['configuration'] != B.CONFIG or
            any(not summary[k] for k in proofs) or summary['optimizer_steps_this_run'] != 0 or
            summary['new_geometry_labels_computed'] or summary['DEV_used']):
        raise ValueError('Closed matching original TRAIN union support required')
    draws = C.R.rows(B.OUT/'draws.jsonl'); rows = C.R.rows(B.OUT/'windows.jsonl')
    if (draws != C.R.rows(C.R.Q.DRAWS) or len(draws) != 1536 or len(rows) != 1491 or
            [r['window']['id'] for r in rows] != list(dict.fromkeys(d['window_id'] for d in draws))):
        raise ValueError('Complete unchanged original TRAIN windows and draws required')
    labels = Counter(); origins = Counter(); replay = Counter()
    records = {r['window']['id']: r for r in rows}
    C.R.validate_feature_summary(json.loads((C.R.FEATURES/'summary.json').read_text()))
    entries = C.R.rows(C.R.PLAN/'feature_plan.jsonl')
    index = {e['window']['id']: e for e in entries if e['window']['id'] in records}
    for ordinal, row in enumerate(rows):
        w = row['window']; ident = w['id']; local = row['candidates']
        if (row['original_ordinal'] != ordinal or w['split'] != 'train' or
                w['seqid'] not in C.R.ALLOW[w['species']] or index[ident]['window'] != w):
            raise ValueError('Original allowed feature-plan geometry differs from A43 support')
        labels.update(str(r['label']) for r in local); origins.update(r['origin'] for r in local)
        if any(row['counts']['label_'+str(y)] != sum(r['label'] == y for r in local) for y in (1, 0, -1)):
            raise ValueError('Actual support labels and per-window counts differ')
    for draw in draws:
        if draw['species'] != records[draw['window_id']]['window']['species']:
            raise ValueError('Original TRAIN draw-window species differs')
        replay.update(records[draw['window_id']]['counts'])
    if (dict(labels) != summary['label_counts_pooled'] or sum(labels.values()) != 19189 or
            origins != {'free': 16321, 'reference_injected': 536, 'A42_natural': 2332} or
            dict(replay) != summary['merged_replayed_draw_counts'] or len(index) != 1491):
        raise ValueError('Complete A43 supervision, source origins or replay counts differ')
    return draws, index, records, summary


def witnesses(records):
    selected = []
    for species in C.R.ALLOW:
        local = [(ident, r) for ident, r in records.items() if r['window']['species'] == species
                 and r['counts']['label_1'] and r['counts']['label_0']]
        used = set()
        for category in ('most_added_known', 'most_added_unknown', 'most_positive_exons', 'most_span_conflicts'):
            def priority(pair):
                rows = pair[1]['candidates']
                if category.startswith('most_added_'):
                    label = 0 if category == 'most_added_known' else -1
                    return sum(r['origin'] == 'A42_natural' and r['label'] == label for r in rows)
                if category == 'most_positive_exons':
                    return sum(len(r['chain']) for r in rows if r['label'] == 1)
                return sum(n['chain'][0][0] < p['chain'][-1][1] and n['chain'][-1][1] > p['chain'][0][0]
                           for p in rows if p['label'] == 1 for n in rows if n['label'] == 0)
            eligible = [pair for pair in local if pair[0] not in used and priority(pair) > 0]
            if not eligible:
                raise ValueError('Missing actual union gradient witness: '+species+' '+category)
            chosen = min(eligible, key=lambda pair: (-priority(pair), pair[1]['original_ordinal']))
            used.add(chosen[0])
            selected.append(dict(window_id=chosen[0], species=species, category=category,
                                 priority=priority(chosen)))
    return selected


def main(support_job):
    draws, index, records, support = population(support_job)
    chosen = witnesses(records)
    model = A36.new_model('cuda').train()
    if sum(p.numel() for p in model.parameters()) != CONFIG['trainable_parameters']:
        raise ValueError('Fresh representation-only parameter count differs')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    weight = {d['window_id']: d['importance_weight_uniform_within_stratum'] for d in draws}
    receipts = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    for witness in chosen:
        ident = witness['window_id']; rec = records[ident]; w = rec['window']; local = rec['candidates']
        started = time.perf_counter()
        cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        features, dna = C.R.Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Detached original cached encoder features and DNA required')
        model.zero_grad(set_to_none=True)
        logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0]
        if logits.dtype != torch.float32:
            raise ValueError('FP32 fresh network required')
        chains = [r['chain'] for r in local]; labels = [r['label'] for r in local]
        gain = trainable_path_scores(logits, chains, w['valid_bases'])
        result = structured_emission_loss(logits, chains, labels, w['valid_bases'])
        weighted = result['loss']*weight[ident]
        if not bool(torch.isfinite(weighted)):
            raise ValueError('Nonfinite actual union TRAIN structured loss')
        # Probe actual joined geometry at the same network potentials. This
        # verifies supervision of new known chains, not prediction filtering.
        probe = gain.detach().requires_grad_(True)
        probe_result = chain_set_loss({'logits': probe, 'null': probe.sum()*0}, chains, labels)
        torch.testing.assert_close(probe_result['loss'], result['loss'], rtol=0, atol=0)
        derivative, = torch.autograd.grad(probe_result['loss']*weight[ident], probe)
        added = [i for i, r in enumerate(local) if r['origin'] == 'A42_natural' and r['label'] == 0]
        protected = [i for i, y in enumerate(labels) if y == -1]
        if (not bool(torch.isfinite(derivative).all()) or
                bool((derivative[added] < 0).any()) or bool((derivative[protected] != 0).any())):
            raise ValueError('Actual new-known suppression or unknown exclusion gradient differs')
        weighted.backward()
        groups = {}
        for name, module in (('project', model.shared.project), ('context', model.shared.context),
            ('dna', model.shared.dna), ('fuse', model.shared.fuse), ('segmentation', model.segmentation)):
            params = list(module.parameters())
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite actual union gradient: '+name)
            groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if groups[name] <= 0:
                raise ValueError('Zero actual union gradient: '+name)
        torch.cuda.synchronize()
        receipts.append(dict(**witness, weighted_loss=float(weighted.detach()), original_importance=weight[ident],
            gradient_norms=groups, actual_added_known=len(added), actual_unknown_excluded=len(protected),
            added_known_positive_gradient_count=int((derivative[added] > 0).sum()),
            maximum_added_known_gain_gradient=float(derivative[added].max()) if added else 0.,
            maximum_unknown_gain_gradient=float(derivative[protected].abs().max()) if protected else 0.,
            load_forward_backward_seconds=time.perf_counter()-started,
            **{k: v for k, v in result.items() if not torch.is_tensor(v)}))
        print(receipts[-1], flush=True)
        del cached, features, dna, logits, gain, result, weighted, probe, probe_result, derivative
    if any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()):
        raise ValueError('Zero-update A43 preflight changed model state')
    C.R.write_summary(OUT, dict(completed=True, experiment='M28-A43-NATURAL-COMPETITOR-PREFLIGHT',
        protocol=PROTOCOL, configuration=CONFIG, source_support_job=support_job,
        windows=1491, original_draws=1536, candidate_window_pairs=19189,
        replayed_original_draw_counts=support['merged_replayed_draw_counts'],
        actual_closed_support_and_original_feature_geometry_reconciled=True,
        actual_new_known_gradient_nonnegative_and_unknown_gradient_zero=True,
        new_model_weights_unchanged=True, fresh_seed0_initialization=True,
        supervised_checkpoint_read=False, generator_forward=False,
        original_TRAIN_catalog_metadata_read=False, original_reference_annotation_file_opened=False,
        new_R4_labelling=False, optimizer_steps_this_run=0, checkpoints_written=False,
        new_candidate_search=False, new_backbone_forward=False,
        DEV_used=False, test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        records=receipts, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Actual finite-support gradient readiness only; no trained new model, full-grammar '
            'training, natural prediction improvement, novelty claim or independent generalization.'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-job', required=True)
    torch.set_num_threads(2)
    main(parser.parse_args().support_job)
