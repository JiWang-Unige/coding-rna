#!/usr/bin/env python3
"""A43 finite TRAIN supervision union from closed saved records, zero updates."""
from collections import Counter
from itertools import islice, zip_longest
import json
import math

import census_natural_competitor_support as S
from src.m28.natural_competitor_pool import augment_training_pool

A42, D = S.A42, S.A42.D
ROOT = D.ROOT
PROTOCOL = S.PROTOCOL
OUT = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-SUPPORT'
A12 = ROOT/'outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC'
A36_FIT = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-FIT'
A36_CPU = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-RESULT-CHECK'
CHECKS = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-SUPPORT-CHECKS'
CONFIG = dict(population='original_1491_TRAIN_windows_1536_unchanged_draws',
    original_supervision='closed_A12_16321_natural_plus_536_positive_injected',
    additions='all_closed_A42_natural_selected_geometry_without_label_or_score_selection',
    identity='complete_exon_geometry_within_original_oriented_window',
    gold='original_A12_positive_geometry_unchanged', unknown='retain_minus1_never_negative',
    new_R4_labels=False, optimizer_steps=0, expected_pairs=19189)


def build_records(draws, training_draws, originals, previous, natural, labelled):
    if len(draws) != len(training_draws):
        raise ValueError('Complete original first-pass draw receipts required')
    for ordinal, (draw, receipt) in enumerate(zip(draws, training_draws)):
        weight = draw['importance_weight_uniform_within_stratum']
        if (receipt['step'] != ordinal+1 or receipt['epoch'] != 1 or receipt['draw_ordinal'] != ordinal or
                receipt['window_id'] != draw['window_id'] or receipt['species'] != draw['species'] or
                receipt['importance_weight'] != weight or not math.isfinite(weight) or weight <= 0):
            raise ValueError('Original draw order, species or importance weight differs')
    ids = list(dict.fromkeys(draw['window_id'] for draw in draws))
    per_species = {s: Counter() for s in D.R.C.R.ALLOW}
    original_per_species = {s: Counter() for s in per_species}
    added_per_species = {s: Counter() for s in per_species}
    injected_positive_matches = Counter({s: 0 for s in per_species})
    original_counts_by_id = {}; merged_counts_by_id = {}; species_by_id = {}; records = []
    for ordinal, (old, prior, new, lab) in enumerate(zip_longest(originals, previous, natural, labelled)):
        if any(value is None for value in (old, prior, new, lab)) or ordinal >= len(ids):
            raise ValueError('Complete original four-way TRAIN population required')
        w = new['window']; s = w['species']; ident = w['id']
        D.validate_target_geometry(old['window'], w)
        if (w != prior['window'] or w != lab['window'] or ident != ids[ordinal] or
                w['split'] != 'train' or w['seqid'] not in D.R.C.R.ALLOW[s] or
                new['original_ordinal'] != ordinal or new['reference_used_in_scoring'] is not False or
                prior['reference_used_in_scoring'] is not False or
                new['selected'] != list(range(len(new['chains'])))):
            raise ValueError('Original allowed ordered reference-free TRAIN window differs')
        local = old['candidates']; old_chains = [r['chain'] for r in local]; labels = [r['label'] for r in local]
        free = prior['variants']['A31_A3']['chains']; nfree = len(free)
        if (old_chains[:nfree] != free or labels[:nfree] != lab['labels']['A31_A3'] or
                any(r['window_id'] != ident or r['species'] != s for r in local) or
                any(r['origin'] != 'free' for r in local[:nfree]) or
                any(r['origin'] != 'reference_injected' or r['label'] != 1 for r in local[nfree:])):
            raise ValueError('Original natural prefix, saved labels or positive injection differs')
        original_counts = dict(windows=1, free_candidates=nfree,
            reference_injected_candidates=len(local)-nfree, zero_known_windows=int(not any(y >= 0 for y in labels)),
            **{'label_'+str(y): labels.count(y) for y in (1, 0, -1)})
        if original_counts != old['counts']:
            raise ValueError('Original A12 per-window counts differ')
        merged = augment_training_pool(old_chains, labels, new['chains'], lab['labels']['A42'])
        gold = [c for c, y in zip(merged['chains'], merged['labels']) if y == 1]
        ordered_gold = sorted(gold, key=lambda c: c[0][0])
        if (len(gold) != lab['eligible_references'] or
                sum(len(c) >= 5 for c in gold) != lab['eligible_five_or_more'] or
                any(a[-1][1] > b[0][0] for a, b in zip(ordered_gold, ordered_gold[1:]))):
            raise ValueError('Original complete compatible gold and long-chain denominators differ')
        old_origin = {A42.key(r['chain']): r['origin'] for r in local}
        for chain, label in zip(new['chains'], lab['labels']['A42']):
            if label == 1 and old_origin[A42.key(chain)] == 'reference_injected':
                injected_positive_matches[s] += 1
        candidates = [dict(chain=c, label=y, origin=r['origin'])
                      for c, y, r in zip(merged['chains'], merged['labels'], local)]
        for i, c, y in zip(merged['added_natural_indices'], merged['chains'][len(local):], merged['labels'][len(local):]):
            candidates.append(dict(chain=c, label=y, origin='A42_natural', natural_index=i))
            added_per_species[s][str(y)] += 1
        counts = dict(windows=1, original_free_candidates=nfree,
            reference_injected_candidates=len(local)-nfree,
            added_A42_natural_candidates=len(merged['added_natural_indices']),
            zero_known_windows=int(not any(y >= 0 for y in merged['labels'])),
            **{'label_'+str(y): merged['labels'].count(y) for y in (1, 0, -1)})
        records.append(dict(window=w, original_ordinal=ordinal, candidates=candidates, counts=counts))
        original_counts_by_id[ident] = original_counts; merged_counts_by_id[ident] = counts
        species_by_id[ident] = s
        original_per_species[s].update(str(y) for y in labels)
        per_species[s].update(str(y) for y in merged['labels'])
    if len(records) != len(ids):
        raise ValueError('Every original unique TRAIN window required')
    original_replay = Counter(); merged_replay = Counter()
    for draw in draws:
        ident = draw['window_id']
        if draw['species'] != species_by_id[ident]:
            raise ValueError('Original draw-window species differs')
        original_replay.update(original_counts_by_id[ident]); merged_replay.update(merged_counts_by_id[ident])
    pooled = Counter(); original_pooled = Counter(); added_pooled = Counter()
    for s in per_species:
        pooled.update(per_species[s]); original_pooled.update(original_per_species[s]); added_pooled.update(added_per_species[s])
    stats = dict(windows=len(records), original_draws=len(draws), candidate_window_pairs=sum(pooled.values()),
        original_candidate_window_pairs=sum(original_pooled.values()),
        original_label_counts_pooled=dict(original_pooled), label_counts_pooled=dict(pooled),
        label_counts_per_species={s: dict(v) for s, v in per_species.items()},
        added_label_counts_pooled=dict(added_pooled),
        added_label_counts_per_species={s: dict(v) for s, v in added_per_species.items()},
        A42_correct_matches_original_reference_injected_per_species=dict(injected_positive_matches),
        original_replayed_draw_counts=dict(original_replay), merged_replayed_draw_counts=dict(merged_replay))
    return records, stats


def main(checks_job):
    jobs = ['13360650', '13408451', '13423961', '13423970', A42.JOBS['readout'], A42.JOBS['CPU'], checks_job]
    D.B.F.A.require_jobs(jobs)
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or
            (CHECKS/'JOBID').read_text().strip() != checks_job):
        raise ValueError('Actual new support adapter checks required')
    report = D.R.completed(A42.READOUT, A42.JOBS['readout']); checked = D.R.completed(A42.CHECK, A42.JOBS['CPU'])
    A42.validate_closed(report, checked)
    source = D.R.completed(A12, '13360650'); old_cpu = D.R.completed(A36_CPU, '13408451')
    membership = D.R.completed(S.OUT, '13423970')
    if (source['source_generator'] != D.R.C.GENERATOR or
            source['source_head'] != 'outputs/M28-A6-CHAIN-QUALITY-FIT/step_004608.pt' or
            source['windows'] != 1491 or source['original_draws'] != 1536 or
            source['unique_candidate_window_pairs'] != 16857 or source['optimizer_steps_this_run'] != 0 or
            source['DEV_forward_used'] or source['test_setaria_used'] or
            not old_cpu['actual_original_training_draw_weights_and_4608_updates_verified'] or
            old_cpu['original_replayed_per_pass_counts'] != source['replayed_original_draw_counts'] or
            membership['protocol'] != PROTOCOL or membership['source_A42_configuration'] != D.CONFIG or
            not membership['actual_original_ordered_support_and_saved_labels_reconciled'] or
            not membership['actual_A42_selected_labels_and_original_denominators_reconciled']):
        raise ValueError('Closed original TRAIN supervision and actual support census required')
    draws = D.R.C.R.rows(D.R.Q.DRAWS)
    with (A36_FIT/'training.jsonl').open() as handle:
        receipts = [json.loads(line) for line in islice(handle, len(draws))]
    records, stats = build_records(draws, receipts, D.R.C.R.rows(A12/'windows.jsonl'),
        D.R.C.R.rows(D.R.READOUT/'natural_scores.jsonl'), D.R.C.R.rows(A42.READOUT/'natural_scores.jsonl'),
        D.R.C.R.rows(A42.READOUT/'windows.jsonl'))
    if (stats['windows'] != 1491 or stats['original_draws'] != 1536 or
            stats['original_candidate_window_pairs'] != 16857 or stats['candidate_window_pairs'] != 19189 or
            stats['original_label_counts_pooled'] != {'1': 2916, '0': 12566, '-1': 1375} or
            stats['label_counts_pooled'] != {'1': 2916, '0': 14429, '-1': 1844} or
            stats['original_replayed_draw_counts'] != source['replayed_original_draw_counts']):
        raise ValueError('Actual unchanged original gold/draws and complete union support differ')
    for s in stats['label_counts_per_species']:
        expected_added = {str(y): membership['per_species'][s][name]['new']
                          for y, name in ((0, 'known_nonmatch'), (-1, 'unknown'))}
        if (stats['added_label_counts_per_species'][s] != expected_added or
                stats['A42_correct_matches_original_reference_injected_per_species'][s] !=
                membership['per_species'][s]['correct']['new']):
            raise ValueError('Actual source-pool novelty and original injected gold do not reconcile')
    for filename, rows in (('windows.jsonl', records), ('draws.jsonl', draws)):
        with (OUT/filename).open('x') as handle:
            for row in rows: handle.write(json.dumps(row, allow_nan=False)+'\n')
    D.R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A43-NATURAL-COMPETITOR-SUPPORT',
        configuration=CONFIG, protocol=PROTOCOL, source_checks_job=checks_job, **stats,
        actual_complete_original_supervision_and_natural_additions_reconciled=True,
        actual_all_A42_correct_geometry_in_original_gold_and_no_positive_added=True,
        actual_original_gold_compatibility_and_reference_denominators_preserved=True,
        actual_original_draw_order_species_and_importance_weights_reconciled=True,
        all_support_files_closed_before_summary=True,
        saved_reference_injected_TRAIN_geometry_opened=True, TRAIN_catalog_opened=False,
        original_reference_annotation_file_opened=False, model_or_feature_tensor_opened=False,
        optimizer_steps_this_run=0, new_geometry_labels_computed=False, new_candidate_search=False,
        model_or_backbone_forward=False, DEV_used=False, independent_species_used=False,
        test_setaria_used=False, A33_outputs_opened=False,
        limitation='Closed finite TRAIN supervision artifact only, not gradient readiness, fit authorization, '
            'model benefit, full-grammar training, raw annotation or generalization.'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-job', required=True)
    main(parser.parse_args().checks_job)
