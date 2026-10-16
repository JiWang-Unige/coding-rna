#!/usr/bin/env python3
"""Complete-core fine supervision eligibility in the existing six-TRAIN pool."""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import json
import os
import time

import numpy as np
import torch
import prepare_full_grammar_fresh_draws as P
import diagnose_full_grammar_native_gap_supervision as G

OUT = P.ROOT/'outputs/M28-FULL-GRAMMAR-POPULATION-SUPERVISION'
CHECKS = P.ROOT/'outputs/M28-FULL-GRAMMAR-POPULATION-SUPERVISION-CHECKS'
PLAN = P.PROTOCOL
CATEGORIES = ('outside_fixed_single_chain_grammar', 'span_exceeds_existing_core_width',
              'no_containing_existing_core', 'containing_cores_without_complete_matching_fine',
              'complete_matching_fine_core_available')


def containing(pool, starts, lo, hi):
    return [w for w in pool[bisect_left(starts, hi-24576):bisect_right(starts, lo)]
            if w['end'] >= hi]


def role_counts(fine, region, roles, window, lo, hi):
    offset = lo-window['start'] if window['strand'] == '+' else window['end']-hi
    observed = fine[offset:offset+len(roles)]
    coarse = region[offset:offset+len(roles)]
    if not 0 <= offset <= len(fine)-len(roles) or len(roles) != hi-lo:
        raise ValueError('Whole reference span must lie in oriented valid core')
    matching = int(np.count_nonzero(observed == roles))
    return dict(matching_fine_rows=matching,
        disagreeing_fine_rows=int(np.count_nonzero((observed >= 0) & (observed != roles))),
        coarse_only_rows=int(np.count_nonzero((observed < 0) & (coarse >= 0))),
        fully_unknown_rows=int(np.count_nonzero((observed < 0) & (coarse < 0))),
        complete_matching_fine=matching == len(roles))


def fixture():
    genome = 'CC'+'ATGCCCTAA'+'C'*9
    expected = np.array([13,13,13,1,2,3,14,14,14])
    for strand, dna, lo, hi, start, end in (
            ('+', genome, 2, 11, 0, 14),
            ('-', P.F.Q.M.rc(genome), 9, 18, 6, 20)):
        tx = dict(id='t', gene_id='g', strand=strand, CDS=[(lo,hi,0)],
                  training_primary=True, CDS_complete_supported=True)
        w = dict(start=start, end=end, width=end-start, valid_bases=end-start, strand=strand,
            overlapping_transcript_ids=['t'], positive_ids=['t'],
            fine_label_unknown_intervals=[], chain_unknown_intervals=[],
            no_annotated_gene_either_strand=False)
        raw = dna[start:end]
        core = raw if strand == '+' else P.F.Q.M.rc(raw)
        target = P.O.build_targets(w, {'t':tx}, core)
        roles = P.E.role_path(P.E.oriented_parts(tx, len(dna)))
        assert np.array_equal(roles, expected)
        counts = role_counts(target['seg'], target['region'], roles, w, lo, hi)
        assert counts == dict(matching_fine_rows=9, disagreeing_fine_rows=0,
            coarse_only_rows=0, fully_unknown_rows=0, complete_matching_fine=True)
        # Real R4 unknown and alternative-isoform masking cannot count as fine supervision.
        masked = {**w, 'fine_label_unknown_intervals':[(2,3)]}
        altered = P.O.build_targets(masked, {'t':tx}, core)
        assert not role_counts(altered['seg'], altered['region'], roles, w, lo, hi)['complete_matching_fine']
        alt = {**tx, 'id':'alt', 'CDS_complete_supported':False}
        alternatives = {**w, 'overlapping_transcript_ids':['t','alt']}
        altered = P.O.build_targets(alternatives, {'t':tx, 'alt':alt}, core)
        assert not role_counts(altered['seg'], altered['region'], roles, w, lo, hi)['complete_matching_fine']
    pool = [dict(id='left', start=0, end=24576),
            dict(id='right', start=12288, end=36864),
            dict(id='tail', start=24576, end=40000)]
    starts = [w['start'] for w in pool]
    assert [w['id'] for w in containing(pool, starts, 12288, 24576)] == ['left','right']
    assert containing(pool, starts, 0, 24577) == []
    assert [w['id'] for w in containing(pool, starts, 39000, 40000)] == ['tail']
    assert containing(pool, starts, 39000, 40001) == []
    return 'plus_minus_role_offsets_R4_unknown_isoform_masks_exact_edges_partial_tail'


def describe(rows):
    eligible = [r for r in rows if r['population_matching_core_window_ids']]
    return dict(reference_chains=len(rows), supported=sum(r['single_chain_supported'] for r in rows),
        categories=dict(Counter(r['population_supervision_category'] for r in rows)),
        population_complete_core_available=len(eligible),
        original_single_full_core=sum(r['original_single_full_core'] for r in rows),
        fresh_single_full_core=sum(r['fresh_single_full_core'] for r in rows),
        old_or_fresh_single_full_core=sum(r['original_single_full_core'] or r['fresh_single_full_core'] for r in rows),
        available_but_absent_from_both_schedules=sum(
            not r['original_single_full_core'] and not r['fresh_single_full_core'] for r in eligible),
        matched_window_reference_pairs=sum(len(r['population_matching_core_window_ids']) for r in rows))


def grouped(rows):
    return dict(all_references=describe(rows),
        fixed_native_groups={name:describe([r for r in rows if G.R.fixed_group(r) == name])
                             for name in G.R.GROUPS},
        all_three_native_exact_anchor_gap=describe([
            r for r in rows if not r['exact_primary_hit'] and all(r['native_exact_hits'].values())]))


def main(fixture_only):
    began = time.perf_counter()
    checked = fixture()
    if fixture_only:
        P.write(CHECKS/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
            fixture=checked, fixture_only=True, research_data_read=False,
            optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
        print(checked, flush=True)
        return
    source = G.D.completed(G.OUT)
    if (source['job_id'] != '13462169' or (G.OUT/'JOBID').read_text().strip() != source['job_id']
            or source['reference_primary_chains'] != 22688 or not source['one_to_one_complete_primary_join']
            or source['optimizer_steps_this_run'] != 0 or source['DEV_used'] or source['sealed_data_used']):
        raise ValueError('Completed fixed supervision source required')
    rows = G.D.read_rows(G.OUT/'references.jsonl')
    if len(G.indexed(rows)) != 22688 or sum(r['single_chain_supported'] for r in rows) != 22479:
        raise ValueError('All 22688 original primary references, including unsupported, required')
    windows, strata, catalogs, original = P.load_population()
    fresh_report = G.D.completed(P.OUT)
    fresh_plan = list(P.rows(P.OUT/'training_plan.jsonl'))
    if (fresh_report['job_id'] != '13453348' or fresh_report['DEV_used'] or fresh_report['sealed_data_used']
            or len(fresh_plan) != 4608 or fresh_report['optimizer_steps_this_run'] != 0):
        raise ValueError('Completed fixed fresh draw schedule required')
    old_ids = {r['window_id'] for r in original}
    fresh_ids = {r['window_id'] for r in fresh_plan}
    if len(old_ids) != 1491 or len(fresh_ids) != 4219 or not (old_ids | fresh_ids) <= set(windows):
        raise ValueError('Original and fresh schedules must remain subsets of the natural pool')
    sequences, descriptors = P.F.P.read_train_DNA()
    tasks, roles_by_id = defaultdict(list), {}
    for row in rows:
        row['population_containing_core_count'] = 0
        row['population_matching_core_window_ids'] = []
        row['population_containing_core_row_counts'] = []
        row['population_supervision_category'] = 'outside_fixed_single_chain_grammar'
        if not row['single_chain_supported']:
            continue
        key = row['species'], row['seqid'], row['strand']
        pool = strata[key]
        lo, hi = row['CDS'][0][0], row['CDS'][-1][1]
        if hi-lo > 24576:
            row['population_supervision_category'] = 'span_exceeds_existing_core_width'
            continue
        possible = containing(pool, [w['start'] for w in pool], lo, hi)
        row['population_containing_core_count'] = len(possible)
        row['population_supervision_category'] = ('containing_cores_without_complete_matching_fine'
                                                 if possible else 'no_containing_existing_core')
        if possible:
            roles_by_id[G.D.record_id(row)] = P.E.role_path(P.E.oriented_parts(row, len(sequences[key[0]][key[1]])))
        for w in possible:
            tasks[w['id']].append(row)
    P.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], plan=PLAN,
        full_existing_TRAIN_windows=len(windows), candidate_containing_windows=len(tasks),
        source_supervision_job=source['job_id'], reference_primary_chains=len(rows)))
    eligibility = {}
    for ordinal, (ident, refs) in enumerate(tasks.items()):
        w = windows[ident]
        chromosome = sequences[w['species']][w['seqid']]
        raw = chromosome[w['start']:w['end']]
        core = raw if w['strand'] == '+' else P.F.Q.M.rc(raw)
        target = P.O.build_targets(w, catalogs[w['species']], core)
        n = w['valid_bases']
        fine, region = target['seg'][:n], target['region'][:n]
        eligibility[ident] = []
        for row in refs:
            key = G.D.record_id(row)
            lo, hi = row['CDS'][0][0], row['CDS'][-1][1]
            counts = role_counts(fine, region, roles_by_id[key], w, lo, hi)
            row['population_containing_core_row_counts'].append(dict(window_id=ident, **counts))
            if counts['complete_matching_fine']:
                row['population_matching_core_window_ids'].append(ident)
                eligibility[ident].append(list(key))
                row['population_supervision_category'] = 'complete_matching_fine_core_available'
        if (ordinal+1) % 512 == 0:
            print(json.dumps(dict(checked_windows=ordinal+1, total_containing_windows=len(tasks),
                                  seconds=time.perf_counter()-began)), flush=True)
    for row in rows:
        ids = set(row['population_matching_core_window_ids'])
        if (len(ids) != len(row['population_matching_core_window_ids'])
                or len(ids & old_ids) != row['original_full_core_unique_views']
                or len(ids & fresh_ids) != row['fresh_full_core_unique_views']):
            raise ValueError('Population eligibility must reproduce every saved original/fresh reference count')
    pooled = grouped(rows)
    per_species = {s:grouped([r for r in rows if r['species'] == s]) for s in P.Q.ALLOW}
    if (pooled['all_references']['old_or_fresh_single_full_core'] != 6007
            or sum(len(v) for v in eligibility.values()) != pooled['all_references']['matched_window_reference_pairs']
            or sum(per_species[s]['all_references']['population_complete_core_available'] for s in per_species)
                != pooled['all_references']['population_complete_core_available']):
        raise ValueError('Complete reference and window incidence accounting must agree')
    with (OUT/'references.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    with (OUT/'window_eligibility.jsonl').open('x') as handle:
        for ident, w in windows.items():
            handle.write(json.dumps(dict(window_id=ident, species=w['species'], seqid=w['seqid'],
                strand=w['strand'], start=w['start'], end=w['end'],
                complete_reference_ids=eligibility.get(ident, [])), allow_nan=False)+'\n')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], plan=PLAN,
        status='complete_existing_original_TRAIN_population_fine_supervision_census_not_model_performance',
        fixture=checked, fixture_only=False, full_existing_TRAIN_windows=len(windows),
        containing_windows_evaluated=len(tasks), reference_primary_chains=len(rows),
        source_supervision_job=source['job_id'], source_fresh_preparation_job=fresh_report['job_id'],
        reference_policy=source['reference_policy'], pooled=pooled, per_species=per_species,
        all_saved_original_and_fresh_unique_core_counts_reproduced=True,
        reference_filters_added=False, labels_changed=False, actual_optimizer_draws_created=False,
        sampling_or_candidate_selected=False, baseline_outcomes_used_to_define_window_eligibility=False,
        model_performance_or_mechanism_causality_proven=False,
        DNA_read=True, existing_TRAIN_R4_labels_read=True, real_reference_GFF_read=False,
        checkpoint_read=False, native_forward=False, feature_tensors_read=False,
        predictions_modified=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        original_TRAIN_DNA_sources=descriptors, seconds=time.perf_counter()-began,
        limitations=['Matching fine labels within one existing core is not a complete model prediction.',
            'Unavailable fine supervision does not mean absent encoder/partition/upstream exposure.',
            'This census does not prove every population fine event feasible in the full grammar.',
            'Native-gap groups are descriptive, never used to choose training windows.',
            'No new draw schedule, fit, DEV candidate, threshold or checkpoint is selected.'])
    P.write(OUT/'summary.json', report)
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-only', action='store_true')
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.fixture_only)
