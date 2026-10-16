#!/usr/bin/env python3
"""Join completed TRAIN supervision and frozen native gaps; no forward or fitting."""
from collections import Counter
import json
import os
from pathlib import Path
import time

import diagnose_native_train_reference as D
import compare_full_grammar_train_recovery as R

ROOT = D.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-GAP-SUPERVISION'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_fresh_draws_plan.md'
SOURCES = {
    'native': (D.OUT, '13461637', 'references.jsonl'),
    'original': (ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-OBSERVATION-DIAGNOSTIC',
                 '13450304', 'reference_observations.jsonl'),
    'views': (ROOT/'outputs/M28-FULL-GRAMMAR-FINAL-TRAIN-VIEWS',
              '13450365', 'reference_view_coverage.jsonl'),
    'fresh': (ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION',
              '13453348', 'reference_coverage.jsonl')}
FULL = 'whole_reference_span_matching_fine_no_fine_disagreement'
NONE = 'no_direct_core_span_observation'


def indexed(rows):
    result = {D.record_id(r): r for r in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate primary reference identity; do not multiply joins')
    return result


def summarize(rows):
    pairs = Counter((r['original_single_full_core'], r['fresh_single_full_core']) for r in rows)
    transitions = Counter((r['original_coverage'], r['fresh_coverage']) for r in rows)
    return dict(reference_chains=len(rows), supported=sum(r['single_chain_supported'] for r in rows),
        original_model_exact_hits=sum(r['exact_primary_hit'] for r in rows),
        original_union_full_fine=sum(r['original_coverage'] == FULL for r in rows),
        fresh_union_full_fine=sum(r['fresh_coverage'] == FULL for r in rows),
        original_single_full_core=sum(r['original_single_full_core'] for r in rows),
        fresh_single_full_core=sum(r['fresh_single_full_core'] for r in rows),
        old_or_fresh_single_full_core=sum(r['original_single_full_core'] or r['fresh_single_full_core'] for r in rows),
        single_full_core_pairs=dict(both=pairs[True, True], fresh_only=pairs[False, True],
            original_only=pairs[True, False], neither=pairs[False, False]),
        original_no_direct_observation=sum(r['original_coverage'] == NONE for r in rows),
        original_no_direct_to_fresh_single_full_core=sum(
            r['original_coverage'] == NONE and r['fresh_single_full_core'] for r in rows),
        original_full_core_model_missed=sum(r['original_single_full_core'] and not r['exact_primary_hit'] for r in rows),
        original_full_core_all_view_MAP_missed=sum(
            r['original_single_full_core'] and r['original_exact_MAP_full_core_views'] == 0 for r in rows),
        coverage_transitions=[dict(original=a, fresh=b, reference_chains=n)
            for (a, b), n in sorted(transitions.items())])


def describe(rows):
    return dict(all_references=summarize(rows),
        fixed_native_groups={g:summarize([r for r in rows if R.fixed_group(r) == g]) for g in R.GROUPS},
        all_three_native_exact_anchor_gap=summarize([
            r for r in rows if not r['exact_primary_hit'] and all(r['native_exact_hits'].values())]))


def fixture():
    # Both schedules covering a chain count once, not twice. Unsupported stays in the denominator.
    rows = [dict(original_single_full_core=a, fresh_single_full_core=b,
        original_coverage=FULL if a else NONE, fresh_coverage=FULL if b else NONE,
        single_chain_supported=True, exact_primary_hit=False, original_exact_MAP_full_core_views=0)
        for a, b in ((True, True), (True, False), (False, True), (False, False))]
    rows.append(dict(rows[-1], single_chain_supported=False,
        original_coverage='unsupported_reference_without_role_join',
        fresh_coverage='unsupported_reference_without_role_join'))
    s = summarize(rows)
    assert s['reference_chains'] == 5 and s['supported'] == 4
    assert s['old_or_fresh_single_full_core'] == 3
    assert s['single_full_core_pairs'] == dict(both=1, fresh_only=1, original_only=1, neither=2)
    assert s['original_no_direct_to_fresh_single_full_core'] == 1
    assert s['original_full_core_all_view_MAP_missed'] == 2
    # Transcript IDs alone are not a join key: chromosome and strand matter.
    ids = [dict(species='a', seqid=q, strand=st, id='same') for q, st in (('q', '+'), ('q', '-'), ('r', '+'))]
    assert len(indexed(ids)) == 3
    try:
        indexed(ids + [ids[0]])
    except ValueError:
        pass
    else:
        raise AssertionError('Duplicate join was accepted')
    return 'single_core_OR_not_sum_unsupported_denominator_composite_ID_duplicate_rejection'


def load_sources():
    reports, tables = {}, {}
    for name, (folder, job, filename) in SOURCES.items():
        s = D.completed(folder)
        if (s['job_id'] != job or (folder/'JOBID').read_text().strip() != job
                or s['DEV_used'] or s['sealed_data_used'] or s['optimizer_steps_this_run'] != 0):
            raise ValueError('Completed fixed original TRAIN source differs: '+name)
        reports[name] = s
        tables[name] = indexed(D.read_rows(folder/filename))
    native = reports['native']
    if (native['source_evaluation_job'] != '13461307'
            or native['source_original_evaluation_job'] != '13443329'
            or native['reference_primary_chains'] != 22688 or native['reference_filters_added']
            or reports['original']['source_fit_job'] != '13434603'
            or reports['fresh']['source_observation_job'] != reports['original']['job_id']
            or reports['views']['checkpoint'] != native['checkpoint']
            or reports['fresh']['original_R4_labels_changed'] or reports['fresh']['primary_reference_changed']):
        raise ValueError('Source lineage, original fit or primary ruler differs')
    keys = set(tables['native'])
    if len(keys) != 22688 or any(set(t) != keys for t in tables.values()):
        raise ValueError('Preserve one-to-one coverage of all22688 primary references')
    refs = list(tables['native'].values())
    if (Counter(r['seqid'] for r in refs) != Counter(D.EXPECTED)
            or len({D.reference_key(r) for r in refs}) != len(refs)
            or D.summarize(refs) != {k:v for k, v in native['pooled'].items() if k != 'strata'}):
        raise ValueError('Fixed geometries or native/model memberships differ')
    return reports, tables


def joined_rows(tables):
    result = []
    for key, ref in tables['native'].items():
        o, v, f = (tables[n][key] for n in ('original', 'views', 'fresh'))
        old, fresh = v['complete_matching_core_views'], f['full_fine_single_core_unique_windows']
        if (not o['single_chain_supported'] == v['single_chain_supported'] == f['single_chain_supported']
                or o['exact_primary_hit'] != ref['exact_primary_hit']
                or o['coverage_class'] != v['original_union_coverage']
                or o['coverage_class'] != f['original_coverage_class']
                or o['category'] != v['whole_contig_category']
                or not 0 <= v['exact_MAP_complete_matching_core_views'] <= old
                or old < 0 or fresh < 0
                or (old > 0 and o['coverage_class'] != FULL)
                or (fresh > 0 and f['coverage_class'] != FULL)
                or (not o['single_chain_supported'] and (old or fresh))):
            raise ValueError('Saved observation/support/view identity differs: '+str(key))
        result.append(dict(ref, single_chain_supported=o['single_chain_supported'],
            original_coverage=o['coverage_class'], fresh_coverage=f['coverage_class'],
            original_single_full_core=old > 0, fresh_single_full_core=fresh > 0,
            original_full_core_unique_views=old, fresh_full_core_unique_views=fresh,
            original_exact_MAP_full_core_views=v['exact_MAP_complete_matching_core_views']))
    return result


def reconcile(rows, sources):
    for name, part in [('pooled', rows)]+[(s, [r for r in rows if r['species'] == s])
                                       for s in sources['native']['per_species']]:
        s = summarize(part)
        old = (sources['original']['pooled'] if name == 'pooled'
               else sources['original']['per_species'][name])
        fresh = sources['fresh']['reference_coverage'][name]
        views = (sources['views']['pooled'] if name == 'pooled'
                 else sources['views']['per_species'][name])['references']
        if (s['reference_chains'] != old['reference'] or s['supported'] != old['supported_role_join']
                or s['original_union_full_fine'] != old['by_coverage'][FULL]['reference']
                or s['fresh_union_full_fine'] != fresh['fresh_coverage'][FULL]
                or s['original_single_full_core'] != views['whole_span_matching_fine_in_at_least_one_original_core']['reference']
                or s['fresh_single_full_core'] != fresh['any_full_fine_single_core_reference']):
            raise ValueError('Full or species counts do not reproduce completed supervision sources')


def main():
    began = time.perf_counter()
    checked = fixture()
    sources, tables = load_sources()
    rows = joined_rows(tables)
    reconcile(rows, sources)
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], plan=PLAN, fixture=checked,
        status='descriptive_original_TRAIN_native_gap_supervision_not_model_performance',
        source_jobs={n:s['job_id'] for n, s in sources.items()},
        reference_primary_chains=len(rows), reference_policy=sources['native']['reference_policy'],
        one_to_one_complete_primary_join=True, full_and_species_counts_reproduce_completed_sources=True,
        pooled=describe(rows), per_species={s:describe([r for r in rows if r['species'] == s])
            for s in sources['native']['per_species']},
        old_or_fresh_means_any_single_core_not_combined_fine_union=True,
        original_core_repeated_passes=3, fresh_schedule_draws=4608,
        prediction_GFF_read=False, DNA_read=False, real_reference_GFF_read=False, checkpoint_read=False,
        saved_reference_profiles_read=True, new_model_inference=False, predictions_modified=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        scientific_winner_or_DEV_candidate_selected=False,
        limitations=[
            'Core supervision is not full encoder, peripheral partition or upstream exposure.',
            'Fresh-only and original-only describe schedule sets, not lost cumulative training knowledge.',
            'Old-or-fresh OR counts one complete matching core; it cannot reconstruct cumulative fine union from partial cores.',
            'Fixed native-gap groups use original prediction outcomes and are descriptive, not causal strata.',
            'Coverage counts are not new prediction recall, F1, generalization or evidence that fresh training succeeds.'],
        seconds=time.perf_counter()-began)
    with (OUT/'references.jsonl').open('x') as handle:
        for r in rows:
            handle.write(json.dumps(r, allow_nan=False)+'\n')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
