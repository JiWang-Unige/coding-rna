#!/usr/bin/env python3
"""Fixed supervision strata: completed LOCAL rehearsal or matched FRESH/STRUCT9216."""
import argparse
import json
import os
import time

import diagnose_full_grammar_native_gap_supervision as G
import compare_full_grammar_train_recovery as R

FRESH = 'M28_FRESH_DRAWS_final9216'
STRATA = ('both_complete', 'fresh_only_complete', 'original_only_complete', 'neither_complete')
CHECKS = R.ROOT/'outputs/M28-FULL-GRAMMAR-SUPERVISION-RECOVERY-CHECKS'
OUT = R.ROOT/'outputs/M28-FULL-GRAMMAR-SAMPLING-SUPERVISION-RECOVERY'


def stratum(row):
    a, b = row['original_single_full_core'], row['fresh_single_full_core']
    return 'both_complete' if a and b else 'original_only_complete' if a else 'fresh_only_complete' if b else 'neither_complete'


def describe(rows, candidate, control):
    consensus = [r for r in rows if not r['exact_primary_hit'] and all(r['native_exact_hits'].values())]
    return dict(all_references=R.pair(rows, candidate, control),
        supervision_strata={k:R.pair([r for r in rows if stratum(r) == k], candidate, control) for k in STRATA},
        fixed_native_groups={g:dict(all_references=R.pair(part, candidate, control),
            supervision_strata={k:R.pair([r for r in part if stratum(r) == k], candidate, control) for k in STRATA})
            for g in R.GROUPS for part in [[r for r in rows if R.fixed_group(r) == g]]},
        all_three_native_exact_anchor_gap=dict(all_references=R.pair(consensus, candidate, control),
            supervision_strata={k:R.pair([r for r in consensus if stratum(r) == k], candidate, control) for k in STRATA}))


def verify_partition(report):
    fields = ('reference_chains', 'control_hits', 'candidate_hits', 'retained', 'gained', 'lost', 'still_missed', 'net_reference_hits')
    for group in [report]+list(report['fixed_native_groups'].values())+[report['all_three_native_exact_anchor_gap']]:
        all_refs = group['all_references']
        for k in fields:
            if sum(r[k] for r in group['supervision_strata'].values()) != all_refs[k]:
                raise ValueError('Supervision strata must partition the full reference/hit counts')


def fixture():
    rows = []
    for old, fresh in ((True, True), (False, True), (True, False), (False, False)):
        for a, b in ((True, True), (False, True), (True, False), (False, False)):
            rows.append(dict(original_single_full_core=old, fresh_single_full_core=fresh,
                exact_primary_hit=a, native_exact_hits={m:b for m in R.D.METHODS},
                model_hits={'control':a, 'candidate':b}))
    report = describe(rows, 'candidate', 'control')
    verify_partition(report)
    assert report['all_references']['reference_chains'] == 16
    assert [report['all_references'][k] for k in ('retained', 'gained', 'lost', 'still_missed')] == [4]*4
    assert all([r[k] for k in ('retained', 'gained', 'lost', 'still_missed')] == [1]*4
               for r in report['supervision_strata'].values())
    assert R.pair([], 'candidate', 'control')['reference_recall'] is None
    return 'four_supervision_strata_four_hit_states_and_empty_groups'


def join_supervision(rows, control):
    s = R.D.completed(G.OUT)
    if (s['job_id'] != '13462169' or (G.OUT/'JOBID').read_text().strip() != s['job_id']
            or s['reference_policy'] != control['reference_policy']
            or not s['one_to_one_complete_primary_join'] or s['DEV_used'] or s['sealed_data_used']
            or s['optimizer_steps_this_run'] != 0):
        raise ValueError('Fixed completed supervision source required')
    saved = G.indexed(R.D.read_rows(G.OUT/'references.jsonl'))
    originals = G.indexed(rows)
    if len(saved) != 22688 or set(saved) != set(originals):
        raise ValueError('Complete one-to-one primary reference supervision join required')
    for key, row in originals.items():
        other = saved[key]
        if (R.D.reference_key(row) != R.D.reference_key(other)
                or row['exact_primary_hit'] != other['exact_primary_hit']
                or row['native_exact_hits'] != other['native_exact_hits']):
            raise ValueError('Frozen reference geometry or original/native outcomes differ')
        for k in ('single_chain_supported', 'original_coverage', 'fresh_coverage',
                  'original_single_full_core', 'fresh_single_full_core', 'original_full_core_unique_views',
                  'fresh_full_core_unique_views', 'original_exact_MAP_full_core_views'):
            row[k] = other[k]
    if G.describe(rows) != s['pooled'] or any(
            G.describe([r for r in rows if r['species'] == sp]) != official for sp, official in s['per_species'].items()):
        raise ValueError('Completed supervision denominators and fixed gaps must reproduce')
    return s


def main(check_local):
    began = time.perf_counter()
    checked = fixture()
    control, candidate = (R.ANCHOR, R.LOCAL) if check_local else (R.CONTROL, FRESH)
    labels = (R.ANCHOR, R.LOCAL) if check_local else (R.ANCHOR, control, candidate)
    # Complete fixed evaluations close before any saved profile/GFF read.
    C, specs, reports = R.load_reports(labels)
    rows, paired = R.source_records(reports[R.ANCHOR], C)
    supervision = join_supervision(rows, reports[R.ANCHOR])
    for label in labels[1:]:
        R.replay(rows, label, reports[label], C)
    pooled = describe(rows, candidate, control)
    per_species = {sp:describe([r for r in rows if r['species'] == sp], candidate, control) for sp in C.SPECIES}
    for report in [pooled]+list(per_species.values()):
        verify_partition(report)
    if pooled['all_references']['net_reference_hits'] != (
            reports[candidate]['pooled']['exact_chain']['tp']-reports[control]['pooled']['exact_chain']['tp']):
        raise ValueError('Paired reference changes must match complete evaluation TP changes')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], plan=G.PLAN,
        status='completed_LOCAL_rehearsal_not_FRESH_vs_STRUCT_result' if check_local else 'completed_matched_FRESH_vs_STRUCT_original_TRAIN_supervision_recovery',
        fixture=checked, fixture_only=False, check_local=check_local,
        formal_sampling_comparison_completed=not check_local,
        control=control, candidate=candidate, source_evaluation_jobs={m:s['job_id'] for m, s in reports.items()},
        source_paired_reference_job=paired['job_id'], source_supervision_job=supervision['job_id'],
        reference_primary_chains=len(rows), reference_policy=reports[R.ANCHOR]['reference_policy'],
        supervision_and_reference_counts_reproduce_completed_sources=True,
        fixed_strata_unchanged_by_new_predictions=True,
        pooled=pooled, per_species=per_species,
        measurements={m:dict(optimizer_steps=specs[m][3],
            **{k:reports[m][k] for k in ('checkpoint', 'pooled', 'per_species',
                'macro_species_exact_chain_F1', 'worst_species_exact_chain_F1')}) for m in (control, candidate)},
        DNA_read=False, real_reference_GFF_read=False, checkpoint_read=False, prediction_GFF_read=True,
        new_model_inference=False, predictions_modified=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, scientific_winner_or_DEV_candidate_selected=False,
        limitations=[
            'Reference strata report recall accounting, not subset prediction precision or F1.',
            'LOCAL rehearsal has extra updates and different objective; future fresh-supervision flags do not describe LOCAL training.',
            'Matched FRESH/STRUCT updates do not isolate coverage from proposal variance and importance effects.',
            'Fixed native-gap groups use original outcomes and remain descriptive, not causal strata.',
            'No direct complete-core supervision is not zero encoder, partition or upstream exposure.',
            'This supplemental readout is not an extra DEV gate and changes no frozen candidate family.'],
        seconds=time.perf_counter()-began)
    out = CHECKS if check_local else OUT
    with (out/'references.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-local', action='store_true', help='Completed LOCAL rehearsal only, not sampling comparison')
    main(parser.parse_args().check_local)
