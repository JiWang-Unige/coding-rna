#!/usr/bin/env python3
"""Post-outcome original-TRAIN reference recovery: dense versus low-rank q/v."""
from collections import Counter
import json
import os
import time
import compare_full_grammar_joint_lora_recovery as R

D = R.D
ROOT = D.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-DENSE-QV-LORA-RECOVERY'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_dense_qv_result.md'
SOURCES = {
    'baseline': ('JOINT-LORA', '13499823', '13499817'),
    'adapted': ('DENSE-QV', '13510351', '13510324'),
}


def main():
    began = time.perf_counter()
    R.checks()
    reports = {}
    for arm, (prefix, evaluation, fit) in SOURCES.items():
        folder = ROOT/('outputs/M28-FULL-GRAMMAR-'+prefix+'-TRAIN-DIAGNOSTIC')
        report = R.completed(folder, evaluation)
        expected_checkpoint = ROOT/('outputs/M28-FULL-GRAMMAR-'+prefix+'-FIT/step_009216.pt')
        if (report['checkpoint'] != str(expected_checkpoint)
                or report['pooled']['exact_chain']['reference'] != 22688
                or set(report['chromosomes']) != set(D.EXPECTED)
                or report['reference_used_in_inference']
                or not report['unsupported_reference_strata_retained']
                or not report['training_eligibility_does_not_filter_reference']
                or report['DEV_used'] or report['sealed_data_used']
                or report['optimizer_steps_this_run']):
            raise ValueError('Two fixed completed final9216 original-TRAIN sources required')
        inputs = json.loads((folder/'paired_input_match.json').read_text())
        if (inputs['job_id'] != evaluation or inputs['adapted_fit_job'] != fit
                or inputs['baseline_fit_job'] != '13494895'
                or not inputs['original_six_TRAIN_inputs_directly_equal']
                or not inputs['comparison_before_reference_access']
                or inputs['DEV_used'] or inputs['sealed_data_used']):
            raise ValueError('Both saved direct-DNA receipts must match the same original TRAIN')
        reports[arm] = report
    baseline = reports['baseline']
    reference = R.completed(D.OUT, '13461637')
    if (reference['reference_policy'] != baseline['reference_policy']
            or reference['reference_primary_chains'] != 22688
            or reference['reference_filters_added']
            or not reference['exact_hit_counts_reproduce_completed_evaluations']
            or reference['DEV_used'] or reference['sealed_data_used']
            or reference['optimizer_steps_this_run']):
        raise ValueError('Completed unchanged saved primary-reference population required')
    rows = D.read_rows(D.OUT/'references.jsonl')
    scope = {(r['species'], seqid):D.EXPECTED[seqid]
             for seqid,r in baseline['chromosomes'].items()}
    if (len(rows) != 22688 or len({D.reference_key(r) for r in rows}) != 22688
            or Counter((r['species'],r['seqid']) for r in rows) != Counter(scope)
            or any(not r['CDS'] or r['CDS_segments'] != len(r['CDS'])
                   or r['span_bp'] != r['CDS'][-1][1]-r['CDS'][0][0] for r in rows)):
        raise ValueError('Saved primary geometry and six-contig population differ')
    for seqid, official in baseline['chromosomes'].items():
        species,length = official['species'],official['length']
        part = [r for r in rows if (r['species'],r['seqid']) == (species,seqid)]
        refs = {D.reference_key(r) for r in part}
        for arm, report in reports.items():
            source = report['chromosomes'][seqid]
            if (report['reference_policy'] != baseline['reference_policy']
                    or (source['species'],source['seqid'],source['length']) != (species,seqid,length)
                    or any(source[k] != official[k] for k in
                           ('reference_primary_records','reference_assessable_isoform_records',
                            'retained_unsupported_reference_strata'))):
                raise ValueError('Paired reference population or chromosome differs')
            prefix = SOURCES[arm][0]
            path = ROOT/source['raw_prediction_path']
            expected = ROOT/('outputs/M28-FULL-GRAMMAR-'+prefix+'-TRAIN-RAW')/seqid/'prediction/full_grammar.gff3'
            if path.resolve() != expected.resolve():
                raise ValueError('Actual fixed source prediction path differs')
            keys = R.checked_predictions(path,species,seqid,length,refs,source)
            for row in part:
                row[arm+'_hit'] = D.reference_key(row) in keys
    views = {'pooled':rows,
             **{s:[r for r in rows if r['species']==s] for s in baseline['per_species']},
             **{q:[r for r in rows if r['seqid']==q] for q in D.EXPECTED}}
    scopes = {name:R.describe(part) for name,part in views.items()}
    for name,value in scopes.items():
        for arm,report in reports.items():
            source = report['pooled'] if name=='pooled' else report['per_species'][name] if name in report['per_species'] else report['chromosomes'][name]
            if (value[arm+'_hits'] != source['exact_chain']['tp']
                    or value['reference_chains'] != source['exact_chain']['reference']):
                raise ValueError('Replay must reconstruct every saved primary-hit count')
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        arm_meaning=dict(baseline='low_rank_qv',adapted='dense_qv'),
        source_evaluation_jobs={a:v[1] for a,v in SOURCES.items()},
        source_fit_jobs={a:v[2] for a,v in SOURCES.items()},
        source_reference_job='13461637',reference_policy=baseline['reference_policy'],
        reference_primary_chains=22688,scopes=scopes,
        artificial_four_way_and_bin_boundary_checks_passed=True,
        exact_hits_unique_predictions_and_duplicates_reproduce_completed_evaluations=True,
        every_reference_partition_reconciled=True,reference_filters_added=False,
        source_DNA_equality_receipts_to_same_original_TRAIN_verified=True,
        DNA_read=False,real_reference_GFF_read=False,saved_primary_reference_profile_read=True,
        prediction_GFF_read=True,checkpoint_read=False,native_forward=False,
        predictions_modified=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,frozen_DEV_family_changed=False,
        post_outcome_exploratory_diagnostic=True,seconds=time.perf_counter()-began,
        limitations=['Post-M-outcome exploratory accounting, not independent confirmation.',
            'Reference intersections and fixed pre-existing bins measure recall, not subgroup precision or F1.',
            'Rank, learning rate, rounding and optimization trajectories differ; not a rank-only causal experiment.',
            'Original TRAIN only; not novel-gene validation, a selected ensemble or frozen-DEV extension.'])
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(completed=True,reference_primary_chains=22688)),flush=True)


if __name__ == '__main__':
    main()
