#!/usr/bin/env python3
"""Complete codon-minus-baseline paired TRAIN endpoints, without selection."""
import argparse
import json
import os
import time
import infer_full_grammar_codon_train as D
import infer_full_grammar_fresh_joint24_train as B
from compare_full_grammar_uniform_coarse_train import checked_prf, SPECIES, SEQIDS, METRICS, ENDPOINTS

ROOT = D.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-PAIRED'


def paired_scopes(reports):
    if reports['codon']['reference_policy'] != reports['baseline']['reference_policy']:
        raise ValueError('Primary reference policy differs')
    for species, seqid, length in D.I.SCOPE:
        a, b = (reports[name]['chromosomes'][seqid] for name in ('codon', 'baseline'))
        if any((r['species'], r['seqid'], r['length']) != (species, seqid, length) for r in (a, b)):
            raise ValueError('Paired chromosome geometry differs')
        if (a['reference_primary_records'] != b['reference_primary_records']
                or a['reference_assessable_isoform_records'] != b['reference_assessable_isoform_records']
                or a['retained_unsupported_reference_strata'] != b['retained_unsupported_reference_strata']):
            raise ValueError('Paired reference population or supported strata changed')
    scopes = {}
    for scope in (*SEQIDS, *SPECIES, 'pooled'):
        arms = {}
        for name, report in reports.items():
            row = report['pooled'] if scope == 'pooled' else report['per_species'][scope] if scope in SPECIES else report['chromosomes'][scope]
            arms[name] = {endpoint:checked_prf(row[endpoint]) for endpoint in ENDPOINTS}
            if scope in SPECIES or scope == 'pooled':
                parts = list(report['chromosomes'].values()) if scope == 'pooled' else [r for r in report['chromosomes'].values() if r['species'] == scope]
                if any(arms[name][e][k] != sum(r[e][k] for r in parts) for e in ENDPOINTS for k in ('tp', 'predicted', 'reference')):
                    raise ValueError('Species and pooled counts must reconstruct chromosomes')
            arms[name]['background'] = row['background']
        if any(arms['codon'][e]['reference'] != arms['baseline'][e]['reference'] for e in ENDPOINTS):
            raise ValueError('Keep common retained reference denominators')
        differences = {e:dict(
            counts={k:arms['codon'][e][k]-arms['baseline'][e][k] for k in ('tp', 'predicted', 'reference')},
            PRF_percentage_points={k:100*(arms['codon'][e][k]-arms['baseline'][e][k]) for k in METRICS}) for e in ENDPOINTS}
        scopes[scope] = dict(arms=arms, codon_minus_baseline=differences)
    for scope in ('macro_species', 'worst_species'):
        arms = {}
        for name, report in reports.items():
            arms[name] = {e:{k:(sum(report['per_species'][s][e][k] for s in SPECIES)/2 if scope == 'macro_species'
                else min(report['per_species'][s][e][k] for s in SPECIES)) for k in METRICS} for e in ENDPOINTS}
        scopes[scope] = dict(arms=arms, codon_minus_baseline_PRF_percentage_points={e:{k:
            100*(arms['codon'][e][k]-arms['baseline'][e][k]) for k in METRICS} for e in ENDPOINTS})
    if any(scopes['pooled']['arms'][name]['exact_chain']['reference'] != 22688 for name in reports):
        raise ValueError('Retain all22688 primary TRAIN chains')
    return scopes


def main(evaluation_job, baseline_evaluation_job):
    began = time.perf_counter()
    reports = {}
    for name, root, job, checkpoint in (('codon', D.EVAL, evaluation_job, D.CHECKPOINT),
                                       ('baseline', B.EVAL, baseline_evaluation_job, B.CHECKPOINT)):
        value = json.loads((root/'summary.json').read_text())
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or (root/'JOBID').read_text().strip() != job
                or not value['completed'] or value['job_id'] != job or value['checkpoint'] != str(checkpoint)
                or value['DEV_used'] or value['sealed_data_used'] or value['optimizer_steps_this_run']
                or value['reference_used_in_inference'] or not value['unsupported_reference_strata_retained']
                or set(value['chromosomes']) != set(SEQIDS) or set(value['per_species']) != set(SPECIES)):
            raise ValueError('Two completed strict original-TRAIN final9216 sources required')
        reports[name] = value
    inputs = json.loads((D.EVAL/'paired_input_match.json').read_text())
    if (inputs['job_id'] != evaluation_job or inputs['protocol'] != D.A.PROTOCOL
            or not inputs['original_six_TRAIN_inputs_directly_equal'] or not inputs['comparison_before_reference_access']
            or inputs['DEV_used'] or inputs['sealed_data_used']):
        raise ValueError('Completed paired DNA match before reference access required')
    D.source(inputs['codon_fit_job'])
    B.source(inputs['baseline_fit_job'])
    scopes = paired_scopes(reports)
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=D.A.PROTOCOL,
        source_jobs=dict(codon=evaluation_job, baseline=baseline_evaluation_job),
        source_fit_jobs=dict(codon=inputs['codon_fit_job'], baseline=inputs['baseline_fit_job']), scopes=scopes,
        all_paired_caller_inputs_and_normalized_original_DNA_directly_equal=True,
        DNA_comparison_before_codon_reference_access=True, counts_and_PRF_reconciled=True,
        reference_policy=reports['codon']['reference_policy'], optimizer_steps_this_run=0,
        native_forward=False, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        limitation='single-fit additional-learning codon-energy contrast including64 new parameters, new AdamW age and joint clipping; not independent repetition, generalization or model selection')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    lines = ['# Original TRAIN complete CDS codon energy comparison', '',
        '| Scope | Codon chain F1 % | Baseline chain F1 % | Difference pp |',
        '| --- | ---: | ---: | ---: |']
    for scope in (*SPECIES, 'pooled'):
        r = scopes[scope]
        values = [100*r['arms'][name]['exact_chain']['f1'] for name in ('codon', 'baseline')]
        lines.append('| '+scope+' | '+' | '.join(f'{v:.4f}' for v in values)+
                     f" | {r['codon_minus_baseline']['exact_chain']['PRF_percentage_points']['f1']:+.4f} |")
    lines.extend(['', report['limitation']])
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(completed=True, job_id=report['job_id'], source_jobs=report['source_jobs'])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation-job', required=True)
    parser.add_argument('--baseline-evaluation-job', required=True)
    args = parser.parse_args()
    main(args.evaluation_job, args.baseline_evaluation_job)
