#!/usr/bin/env python3
"""All fixed paired coarse-information additional-fit TRAIN endpoints."""
import argparse
import json
import math
import os
import time
import infer_full_grammar_uniform_coarse_train as U
import infer_full_grammar_fresh_joint24_train as B

ROOT = U.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-PAIRED'
SPECIES = ('arabidopsis_thaliana','oryza_sativa')
SEQIDS = tuple(x[1] for x in U.I.SCOPE)
METRICS = ('precision','recall','f1')
ENDPOINTS = ('exact_chain','CDS_base_unstranded_union')


def checked_prf(value):
    tp,predicted,reference = (value[k] for k in ('tp','predicted','reference'))
    if not 0 <= tp <= min(predicted,reference):
        raise ValueError('Invalid endpoint counts')
    expected = dict(precision=tp/predicted if predicted else 0.,recall=tp/reference if reference else 0.,
        f1=2*tp/(predicted+reference) if predicted+reference else 0.)
    if any(not math.isclose(value[k],v,rel_tol=1e-12,abs_tol=1e-12) for k,v in expected.items()):
        raise ValueError('Endpoint PRF must reconstruct exact counts')
    return dict(value)


def main(evaluation_job, informative_evaluation_job):
    began = time.perf_counter()
    reports = {}
    for name,root,job,checkpoint in (('informative',B.EVAL,informative_evaluation_job,B.CHECKPOINT),
                                    ('uniform',U.EVAL,evaluation_job,U.CHECKPOINT)):
        value = json.loads((root/'summary.json').read_text())
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or (root/'JOBID').read_text().strip() != job
                or not value['completed'] or value['job_id'] != job or value['checkpoint'] != str(checkpoint)
                or value['DEV_used'] or value['sealed_data_used'] or value['optimizer_steps_this_run']
                or value['reference_used_in_inference'] or not value['unsupported_reference_strata_retained']
                or set(value['chromosomes']) != set(SEQIDS) or set(value['per_species']) != set(SPECIES)):
            raise ValueError('Two completed strict original-TRAIN final9216 sources required')
        reports[name] = value
    inputs = json.loads((U.EVAL/'paired_input_match.json').read_text())
    if (inputs['job_id'] != evaluation_job or inputs['protocol'] != U.A.PROTOCOL
            or not inputs['original_six_TRAIN_inputs_directly_equal']
            or not inputs['comparison_before_reference_access'] or inputs['DEV_used'] or inputs['sealed_data_used']):
        raise ValueError('Completed paired DNA match before reference access required')
    U.source(inputs['uniform_fit_job'])
    B.source(inputs['informative_fit_job'])
    if reports['uniform']['reference_policy'] != reports['informative']['reference_policy']:
        raise ValueError('Primary reference policy differs')
    for species,seqid,length in U.I.SCOPE:
        a,b = (reports[name]['chromosomes'][seqid] for name in ('informative','uniform'))
        if any((r['species'],r['seqid'],r['length']) != (species,seqid,length) for r in (a,b)):
            raise ValueError('Paired chromosome geometry differs')
        if (a['reference_primary_records'] != b['reference_primary_records']
                or a['reference_assessable_isoform_records'] != b['reference_assessable_isoform_records']
                or a['retained_unsupported_reference_strata'] != b['retained_unsupported_reference_strata']):
            raise ValueError('Paired reference population or supported strata changed')
    scopes = {}
    for scope in (*SEQIDS,*SPECIES,'pooled'):
        arms = {}
        for name,report in reports.items():
            row = report['pooled'] if scope=='pooled' else report['per_species'][scope] if scope in SPECIES else report['chromosomes'][scope]
            arms[name] = {endpoint:checked_prf(row[endpoint]) for endpoint in ENDPOINTS}
            if scope in SPECIES or scope=='pooled':
                parts = list(report['chromosomes'].values()) if scope=='pooled' else [r for r in report['chromosomes'].values() if r['species']==scope]
                if any(arms[name][e][k] != sum(r[e][k] for r in parts) for e in ENDPOINTS for k in ('tp','predicted','reference')):
                    raise ValueError('Species and pooled endpoint counts must reconstruct chromosomes')
            arms[name]['background'] = row['background']
        if any(arms['informative'][e]['reference'] != arms['uniform'][e]['reference'] for e in ENDPOINTS):
            raise ValueError('Keep common retained reference denominators')
        differences = {e:dict(
            counts={k:arms['informative'][e][k]-arms['uniform'][e][k] for k in ('tp','predicted','reference')},
            PRF_percentage_points={k:100*(arms['informative'][e][k]-arms['uniform'][e][k]) for k in METRICS}) for e in ENDPOINTS}
        scopes[scope] = dict(arms=arms,informative_minus_uniform=differences)
    for scope in ('macro_species','worst_species'):
        arms = {}
        for name,report in reports.items():
            arms[name] = {e:{k:(sum(report['per_species'][s][e][k] for s in SPECIES)/2 if scope=='macro_species'
                else min(report['per_species'][s][e][k] for s in SPECIES)) for k in METRICS} for e in ENDPOINTS}
        scopes[scope] = dict(arms=arms,informative_minus_uniform_PRF_percentage_points={e:{k:
            100*(arms['informative'][e][k]-arms['uniform'][e][k]) for k in METRICS} for e in ENDPOINTS})
    if any(scopes['pooled']['arms'][name]['exact_chain']['reference'] != 22688 for name in reports):
        raise ValueError('Retain all22688 primary TRAIN chains')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=U.A.PROTOCOL,
        source_jobs=dict(informative=informative_evaluation_job,uniform=evaluation_job),
        source_fit_jobs=dict(informative=inputs['informative_fit_job'],uniform=inputs['uniform_fit_job']),
        scopes=scopes,all_paired_caller_inputs_and_normalized_original_DNA_directly_equal=True,
        DNA_comparison_before_uniform_reference_access=True,counts_and_PRF_reconciled=True,
        reference_policy=reports['uniform']['reference_policy'],optimizer_steps_this_run=0,
        native_forward=False,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began,
        limitation='single-fit additional-learning contrast from already informative-trained source4608; not full-training coarse necessity, independent generalization or model selection')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False);handle.write('\n')
    lines = ['# Original TRAIN coarse information additional learning comparison','',
        '| Scope | Informative chain F1 % | Uniform chain F1 % | Difference pp |',
        '| --- | ---: | ---: | ---: |']
    for scope in (*SPECIES,'pooled'):
        r = scopes[scope]
        values = [100*r['arms'][name]['exact_chain']['f1'] for name in ('informative','uniform')]
        lines.append('| '+scope+' | '+' | '.join(f'{v:.4f}' for v in values)+
                     f" | {r['informative_minus_uniform']['exact_chain']['PRF_percentage_points']['f1']:+.4f} |")
    lines.extend(['',report['limitation']])
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(completed=True,job_id=report['job_id'],source_jobs=report['source_jobs'])),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation-job',required=True)
    parser.add_argument('--informative-evaluation-job',required=True)
    args = parser.parse_args()
    main(args.evaluation_job,args.informative_evaluation_job)
