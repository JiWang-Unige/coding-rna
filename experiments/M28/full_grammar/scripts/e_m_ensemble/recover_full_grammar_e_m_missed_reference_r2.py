#!/usr/bin/env python3
"""Finish two original TRAIN scopes and reconcile four completed timeout scopes."""
import argparse
import json
import os
import time

import diagnose_full_grammar_e_m_missed_reference as D
import readout_full_grammar_gene_posterior_train_r1 as R1

P = D.P
ROOT = P.ROOT
OLD = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE-R1'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE-R2'
CHECKS_JOB = '13549411'
FAILED_JOB = '13554511'
PARENT_JOB = '13554510'
REUSED = tuple(q for _,q,_ in P.E.I.SCOPE[:4])


def guard():
    started = P.read(OLD/'started.json')
    checked = P.read(D.CHECKS/CHECKS_JOB/'summary.json')
    if (started['job_id'] != FAILED_JOB or started['protocol'] != D.PROTOCOL
            or started['checks_job'] != CHECKS_JOB or started['source_posterior_evaluation_job'] != PARENT_JOB
            or started['source_profile_job'] != '13449803' or started['source_evaluation_jobs'] != D.EVALUATIONS
            or any(started[k] for k in ('checkpoint_read','native_forward','predictions_modified',
                'optimizer_steps_this_run','DEV_used','sealed_data_used'))
            or not checked['completed'] or checked['job_id'] != CHECKS_JOB
            or checked['tests_run'] != 4 or checked['failures'] or checked['errors']):
        raise ValueError('Exact unchanged source/qualification of the preserved timeout attempt required')
    parent = R1.OUT/'evaluation'
    if ((parent/'STATUS').read_text().strip() != 'COMPLETED'
            or (parent/'JOBID').read_text().strip() != PARENT_JOB or not (parent/'summary.json').is_file()):
        raise ValueError('Completed original full E/M posterior evaluation required')
    return D.sources()


def reference_keys(seqid):
    refs = []
    with (D.D.PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['seqid'] == seqid: refs.append(row)
    return {D.R.key(dict(r,reference_CDS=[(a,b) for a,b,_ in r['CDS']])) for r in refs}, len(refs)


def read_block(folder, seqid, name, expected, evaluations):
    summary = P.read(folder/'summary.json')
    with (folder/'reference_energy.jsonl').open() as handle:
        rows = [json.loads(line) for line in handle]
    keys = {D.R.key(r) for r in rows}
    rk, n = reference_keys(seqid)
    measured = D.counts(rows)
    if (n != expected or len(rk) != n or len(rows) != expected or len(keys) != len(rows) or keys != rk
            or any(r['seqid'] != seqid for r in rows) or measured != summary['counts']
            or measured['exact_primary_hits'] != evaluations[name]['chromosomes'][seqid]['exact_chain']['tp']
            or len(summary['scopes']) != 2 or {s['strand'] for s in summary['scopes']} != {'+','-'}
            or any(not s['selected_set_gain_replayed'] for s in summary['scopes'])
            or measured['replacement_improvement_counterexamples']):
        raise ValueError('Only complete exact original primary rows and replayed non-counterexample scopes can be reused')
    if D.T.summarize(rows) != summary['role_terms']:
        raise ValueError('Saved complete role-term summary must reconcile with its rows')
    return rows


def validate_reuse():
    began=time.perf_counter(); profile,evaluations=guard()
    out=OUT/'reuse-checks'/os.environ['SLURM_JOB_ID']; verified={}
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=D.PROTOCOL,
        source_failed_job=FAILED_JOB,reused_seqids=REUSED,source_posterior_evaluation_job=PARENT_JOB))
    for species,seqid,length in P.E.I.SCOPE[:4]:
        blocks={name:read_block(OLD/seqid/name,seqid,name,
            profile['chromosomes'][seqid]['exact_chain']['reference'],evaluations) for name in P.MODELS}
        verified[seqid]=D.paired_rows(blocks['E_RANK4'],blocks['M_DENSE'])
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=D.PROTOCOL,
        source_failed_job=FAILED_JOB,source_posterior_evaluation_job=PARENT_JOB,source_profile_job='13449803',
        reused_seqids=REUSED,verified_scopes=verified,all_eight_complete_blocks_reconciled=True,
        incomplete_fifth_scope_not_reused=True,research_DNA_read=False,role_arrays_read=False,
        checkpoint_read=False,native_forward=False,predictions_modified=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,scientific_threshold_changed=False,seconds=time.perf_counter()-began))


def reuse_checked(job):
    folder=OUT/'reuse-checks'/job; report=P.read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip()!='COMPLETED' or not report['completed']
            or report['job_id']!=job or report['protocol']!=D.PROTOCOL
            or report['source_failed_job']!=FAILED_JOB or report['source_posterior_evaluation_job']!=PARENT_JOB
            or tuple(report['reused_seqids'])!=REUSED or not report['all_eight_complete_blocks_reconciled']):
        raise ValueError('Completed actual four-scope partial-output reconciliation is required')


def component(task, reuse_check_job):
    reuse_checked(reuse_check_job)
    began = time.perf_counter(); profile,evaluations = guard()
    species,seqid,length = P.E.I.SCOPE[task]; out = OUT/seqid
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=D.PROTOCOL,
        task=task,species=species,seqid=seqid,length=length,source_failed_job=FAILED_JOB,
        checks_job=CHECKS_JOB,reuse_check_job=reuse_check_job,source_posterior_evaluation_job=PARENT_JOB,source_profile_job='13449803',
        original_scientific_computation_unchanged=True,DEV_used=False,sealed_data_used=False))
    chromosome, blocks = D.diagnose_chromosome(species,seqid,length,
        profile['chromosomes'][seqid]['exact_chain']['reference'],evaluations,OUT)
    counterexamples = sum(D.counts(rows)['replacement_improvement_counterexamples'] for rows in blocks.values())
    P.write(out/'summary.json',dict(completed=not counterexamples,job_id=os.environ['SLURM_JOB_ID'],
        protocol=D.PROTOCOL,task=task,species=species,seqid=seqid,length=length,source_failed_job=FAILED_JOB,
        checks_job=CHECKS_JOB,reuse_check_job=reuse_check_job,source_posterior_evaluation_job=PARENT_JOB,source_profile_job='13449803',
        chromosome=chromosome,original_scientific_computation_unchanged=True,
        unexpected_improving_single_reference_replacements=counterexamples,
        checkpoint_read=False,native_forward=False,predictions_modified=False,decoder_rerun=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,scientific_threshold_changed=False,
        seconds=time.perf_counter()-began))
    if counterexamples: raise ValueError('Improving replacement requires MAP/replay investigation, not recovery interpretation')


def aggregate(jobs, reuse_check_job):
    reuse_checked(reuse_check_job)
    began=time.perf_counter(); profile,evaluations=guard()
    if len(jobs)!=2 or len(set(jobs))!=2: raise ValueError('Two distinct actual remaining-scope job IDs required')
    for task,job in zip((4,5),jobs):
        seqid=P.E.I.SCOPE[task][1]; folder=OUT/seqid; report=P.read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip()!='COMPLETED' or (folder/'JOBID').read_text().strip()!=job
                or not report['completed'] or report['job_id']!=job or report['task']!=task
                or report['seqid']!=seqid or report['protocol']!=D.PROTOCOL
                or report['source_failed_job']!=FAILED_JOB or report['reuse_check_job']!=reuse_check_job
                or not report['original_scientific_computation_unchanged']
                or report['unexpected_improving_single_reference_replacements']):
            raise ValueError('Both actual recovery computations must complete before aggregation')
    out=OUT/'evaluation'
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=D.PROTOCOL,
        component_jobs=jobs,reuse_check_job=reuse_check_job,source_failed_job=FAILED_JOB,reused_seqids=REUSED,
        source_posterior_evaluation_job=PARENT_JOB,source_profile_job='13449803'))
    all_rows={name:[] for name in P.MODELS}; chromosomes={}
    for species,seqid,length in P.E.I.SCOPE:
        root=OLD if seqid in REUSED else OUT
        blocks={name:read_block(root/seqid/name,seqid,name,
            profile['chromosomes'][seqid]['exact_chain']['reference'],evaluations) for name in P.MODELS}
        for name,rows in blocks.items(): all_rows[name].extend(rows)
        chromosomes[seqid]=dict(species=species,length=length,
            models={name:D.counts(rows) for name,rows in blocks.items()},
            paired=D.paired_rows(blocks['E_RANK4'],blocks['M_DENSE']))
    scopes={}
    for scope in ('arabidopsis_thaliana','oryza_sativa','pooled'):
        blocks={name:[r for r in rows if scope=='pooled' or r['species']==scope] for name,rows in all_rows.items()}
        scopes[scope]=dict(models={name:dict(counts=D.counts(rows),role_terms=D.T.summarize(rows))
            for name,rows in blocks.items()},paired=D.paired_rows(blocks['E_RANK4'],blocks['M_DENSE']))
        for name,report in evaluations.items():
            expected=(report['pooled'] if scope=='pooled' else report['per_species'][scope])['exact_chain']
            measured=scopes[scope]['models'][name]['counts']
            if (measured['reference'],measured['exact_primary_hits'])!=(expected['reference'],expected['tp']):
                raise ValueError('Entire retained reference and exact endpoints must match original E/M results')
    if any(len(rows)!=22688 for rows in all_rows.values()): raise ValueError('Retain all22688 original primary references')
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=D.PROTOCOL,
        checks_job=CHECKS_JOB,reuse_check_job=reuse_check_job,component_jobs=jobs,source_failed_job=FAILED_JOB,
        source_failed_job_terminal='TIMEOUT_batch_CANCELLED_0:15_extern_COMPLETED_0:0',
        reused_completed_scope_folders={seqid:str(OLD/seqid) for seqid in REUSED},
        source_posterior_evaluation_job=PARENT_JOB,source_profile_job='13449803',source_evaluation_jobs=D.EVALUATIONS,
        chromosomes=chromosomes,scopes=scopes,all_original_primary_references_retained=True,
        fixed_grammar_and_reference_role_rows_shared=True,original_exact_chain_counts_replayed=True,
        all_reused_reference_rows_and_role_terms_reconciled=True,original_scientific_computation_unchanged=True,
        unexpected_improving_single_reference_replacements=0,checkpoint_read=False,native_forward=False,
        real_reference_GFF_read=False,reference_profile_rows_read=True,reference_used_for_posthoc_scoring_only=True,
        predictions_modified=False,decoder_rerun=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        scientific_threshold_changed=False,numerical_near_zero_tolerance=D.D.EPSILON,
        seconds=time.perf_counter()-began,
        interpretation='complete_original_TRAIN_score_diagnosis_with_preserved_timeout_not_capacity_LR_causal_or_generalization'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('validate-reuse','component','aggregate'))
    parser.add_argument('--task',type=int,choices=(4,5))
    parser.add_argument('--component-jobs',nargs=2)
    parser.add_argument('--reuse-check-job')
    args=parser.parse_args()
    D.torch.set_num_threads(2)
    if args.mode=='validate-reuse': validate_reuse()
    elif args.mode=='component': component(args.task,args.reuse_check_job)
    else: aggregate(args.component_jobs,args.reuse_check_job)
