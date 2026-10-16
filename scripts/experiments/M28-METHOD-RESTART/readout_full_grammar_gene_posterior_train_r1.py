#!/usr/bin/env python3
"""Same immutable TRAIN posterior/ruler with qualified numerical kernels only."""
import argparse
import math
import os
import time
from unittest.mock import patch

import numpy as np
import torch

import readout_full_grammar_gene_posterior_train as P
import run_full_grammar_e_m_energy_pool_train as Q
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries

OUT = P.ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-R1'
CHECKS_JOB = '13548771'  # Unchanged completed coordinate/ranking/source qualification.
NUMERICAL_CHECKS_JOB = '13552225'
NATIVE_JOB = '13552583'


def qualified_native():
    P.qualification(CHECKS_JOB)
    checks = P.ROOT/'outputs/M28-FULL-GRAMMAR-SCALED-MASS-CHECKS'/NUMERICAL_CHECKS_JOB
    c = P.read(checks/'summary.json')
    if ((checks/'STATUS').read_text().strip() != 'COMPLETED' or not c['completed']
            or c['job_id'] != NUMERICAL_CHECKS_JOB or c['tests_run'] != 8 or c['failures'] or c['errors']
            or not c['original_probability_and_mass_tolerances_unchanged'] or c['production_algorithm_modified']):
        raise ValueError('Completed unchanged-tolerance eight-test numerical qualification required')
    folder = P.ROOT/'outputs/M28-FULL-GRAMMAR-SCALED-MASS-NATIVE'
    n = P.read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not n['completed']
            or n['job_id'] != NATIVE_JOB or n['checks_job'] != NUMERICAL_CHECKS_JOB
            or n['protocol'] != P.PROTOCOL or n['candidate'] != c['candidate']
            or (n['model'],n['task'],n['strand'],n['seqid'],n['length'])
                != ('E_RANK4',0,'+','NC_003070.9',30427671)
            or n['source_emissions_job'] != '13548914' or n['source_fit_job'] != '13499817'
            or not n['original_forward_backward_check_passed']
            or not n['all_queried_genes_literal_DNA_legal'] or not n['original_MAP_GFF_gains_preserved']
            or not n['original_probability_and_mass_tolerances_unchanged']
            or not n['denominator_all_legal_gene_sets_not_selected_pool']
            or any(n[k] for k in ('production_algorithm_modified','reference_used','checkpoint_read',
                'native_forward','predictions_modified','weights_modified','optimizer_steps_this_run','DEV_used','sealed_data_used'))):
        raise ValueError('Exact completed original E0 plus whole-chromosome qualification required')
    return n


def checks():
    began = time.perf_counter(); native = qualified_native()
    for name in P.MODELS:
        module,_,reports = P.source(name)
        for task,(_,seqid,_) in enumerate(P.E.I.SCOPE):
            Q.emissions(name,task,reports[seqid])
    module,_,reports = P.source('E_RANK4'); _,seqid,length = P.E.I.SCOPE[0]
    ordered = sorted((t for t in P.predictions(module,seqid,length,reports[seqid]) if t['strand']=='+'),
        key=lambda t:t['CDS'][0][0])
    if [t['id'] for t in ordered] != [r['id'] for r in native['genes']]:
        raise ValueError('Actual native reuse must retain every original E0 plus query in order')
    for tx,row in zip(ordered,native['genes']):
        if (not math.isfinite(row['log_probability']) or row['log_probability']>0
                or not math.isclose(row['gain'],tx['score'],rel_tol=6e-8,abs_tol=1e-5)):
            raise ValueError('Actual native reuse must preserve finite supported logs and original GFF gains')
    out = P.ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-R1-BINDING'/os.environ['SLURM_JOB_ID']
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        numerical_checks_job=NUMERICAL_CHECKS_JOB,native_numerical_qualification_job=NATIVE_JOB,
        all_twelve_completed_emission_bindings_verified=True,actual_E0_plus_reuse_queries_and_gains_verified=True,
        research_DNA_read=False,reference_used=False,role_arrays_read=False,checkpoint_read=False,native_forward=False,
        production_algorithm_modified=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began))


def posterior(name, task):
    native = qualified_native(); module,fit_job,reports = P.source(name)
    species,seqid,length = P.E.I.SCOPE[task]
    emitted = Q.emissions(name,task,reports[seqid]); out = OUT/name/seqid/'posterior'
    began = time.perf_counter()
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        model=name,task=task,species=species,seqid=seqid,length=length,checks_job=CHECKS_JOB,
        numerical_checks_job=NUMERICAL_CHECKS_JOB,native_numerical_qualification_job=NATIVE_JOB,
        source_emissions_job=Q.EMISSION_JOBS[name][task],source_fit_job=fit_job,
        reference_used=False,DEV_used=False,sealed_data_used=False,native_forward=False))
    sequence = P.C.M.read_allowed(module.OUT/seqid/'input.fasta',{seqid:length})[seqid]
    txs = P.predictions(module,seqid,length,reports[seqid]); rows,records = [],[]
    for strand,filename in (('+','plus.npy'),('-','minus.npy')):
        ordered = sorted((t for t in txs if t['strand']==strand),
            key=lambda t:t['CDS'][0][0] if strand=='+' else -t['CDS'][-1][1])
        chains = [P.oriented_chain([(a,b) for a,b,_ in t['CDS']],length,strand) for t in ordered]
        reused = name=='E_RANK4' and task==0 and strand=='+'
        started = time.perf_counter()
        if reused:
            if [t['id'] for t in ordered] != [r['id'] for r in native['genes']]:
                raise ValueError('Qualification queries must exactly match every original E0 plus query in order')
            gains = [r['gain'] for r in native['genes']]; logs = [r['log_probability'] for r in native['genes']]
            result = dict(log_partition=native['log_partition'],reverse_log_partition=native['reverse_log_partition'],
                queried_set_log_probability=native['queried_set_log_probability'],
                empty_set_log_probability=native['empty_set_log_probability'],state_count=native['state_count'],
                stored_boundary_scalars=len({0,*(c[0][0] for c in chains)})+len({length,*(c[-1][1] for c in chains)}))
        else:
            dna = sequence if strand=='+' else P.C.M.rc(sequence)
            array = np.load(emitted/filename,mmap_mode='c',allow_pickle=False)
            if array.shape != (length,15) or array.dtype != np.float32:
                raise ValueError('Original full chromosome FP32 fifteen-role array required')
            # Only this readout process substitutes the qualified numerical messages.
            with patch.multiple(G,_prefix_queries=_scaled_prefix_queries,_suffix_queries=_scaled_suffix_queries):
                result = G.complete_gene_set_posteriors(dna,torch.from_numpy(array),chains,128)
            gains,logs = result['chain_gains'],result['gene_log_probabilities']
            del array,dna
        elapsed = time.perf_counter()-started
        if len(gains)!=len(ordered) or len(logs)!=len(ordered): raise ValueError('Every immutable query retained')
        for tx,chain,gain,logp in zip(ordered,chains,gains,logs):
            if not math.isclose(gain,float(tx['score']),rel_tol=6e-8,abs_tol=1e-5):
                raise ValueError('Original rounded MAP GFF gain must remain unchanged')
            rows.append(dict(id=tx['id'],strand=strand,genomic_CDS=[(a,b) for a,b,_ in tx['CDS']],
                oriented_CDS=chain,gain=gain,original_rounded_GFF_gain=float(tx['score']),log_probability=logp))
        records.append(dict(strand=strand,supplied_bases=length,queried_genes=len(chains),
            **{k:result[k] for k in ('log_partition','reverse_log_partition','queried_set_log_probability',
                'empty_set_log_probability','state_count','stored_boundary_scalars')},
            posterior_seconds=elapsed,native_qualification_reused=reused))
        print(P.json.dumps(dict(stage='whole_chromosome_posterior_R1',model=name,seqid=seqid,strand=strand,
            native_qualification_reused=reused,seconds=elapsed)),flush=True)
        del result
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        model=name,task=task,species=species,seqid=seqid,length=length,checks_job=CHECKS_JOB,
        numerical_checks_job=NUMERICAL_CHECKS_JOB,native_numerical_qualification_job=NATIVE_JOB,
        numerical_candidate=native['candidate'],source_emissions_job=Q.EMISSION_JOBS[name][task],
        source_fit_job=fit_job,checkpoint=str(module.CHECKPOINT),genes=rows,records=records,
        source_DNA_equality_receipt_verified=True,full_chromosome_conditional_posterior_measured=True,
        denominator_all_legal_gene_sets_not_selected_pool=True,original_MAP_scores_and_exported_chains_preserved=True,
        original_probability_and_mass_tolerances_unchanged=True,production_algorithm_modified=False,
        predictions_modified=False,reference_used=False,labels_or_feature_cache_used=False,native_forward=False,
        checkpoint_read=False,DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0,
        posterior_threshold_selected=False,fitted_temperature=False,seconds=time.perf_counter()-began,
        interpretation='model_structural_support_not_empirically_calibrated_biological_correctness'))


def evaluate(jobs):
    native = qualified_native()
    if len(jobs)!=12 or len(set(jobs))!=12: raise ValueError('Twelve distinct actual R1 posterior job IDs required')
    for i,(name,task) in enumerate((n,t) for n in P.MODELS for t in range(6)):
        seqid = P.E.I.SCOPE[task][1]; folder = OUT/name/seqid/'posterior'; report = P.read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip()!='COMPLETED' or (folder/'JOBID').read_text().strip()!=jobs[i]
                or report['job_id']!=jobs[i] or report['numerical_checks_job']!=NUMERICAL_CHECKS_JOB
                or report['native_numerical_qualification_job']!=NATIVE_JOB
                or report['numerical_candidate']!=native['candidate']
                or not report['original_probability_and_mass_tolerances_unchanged']
                or report['production_algorithm_modified']):
            raise ValueError('All exact completed R1 numerical identities required before reference access')
    # Reuse the unchanged ruler, ranking/bins, full-reference denominators and receipt guards.
    P.OUT = OUT
    P.evaluate(CHECKS_JOB)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('checks','posterior','evaluate'))
    parser.add_argument('--model',choices=tuple(P.MODELS))
    parser.add_argument('--task',type=int,choices=range(6))
    parser.add_argument('--posterior-jobs',nargs=12)
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode=='checks': checks()
    elif args.mode=='posterior': posterior(args.model,args.task)
    else: evaluate(args.posterior_jobs)
