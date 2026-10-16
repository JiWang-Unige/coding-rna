#!/usr/bin/env python3
"""Immutable equal-energy-pool whole-TRAIN event probabilities and fixed diagnostics."""
import argparse
import math
import os
import time
import unittest
from unittest.mock import patch

import numpy as np
import torch
import readout_full_grammar_gene_posterior_train as P
import readout_full_grammar_gene_posterior_train_r1 as R1
import run_full_grammar_e_m_energy_pool_train as Q
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_energy_pool import equal_role_energy_pool
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries

ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-TRAIN'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_energy_pool_posterior_train.md'
RAW_JOBS = ('13551143','13551145','13551146','13551147','13551148','13551149')


def sources():
    numerical = R1.qualified_native()
    Q.qualification('13551105')
    folder = Q.OUT/'evaluation'; pooled = P.read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not pooled['completed']
            or pooled['job_id'] != '13551150' or pooled['identity'] != Q.IDENTITY
            or tuple(pooled['raw_jobs']) != RAW_JOBS or not pooled['all_original_primary_references_retained']
            or not pooled['unsupported_reference_strata_retained']
            or any(pooled[k] for k in ('original_predictions_modified','weights_modified',
                'reference_used_in_inference','optimizer_steps_this_run','DEV_used','sealed_data_used',
                'threshold_selected','frozen_DEV_family_changed'))):
        raise ValueError('The completed unchanged fixed E/M pool and all primary references are required')
    native_folder = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-RAW-NATIVE'
    native = P.read(native_folder/'summary.json')
    if ((native_folder/'STATUS').read_text().strip() != 'COMPLETED' or not native['completed']
            or native['job_id'] != '13556298' or native['identity'] != Q.IDENTITY
            or not native['all_six_native_raw_outputs_match_independent_mean_energy_decode']):
        raise ValueError('Completed actual raw-native parity of the identical pool rule is required')
    return pooled, numerical


def checked(job):
    folder = OUT/'checks'/job; report = P.read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
            or report['job_id'] != job or report['protocol'] != PROTOCOL or report['identity'] != Q.IDENTITY
            or report['tests_run'] != 2 or report['failures'] or report['errors']
            or not report['all_twelve_emission_sources_and_six_pool_outputs_bound']):
        raise ValueError('New completed FP64 mean-energy/event qualification is required')


def checks():
    import test_m28_energy_pool_gene_posterior as T
    began = time.perf_counter(); out = OUT/'checks'/os.environ['SLURM_JOB_ID']
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,identity=Q.IDENTITY,
        research_DNA_read=False,reference_used=False,role_arrays_read=False,checkpoint_read=False,native_forward=False))
    tested = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(T))
    if not tested.wasSuccessful():
        raise RuntimeError('Pooled FP64 gene-event or complete-set accounting failed')
    pooled, numerical = sources()
    for name in P.MODELS:
        _,_,reports = P.source(name)
        for task,(_,seqid,_) in enumerate(P.E.I.SCOPE):
            Q.emissions(name,task,reports[seqid])
    for task,(_,seqid,length) in enumerate(P.E.I.SCOPE):
        folder = Q.OUT/seqid; raw = P.read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not raw['completed']
                or (folder/'JOBID').read_text().strip() != RAW_JOBS[task] or raw['job_id'] != RAW_JOBS[task]
                or raw['identity'] != Q.IDENTITY or (raw['seqid'],raw['length']) != (seqid,length)
                or raw['primary_CDS_chain_predictions'] != pooled['energy_pool']['chromosomes'][seqid]['exact_chain']['predicted']):
            raise ValueError('All six exact completed raw-pool outputs must remain bound')
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        identity=Q.IDENTITY,tests_run=tested.testsRun,failures=len(tested.failures),errors=len(tested.errors),
        all_twelve_emission_sources_and_six_pool_outputs_bound=True,source_pool_evaluation_job='13551150',
        source_raw_native_parity_job='13556298',numerical_candidate=numerical['candidate'],
        original_probability_and_mass_tolerances_unchanged=True,research_DNA_read=False,reference_used=False,
        role_arrays_read=False,checkpoint_read=False,native_forward=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began))


def predictions(seqid, length, expected_count):
    path = Q.OUT/seqid/'full_grammar.gff3'
    txs = P.C.S.E.primary_transcripts(P.C.S.E.parse_annotation(path,{seqid:length}),complete_only=False)
    scores = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip(): continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] == 'mRNA': scores[P.C.S.E.parse_attrs(fields[8])['ID']] = float(fields[5])
    if (len(txs) != expected_count or len({P.C.S.key(t) for t in txs}) != len(txs)
            or set(scores) != {t['id'] for t in txs} or any(not math.isfinite(v) for v in scores.values())):
        raise ValueError('Preserve every unique immutable pool query and finite eight-digit GFF gain')
    for tx in txs: tx['score'] = scores[tx['id']]
    return txs


def posterior(task, checks_job):
    checked(checks_job); pooled,numerical = sources()
    species,seqid,length = P.E.I.SCOPE[task]; out = OUT/seqid/'posterior'
    began = time.perf_counter()
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,identity=Q.IDENTITY,
        task=task,species=species,seqid=seqid,length=length,checks_job=checks_job,
        source_pool_raw_job=RAW_JOBS[task],source_emission_jobs={n:j[task] for n,j in Q.EMISSION_JOBS.items()},
        reference_used=False,DEV_used=False,sealed_data_used=False,native_forward=False))
    folders,sequences = {},{}
    for name in P.MODELS:
        module,_,raw = P.source(name)
        folders[name] = Q.emissions(name,task,raw[seqid])
        sequences[name] = P.C.M.read_allowed(module.OUT/seqid/'input.fasta',{seqid:length})[seqid]
    if sequences['E_RANK4'] != sequences['M_DENSE']:
        raise ValueError('Literal identical original whole-TRAIN DNA is required for pooling')
    sequence = sequences['E_RANK4']; del sequences
    txs = predictions(seqid,length,pooled['energy_pool']['chromosomes'][seqid]['exact_chain']['predicted'])
    rows,records = [],[]
    for strand,filename in (('+','plus.npy'),('-','minus.npy')):
        ordered = sorted((t for t in txs if t['strand']==strand),
            key=lambda t:t['CDS'][0][0] if strand=='+' else -t['CDS'][-1][1])
        chains = [P.oriented_chain([(a,b) for a,b,_ in t['CDS']],length,strand) for t in ordered]
        dna = sequence if strand=='+' else P.C.M.rc(sequence)
        arrays = [np.load(folders[name]/filename,mmap_mode='c',allow_pickle=False) for name in P.MODELS]
        if any(a.shape!=(length,15) or a.dtype!=np.float32 for a in arrays):
            raise ValueError('Aligned original full-chromosome FP32 fifteen-role sources required')
        values = equal_role_energy_pool(*(torch.from_numpy(a) for a in arrays)); del arrays
        started = time.perf_counter()
        with patch.multiple(G,_prefix_queries=_scaled_prefix_queries,_suffix_queries=_scaled_suffix_queries):
            result = G.complete_gene_set_posteriors(dna,values,chains,128)
        elapsed = time.perf_counter()-started; del values,dna
        if len(result['chain_gains']) != len(ordered) or len(result['gene_log_probabilities']) != len(ordered):
            raise ValueError('Every existing pool query must remain present')
        for tx,chain,gain,logp in zip(ordered,chains,result['chain_gains'],result['gene_log_probabilities']):
            if not math.isclose(gain,float(tx['score']),rel_tol=6e-8,abs_tol=1e-5):
                raise ValueError('Pooled gene-event readout must preserve the immutable MAP GFF gain')
            rows.append(dict(id=tx['id'],strand=strand,genomic_CDS=[(a,b) for a,b,_ in tx['CDS']],
                oriented_CDS=chain,gain=gain,original_rounded_GFF_gain=float(tx['score']),log_probability=logp))
        records.append(dict(strand=strand,supplied_bases=length,queried_genes=len(chains),
            **{k:result[k] for k in ('log_partition','reverse_log_partition','queried_set_log_probability',
                'empty_set_log_probability','state_count','stored_boundary_scalars')},posterior_seconds=elapsed))
        print(P.json.dumps(dict(stage='whole_pool_gene_posterior',seqid=seqid,strand=strand,seconds=elapsed)),flush=True)
        del result
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        identity=Q.IDENTITY,task=task,species=species,seqid=seqid,length=length,checks_job=checks_job,
        source_pool_raw_job=RAW_JOBS[task],source_emission_jobs={n:j[task] for n,j in Q.EMISSION_JOBS.items()},
        numerical_candidate=numerical['candidate'],numerical_checks_job=R1.NUMERICAL_CHECKS_JOB,
        genes=rows,records=records,source_DNA_equality_before_pooling=True,
        full_chromosome_conditional_posterior_measured=True,denominator_all_legal_gene_sets_not_selected_pool=True,
        original_MAP_scores_and_exported_chains_preserved=True,original_probability_and_mass_tolerances_unchanged=True,
        production_algorithm_modified=False,predictions_modified=False,weights_modified=False,pooled_role_arrays_saved=False,
        reference_used=False,labels_or_feature_cache_used=False,native_forward=False,checkpoint_read=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,threshold_selected=False,
        fitted_temperature=False,frozen_DEV_family_changed=False,seconds=time.perf_counter()-began,
        interpretation='model_structural_inclusion_not_empirically_calibrated_biological_correctness'))


def diagnostics(rows, reference_count):
    result = P.diagnostics(rows,reference_count)
    # The full tie-group curve is used to calculate every metric, not duplicated
    # at three aggregation levels. Queries/logs remain in the six immutable readouts.
    for label in result.values():
        for rank in ('posterior','internal_gain'):
            del label[rank]['grouped_precision_recall_curve']
    return result


def evaluate(checks_job, jobs):
    checked(checks_job); pooled,numerical = sources(); began=time.perf_counter()
    if len(jobs)!=6 or len(set(jobs))!=6: raise ValueError('Six distinct actual pool posterior job IDs required')
    raw = {}
    for task,(species,seqid,length) in enumerate(P.E.I.SCOPE):
        folder=OUT/seqid/'posterior'; r=P.read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip()!='COMPLETED' or (folder/'JOBID').read_text().strip()!=jobs[task]
                or not r['completed'] or r['job_id']!=jobs[task] or r['protocol']!=PROTOCOL
                or r['identity']!=Q.IDENTITY or r['checks_job']!=checks_job
                or (r['task'],r['species'],r['seqid'],r['length'])!=(task,species,seqid,length)
                or r['numerical_candidate']!=numerical['candidate'] or len(r['records'])!=2
                or not r['full_chromosome_conditional_posterior_measured']
                or not r['original_MAP_scores_and_exported_chains_preserved']
                or not r['original_probability_and_mass_tolerances_unchanged']
                or any(r[k] for k in ('reference_used','predictions_modified','weights_modified','DEV_used',
                    'sealed_data_used','optimizer_steps_this_run','threshold_selected','fitted_temperature'))):
            raise ValueError('All six completed immutable full-strand pool posteriors required before reference access')
        raw[seqid]=r
    parent_folder=R1.OUT/'evaluation'; parent=P.read(parent_folder/'summary.json')
    if ((parent_folder/'STATUS').read_text().strip()!='COMPLETED' or not parent['completed']
            or parent['job_id']!='13554510' or parent['protocol']!=P.PROTOCOL
            or not parent['all_twelve_posteriors_verified_before_reference_access']
            or not parent['original_exact_chain_counts_replayed'] or parent['DEV_used'] or parent['sealed_data_used']):
        raise ValueError('Completed fixed E/M R1 ranking comparison is required')
    out=OUT/'evaluation'
    P.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,identity=Q.IDENTITY,
        checks_job=checks_job,posterior_jobs=jobs,source_parent_evaluation_job='13554510',
        all_six_pool_posteriors_complete_before_reference_access=True))
    chromosomes,by_species,reference_counts={},{},{}
    for species,seqid,length in P.E.I.SCOPE:
        sequence=P.C.M.read_allowed(P.E.OUT/seqid/'input.fasta',{seqid:length})
        data=ROOT/'data/m1_screen'/species
        coding=P.C.S.E.parse_annotation(data/'reference.gff3',{seqid:length},protein_coding_only=True)
        with (data/'reference.gff3').open() as handle: _,_,raw_cds,_=P.C.P.raw_features(handle,seqid)
        refs,isoforms=P.C.select_reference(coding,raw_cds,sequence)
        if any(not P.C.R7.phase_chain_consistent(t) for t in refs): raise ValueError('Retained reference phase mismatch')
        primary,all_iso=({P.C.S.key(t) for t in group} for group in (refs,isoforms))
        rows=raw[seqid]['genes']
        for row in rows:
            key=P.C.S.key(dict(seqid=seqid,strand=row['strand'],CDS=[(a,b,'0') for a,b in row['genomic_CDS']]))
            row['exact_primary']=key in primary; row['exact_any_assessable_isoform']=key in all_iso
        old=pooled['energy_pool']['chromosomes'][seqid]['exact_chain']
        if ((len(rows),sum(r['exact_primary'] for r in rows),len(primary))!=(old['predicted'],old['tp'],old['reference'])
                or len(primary)!=parent['retained_primary_reference_counts'][seqid]):
            raise ValueError('Full immutable pool endpoints and parent reference denominators must replay')
        reference_counts[seqid]=len(primary); chromosomes[seqid]=diagnostics(rows,len(primary))
        by_species.setdefault(species,[]).extend(rows)
    per_species={species:diagnostics(rows,sum(reference_counts[q] for s,q,_ in P.E.I.SCOPE if s==species))
                 for species,rows in by_species.items()}
    all_rows=[r for group in by_species.values() for r in group]
    if sum(reference_counts.values())!=22688: raise ValueError('Retain all22688 original primary references')
    combined=diagnostics(all_rows,22688)
    comparisons={name:{label:dict(
        prediction_pool_AP_difference=combined[label]['posterior']['prediction_pool_average_precision']-model['pooled'][label]['posterior']['prediction_pool_average_precision'],
        prediction_pool_AUROC_difference=combined[label]['posterior']['prediction_pool_AUROC']-model['pooled'][label]['posterior']['prediction_pool_AUROC'],
        reference_relative_Brier_difference=combined[label]['probability']['reference_relative_Brier']-model['pooled'][label]['probability']['reference_relative_Brier'])
        for label in combined} for name,model in parent['models'].items()}
    P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        identity=Q.IDENTITY,checks_job=checks_job,posterior_jobs=jobs,source_parent_evaluation_job='13554510',
        energy_pool=dict(chromosomes=chromosomes,per_species=per_species,pooled=combined),pool_minus_parent=comparisons,
        retained_primary_reference_counts=reference_counts,reference_policy=pooled['reference_policy'],
        all_six_pool_posteriors_complete_before_reference_access=True,original_exact_chain_counts_replayed=True,
        all_original_primary_references_retained=True,full_tie_curves_used_in_AP_AUC_and_fixed_fractions=True,
        full_grouped_curves_duplicated_in_summary=False,reference_used_only_after_readout=True,
        predictions_modified=False,weights_modified=False,native_forward=False,checkpoint_read=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,threshold_selected=False,
        calibration_fitted=False,frozen_DEV_family_changed=False,seconds=time.perf_counter()-began,
        interpretation='post_outcome_TRAIN_selected_pool_reference_relative_ranking_not_independent_detection_AP_or_calibration'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('checks','posterior','evaluate'))
    parser.add_argument('--task',type=int,choices=range(6))
    parser.add_argument('--checks-job')
    parser.add_argument('--posterior-jobs',nargs=6)
    args=parser.parse_args(); torch.set_num_threads(2)
    if args.mode=='checks': checks()
    elif args.mode=='posterior': posterior(args.task,args.checks_job)
    else: evaluate(args.checks_job,args.posterior_jobs)
