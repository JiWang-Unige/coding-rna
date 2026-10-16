#!/usr/bin/env python3
"""Fixed local energy-disagreement penalty on original TRAIN, same common ruler."""
import argparse
import json
import os
import time

import numpy as np
import torch

import annotate_full_grammar_coarse_channel as R
import common_ruler as C
import evaluate_full_grammar_original_train as F
import evaluate_native_train_baselines as N
import readout_full_grammar_gene_posterior_train as P
from src.m28.contig_emission_lattice import context_segments
from src.m28.full_grammar_energy_pool import minimum_centered_role_energy_pool
from src.m28.full_grammar_skyline import decode_full_grammar_skyline

ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-TRAIN'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_minimum_centered_energy_train_plan.md'
EMISSION_JOBS = {
    'E_RANK4': ('13548914','13548916','13548918','13548920','13548922','13548924'),
    'M_DENSE': ('13548926','13548928','13548930','13548933','13548935','13548938')}
IDENTITY = dict(max_exons=128,
    pool='per_base_role_minimum_of_FP64_role_minus_own_branch_background',
    absolute_disagreement_coefficient=0.5, whole_gene_minimum_used=False,
    score_accumulation='FP64_per_base_forward_order', context_policy=R.CONTEXT,
    decoder='exact_exon_budget_skyline', null=0, role_weights=1,
    selected_chain_pool_used=False, probability_calibration_used=False)


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def emissions(name, task, raw):
    module, fit_job = P.MODELS[name]
    species, seqid, length = P.E.I.SCOPE[task]
    folder = P.OUT/name/seqid/'emissions'; report = read(folder/'summary.json')
    expected = [[a,min(a+C.M.WIDTH,length)] for a in C.M.starts(length)]
    segments = [list(x) for x in context_segments(length, expected)]
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED'
            or (folder/'JOBID').read_text().strip() != EMISSION_JOBS[name][task]
            or not report['completed'] or report['job_id'] != EMISSION_JOBS[name][task]
            or report['protocol'] != P.PROTOCOL or report['checks_job'] != '13548771'
            or (report['model'],report['task'],report['species'],report['seqid'],report['length'])
                != (name,task,species,seqid,length)
            or report['source_fit_job'] != fit_job or report['checkpoint'] != str(module.CHECKPOINT)
            or report['windows'] != expected or not report['source_DNA_equality_receipt_verified']
            or not report['whole_original_TRAIN_contig_both_strands']
            or any(report[k] for k in ('reference_used','prediction_GFF_content_read',
                'labels_or_feature_cache_used','DEV_used','sealed_data_used','optimizer_steps_this_run',
                'weights_modified','original_predictions_modified'))
            or len(report['records']) != 2):
        raise ValueError('Exact completed original final E/M whole-chromosome emission source required')
    for record, strand, filename in zip(report['records'], ('+','-'), ('plus.npy','minus.npy')):
        if (record['strand'] != strand or record['filename'] != filename
                or record['shape'] != [length,15] or record['dtype'] != 'torch.float32'
                or record['genomic_context_segments'] != segments
                or any(s['windows'] != expected or s['genomic_context_segments'] != segments
                       for s in raw['scopes'])):
            raise ValueError('Identical oriented row ownership and fifteen-role layout required')
    return folder


def qualification(job):
    folder = CHECKS/job; result = read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or job != '13567792'
            or result['job_id'] != job or not result['completed'] or result['protocol'] != PROTOCOL
            or result['tests_run'] != 3 or result['failures'] or result['errors']
            or not result['all_twelve_completed_emission_sources_verified']
            or result['energy'] != 'minimum_of_FP64_per_branch_role_minus_own_background'
            or result['equivalent_fixed_penalty'] != 'mean_centered_energy_minus_half_absolute_disagreement'
            or any(result[k] for k in ('whole_gene_minimum_used','research_DNA_read','research_reference_read',
                'role_arrays_read','checkpoint_read','native_forward','optimizer_steps_this_run','DEV_used','sealed_data_used'))):
        raise ValueError('Completed fixed local-minimum math/source qualification required')


COMPARATOR_JOBS = dict(E_RANK4='13499823', M_DENSE='13510351', NATIVE='13461307',
    CONSENSUS='13550067', EQUAL_ENERGY_POOL='13551150', CARBON_A='13503419')


def comparator_sources():
    folders = {name:module.EVAL for name,(module,_) in P.MODELS.items()}
    folders.update(NATIVE=N.OUT, CONSENSUS=ROOT/'outputs/M28-FULL-GRAMMAR-E-M-CONSENSUS-TRAIN',
        EQUAL_ENERGY_POOL=ROOT/'outputs/M28-FULL-GRAMMAR-E-M-ENERGY-POOL-TRAIN/evaluation',
        CARBON_A=ROOT/'outputs/M28-CARBON-A-TRAIN-EVAL')
    reports = {}
    for name,folder in folders.items():
        report=read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip()!='COMPLETED' or not report['completed']
                or report['job_id']!=COMPARATOR_JOBS[name] or report['DEV_used'] or report['sealed_data_used']):
            raise ValueError('Exact completed original TRAIN comparator required')
        reports[name]=report
    targets={name:reports[name] for name in P.MODELS}
    targets.update(reports['NATIVE']['methods'], CONSENSUS=reports['CONSENSUS']['consensus'],
        EQUAL_ENERGY_POOL=reports['EQUAL_ENERGY_POOL']['energy_pool'], **reports['CARBON_A']['methods'])
    if set(reports['CARBON_A']['methods']) != {'CarbonA_confidence_ge_0.1_or_NA','CarbonA_unfiltered'}:
        raise ValueError('Both unchanged declared native Carbon-A operating points required')
    reference_policy=reports['NATIVE']['reference_policy']
    if any(r['reference_policy']!=reference_policy for r in reports.values()):
        raise ValueError('All comparator reference policies must match')
    for target in targets.values():
        if (set(target['chromosomes'])!=set(N.EXPECTED) or set(target['per_species'])!=
                {'arabidopsis_thaliana','oryza_sativa'} or target['pooled']['exact_chain']['reference']!=22688):
            raise ValueError('Complete six-scope, two-species and full retained comparator denominator required')
        for seqid in N.EXPECTED:
            if N.reference_metadata(target['chromosomes'][seqid])!=N.reference_metadata(targets['ANNEVO']['chromosomes'][seqid]):
                raise ValueError('Comparator reference strata/base/background population mismatch')
    return targets, reference_policy


def bindings(checks_job):
    qualification(checks_job); began=time.perf_counter()
    targets,reference_policy=comparator_sources()
    out=ROOT/'outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-BINDINGS'/os.environ['SLURM_JOB_ID']
    write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        checks_job=checks_job,identity=IDENTITY,source_comparator_jobs=COMPARATOR_JOBS,
        all_comparator_headers_and_reference_populations_reconciled=True,
        comparator_names=list(targets),reference_policy=reference_policy,
        research_DNA_read=False,research_reference_read=False,role_arrays_read=False,
        checkpoint_read=False,native_forward=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began))


def bound(job):
    folder=ROOT/'outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-BINDINGS'/job
    report=read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip()!='COMPLETED' or not report['completed']
            or report['job_id']!=job or report['protocol']!=PROTOCOL or report['checks_job']!='13567792'
            or report['identity']!=IDENTITY or report['source_comparator_jobs']!=COMPARATOR_JOBS
            or not report['all_comparator_headers_and_reference_populations_reconciled']):
        raise ValueError('Actual completed new driver/comparator interface qualification required')


def raw(task, checks_job, binding_job):
    qualification(checks_job); bound(binding_job)
    began = time.perf_counter(); species, seqid, length = P.E.I.SCOPE[task]; out = OUT/seqid
    folders, sequences = {}, {}
    for name in P.MODELS:
        module, _, reports = P.source(name)
        folders[name] = emissions(name, task, reports[seqid])
        sequences[name] = C.M.read_allowed(module.OUT/seqid/'input.fasta', {seqid:length})[seqid]
    if sequences['E_RANK4'] != sequences['M_DENSE']:
        raise ValueError('Source normalized DNA must be directly identical before pooling')
    sequence = sequences['E_RANK4']; del sequences
    write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        task=task, species=species, seqid=seqid, length=length, checks_job=checks_job, binding_job=binding_job, identity=IDENTITY,
        source_emission_jobs={name:jobs[task] for name,jobs in EMISSION_JOBS.items()},
        reference_used=False, DEV_used=False, sealed_data_used=False, native_forward=False))
    transcripts, scopes = [], []
    for strand, filename in (('+','plus.npy'), ('-','minus.npy')):
        started = time.perf_counter()
        arrays = [np.load(folders[name]/filename, mmap_mode='c', allow_pickle=False) for name in P.MODELS]
        if any(a.shape != (length,15) or a.dtype != np.float32 for a in arrays):
            raise ValueError('Full aligned chromosome FP32 fifteen-role source arrays required')
        values = minimum_centered_role_energy_pool(*(torch.from_numpy(a) for a in arrays))
        pooling_seconds = time.perf_counter()-started; del arrays
        dna = sequence if strand == '+' else C.M.rc(sequence)
        started = time.perf_counter()
        result = decode_full_grammar_skyline(dna, values, 128)
        decode_seconds = time.perf_counter()-started
        for i,(chain,gain) in enumerate(zip(result['chains'],result['chain_gains'])):
            genomic = R.genomic_chain(chain,0,length,strand)
            tx = R.as_transcript(genomic,seqid,strand,'M28_MINROLE_'+str(task)+'_'+strand.replace('+','p').replace('-','m')+'_'+str(i),gain)
            failures,_ = C.M.chain_checks(tx,sequence)
            if failures or not 1 <= len(tx['CDS']) <= 128:
                raise ValueError('Pooled raw chain DNA, phase or cap failure')
            transcripts.append(tx)
        scopes.append(dict(strand=strand, gain=result['gain'], selected_complete_CDS_chains=len(result['chains']),
            counts=result['counts'], tie_policy=result['tie_policy'], pruning=result['pruning'],
            score_accumulation=result['score_accumulation'], pooling_seconds=pooling_seconds,
            decode_seconds=decode_seconds, whole_oriented_contig=True, reference_used=False))
        print(json.dumps(dict(stage='pooled_strand_complete',seqid=seqid,strand=strand,
            predictions=len(result['chains']),decode_seconds=decode_seconds)),flush=True)
        del values, result, dna
    R.write_gff3(out/'full_grammar.gff3',transcripts,{seqid:length},'M28_MINROLE')
    artifact_bytes = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    if artifact_bytes > 64*1024**2:
        raise ValueError('Declared small-GFF output budget exceeded')
    write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        task=task,species=species,seqid=seqid,length=length,checks_job=checks_job,binding_job=binding_job,identity=IDENTITY,
        source_emission_jobs={name:jobs[task] for name,jobs in EMISSION_JOBS.items()},
        input_fasta=str(P.E.OUT/seqid/'input.fasta'), direct_paired_DNA_equality_before_reference_access=True,
        primary_CDS_chain_predictions=len(transcripts),scopes=scopes,new_GFF_written=True,
        score_kind='fixed_minimum_centered_fifteen_role_whole_contig_gain',
        pooled_role_arrays_saved=False, original_predictions_modified=False, weights_modified=False,
        reference_used=False, prediction_GFF_content_read=False, labels_or_feature_cache_used=False,
        observation_model_used_at_inference=False,native_forward=False,checkpoint_read=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,threshold_selected=False,
        wall_loop_seconds=time.perf_counter()-began,artifact_bytes_before_summary=artifact_bytes))


def evaluate(checks_job, raw_jobs, binding_job):
    qualification(checks_job); bound(binding_job); began = time.perf_counter(); runs = {}
    if len(raw_jobs) != 6 or len(set(raw_jobs)) != 6:
        raise ValueError('Six distinct actual new raw job IDs required')
    for task, (species,seqid,length) in enumerate(P.E.I.SCOPE):
        folder = OUT/seqid; report = read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED'
                or (folder/'JOBID').read_text().strip() != raw_jobs[task]
                or not report['completed'] or report['job_id'] != raw_jobs[task]
                or report['protocol'] != PROTOCOL or report['checks_job'] != checks_job or report['binding_job'] != binding_job
                or report['identity'] != IDENTITY
                or (report['task'],report['species'],report['seqid'],report['length']) != (task,species,seqid,length)
                or not report['direct_paired_DNA_equality_before_reference_access']
                or any(report[k] for k in ('reference_used','prediction_GFF_content_read',
                    'labels_or_feature_cache_used','observation_model_used_at_inference','native_forward',
                    'checkpoint_read','optimizer_steps_this_run','DEV_used','sealed_data_used',
                    'threshold_selected','weights_modified','original_predictions_modified','pooled_role_arrays_saved'))
                or len(report['scopes']) != 2 or {s['strand'] for s in report['scopes']} != {'+','-'}
                or any(not s['whole_oriented_contig'] or s['reference_used']
                       or s['score_accumulation'] != IDENTITY['score_accumulation'] for s in report['scopes'])):
            raise ValueError('All six fixed raw outputs must complete before any reference access')
        runs[seqid] = report
    targets, reference_policy = comparator_sources()
    write(OUT/'evaluation/started.json', dict(job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        checks_job=checks_job,binding_job=binding_job,raw_jobs=raw_jobs,identity=IDENTITY,source_comparator_jobs=COMPARATOR_JOBS,all_raw_complete_before_reference_access=True))
    chromosomes, novelty = {}, {}
    for species,seqid,length in P.E.I.SCOPE:
        sequence = C.M.read_allowed(P.E.OUT/seqid/'input.fasta',{seqid:length})[seqid]
        other = C.M.read_allowed(P.M.OUT/seqid/'input.fasta',{seqid:length})[seqid]
        if sequence != other: raise ValueError('Exact paired input DNA required before reference access')
        del other
        txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(OUT/seqid/'full_grammar.gff3',
            {seqid:length}),complete_only=False)
        pk = {C.S.key(t) for t in txs}
        if len(txs) != len(pk) or len(txs) != runs[seqid]['primary_CDS_chain_predictions']:
            raise ValueError('New GFF3 roundtrip must preserve unique output chains')
        for tx in txs:
            failures,_ = C.M.chain_checks(tx,sequence)
            if failures or len(tx['CDS']) > 128: raise ValueError('Exported pooled DNA/phase/cap failure')
        for strand in ('+','-'):
            spans = sorted((t['CDS'][0][0],t['CDS'][-1][1]) for t in txs if t['strand']==strand)
            if any(b>c for (_,b),(c,_) in zip(spans,spans[1:])):
                raise ValueError('New same-strand selected gene spans must remain nonoverlapping')
        keys = {name:{C.S.key(t) for t in C.S.E.primary_transcripts(C.S.E.parse_annotation(
            module.OUT/seqid/'prediction/full_grammar.gff3',{seqid:length}),complete_only=False)}
            for name,(module,_) in P.MODELS.items()}
        groups = {'matches_both_parents':pk & keys['E_RANK4'] & keys['M_DENSE'],
            'matches_E_only':(pk & keys['E_RANK4'])-keys['M_DENSE'],
            'matches_M_only':(pk & keys['M_DENSE'])-keys['E_RANK4'],
            'new_chain_neither_parent_emitted':pk-keys['E_RANK4']-keys['M_DENSE']}
        data = ROOT/'data/m1_screen'/species
        coding = C.S.E.parse_annotation(data/'reference.gff3',{seqid:length},protein_coding_only=True)
        annotation = C.S.E.parse_annotation(data/'reference.gff3',{seqid:length})
        with (data/'reference.gff3').open() as handle: _,_,raw_cds,_ = C.P.raw_features(handle,seqid)
        refs,isoforms = C.select_reference(coding,raw_cds,{seqid:sequence})
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('All assessable primary references must retain consistent phase')
        overlap = C.M.C.overlapping_ids(refs); unsupported = C.Counter()
        for tx in refs:
            failures,_ = C.M.chain_checks(tx,sequence)
            if tx['id'] in overlap: failures.append('same_strand_primary_overlap')
            if len(tx['CDS']) > 128: failures.append('over_128_exons')
            if tx['CDS'][-1][1]-tx['CDS'][0][0] > C.M.WIDTH: failures.append('longer_than_native_window')
            unsupported.update(set(failures))
        row = C.score(txs,refs,isoforms,length,annotation)
        row['prediction_nonexact_primary'] = row.pop('reference_nonexact_primary')
        row.update(species=species,seqid=seqid,length=length,reference_primary_records=len(refs),
            reference_assessable_isoform_records=len(isoforms),retained_unsupported_reference_strata=dict(unsupported),
            phase_checked_prediction_records=len(txs),raw_inference_loop_seconds=runs[seqid]['wall_loop_seconds'])
        if any(N.reference_metadata(row) != N.reference_metadata(old['chromosomes'][seqid])
               for old in targets.values()):
            raise ValueError('Full reference population, strata, base denominator and background must match')
        chromosomes[seqid] = row
        rk,ik = ({C.S.key(t) for t in block} for block in (refs,isoforms))
        novelty[seqid] = {g:dict(unique_predictions=len(k),exact_primary=len(k & rk),
            exact_other_assessable_isoform=len((k & ik)-rk),
            neither_exact_primary_nor_any_assessable_isoform=len(k-ik)) for g,k in groups.items()}
        print(json.dumps(dict(completed_pooled_TRAIN_evaluation=seqid)),flush=True)
    per_species = {s:F.aggregate([r for r in chromosomes.values() if r['species']==s])
                   for s in ('arabidopsis_thaliana','oryza_sativa')}
    pooled = F.aggregate(list(chromosomes.values())); f1s = [r['exact_chain']['f1'] for r in per_species.values()]
    if pooled['exact_chain']['reference'] != 22688: raise ValueError('All22688 primary references retained')
    method = dict(chromosomes=chromosomes,per_species=per_species,pooled=pooled,
        macro_species_exact_chain_F1=sum(f1s)/2,worst_species_exact_chain_F1=min(f1s))
    comparisons = {name:dict(pooled_F1_pp=100*(pooled['exact_chain']['f1']-t['pooled']['exact_chain']['f1']),
        per_species_F1_pp={s:100*(per_species[s]['exact_chain']['f1']-t['per_species'][s]['exact_chain']['f1'])
                          for s in per_species},
        macro_species_F1_pp=100*(method['macro_species_exact_chain_F1']-sum(t['per_species'][s]['exact_chain']['f1'] for s in per_species)/2),
        worst_species_F1_pp=100*(method['worst_species_exact_chain_F1']-min(t['per_species'][s]['exact_chain']['f1'] for s in per_species))) for name,t in targets.items()}
    artifact_bytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
    if artifact_bytes>256*1024**2: raise ValueError('Declared total new output budget exceeded')
    write(OUT/'evaluation/summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        checks_job=checks_job,binding_job=binding_job,raw_jobs=raw_jobs,identity=IDENTITY,
        source_comparator_jobs=COMPARATOR_JOBS,energy_pool=method,novelty=novelty,
        energy_pool_minus_comparator=comparisons,reference_policy=reference_policy,
        all_raw_complete_before_reference_access=True,all_original_primary_references_retained=True,
        full_reference_metadata_matches_all_comparators=True,unsupported_reference_strata_retained=True,
        reference_used_in_inference=False,original_predictions_modified=False,weights_modified=False,
        native_forward=False,checkpoint_read=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        threshold_selected=False,frozen_DEV_family_changed=False,
        artifact_bytes_before_summary=artifact_bytes,
        interpretation='post_outcome_TRAIN_fixed_cached_role_energy_pool_not_independent_validation_or_standalone_raw_caller',
        seconds=time.perf_counter()-began))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('bindings','raw','evaluate'))
    parser.add_argument('--task',type=int,choices=range(6))
    parser.add_argument('--checks-job')
    parser.add_argument('--raw-jobs',nargs=6)
    parser.add_argument('--binding-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'bindings': bindings(args.checks_job)
    elif args.mode == 'raw': raw(args.task,args.checks_job,args.binding_job)
    else: evaluate(args.checks_job,args.raw_jobs,args.binding_job)
