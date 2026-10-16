#!/usr/bin/env python3
"""Whole original-six TRAIN structural posteriors; immutable E/M predictions."""
import argparse
import json
import math
import os
import time

import numpy as np
import torch

import common_ruler as C
import infer_full_grammar_joint_lora_train as E
import infer_full_grammar_dense_qv_train as M
from src.m28.full_grammar_gene_posterior import complete_gene_set_posteriors
from src.m28.gene_posterior_ranking import oriented_chain, ranking_summary, probability_diagnostics

ROOT = E.ROOT
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_gene_posterior_train_plan.md'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-CHECKS'
MODELS = {'E_RANK4': (E, '13499817'), 'M_DENSE': (M, '13510324')}


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def source(name):
    module, fit_job = MODELS[name]
    module.source(fit_job)
    reports = module.reconcile_inputs()
    receipt = read(module.EVAL/'paired_input_match.json')
    if (receipt['protocol'] != module.A.PROTOCOL or receipt['adapted_fit_job'] != fit_job
            or not receipt['original_six_TRAIN_inputs_directly_equal']
            or not receipt['comparison_before_reference_access'] or receipt['DEV_used'] or receipt['sealed_data_used']):
        raise ValueError('Existing exact original-six normalized DNA receipt required')
    return module, fit_job, reports


def qualification(job):
    folder = CHECKS/job
    result = read(folder/'summary.json')
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or result['job_id'] != job
            or result['protocol'] != PROTOCOL or not result['completed']
            or result['tests_run'] != 4 or result['failures'] or result['errors']
            or not result['old_E_M_raw_source_receipts_verified']):
        raise ValueError('New coordinate/ranking/source qualification required')


def checks():
    import unittest
    import test_m28_gene_posterior_ranking as T
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(T))
    if not result.wasSuccessful():
        raise RuntimeError('New coordinate/ranking arithmetic checks failed')
    for name in MODELS:
        module, _, reports = source(name)
        for _, seqid, length in E.I.SCOPE:
            predictions(module, seqid, length, reports[seqid])
    old = read(ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-PREFLIGHT/summary.json')
    if (not old['completed'] or old['job_id'] != '13547015' or old['fixed_native_cases'] != 8
            or not old['denominator_all_legal_gene_sets_not_selected_pool']
            or any(not r['original_partition_and_direct_event_match'] or not r['original_MAP_scores_preserved']
                   for r in old['records'])):
        raise ValueError('Completed previous eight native event cases required')
    write(CHECKS/os.environ['SLURM_JOB_ID']/'summary.json', dict(completed=True,
        job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors), old_E_M_raw_source_receipts_verified=True,
        original_native_preflight_job='13547015', research_DNA_read=False, research_reference_read=False,
        prediction_GFF_content_read=True, checkpoint_read=False, native_model_read=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began))


def encode(name, task, checks_job):
    qualification(checks_job)
    module, fit_job, reports = source(name)
    species, seqid, length = E.I.SCOPE[task]
    out = OUT/name/seqid/'emissions'
    began = time.perf_counter()
    write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        model=name, task=task, species=species, seqid=seqid, length=length, checks_job=checks_job,
        source_fit_job=fit_job, checkpoint=str(module.CHECKPOINT), reference_used=False, DEV_used=False))
    sequence = module.I.R.M.read_allowed(module.OUT/seqid/'input.fasta', {seqid:length})[seqid]
    if set(sequence)-set('ACGTN'):
        raise ValueError('Previously normalized strict raw TRAIN DNA required')
    locations = [(a, min(a+module.A.R.M.WIDTH, length)) for a in module.A.R.M.starts(length)]
    if any(s['windows'] != [list(w) for w in locations] for s in reports[seqid]['scopes']):
        raise ValueError('Exact original whole-contig window geometry required')
    state = torch.load(module.CHECKPOINT, map_location='cpu', weights_only=True)
    model = module.A.model_from_state(state, 'cuda')
    backbone, tokenizer = module.A.load_native_adapted(state)
    del state
    torch.cuda.reset_peak_memory_stats()
    records = []
    for strand, filename in (('+', 'plus.npy'), ('-', 'minus.npy')):
        started = time.perf_counter()
        with torch.inference_mode():
            joined = module.A.J.contig_emissions(model, backbone, tokenizer, sequence, locations, strand, 'cuda')
        values = joined['local_potentials']
        with (out/filename).open('xb') as handle:
            np.save(handle, values.numpy(), allow_pickle=False)
        records.append(dict(strand=strand, filename=filename, shape=list(values.shape), dtype=str(values.dtype),
            genomic_context_segments=joined['genomic_context_segments'], encoding_and_write_seconds=time.perf_counter()-started))
        del joined, values
    artifact_bytes = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    if artifact_bytes > length*15*4*2+2*1024*1024:
        raise ValueError('Declared two FP32 role arrays plus small metadata budget exceeded')
    write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        model=name, task=task, species=species, seqid=seqid, length=length, checks_job=checks_job,
        source_fit_job=fit_job, checkpoint=str(module.CHECKPOINT), windows=locations, records=records,
        source_DNA_equality_receipt_verified=True, whole_original_TRAIN_contig_both_strands=True,
        reference_used=False, prediction_GFF_content_read=False, labels_or_feature_cache_used=False,
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0, weights_modified=False,
        original_predictions_modified=False, artifact_bytes_before_summary=artifact_bytes,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(), seconds=time.perf_counter()-began))


def predictions(module, seqid, length, report):
    path = module.OUT/seqid/'prediction/full_grammar.gff3'
    txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(path,
        {seqid:length}), complete_only=False)
    if len(txs) != report['primary_CDS_chain_predictions'] or len({C.S.key(t) for t in txs}) != len(txs):
        raise ValueError('Original exported prediction count/unique chains must be preserved')
    # The common evaluation parser intentionally discards GFF's score column.
    # Recover only the original export's mRNA scores; do not alter that parser.
    scores = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] == 'mRNA':
                scores[C.S.E.parse_attrs(fields[8])['ID']] = float(fields[5])
    if set(scores) != {t['id'] for t in txs} or any(not math.isfinite(v) for v in scores.values()):
        raise ValueError('Original finite exported mRNA gain required for every unique query')
    for tx in txs:
        tx['score'] = scores[tx['id']]
    return txs


def posterior(name, task, checks_job, emissions_job):
    qualification(checks_job)
    module, fit_job, reports = source(name)
    species, seqid, length = E.I.SCOPE[task]
    emitted = OUT/name/seqid/'emissions'; out = OUT/name/seqid/'posterior'
    report = read(emitted/'summary.json')
    if ((emitted/'STATUS').read_text().strip() != 'COMPLETED' or report['job_id'] != emissions_job
            or (report['model'], report['task'], report['seqid'], report['length']) != (name, task, seqid, length)
            or report['checkpoint'] != str(module.CHECKPOINT) or report['checks_job'] != checks_job
            or report['protocol'] != PROTOCOL or not report['completed'] or report['reference_used']
            or report['DEV_used'] or report['sealed_data_used'] or report['optimizer_steps_this_run']):
        raise ValueError('This exact completed model/chromosome emission task required')
    began = time.perf_counter()
    write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL, model=name,
        task=task, species=species, seqid=seqid, length=length, source_emissions_job=emissions_job,
        checks_job=checks_job, reference_used=False, native_forward=False, DEV_used=False))
    sequence = module.I.R.M.read_allowed(module.OUT/seqid/'input.fasta', {seqid:length})[seqid]
    txs = predictions(module, seqid, length, reports[seqid])
    rows, records = [], []
    for strand, filename in (('+', 'plus.npy'), ('-', 'minus.npy')):
        ordered = sorted((t for t in txs if t['strand'] == strand),
            key=lambda t: t['CDS'][0][0] if strand == '+' else -t['CDS'][-1][1])
        chains = [oriented_chain([(a,b) for a,b,_ in t['CDS']], length, strand) for t in ordered]
        dna = sequence if strand == '+' else module.I.R.M.rc(sequence)
        # Read-only mmap; posterior never modifies inputs and scans one orientation at a time.
        array = np.load(emitted/filename, mmap_mode='c', allow_pickle=False)
        if array.shape != (length, 15) or array.dtype != np.float32:
            raise ValueError('Full chromosome original fifteen-role FP32 array required')
        values = torch.from_numpy(array)
        started = time.perf_counter()
        result = complete_gene_set_posteriors(dna, values, chains, 128)
        elapsed = time.perf_counter()-started
        # GFF scores have eight significant digits, unlike the retained FP64 energies.
        for tx, chain, gain, logp in zip(ordered, chains, result['chain_gains'], result['gene_log_probabilities']):
            if not math.isclose(gain, float(tx['score']), rel_tol=6e-8, abs_tol=1e-5):
                raise ValueError('Fresh whole-chromosome gene gain differs from immutable original MAP GFF score')
            rows.append(dict(id=tx['id'], strand=strand, genomic_CDS=[(a,b) for a,b,_ in tx['CDS']],
                oriented_CDS=chain, gain=gain, original_rounded_GFF_gain=float(tx['score']), log_probability=logp))
        records.append(dict(strand=strand, supplied_bases=len(dna), queried_genes=len(chains),
            log_partition=result['log_partition'], reverse_log_partition=result['reverse_log_partition'],
            queried_set_log_probability=result['queried_set_log_probability'],
            empty_set_log_probability=result['empty_set_log_probability'],
            state_count=result['state_count'], stored_boundary_scalars=result['stored_boundary_scalars'], posterior_seconds=elapsed))
        print(json.dumps(dict(stage='whole_chromosome_posterior', model=name, seqid=seqid, strand=strand, seconds=elapsed)), flush=True)
        del result, values, array, dna
    write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        model=name, task=task, species=species, seqid=seqid, length=length, checks_job=checks_job,
        source_emissions_job=emissions_job, source_fit_job=fit_job, checkpoint=str(module.CHECKPOINT),
        genes=rows, records=records, source_DNA_equality_receipt_verified=True,
        full_chromosome_conditional_posterior_measured=True, denominator_all_legal_gene_sets_not_selected_pool=True,
        original_MAP_scores_and_exported_chains_preserved=True, predictions_modified=False,
        reference_used=False, labels_or_feature_cache_used=False, native_forward=False, checkpoint_read=False,
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0,
        posterior_threshold_selected=False, fitted_temperature=False, seconds=time.perf_counter()-began,
        interpretation='model_structural_support_not_empirically_calibrated_biological_correctness'))


def diagnostics(rows, reference_count):
    return {label:dict(posterior=ranking_summary(rows, 'log_probability', label,
                    reference_count if label == 'exact_primary' else None),
                internal_gain=ranking_summary(rows, 'gain', label,
                    reference_count if label == 'exact_primary' else None),
                probability=probability_diagnostics(rows, label))
            for label in ('exact_primary', 'exact_any_assessable_isoform')}


def evaluate(checks_job):
    qualification(checks_job)
    began = time.perf_counter()
    raw, source_modules = {}, {}
    # All twelve posterior tasks qualify before any real reference is opened.
    for name in MODELS:
        module, _, reports = source(name); source_modules[name] = module
        for task, (species, seqid, length) in enumerate(E.I.SCOPE):
            folder = OUT/name/seqid/'posterior'; report = read(folder/'summary.json')
            if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                    or report['protocol'] != PROTOCOL or report['checks_job'] != checks_job
                    or (report['model'], report['task'], report['species'], report['seqid'], report['length']) != (name, task, species, seqid, length)
                    or report['checkpoint'] != str(module.CHECKPOINT) or report['reference_used']
                    or not report['full_chromosome_conditional_posterior_measured']
                    or not report['original_MAP_scores_and_exported_chains_preserved']
                    or report['DEV_used'] or report['sealed_data_used'] or report['predictions_modified']
                    or len(report['genes']) != reports[seqid]['primary_CDS_chain_predictions']):
                raise ValueError('All twelve fixed whole-chromosome posterior receipts required before truth')
            raw[name, seqid] = report
    models = {name:dict(chromosomes={}, per_species={}) for name in MODELS}
    by_species = {name:{} for name in MODELS}; reference_counts = {}
    for species, seqid, length in E.I.SCOPE:
        sequences = C.M.read_allowed(E.OUT/seqid/'input.fasta', {seqid:length})
        data = ROOT/'data/m1_screen'/species
        coding = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length}, protein_coding_only=True)
        with (data/'reference.gff3').open() as handle:
            _, _, raw_cds, _ = C.P.raw_features(handle, seqid)
        refs, isoforms = C.select_reference(coding, raw_cds, sequences)
        primary, all_iso = ({C.S.key(t) for t in ts} for ts in (refs, isoforms))
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Retained primary reference phase must remain unchanged')
        reference_counts[seqid] = len(primary)
        for name, module in source_modules.items():
            rows = raw[name, seqid]['genes']
            for row in rows:
                key = C.S.key(dict(seqid=seqid, strand=row['strand'], CDS=[(a,b,'0') for a,b in row['genomic_CDS']]))
                row['exact_primary'] = key in primary
                row['exact_any_assessable_isoform'] = key in all_iso
            old = read(module.EVAL/'summary.json')['chromosomes'][seqid]['exact_chain']
            if (len(rows), sum(r['exact_primary'] for r in rows), len(primary)) != (old['predicted'], old['tp'], old['reference']):
                raise ValueError('Immutable predictions and full retained primary endpoints must replay')
            models[name]['chromosomes'][seqid] = diagnostics(rows, len(primary))
            by_species[name].setdefault(species, []).extend(rows)
    for name in MODELS:
        pooled = []
        for species, rows in by_species[name].items():
            count = sum(reference_counts[q] for s,q,_ in E.I.SCOPE if s == species)
            models[name]['per_species'][species] = diagnostics(rows, count); pooled.extend(rows)
        models[name]['pooled'] = diagnostics(pooled, sum(reference_counts.values()))
    write(OUT/'evaluation/summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        protocol=PROTOCOL, checks_job=checks_job, models=models, retained_primary_reference_counts=reference_counts,
        all_twelve_posteriors_verified_before_reference_access=True, original_exact_chain_counts_replayed=True,
        reference_used_only_after_readout=True, native_forward=False, checkpoint_read=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, predictions_modified=False,
        threshold_selected=False, calibration_fitted=False, seconds=time.perf_counter()-began,
        interpretation='TRAIN_selected_prediction_reference_relative_ranking_not_generalization_or_biological_calibration'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks', 'encode', 'posterior', 'evaluate'))
    parser.add_argument('--model', choices=tuple(MODELS))
    parser.add_argument('--task', type=int, choices=range(6))
    parser.add_argument('--checks-job')
    parser.add_argument('--emissions-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'checks':
        checks()
    elif args.checks_job is None:
        parser.error('--checks-job required')
    elif args.mode == 'evaluate':
        evaluate(args.checks_job)
    elif args.model is None or args.task is None:
        parser.error('--model and --task required')
    elif args.mode == 'encode':
        encode(args.model, args.task, args.checks_job)
    elif args.emissions_job is None:
        parser.error('--emissions-job required')
    else:
        posterior(args.model, args.task, args.checks_job, args.emissions_job)
