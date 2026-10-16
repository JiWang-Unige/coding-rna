#!/usr/bin/env python3
"""New posterior math checks and eight fixed native TRAIN-view feasibility cases."""
import argparse
import json
import math
import os
from pathlib import Path
import time

import torch

from src.m28.full_grammar_gene_posterior import complete_gene_set_posteriors
from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition
from src.m28.full_grammar_context_view import context_geometry
from src.m28.phase_evidence_chain import candidate_roles
import infer_full_grammar_joint_lora_train as E
import infer_full_grammar_dense_qv_train as M

ROOT = E.ROOT
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_gene_posterior_preflight.md'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-CHECKS'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-PREFLIGHT'
VIEWS = (('arabidopsis_thaliana', 'NC_003070.9', 30427671, 14389248),
         ('oryza_sativa', 'NC_089046.1', 28271460, 14462976))
SOURCES = (('E_RANK4', E, '13499817'), ('M_DENSE', M, '13510324'))


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def checks():
    import unittest
    import test_m28_full_grammar_gene_posterior as T
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(T))
    write(CHECKS/os.environ['SLURM_JOB_ID']/'summary.json', dict(completed=result.wasSuccessful(),
        job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors), research_DNA_read=False,
        research_reference_read=False, checkpoint_read=False, native_model_read=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began))
    if not result.wasSuccessful():
        raise RuntimeError('New complete-gene posterior checks failed')


def native(checks_job):
    began = time.perf_counter()
    qualification = json.loads((CHECKS/checks_job/'summary.json').read_text())
    if ((CHECKS/checks_job/'STATUS').read_text().strip() != 'COMPLETED'
            or qualification['job_id'] != checks_job or qualification['protocol'] != PROTOCOL
            or not qualification['completed'] or qualification['tests_run'] != 4
            or qualification['errors'] or qualification['failures']):
        raise ValueError('Four completed independent posterior-accounting groups required')
    # Existing TRAIN callers already directly matched these normalized input files
    # to original DNA before reference access. Only source metadata is read here.
    input_receipt = json.loads((E.EVAL/'paired_input_match.json').read_text())
    if (input_receipt['adapted_fit_job'] != '13499817' or input_receipt['baseline_fit_job'] != '13494895'
            or not input_receipt['original_six_TRAIN_inputs_directly_equal']
            or not input_receipt['comparison_before_reference_access']
            or input_receipt['DEV_used'] or input_receipt['sealed_data_used']):
        raise ValueError('Completed original-TRAIN caller DNA equality receipt required')
    for _, source, fit_job in SOURCES:
        source.source(fit_job)
    write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL, checks_job=checks_job,
        source_fit_jobs={name:job for name, _, job in SOURCES}, source_head_steps=9216,
        source_native_steps=4608, views=VIEWS, strands=['+', '-'], max_exons=128,
        optimizer_steps_this_run=0, reference_used=False, labels_or_feature_cache_used=False,
        DEV_used=False, sealed_data_used=False))
    sequences = {}
    for species, seqid, length, start in VIEWS:
        folder = E.OUT/seqid
        receipt = json.loads((folder/'summary.json').read_text())
        if (not receipt['inference_complete'] or receipt['source_fit_job'] != '13499817'
                or (receipt['species'], receipt['seqid'], receipt['length']) != (species, seqid, length)
                or receipt['reference_used'] or receipt['labels_or_feature_cache_used']
                or receipt['DEV_used'] or receipt['sealed_data_used']):
            raise ValueError('Two declared completed original-TRAIN FASTA sources required')
        sequences[seqid] = E.I.R.M.read_allowed(folder/'input.fasta', {seqid:length})[seqid]
    records = []
    for name, source, fit_job in SOURCES:
        state = torch.load(source.CHECKPOINT, map_location='cpu', weights_only=True)
        model = source.A.model_from_state(state, 'cuda')
        backbone, tokenizer = source.A.load_native_adapted(state)
        del state
        if any(p.requires_grad for p in [*model.parameters(), *backbone.parameters()]):
            raise ValueError('Native posterior readout must remain gradient-free')
        torch.cuda.reset_peak_memory_stats()
        for species, seqid, length, start in VIEWS:
            chromosome = sequences[seqid]
            geometry = context_geometry(start, start+24576, length)
            a, b = geometry['genomic_context']
            sequence = chromosome[a:b]
            if len(sequence) != 49152 or set(sequence)-set('ACGTN'):
                raise ValueError('Two fixed original normalized49152bp contexts required')
            for strand in ('+', '-'):
                locations = list(geometry['local_windows'])
                encoding_started = time.perf_counter()
                with torch.inference_mode():
                    joined = source.A.J.contig_emissions(model, backbone, tokenizer, sequence, locations, strand, 'cuda')
                    dna, values = joined['sequence'], joined['local_potentials']
                    decoded = source.A.R.decode_full_grammar_skyline(dna, values, 128)
                encoding_and_MAP_seconds = time.perf_counter()-encoding_started
                posterior_started = time.perf_counter()
                result = complete_gene_set_posteriors(dna, values, decoded['chains'], 128)
                posterior_seconds = time.perf_counter()-posterior_started
                if any(not math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-8)
                       for x, y in zip(result['chain_gains'], decoded['chain_gains'])):
                    raise ValueError('Readout changed the original MAP gene energies')
                all_z = float(compiled_boundary_full_grammar_log_partition(dna, values, 128)['log_partition'])
                if not math.isclose(all_z, result['log_partition'], rel_tol=1e-12, abs_tol=1e-8):
                    raise ValueError('New forward/backward denominator differs from original training grammar')
                allowed = torch.ones(len(dna), 15, dtype=torch.bool)
                if decoded['chains']:
                    first = decoded['chains'][0]
                    left, right = first[0][0], first[-1][1]
                    query_roles = candidate_roles(first, len(dna))
                    allowed[left:right] = torch.nn.functional.one_hot(query_roles[left:right], num_classes=15).bool()
                    expected_log_p = result['gene_log_probabilities'][0]
                else:
                    allowed[:] = False; allowed[:, 0] = True
                    expected_log_p = result['empty_set_log_probability']
                event_z = float(compiled_boundary_full_grammar_log_partition(dna, values, 128, allowed)['log_partition'])
                if not math.isclose(event_z-all_z, expected_log_p, rel_tol=1e-10, abs_tol=1e-7):
                    raise ValueError('Event marginal differs from direct forced-prediction-role numerator')
                row = dict(model=name, source_fit_job=fit_job, checkpoint=str(source.CHECKPOINT),
                    species=species, seqid=seqid, genomic_context=[a, b], core_start=start,
                    strand=strand, supplied_view_bases=len(dna), native_windows=locations,
                    posterior=result, encoding_and_MAP_seconds=encoding_and_MAP_seconds,
                    posterior_seconds=posterior_seconds, original_partition_and_direct_event_match=True,
                    original_MAP_scores_preserved=True, complete_gene_event_checked=bool(decoded['chains']),
                    first_event_log_probability=expected_log_p,
                    reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
                    DEV_used=False, sealed_data_used=False,
                    peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
                    peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved())
                records.append(row)
                print(json.dumps(dict(model=name, seqid=seqid, strand=strand,
                    selected_complete_CDS_chains=len(decoded['chains']), posterior_seconds=posterior_seconds)), flush=True)
                del joined, values, decoded, result, allowed
        del model, backbone, tokenizer
        torch.cuda.empty_cache()
    write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        checks_job=checks_job, source_fit_jobs={name:job for name, _, job in SOURCES}, records=records,
        fixed_native_cases=8, source_DNA_equality_receipt_verified=True,
        denominator_all_legal_gene_sets_not_selected_pool=True, original_predictions_modified=False,
        checkpoints_written=False, weights_or_optimizer_updated=False, optimizer_steps_this_run=0,
        reference_used=False, labels_or_feature_cache_used=False, DEV_used=False, sealed_data_used=False,
        posterior_threshold_selected=False, calibration_or_accuracy_measured=False,
        full_chromosome_conditional_posterior_measured=False, seconds=time.perf_counter()-began,
        interpretation='fixed_TRAIN_view_native_math_and_cost_feasibility_not_biological_calibration_or_generalization'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks', 'native'))
    parser.add_argument('--checks-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'checks':
        checks()
    elif args.checks_job is None:
        parser.error('--checks-job required')
    else:
        native(args.checks_job)
