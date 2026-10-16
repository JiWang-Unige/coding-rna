#!/usr/bin/env python3
"""Cold full-length synthetic boundary-message resource gate; no research input."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import time

import torch

import profile_emission_lattice_synthetic as S
from src.m28.full_grammar_boundary_compiled import (
    _forward_messages, _reverse_messages, compiled_boundary_full_grammar_log_partition)

LENGTH = 24576
CHECKS = Path('outputs/M28-FULL-GRAMMAR-BOUNDARY-COMPILED-CHECKS')
CHECKS_JOB = '13428495'
BASELINE = Path('outputs/M28-FULL-GRAMMAR-COMPILED-REVERSE-WINDOW-PROBE')
BASELINE_JOB = '13427267'


def save(path, value):
    with path.open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False)


def main(out):
    torch.set_num_threads(2)
    checked = json.loads((CHECKS/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or (CHECKS/'JOBID').read_text().strip() != CHECKS_JOB
            or checked['job_id'] != CHECKS_JOB or not checked['completed'] or checked['tests_run'] != 10
            or checked['failures'] or checked['errors'] or not checked['synthetic_only']
            or not checked['complete_gene_set_logZ_and_all_gradients_match_three_oracles']
            or not checked['every_supported_output_and_incoming_message_gradient_verified']
            or not checked['checkpoint_recomputation_invariance_verified'] or checked['optimizer_steps_this_run'] != 0):
        raise ValueError('Matching completed compiled-message mathematics required')
    previous = json.loads((BASELINE/'summary.json').read_text())
    if ((BASELINE/'STATUS').read_text().strip() != 'COMPLETED' or (BASELINE/'JOBID').read_text().strip() != BASELINE_JOB
            or previous['job_id'] != BASELINE_JOB or not previous['completed'] or not previous['synthetic_only']
            or not previous['full_window_forward_and_backward_completed'] or previous['optimizer_steps_this_run'] != 0
            or (previous['length'],previous['DNA_generation_length'],previous['DNA_seed'],previous['potential_seed'],
                previous['max_exons'],previous['role_weights'],previous['null'],previous['score_accumulation'])
                != (LENGTH,49152,0,41,128,1,0,'FP64')
            or previous['forward']['log_partition'] != 10992.597775197577):
        raise ValueError('Original completed full synthetic input identity differs')
    signatures_before = (len(_forward_messages.nopython_signatures),len(_reverse_messages.nopython_signatures))
    if signatures_before != (0,0): raise ValueError('Fresh process must not prewarm either kernel')
    dna_generator = torch.Generator().manual_seed(0)
    sequence = ''.join('ACGT'[i] for i in torch.randint(0,4,(49152,),generator=dna_generator).tolist())[:LENGTH]
    potential_generator = torch.Generator().manual_seed(41)
    logits = torch.randint(-3,4,(LENGTH,15),generator=potential_generator).float().requires_grad_(True)
    motif_counts = S.motif_counts(sequence)
    if motif_counts != previous['motif_counts']: raise ValueError('Full synthetic DNA motif counts differ')
    started = dict(job_id=os.environ['SLURM_JOB_ID'],source_mathematical_checks_job=CHECKS_JOB,
        full_synthetic_baseline_job=BASELINE_JOB,implementation='checkpointed_nopython_complete_gene_boundary_messages',
        length=LENGTH,strand='+',DNA_seed=0,DNA_generation_length=49152,potential_seed=41,
        potential_definition='first_window_CPU_FP32_integers_minus3_to3_all15_roles',max_exons=128,
        role_weights=1,null=0,score_accumulation='FP64',checkpoint_bp=256,threads=2,synthetic_only=True,
        motif_counts=motif_counts,torch_version=str(torch.__version__),started_utc=S.utc_now(),
        forward_signatures_before=signatures_before[0],reverse_signatures_before=signatures_before[1],
        full_logZ_comparison_rtol=0,full_logZ_comparison_atol=1e-8)
    save(out/'started.json',started)
    print(json.dumps(dict(stage='full_window_forward_started',**started)),flush=True)
    began = time.perf_counter()
    result = compiled_boundary_full_grammar_log_partition(sequence,logits,128,checkpoint_bp=256)
    forward_seconds = time.perf_counter()-began
    value = result['log_partition']; actual = float(value.detach())
    difference = actual-previous['forward']['log_partition']
    forward = dict(forward_seconds=forward_seconds,log_partition=actual,baseline_logZ_difference=difference,
        statistics=[dict(item) for item in result['statistics']],
        process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        forward_signatures_after=len(_forward_messages.nopython_signatures),forward_completed_utc=S.utc_now())
    save(out/'forward.json',forward)
    if not math.isfinite(actual) or actual < 0 or not math.isclose(actual,previous['forward']['log_partition'],rel_tol=0,abs_tol=1e-8):
        raise ValueError('Complete fixed-input null-inclusive logZ differs beyond the declared FP64 tolerance')
    print(json.dumps(dict(stage='full_window_backward_started',**forward)),flush=True)
    began = time.perf_counter(); value.backward(); backward_seconds = time.perf_counter()-began
    gradient = logits.grad
    if (gradient is None or tuple(gradient.shape) != (LENGTH,15) or gradient.dtype != torch.float32
            or not bool(torch.isfinite(gradient).all())):
        raise ValueError('Every original FP32 potential row must receive finite gradients')
    row_shift = float(gradient.sum(1).abs().max())
    minimum = float(gradient[:,1:].min()); maximum = float(gradient[:,1:].sum(1).max())
    if row_shift > 1e-5 or minimum < -1e-6 or maximum > 1+1e-6:
        raise ValueError('All-input role occupancy or common-row-shift derivative exceeds the original tolerance')
    if len(_forward_messages.nopython_signatures) < 1 or len(_reverse_messages.nopython_signatures) < 1:
        raise ValueError('Both full-length kernels must actually compile in nopython mode')
    summary = dict(**started,completed=True,completed_utc=S.utc_now(),
        protocol='reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md',forward=forward,
        backward_seconds=backward_seconds,statistics_after_backward=[dict(item) for item in result['statistics']],
        process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        forward_signatures_after=len(_forward_messages.nopython_signatures),reverse_signatures_after=len(_reverse_messages.nopython_signatures),
        cold_forward_includes_state_table_and_JIT=True,cold_backward_includes_JIT=True,
        complete_fixed_input_logZ_matches_original_with_declared_tolerance=True,
        common_row_shift_gradient_max_abs=row_shift,minimum_nonintergenic_role_gradient=minimum,
        maximum_nonintergenic_occupancy_probability=maximum,all_input_role_posterior_bounds_verified=True,
        full_window_forward_and_backward_completed=True,all_input_gradients_finite=True,
        baseline_full_gradient_equality_checked=False,old_recurrence_or_tape_counts_compared=False,
        optimizer_steps_this_run=0,boundary_prior_chosen=False,context_extent_or_padding_rule_selected=False,
        loss_normalizer_or_fit_coefficients_selected=False,model_or_backbone_forward=False,
        research_inputs_opened=False,labels_or_catalog_opened=False,DEV_used=False,independent_species_used=False,
        A33_outputs_opened=False,observed_full_window_NLL_executed=False,original_R4_policy_changed=False,
        observation_rows_deleted=False,training_or_inference_caller_changed=False,
        interpretation='single_full_artificial_CPU_resource_gate_not_real_context_observed_NLL_fit_GPU_or_performance')
    save(out/'summary.json',summary)
    print(json.dumps(summary),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    main(parser.parse_args().out)
