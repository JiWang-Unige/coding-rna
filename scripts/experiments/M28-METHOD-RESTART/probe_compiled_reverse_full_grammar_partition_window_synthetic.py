#!/usr/bin/env python3
"""Cold compiled-reverse full artificial window; unchanged original forward."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import torch

import profile_emission_lattice_synthetic as S
from src.m28.full_grammar_partition_compiled_reverse import (compiled_reverse_full_grammar_log_partition,
    _accumulate_reverse_nodes)


LENGTH = 24576
CHECKS = Path('outputs/M28-FULL-GRAMMAR-COMPILED-REVERSE-CHECKS')
CHECKS_JOB = '13427223'


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)


def main(out):
    torch.set_num_threads(2)
    checked = json.loads((CHECKS/'summary.json').read_text())
    if (CHECKS.joinpath('STATUS').read_text().strip() != 'COMPLETED'
            or CHECKS.joinpath('JOBID').read_text().strip() != CHECKS_JOB
            or checked['job_id'] != CHECKS_JOB or checked['completed'] is not True
            or checked['tests_run'] != 10 or checked['failures'] != 0 or checked['errors'] != 0
            or checked['mathematical_partition_and_gradient_checks_passed'] is not True
            or checked['synthetic_only'] is not True or checked['optimizer_steps_this_run'] != 0
            or checked['implementation'] != 'sequential_nopython_reverse_only_original_forward'
            or checked['original_mathematical_checks_job'] != '13426777'
            or checked['indexed_mathematical_checks_job'] != '13427003'
            or checked['higher_order_gradients_supported'] is not False
            or checked['GPU_supported'] is not False
            or checked['explicit_reverse_checks_job'] != '13427085'
            or checked['actual_nopython_signature_count'] < 1):
        raise ValueError('The assigned actual mathematical checks must be closed')
    if len(_accumulate_reverse_nodes.nopython_signatures) != 0:
        raise ValueError('This new process must not prewarm the compiled reverse')
    previous_root = Path('outputs/M28-FULL-GRAMMAR-REVERSE-WINDOW-PROBE')
    previous_started = json.loads((previous_root/'started.json').read_text())
    if (previous_root.joinpath('JOBID').read_text().strip() != '13427113'
            or previous_started['job_id'] != '13427113'
            or (previous_started['length'], previous_started['max_exons'], previous_started['DNA_seed'],
                previous_started['potential_seed']) != (LENGTH, 128, 0, 41)):
        raise ValueError('Original complete artificial forward identity differs')
    previous_forward = json.loads((previous_root/'forward.json').read_text())
    dna_generator = torch.Generator().manual_seed(0)
    # Match the original artificial DNA generation, then take its first window.
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (49152,),
        generator=dna_generator).tolist())[:LENGTH]
    generator = torch.Generator().manual_seed(41)
    logits = torch.randint(-3, 4, (LENGTH, 15), generator=generator).float().requires_grad_(True)
    started = dict(job_id=os.environ['SLURM_JOB_ID'], source_mathematical_checks_job=CHECKS_JOB,
        implementation='sequential_nopython_reverse_only_original_forward', original_naive_probe_job='13426855', original_indexed_probe_job='13427015', original_explicit_reverse_probe_job='13427113',
        nopython_signatures_before_backward=len(_accumulate_reverse_nodes.nopython_signatures),
        length=LENGTH, strand='+', DNA_seed=0, DNA_generation_length=49152,
        potential_seed=41, potential_definition='first_window_CPU_FP32_integers_minus3_to3_all15_roles',
        max_exons=128, null=0, role_weights=1, score_accumulation='FP64',
        torch_version=str(torch.__version__), threads=2, synthetic_only=True,
        motif_counts=S.motif_counts(sequence), started_utc=S.utc_now())
    save(out/'started.json', started)
    print(json.dumps(dict(stage='full_window_forward_started', **started)), flush=True)
    began = time.perf_counter()
    result = compiled_reverse_full_grammar_log_partition(sequence, logits, 128)
    forward_seconds = time.perf_counter()-began
    value = result['log_partition']
    if not bool(torch.isfinite(value)) or float(value.detach()) < 0:
        raise ValueError('The complete null-inclusive partition must be finite and nonnegative')
    forward = dict(forward_seconds=forward_seconds, log_partition=float(value.detach()),
        counts=result['counts'], process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        forward_completed_utc=S.utc_now())
    if (forward['log_partition'] != previous_forward['log_partition']
            or forward['counts'] != previous_forward['counts']):
        raise ValueError('Full original logZ and all recurrence/tape counts must be identical')
    save(out/'forward.json', forward)
    print(json.dumps(dict(stage='full_window_backward_started', **forward)), flush=True)
    began = time.perf_counter()
    value.backward()
    backward_seconds = time.perf_counter()-began
    gradient = logits.grad
    if (gradient is None or tuple(gradient.shape) != (LENGTH, 15)
            or gradient.dtype != torch.float32 or not bool(torch.isfinite(gradient).all())):
        raise ValueError('All original FP32 input positions must receive finite gradients')
    row_shift_residual = float(gradient.sum(1).abs().max())
    if row_shift_residual > 1e-5:
        raise ValueError('Full-window common-row logit-shift gradient exceeds FP32 tolerance')
    minimum_role_gradient = float(gradient[:, 1:].min())
    maximum_covered_probability = float(gradient[:, 1:].sum(1).max())
    if minimum_role_gradient < -1e-6 or maximum_covered_probability > 1 + 1e-6:
        raise ValueError('All-path unit-role derivatives must be nonnegative occupancy probabilities at most one')
    summary = dict(**started, completed=True, completed_utc=S.utc_now(),
        protocol='full_grammar_partition_preparation.md', forward=forward,
        backward_seconds=backward_seconds, common_row_shift_gradient_max_abs=row_shift_residual,
        process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        actual_complete_forward_logZ_and_all_counts_equal_original=True,
        nopython_signatures_after_backward=len(_accumulate_reverse_nodes.nopython_signatures),
        cold_nopython_compilation_in_backward=True,
        minimum_nonintergenic_role_gradient=minimum_role_gradient,
        maximum_nonintergenic_occupancy_probability=maximum_covered_probability,
        all_input_role_posterior_bounds_verified=True,
        full_window_forward_and_backward_completed=True, all_input_gradients_finite=True,
        optimizer_steps_this_run=0, model_or_backbone_forward=False,
        research_inputs_opened=False, labels_or_catalog_opened=False, DEV_used=False,
        independent_species_used=False, observed_mask_or_gold_likelihood_implemented=False,
        training_or_inference_caller_changed=False,
        interpretation='single_artificial_full_window_CPU_resource_probe_not_fit_GPU_feasibility_or_performance')
    save(out/'summary.json', summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    main(parser.parse_args().out)
