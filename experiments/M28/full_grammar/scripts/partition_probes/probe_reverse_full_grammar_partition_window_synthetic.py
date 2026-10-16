#!/usr/bin/env python3
"""One full-length artificial probe of compact explicit first-order reverse."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import torch

import profile_emission_lattice_synthetic as S
from src.m28.full_grammar_partition_reverse import reverse_full_grammar_log_partition


LENGTH = 24576
CHECKS = Path('outputs/M28-FULL-GRAMMAR-REVERSE-CHECKS')
CHECKS_JOB = '13427085'


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)


def main(out):
    torch.set_num_threads(2)
    checked = json.loads((CHECKS/'summary.json').read_text())
    if (CHECKS.joinpath('STATUS').read_text().strip() != 'COMPLETED'
            or CHECKS.joinpath('JOBID').read_text().strip() != CHECKS_JOB
            or checked['job_id'] != CHECKS_JOB or checked['completed'] is not True
            or checked['tests_run'] != 9 or checked['failures'] != 0 or checked['errors'] != 0
            or checked['mathematical_partition_and_gradient_checks_passed'] is not True
            or checked['synthetic_only'] is not True or checked['optimizer_steps_this_run'] != 0
            or checked['implementation'] != 'compact_grammar_specific_explicit_first_order_CPU_reverse'
            or checked['original_mathematical_checks_job'] != '13426777'
            or checked['indexed_mathematical_checks_job'] != '13427003'
            or checked['higher_order_gradients_supported'] is not False
            or checked['GPU_supported'] is not False):
        raise ValueError('The assigned actual mathematical checks must be closed')
    dna_generator = torch.Generator().manual_seed(0)
    # Match the original artificial DNA generation, then take its first window.
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (49152,),
        generator=dna_generator).tolist())[:LENGTH]
    generator = torch.Generator().manual_seed(41)
    logits = torch.randint(-3, 4, (LENGTH, 15), generator=generator).float().requires_grad_(True)
    started = dict(job_id=os.environ['SLURM_JOB_ID'], source_mathematical_checks_job=CHECKS_JOB,
        implementation='compact_grammar_specific_explicit_first_order_CPU_reverse', original_naive_probe_job='13426855', original_indexed_probe_job='13427015',
        length=LENGTH, strand='+', DNA_seed=0, DNA_generation_length=49152,
        potential_seed=41, potential_definition='first_window_CPU_FP32_integers_minus3_to3_all15_roles',
        max_exons=128, null=0, role_weights=1, score_accumulation='FP64',
        torch_version=str(torch.__version__), threads=2, synthetic_only=True,
        motif_counts=S.motif_counts(sequence), started_utc=S.utc_now())
    save(out/'started.json', started)
    print(json.dumps(dict(stage='full_window_forward_started', **started)), flush=True)
    began = time.perf_counter()
    result = reverse_full_grammar_log_partition(sequence, logits, 128)
    forward_seconds = time.perf_counter()-began
    value = result['log_partition']
    if not bool(torch.isfinite(value)) or float(value.detach()) < 0:
        raise ValueError('The complete null-inclusive partition must be finite and nonnegative')
    forward = dict(forward_seconds=forward_seconds, log_partition=float(value.detach()),
        counts=result['counts'], process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        forward_completed_utc=S.utc_now())
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
    summary = dict(**started, completed=True, completed_utc=S.utc_now(),
        protocol='full_grammar_partition_preparation.md', forward=forward,
        backward_seconds=backward_seconds, common_row_shift_gradient_max_abs=row_shift_residual,
        process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
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
