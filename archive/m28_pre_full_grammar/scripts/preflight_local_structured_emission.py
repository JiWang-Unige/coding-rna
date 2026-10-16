#!/usr/bin/env python3
"""Original R4 local/finite supervision through a fresh model, zero updates."""
from collections import Counter
import json
import math
import time
from pathlib import Path
import torch

import preflight_natural_competitor as A
from src.m28.labels import build_targets
from src.m28.local_structured_emission import (
    local_structured_emission_loss, reconcile_local_supervision)

ROOT, Q = A.ROOT, A.C.R.Q
OUT = ROOT/'outputs/M28-A44-LOCAL-STRUCTURED-PREFLIGHT'
CHECKS = ROOT/'outputs/M28-A44-LOCAL-STRUCTURED-SOURCE-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/a44_local_structured_protocol.md'
CONFIG = dict(model=A.CONFIG['model'], trainable_parameters=366191,
    feature_dim=2048, hidden=128, seed=0, precision='FP32',
    population='original_1491_unique_TRAIN_windows_1536_unchanged_draws',
    local_source='original_R4_TRAIN_windows_and_all_overlapping_transcripts',
    local_targets='original_build_targets_from_original_oriented_cached_DNA_not_old_cache_labels',
    local_normalizer='actual_fine_or_grouped_region_supervised_positions_per_window',
    structured_source_configuration=A.CONFIG['source_support_configuration'],
    structured_normalizer='original_compatible_known_candidate_count_per_window',
    original_importance='multiply_total_once',
    probe_coefficients=[1., 1.], probe_coefficients_are_fit_selection=False,
    optimizer_steps=0, pretrained_encoder_frozen=True,
    unknown='preserve_original_local_masks_and_separate_minus1_chain_exclusion',
    source_supervised_weights=False, new_R4_policy=False)


def population(support_job):
    draws, index, records, support = A.population(support_job)
    local_draws, local_index, windows, catalogs = Q.load_train_inputs()
    if local_draws != draws or set(windows) != set(records):
        raise ValueError('Original local and finite TRAIN draw exposure differs')
    counts = {s: Counter() for s in A.C.R.ALLOW}
    for ident, rec in records.items():
        w = windows[ident]
        reconcile_local_supervision(rec, w, catalogs[w['species']],
            index[ident]['artifact'], local_index[ident]['artifact'])
        counts[w['species']].update(references=len(w['positive_ids']),
            long_references=sum(len(catalogs[w['species']][tid]['CDS']) >= 5 for tid in w['positive_ids']))
    denominators = {s: dict(v) for s, v in counts.items()}
    if denominators != {
            'arabidopsis_thaliana': {'references': 2172, 'long_references': 854},
            'oryza_sativa': {'references': 744, 'long_references': 232}}:
        raise ValueError('Complete original gold and long-chain denominators differ')
    return draws, index, records, support, windows, catalogs, denominators


def gradient_norms(model, gradients):
    values = dict(zip((name for name, _ in model.named_parameters()), gradients))
    groups = {}
    for group in ('shared.project', 'shared.context', 'shared.dna', 'shared.fuse', 'segmentation'):
        local = [v for name, v in values.items() if name.startswith(group+'.')]
        if not local or any(not bool(torch.isfinite(v).all()) for v in local):
            raise ValueError('Missing or nonfinite actual supervision gradient: '+group)
        groups[group] = math.sqrt(sum(float(v.double().square().sum()) for v in local))
        if groups[group] <= 0:
            raise ValueError('Zero actual supervision gradient: '+group)
    return groups


def main(support_job, checks_job, out=OUT):
    A.B.D.B.F.A.require_jobs([checks_job])
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED'
            or (CHECKS/'JOBID').read_text().strip() != checks_job):
        raise ValueError('Actual original-source join checks required')
    draws, index, records, support, windows, catalogs, denominators = population(support_job)
    chosen = A.witnesses(records)
    model = A.A36.new_model('cuda').train()
    paired = A.A36.new_model('cuda').train()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if (sum(p.numel() for p in model.parameters()) != CONFIG['trainable_parameters']
            or any(not torch.equal(v, paired.state_dict()[k]) for k, v in before.items())):
        raise ValueError('Same-capacity fresh seed0 paired initialization differs')
    del paired
    params = list(model.parameters())
    weights = {d['window_id']: d['importance_weight_uniform_within_stratum'] for d in draws}
    receipts = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    for witness in chosen:
        started = time.perf_counter(); ident = witness['window_id']
        w = windows[ident]; rec = records[ident]
        cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        features, dna = Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Original cached inputs must be detached')
        targets = build_targets(w, catalogs[w['species']], cached['sequence'])
        logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0]
        if logits.dtype != torch.float32:
            raise ValueError('Original FP32 representation required')
        result = local_structured_emission_loss(logits,
            [r['chain'] for r in rec['candidates']], [r['label'] for r in rec['candidates']],
            targets, w['valid_bases'], local_coefficient=1., structured_coefficient=1.)
        weighted = {key: result[key]*weights[ident]
                    for key in ('local_loss', 'structured_loss', 'loss')}
        if any(not bool(torch.isfinite(v)) for v in weighted.values()):
            raise ValueError('Nonfinite actual local or finite supervision loss')
        fine = torch.as_tensor(targets['seg'], device=logits.device)
        region = torch.as_tensor(targets['region'], device=logits.device)
        fine_known = fine >= 0
        grouped = (fine < 0) & (region >= 0) & (region < 3)
        supervised = fine_known | grouped
        if int(supervised.sum()) != result['local_supervised_positions']:
            raise ValueError('Actual original local position denominator differs')
        emission_gradient, = torch.autograd.grad(weighted['local_loss'], logits, retain_graph=True)
        if (not bool(torch.isfinite(emission_gradient).all())
                or bool((emission_gradient[~supervised] != 0).any())):
            raise ValueError('Original unknown/N/padded local positions have a gradient')
        structured_emission_gradient, = torch.autograd.grad(weighted['structured_loss'], logits, retain_graph=True)
        joint_emission_gradient, = torch.autograd.grad(weighted['loss'], logits, retain_graph=True)
        # Additivity belongs at the newly composed loss interface. Replaying
        # FP32/TF32 upstream kernels three times is not bitwise linear, even
        # for a mathematically linear vector-Jacobian product (13425702).
        expected_emission_gradient = emission_gradient+structured_emission_gradient
        torch.testing.assert_close(joint_emission_gradient, expected_emission_gradient,
                                   rtol=1e-4, atol=1e-6)
        emission_error = float((joint_emission_gradient-expected_emission_gradient).abs().max())
        gradients = {}
        for key in ('local_loss', 'structured_loss', 'loss'):
            gradients[key] = torch.autograd.grad(weighted[key], params, retain_graph=key != 'loss')
        max_error = 0.
        for local, structured, total in zip(gradients['local_loss'], gradients['structured_loss'], gradients['loss']):
            summed = local+structured
            max_error = max(max_error, float((total-summed).abs().max()))
        norms = {k: gradient_norms(model, g) for k, g in gradients.items()}
        magnitude = {k: math.sqrt(sum(float(v.double().square().sum()) for v in g))
                     for k, g in gradients.items()}
        dot = sum(float((a.double()*b.double()).sum())
                  for a, b in zip(gradients['local_loss'], gradients['structured_loss']))
        if any(p.grad is not None for p in params):
            raise ValueError('Zero-update probe unexpectedly accumulated parameter gradients')
        torch.cuda.synchronize()
        receipts.append(dict(**witness, original_importance=weights[ident],
            weighted_losses={k: float(v.detach()) for k, v in weighted.items()},
            fine_positions=int(fine_known.sum()), grouped_region_positions=int(grouped.sum()),
            local_unknown_or_uncallable_or_padded_positions=int((~supervised).sum()),
            local_supervised_positions=result['local_supervised_positions'],
            known_candidates=result['known_candidates'], unknown_chains_excluded=result['unknown_excluded'],
            gradient_norms=norms, total_gradient_norms=magnitude,
            structured_to_local_gradient_norm_ratio=magnitude['structured_loss']/magnitude['local_loss'],
            local_structured_gradient_cosine=dot/(magnitude['structured_loss']*magnitude['local_loss']),
            maximum_loss_interface_gradient_additivity_error=emission_error,
            maximum_upstream_gradient_additivity_error_diagnostic=max_error,
            masked_local_emission_gradient_zero=True,
            load_forward_three_gradient_seconds=time.perf_counter()-started))
        print(json.dumps(receipts[-1], allow_nan=False), flush=True)
        del cached, features, dna, targets, logits, result, weighted, gradients, emission_gradient
        del structured_emission_gradient, joint_emission_gradient, expected_emission_gradient
    if any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()):
        raise ValueError('Zero-update local/structured preflight changed model state')
    report = dict(completed=True, experiment='M28-A44-LOCAL-STRUCTURED-PREFLIGHT',
        protocol=PROTOCOL, configuration=CONFIG, source_support_job=support_job, source_checks_job=checks_job,
        feature_revision=A.C.I.F.REVISION, windows=1491, original_draws=1536,
        runtime_precision_flags=dict(cudnn_allow_tf32=torch.backends.cudnn.allow_tf32,
            matmul_allow_tf32=torch.backends.cuda.matmul.allow_tf32,
            float32_matmul_precision=torch.get_float32_matmul_precision()),
        candidate_window_pairs=19189, eligible_reference_pairs=denominators,
        replayed_original_draw_counts=support['merged_replayed_draw_counts'], records=receipts,
        actual_original_R4_and_complete_support_geometry_cache_gold_reconciled=True,
        actual_local_unknown_callable_padding_masks_and_position_denominators_verified=True,
        actual_separate_and_joint_shared_parameter_gradients_verified=True,
        actual_loss_interface_gradient_additivity_verified=True,
        upstream_parameter_additivity_is_diagnostic_not_pass_gate=True,
        original_parameter_guard_failed_job='13425634', arithmetic_diagnosis_job='13425702',
        same_capacity_fresh_seed0_paired_initialization_verified=True,
        new_model_weights_unchanged=True, optimizer_steps_this_run=0, optimizer_constructed=False,
        checkpoints_written=False, supervised_checkpoint_read=False, original_TRAIN_catalog_metadata_read=True,
        original_reference_annotation_file_opened=False, original_R4_local_targets_regenerated=True,
        old_cached_label_fields_used=False, new_R4_policy_or_chain_labels=False,
        new_candidate_search=False, new_backbone_forward=False, DEV_used=False,
        test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        fit_coefficients_selected=False, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Eight original TRAIN zero-update wiring witnesses only; 1/1 is an explicit probe, '
            'not fit coefficient selection, model benefit, matched gradient scale, novelty or generalization.')
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'records'}, allow_nan=False), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-job', required=True)
    parser.add_argument('--checks-job', required=True)
    parser.add_argument('--output-root', type=Path, default=OUT)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.support_job, args.checks_job, ROOT/args.output_root)
