#!/usr/bin/env python3
"""Fixed A44 local-only or local plus per-base finite-set supervision fit."""
from collections import Counter
import json
import time
import torch

import preflight_local_structured_emission as P
from src.m28.labels import build_targets
from src.m28.local_structured_emission import position_normalized_local_structured_loss

ROOT, A = P.ROOT, P.A
PROTOCOL = 'reports/M28-METHOD-RESTART/a44_paired_fit_protocol.md'
CHECKS = ROOT/'outputs/M28-A44-LOCAL-STRUCTURED-POSITION-CHECKS'
CONFIG = {**A.CONFIG,
    'loss': 'original_local_position_mean_plus_finite_set_NLL_per_valid_base',
    'local_cross_entropy': True,
    'local_source': P.CONFIG['local_source'], 'local_targets': P.CONFIG['local_targets'],
    'local_normalizer': P.CONFIG['local_normalizer'],
    'structured_normalizer': 'actual_valid_bases_per_window_not_known_candidate_count',
    'comparison': 'same_capacity_seed0_local_only_and_frozen_A31_A3_full_grammar_TRAIN',
    'loss_coefficients_by_arm': {'A44_LOCAL': [1., 0.], 'A44_JOINT': [1., 1.]},
    'coefficient_selection': 'fixed_once_by_loss_units_no_gradient_balancing_or_search',
    'source_zero_update_preflight_job': '13425772',
    'main_checkpoint': 'last_step_004608_only'}
ARMS = ('A44_LOCAL', 'A44_JOINT')


def output_root(arm):
    if arm not in ARMS:
        raise ValueError('Only the two fixed A44 paired arms allowed')
    return ROOT/('outputs/M28-'+arm.replace('_', '-')+'-FIT')


def ready(checks_job):
    A.B.D.B.F.A.require_jobs(['13425772', checks_job])
    closed = ROOT/'outputs/M28-A44-LOCAL-STRUCTURED-PREFLIGHT-R1'
    value = A.B.D.R.completed(closed, '13425772')
    proofs = ('actual_original_R4_and_complete_support_geometry_cache_gold_reconciled',
        'actual_local_unknown_callable_padding_masks_and_position_denominators_verified',
        'actual_separate_and_joint_shared_parameter_gradients_verified',
        'actual_loss_interface_gradient_additivity_verified',
        'same_capacity_fresh_seed0_paired_initialization_verified', 'new_model_weights_unchanged')
    if (value['configuration'] != P.CONFIG or value['protocol'] != P.PROTOCOL
            or value['source_support_job'] != '13424089' or value['source_checks_job'] != '13425615'
            or any(not value[k] for k in proofs) or len(value['records']) != 8
            or value['optimizer_steps_this_run'] != 0 or value['fit_coefficients_selected']
            or value['DEV_used'] or value['test_setaria_used']):
        raise ValueError('Actual closed original local/finite gradient readiness required')
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED'
            or (CHECKS/'JOBID').read_text().strip() != checks_job):
        raise ValueError('Actual closed new position-unit checks required')


def main(arm, checks_job):
    out = output_root(arm); ready(checks_job)
    draws, index, records, support, windows, catalogs, denominators = P.population('13424089')
    model = A.A36.new_model('cuda').train(); paired = A.A36.new_model('cuda').train()
    if (sum(p.numel() for p in model.parameters()) != 366191
            or any(not torch.equal(v, paired.state_dict()[k]) for k, v in model.state_dict().items())):
        raise ValueError('Same-capacity fresh seed0 paired initialization differs')
    del paired
    params = list(model.parameters())
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    coefficients = CONFIG['loss_coefficients_by_arm'][arm]
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step = 0; epochs = []
    with (out/'training.jsonl').open('x') as handle:
        for epoch in range(3):
            started = time.perf_counter(); total = 0.; counts = Counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; rec = records[ident]; w = windows[ident]
                cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
                features, dna = P.Q.window_inputs(cached, w, 'cuda')
                targets = build_targets(w, catalogs[w['species']], cached['sequence'])
                optimizer.zero_grad(set_to_none=True)
                logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0]
                if logits.dtype != torch.float32 or features.requires_grad or dna.requires_grad:
                    raise ValueError('Original FP32 model and detached cached inputs required')
                candidates = rec['candidates'] if arm == 'A44_JOINT' else []
                result = position_normalized_local_structured_loss(logits,
                    [r['chain'] for r in candidates], [r['label'] for r in candidates], targets,
                    w['valid_bases'], local_coefficient=coefficients[0],
                    structured_coefficient=coefficients[1])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A44 paired TRAIN loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident,
                    species=w['species'], arm=arm, valid_bases=w['valid_bases'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()), local_loss=float(result['local_loss'].detach()),
                    structured_loss=float(result['structured_loss'].detach()),
                    structured_candidate_mean_loss=float(result['structured_candidate_mean_loss'].detach()),
                    log_partition=float(result['log_partition'].detach()), gold_gain=float(result['gold_gain'].detach()),
                    local_supervised_positions=result['local_supervised_positions'],
                    structured_known_candidates_used=result['known_candidates'],
                    structured_unknown_chains_excluded=result['unknown_excluded'],
                    loss_coefficients=coefficients, gradient_norm_before_clip=float(norm), **rec['counts'])
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                total += row['weighted_loss']; counts.update(rec['counts'])
                if step % 64 == 0:
                    handle.flush(); print(json.dumps(dict(arm=arm, step=step, total=4608,
                        seconds=time.perf_counter()-began)), flush=True)
                del cached, features, dna, targets, logits, result, weighted
            if dict(counts) != support['merged_replayed_draw_counts']:
                raise ValueError('A44 pass differs from unchanged original support exposure')
            state = model.state_dict()
            if any(not bool(torch.isfinite(v).all()) for v in state.values()):
                raise ValueError('Nonfinite completed-pass A44 model state')
            checkpoint = out/('step_'+str(step).zfill(6)+'.pt')
            with checkpoint.open('xb') as file:
                torch.save(dict(arm=arm, configuration=CONFIG, protocol=PROTOCOL, step=step,
                    revision=A.C.I.F.REVISION, model=state, optimizer=optimizer.state_dict()), file)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                mean_weighted_loss=total/1536, seconds=time.perf_counter()-started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 4608:
        raise ValueError('Incomplete fixed A44 original draw exposure')
    A.C.R.write_summary(out, dict(completed=True, arm=arm, configuration=CONFIG, protocol=PROTOCOL,
        source_support_job='13424089', source_preflight_job='13425772', position_checks_job=checks_job,
        feature_revision=A.C.I.F.REVISION, feature_spec=A.C.R.FEATURE_SPEC,
        optimizer_steps=step, epochs=epochs, trainable_parameters=366191,
        final_checkpoint=str((out/'step_004608.pt').relative_to(ROOT)),
        initialization='fresh_seed0_no_supervised_weight_inheritance',
        same_capacity_fresh_seed0_paired_initialization_verified=True,
        loss_coefficients=coefficients, eligible_reference_pairs=denominators,
        runtime_precision_flags=dict(cudnn_allow_tf32=torch.backends.cudnn.allow_tf32,
            matmul_allow_tf32=torch.backends.cuda.matmul.allow_tf32,
            float32_matmul_precision=torch.get_float32_matmul_precision()),
        supervised_checkpoint_read=False, original_TRAIN_catalog_metadata_read=True,
        original_reference_annotation_file_opened=False, original_R4_local_targets_regenerated=True,
        old_cached_label_fields_used=False, R4_policy_or_catalog_changed=False,
        structured_candidate_supervision_used=arm == 'A44_JOINT',
        pretrained_encoder_frozen=True, generator_forward=False, generator_updated=False,
        new_candidate_search=False, new_backbone_forward=False, DEV_used_for_selection=False,
        test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One fixed paired TRAIN arm, not full-grammar training, matched total label information '
            'or gradient size, annotation improvement, raw-assembly performance, novelty or generalization.'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', required=True, choices=ARMS)
    parser.add_argument('--checks-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.arm, args.checks_job)
