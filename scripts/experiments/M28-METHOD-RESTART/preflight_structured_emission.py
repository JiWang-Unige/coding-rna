#!/usr/bin/env python3
"""A36 fresh representation + whole-chain NLL wiring, zero optimizer updates."""
from collections import Counter
import math
import time
import torch
import run_chain_set as A14
from diagnose_missing_path_emission_train import require_jobs
from src.m28.core import M28Core
from src.m28.local_path_score import local_path_scores
from src.m28.structured_emission import trainable_path_scores, structured_emission_loss

P, C, ROOT = A14.P, A14.C, A14.ROOT
OUT = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a36_structured_emission_protocol.md'
CONFIG = dict(feature_dim=2048, hidden=128, seed=0, precision='FP32',
    model='fresh_M28Core_C0_SharedContext_and_15_class_segmentation_only',
    trainable_parameters=366191, pretrained_encoder_frozen=True,
    source_natural_generator=C.GENERATOR, generator_optimizer_steps=4608,
    candidate_pool='original_A12_16321_natural_plus_536_TRAIN_positive_injections',
    loss='exact_compatible_known_candidate_set_NLL_divided_by_known_count',
    score='unit_role_logit_minus_intergenic_logit_sum_FP64', null=0.,
    threshold=0., threshold_scan=False, local_cross_entropy=False,
    loss_unknown='exclude_-1_from_partition_and_gold',
    draw_order_and_importance_weights='original_1536_draws_1491_unique_TRAIN',
    proposed_passes=3, proposed_optimizer_steps=4608,
    optimizer='AdamW', lr=3e-4, weight_decay=.01, betas=[.9,.999], eps=1e-8, clip_norm=1.,
    proposed_main_checkpoint='last_step_004608_only',
    comparison='frozen_A31_A3_same_natural_pool_and_same_null0')


def new_model(device):
    torch.manual_seed(CONFIG['seed'])
    return M28Core('C0', CONFIG['feature_dim'], CONFIG['hidden']).to(device)


def main():
    require_jobs(['13403916', '13403920', '13404094'])
    # Reuse the original full-cohort identity/label and compatible-gold checks.
    draws, index, windows, records, old_rows = A14.population()
    if (len(draws) != 1536 or len(records) != 1491 or len(old_rows) != 16857 or
        Counter(r['origin'] for r in old_rows) != {'free': 16321, 'reference_injected': 536} or
        Counter(r['label'] for r in old_rows) != {1: 2916, 0: 12566, -1: 1375}):
        raise ValueError('Only complete unchanged original TRAIN candidate supervision is allowed')
    witnesses = P.select_gradient_windows(records)
    model = new_model('cuda').train()
    if sum(p.numel() for p in model.parameters()) != CONFIG['trainable_parameters']:
        raise ValueError('Fresh representation-only parameter count changed')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    weight = {d['window_id']: d['importance_weight_uniform_within_stratum'] for d in draws}
    receipts = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    for witness in witnesses:
        ident = witness['window_id']; w = windows[ident]; local = records[ident]['candidates']
        started = time.perf_counter()
        cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        features, dna = C.R.Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Cached encoder features and raw one-hot DNA must be detached')
        model.zero_grad(set_to_none=True)
        logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0]
        if logits.dtype != torch.float32:
            raise ValueError('FP32 network required')
        chains, labels = [r['chain'] for r in local], [r['label'] for r in local]
        gain = trainable_path_scores(logits, chains, w['valid_bases'])
        prior = local_path_scores(logits, chains, w['valid_bases'])
        torch.testing.assert_close(gain, prior, rtol=0, atol=0)
        result = structured_emission_loss(logits, chains, labels, w['valid_bases'])
        weighted = result['loss']*weight[ident]
        if not bool(torch.isfinite(weighted)):
            raise ValueError('Nonfinite actual TRAIN structured loss')
        weighted.backward()
        groups = {}
        for name, module in (('project', model.shared.project), ('context', model.shared.context),
            ('dna', model.shared.dna), ('fuse', model.shared.fuse), ('segmentation', model.segmentation)):
            params = list(module.parameters())
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite actual gradient: '+name)
            groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if groups[name] <= 0:
                raise ValueError('Zero actual gradient: '+name)
        torch.cuda.synchronize()
        receipts.append(dict(**witness, weighted_loss=float(weighted.detach()),
            gradient_norms=groups, original_A31_potential_error=float((gain-prior).abs().max()),
            load_forward_backward_seconds=time.perf_counter()-started,
            **{k: v for k, v in result.items() if not torch.is_tensor(v)}))
        print(receipts[-1], flush=True)
        del cached, features, dna, logits, gain, prior, result, weighted
    if any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()):
        raise ValueError('Zero-update preflight changed model state')
    replay = Counter()
    for draw in draws:
        replay.update(records[draw['window_id']]['counts'])
    C.R.write_summary(OUT, dict(completed=True, experiment='M28-A36-STRUCTURED-EMISSION-PREFLIGHT',
        protocol=PROTOCOL, configuration=CONFIG, windows=1491, original_draws=1536,
        candidate_window_pairs=16857, natural_candidates=16321, reference_injections=536,
        replayed_original_draw_counts=dict(replay), all_original_gold_sets_compatible=True,
        new_model_weights_unchanged=True, fresh_seed0_initialization=True,
        supervised_checkpoint_read=False, generator_forward=False,
        original_TRAIN_catalog_metadata_read=True, new_R4_labelling=False,
        optimizer_steps_this_run=0, checkpoints_written=False, new_candidate_search=False,
        new_backbone_forward=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        records=receipts, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Zero-update whole-chain gradient/identity checks only; no trained new model, '
            'natural prediction improvement, novelty claim or generalization.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
