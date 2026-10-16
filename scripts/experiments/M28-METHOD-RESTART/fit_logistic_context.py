#!/usr/bin/env python3
"""One fixed A38 natural-domain joint-context fit, no model-selection retry."""
from collections import Counter
import argparse
import json
import time
import torch
import preflight_logistic_context as A
from src.m28.core import chain_loss

ROOT, CONFIG, PROTOCOL = A.ROOT, A.CONFIG, A.PROTOCOL
OUT = ROOT/'outputs/M28-A38-LOGISTIC-CONTEXT-FIT'
PREFLIGHT_JOB = '13412601'


def ready():
    A.require_jobs([PREFLIGHT_JOB])
    report = A.B.R.completed(A.OUT, PREFLIGHT_JOB)
    if (report['protocol'] != PROTOCOL or report['configuration'] != CONFIG or
            report['windows'] != 1491 or report['original_draws'] != 1536 or
            report['natural_candidates'] != 16321 or len(report['records']) != 8 or
            report['eligible_reference_pairs'] != {'eligible':2916, 'five_or_more':1086} or
            report['optimizer_steps_this_run'] != 0 or report['reference_injected_training_candidates'] != 0 or
            any(report[k] is not True for k in ('fresh_weights_unchanged',
                'all_natural_gold_sets_compatible', 'both_species_known_junction_gradient_coverage',
                'actual_mean_BCE_gain_gradient_verified', 'initialization_identical_constructor_and_seed_to_A37',
                'original_natural_identity_labels_and_stronger_baseline_reconciled')) or
            any(report[k] is not False for k in ('supervised_checkpoint_read', 'generator_forward',
                'checkpoints_written', 'new_candidate_search', 'new_backbone_forward',
                'new_R4_labelling', 'DEV_used', 'test_setaria_used', 'A33_outputs_opened'))):
        raise ValueError('Actual unchanged zero-update joint-context preparation required')


def main():
    ready()
    draws, index, windows, records, _, _, _ = A.population()
    model = A.new_model('cuda').train(); params = list(model.parameters())
    if sum(p.numel() for p in params) != CONFIG['trainable_parameters']:
        raise ValueError('Fixed fresh joint-context structure differs')
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step = 0; epochs = []
    with (OUT/'training.jsonl').open('x') as handle:
        for epoch in range(CONFIG['proposed_passes']):
            started = time.perf_counter(); counts = Counter(); total = 0.
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; w = windows[ident]; rec = records[ident]
                cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
                features, dna = A.R.C.R.Q.window_inputs(cached, w, 'cuda')
                if features.requires_grad or dna.requires_grad:
                    raise ValueError('Original encoder features and DNA must remain detached')
                optimizer.zero_grad(set_to_none=True)
                score = model(features, dna, [w['valid_bases']], [r['chain'] for r in rec['candidates']])
                if any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
                    raise ValueError('Finite original FP32 joint scores required')
                loss = chain_loss(score, [r['label'] for r in rec['candidates']])
                weighted = loss*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A38 actual mean-known BCE')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident,
                    species=w['species'], importance_weight=draw['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()), unweighted_mean_known_BCE=float(loss.detach()), gradient_norm_before_clip=float(norm), **rec['counts'])
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                total += row['weighted_loss']; counts.update(rec['counts'])
                if step % 64 == 0:
                    handle.flush(); print(json.dumps(dict(step=step, total=4608,
                        seconds=time.perf_counter()-began)), flush=True)
                del cached, features, dna, score, loss, weighted
            state = model.state_dict()
            if any(not bool(torch.isfinite(v).all()) for v in state.values()):
                raise ValueError('Nonfinite completed-pass parameters')
            checkpoint = OUT/('step_'+str(step).zfill(6)+'.pt')
            with checkpoint.open('xb') as file:
                torch.save(dict(arm='A38', configuration=CONFIG, protocol=PROTOCOL, step=step,
                    revision=A.R.C.I.F.REVISION, model=state, optimizer=optimizer.state_dict()), file)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                mean_weighted_loss=total/1536, seconds=time.perf_counter()-started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 4608 or step != CONFIG['proposed_optimizer_steps']:
        raise ValueError('Incomplete fixed A38 exposure')
    exposure = {name: int(optimizer.state[p]['step']) if p in optimizer.state else 0
                for name, p in model.named_parameters()}
    A.R.C.R.write_summary(OUT, dict(completed=True, arm='A38', configuration=CONFIG, protocol=PROTOCOL,
        preflight_job=PREFLIGHT_JOB, feature_revision=A.R.C.I.F.REVISION,
        feature_spec=A.R.C.R.FEATURE_SPEC, optimizer_steps=step, optimizer_steps_per_parameter=exposure,
        epochs=epochs, trainable_parameters=CONFIG['trainable_parameters'],
        final_checkpoint=str((OUT/'step_004608.pt').relative_to(ROOT)),
        main_evaluation='last_step_004608_only_same_natural_TRAIN_pool_stronger_A31_A3',
        initialization=CONFIG['initialization'], supervised_checkpoint_read=False,
        original_TRAIN_catalog_metadata_read=True, pretrained_encoder_frozen=True,
        reference_injected_training_candidates=0, generator_forward=False, generator_updated=False,
        new_candidate_search=False, new_backbone_forward=False, new_R4_labelling=False,
        DEV_used_for_selection=False, test_setaria_used=False, A33_outputs_opened=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One fixed joint-context TRAIN fit only; no natural readout gate, '
            'matched-capacity causal result, raw annotation, independent generalization or scientific PASS.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight-job', required=True)
    args = parser.parse_args()
    if args.preflight_job != PREFLIGHT_JOB:
        parser.error('Only the actual declared A38 preflight job is allowed')
    torch.set_num_threads(2); main()
