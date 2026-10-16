#!/usr/bin/env python3
"""One A43 fresh finite-support fit, fixed original draws and last checkpoint."""
from collections import Counter
import json
import time
import torch

import preflight_natural_competitor as A
from src.m28.structured_emission import structured_emission_loss

ROOT, CONFIG, PROTOCOL = A.ROOT, A.CONFIG, A.PROTOCOL
OUT = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-FIT'


def ready(support_job, preflight_job):
    A.B.D.B.F.A.require_jobs([preflight_job])
    report = A.B.D.R.completed(A.OUT, preflight_job)
    if (report['protocol'] != PROTOCOL or report['configuration'] != CONFIG or
            report['source_support_job'] != support_job or
            not report['actual_closed_support_and_original_feature_geometry_reconciled'] or
            not report['actual_new_known_gradient_nonnegative_and_unknown_gradient_zero'] or
            not report['new_model_weights_unchanged'] or not report['fresh_seed0_initialization'] or
            report['supervised_checkpoint_read'] or report['checkpoints_written'] or
            report['optimizer_steps_this_run'] != 0 or len(report['records']) != 8 or
            report['DEV_used'] or report['test_setaria_used']):
        raise ValueError('Actual matching closed zero-update union gradient readiness required')


def main(support_job, preflight_job):
    ready(support_job, preflight_job)
    draws, index, records, support = A.population(support_job)
    model = A.A36.new_model('cuda').train()
    params = list(model.parameters())
    if sum(p.numel() for p in params) != CONFIG['trainable_parameters']:
        raise ValueError('Fixed fresh representation-only model differs')
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    step = 0; epochs = []
    with (OUT/'training.jsonl').open('x') as handle:
        for epoch in range(CONFIG['proposed_passes']):
            started = time.perf_counter(); total = 0.; counts = Counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; rec = records[ident]; w = rec['window']
                cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
                features, dna = A.C.R.Q.window_inputs(cached, w, 'cuda')
                optimizer.zero_grad(set_to_none=True)
                logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0]
                if logits.dtype != torch.float32 or features.requires_grad or dna.requires_grad:
                    raise ValueError('FP32 fresh network and detached original cached inputs required')
                result = structured_emission_loss(logits, [r['chain'] for r in rec['candidates']],
                    [r['label'] for r in rec['candidates']], w['valid_bases'])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A43 finite-support TRAIN loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident,
                    species=w['species'], importance_weight=draw['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()), log_partition=float(result['log_partition'].detach()),
                    gold_gain=float(result['gold_gain'].detach()), gradient_norm_before_clip=float(norm), **rec['counts'])
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                total += row['weighted_loss']; counts.update(rec['counts'])
                if step % 64 == 0:
                    handle.flush(); print(json.dumps(dict(step=step, total=4608,
                        seconds=time.perf_counter()-began)), flush=True)
                del cached, features, dna, logits, result, weighted
            if dict(counts) != support['merged_replayed_draw_counts']:
                raise ValueError('Completed pass differs from original ordered support exposure')
            state = model.state_dict()
            if any(not bool(torch.isfinite(v).all()) for v in state.values()):
                raise ValueError('Nonfinite completed-pass A43 parameters')
            checkpoint = OUT/('step_'+str(step).zfill(6)+'.pt')
            with checkpoint.open('xb') as file:
                torch.save(dict(arm='A43', configuration=CONFIG, protocol=PROTOCOL, step=step,
                    revision=A.C.I.F.REVISION, model=state, optimizer=optimizer.state_dict()), file)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                mean_weighted_loss=total/1536, seconds=time.perf_counter()-started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != CONFIG['proposed_optimizer_steps'] or step != 4608:
        raise ValueError('Incomplete fixed A43 original draw exposure')
    A.C.R.write_summary(OUT, dict(completed=True, arm='A43', configuration=CONFIG, protocol=PROTOCOL,
        source_support_job=support_job, preflight_job=preflight_job,
        feature_revision=A.C.I.F.REVISION, feature_spec=A.C.R.FEATURE_SPEC,
        optimizer_steps=step, epochs=epochs, trainable_parameters=CONFIG['trainable_parameters'],
        final_checkpoint=str((OUT/'step_004608.pt').relative_to(ROOT)),
        main_evaluation='last_step_004608_only_full_grammar_natural_TRAIN_then_separate_raw_protocol',
        initialization='fresh_seed0_no_supervised_weight_inheritance', supervised_checkpoint_read=False,
        original_TRAIN_catalog_metadata_read=False, original_reference_annotation_file_opened=False,
        new_R4_labelling=False, pretrained_encoder_frozen=True, generator_forward=False, generator_updated=False,
        new_candidate_search=False, new_backbone_forward=False, DEV_used_for_selection=False,
        test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One finite offline-competitor TRAIN fit; not full-grammar training, checkpoint selection, '
            'natural annotation benefit, matched supervision density or independent generalization.'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-job', required=True)
    parser.add_argument('--preflight-job', required=True)
    torch.set_num_threads(2)
    args = parser.parse_args()
    main(args.support_job, args.preflight_job)
