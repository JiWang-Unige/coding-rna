#!/usr/bin/env python3
"""One fixed A36 TRAIN fit; actual closed preflight required, no automatic retry."""
from collections import Counter
import time
import json
import torch
import preflight_structured_emission as A
from src.m28.structured_emission import structured_emission_loss

ROOT, CONFIG, PROTOCOL = A.ROOT, A.CONFIG, A.PROTOCOL
OUT = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-FIT'
PREFLIGHT_JOB = '13404309'


def ready():
    # Computational readiness only: these artifacts do not grant permissions.
    A.require_jobs([PREFLIGHT_JOB])
    if ((A.OUT/'STATUS').read_text().strip() != 'COMPLETED' or
            (A.OUT/'JOBID').read_text().strip() != PREFLIGHT_JOB):
        raise ValueError('Actual matching preflight receipt required')
    report = json.loads((A.OUT/'summary.json').read_text())
    if (not report['completed'] or report['protocol'] != PROTOCOL or report['configuration'] != CONFIG or
            not report['new_model_weights_unchanged'] or not report['fresh_seed0_initialization'] or
            not report['all_original_gold_sets_compatible'] or report['optimizer_steps_this_run'] != 0 or
            report['supervised_checkpoint_read'] or report['checkpoints_written'] or len(report['records']) != 8):
        raise ValueError('Closed unchanged zero-update whole-chain gradient preparation required')


def main():
    ready()
    draws, index, windows, records, old_rows = A.A14.population()
    if len(draws) != 1536 or len(records) != 1491 or len(old_rows) != 16857:
        raise ValueError('Unchanged complete original TRAIN population required')
    model = A.new_model('cuda').train()
    params = list(model.parameters())
    if sum(p.numel() for p in params) != CONFIG['trainable_parameters']:
        raise ValueError('Fixed representation model changed')
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    step = 0; epochs = []
    with (OUT/'training.jsonl').open('x') as handle:
        for epoch in range(CONFIG['proposed_passes']):
            started = time.perf_counter(); total = 0.; counts = Counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; w = windows[ident]; rec = records[ident]
                cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
                features, dna = A.C.R.Q.window_inputs(cached, w, 'cuda')
                optimizer.zero_grad(set_to_none=True)
                logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0]
                if logits.dtype != torch.float32 or features.requires_grad or dna.requires_grad:
                    raise ValueError('FP32 representation and detached encoder/DNA inputs required')
                result = structured_emission_loss(logits, [r['chain'] for r in rec['candidates']],
                    [r['label'] for r in rec['candidates']], w['valid_bases'])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A36 structured TRAIN loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident,
                    species=w['species'], importance_weight=draw['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()), log_partition=float(result['log_partition'].detach()),
                    gold_gain=float(result['gold_gain'].detach()), gradient_norm_before_clip=float(norm),
                    **rec['counts'])
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                total += row['weighted_loss']; counts.update(rec['counts'])
                if step % 64 == 0:
                    handle.flush(); print(json.dumps(dict(step=step, total=4608,
                        seconds=time.perf_counter()-began)), flush=True)
                del cached, features, dna, logits, result, weighted
            state = model.state_dict()
            if any(not bool(torch.isfinite(v).all()) for v in state.values()):
                raise ValueError('Nonfinite completed-pass parameters')
            checkpoint = OUT/('step_'+str(step).zfill(6)+'.pt')
            with checkpoint.open('xb') as file:
                torch.save(dict(arm='A36', configuration=CONFIG, protocol=PROTOCOL, step=step,
                    revision=A.C.I.F.REVISION, model=state, optimizer=optimizer.state_dict()), file)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                mean_weighted_loss=total/1536, seconds=time.perf_counter()-started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != CONFIG['proposed_optimizer_steps'] or step != 4608:
        raise ValueError('Incomplete fixed A36 exposure')
    A.C.R.write_summary(OUT, dict(completed=True, arm='A36', configuration=CONFIG, protocol=PROTOCOL,
        preflight_job=PREFLIGHT_JOB, feature_revision=A.C.I.F.REVISION, feature_spec=A.C.R.FEATURE_SPEC,
        optimizer_steps=step, epochs=epochs, trainable_parameters=CONFIG['trainable_parameters'],
        final_checkpoint=str((OUT/'step_004608.pt').relative_to(ROOT)),
        main_evaluation='last_step_004608_only_same_natural_pool_TRAIN_then_separate_raw_FASTA_protocol',
        initialization='fresh_seed0_no_supervised_weight_inheritance', supervised_checkpoint_read=False,
        original_TRAIN_catalog_metadata_read=True, new_R4_labelling=False,
        pretrained_encoder_frozen=True, generator_forward=False, generator_updated=False,
        new_candidate_search=False, new_backbone_forward=False, DEV_used_for_selection=False,
        test_setaria_used=False, A33_outputs_opened=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One original-TRAIN representation fit, not model selection, free annotation, '
            'matched total exposure, independent generalization or a scientific PASS.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
