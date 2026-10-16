#!/usr/bin/env python3
"""A17 one engineering fresh restart; original compute/loss, earlier saves."""
import argparse
from collections import Counter
import json
import time
import torch
import run_query_local_context as R
import verify_query_local_context_train as V

FAILED = R.FIT
R.FIT = R.ROOT/'outputs/M28-A17-LOCAL-CONTEXT-FIT-RECOVERY-R1'
R.READOUT = R.ROOT/'outputs/M28-A17-LOCAL-CONTEXT-TRAIN-READOUT-RECOVERY-R1'
V.OUT = R.ROOT/'outputs/M28-A17-LOCAL-CONTEXT-RESULT-CHECK-RECOVERY-R1'


def fit():
    assert (FAILED/'JOBID').read_text().strip() == '13367504'
    assert not list(FAILED.glob('step_*.pt')) and not (FAILED/'summary.json').exists()
    ready = json.loads((R.PREFLIGHT/'summary.json').read_text())
    assert (R.PREFLIGHT/'STATUS').read_text().strip() == 'COMPLETED'
    assert ready['completed'] and ready['configuration'] == R.CONFIG
    assert ready['exact_A3_fresh_seed0_parameter_identity'] and ready['model_weights_unchanged']
    assert ready['optimizer_steps_this_run'] == 0
    draws, index, windows, catalogs = R.Q.load_train_inputs()
    model = R.new_model().train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step = 0; epochs = []
    with (R.FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            sums = Counter(); counts = {s: Counter() for s in R.A.ALLOW}; epoch_began = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                w = windows[draw['window_id']]
                cached = torch.load(R.ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
                features, dna = R.Q.window_inputs(cached, w, 'cuda')
                optimizer.zero_grad(set_to_none=True)
                outputs = model(features, dna, [w['valid_bases']])
                result = R.A.supervised_readout(model, outputs, w, catalogs[w['species']], cached['sequence'],
                                              draw['importance_weight_uniform_within_stratum'])
                result['total'].backward()
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1., error_if_nonfinite=True)
                optimizer.step(); step += 1
                labels = result['query']['query_labels']; values = {k: float(v.detach()) for k, v in result['losses'].items()}
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=w['id'], species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'], losses=values,
                    weighted_total=float(result['total'].detach()), gradient_norm_before_clip=float(norm),
                    query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
                    protection_reasons=result['protection_reasons'], gene_to_query=result['query']['gene_to_query'])
                log.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update(dict(values, weighted_total=row['weighted_total']))
                counts[w['species']].update({'query_'+k: v for k, v in row['query_label_counts'].items()})
                counts[w['species']].update(result['protection_reasons'])
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608)), flush=True)
                if step % 256 == 0:
                    checkpoint = R.FIT/('step_'+str(step).zfill(6)+'.pt')
                    with checkpoint.open('xb') as handle:
                        torch.save(dict(arm='A17', configuration=R.CONFIG, protocol=R.PROTOCOL, step=step,
                                        revision=R.A.F.REVISION, model=model.state_dict(), optimizer=optimizer.state_dict()), handle)
                del cached, features, dna, outputs, result
            epochs.append(dict(epoch=epoch+1, steps=1536, supervision_by_species={s: dict(c) for s, c in counts.items()},
                mean_losses={k: v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-epoch_began,
                checkpoint=str(checkpoint.relative_to(R.ROOT))))
    assert step == 4608
    R.write(R.FIT, dict(arm='A17', completed=True, feature_revision=R.A.F.REVISION, optimizer_steps=step, epochs=epochs,
        final_checkpoint=str((R.FIT/'step_004608.pt').relative_to(R.ROOT)), main_evaluation='last_step_004608_only',
        trainable_parameters=649457, supervised_checkpoint_inheritance=False, new_backbone_forward=False,
        DEV_used_for_selection=False, test_setaria_used=False, wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), checkpoint_interval=256,
        engineering_recovery='fresh_seed0_restart_after_timeout_original_compute_path',
        recovery_of_failed_job='13367504', discarded_attempt_persisted_updates=576,
        discarded_attempt_unflushed_updates_unknown=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('fit', 'readout', 'check'), required=True)
    parser.add_argument('--fit-job'); parser.add_argument('--readout-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'fit':
        fit()
    elif args.mode == 'readout':
        R.readout()
    else:
        assert args.fit_job and args.readout_job
        V.main(args.fit_job, args.readout_job)
