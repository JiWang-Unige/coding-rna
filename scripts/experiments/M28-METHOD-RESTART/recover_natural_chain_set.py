#!/usr/bin/env python3
"""One A15 engineering restart: unchanged science, immutable GPU hidden reuse."""
import argparse
from collections import Counter
import json
import time
import torch
import run_natural_chain_set as R
import verify_natural_chain_set_train as V

FAILED = R.FIT
R.FIT = R.ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-FIT-RECOVERY-R1'
R.READOUT = R.ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-TRAIN-READOUT-RECOVERY-R1'
V.OUT = R.ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-RESULT-CHECK-RECOVERY-R1'
TIMING = R.ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-TIMING'


def fit():
    assert (FAILED/'JOBID').read_text().strip() == '13367306'
    assert not list(FAILED.glob('step_*.pt')) and not (FAILED/'summary.json').exists()
    ready = json.loads((R.PREFLIGHT/'summary.json').read_text())
    assert (R.PREFLIGHT/'STATUS').read_text().strip() == 'COMPLETED'
    assert ready['completed'] and ready['configuration'] == R.CONFIG
    assert ready['all_natural_gold_sets_compatible'] and ready['head_weights_unchanged']
    assert ready['optimizer_steps_this_run'] == 0
    timing = json.loads((TIMING/'summary.json').read_text())
    assert (TIMING/'STATUS').read_text().strip() == 'COMPLETED'
    assert (TIMING/'JOBID').read_text().strip() == '13367598'
    assert timing['completed'] and timing['configuration'] == R.CONFIG
    assert timing['optimizer_steps_this_run'] == 0 and timing['head_weights_unchanged']
    assert timing['mean_stage_seconds']['cache_read_seconds'] > sum(
        v for k, v in timing['mean_stage_seconds'].items() if k != 'cache_read_seconds')
    began = time.perf_counter(); draws, index, windows, records, _ = R.population()
    generator = R.C.load_generator()
    required = sum(w['width']*128*4 for w in windows.values())
    assert required == timing['predicted_FP32_all_unique_hidden_bytes']
    assert required + timing['peak_GPU_allocated_bytes'] < torch.cuda.mem_get_info()[0]
    hidden_by_id = {}; cache_began = time.perf_counter()
    for ordinal, (ident, w) in enumerate(windows.items()):
        hidden = R.P.inputs(generator, index, w).detach()
        assert hidden.dtype == torch.float32 and not hidden.requires_grad
        assert list(hidden.shape) == [w['width'], 128] and bool(torch.isfinite(hidden).all())
        hidden_by_id[ident] = hidden
        if (ordinal+1) % 64 == 0:
            print(json.dumps(dict(phase='immutable_hidden', windows=ordinal+1, total=1491)), flush=True)
    cache_seconds = time.perf_counter()-cache_began
    assert len(hidden_by_id) == 1491 and sum(h.numel()*h.element_size() for h in hidden_by_id.values()) == required
    head = R.C.new_quality('cuda').train()
    params = [p for p in head.parameters() if p.requires_grad]
    assert sum(p.numel() for p in params) == 165635
    initial = {k: v.detach().clone() for k, v in head.state_dict().items()}
    witnesses = []
    for witness in R.P.select_gradient_windows(records):
        ident = witness['window_id']; w = windows[ident]; local = records[ident]['candidates']
        live = R.P.inputs(generator, index, w)
        torch.testing.assert_close(hidden_by_id[ident], live, rtol=1e-6, atol=1e-6)
        prior_gradients = None; prior_values = None; prior_scores = None
        for hidden in (live, hidden_by_id[ident]):
            head.zero_grad(set_to_none=True)
            score = head(hidden, [r['chain'] for r in local], w['valid_bases'])
            result = R.chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
            result['loss'].backward()
            norm = torch.nn.utils.clip_grad_norm_(params, 1., error_if_nonfinite=True)
            gradients = [p.grad.detach().clone() for p in params]
            values = torch.stack([result['loss'].detach(), result['log_partition'].detach(), result['gold_gain'].detach(), norm])
            if prior_gradients is not None:
                torch.testing.assert_close(score, prior_scores, rtol=1e-6, atol=1e-6)
                torch.testing.assert_close(values, prior_values, rtol=1e-6, atol=1e-6)
                for a, b in zip(gradients, prior_gradients):
                    torch.testing.assert_close(a, b, rtol=1e-6, atol=1e-6)
            prior_gradients, prior_values = gradients, values
            prior_scores = {k: v.detach().clone() for k, v in score.items()}
        witnesses.append(dict(**witness, hidden_max_abs_difference=float((live-hidden_by_id[ident]).abs().max()),
                              score_loss_clipped_gradient_equivalence=True))
        del live, score, result, gradients, prior_gradients, prior_values, prior_scores
    assert all(torch.equal(initial[k], v) for k, v in head.state_dict().items())
    assert not any(p.grad is not None for p in generator.parameters())
    assert not any(p.grad is not None for p in head.link.parameters())
    head.zero_grad(set_to_none=True)
    receipt = dict(completed=True, windows=1491, hidden_bytes=required, cache_seconds=cache_seconds,
                   live_cached_equivalence_witnesses=witnesses, optimizer_steps_this_check=0,
                   head_weights_unchanged=True, generator_frozen=True, precision='FP32',
                   source_timing_job='13367598', recovery_of_failed_job='13367306')
    with (R.FIT/'cache_receipt.json').open('x') as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False); handle.write('\n')
    optimizer = torch.optim.AdamW(params, lr=R.CONFIG['lr'], weight_decay=R.CONFIG['weight_decay'],
                                 betas=tuple(R.CONFIG['betas']), eps=R.CONFIG['eps'])
    step = 0; epochs = []; torch.cuda.reset_peak_memory_stats()
    with (R.FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts = Counter(); total = 0.; started = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; w = windows[ident]; local = records[ident]['candidates']
                optimizer.zero_grad(set_to_none=True)
                score = head(hidden_by_id[ident], [r['chain'] for r in local], w['valid_bases'])
                result = R.chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite unchanged A15 score/loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, R.CONFIG['clip_norm'], error_if_nonfinite=True)
                assert not any(p.grad is not None for p in generator.parameters())
                assert not any(p.grad is not None for p in head.link.parameters())
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident, species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'], weighted_loss=float(weighted.detach()),
                    log_partition=float(result['log_partition'].detach()), gold_gain=float(result['gold_gain'].detach()),
                    gradient_norm_before_clip=float(norm), **records[ident]['counts'])
                log.write(json.dumps(row, allow_nan=False)+'\n'); counts.update(records[ident]['counts']); total += row['weighted_loss']
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608)), flush=True)
                if step % 256 == 0:
                    with (R.FIT/('step_'+str(step).zfill(6)+'.pt')).open('xb') as handle:
                        torch.save(dict(arm='A15', configuration=R.CONFIG, protocol=R.PROTOCOL, step=step,
                            revision=R.C.I.F.REVISION, head=head.state_dict(), optimizer=optimizer.state_dict()), handle)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                               mean_weighted_loss=total/1536, seconds=time.perf_counter()-started))
    assert step == 4608
    R.C.R.write_summary(R.FIT, dict(arm='A15', completed=True, configuration=R.CONFIG, protocol=R.PROTOCOL,
        feature_revision=R.C.I.F.REVISION, feature_spec=R.C.R.FEATURE_SPEC, optimizer_steps=step,
        final_checkpoint=str((R.FIT/'step_004608.pt').relative_to(R.ROOT)), epochs=epochs,
        source_free_pool=str((R.C.FIT/'free_pool.jsonl').relative_to(R.ROOT)), source_TRAIN_cohort=str(R.D.OUT.relative_to(R.ROOT)),
        generator_trainable_parameters=0, quality_trainable_parameters=165635, DEV_used_for_selection=False,
        test_setaria_used=False, new_candidate_search=False, new_backbone_forward=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        engineering_recovery='fresh_seed0_restart_after_timeout_no_weights_inherited',
        recovery_of_failed_job='13367306', discarded_attempt_persisted_updates=896,
        discarded_attempt_unflushed_updates_unknown=True, checkpoint_interval=256,
        immutable_hidden_cache_receipt=str((R.FIT/'cache_receipt.json').relative_to(R.ROOT)),
        limitation='Same final-model 4608 updates; discarded failed-attempt compute is additional, not exact resume.'))


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
        receipt = json.loads((R.FIT/'cache_receipt.json').read_text())
        assert receipt['completed'] and receipt['windows'] == 1491 and receipt['optimizer_steps_this_check'] == 0
        assert receipt['head_weights_unchanged'] and receipt['generator_frozen']
        assert len(receipt['live_cached_equivalence_witnesses']) == 8
        assert all(w['score_loss_clipped_gradient_equivalence'] for w in receipt['live_cached_equivalence_witnesses'])
        V.main(args.fit_job, args.readout_job)
