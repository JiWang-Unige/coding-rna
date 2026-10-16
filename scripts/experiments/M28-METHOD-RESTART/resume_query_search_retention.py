#!/usr/bin/env python3
"""One A30 own-trajectory continuation after actual TIMEOUT; no fresh restart."""
import argparse
from collections import Counter
import json
import os
import shutil
import time
import torch
import inspect_search_retention_checkpoint as K
import verify_query_search_retention_train as V

R = K.R
FIT = R.ROOT/'outputs/M28-A30-SEARCH-RETENTION-FIT-RECOVERY-R1'
READOUT = R.ROOT/'outputs/M28-A30-SEARCH-RETENTION-TRAIN-READOUT-RECOVERY-R1'
CHECK = R.ROOT/'outputs/M28-A30-SEARCH-RETENTION-RESULT-CHECK-RECOVERY-R1'


def bind_outputs():
    # The source Path in K was captured before this one-off binding.
    R.FIT, R.READOUT, V.OUT = FIT, READOUT, CHECK


def accumulator():
    return dict(losses=[Counter() for _ in range(3)],
        counts=[{s: Counter() for s in R.A.ALLOW} for _ in range(3)],
        steps=[0]*3, seconds=[0.]*3, margin_steps=Counter(), total=0)


def accumulate(saved, row):
    i = saved['total']; epoch, ordinal = divmod(i, 1536)
    if (row['step'], row['epoch'], row['draw_ordinal']) != (i+1, epoch+1, ordinal):
        raise ValueError('Copied and continued rows must form one ordered 4608-step trajectory')
    species = row['species']; saved['total'] += 1; saved['steps'][epoch] += 1
    saved['losses'][epoch].update(dict(row['losses'],
        search_retention_margin=row['search_retention_margin'], weighted_total=row['weighted_total']))
    saved['counts'][epoch][species].update({'query_'+k: v for k, v in row['query_label_counts'].items()})
    saved['counts'][epoch][species].update(row['protection_reasons'])
    saved['margin_steps'][species] += int(row['search_retention_margin'] > 0)
    saved['seconds'][epoch] += row['load_forward_supervision_search_backward_optimizer_seconds']


def epochs(saved):
    if saved['total'] != 4608 or saved['steps'] != [1536]*3:
        raise ValueError('Exactly three original 1536-draw epochs required')
    return [dict(epoch=i+1, steps=1536,
        supervision_by_species={s: dict(c) for s, c in saved['counts'][i].items()},
        mean_losses={k: v/1536 for k, v in saved['losses'][i].items()}, seconds=saved['seconds'][i],
        checkpoint=str((FIT/('step_'+str((i+1)*1536).zfill(6)+'.pt')).relative_to(R.ROOT))) for i in range(3)]


def current_update(model, optimizer, draw, w, cached, catalog):
    """Unchanged original four losses, current safe search term and clip1 update."""
    started = time.perf_counter()
    features, dna = R.Q.window_inputs(cached, w, 'cuda')
    optimizer.zero_grad(set_to_none=True)
    outputs = model(features, dna, [w['valid_bases']])
    weight = draw['importance_weight_uniform_within_stratum']
    base = R.A.supervised_readout(model, outputs, w, catalog, cached['sequence'], weight)
    targets = R.Q.build_targets(w, catalog, cached['sequence'])
    window = model.query_window(outputs['features'][0], w['valid_bases'])
    t = time.perf_counter()
    extra = R.search_retention_loss(window, cached['sequence'], base['query']['positive_chains'],
                                   targets, R.CONFIG['search_retention']['margin_log_score_units'])
    torch.cuda.synchronize(); collection_seconds = time.perf_counter()-t
    total = base['total']+weight*R.CONFIG['search_retention']['coefficient']*extra['loss']
    if not bool(torch.isfinite(total)):
        raise ValueError('Nonfinite fixed A30 objective')
    total.backward()
    norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1., error_if_nonfinite=True)
    optimizer.step()
    labels = base['query']['query_labels']; values = {k: float(v.detach()) for k, v in base['losses'].items()}
    torch.cuda.synchronize()
    return dict(window_id=w['id'], species=w['species'], importance_weight=weight,
        losses=values, weighted_original_total=float(base['total'].detach()),
        search_retention_margin=float(extra['loss'].detach()), weighted_total=float(total.detach()),
        gradient_norm_before_clip=float(norm), positive_queries=extra['positive_queries'],
        collection_records=extra['records'], reference_injected_candidates=extra['reference_injected_candidates'],
        query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
        protection_reasons=base['protection_reasons'], gene_to_query=base['query']['gene_to_query'],
        current_search_collection_and_live_score_seconds=collection_seconds,
        load_forward_supervision_search_backward_optimizer_seconds=time.perf_counter()-started)


def fit():
    # Authoritative source gate before weights, TRAIN inputs, CUDA or output changes.
    checkpoint = K.latest_timeout_checkpoint()
    model, optimizer, start_step = K.restore_checkpoint(checkpoint)
    draws, index, windows, catalogs = R.Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population required')
    # Validate prefix before creating output or initializing GPU work.
    saved = accumulator()
    for row in K.validated_prefix(K.SOURCE/'training.jsonl', start_step, draws, windows):
        accumulate(saved, row)
    job = os.environ['SLURM_JOB_ID']
    bind_outputs(); FIT.mkdir(parents=True, exist_ok=False)
    (FIT/'STATUS').write_text('RUNNING\n'); (FIT/'JOBID').write_text(job+'\n')
    for step in range(256, start_step+1, 256):
        path = K.SOURCE/('step_'+str(step).zfill(6)+'.pt')
        if path.stat().st_size == 0:
            raise ValueError('Original engineering checkpoint prefix missing')
        shutil.copyfile(path, FIT/path.name)
    began = time.perf_counter()
    # Reload after moving the model: PyTorch maps moments to parameter devices
    # while retaining the original per-parameter AdamW step counters.
    moments = optimizer.state_dict(); model.to('cuda').train()
    optimizer = K.load_optimizer(model, moments, start_step); del moments
    torch.cuda.reset_peak_memory_stats()
    with (FIT/'training.jsonl').open('x') as log:
        # Copy source rows verbatim; no discarded post-checkpoint row is used.
        with (K.SOURCE/'training.jsonl').open() as source:
            for _ in range(start_step):
                log.write(source.readline())
        log.flush()
        for i in range(start_step, 4608):
            started = time.perf_counter()
            epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; w = windows[draw['window_id']]
            cached = torch.load(R.ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
            row = current_update(model, optimizer, draw, w, cached, catalogs[w['species']])
            row['load_forward_supervision_search_backward_optimizer_seconds'] = time.perf_counter()-started
            step = i+1; row.update(step=step, epoch=epoch+1, draw_ordinal=ordinal)
            log.write(json.dumps(row, allow_nan=False)+'\n'); accumulate(saved, row)
            if step % 64 == 0:
                log.flush(); print(json.dumps(dict(step=step, total=4608, resumed_from=start_step,
                    margin=row['search_retention_margin'], seconds=time.perf_counter()-began)), flush=True)
            if step % 256 == 0:
                with (FIT/('step_'+str(step).zfill(6)+'.pt')).open('xb') as handle:
                    torch.save(dict(arm='A30', configuration=R.CONFIG, protocol=R.PROTOCOL, step=step,
                        revision=R.A.F.REVISION, model=model.state_dict(), optimizer=optimizer.state_dict()), handle)
            del cached, row
    actual_epochs = epochs(saved)
    if any(saved['margin_steps'][s] == 0 for s in R.A.ALLOW):
        raise ValueError('Both-species current-search signal missing')
    R.write(FIT, dict(arm='A30', completed=True, feature_revision=R.A.F.REVISION, optimizer_steps=4608,
        epochs=actual_epochs, final_checkpoint=str((FIT/'step_004608.pt').relative_to(R.ROOT)),
        main_evaluation='last_step_004608_only', trainable_parameters=649457,
        supervised_checkpoint_inheritance=False, supervised_checkpoint_from_other_fit=False,
        own_interrupted_checkpoint_restored=True, new_backbone_forward=False,
        reference_injected_candidates=0, DEV_used_for_selection=False, test_setaria_used=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        checkpoint_interval=256, positive_search_margin_steps_by_species=dict(saved['margin_steps']),
        engineering_recovery='same_seed0_A30_own_model_and_AdamW_trajectory_continuation_after_actual_TIMEOUT',
        recovery_of_failed_job=K.SOURCE_JOB, source_terminal_TIMEOUT_verified=True,
        restored_source_checkpoint=str(checkpoint.relative_to(R.ROOT)), restored_logical_step=start_step,
        reused_source_prefix_records=start_step, optimizer_updates_this_recovery=4608-start_step,
        source_artifacts_modified=False, abandoned_source_tail_not_in_final_optimizer_trajectory=True,
        discarded_source_tail_updates_not_assumed_known=True,
        epoch_seconds_semantics='sum_saved_step_wall_seconds_including_source_prefix',
        source_RNG_state_not_saved=True, current_training_modules_have_no_enabled_dropout=True,
        limitation='Same logical 4608-update trajectory, not a fresh repeat or proof of bitwise uninterrupted GPU execution.'))
    (FIT/'STATUS').write_text('COMPLETED\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('fit', 'readout', 'check'), required=True)
    parser.add_argument('--fit-job'); parser.add_argument('--readout-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'fit':
        fit()
    else:
        bind_outputs()
        if args.mode == 'readout': R.readout()
        else:
            if not args.fit_job or not args.readout_job:
                raise ValueError('Actual recovery fit/readout job IDs required')
            V.main(args.fit_job, args.readout_job)
