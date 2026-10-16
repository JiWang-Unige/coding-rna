#!/usr/bin/env python3
"""Verify the completed fixed fit and summarize online TRAIN losses, no updates."""
from collections import Counter
import json
import math
import os
import statistics
import time

import torch
import fit_full_grammar_coarse_channel as F

OUT = F.ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIT-VERIFICATION-R1'
LOSS_KEYS = ('raw_NLL', 'normalized_NLL', 'weighted_total')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def write(name, value):
    with (OUT/name).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def losses(rows):
    return dict(draws=len(rows), mean={k:statistics.mean(r[k] for r in rows) for k in LOSS_KEYS})


def main():
    began = time.perf_counter()
    job = os.environ['SLURM_JOB_ID']
    write('started.json', dict(job_id=job, source_fit_job='13434603', configuration=F.CONFIG,
        scope='completed_original_fit_logs_and_three_checkpoint_states_only', optimizer_steps_this_run=0))
    summary = json.loads((F.FIT/'summary.json').read_text())
    started = json.loads((F.FIT/'started.json').read_text())
    require((F.FIT/'STATUS').read_text().strip() == 'COMPLETED'
        and (F.FIT/'JOBID').read_text().strip() == started['job_id'] == summary['job_id'] == '13434603'
        and summary['completed'] and summary['configuration'] == started['configuration'] == F.CONFIG
        and summary['protocol'] == F.PROTOCOL and summary['preparation_job'] == '13429756',
        'Completed unchanged fixed fit identity required')
    stdout = F.ROOT/'outputs/M28-PILOT-R5/logs/M28COARSEFIRSTFIT_13434603.out'
    require(json.loads(stdout.read_text().splitlines()[-1]) == summary
        and stdout.with_suffix('.err').stat().st_size == 0, 'Fit stdout/summary mismatch or stderr')
    _, _, windows, _, plan, _ = F.original_plan()
    require(plan == [json.loads(x) for x in (F.READY/'training_plan.jsonl').read_text().splitlines()],
        'Original preparation draw plan drifted')
    rows = [json.loads(x) for x in (F.FIT/'training.jsonl').read_text().splitlines()]
    require(len(rows) == 4608 and len(windows) == 1491, 'Full three-pass exposure required')
    numeric = LOSS_KEYS + ('gradient_norm_before_clip', 'seconds')
    unknown = []
    for ordinal, row in enumerate(rows):
        expected = plan[ordinal % 1536]
        require(row['step'] == ordinal+1 and row['epoch'] == ordinal//1536+1,
            'Noncontinuous optimizer step or epoch')
        for key in ('draw_ordinal', 'window_id', 'species', 'importance_weight',
                    'known_observation_denominator', 'fine_positions', 'coarse_positions', 'fully_unknown_positions'):
            require(row[key] == expected[key], 'Original draw/weight/observation drift: '+key)
        require(all(math.isfinite(row[k]) for k in numeric) and row['raw_NLL'] >= -1e-8
            and row['gradient_norm_before_clip'] >= 0 and row['seconds'] >= 0,
            'Nonfinite or invalid training numeric record')
        require(math.isclose(row['normalized_NLL'], row['raw_NLL']/expected['known_observation_denominator'],
                rel_tol=1e-12, abs_tol=1e-12)
            and math.isclose(row['weighted_total'], row['raw_NLL']*expected['objective_multiplier'],
                rel_tol=1e-12, abs_tol=1e-12), 'Recorded loss normalization or importance mismatch')
        if expected['fine_positions']+expected['coarse_positions'] == 0:
            require(abs(row['raw_NLL']) <= 1e-8 and row['gradient_norm_before_clip'] <= 1e-6,
                'All-unknown exposure acquired an observation gradient')
            unknown.append({k:row[k] for k in ('step', 'draw_ordinal', 'window_id', 'raw_NLL', 'gradient_norm_before_clip')})
    epochs = []
    for epoch in range(1, 4):
        part = rows[(epoch-1)*1536:epoch*1536]
        stats = losses(part)
        for key in LOSS_KEYS:
            require(math.isclose(stats['mean'][key], summary['epochs'][epoch-1]['mean_losses'][key],
                rel_tol=1e-12, abs_tol=1e-12), 'Summary epoch loss differs from complete records')
        species = {s:losses([r for r in part if r['species'] == s]) for s in sorted(F.Q.ALLOW)}
        require(dict(Counter(r['species'] for r in part)) == summary['epochs'][epoch-1]['by_species'],
            'Species exposure mismatch')
        strata = {s:losses([r for r, p in zip(part, plan) if json.dumps(p['stratum']) == s])
            for s in sorted({json.dumps(p['stratum']) for p in plan})}
        negatives = [r for r, p in zip(part, plan) if not p['original_hard_event_support_feasible']]
        epochs.append(dict(epoch=epoch, **stats, by_species=species, by_stratum=strata,
            original_hard_negative_draws=losses(negatives)))

    torch.manual_seed(0)
    model = F.G.M28Core('C0', 2048, 128)
    initial = {k:v.detach().clone() for k, v in model.state_dict().items()}
    require(sum(p.numel() for p in model.parameters()) == 366191, 'Fresh C0 capacity drifted')
    previous = initial
    checkpoints = []
    for epoch in range(1, 4):
        step = epoch*1536
        path = F.FIT/('step_%06d.pt' % step)
        saved = torch.load(path, map_location='cpu', weights_only=True)
        require(saved['configuration'] == F.CONFIG and saved['protocol'] == F.PROTOCOL
            and saved['encoder_revision'] == F.Q.F.REVISION and saved['step'] == step
            and saved['epoch'] == epoch, 'Checkpoint configuration or step mismatch')
        model.load_state_dict(saved['model'], strict=True)
        state = saved['model']
        require(all(v.dtype == torch.float32 and torch.isfinite(v).all() for v in state.values()),
            'Checkpoint has nonfinite or non-FP32 C0 weights')
        changed = [k for k in state if not torch.equal(state[k], initial[k])]
        between = [k for k in state if not torch.equal(state[k], previous[k])]
        require(changed and between, 'Fresh C0 did not change or epoch updates made no weight change')
        parameters = list(model.parameters())
        opt = saved['optimizer']
        require(len(opt['param_groups']) == 1, 'Original single AdamW parameter group required')
        group = opt['param_groups'][0]
        for key, value in dict(lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8).items():
            require(group[key] == value, 'Optimizer hyperparameter drift: '+key)
        identifiers = group['params']
        require(len(identifiers) == len(parameters) and len(set(identifiers)) == len(parameters)
            and set(opt['state']) == set(identifiers), 'Optimizer parameter/state coverage mismatch')
        for ident, parameter in zip(identifiers, parameters):
            item = opt['state'][ident]
            require(float(item['step']) == step, 'AdamW step does not match all original exposures')
            for key in ('exp_avg', 'exp_avg_sq'):
                value = item[key]
                require(value.shape == parameter.shape and value.dtype == torch.float32
                    and torch.isfinite(value).all(), 'Invalid optimizer moment: '+key)
            require((item['exp_avg_sq'] >= 0).all(), 'Negative second moment')
        delta = math.sqrt(sum(float((state[k].double()-initial[k].double()).square().sum()) for k in state))
        checkpoints.append(dict(epoch=epoch, step=step, file_bytes=path.stat().st_size,
            state_tensors=len(state), changed_tensors_from_seed0=len(changed), changed_tensors_from_previous=len(between),
            parameter_change_L2_from_seed0=delta, finite_FP32_weights_and_optimizer_moments=True,
            all_optimizer_parameter_steps=step, optimizer_parameter_states=len(identifiers), strict_C0_load=True))
        previous = state

    before, after = rows[:1536], rows[3072:]
    pairs = [(a, b) for a, b in zip(before, after) if a['fine_positions']+a['coarse_positions'] > 0]
    improved = sum(b['weighted_total'] < a['weighted_total']-1e-8 for a, b in pairs)
    worsened = sum(b['weighted_total'] > a['weighted_total']+1e-8 for a, b in pairs)
    sizes = {p.name:p.stat().st_size for p in F.FIT.iterdir() if p.is_file()}
    require(sum(sizes.values()) <= F.CONFIG['maximum_new_fit_bytes'], 'Fit output ceiling exceeded')
    result = dict(completed=True, job_id=job, source_fit_job='13434603', configuration=F.CONFIG,
        original_draw_plan_and_all_numeric_records_verified=True, training_rows=len(rows),
        original_unique_windows=len(windows), fit_artifact_bytes=sizes, epochs=epochs, checkpoints=checkpoints,
        no_known_observation_records=unknown,
        online_weighted_mean_reduction_first_to_third=1-epochs[-1]['mean']['weighted_total']/epochs[0]['mean']['weighted_total'],
        paired_online_known_draws=dict(count=len(pairs), improved_at_absolute_1e_minus8=improved,
            worsened_at_absolute_1e_minus8=worsened, unchanged_at_absolute_1e_minus8=len(pairs)-improved-worsened),
        final_checkpoint=summary['final_checkpoint'], checkpoint_selection_changed=False,
        actual_research_DNA_or_features_read_this_run=False, native_GLM_forward_this_run=False,
        optimizer_steps_this_run=0, DEV_used=False, independent_species_used=False,
        seconds=time.perf_counter()-began,
        interpretation='completed_numerically_valid_original_TRAIN_online_fit_not_fixed_checkpoint_evaluation_or_accuracy')
    write('summary.json', result)
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
