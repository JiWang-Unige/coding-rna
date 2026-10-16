#!/usr/bin/env python3
"""Reconcile original exposures, logged objective arithmetic and checkpoints."""
import argparse
from collections import Counter
import json
import math
import os
import time
import torch
import fit_full_grammar_objective_contrast as C
import annotate_full_grammar_objective_contrast as A

OUT = C.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-VERIFICATION'


def require(value, message):
    if not value:
        raise ValueError(message)


def main(arm, fit_job):
    began = time.perf_counter()
    fit, out = C.OUT/arm, OUT/arm
    summary = json.loads((fit/'summary.json').read_text())
    started = json.loads((fit/'started.json').read_text())
    configuration = C.arm_config(arm)
    require((fit/'STATUS').read_text().strip() == 'COMPLETED'
        and started['job_id'] == summary['job_id'] == fit_job and summary['completed']
        and started['configuration'] == summary['configuration'] == configuration
        and summary['protocol'] == C.PROTOCOL and summary['preflight_job'] == '13451086'
        and summary['optimizer_steps_this_run'] == 4608 and summary['cumulative_optimizer_steps'] == 9216
        and not summary['DEV_used'] and not summary['sealed_data_used'], 'Completed declared continuation required')
    plan = [json.loads(x) for x in (C.F.READY/'training_plan.jsonl').read_text().splitlines()]
    require(len(plan) == 1536 and len({x['window_id'] for x in plan}) == 1491,
            'Original exposure plan size differs')
    rows = [json.loads(x) for x in (fit/'training.jsonl').read_text().splitlines()]
    require(len(rows) == 4608, 'Every additional optimizer update must be recorded')
    component_keys = {'structured_NLL'} if arm == C.ARMS[0] else {'local_NLL'} if arm == C.ARMS[1] else {'structured_NLL', 'local_NLL'}
    for i, row in enumerate(rows):
        expected = plan[i % 1536]
        require((row['step'], row['additional_step'], row['epoch']) == (4609+i, 1+i, 4+i//1536),
                'Step or epoch discontinuity')
        for key in ('draw_ordinal', 'window_id', 'species', 'importance_weight', 'known_observation_denominator',
                    'fine_positions', 'coarse_positions', 'fully_unknown_positions'):
            require(row[key] == expected[key], 'Original exposure differs: '+key)
        require(set(row['components']) == component_keys, 'Wrong training objective components')
        values = [row[k] for k in ('raw_NLL', 'normalized_NLL', 'weighted_total', 'gradient_norm_before_clip', 'seconds')]
        values.extend(row['components'].values())
        require(all(math.isfinite(v) for v in values) and min(values) >= -1e-8, 'Invalid logged numeric result')
        raw = sum(row['components'].values())*(.5 if arm == C.ARMS[2] else 1.)
        require(math.isclose(row['raw_NLL'], raw, rel_tol=1e-12, abs_tol=1e-10)
            and math.isclose(row['normalized_NLL'], raw/expected['known_observation_denominator'], rel_tol=1e-12, abs_tol=1e-12)
            and math.isclose(row['weighted_total'], raw*expected['objective_multiplier'], rel_tol=1e-12, abs_tol=1e-12),
            'Objective mixture, normalization or importance differs')
        if expected['fine_positions']+expected['coarse_positions'] == 0:
            require(abs(raw) <= 1e-8 and row['gradient_norm_before_clip'] <= 1e-6,
                    'All-unknown core gained an observation gradient')
    previous = C.starting_state()['model']
    checkpoints = []
    for index in range(3):
        step = 6144+index*1536
        path = fit/('step_%06d.pt' % step)
        state = torch.load(path, map_location='cpu', weights_only=True)
        require(state['configuration'] == configuration and state['protocol'] == C.PROTOCOL
            and state['step'] == step and state['epoch'] == index+4
            and state['encoder_revision'] == C.Q.F.REVISION, 'Stored checkpoint identity differs')
        model = C.G.M28Core('C0', 2048, 128)
        model.load_state_dict(state['model'], strict=True)
        require(all(v.dtype == torch.float32 and bool(torch.isfinite(v).all()) for v in state['model'].values()),
                'Nonfinite or non-FP32 checkpoint weights')
        require(any(not torch.equal(v, previous[k]) for k, v in state['model'].items()), 'Pass did not change model weights')
        optimizer = state['optimizer']
        ids = [p for group in optimizer['param_groups'] for p in group['params']]
        parameters = list(model.parameters())
        require(len(ids) == len(parameters) == len(optimizer['state']), 'Incomplete AdamW state')
        for ident, parameter in zip(ids, parameters):
            saved = optimizer['state'][ident]
            require(int(saved['step']) == step, 'AdamW state step differs')
            for key in ('exp_avg', 'exp_avg_sq'):
                require(saved[key].shape == parameter.shape and bool(torch.isfinite(saved[key]).all()),
                        'Invalid AdamW moment tensor')
        part = rows[index*1536:(index+1)*1536]
        require(dict(Counter(x['species'] for x in part)) == summary['epochs'][index]['by_species'],
                'Summary species exposure differs')
        for key in ('raw_NLL', 'normalized_NLL', 'weighted_total'):
            require(math.isclose(sum(x[key] for x in part)/1536, summary['epochs'][index]['mean_losses'][key],
                    rel_tol=1e-12, abs_tol=1e-12), 'Summary loss differs from complete log')
        checkpoints.append(dict(step=step, epoch=index+4, file=str(path.relative_to(C.ROOT))))
        previous = state['model']
    A.validate_checkpoint(state, arm)
    C.write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], arm=arm,
        source_fit_job=fit_job, configuration=configuration, protocol=C.PROTOCOL,
        final_checkpoint=summary['final_checkpoint'], checkpoints=checkpoints,
        verified_additional_optimizer_steps=4608, verified_cumulative_optimizer_steps=9216,
        original_draw_order_observations_and_importance_match=True, logged_objective_arithmetic_match=True,
        final_checkpoint_strict_load=True, all_checkpoint_weights_and_optimizer_moments_finite=True,
        optimizer_steps_this_run=0, native_forward=False, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began, interpretation='training_record_integrity_not_chain_accuracy'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', choices=C.ARMS, required=True)
    parser.add_argument('--fit-job', required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.arm, args.fit_job)
