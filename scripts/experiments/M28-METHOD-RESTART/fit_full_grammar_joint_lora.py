#!/usr/bin/env python3
"""One fixed same-start fresh/joint adapter continuation and CPU record verification."""
import argparse
from collections import Counter
import json
import math
import os
import time
import torch
import preflight_full_grammar_joint_lora as N
import annotate_full_grammar_joint_lora as A

B, C, G, Q, ROOT, L = N.B, N.C, N.G, N.Q, N.ROOT, N.L
PROTOCOL = A.PROTOCOL
CONFIG = dict(B.CONFIG, **A.IDENTITY)
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-LORA-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-LORA-VERIFICATION'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-LORA-TRAIN-CHECKS'


def verify_optimizer(saved, model, adapters, head_step, adapter_step):
    names = list(A.adapter_shapes())
    parameters = list(model.parameters()) + [adapters[k] for k in names]
    groups = saved['param_groups']
    if (len(groups) != 2 or len(groups[0]['params']) != len(list(model.parameters()))
            or groups[1].get('param_names') != names or len(groups[1]['params']) != 16):
        raise ValueError('Original head group and sixteen canonically named adapter parameters required')
    ids = [ident for group in groups for ident in group['params']]
    if (len(ids) != len(parameters) or len(set(ids)) != len(ids)
            or set(ids) != set(saved['state'])):
        raise ValueError('Complete one-to-one AdamW parameter registration required')
    for group in groups:
        if (group['lr'] != CONFIG['lr'] or group['weight_decay'] != CONFIG['weight_decay']
                or list(group['betas']) != CONFIG['betas'] or group['eps'] != CONFIG['eps']):
            raise ValueError('Head or adapter optimizer hyperparameters differ')
    for ordinal, (ident, parameter) in enumerate(zip(ids, parameters)):
        state = saved['state'][ident]
        age = head_step if ordinal < len(list(model.parameters())) else adapter_step
        if (int(state['step']) != age or any(state[k].shape != parameter.shape
                or state[k].dtype != torch.float32 or not bool(torch.isfinite(state[k]).all())
                for k in ('exp_avg', 'exp_avg_sq'))):
            raise ValueError('Misregistered, nonfinite or wrong-age AdamW moments')


def update(model, adapters, optimizer, backbone, raw, multiplier):
    total = raw * multiplier
    if not bool(torch.isfinite(total)) or float(raw.detach()) < -1e-8:
        raise ValueError('Finite nonnegative observation NLL required')
    total.backward()
    parameters = list(model.parameters()) + L.adapter_parameters(adapters)
    if any(p.grad is None for p in parameters):
        raise ValueError('Missing head or native adapter gradient')
    norm = torch.nn.utils.clip_grad_norm_(parameters, 1., error_if_nonfinite=True)
    adapter_ids = {id(p) for p in L.adapter_parameters(adapters)}
    if any(p.requires_grad or p.grad is not None for p in backbone.parameters() if id(p) not in adapter_ids):
        raise ValueError('All released native encoder base weights must remain frozen and gradient-free')
    optimizer.step()
    if any(not bool(torch.isfinite(p).all()) for p in parameters):
        raise ValueError('Nonfinite updated head or adapter parameter')
    torch.cuda.synchronize()
    return float(total.detach()), float(norm)


def fit(preflight_job, checks_job):
    began = time.perf_counter()
    ready = json.loads((N.OUT/'summary.json').read_text())
    if ((N.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != preflight_job or ready['configuration'] != N.CONFIG
            or len(ready['records']) != 2 or ready['optimizer_steps_this_run'] != 2
            or ready['checkpoints_written'] or ready['DEV_used'] or ready['sealed_data_used']
            or not all(r['zero_adapter_native_potentials_NLL_and_role_VJP_parent_exact']
                and r['actual_updated_native_logits_raw_exact'] and r['actual_updated_context_MAP_raw_exact']
                for r in ready['records'])):
        raise ValueError('Completed matching disposable native adapter preflight required')
    checks = CHECKS/checks_job
    checked = json.loads((checks/'summary.json').read_text())
    if ((checks/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['job_id'] != checks_job or checked['configuration'] != CONFIG
            or checked['tests_run'] != 11 or checked['errors'] or checked['failures']):
        raise ValueError('Completed strict adapter checkpoint, endpoint and arithmetic checks required')
    plan, windows, catalogs = B.plan_and_population()
    C.write(FIT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, preflight_job=preflight_job, checks_job=checks_job,
        starting_fit_job='13434603', support_job='13453408'))
    source = C.starting_state()
    model, optimizer = C.load_model_optimizer(source)
    sequences, descriptors = C.P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    adapters = L.attach_adapters(backbone)
    backbone.eval()
    parameters = L.adapter_parameters(adapters)
    optimizer.add_param_group({'params': parameters, 'param_names': list(A.adapter_shapes())})
    if sum(p.numel() for p in parameters) != L.PARAMETERS:
        raise ValueError('Exactly102400 fixed native adapter parameters required')
    torch.cuda.reset_peak_memory_stats()
    epochs = []
    with (FIT/'training.jsonl').open('x', buffering=512*1024) as handle:
        for pass_index in range(3):
            pass_started, sums, counts = time.perf_counter(), Counter(), Counter()
            for expected in plan[pass_index*1536:(pass_index+1)*1536]:
                started = time.perf_counter()
                w = windows[expected['window_id']]
                chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
                optimizer.zero_grad(set_to_none=True)
                joined = N.potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
                raw, _ = C.objective('STRUCTURED_CONTINUE', context, joined, offset, w['valid_bases'], targets)
                total, norm = update(model, adapters, optimizer, backbone, raw, expected['objective_multiplier'])
                step = 4608 + expected['additional_step']
                row = dict(step=step, epoch=pass_index+4, adapter_optimizer_step=expected['additional_step'],
                    **{k:expected[k] for k in ('additional_step', 'pass_index', 'draw_ordinal', 'window_id', 'species',
                        'importance_weight', 'known_observation_denominator', 'fine_positions', 'coarse_positions',
                        'fully_unknown_positions')}, raw_NLL=float(raw.detach()),
                    normalized_NLL=float(raw.detach())/expected['known_observation_denominator'],
                    weighted_total=total, gradient_norm_before_clip=norm, seconds=time.perf_counter()-started)
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update({k:row[k] for k in ('raw_NLL', 'normalized_NLL', 'weighted_total')})
                counts[w['species']] += 1
                if expected['additional_step'] <= 2 or step % 64 == 0:
                    handle.flush()
                    print(json.dumps(dict(step=step, final_step=9216, recent=row,
                        elapsed_seconds=time.perf_counter()-began,
                        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
                        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved())), flush=True)
                del raw, joined, targets
            checkpoint = FIT/('step_%06d.pt' % step)
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, epoch=pass_index+4, step=step,
                adapter_optimizer_step=(pass_index+1)*1536, encoder_revision=Q.F.REVISION,
                model=model.state_dict(), optimizer=optimizer.state_dict(), adapters=L.adapter_state(adapters)), checkpoint)
            epochs.append(dict(epoch=pass_index+4, steps=1536, by_species=dict(counts),
                mean_losses={k:v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if (step != 9216 or any(int(optimizer.state[p]['step']) != 9216 for p in model.parameters())
            or any(int(optimizer.state[p]['step']) != 4608 for p in parameters)):
        raise ValueError('Complete original head9216 and all native adapters4608 required')
    if sum(p.stat().st_size for p in FIT.iterdir() if p.is_file()) > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared fit artifact ceiling exceeded')
    C.write(FIT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG, protocol=PROTOCOL, preflight_job=preflight_job, checks_job=checks_job,
        source_fit_job='13434603', support_job='13453408', cumulative_head_optimizer_steps=9216,
        adapter_optimizer_steps=4608, optimizer_steps_this_run=4608, head_parameters=366191,
        native_adapter_parameters=102400, trained_parameters=468591, fresh_draws=4608,
        fresh_unique_windows=4219, draws_per_pass=1536, epochs=epochs,
        final_checkpoint=str(checkpoint.relative_to(ROOT)), original_TRAIN_DNA_sources=descriptors,
        fixed_fresh_draw_order_and_unit_importance_preserved=True, all_fixed_draws_retained=True,
        original_targets_changed=False, context_or_structural_grammar_changed=False,
        encoder_condition=A.CONDITION, pretrained_encoder_base_weights_frozen_no_grad=True,
        saved_source_weights_and_optimizer_unchanged=True, cached_features_used=False,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        interpretation='same_start_encoder_capacity_TRAIN_continuation_not_generalization'))


def verify(fit_job):
    began = time.perf_counter()
    summary = json.loads((FIT/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed']
            or summary['job_id'] != fit_job or summary['configuration'] != CONFIG
            or summary['optimizer_steps_this_run'] != 4608 or summary['cumulative_head_optimizer_steps'] != 9216
            or summary['adapter_optimizer_steps'] != 4608 or summary['DEV_used'] or summary['sealed_data_used']):
        raise ValueError('Completed matching final joint native-adapter fit required')
    plan, _, _ = B.plan_and_population()
    rows = list(B.P.rows(FIT/'training.jsonl'))
    if len(rows) != 4608:
        raise ValueError('Complete fixed fresh exposure log required')
    for i, (row, expected) in enumerate(zip(rows, plan)):
        if (row['step'] != 4609+i or row['adapter_optimizer_step'] != i+1 or row['epoch'] != 4+i//1536
                or any(row[k] != expected[k] for k in ('additional_step', 'pass_index', 'draw_ordinal',
                    'window_id', 'species', 'importance_weight', 'known_observation_denominator',
                    'fine_positions', 'coarse_positions', 'fully_unknown_positions'))):
            raise ValueError('Actual draw or head/adapter age differs')
        raw = row['raw_NLL']
        if (not all(math.isfinite(row[k]) and row[k] >= -1e-8 for k in
                    ('raw_NLL', 'normalized_NLL', 'weighted_total', 'gradient_norm_before_clip', 'seconds'))
                or not math.isclose(row['normalized_NLL'], raw/expected['known_observation_denominator'], rel_tol=1e-12, abs_tol=1e-12)
                or not math.isclose(row['weighted_total'], raw*expected['objective_multiplier'], rel_tol=1e-12, abs_tol=1e-12)):
            raise ValueError('Logged structured objective arithmetic differs')
        if expected['fine_positions'] + expected['coarse_positions'] == 0:
            if abs(raw) > 1e-8 or row['gradient_norm_before_clip'] > 1e-6:
                raise ValueError('All-unknown core gained observation gradient')
    previous = C.starting_state()['model']
    previous_adapters = None
    for index in range(3):
        step, adapter_step = 6144+index*1536, (index+1)*1536
        state = torch.load(FIT/('step_%06d.pt' % step), map_location='cpu', weights_only=True)
        if (state['configuration'] != CONFIG or state['protocol'] != PROTOCOL or state['step'] != step
                or state['adapter_optimizer_step'] != adapter_step or state['epoch'] != index+4
                or state['encoder_revision'] != Q.F.REVISION):
            raise ValueError('Checkpoint identity differs')
        model = G.M28Core('C0', 2048, 128)
        model.load_state_dict(state['model'], strict=True)
        A.validate_adapter_state(state['adapters'])
        if (not all(v.dtype == torch.float32 and bool(torch.isfinite(v).all()) for v in state['model'].values())
                or all(torch.equal(v, previous[k]) for k, v in state['model'].items())
                or all(torch.count_nonzero(v) == 0 for k, v in state['adapters'].items() if k.endswith('.B'))
                or previous_adapters is not None and all(torch.equal(v, previous_adapters[k]) for k, v in state['adapters'].items())):
            raise ValueError('Nonfinite, wrong-dtype or unchanged head/adapter checkpoint')
        verify_optimizer(state['optimizer'], model, state['adapters'], step, adapter_step)
        part = rows[index*1536:(index+1)*1536]
        if (dict(Counter(r['species'] for r in part)) != summary['epochs'][index]['by_species']
                or any(not math.isclose(sum(r[k] for r in part)/1536, summary['epochs'][index]['mean_losses'][k],
                    rel_tol=1e-12, abs_tol=1e-12) for k in ('raw_NLL', 'normalized_NLL', 'weighted_total'))):
            raise ValueError('Complete pass summary differs from log')
        previous, previous_adapters = state['model'], state['adapters']
    A.model_from_state(state, 'cpu')
    C.write(VERIFY/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], source_fit_job=fit_job,
        configuration=CONFIG, protocol=PROTOCOL, final_checkpoint=summary['final_checkpoint'],
        verified_additional_optimizer_steps=4608, verified_final_head_steps=9216, verified_adapter_steps=4608,
        fixed_draws_original_observations_and_objective_arithmetic_match=True,
        all_checkpoint_parameters_and_moments_finite_named_and_correct_age=True,
        final_checkpoint_strict_load=True, native_forward=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        interpretation='record_integrity_not_chain_accuracy'))


def checks():
    import unittest
    import test_m28_full_grammar_joint_lora_train as T
    import test_m28_joint_encoder_lora as V
    began = time.perf_counter()
    suite = unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromModule(module) for module in (T, V)])
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    C.write(CHECKS/os.environ['SLURM_JOB_ID']/'summary.json', dict(completed=result.wasSuccessful(),
        job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors), research_DNA_read=False,
        research_reference_read=False, checkpoint_read=False, native_model_read=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began))
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('fit', 'verify', 'checks'))
    parser.add_argument('--source-job')
    parser.add_argument('--checks-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'checks':
        checks()
    elif args.source_job is None or args.mode == 'fit' and args.checks_job is None:
        parser.error('--source-job and, for fit, --checks-job required')
    elif args.mode == 'fit':
        fit(args.source_job, args.checks_job)
    else:
        verify(args.source_job)
