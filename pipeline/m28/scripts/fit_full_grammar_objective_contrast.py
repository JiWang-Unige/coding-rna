#!/usr/bin/env python3
"""Three predeclared original-TRAIN continuations from final4608 plus AdamW."""
import argparse
from collections import Counter
import copy
import json
import math
import os
import re
import time
import torch
import fit_full_grammar_coarse_channel as F
import annotate_full_grammar_coarse_channel as R
from src.m28.full_grammar_coarse_channel import coarse_channel_nll
from src.m28.local_coarse_channel import local_coarse_channel_nll

ROOT, G, Q, P = F.ROOT, F.G, F.Q, F.P
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST'
PREFLIGHT = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-PREFLIGHT'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-CHECKS'
CHECKPOINT = F.FIT/'step_004608.pt'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_objective_contrast_contract.md'
ARMS = ('STRUCTURED_CONTINUE', 'LOCAL_CHANNEL', 'HALF_STRUCTURED_HALF_LOCAL')
CONFIG = dict(F.CONFIG, model='final4608_M28Core_C0_2048_128_366191_parameters_FP32_continuation',
    starting_checkpoint=str(CHECKPOINT.relative_to(ROOT)), starting_fit_job='13434603',
    starting_optimizer_steps=4608, optimizer_state='restore_original_final4608_AdamW_moments_and_step',
    passes=3, optimizer_steps=9216, additional_optimizer_steps=4608,
    formal_checkpoint='last_cumulative_step_009216_only_no_DEV_selection',
    loss='arm_specific_same_fine_coarse_channel_same_original_importance_and_known_denominator',
    arms=list(ARMS), mixed_structured_coefficient=0.5, mixed_local_coefficient=0.5,
    auxiliary_losses=False, maximum_new_fit_bytes=96*1024*1024)


def arm_config(arm):
    if arm not in ARMS:
        raise ValueError('One of the three predeclared objectives required')
    return dict(CONFIG, arm=arm)


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def starting_state():
    summary = json.loads((F.FIT/'summary.json').read_text())
    if ((F.FIT/'STATUS').read_text().strip() != 'COMPLETED' or summary['job_id'] != '13434603'
            or summary['configuration'] != F.CONFIG or not summary['completed']
            or summary['optimizer_steps'] != 4608 or summary['DEV_used']):
        raise ValueError('Completed unchanged original final4608 fit required')
    state = torch.load(CHECKPOINT, map_location='cpu', weights_only=True)
    R.validate_checkpoint(state)
    if state['configuration'] != F.CONFIG:
        raise ValueError('Original full fit configuration differs')
    return state


def load_model_optimizer(state):
    torch.manual_seed(0)
    model = G.M28Core('C0', 2048, 128).to('cuda').train()
    model.load_state_dict(state['model'], strict=True)
    parameters = list(model.parameters())
    if sum(p.numel() for p in parameters) != 366191:
        raise ValueError('C0 capacity differs')
    optimizer = torch.optim.AdamW(parameters, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
        betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    # load_state_dict may share CPU moment tensors; each arm gets its own copy.
    optimizer.load_state_dict(copy.deepcopy(state['optimizer']))
    for group in optimizer.param_groups:
        if (group['lr'] != CONFIG['lr'] or group['weight_decay'] != CONFIG['weight_decay']
                or list(group['betas']) != CONFIG['betas'] or group['eps'] != CONFIG['eps']):
            raise ValueError('Restored optimizer hyperparameters differ')
    if (len(optimizer.state) != len(parameters)
            or any(int(s['step']) != 4608 for s in optimizer.state.values())):
        raise ValueError('Every original AdamW parameter must start at step4608')
    return model, optimizer


def example(w, expected, sequences, catalogs):
    chromosome = sequences[w['species']][w['seqid']]
    geometry = G.context_geometry(w['start'], w['end'], len(chromosome))
    a, b = geometry['genomic_context']
    raw = chromosome[w['start']:w['end']]
    core = raw if w['strand'] == '+' else Q.M.rc(raw)
    targets = P.P.build_targets(w, catalogs[w['species']], core)
    n = w['valid_bases']
    fine, region = (torch.as_tensor(targets[k])[:n] for k in ('seg', 'region'))
    counts = (int((fine >= 0).sum()), int(((fine < 0) & (region >= 0)).sum()),
              int(((fine < 0) & (region < 0)).sum()))
    offset = w['start']-a if w['strand'] == '+' else b-w['end']
    if (counts != tuple(expected[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions'))
            or [a, b] != expected['genomic_context'] or [offset, offset+n] != expected['core_oriented_interval']):
        raise ValueError('Original observations or fixed context changed')
    context = re.sub('[^ACGT]', 'N', chromosome[a:b])
    if w['strand'] == '-':
        context = Q.M.rc(context)
    return chromosome, geometry, context, offset, targets


def potentials(model, backbone, tokenizer, chromosome, geometry, strand):
    values = []
    for left, right in geometry['genomic_windows']:
        sequence = re.sub('[^ACGT]', 'N', chromosome[left:right])
        if strand == '-':
            sequence = Q.M.rc(sequence)
        features = G.encode_window(backbone, tokenizer, sequence, 'cuda').clone()
        dna = G.one_hot(sequence, Q.M.WIDTH)[None].to('cuda')
        values.append(model(features, dna, [right-left])['segmentation_logits'][0, :right-left])
    return G.owned_context_logits(values, geometry, strand).cpu()


def objective(arm, context, logits, offset, n, targets):
    components = {}
    if arm != 'LOCAL_CHANNEL':
        cuts = tuple(x for x in (offset, offset+n) if 0 < x < len(context))
        components['structured_NLL'] = coarse_channel_nll(context, logits, 128, offset, offset+n,
            targets, .1, cuts=cuts, checkpoint_bp=256)['loss']
    if arm != 'STRUCTURED_CONTINUE':
        components['local_NLL'] = local_coarse_channel_nll(logits[offset:offset+n], targets, .1)['loss']
    raw = sum(components.values()) * (0.5 if arm == 'HALF_STRUCTURED_HALF_LOCAL' else 1.)
    return raw, components


def update(model, optimizer, backbone, raw, multiplier):
    total = raw * multiplier
    if not bool(torch.isfinite(total)) or float(raw.detach()) < -1e-8:
        raise ValueError('Nonfinite or negative observation NLL')
    total.backward()
    parameters = list(model.parameters())
    if any(p.grad is None for p in parameters):
        raise ValueError('Missing C0 gradient')
    norm = torch.nn.utils.clip_grad_norm_(parameters, 1., error_if_nonfinite=True)
    if any(p.requires_grad or p.grad is not None for p in backbone.parameters()):
        raise ValueError('Native pretrained encoder must remain frozen')
    optimizer.step()
    if any(not bool(torch.isfinite(p).all()) for p in parameters):
        raise ValueError('Nonfinite updated C0 parameter')
    torch.cuda.synchronize()
    return float(total.detach()), float(norm)


def preflight(checks_job):
    began = time.perf_counter()
    checks = json.loads((CHECKS/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not checks['completed']
            or checks['job_id'] != checks_job or checks['configuration'] != CONFIG):
        raise ValueError('Successful matching likelihood and gradient arithmetic required')
    state = starting_state()
    draws, _, windows, catalogs, plan, _ = F.original_plan()
    if plan != [json.loads(x) for x in (F.READY/'training_plan.jsonl').read_text().splitlines()]:
        raise ValueError('Original preparation plan differs')
    # Fixed before reading losses: first original draw of each species.
    selected = [next(i for i, d in enumerate(draws) if d['species'] == species) for species in Q.ALLOW]
    write(PREFLIGHT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        checks_job=checks_job, selected_draw_ordinals=selected, protocol=PROTOCOL))
    sequences, descriptors = P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    records = []
    for ordinal in selected:
        expected = plan[ordinal]
        w = windows[draws[ordinal]['window_id']]
        chromosome, geometry, context, offset, targets = example(w, expected, sequences, catalogs)
        for arm in ARMS:
            started = time.perf_counter()
            model, optimizer = load_model_optimizer(state)
            joined = potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
            if joined.shape != (len(context), 15) or not bool(torch.isfinite(joined).all()):
                raise ValueError('Finite native fifteen-role potentials required')
            raw, components = objective(arm, context, joined, offset, w['valid_bases'], targets)
            values = {k:float(v.detach()) for k, v in components.items()}
            total, norm = update(model, optimizer, backbone, raw, expected['objective_multiplier'])
            if (any(int(s['step']) != 4609 for s in optimizer.state.values())
                    or all(torch.equal(p.cpu(), state['model'][k]) for k, p in model.state_dict().items())):
                raise ValueError('Restored optimizer did not perform an independent step4609')
            row = dict(arm=arm, draw_ordinal=ordinal, window_id=w['id'], species=w['species'],
                raw_NLL=float(raw.detach()), components=values, weighted_total=total,
                gradient_norm_before_clip=norm, seconds=time.perf_counter()-started,
                initial_optimizer_step=4608, final_disposable_optimizer_step=4609)
            records.append(row)
            print(json.dumps(row, allow_nan=False), flush=True)
            del model, optimizer, raw, components, joined
    write(PREFLIGHT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG, protocol=PROTOCOL, checks_job=checks_job, records=records,
        original_TRAIN_DNA_sources=descriptors, formal_training_started=False,
        optimizer_steps_this_run=len(records), disposable_updates_only=True, checkpoints_written=False,
        saved_source_weights_and_optimizer_unchanged=True, pretrained_encoder_frozen_no_grad=True,
        cached_features_used=False, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began, interpretation='native_training_wiring_not_performance_or_arm_selection'))


def fit(arm, preflight_job):
    began = time.perf_counter()
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != preflight_job or ready['configuration'] != CONFIG
            or len(ready['records']) != 6 or ready['DEV_used'] or ready['checkpoints_written']):
        raise ValueError('Complete matching six-update native preflight required')
    configuration = arm_config(arm)
    out = OUT/arm
    write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=configuration,
        protocol=PROTOCOL, preflight_job=preflight_job, starting_fit_job='13434603'))
    state = starting_state()
    model, optimizer = load_model_optimizer(state)
    draws, _, windows, catalogs, plan, _ = F.original_plan()
    if plan != [json.loads(x) for x in (F.READY/'training_plan.jsonl').read_text().splitlines()]:
        raise ValueError('Original exposure plan changed')
    sequences, descriptors = P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    step = 4608
    epochs = []
    with (out/'training.jsonl').open('x') as handle:
        for pass_index in range(3):
            pass_started = time.perf_counter()
            sums = Counter()
            species_counts = Counter()
            for d, expected in zip(draws, plan):
                started = time.perf_counter()
                w = windows[d['window_id']]
                chromosome, geometry, context, offset, targets = example(w, expected, sequences, catalogs)
                optimizer.zero_grad(set_to_none=True)
                joined = potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
                raw, components = objective(arm, context, joined, offset, w['valid_bases'], targets)
                component_values = {k:float(v.detach()) for k, v in components.items()}
                total, norm = update(model, optimizer, backbone, raw, expected['objective_multiplier'])
                step += 1
                row = dict(step=step, additional_step=step-4608, epoch=pass_index+4,
                    **{k:expected[k] for k in ('draw_ordinal', 'window_id', 'species', 'importance_weight',
                        'known_observation_denominator', 'fine_positions', 'coarse_positions', 'fully_unknown_positions')},
                    raw_NLL=float(raw.detach()), normalized_NLL=float(raw.detach())/expected['known_observation_denominator'],
                    weighted_total=total, gradient_norm_before_clip=norm, components=component_values,
                    seconds=time.perf_counter()-started)
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update({k:row[k] for k in ('raw_NLL', 'normalized_NLL', 'weighted_total')})
                species_counts[w['species']] += 1
                if step % 64 == 0:
                    handle.flush()
                    print(json.dumps(dict(arm=arm, step=step, final_step=9216, recent=row,
                        elapsed_seconds=time.perf_counter()-began)), flush=True)
                del raw, components, joined, targets
            checkpoint = out/('step_%06d.pt' % step)
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(configuration=configuration, protocol=PROTOCOL, epoch=pass_index+4, step=step,
                encoder_revision=Q.F.REVISION, model=model.state_dict(), optimizer=optimizer.state_dict()), checkpoint)
            epochs.append(dict(epoch=pass_index+4, steps=1536, by_species=dict(species_counts),
                mean_losses={k:v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 9216 or any(int(s['step']) != 9216 for s in optimizer.state.values()):
        raise ValueError('Complete cumulative9216 continuation required')
    if all(torch.equal(p.cpu(), state['model'][k]) for k, p in model.state_dict().items()):
        raise ValueError('Continuation did not change original model')
    size = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    if size > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared new per-arm artifact ceiling exceeded')
    write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        configuration=configuration, protocol=PROTOCOL, preflight_job=preflight_job,
        starting_fit_job='13434603', cumulative_optimizer_steps=step, optimizer_steps_this_run=4608,
        epochs=epochs, final_checkpoint=str(checkpoint.relative_to(ROOT)), trained_parameters=366191,
        original_draws_per_pass=len(draws), original_unique_windows=len(windows),
        original_TRAIN_DNA_sources=descriptors, original_draw_order_and_importance_preserved=True,
        original_targets_changed=False, context_or_inference_grammar_changed=False,
        pretrained_encoder_frozen_no_grad=True, saved_source_weights_and_optimizer_unchanged=True,
        cached_features_used=False, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began,
        interpretation='original_TRAIN_objective_contrast_not_generalization_competitiveness_or_publication_success'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('preflight', 'fit'))
    parser.add_argument('--arm', choices=ARMS)
    parser.add_argument('--source-job', required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'preflight':
        preflight(args.source_job)
    else:
        if args.arm is None:
            parser.error('--arm required for fit')
        fit(args.arm, args.source_job)
