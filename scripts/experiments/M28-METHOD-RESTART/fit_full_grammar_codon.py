#!/usr/bin/env python3
"""Fixed fresh/joint TRAIN codon-factor continuation, then record verification."""
import argparse
from collections import Counter
import json
import math
import os
import time
import torch
import preflight_full_grammar_codon as N
from src.m28.full_grammar_codon import codon_coarse_channel_nll, centered_codon_scores

B, C, J, G, Q, ROOT = N.B, N.C, N.J, N.G, N.Q, N.ROOT
A, CONFIG, PROTOCOL = N.A, N.CONFIG, N.A.PROTOCOL
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-VERIFICATION'


def update(model, raw_scores, optimizer, backbone, loss, multiplier):
    total = loss * multiplier
    if not bool(torch.isfinite(total)) or float(loss.detach()) < -1e-8:
        raise ValueError('Nonfinite or negative observation NLL')
    total.backward()
    parameters = list(model.parameters()) + [raw_scores]
    if any(p.grad is None for p in parameters):
        raise ValueError('Missing head or codon gradient')
    norm = torch.nn.utils.clip_grad_norm_(parameters, 1., error_if_nonfinite=True)
    if any(p.requires_grad or p.grad is not None for p in backbone.parameters()):
        raise ValueError('Pretrained encoder must remain frozen')
    optimizer.step()
    if any(not bool(torch.isfinite(p).all()) for p in parameters):
        raise ValueError('Nonfinite updated head or codon coefficient')
    torch.cuda.synchronize()
    return float(total.detach()), float(norm)


def fit(preflight_job):
    began = time.perf_counter()
    ready = json.loads((N.OUT/'summary.json').read_text())
    if ((N.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != preflight_job or ready['configuration'] != CONFIG
            or len(ready['records']) != 2 or ready['optimizer_steps_this_run'] != 2
            or ready['checkpoints_written'] or ready['DEV_used'] or ready['sealed_data_used']
            or not all(r['zero_codon_native_NLL_parent_exact']
                and r['actual_native_logits_raw_adapter_exact'] and r['actual_codon_scores_raw_adapter_exact']
                and r['actual_updated_context_MAP_raw_adapter_exact'] for r in ready['records'])):
        raise ValueError('Completed matching disposable native codon preflight required')
    plan, windows, catalogs = B.plan_and_population()
    C.write(FIT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, preflight_job=preflight_job, starting_fit_job='13434603', support_job='13453408'))
    source = C.starting_state()
    model, optimizer = C.load_model_optimizer(source)
    raw_scores = torch.nn.Parameter(torch.zeros(64, dtype=torch.float32, device='cpu'))
    optimizer.add_param_group({'params': [raw_scores]})
    sequences, descriptors = C.P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    epochs = []
    with (FIT/'training.jsonl').open('x', buffering=512*1024) as handle:
        for pass_index in range(3):
            pass_started, sums, counts = time.perf_counter(), Counter(), Counter()
            for expected in plan[pass_index*1536:(pass_index+1)*1536]:
                started = time.perf_counter()
                w = windows[expected['window_id']]
                chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
                optimizer.zero_grad(set_to_none=True)
                joined = J.potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
                cuts = tuple(x for x in (offset, offset+w['valid_bases']) if 0 < x < len(context))
                result = codon_coarse_channel_nll(context, joined, 128, raw_scores, offset,
                    offset+w['valid_bases'], targets, .1, cuts=cuts, checkpoint_bp=256)
                raw = result['loss']
                total, norm = update(model, raw_scores, optimizer, backbone, raw, expected['objective_multiplier'])
                step = 4608 + expected['additional_step']
                row = dict(step=step, epoch=pass_index+4, codon_optimizer_step=expected['additional_step'],
                    **{k:expected[k] for k in ('additional_step', 'pass_index', 'draw_ordinal', 'window_id', 'species',
                        'importance_weight', 'known_observation_denominator', 'fine_positions', 'coarse_positions',
                        'fully_unknown_positions')}, raw_NLL=float(raw.detach()),
                    normalized_NLL=float(raw.detach())/expected['known_observation_denominator'],
                    weighted_total=total, gradient_norm_before_clip=norm, seconds=time.perf_counter()-started)
                handle.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update({k:row[k] for k in ('raw_NLL', 'normalized_NLL', 'weighted_total')})
                counts[w['species']] += 1
                if step % 64 == 0:
                    handle.flush()
                    print(json.dumps(dict(step=step, final_step=9216, recent=row,
                        elapsed_seconds=time.perf_counter()-began)), flush=True)
                del result, raw, joined, targets
            checkpoint = FIT/('step_%06d.pt' % step)
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, epoch=pass_index+4, step=step,
                encoder_revision=Q.F.REVISION, model=model.state_dict(), optimizer=optimizer.state_dict(),
                raw_codon_scores=raw_scores.detach().clone()), checkpoint)
            epochs.append(dict(epoch=pass_index+4, steps=1536, by_species=dict(counts),
                mean_losses={k:v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if (step != 9216 or any(int(optimizer.state[p]['step']) != 9216 for p in model.parameters())
            or int(optimizer.state[raw_scores]['step']) != 4608):
        raise ValueError('Original head9216 and new codon4608 optimizer ages required')
    if sum(p.stat().st_size for p in FIT.iterdir() if p.is_file()) > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared fit artifact ceiling exceeded')
    C.write(FIT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG, protocol=PROTOCOL, preflight_job=preflight_job, source_fit_job='13434603',
        support_job='13453408', cumulative_head_optimizer_steps=9216, optimizer_steps_this_run=4608,
        codon_optimizer_steps=4608, head_parameters=366191, raw_codon_parameters=64, trained_parameters=366255,
        fresh_draws=4608, fresh_unique_windows=4219, draws_per_pass=1536, epochs=epochs,
        final_checkpoint=str(checkpoint.relative_to(ROOT)), original_TRAIN_DNA_sources=descriptors,
        centered_codon_scores=centered_codon_scores(raw_scores.detach()).tolist(),
        fixed_fresh_draw_order_and_unit_importance_preserved=True, all_fixed_draws_retained=True,
        original_targets_changed=False, context_or_structural_grammar_changed=False,
        encoder_condition='JOINT24576', native_CDS_priors_used=False, pretrained_encoder_frozen_no_grad=True,
        saved_source_weights_and_optimizer_unchanged=True, cached_features_used=False, DEV_used=False,
        sealed_data_used=False, seconds=time.perf_counter()-began,
        interpretation='complete_codon_energy_TRAIN_continuation_not_generalization'))


def verify(fit_job):
    began = time.perf_counter()
    summary = json.loads((FIT/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed']
            or summary['job_id'] != fit_job or summary['configuration'] != CONFIG
            or summary['optimizer_steps_this_run'] != 4608 or summary['cumulative_head_optimizer_steps'] != 9216
            or summary['codon_optimizer_steps'] != 4608 or summary['DEV_used'] or summary['sealed_data_used']):
        raise ValueError('Completed matching final codon fit required')
    plan, _, _ = B.plan_and_population()
    rows = list(B.P.rows(FIT/'training.jsonl'))
    if len(rows) != 4608:
        raise ValueError('Complete fresh exposure log required')
    for i, (row, expected) in enumerate(zip(rows, plan)):
        if (row['step'] != 4609+i or row['codon_optimizer_step'] != i+1 or row['epoch'] != 4+i//1536
                or any(row[k] != expected[k] for k in ('additional_step', 'pass_index', 'draw_ordinal',
                    'window_id', 'species', 'importance_weight', 'known_observation_denominator',
                    'fine_positions', 'coarse_positions', 'fully_unknown_positions'))):
            raise ValueError('Actual draw, head step or codon step differs')
        raw = row['raw_NLL']
        if (not all(math.isfinite(row[k]) and row[k] >= -1e-8 for k in
                    ('raw_NLL', 'normalized_NLL', 'weighted_total', 'gradient_norm_before_clip', 'seconds'))
                or not math.isclose(row['normalized_NLL'], raw/expected['known_observation_denominator'], rel_tol=1e-12, abs_tol=1e-12)
                or not math.isclose(row['weighted_total'], raw*expected['objective_multiplier'], rel_tol=1e-12, abs_tol=1e-12)):
            raise ValueError('Logged original structured objective arithmetic differs')
        if expected['fine_positions'] + expected['coarse_positions'] == 0:
            if abs(raw) > 1e-8 or row['gradient_norm_before_clip'] > 1e-6:
                raise ValueError('All-unknown core gained observation gradient')
    previous = C.starting_state()['model']
    for index in range(3):
        step, codon_step = 6144+index*1536, (index+1)*1536
        checkpoint = FIT/('step_%06d.pt' % step)
        state = torch.load(checkpoint, map_location='cpu', weights_only=True)
        if (state['configuration'] != CONFIG or state['protocol'] != PROTOCOL or state['step'] != step
                or state['epoch'] != index+4 or state['encoder_revision'] != Q.F.REVISION):
            raise ValueError('Checkpoint identity differs')
        model = G.M28Core('C0', 2048, 128)
        model.load_state_dict(state['model'], strict=True)
        raw_scores = state['raw_codon_scores']
        centered_codon_scores(raw_scores)
        if (raw_scores.dtype != torch.float32
                or not all(v.dtype == torch.float32 and bool(torch.isfinite(v).all()) for v in state['model'].values())
                or all(torch.equal(v, previous[k]) for k, v in state['model'].items())):
            raise ValueError('Nonfinite, wrong-dtype or unchanged head checkpoint')
        parameters = list(model.parameters()) + [raw_scores]
        saved_optimizer = state['optimizer']
        groups = saved_optimizer['param_groups']
        if len(groups) != 2 or len(groups[1]['params']) != 1:
            raise ValueError('Original head and one new codon parameter group required')
        ids = [ident for group in groups for ident in group['params']]
        if len(ids) != len(parameters) or len(saved_optimizer['state']) != len(parameters):
            raise ValueError('Incomplete AdamW state')
        for group in groups:
            if (group['lr'] != CONFIG['lr'] or group['weight_decay'] != CONFIG['weight_decay']
                    or list(group['betas']) != CONFIG['betas'] or group['eps'] != CONFIG['eps']):
                raise ValueError('Saved head or codon optimizer hyperparameters differ')
        for ordinal, (ident, parameter) in enumerate(zip(ids, parameters)):
            saved = saved_optimizer['state'][ident]
            age = codon_step if ordinal == len(parameters)-1 else step
            if (int(saved['step']) != age or any(saved[k].shape != parameter.shape
                    or not bool(torch.isfinite(saved[k]).all()) for k in ('exp_avg', 'exp_avg_sq'))):
                raise ValueError('Nonfinite, misregistered or wrong-age AdamW moments')
        part = rows[index*1536:(index+1)*1536]
        if (dict(Counter(r['species'] for r in part)) != summary['epochs'][index]['by_species']
                or any(not math.isclose(sum(r[k] for r in part)/1536, summary['epochs'][index]['mean_losses'][k],
                    rel_tol=1e-12, abs_tol=1e-12) for k in ('raw_NLL', 'normalized_NLL', 'weighted_total'))):
            raise ValueError('Complete pass summary differs from log')
        previous = state['model']
    A.validate_checkpoint(state)
    if summary['centered_codon_scores'] != centered_codon_scores(state['raw_codon_scores']).tolist():
        raise ValueError('Saved final codon score readout differs')
    C.write(VERIFY/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], source_fit_job=fit_job,
        configuration=CONFIG, protocol=PROTOCOL, final_checkpoint=summary['final_checkpoint'],
        verified_additional_optimizer_steps=4608, verified_final_head_steps=9216, verified_codon_steps=4608,
        fixed_draws_original_observations_and_objective_arithmetic_match=True,
        all_checkpoint_parameters_and_moments_finite_and_correct_age=True,
        final_checkpoint_strict_load=True, native_forward=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        interpretation='record_integrity_not_chain_accuracy'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('fit', 'verify'))
    parser.add_argument('--source-job', required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'fit':
        fit(args.source_job)
    else:
        verify(args.source_job)
