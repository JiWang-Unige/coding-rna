#!/usr/bin/env python3
"""Matched original-TRAIN structural continuation, adding only28 native priors."""
import argparse
from collections import Counter
import copy
import io
import json
import math
import os
import time

import torch
import fit_full_grammar_objective_contrast as C
import annotate_full_grammar_native_cds_prior as A
import preflight_full_grammar_native_cds_prior as V
from native_cds_role_prior import NativeCDSRolePrior, native_cds_log_odds

ROOT, F, G, Q = C.ROOT, C.F, C.G, C.Q
PREFLIGHT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-UPDATE-PREFLIGHT'
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-FIT'
PROTOCOL = A.PROTOCOL
CONFIG = dict(F.CONFIG, **{k:v for k, v in A.IDENTITY.items() if k not in F.CONFIG})
CONFIG.update(A.IDENTITY)
CONFIG.update(new_prior_optimizer_steps=4608, old_C0_optimizer_steps=9216,
    original_C0_core_or_DNA_branch_changed=False, maximum_new_fit_bytes=96*1024*1024,
    matched_control='OBJECTIVE-CONTRAST/STRUCTURED_CONTINUE_final9216',
    gradient_clip_includes_new28=True)


def load_model_optimizer(state):
    torch.manual_seed(0)
    model = NativeCDSRolePrior().to('cuda').train()
    model.load_c0_weights(state['model'])
    old = [p for name, p in model.named_parameters() if name != 'native_cds_role_prior.weight']
    if (sum(p.numel() for p in old) != 366191 or sum(p.numel() for p in model.parameters()) != 366219
            or [name for name, _ in model.named_parameters() if name != 'native_cds_role_prior.weight']
               != list(state['model'])):
        raise ValueError('Original named C0 parameter order and exactly28 new coefficients required')
    optimizer = torch.optim.AdamW(old, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
        betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    optimizer.load_state_dict(copy.deepcopy(state['optimizer']))
    if (len(optimizer.state) != len(old) or any(int(optimizer.state[p]['step']) != 4608 for p in old)
            or any(optimizer.state[p][k].shape != p.shape for p in old for k in ('exp_avg', 'exp_avg_sq'))):
        raise ValueError('Every old parameter must restore its own final4608 AdamW moments')
    optimizer.add_param_group(dict(params=[model.native_cds_role_prior.weight]))
    for group in optimizer.param_groups:
        if any(group[k] != CONFIG[k] for k in ('lr', 'weight_decay', 'eps')) or list(group['betas']) != CONFIG['betas']:
            raise ValueError('Both groups must retain the original AdamW hyperparameters')
    if model.native_cds_role_prior.weight in optimizer.state:
        raise ValueError('New prior has no inherited moments or optimizer step')
    return model, optimizer


def require_gradient_preflight(job):
    value = json.loads((V.OUT/'summary.json').read_text())
    if ((V.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
            or value['job_id'] != job or value['configuration'] != V.CONFIG or value['actual_witnesses'] != 8
            or not value['all_original_eval_predictions_replayed_exact']
            or not value['both_native_core_margins_replayed_exact']
            or not value['all_real_prior_VJPs_match_analytic'] or value['optimizer_steps_this_run'] != 0):
        raise ValueError('Complete exact-replay and actual-gradient initialization evidence required')


def plan_inputs():
    draws, _, windows, catalogs, plan, _ = F.original_plan()
    if plan != [json.loads(x) for x in (F.READY/'training_plan.jsonl').read_text().splitlines()]:
        raise ValueError('Original draw order, labels or likelihood multipliers changed')
    sequences, descriptors = C.P.read_train_DNA()
    return draws, windows, catalogs, plan, sequences, descriptors


def preflight(gradient_job):
    began = time.perf_counter(); require_gradient_preflight(gradient_job)
    state = C.starting_state()
    draws, windows, catalogs, plan, sequences, descriptors = plan_inputs()
    selected = [next(i for i, d in enumerate(draws) if d['species'] == s) for s in Q.ALLOW]
    C.write(PREFLIGHT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, gradient_preflight_job=gradient_job, selected_draw_ordinals=selected))
    full, tokenizer = A.load_native(); records = []
    for ordinal in selected:
        began_row = time.perf_counter(); w = windows[draws[ordinal]['window_id']]
        chromosome, geometry, context, offset, targets = C.example(w, plan[ordinal], sequences, catalogs)
        core, optimizer = load_model_optimizer(state)
        adapter = A.NativeCDSHeadAdapter(core, full.score)
        joined = C.potentials(adapter, full.model, tokenizer, chromosome, geometry, w['strand'])
        raw, components = C.objective('STRUCTURED_CONTINUE', context, joined, offset, w['valid_bases'], targets)
        total, norm = C.update(core, optimizer, full, raw, plan[ordinal]['objective_multiplier'])
        prior = core.native_cds_role_prior.weight
        old = [p for name, p in core.named_parameters() if name != 'native_cds_role_prior.weight']
        if (any(int(optimizer.state[p]['step']) != 4609 for p in old)
                or int(optimizer.state[prior]['step']) != 1 or not bool(prior.detach().count_nonzero())):
            raise ValueError('Actual disposable update must advance old4608->4609 and new0->1')
        # Round-trip only the updated neural weights in RAM. No scientific
        # checkpoint file is created and the strict raw caller still rejects4609.
        buffer = io.BytesIO(); torch.save(core.state_dict(), buffer); buffer.seek(0)
        restored = NativeCDSRolePrior().to('cuda').eval()
        restored.load_state_dict(torch.load(buffer, map_location='cuda', weights_only=True), strict=True)
        core.eval()
        left, right = geometry['genomic_windows'][0]
        dna_string = C.re.sub('[^ACGT]', 'N', chromosome[left:right])
        if w['strand'] == '-': dna_string = Q.M.rc(dna_string)
        features = G.encode_window(full.model, tokenizer, dna_string, 'cuda')
        dna = G.one_hot(dna_string, Q.M.WIDTH)[None].to('cuda')
        with torch.inference_mode():
            direct = core(features, dna, [right-left], native_cds_log_odds(full.score, features))['segmentation_logits']
            replay = A.NativeCDSHeadAdapter(restored, full.score)(features, dna, [right-left])['segmentation_logits']
        if not torch.equal(direct, replay):
            raise ValueError('Strictly restored updated weights must replay the real three-input raw adapter')
        rejection = dict(configuration=CONFIG, protocol=PROTOCOL, step=4609, epoch=4, encoder_revision=Q.F.REVISION)
        try: A.validate_checkpoint(rejection)
        except ValueError: pass
        else: raise ValueError('Formal raw caller accepted a disposable nonfinal checkpoint')
        row = dict(draw_ordinal=ordinal, window_id=w['id'], species=w['species'],
            raw_NLL=float(raw.detach()), weighted_total=total, gradient_norm_before_clip=norm,
            old_optimizer_step=4609, new_prior_optimizer_step=1,
            prior_weights=prior.detach().cpu().tolist(),
            strict_neural_weight_in_memory_roundtrip_and_raw_adapter_exact=True,
            nonfinal_checkpoint_rejected=True, seconds=time.perf_counter()-began_row)
        records.append(row); print(json.dumps(row, allow_nan=False), flush=True)
        del core, optimizer, adapter, restored, raw, components, joined, direct, replay, buffer
    C.write(PREFLIGHT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, gradient_preflight_job=gradient_job, records=records,
        original_TRAIN_DNA_sources=descriptors, optimizer_steps_this_run=2, disposable_updates_only=True,
        checkpoint_created=False, formal_fit_started=False, saved_source_weights_and_optimizer_unchanged=True,
        pretrained_encoder_and_heads_frozen_no_grad=True, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began, interpretation='matched_optimizer_and_raw_interface_not_annotation_accuracy'))


def fit(preflight_job):
    began = time.perf_counter(); value = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
            or value['job_id'] != preflight_job or value['configuration'] != CONFIG
            or len(value['records']) != 2 or value['checkpoint_created'] or not value['disposable_updates_only']):
        raise ValueError('Complete matching two-update raw-interface and optimizer preflight required')
    C.write(FIT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, preflight_job=preflight_job, starting_fit_job='13434603'))
    state = C.starting_state(); core, optimizer = load_model_optimizer(state)
    draws, windows, catalogs, plan, sequences, descriptors = plan_inputs()
    full, tokenizer = A.load_native(); adapter = A.NativeCDSHeadAdapter(core, full.score)
    step = 4608; epochs = []; torch.cuda.reset_peak_memory_stats()
    with (FIT/'training.jsonl').open('x') as log:
        for pass_index in range(3):
            pass_started = time.perf_counter(); sums = Counter(); species_counts = Counter()
            for d, expected in zip(draws, plan):
                began_row = time.perf_counter(); w = windows[d['window_id']]
                chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
                optimizer.zero_grad(set_to_none=True)
                joined = C.potentials(adapter, full.model, tokenizer, chromosome, geometry, w['strand'])
                raw, components = C.objective('STRUCTURED_CONTINUE', context, joined, offset, w['valid_bases'], targets)
                total, norm = C.update(core, optimizer, full, raw, expected['objective_multiplier']); step += 1
                row = dict(step=step, additional_step=step-4608, epoch=pass_index+4,
                    **{k:expected[k] for k in ('draw_ordinal', 'window_id', 'species', 'importance_weight',
                        'known_observation_denominator', 'fine_positions', 'coarse_positions', 'fully_unknown_positions')},
                    raw_NLL=float(raw.detach()), normalized_NLL=float(raw.detach())/expected['known_observation_denominator'],
                    weighted_total=total, gradient_norm_before_clip=norm,
                    prior_weight_l2=float(core.native_cds_role_prior.weight.detach().norm()),
                    seconds=time.perf_counter()-began_row)
                log.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update({k:row[k] for k in ('raw_NLL', 'normalized_NLL', 'weighted_total')}); species_counts[w['species']] += 1
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, final_step=9216, recent=row,
                        elapsed_seconds=time.perf_counter()-began)), flush=True)
                del raw, components, joined, targets
            checkpoint = FIT/('step_%06d.pt' % step)
            if checkpoint.exists(): raise FileExistsError(checkpoint)
            old = [p for name, p in core.named_parameters() if name != 'native_cds_role_prior.weight']
            if (any(int(optimizer.state[p]['step']) != step for p in old)
                    or int(optimizer.state[core.native_cds_role_prior.weight]['step']) != step-4608):
                raise ValueError('Every saved old/new AdamW exposure must match its own starting step')
            torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, epoch=pass_index+4, step=step,
                encoder_revision=Q.F.REVISION, model=core.state_dict(), optimizer=optimizer.state_dict()), checkpoint)
            epochs.append(dict(epoch=pass_index+4, steps=1536, by_species=dict(species_counts),
                mean_losses={k:v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 9216: raise ValueError('All4608 new optimizer updates required')
    if sum(p.stat().st_size for p in FIT.iterdir() if p.is_file()) > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared new fit artifact ceiling exceeded')
    C.write(FIT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, preflight_job=preflight_job, starting_fit_job='13434603',
        cumulative_optimizer_steps=step, optimizer_steps_this_run=4608, new_prior_optimizer_steps=4608,
        old_C0_optimizer_steps=9216, epochs=epochs, final_checkpoint=str(checkpoint.relative_to(ROOT)),
        trained_parameters=366219, original_draws_per_pass=len(draws), original_unique_windows=len(windows),
        original_TRAIN_DNA_sources=descriptors, original_draw_order_and_importance_preserved=True,
        original_targets_changed=False, context_or_inference_grammar_changed=False,
        pretrained_encoder_and_heads_frozen_no_grad=True, saved_source_weights_and_optimizer_unchanged=True,
        cached_features_used=False, DEV_used=False, sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began,
        interpretation='native_prior_original_TRAIN_continuation_not_generalization_or_publication_success'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('preflight', 'fit')); parser.add_argument('--source-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'preflight': preflight(args.source_job)
    else: fit(args.source_job)
