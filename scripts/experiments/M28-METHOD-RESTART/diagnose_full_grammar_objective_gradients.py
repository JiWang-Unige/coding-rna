#!/usr/bin/env python3
"""Fixed final4608 TRAIN parameter-gradient geometry, without optimizer updates."""
import json
import math
import os
import statistics
import time
import torch
import fit_full_grammar_objective_contrast as O

OUT = O.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-GRADIENTS'
CONFIG = dict(starting_checkpoint=str(O.CHECKPOINT.relative_to(O.ROOT)),
    starting_optimizer_steps=4608, selection='first_12_original_draws_per_species_in_original_order',
    draws_per_species=12, species=list(O.Q.ALLOW), seed=0,
    objective_configuration=O.CONFIG, parameter_gradients='same_native_potentials_for_both_components',
    normalization='original_importance_and_known_observation_denominator',
    mixed_coefficients=[0.5, 0.5], gradient_measurement='FP64_vector_norms_and_cosines_of_FP32_C0_gradients',
    optimizer_steps_this_run=0, cached_features_used=False, DEV_used=False, sealed_data_used=False)


def cosine(a, b):
    denominator = float(torch.linalg.vector_norm(a)*torch.linalg.vector_norm(b))
    return None if denominator == 0 else float(torch.dot(a, b))/denominator


def geometry(structured, local):
    mixed = (structured+local)*0.5
    sn, ln, mn = (float(torch.linalg.vector_norm(v)) for v in (structured, local, mixed))
    return dict(structured_norm=sn, local_norm=ln, mixed_norm=mn,
        local_over_structured_norm=None if sn == 0 else ln/sn,
        structured_local_cosine=cosine(structured, local),
        mixed_local_cosine=cosine(mixed, local), mixed_structured_cosine=cosine(mixed, structured))


def describe(records):
    result = dict(draws=len(records), local_norm_larger=sum(r['local_norm'] > r['structured_norm'] for r in records),
        zero_structured_gradients=sum(r['structured_norm'] == 0 for r in records),
        zero_local_gradients=sum(r['local_norm'] == 0 for r in records),
        opposing_component_gradients=sum(r['structured_local_cosine'] is not None and
            r['structured_local_cosine'] < 0 for r in records))
    for field in ('local_over_structured_norm', 'structured_local_cosine',
                  'mixed_local_cosine', 'mixed_structured_cosine'):
        values = [r[field] for r in records if r[field] is not None]
        result[field] = dict(defined=len(values), minimum=min(values) if values else None,
            median=statistics.median(values) if values else None, maximum=max(values) if values else None)
    return result


def run():
    began = time.perf_counter()
    torch.set_num_threads(2)
    state = O.starting_state()
    draws, _, windows, catalogs, plan, _ = O.F.original_plan()
    if plan != [json.loads(line) for line in (O.F.READY/'training_plan.jsonl').read_text().splitlines()]:
        raise ValueError('Original preparation plan differs')
    selected = sorted(i for species in O.Q.ALLOW for i in
        [j for j, draw in enumerate(draws) if draw['species'] == species][:12])
    if any(sum(draws[i]['species'] == s for i in selected) != 12 for s in O.Q.ALLOW):
        raise ValueError('Twelve original draws of each species required')
    O.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        selected_draw_ordinals=selected, protocol=O.PROTOCOL))
    torch.manual_seed(0)
    model = O.G.M28Core('C0', 2048, 128).to('cuda').train()
    model.load_state_dict(state['model'], strict=True)
    parameters = tuple(model.parameters())
    sequences, descriptors = O.P.read_train_DNA()
    backbone, tokenizer = O.G.load_native()
    preflight = json.loads((O.PREFLIGHT/'summary.json').read_text())
    if not preflight['completed'] or preflight['job_id'] != '13451086':
        raise ValueError('Completed original objective preflight required')
    prior = {(r['draw_ordinal'], r['arm']):r for r in preflight['records']}
    records = []
    matched_preflight_records = 0
    with (OUT/'gradients.jsonl').open('x') as handle:
        for ordinal in selected:
            started = time.perf_counter()
            expected = plan[ordinal]
            w = windows[draws[ordinal]['window_id']]
            chromosome, context_geometry, context, offset, targets = O.example(w, expected, sequences, catalogs)
            joined = O.potentials(model, backbone, tokenizer, chromosome, context_geometry, w['strand'])
            raw, components = O.objective('HALF_STRUCTURED_HALF_LOCAL', context, joined, offset,
                w['valid_bases'], targets)
            losses = [components[key]*expected['objective_multiplier'] for key in ('structured_NLL', 'local_NLL')]
            if any(not bool(torch.isfinite(loss)) or float(loss.detach()) < -1e-8 for loss in losses):
                raise ValueError('Finite nonnegative component losses required')
            gradients = []
            for index, loss in enumerate(losses):
                values = torch.autograd.grad(loss, parameters, retain_graph=index == 0)
                vector = torch.cat([value.detach().reshape(-1).cpu().double() for value in values])
                if not bool(torch.isfinite(vector).all()):
                    raise ValueError('Nonfinite parameter gradient')
                gradients.append(vector)
            row = dict(draw_ordinal=ordinal, window_id=w['id'], species=w['species'],
                fine_positions=expected['fine_positions'], coarse_positions=expected['coarse_positions'],
                fully_unknown_positions=expected['fully_unknown_positions'],
                importance_weight=expected['importance_weight'],
                known_observation_denominator=expected['known_observation_denominator'],
                objective_multiplier=expected['objective_multiplier'],
                raw_components={k:float(v.detach()) for k,v in components.items()},
                weighted_components=[float(loss.detach()) for loss in losses], **geometry(*gradients))
            # Existing independent native preflight anchors the first draw of each species.
            for arm, field in zip(O.ARMS, ('structured_norm', 'local_norm', 'mixed_norm')):
                if (ordinal, arm) in prior:
                    old = prior[ordinal, arm]
                    weighted = sum(row['weighted_components'])*0.5 if arm == O.ARMS[2] else \
                        row['weighted_components'][0 if arm == O.ARMS[0] else 1]
                    if (not math.isclose(weighted, old['weighted_total'], rel_tol=1e-6, abs_tol=1e-10) or
                            not math.isclose(row[field], old['gradient_norm_before_clip'], rel_tol=2e-4, abs_tol=1e-10)):
                        raise ValueError('First-draw loss or gradient differs from original native preflight')
                    matched_preflight_records += 1
            torch.cuda.synchronize()
            row['seconds'] = time.perf_counter()-started
            records.append(row)
            handle.write(json.dumps(row, allow_nan=False)+'\n')
            handle.flush()
            print(json.dumps(row, allow_nan=False), flush=True)
            del raw, joined, components, losses, gradients, values, vector, targets
    if (matched_preflight_records != 6 or any(p.grad is not None for p in parameters) or
            any(p.requires_grad or p.grad is not None for p in backbone.parameters()) or
            any(not torch.equal(value.cpu(), state['model'][name]) for name,value in model.state_dict().items())):
        raise ValueError('Original weights, zero updates, frozen encoder, or preflight binding changed')
    O.write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        selected_draw_ordinals=selected, matched_preflight_records=matched_preflight_records,
        species={s:describe([r for r in records if r['species'] == s]) for s in O.Q.ALLOW},
        pooled=describe(records), original_TRAIN_DNA_sources=descriptors,
        research_data_read=True, models_or_weights_read=True, native_forward_run=True,
        optimizer_constructed=False, optimizer_steps_this_run=0, checkpoints_written=False,
        original_C0_weights_unchanged=True, pretrained_encoder_frozen_no_grad=True,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        interpretation='fixed_checkpoint_diagnostic_not_AdamW_update_trajectory_or_arm_selection'))


if __name__ == '__main__':
    run()
