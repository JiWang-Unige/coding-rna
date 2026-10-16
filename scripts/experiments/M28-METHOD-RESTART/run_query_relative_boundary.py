#!/usr/bin/env python3
"""A23 fixed fresh fit and reference-free natural TRAIN pool generation."""
import argparse
from collections import Counter
import json
import math
import time
import torch
import preflight_query_relative_boundary as P
import diagnose_quality_train as D
from src.m28.labels import build_targets, candidate_targets

A, Q, C, ROOT = P.A, P.Q, P.C, P.ROOT
CONFIG, PROTOCOL, PREFLIGHT = P.CONFIG, P.PROTOCOL, P.OUT
FIT = ROOT/'outputs/M28-A23-RELATIVE-BOUNDARY-FIT'
READOUT = ROOT/'outputs/M28-A23-RELATIVE-BOUNDARY-TRAIN-READOUT'
new_model = P.new_model


def write(out, result):
    C.R.write_summary(out, dict(protocol=PROTOCOL, configuration=CONFIG, **result))


def validate_fit(value, state=False):
    if (value['arm'] != 'A23' or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
        value['step' if state else 'optimizer_steps'] != 4608 or
        value['revision' if state else 'feature_revision'] != A.F.REVISION):
        raise ValueError('Only fixed A23 last generator is permitted')
    if not state and (not value['completed'] or value['DEV_used_for_selection'] or value['test_setaria_used'] or
                      value['final_checkpoint'] != str((FIT/'step_004608.pt').relative_to(ROOT))):
        raise ValueError('No incomplete or selected A23 checkpoint')


def relative_gradients(model):
    return {name: math.sqrt(sum(float(p.grad.double().square().sum()) for p in layer.parameters()
                               if p.grad is not None))
            for name, layer in (('input', model.query.relative_boundary[0]),
                                ('output', model.query.relative_boundary[-1]))}


def fit():
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed'] or
        ready['configuration'] != CONFIG or not ready['exact_A3_common_fresh_parameters'] or
        not ready['zero_delta_actual_losses_and_assignment_verified'] or not ready['model_weights_unchanged'] or
        ready['optimizer_steps_this_run'] != 0 or ready['trainable_parameters'] != 649617):
        raise ValueError('Matching completed zero-update relative-boundary preflight required')
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population required')
    model = new_model().train()
    fresh_relative = {k: v.detach().clone() for k, v in model.query.relative_boundary.state_dict().items()}
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step = 0; epochs = []
    gradient_steps = Counter()
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            sums = Counter(); counts = {s: Counter() for s in A.ALLOW}; epoch_began = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                w = windows[draw['window_id']]
                cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
                features, dna = Q.window_inputs(cached, w, 'cuda')
                optimizer.zero_grad(set_to_none=True)
                outputs = model(features, dna, [w['valid_bases']])
                result = A.supervised_readout(model, outputs, w, catalogs[w['species']], cached['sequence'],
                                              draw['importance_weight_uniform_within_stratum'])
                result['total'].backward()
                relative = relative_gradients(model)
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1., error_if_nonfinite=True)
                optimizer.step(); step += 1
                gradient_steps.update({k: int(v > 0) for k, v in relative.items()})
                labels = result['query']['query_labels']; values = {k: float(v.detach()) for k, v in result['losses'].items()}
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=w['id'], species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'], losses=values,
                    weighted_total=float(result['total'].detach()), gradient_norm_before_clip=float(norm),
                    relative_layer_gradient_norms_before_clip=relative,
                    query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
                    protection_reasons=result['protection_reasons'], gene_to_query=result['query']['gene_to_query'])
                log.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update(dict(values, weighted_total=row['weighted_total']))
                counts[w['species']].update({'query_'+k: v for k, v in row['query_label_counts'].items()})
                counts[w['species']].update(result['protection_reasons'])
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608)), flush=True)
                if step % 256 == 0:
                    checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
                    with checkpoint.open('xb') as handle:
                        torch.save(dict(arm='A23', configuration=CONFIG, protocol=PROTOCOL, step=step,
                                        revision=A.F.REVISION, model=model.state_dict(), optimizer=optimizer.state_dict()), handle)
                del cached, features, dna, outputs, result
            epochs.append(dict(epoch=epoch+1, steps=1536, supervision_by_species={s: dict(c) for s, c in counts.items()},
                mean_losses={k: v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-epoch_began,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 4608 or any(gradient_steps[k] == 0 for k in ('input', 'output')):
        raise ValueError('Incomplete exposure or an inactive relative-boundary layer')
    changes = {k: float((v-fresh_relative[k]).double().norm())
               for k, v in model.query.relative_boundary.state_dict().items()}
    write(FIT, dict(arm='A23', completed=True, feature_revision=A.F.REVISION, optimizer_steps=step, epochs=epochs,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(ROOT)), main_evaluation='last_step_004608_only',
        trainable_parameters=649617, supervised_checkpoint_inheritance=False, new_backbone_forward=False,
        DEV_used_for_selection=False, test_setaria_used=False, wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), checkpoint_interval=256,
        relative_layer_nonzero_gradient_steps=dict(gradient_steps), relative_parameter_change_norms=changes))


def counts(rows):
    return dict(free_candidates=len(rows), correct=sum(r['label'] == 1 for r in rows),
        known_nonmatch=sum(r['label'] == 0 for r in rows), unknown=sum(r['label'] == -1 for r in rows),
        five_or_more_correct=sum(r['label'] == 1 and len(r['chain']) >= 5 for r in rows))


def gate(per_species):
    checks = {s+'_free_correct_increased': v['A23']['correct'] > v['A3']['correct']
              for s, v in per_species.items()}
    checks['pooled_free_five_or_more_correct_increased'] = (
        sum(v['A23']['five_or_more_correct'] for v in per_species.values()) >
        sum(v['A3']['five_or_more_correct'] for v in per_species.values()))
    return checks


@torch.inference_mode()
def readout():
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    if (FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal completed A23 fit required')
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True); validate_fit(state, state=True)
    model = new_model().eval().requires_grad_(False); model.load_state_dict(state['model'], strict=True); del state
    draws = C.R.rows(Q.DRAWS); ids = list(dict.fromkeys(d['window_id'] for d in draws))
    index = {r['window_id']: r for r in Q.train_records(Q.FEATURES/'index.jsonl') if r['window_id'] in ids}
    prior_fit = json.loads((C.FIT/'summary.json').read_text()); C.validate_head(prior_fit)
    if (C.FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Closed original A3 natural pools required')
    original_free = {r['window']['id']: r for r in C.R.rows(C.FIT/'free_pool.jsonl')}
    if len(draws) != 1536 or len(ids) != 1491 or set(index) != set(ids) or set(original_free) != set(ids):
        raise ValueError('Exactly original ordered natural TRAIN geometry required')
    free_records = {}; sequences = {}; began = time.perf_counter()
    # Labels are loaded only after the full raw-DNA candidate generation file closes.
    with (READOUT/'free_pool.jsonl').open('x') as log:
        for ordinal, ident in enumerate(ids):
            w = {k: original_free[ident]['window'][k] for k in C.R.GEOMETRY}
            if w['split'] != 'train' or w['seqid'] not in C.R.ALLOW[w['species']]:
                raise ValueError('Original allowed TRAIN geometry only')
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            features, dna = Q.window_inputs(cached, w, 'cuda')
            hidden = model(features, dna, [w['valid_bases']])['features'][0]
            chains, sources = C.free_pool(model, hidden, cached['sequence'])
            record = dict(window=w, chains=chains, decoder_sources=sources, reference_used_in_pool_generation=False)
            log.write(json.dumps(record, allow_nan=False)+'\n'); free_records[ident] = record
            sequences[ident] = cached['sequence']
            if (ordinal+1) % 32 == 0:
                log.flush(); print(json.dumps(dict(natural_TRAIN_windows=ordinal+1, total=1491)), flush=True)
            del cached, features, dna, hidden
    pool_seconds = time.perf_counter()-began
    loaded_draws, _, windows, catalogs = Q.load_train_inputs()
    if draws != loaded_draws:
        raise ValueError('Natural loop and target population differ')
    prior = json.loads((D.OUT/'summary.json').read_text())
    if ((D.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not prior['completed'] or
        prior['source_generator'] != C.GENERATOR or prior['windows'] != 1491):
        raise ValueError('Original A12 closed natural-label baseline required')
    old_records = {r['window']['id']: r for r in C.R.rows(D.OUT/'windows.jsonl')}
    totals = {s: {a: Counter() for a in ('A3', 'A23')} for s in A.ALLOW}
    with (READOUT/'windows.jsonl').open('x') as log:
        for ident in ids:
            w, free = windows[ident], free_records[ident]
            if any(w[k] != free['window'][k] or w[k] != old_records[ident]['window'][k] for k in C.R.GEOMETRY):
                raise ValueError('Natural/label/source geometry mismatch')
            targets = build_targets(w, catalogs[w['species']], sequences[ident])
            baseline = [r for r in old_records[ident]['candidates'] if r['origin'] == 'free']
            if ([r['chain'] for r in baseline] != original_free[ident]['chains'] or
                candidate_targets([r['chain'] for r in baseline], targets) != [r['label'] for r in baseline]):
                raise ValueError('Original actual natural candidates and labels changed')
            local = [dict(chain=chain, label=label) for chain, label in
                     zip(free['chains'], candidate_targets(free['chains'], targets))]
            values = dict(A3=counts(baseline), A23=counts(local))
            log.write(json.dumps(dict(window=free['window'], candidates=local, natural_counts=values), allow_nan=False)+'\n')
            for arm in ('A3', 'A23'):
                totals[w['species']][arm].update(values[arm])
    actual = {s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()}; checks = gate(actual)
    write(READOUT, dict(completed=True, arm='A23', checkpoint=fitted['final_checkpoint'], windows=1491,
        natural_counts_per_species=actual, original_A3_actual_natural_pool_labels_recomputed=True,
        predeclared_TRAIN_generator_gate=checks, TRAIN_generator_supported=all(checks.values()),
        reference_used_in_pool_generation=False, reference_injected_candidates=0, new_candidate_search=True,
        optimizer_steps_this_run=0, new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        pool_seconds=pool_seconds, seconds=time.perf_counter()-began,
        limitation='In-sample candidate/window support, not selected genes, F1 or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('fit', 'readout'), required=True)
    torch.set_num_threads(2)
    {'fit': fit, 'readout': readout}[parser.parse_args().mode]()
