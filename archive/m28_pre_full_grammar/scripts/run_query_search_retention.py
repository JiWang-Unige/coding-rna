#!/usr/bin/env python3
"""A30 one fresh fixed fit and original reference-free natural TRAIN pools."""
import argparse
from collections import Counter
import json
import time
import torch
import preflight_query_search_retention as P
import diagnose_quality_train as D
from src.m28.query_collaboration import CollaborativeGeneQueryModel
from src.m28.labels import build_targets, candidate_targets
from query_search_retention_training import search_retention_loss

A, Q, C, ROOT = P.A, P.Q, P.C, P.ROOT
CONFIG, PROTOCOL, PREFLIGHT = P.CONFIG, P.PROTOCOL, P.OUT
FIT = ROOT/'outputs/M28-A30-SEARCH-RETENTION-FIT'
READOUT = ROOT/'outputs/M28-A30-SEARCH-RETENTION-TRAIN-READOUT'


def new_model(device='cuda'):
    torch.manual_seed(CONFIG['seed'])
    return CollaborativeGeneQueryModel(CONFIG['feature_dim'], CONFIG['num_queries'], CONFIG['hidden']).to(device)


def write(out, result):
    C.R.write_summary(out, dict(protocol=PROTOCOL, configuration=CONFIG, **result))


def validate_fit(value, state=False):
    if (value['arm'] != 'A30' or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
        value['step' if state else 'optimizer_steps'] != 4608 or
        value['revision' if state else 'feature_revision'] != A.F.REVISION):
        raise ValueError('Only fixed A30 last generator is permitted')
    if not state and (not value['completed'] or value['DEV_used_for_selection'] or value['test_setaria_used'] or
                      value['final_checkpoint'] != str((FIT/'step_004608.pt').relative_to(ROOT))):
        raise ValueError('No incomplete or selected A30 checkpoint')


def fit():
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    flags = ('exact_complete_A3_fresh_parameters', 'original_four_losses_assignment_and_protection_verified',
             'original_search_geometry_order_scores_and_first_loss_verified',
             'finite_actual_combined_and_extra_gradients',
             'both_species_actual_shared_and_query_extra_gradients_verified', 'model_weights_unchanged')
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or
        (PREFLIGHT/'JOBID').read_text().strip() != '13398207' or not ready['completed'] or
        ready['configuration'] != CONFIG or ready['protocol'] != PROTOCOL or
        not all(ready[k] for k in flags) or ready['optimizer_steps_this_run'] != 0 or
        ready['trainable_parameters'] != 649457 or ready['cases'] != 8):
        raise ValueError('Matching completed original eight-case zero-update preflight required')
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population required')
    model = new_model().train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step = 0; epochs = []
    margin_steps = Counter()
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            sums = Counter(); counts = {s: Counter() for s in A.ALLOW}; epoch_began = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                started = time.perf_counter(); w = windows[draw['window_id']]
                cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
                features, dna = Q.window_inputs(cached, w, 'cuda')
                optimizer.zero_grad(set_to_none=True)
                outputs = model(features, dna, [w['valid_bases']])
                weight = draw['importance_weight_uniform_within_stratum']
                base = A.supervised_readout(model, outputs, w, catalogs[w['species']], cached['sequence'], weight)
                targets = Q.build_targets(w, catalogs[w['species']], cached['sequence'])
                window = model.query_window(outputs['features'][0], w['valid_bases'])
                t = time.perf_counter()
                extra = search_retention_loss(window, cached['sequence'], base['query']['positive_chains'],
                    targets, CONFIG['search_retention']['margin_log_score_units'])
                torch.cuda.synchronize(); collection_seconds = time.perf_counter()-t
                total = base['total']+weight*CONFIG['search_retention']['coefficient']*extra['loss']
                if not bool(torch.isfinite(total)):
                    raise ValueError('Nonfinite fixed A30 objective')
                total.backward()
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1., error_if_nonfinite=True)
                optimizer.step(); step += 1
                labels = base['query']['query_labels']; values = {k: float(v.detach()) for k, v in base['losses'].items()}
                margin = float(extra['loss'].detach()); margin_steps[w['species']] += int(margin > 0)
                torch.cuda.synchronize()
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=w['id'], species=w['species'],
                    importance_weight=weight, losses=values, weighted_original_total=float(base['total'].detach()),
                    search_retention_margin=margin, weighted_total=float(total.detach()),
                    gradient_norm_before_clip=float(norm), positive_queries=extra['positive_queries'],
                    collection_records=extra['records'], reference_injected_candidates=extra['reference_injected_candidates'],
                    query_label_counts={str(y): labels.count(y) for y in (-1, 0, 1)},
                    protection_reasons=base['protection_reasons'], gene_to_query=base['query']['gene_to_query'],
                    current_search_collection_and_live_score_seconds=collection_seconds,
                    load_forward_supervision_search_backward_optimizer_seconds=time.perf_counter()-started)
                log.write(json.dumps(row, allow_nan=False)+'\n')
                sums.update(dict(values, search_retention_margin=margin, weighted_total=row['weighted_total']))
                counts[w['species']].update({'query_'+k: v for k, v in row['query_label_counts'].items()})
                counts[w['species']].update(base['protection_reasons'])
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608, margin=margin,
                                                       seconds=time.perf_counter()-began)), flush=True)
                if step % 256 == 0:
                    checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
                    with checkpoint.open('xb') as handle:
                        torch.save(dict(arm='A30', configuration=CONFIG, protocol=PROTOCOL, step=step,
                            revision=A.F.REVISION, model=model.state_dict(), optimizer=optimizer.state_dict()), handle)
                del cached, features, dna, outputs, base, targets, window, extra, total
            epochs.append(dict(epoch=epoch+1, steps=1536, supervision_by_species={s: dict(c) for s, c in counts.items()},
                mean_losses={k: v/1536 for k, v in sums.items()}, seconds=time.perf_counter()-epoch_began,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 4608 or any(margin_steps[s] == 0 for s in A.ALLOW):
        raise ValueError('Incomplete exposure or both-species current-search signal missing')
    write(FIT, dict(arm='A30', completed=True, feature_revision=A.F.REVISION, optimizer_steps=step, epochs=epochs,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(ROOT)), main_evaluation='last_step_004608_only',
        trainable_parameters=649457, supervised_checkpoint_inheritance=False, new_backbone_forward=False,
        reference_injected_candidates=0, DEV_used_for_selection=False, test_setaria_used=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        checkpoint_interval=256, positive_search_margin_steps_by_species=dict(margin_steps)))


def counts(rows):
    return dict(free_candidates=len(rows), correct=sum(r['label'] == 1 for r in rows),
        known_nonmatch=sum(r['label'] == 0 for r in rows), unknown=sum(r['label'] == -1 for r in rows),
        five_or_more_correct=sum(r['label'] == 1 and len(r['chain']) >= 5 for r in rows))


def gate(per_species):
    checks = {s+'_free_correct_increased': v['A30']['correct'] > v['A3']['correct']
              for s, v in per_species.items()}
    checks['pooled_free_five_or_more_correct_increased'] = (
        sum(v['A30']['five_or_more_correct'] for v in per_species.values()) >
        sum(v['A3']['five_or_more_correct'] for v in per_species.values()))
    checks.update({s+'_free_known_nonmatch_not_increased': v['A30']['known_nonmatch'] <= v['A3']['known_nonmatch']
                   for s, v in per_species.items()})
    return checks


@torch.inference_mode()
def readout():
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    if (FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal completed A30 fit required')
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
    # Complete the reference-free candidate file before loading any TRAIN catalog.
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
    totals = {s: {a: Counter() for a in ('A3', 'A30')} for s in A.ALLOW}
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
            values = dict(A3=counts(baseline), A30=counts(local))
            log.write(json.dumps(dict(window=free['window'], candidates=local, natural_counts=values), allow_nan=False)+'\n')
            for arm in ('A3', 'A30'):
                totals[w['species']][arm].update(values[arm])
    actual = {s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()}; checks = gate(actual)
    write(READOUT, dict(completed=True, arm='A30', checkpoint=fitted['final_checkpoint'], windows=1491,
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
