#!/usr/bin/env python3
"""A14 fixed candidate-set objective and natural-pool TRAIN selection readout."""
import argparse
from collections import Counter
import json
import math
import time
import torch
import run_chain_competition as P
from src.m28.chain_set import chain_set_loss
from src.m28.core import select_nonoverlapping

C, D, ROOT = P.C, P.D, P.ROOT
PREFLIGHT = ROOT/'outputs/M28-A14-CHAIN-SET-PREFLIGHT'
FIT = ROOT/'outputs/M28-A14-CHAIN-SET-FIT'
READOUT = ROOT/'outputs/M28-A14-CHAIN-SET-TRAIN-READOUT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a14_chain_set_protocol.md'
CONFIG = {**C.CONFIG, 'loss': 'exact_nonoverlap_candidate_set_NLL_divided_by_known_candidate_count',
          'pool': 'same_immutable_A6_free_pool_plus_same_TRAIN_positive_injection',
          'gold': 'all_known_positive_chains_as_one_compatible_set',
          'unknown': 'excluded_from_partition_and_gold_never_a_negative',
          'initialization': 'same_A6_seed0_from_scratch',
          'main_evaluation': 'last_head_step_004608_only_natural_pool_TRAIN_selection'}


def population():
    draws, index, windows, records, rows = P.load_cohort()
    # Reject a genuinely incompatible supervision domain, never drop positives
    # or relabel protected chains merely to make the structured target feasible.
    for rec in records.values():
        local = rec['candidates']
        chain_set_loss({'logits': torch.zeros(len(local), dtype=torch.float64), 'null': torch.tensor(0.)},
                       [r['chain'] for r in local], [r['label'] for r in local])
    return draws, index, windows, records, rows


def natural_selection(rows):
    free = [r for r in rows if r['origin'] == 'free']
    chosen = select_nonoverlapping([r['chain'] for r in free], [r['gain'] for r in free])
    selected = [free[i] for i in chosen]
    return dict(free_candidates=len(free), free_correct_candidates=sum(r['label'] == 1 for r in free),
        free_known_nonmatch_candidates=sum(r['label'] == 0 for r in free),
        free_unknown_candidates=sum(r['label'] == -1 for r in free), selected=len(selected),
        selected_correct=sum(r['label'] == 1 for r in selected),
        selected_known_nonmatch=sum(r['label'] == 0 for r in selected),
        selected_unknown=sum(r['label'] == -1 for r in selected),
        selected_five_or_more_correct=sum(r['label'] == 1 and len(r['chain']) >= 5 for r in selected))


def preflight():
    draws, index, windows, records, old_rows = population()
    baseline = {s: Counter() for s in C.R.ALLOW}
    for rec in records.values():
        baseline[rec['window']['species']].update(natural_selection(rec['candidates']))
    generator = C.load_generator(); head = C.new_quality('cuda').train()
    initial = {k: v.detach().clone() for k, v in head.named_parameters()}
    witnesses = P.select_gradient_windows(records)
    began = time.perf_counter(); receipts = []
    for witness in witnesses:
        ident = witness['window_id']; w = windows[ident]; local = records[ident]['candidates']
        hidden = P.inputs(generator, index, w)
        head.zero_grad(set_to_none=True)
        score = head(hidden, [r['chain'] for r in local], w['valid_bases'])
        result = chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
        if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(result['loss'])):
            raise ValueError('Nonfinite actual TRAIN score/structured loss')
        result['loss'].backward()
        groups = {}
        for name in ('exon', 'sequence', 'additive', 'nonadditive', 'null'):
            params = list(getattr(head, name).parameters())
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite actual gradient: '+name)
            groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if groups[name] <= 0:
                raise ValueError('Zero actual gradient: '+name)
        if any(p.grad is not None for p in generator.parameters()) or any(p.grad is not None for p in head.link.parameters()):
            raise ValueError('Frozen generator/link received gradient')
        receipts.append(dict(**witness, gradient_norms=groups, loss=float(result['loss'].detach()),
                             **{k: v for k, v in result.items() if not torch.is_tensor(v)}))
    if any(not torch.equal(p.detach(), initial[k]) for k, p in head.named_parameters()):
        raise ValueError('Zero-update preflight changed head weights')
    C.R.write_summary(PREFLIGHT, dict(experiment='M28-A14-CHAIN-SET-PREFLIGHT', completed=True,
        configuration=CONFIG, protocol=PROTOCOL, original_TRAIN_windows=len(records),
        candidate_window_pairs=len(old_rows), original_draws=len(draws),
        all_gold_sets_compatible=True, original_A12_actual_records_verified=True,
        A6_natural_selection_per_species={s: dict(v) for s, v in baseline.items()},
        quality_trainable_parameters=sum(p.numel() for p in head.parameters() if p.requires_grad),
        optimizer_steps_this_run=0, head_weights_unchanged=True, gradient_witnesses=receipts,
        new_candidate_search=False, new_backbone_forward=False, DEV_forward_used=False,
        test_setaria_used=False, seconds=time.perf_counter()-began))


def validate_fit(value, state=False):
    if (value['arm'] != 'A14' or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
        value['step' if state else 'optimizer_steps'] != 4608 or
        value['revision' if state else 'feature_revision'] != C.I.F.REVISION):
        raise ValueError('Only fixed A14 last head is permitted')
    if not state and (not value['completed'] or value['DEV_used_for_selection'] or value['test_setaria_used'] or
                     value['final_checkpoint'] != str((FIT/'step_004608.pt').relative_to(ROOT))):
        raise ValueError('No selected or incomplete A14 checkpoint')


def fit():
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed'] or
        ready['configuration'] != CONFIG or not ready['all_gold_sets_compatible'] or
        not ready['head_weights_unchanged'] or ready['optimizer_steps_this_run'] != 0):
        raise ValueError('Completed unchanged actual-data structured preflight required')
    draws, index, windows, records, _ = population()
    generator = C.load_generator(); head = C.new_quality('cuda').train()
    params = [p for p in head.parameters() if p.requires_grad]
    if sum(p.numel() for p in params) != 165635:
        raise ValueError('Same A6 head parameter count required')
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    step = 0; epochs = []
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts = Counter(); total = 0.; started = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; w = windows[ident]; rec = records[ident]
                local = rec['candidates']; hidden = P.inputs(generator, index, w)
                optimizer.zero_grad(set_to_none=True)
                score = head(hidden, [r['chain'] for r in local], w['valid_bases'])
                result = chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A14 score/loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                if any(p.grad is not None for p in generator.parameters()) or any(p.grad is not None for p in head.link.parameters()):
                    raise ValueError('Frozen generator/link received gradient')
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident, species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'], weighted_loss=float(weighted.detach()),
                    log_partition=float(result['log_partition'].detach()), gold_gain=float(result['gold_gain'].detach()),
                    gradient_norm_before_clip=float(norm), **rec['counts'])
                log.write(json.dumps(row, allow_nan=False)+'\n'); counts.update(rec['counts']); total += row['weighted_loss']
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608)), flush=True)
            checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(arm='A14', configuration=CONFIG, protocol=PROTOCOL, step=step,
                            revision=C.I.F.REVISION, head=head.state_dict(), optimizer=optimizer.state_dict()), checkpoint)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                               mean_weighted_loss=total/1536, seconds=time.perf_counter()-started))
    if step != 4608:
        raise ValueError('Incomplete fixed A14 exposure')
    C.R.write_summary(FIT, dict(arm='A14', completed=True, configuration=CONFIG, protocol=PROTOCOL,
        feature_revision=C.I.F.REVISION, feature_spec=C.R.FEATURE_SPEC, optimizer_steps=step,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(ROOT)), epochs=epochs,
        source_free_pool=str((C.FIT/'free_pool.jsonl').relative_to(ROOT)), source_TRAIN_cohort=str(D.OUT.relative_to(ROOT)),
        generator_trainable_parameters=0, quality_trainable_parameters=165635, DEV_used_for_selection=False,
        test_setaria_used=False, new_candidate_search=False, new_backbone_forward=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One fixed TRAIN mechanism candidate, not final raw annotation or independent generalization.'))


@torch.inference_mode()
def readout():
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    if (FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal fit required')
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    head = C.new_quality('cuda').eval().requires_grad_(False); head.load_state_dict(state['head'], strict=True)
    del state
    generator = C.load_generator(); _, index, windows, records, old_rows = population()
    baseline = {s: Counter() for s in C.R.ALLOW}; candidate = {s: Counter() for s in C.R.ALLOW}
    new_rows = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    with (READOUT/'windows.jsonl').open('x') as log:
        for ordinal, (ident, rec) in enumerate(records.items()):
            w = windows[ident]; hidden = P.inputs(generator, index, w)
            score = head(hidden, [r['chain'] for r in rec['candidates']], w['valid_bases'])
            values = {k: score[k].float().cpu().tolist() for k in ('logits', 'additive', 'nonadditive')}
            null = float(score['null']); local = []
            for i, prior in enumerate(rec['candidates']):
                row = {**prior, 'gain': values['logits'][i], 'additive': values['additive'][i],
                       'nonadditive': values['nonadditive'][i], 'null': null}
                if not all(math.isfinite(row[k]) for k in D.COMPONENTS):
                    raise ValueError('Nonfinite last-head TRAIN score')
                local.append(row)
            for row in local:
                row['max_overlapping_known_nonmatch_gain'] = D.overlapping_max(row['chain'], local) if row['label'] == 1 else None
            a, b = natural_selection(rec['candidates']), natural_selection(local)
            baseline[w['species']].update(a); candidate[w['species']].update(b)
            log.write(json.dumps(dict(window=w, counts=rec['counts'], candidates=local,
                natural_selection={'A6': a, 'A14': b}), allow_nan=False)+'\n'); new_rows.extend(local)
            if (ordinal+1) % 64 == 0:
                log.flush(); print(json.dumps(dict(scored_TRAIN_windows=ordinal+1, total=1491)), flush=True)
    checks = {}
    for species in baseline:
        a, b = baseline[species], candidate[species]
        for name in ('free_candidates', 'free_correct_candidates', 'free_known_nonmatch_candidates', 'free_unknown_candidates'):
            if a[name] != b[name]:
                raise ValueError('Same natural population required')
        checks[species+'_natural_selected_correct_increased'] = b['selected_correct'] > a['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = b['selected_known_nonmatch'] <= a['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_increased'] = (
        sum(v['selected_five_or_more_correct'] for v in candidate.values()) >
        sum(v['selected_five_or_more_correct'] for v in baseline.values()))
    C.R.write_summary(READOUT, dict(experiment='M28-A14-CHAIN-SET-TRAIN-READOUT', completed=True,
        configuration=CONFIG, protocol=PROTOCOL, checkpoint=fitted['final_checkpoint'], windows=1491,
        candidate_window_pairs=len(new_rows), source_generator=C.GENERATOR,
        natural_selection_per_species={s: {'A6': dict(baseline[s]), 'A14': dict(candidate[s])} for s in baseline},
        score_profiles={'A6': P.compact(old_rows), 'A14': P.compact(new_rows)},
        predeclared_TRAIN_preparation_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        candidate_geometry_labels_and_origins_unchanged=True, optimizer_steps_this_run=0,
        new_candidate_search=False, new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        threshold=0., threshold_scan=False, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Natural selection is per TRAIN window, no injected chain or label filter in selection. Counts are window/candidate pairs, not independent gene F1, final genome background or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preflight', 'fit', 'readout'), required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    {'preflight': preflight, 'fit': fit, 'readout': readout}[args.mode]()
