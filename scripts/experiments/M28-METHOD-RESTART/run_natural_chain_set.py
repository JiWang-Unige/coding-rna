#!/usr/bin/env python3
"""A15 fixed natural-domain set fit and matched A6/A14 TRAIN readout."""
import argparse
from collections import Counter
import json
import math
import time
import torch
import run_chain_set as S
from src.m28.chain_set import chain_set_loss
from src.m28.natural_chain_set import natural_record

P, C, D, ROOT = S.P, S.C, S.D, S.ROOT
PREFLIGHT = ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-PREFLIGHT'
FIT = ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-FIT'
READOUT = ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-TRAIN-READOUT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a15_natural_chain_set_protocol.md'
CONFIG = {**S.CONFIG, 'pool': 'same_immutable_A6_free_pool_without_reference_injection',
          'gold': 'all_natural_known_positive_chains_as_one_compatible_set',
          'main_evaluation': 'last_head_step_004608_only_natural_TRAIN_A6_gate_A14_pair'}
IDENTITY = ('window_id', 'species', 'chain', 'origin', 'label', 'exon_class')


def population():
    draws, index, windows, source, _ = P.load_cohort()
    records = {ident: natural_record(rec) for ident, rec in source.items()}
    for rec in records.values():
        local = rec['candidates']
        chain_set_loss(dict(logits=torch.zeros(len(local), dtype=torch.float64), null=torch.tensor(0.)),
                       [r['chain'] for r in local], [r['label'] for r in local])
    return draws, index, windows, records, source


def totals(records, ids):
    counts = Counter()
    for ident in ids:
        counts.update(records[ident]['counts'])
    return dict(counts)


def preflight():
    draws, index, windows, records, source = population()
    generator = C.load_generator(); head = C.new_quality('cuda').train()
    initial = {k: v.detach().clone() for k, v in head.named_parameters()}
    receipts = []; began = time.perf_counter()
    for witness in P.select_gradient_windows(records):
        ident = witness['window_id']; w = windows[ident]; local = records[ident]['candidates']
        head.zero_grad(set_to_none=True)
        score = head(P.inputs(generator, index, w), [r['chain'] for r in local], w['valid_bases'])
        result = chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
        if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(result['loss'])):
            raise ValueError('Nonfinite natural-domain preflight score/loss')
        result['loss'].backward(); groups = {}
        for name in ('exon', 'sequence', 'additive', 'nonadditive', 'null'):
            params = list(getattr(head, name).parameters())
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite natural-domain gradient: '+name)
            groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if groups[name] <= 0:
                raise ValueError('Zero natural-domain gradient: '+name)
        if any(p.grad is not None for p in generator.parameters()) or any(p.grad is not None for p in head.link.parameters()):
            raise ValueError('Frozen generator/link received gradient')
        receipts.append(dict(**witness, gradient_norms=groups, loss=float(result['loss'].detach()),
                             **{k: v for k, v in result.items() if not torch.is_tensor(v)}))
    if any(not torch.equal(p.detach(), initial[k]) for k, p in head.named_parameters()):
        raise ValueError('Preflight updated the head')
    C.R.write_summary(PREFLIGHT, dict(experiment='M28-A15-NATURAL-CHAIN-SET-PREFLIGHT', completed=True,
        configuration=CONFIG, protocol=PROTOCOL, original_A12_actual_records_verified=True,
        all_natural_gold_sets_compatible=True, original_draws=len(draws), windows=len(records),
        source_counts=totals(source, source), natural_counts=totals(records, records),
        replayed_draw_counts=totals(records, [d['window_id'] for d in draws]),
        windows_losing_all_gold=sum(bool(source[i]['counts']['label_1']) and
                                   not records[i]['counts']['label_1'] for i in records),
        gradient_witnesses=receipts, optimizer_steps_this_run=0, head_weights_unchanged=True,
        quality_trainable_parameters=sum(p.numel() for p in head.parameters() if p.requires_grad),
        new_candidate_search=False, new_backbone_forward=False, DEV_forward_used=False,
        test_setaria_used=False, seconds=time.perf_counter()-began))


def validate_fit(value, state=False):
    if (value['arm'] != 'A15' or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
        value['step' if state else 'optimizer_steps'] != 4608 or
        value['revision' if state else 'feature_revision'] != C.I.F.REVISION):
        raise ValueError('Only the fixed A15 last head is permitted')
    if not state and (not value['completed'] or value['DEV_used_for_selection'] or value['test_setaria_used'] or
                     value['final_checkpoint'] != str((FIT/'step_004608.pt').relative_to(ROOT))):
        raise ValueError('No selected or incomplete A15 checkpoint')


def fit():
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed'] or
        ready['configuration'] != CONFIG or not ready['all_natural_gold_sets_compatible'] or
        not ready['head_weights_unchanged'] or ready['optimizer_steps_this_run'] != 0):
        raise ValueError('Completed unchanged natural-domain preflight required')
    draws, index, windows, records, _ = population()
    generator = C.load_generator(); head = C.new_quality('cuda').train()
    params = [p for p in head.parameters() if p.requires_grad]
    if sum(p.numel() for p in params) != 165635:
        raise ValueError('Same head parameter count required')
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step = 0; epochs = []
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts = Counter(); total = 0.; started = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; w = windows[ident]; local = records[ident]['candidates']
                optimizer.zero_grad(set_to_none=True)
                score = head(P.inputs(generator, index, w), [r['chain'] for r in local], w['valid_bases'])
                result = chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A15 score/loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                if any(p.grad is not None for p in generator.parameters()) or any(p.grad is not None for p in head.link.parameters()):
                    raise ValueError('Frozen generator/link received gradient')
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident, species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'], weighted_loss=float(weighted.detach()),
                    log_partition=float(result['log_partition'].detach()), gold_gain=float(result['gold_gain'].detach()),
                    gradient_norm_before_clip=float(norm), **records[ident]['counts'])
                log.write(json.dumps(row, allow_nan=False)+'\n'); counts.update(records[ident]['counts']); total += row['weighted_loss']
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608)), flush=True)
            checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(arm='A15', configuration=CONFIG, protocol=PROTOCOL, step=step,
                            revision=C.I.F.REVISION, head=head.state_dict(), optimizer=optimizer.state_dict()), checkpoint)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                               mean_weighted_loss=total/1536, seconds=time.perf_counter()-started))
    if step != 4608:
        raise ValueError('Incomplete fixed exposure')
    C.R.write_summary(FIT, dict(arm='A15', completed=True, configuration=CONFIG, protocol=PROTOCOL,
        feature_revision=C.I.F.REVISION, feature_spec=C.R.FEATURE_SPEC, optimizer_steps=step,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(ROOT)), epochs=epochs,
        source_free_pool=str((C.FIT/'free_pool.jsonl').relative_to(ROOT)), source_TRAIN_cohort=str(D.OUT.relative_to(ROOT)),
        generator_trainable_parameters=0, quality_trainable_parameters=165635, DEV_used_for_selection=False,
        test_setaria_used=False, new_candidate_search=False, new_backbone_forward=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='A14-matched draw/optimizer exposure, not matched positive exposure; natural gold only.'))


def a14_comparator(records):
    summary = json.loads((S.READOUT/'summary.json').read_text())
    checked = ROOT/'outputs/M28-A14-CHAIN-SET-RESULT-CHECK'
    if ((S.READOUT/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed'] or
        summary['configuration'] != S.CONFIG or (checked/'STATUS').read_text().strip() != 'COMPLETED' or
        (checked/'JOBID').read_text().strip() != '13367225'):
        raise ValueError('Closed verified A14 paired comparator required')
    old = C.R.rows(S.READOUT/'windows.jsonl')
    if [r['window']['id'] for r in old] != list(records):
        raise ValueError('A14 ordered window population differs')
    paired = {}
    for rec in old:
        ident = rec['window']['id']; local = natural_record(rec)['candidates']
        if ([{k: r[k] for k in IDENTITY} for r in local] !=
            [{k: r[k] for k in IDENTITY} for r in records[ident]['candidates']]):
            raise ValueError('A14 natural geometry/labels/origins differ')
        paired[ident] = local
    return paired


@torch.inference_mode()
def readout():
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    if (FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal fit required')
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    head = C.new_quality('cuda').eval().requires_grad_(False); head.load_state_dict(state['head'], strict=True); del state
    generator = C.load_generator(); _, index, windows, records, _ = population()
    paired = a14_comparator(records)
    selected = {s: {arm: Counter() for arm in ('A6', 'A14', 'A15')} for s in C.R.ALLOW}
    new_rows = []; old_rows = []; a14_rows = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    with (READOUT/'windows.jsonl').open('x') as log:
        for ordinal, (ident, rec) in enumerate(records.items()):
            w = windows[ident]; score = head(P.inputs(generator, index, w), [r['chain'] for r in rec['candidates']], w['valid_bases'])
            values = {k: score[k].float().cpu().tolist() for k in ('logits', 'additive', 'nonadditive')}
            null = float(score['null']); local = []
            for i, prior in enumerate(rec['candidates']):
                row = {**prior, 'gain': values['logits'][i], 'additive': values['additive'][i],
                       'nonadditive': values['nonadditive'][i], 'null': null}
                if (not all(math.isfinite(row[k]) for k in D.COMPONENTS) or
                    abs(row['gain']-(row['additive']+row['nonadditive']-row['null'])) > 1e-4):
                    raise ValueError('Nonfinite/inconsistent last-head score components')
                local.append(row)
            for row in local:
                row['max_overlapping_known_nonmatch_gain'] = D.overlapping_max(row['chain'], local) if row['label'] == 1 else None
            comparison = {arm: S.natural_selection(rows) for arm, rows in
                          (('A6', rec['candidates']), ('A14', paired[ident]), ('A15', local))}
            for arm, counts in comparison.items():
                selected[w['species']][arm].update(counts)
            log.write(json.dumps(dict(window=w, counts=rec['counts'], candidates=local,
                natural_selection=comparison), allow_nan=False)+'\n')
            new_rows.extend(local); old_rows.extend(rec['candidates']); a14_rows.extend(paired[ident])
            if (ordinal+1) % 64 == 0:
                log.flush(); print(json.dumps(dict(scored_TRAIN_windows=ordinal+1, total=1491)), flush=True)
    checks = {}
    for species, arms in selected.items():
        a, b = arms['A6'], arms['A15']
        for name in ('free_candidates', 'free_correct_candidates', 'free_known_nonmatch_candidates', 'free_unknown_candidates'):
            if not a[name] == arms['A14'][name] == b[name]:
                raise ValueError('Same natural population required')
        checks[species+'_natural_selected_correct_increased'] = b['selected_correct'] > a['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = b['selected_known_nonmatch'] <= a['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_increased'] = (
        sum(v['A15']['selected_five_or_more_correct'] for v in selected.values()) >
        sum(v['A6']['selected_five_or_more_correct'] for v in selected.values()))
    C.R.write_summary(READOUT, dict(experiment='M28-A15-NATURAL-CHAIN-SET-TRAIN-READOUT', completed=True,
        configuration=CONFIG, protocol=PROTOCOL, checkpoint=fitted['final_checkpoint'], windows=len(records),
        candidate_window_pairs=len(new_rows), source_generator=C.GENERATOR,
        natural_selection_per_species={s: {a: dict(c) for a, c in v.items()} for s, v in selected.items()},
        score_profiles={'A6': P.compact(old_rows), 'A14': P.compact(a14_rows), 'A15': P.compact(new_rows)},
        predeclared_TRAIN_preparation_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        natural_candidate_geometry_labels_and_origins_unchanged=True, reference_injected_candidates=0,
        optimizer_steps_this_run=0, new_candidate_search=False, new_backbone_forward=False,
        DEV_forward_used=False, test_setaria_used=False, threshold=0., threshold_scan=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='TRAIN window/candidate pairs, not gene F1, raw-genome background, generalization or causal proof.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preflight', 'fit', 'readout'), required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    {'preflight': preflight, 'fit': fit, 'readout': readout}[args.mode]()
