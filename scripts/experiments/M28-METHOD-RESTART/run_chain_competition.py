#!/usr/bin/env python3
"""A13 fixed matched-head competition fit and TRAIN-only last-head readout."""
import argparse
import json
import math
import time
from collections import Counter, defaultdict
import torch
import run_chain_quality as C
import diagnose_quality_train as D
from src.m28.chain_competition import chain_competition_loss

ROOT = C.ROOT
PREFLIGHT = ROOT/'outputs/M28-A13-CHAIN-COMPETITION-PREFLIGHT'
FIT = ROOT/'outputs/M28-A13-CHAIN-COMPETITION-FIT'
READOUT = ROOT/'outputs/M28-A13-CHAIN-COMPETITION-TRAIN-READOUT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a13_chain_competition_protocol.md'
CONFIG = {**C.CONFIG, 'loss': 'local_positive_vs_null_and_known_span_conflicts_plus_negative_null_anchor',
          'pool': 'reuse_immutable_A6_free_pool_plus_same_TRAIN_positive_injection',
          'positive_negative_group_mass': 'equal_over_present_groups_per_window',
          'conflicts': 'strict_span_overlap_known_label0_only',
          'temperature': 1., 'margin': 0., 'initialization': 'same_A6_seed0_from_scratch',
          'main_evaluation': 'last_head_step_004608_only_TRAIN_then_separate_DEV_protocol'}


def compact(rows):
    return {name: D.profile(values) for name, values in {
        'positive': [r for r in rows if r['label'] == 1],
        'free_positive': [r for r in rows if r['label'] == 1 and r['origin'] == 'free'],
        'known_nonmatch': [r for r in rows if r['label'] == 0],
        'five_or_more_positive': [r for r in rows if r['label'] == 1 and r['exon_class'] == 'five_or_more'],
        'unknown': [r for r in rows if r['label'] == -1]}.items()}


def load_cohort():
    source = json.loads((D.OUT/'summary.json').read_text())
    if ((D.OUT/'STATUS').read_text().strip() != 'COMPLETED' or
        (D.OUT/'JOBID').read_text().strip() != '13360650' or
        not source['completed'] or source['source_generator'] != C.GENERATOR or
        source['source_head'] != 'outputs/M28-A6-CHAIN-QUALITY-FIT/step_004608.pt' or
        source['windows'] != 1491 or source['original_draws'] != 1536 or
        source['optimizer_steps_this_run'] != 0 or source['DEV_forward_used'] or source['test_setaria_used']):
        raise ValueError('Completed original zero-update TRAIN diagnosis required')
    draws, index, windows, catalogs = C.R.Q.load_train_inputs()
    records = C.R.rows(D.OUT/'windows.jsonl')
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    if [r['window']['id'] for r in records] != ids:
        raise ValueError('Exactly original 1491 ordered unique TRAIN windows required')
    free = {r['window']['id']: r for r in C.R.rows(C.FIT/'free_pool.jsonl')}
    rows = []
    by_id = {}
    replay = Counter()
    for rec in records:
        ident = rec['window']['id']; w = windows[ident]; local = rec['candidates']
        if any(w[k] != rec['window'][k] or w[k] != free[ident]['window'][k] for k in C.R.GEOMETRY):
            raise ValueError('Original feature/label/free-pool geometry changed')
        nfree = len(free[ident]['chains'])
        if ([r['chain'] for r in local[:nfree]] != free[ident]['chains'] or
            any(r['origin'] != 'free' for r in local[:nfree]) or
            any(r['origin'] != 'reference_injected' or r['label'] != 1 for r in local[nfree:])):
            raise ValueError('Same free/injected TRAIN candidate population required')
        for r in local:
            if (r['window_id'] != ident or r['species'] != w['species'] or
                r['label'] not in (-1, 0, 1) or r['exon_class'] != D.exon_class(r['chain']) or
                not all(math.isfinite(r[k]) for k in D.COMPONENTS) or
                abs(r['gain']-(r['additive']+r['nonadditive']-r['null'])) > 1e-4 or
                r['max_overlapping_known_nonmatch_gain'] !=
                    (D.overlapping_max(r['chain'], local) if r['label'] == 1 else None)):
                raise ValueError('Original A12 candidate/component/overlap identity changed')
        labels = [r['label'] for r in local]
        counts = dict(windows=1, free_candidates=nfree, reference_injected_candidates=len(local)-nfree,
                      zero_known_windows=int(not any(y >= 0 for y in labels)),
                      **{'label_'+str(y): labels.count(y) for y in (1, 0, -1)})
        if counts != rec['counts']:
            raise ValueError('A12 actual record count mismatch')
        by_id[ident] = rec; rows.extend(local)
    for d in draws:
        replay.update(by_id[d['window_id']]['counts'])
    if dict(replay) != source['replayed_original_draw_counts']:
        raise ValueError('Do not change original A6 draw/label counts')
    for y in (-1, 0, 1):
        actual = D.profile([r for r in rows if r['label'] == y])
        if actual != source['profiles']['label:'+str(y)]:
            raise ValueError('A12 primary label profile must match actual closed records')
    return draws, index, windows, by_id, rows


def inputs(generator, index, w):
    cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
    features, dna = C.R.Q.window_inputs(cached, w, 'cuda')
    with torch.no_grad():
        hidden = generator(features, dna, [w['valid_bases']])['features'][0]
    return hidden


def select_gradient_windows(records):
    selected = []
    for species in C.R.ALLOW:
        local = [(ordinal, ident, r) for ordinal, (ident, r) in enumerate(records.items())
                 if r['window']['species'] == species and r['counts']['label_1'] and r['counts']['label_0']]
        used = set()
        for category in ('most_positive', 'most_positive_exons', 'most_span_conflicts', 'positive_protected'):
            eligible = [p for p in local if p[1] not in used and
                        (category != 'positive_protected' or p[2]['counts']['label_-1'])]
            def priority(p):
                rows = p[2]['candidates']
                if category == 'most_positive':
                    return p[2]['counts']['label_1']
                if category == 'most_positive_exons':
                    return sum(len(r['chain']) for r in rows if r['label'] == 1)
                if category == 'most_span_conflicts':
                    return sum(n['chain'][0][0] < r['chain'][-1][1] and
                               n['chain'][-1][1] > r['chain'][0][0]
                               for r in rows if r['label'] == 1 for n in rows if n['label'] == 0)
                return 0
            if not eligible:
                raise ValueError('Missing original-draw gradient category: '+species+' '+category)
            chosen = min(eligible, key=lambda p: (-priority(p), p[0]))
            used.add(chosen[1])
            selected.append(dict(window_id=chosen[1], category=category, species=species))
    return selected


def preflight():
    draws, index, windows, records, old_rows = load_cohort()
    generator = C.load_generator()
    head = C.new_quality('cuda').train()
    initial = {k: v.detach().clone() for k, v in head.named_parameters()}
    witnesses = select_gradient_windows(records)
    began = time.perf_counter(); receipts = []
    for witness in witnesses:
        ident = witness['window_id']; w = windows[ident]; rec = records[ident]
        hidden = inputs(generator, index, w)
        chains = [r['chain'] for r in rec['candidates']]
        labels = [r['label'] for r in rec['candidates']]
        head.zero_grad(set_to_none=True)
        score = head(hidden, chains, w['valid_bases'])
        result = chain_competition_loss(score, chains, labels)
        if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(result['loss'])):
            raise ValueError('Nonfinite real-data score/loss')
        result['loss'].backward()
        groups = {}
        for name in ('exon', 'sequence', 'additive', 'nonadditive', 'null'):
            ps = list(getattr(head, name).parameters())
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in ps):
                raise ValueError('Missing/nonfinite gradient: '+name)
            groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in ps))
            if groups[name] <= 0:
                raise ValueError('Zero real-data gradient: '+name)
        if any(p.grad is not None for p in generator.parameters()) or any(p.grad is not None for p in head.link.parameters()):
            raise ValueError('Frozen generator/link received gradients')
        receipts.append(dict(window_id=ident, category=witness['category'], gradient_norms=groups,
                             **{k: v for k, v in result.items() if not torch.is_tensor(v)}))
    if any(not torch.equal(v.detach(), initial[k]) for k, v in head.named_parameters()):
        raise ValueError('Preflight must not update head weights')
    C.R.write_summary(PREFLIGHT, dict(experiment='M28-A13-CHAIN-COMPETITION-PREFLIGHT', completed=True,
        configuration=CONFIG, protocol=PROTOCOL, windows=8, A12_actual_records_verified=True,
        A12_candidate_window_pairs=len(old_rows), A12_replayed_original_draw_counts=
            {k: sum(records[d['window_id']]['counts'][k] for d in draws) for k in records[draws[0]['window_id']]['counts']},
        quality_trainable_parameters=sum(p.numel() for p in head.parameters() if p.requires_grad),
        optimizer_steps_this_run=0, head_weights_unchanged=True, generator_gradients=False,
        new_candidate_search=False, new_backbone_forward=False, DEV_forward_used=False,
        test_setaria_used=False, records=receipts, seconds=time.perf_counter()-began))


def validate_fit(summary, state=False):
    if (summary['arm'] != 'A13' or summary['configuration'] != CONFIG or summary['protocol'] != PROTOCOL or
        summary['step' if state else 'optimizer_steps'] != 4608 or
        summary['revision' if state else 'feature_revision'] != C.I.F.REVISION):
        raise ValueError('Only fixed A13 last head is permitted')
    if not state and (not summary['completed'] or summary['DEV_used_for_selection'] or
                      summary['test_setaria_used'] or summary['final_checkpoint'] !=
                      str((FIT/'step_004608.pt').relative_to(ROOT))):
        raise ValueError('No early/selected/mismatched A13 checkpoint')


def fit():
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed'] or
        ready['configuration'] != CONFIG or not ready['A12_actual_records_verified'] or
        not ready['head_weights_unchanged'] or ready['optimizer_steps_this_run'] != 0):
        raise ValueError('Matching completed zero-update actual-data preflight required')
    draws, index, windows, records, old_rows = load_cohort()
    generator = C.load_generator()
    head = C.new_quality('cuda').train()
    params = [p for p in head.parameters() if p.requires_grad]
    if sum(p.numel() for p in params) != 165635:
        raise ValueError('Same head parameter count required')
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    step = 0; epochs = []
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts = Counter(); total = 0.; started = time.perf_counter()
            for ordinal, d in enumerate(draws):
                ident = d['window_id']; rec = records[ident]; w = windows[ident]
                hidden = inputs(generator, index, w)
                chains = [r['chain'] for r in rec['candidates']]
                labels = [r['label'] for r in rec['candidates']]
                optimizer.zero_grad(set_to_none=True)
                score = head(hidden, chains, w['valid_bases'])
                result = chain_competition_loss(score, chains, labels)
                weighted = result['loss']*d['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(score['logits']).all()) or not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite A13 score/loss')
                weighted.backward()
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                if any(p.grad is not None for p in generator.parameters()):
                    raise ValueError('Frozen generator received gradient')
                optimizer.step(); step += 1
                record = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident,
                    species=w['species'], importance_weight=d['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()), gradient_norm_before_clip=float(norm),
                    positive_competition=float(result['positive_competition'].detach()),
                    negative_null_anchor=float(result['negative_null_anchor'].detach()),
                    conflict_pairs=result['conflict_pairs'], **rec['counts'])
                log.write(json.dumps(record, allow_nan=False)+'\n')
                counts.update(rec['counts']); total += record['weighted_loss']
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(step=step, total=4608)), flush=True)
            checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(arm='A13', configuration=CONFIG, protocol=PROTOCOL, step=step,
                            revision=C.I.F.REVISION, head=head.state_dict(), optimizer=optimizer.state_dict()), checkpoint)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                               mean_weighted_loss=total/1536, seconds=time.perf_counter()-started))
    C.R.write_summary(FIT, dict(arm='A13', completed=True, configuration=CONFIG, protocol=PROTOCOL,
        feature_revision=C.I.F.REVISION, feature_spec=C.R.FEATURE_SPEC, optimizer_steps=step,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(ROOT)), epochs=epochs,
        source_free_pool=str((C.FIT/'free_pool.jsonl').relative_to(ROOT)),
        source_TRAIN_cohort=str(D.OUT.relative_to(ROOT)), generator_trainable_parameters=0,
        quality_trainable_parameters=165635, DEV_used_for_selection=False, test_setaria_used=False,
        new_candidate_search=False, new_backbone_forward=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Matched A6 head initialization/exposure/pool; objective changes both local competition and group mass, not a component-isolated causal comparison.'))


@torch.inference_mode()
def readout():
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    if (FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal fit required')
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    head = C.new_quality('cuda').eval().requires_grad_(False)
    head.load_state_dict(state['head'], strict=True); del state
    generator = C.load_generator()
    draws, index, windows, records, old_rows = load_cohort()
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    rows = []
    with (READOUT/'windows.jsonl').open('x') as log:
        for ordinal, ident in enumerate(records):
            w = windows[ident]; old = records[ident]
            hidden = inputs(generator, index, w)
            score = head(hidden, [r['chain'] for r in old['candidates']], w['valid_bases'])
            vals = {k: score[k].float().cpu().tolist() for k in ('logits', 'additive', 'nonadditive')}
            null = float(score['null']); local = []
            for i, prior in enumerate(old['candidates']):
                r = {**prior, 'gain': vals['logits'][i], 'additive': vals['additive'][i],
                     'nonadditive': vals['nonadditive'][i], 'null': null}
                if not all(math.isfinite(r[k]) for k in D.COMPONENTS):
                    raise ValueError('Nonfinite last-head readout')
                local.append(r)
            for r in local:
                r['max_overlapping_known_nonmatch_gain'] = D.overlapping_max(r['chain'], local) if r['label'] == 1 else None
            log.write(json.dumps(dict(window=w, counts=old['counts'], candidates=local), allow_nan=False)+'\n')
            rows.extend(local)
            if (ordinal+1) % 64 == 0:
                log.flush(); print(json.dumps(dict(scored_TRAIN_windows=ordinal+1, total=1491)), flush=True)
    baseline, candidate = compact(old_rows), compact(rows)
    per_species = {s: {'A6': compact([r for r in old_rows if r['species'] == s]),
                       'A13': compact([r for r in rows if r['species'] == s])} for s in C.R.ALLOW}
    checks = {}
    for s, value in per_species.items():
        a, b = value['A6'], value['A13']
        checks[s+'_free_correct_positive_gain_increased'] = b['free_positive']['positive_gain'] > a['free_positive']['positive_gain']
        checks[s+'_all_correct_positive_gain_increased'] = b['positive']['positive_gain'] > a['positive']['positive_gain']
        checks[s+'_known_nonmatch_positive_gain_not_increased'] = b['known_nonmatch']['positive_gain'] <= a['known_nonmatch']['positive_gain']
    checks['pooled_five_or_more_correct_positive_gain_increased'] = candidate['five_or_more_positive']['positive_gain'] > baseline['five_or_more_positive']['positive_gain']
    checks['pooled_correct_not_above_known_conflict_not_increased'] = candidate['positive']['positive_candidates_not_above_max_known_overlap'] <= baseline['positive']['positive_candidates_not_above_max_known_overlap']
    C.R.write_summary(READOUT, dict(experiment='M28-A13-CHAIN-COMPETITION-TRAIN-READOUT', completed=True,
        configuration=CONFIG, protocol=PROTOCOL, checkpoint=fitted['final_checkpoint'],
        source_generator=C.GENERATOR, source_A6=str(D.OUT.relative_to(ROOT)), windows=1491,
        original_draws=1536, candidate_window_pairs=len(rows), per_species=per_species,
        pooled={'A6': baseline, 'A13': candidate}, predeclared_TRAIN_preparation_gate=checks,
        TRAIN_preparation_supported=all(checks.values()), optimizer_steps_this_run=0,
        candidate_geometry_labels_and_origins_unchanged=True, new_candidate_search=False,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        threshold=0., threshold_scan=False, head_weights_updated=False, checkpoints_written=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='In-sample TRAIN continuation gate only, not final gene F1/background, independent generalization or publication success.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preflight', 'fit', 'readout'), required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    {'preflight': preflight, 'fit': fit, 'readout': readout}[args.mode]()
