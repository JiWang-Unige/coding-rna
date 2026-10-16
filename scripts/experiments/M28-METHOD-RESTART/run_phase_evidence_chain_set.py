#!/usr/bin/env python3
"""A26 frozen candidate-phase evidence fit, final readout and reconciliation."""
import argparse
from collections import Counter
import json
import math
import time
import torch
import preflight_phase_evidence_chain as P

R, ROOT = P.R, P.ROOT
PREFLIGHT, PROTOCOL, CONFIG = P.OUT, P.PROTOCOL, P.CONFIG
FIT = ROOT/'outputs/M28-A26-PHASE-EVIDENCE-FIT'
READOUT = ROOT/'outputs/M28-A26-PHASE-EVIDENCE-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A26-PHASE-EVIDENCE-RESULT-CHECK-REPAIR-R1'
A21_READOUT, A21_CHECK = P.J.READOUT, P.J.CHECK
COMPONENTS = P.J.COMPONENTS+('recurrent_nonadditive', 'phase_evidence_term')
new_head, score_loss = P.new_head, P.score_loss

def paired_a21(records):
    for folder, job in ((A21_READOUT, '13390188'), (A21_CHECK, '13390189')):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
    checked = json.loads((A21_CHECK/'summary.json').read_text())
    assert checked['completed'] and checked['actual_last_head_shapes_finite_tensors_components_selection_profiles_and_gate_verified']
    saved = R.C.R.rows(A21_READOUT/'windows.jsonl')
    assert [r['window']['id'] for r in saved] == list(records)
    result = {}
    for rec in saved:
        ident = rec['window']['id']; original = records[ident]
        assert rec['window'] == original['window'] and rec['counts'] == original['counts']
        assert [{k: r[k] for k in R.IDENTITY} for r in rec['candidates']] == [
            {k: r[k] for k in R.IDENTITY} for r in original['candidates']]
        result[ident] = rec['candidates']
    return result


def validate_fit(value, state=False):
    assert value['arm'] == 'A26' and value['configuration'] == CONFIG and value['protocol'] == PROTOCOL
    assert value['step' if state else 'optimizer_steps'] == 4608
    assert value['revision' if state else 'feature_revision'] == R.C.I.F.REVISION
    if not state:
        assert value['completed'] and not value['DEV_used_for_selection'] and not value['test_setaria_used']
        assert value['final_checkpoint'] == str((FIT/'step_004608.pt').relative_to(ROOT))


def gate(selected):
    checks = {}
    for species, values in selected.items():
        old, new = values['A6'], values['A26']
        checks[species+'_natural_selected_correct_increased'] = new['selected_correct'] > old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = new['selected_known_nonmatch'] <= old['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_increased'] = (
        sum(v['A26']['selected_five_or_more_correct'] for v in selected.values()) >
        sum(v['A6']['selected_five_or_more_correct'] for v in selected.values()))
    return checks


@torch.inference_mode()
def readout(hidden_by_id, evidence_by_id, windows, records):
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    assert (FIT/'STATUS').read_text().strip() == 'COMPLETED'
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    head = new_head('cuda').eval().requires_grad_(False); head.load_state_dict(state['head'], strict=True)
    del state
    paired = paired_a21(records); began = time.perf_counter()
    selected = {s: {a: Counter() for a in ('A6', 'A21', 'A26')} for s in R.C.R.ALLOW}
    all_rows = {a: [] for a in ('A6', 'A21', 'A26')}
    with (READOUT/'windows.jsonl').open('x') as log:
        for ordinal, (ident, rec) in enumerate(records.items()):
            scores, _ = score_loss(head, hidden_by_id[ident], evidence_by_id[ident], windows[ident], rec)
            values = {k: scores[k].float().cpu().tolist() for k in
                      ('logits', 'additive', 'nonadditive', 'exon_additive', 'junction_additive', 'recurrent_nonadditive', 'phase_evidence_term')}
            local = []
            for i, prior in enumerate(rec['candidates']):
                row = {**prior, 'gain': values['logits'][i], 'null': float(scores['null']),
                    **{k: values[k][i] for k in ('additive', 'nonadditive', 'exon_additive', 'junction_additive', 'recurrent_nonadditive', 'phase_evidence_term')}}
                assert abs(row['gain']-(row['additive']+row['nonadditive']-row['null'])) <= 1e-4
                assert abs(row['additive']-row['exon_additive']-row['junction_additive']) <= 1e-4
                local.append(row)
            for row in local:
                row['max_overlapping_known_nonmatch_gain'] = R.D.overlapping_max(row['chain'], local) if row['label'] == 1 else None
            candidates = dict(A6=rec['candidates'], A21=paired[ident], A26=local)
            comparison = {a: R.S.natural_selection(rows) for a, rows in candidates.items()}
            for arm, counts in comparison.items():
                selected[windows[ident]['species']][arm].update(counts); all_rows[arm].extend(candidates[arm])
            log.write(json.dumps(dict(window=rec['window'], counts=rec['counts'], candidates=local,
                                     natural_selection=comparison), allow_nan=False)+'\n')
            if (ordinal+1) % 64 == 0:
                log.flush(); print(json.dumps(dict(phase='TRAIN_readout', windows=ordinal+1, total=1491)), flush=True)
    checks = gate(selected)
    R.C.R.write_summary(READOUT, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        checkpoint=fitted['final_checkpoint'], windows=len(records), candidate_window_pairs=len(all_rows['A26']),
        natural_selection_per_species={s: {a: dict(c) for a, c in v.items()} for s, v in selected.items()},
        score_profiles={a: R.P.compact(rows) for a, rows in all_rows.items()},
        predeclared_TRAIN_preparation_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        optimizer_steps_this_run=0, new_candidate_search=False, new_backbone_forward=False,
        DEV_forward_used=False, test_setaria_used=False, threshold=0., threshold_scan=False,
        natural_candidate_geometry_labels_and_origins_unchanged=True, reference_injected_candidates=0,
        same_fit_immutable_hidden_used=True, same_fit_frozen_phase_evidence_used=True,
        reference_used_in_evidence=False, seconds=time.perf_counter()-began,
        limitation='TRAIN candidate/window pairs, increased head capacity; not causal isolation, gene F1 or generalization.'))


def fit_readout():
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    assert (PREFLIGHT/'STATUS').read_text().strip() == 'COMPLETED'
    assert ready['completed'] and ready['configuration'] == CONFIG and ready['protocol'] == PROTOCOL
    assert ready['head_weights_unchanged'] and ready['optimizer_steps_this_run'] == 0
    assert ready['exact_A21_common_fresh_parameters'] and ready['zero_prior_actual_scores_and_loss_verified']
    assert ready['generator_frozen'] and not ready['reference_used_in_evidence']
    assert len(ready['gradient_witnesses']) == 8
    assert all(w['unknown_direct_gain_gradient_zero'] for w in ready['gradient_witnesses'])
    began = time.perf_counter()
    draws, index, windows, records, _ = R.population()
    assert ready['natural_counts'] == R.totals(records, records)
    assert ready['replayed_draw_counts'] == R.totals(records, [d['window_id'] for d in draws])
    generator = R.C.load_generator(); head = new_head('cuda').train()
    hidden_bytes = sum(w['width']*128*4 for w in windows.values())
    assert hidden_bytes == 18761121792
    timing = json.loads((ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-TIMING/summary.json').read_text())
    assert timing['completed'] and timing['optimizer_steps_this_run'] == 0
    evidence_bytes = sum(len(r['candidates'])*28*4 for r in records.values())
    assert evidence_bytes == 1827952
    # Measured old pipeline peak plus 64 MiB covers new head/transients.
    assert hidden_bytes+evidence_bytes+timing['peak_GPU_allocated_bytes']+64*1024**2 < torch.cuda.mem_get_info()[0]
    hidden_by_id, evidence_by_id = {}, {}; cache_began = time.perf_counter()
    torch.cuda.reset_peak_memory_stats()
    with (FIT/'phase_evidence.jsonl').open('x') as evidence_log:
        for ordinal, (ident, window) in enumerate(windows.items()):
            hidden = R.P.inputs(generator, index, window).detach()
            assert hidden.dtype == torch.float32 and list(hidden.shape) == [window['width'], 128]
            assert not hidden.requires_grad and bool(torch.isfinite(hidden).all())
            evidence = P.evidence_inputs(generator, hidden, window, records[ident])
            assert evidence.dtype == torch.float32 and list(evidence.shape) == [len(records[ident]['candidates']), 28]
            hidden_by_id[ident], evidence_by_id[ident] = hidden, evidence
            evidence_log.write(json.dumps(dict(window_id=ident,
                chains=[r['chain'] for r in records[ident]['candidates']],
                vectors=evidence.cpu().tolist()), allow_nan=False)+'\n')
            if (ordinal+1) % 64 == 0:
                evidence_log.flush()
                print(json.dumps(dict(phase='immutable_hidden_and_state_evidence', windows=ordinal+1, total=1491)), flush=True)
    cache_seconds = time.perf_counter()-cache_began
    assert len(hidden_by_id) == len(evidence_by_id) == 1491
    assert sum(h.numel()*h.element_size() for h in hidden_by_id.values()) == hidden_bytes
    assert sum(e.numel()*e.element_size() for e in evidence_by_id.values()) == evidence_bytes
    witnesses = P.zero_update_witnesses(head, generator, index, windows, records, hidden_by_id, evidence_by_id)
    with (FIT/'cache_receipt.json').open('x') as handle:
        json.dump(dict(completed=True, windows=1491, hidden_bytes=hidden_bytes, precision='FP32',
            evidence_bytes=evidence_bytes, evidence_dimensions=28, candidate_window_pairs=16321,
            cache_seconds=cache_seconds, witnesses=witnesses, optimizer_steps_this_check=0,
            head_weights_unchanged=True, generator_frozen=True, reference_used_in_evidence=False,
            materialization_peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated()), handle, indent=2, allow_nan=False)
        handle.write('\n')
    params = list(head.parameters())
    optimizer = torch.optim.AdamW(params, lr=CONFIG['lr'], weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']), eps=CONFIG['eps'])
    step = 0; epochs = []; evidence_gradient_steps = 0
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts = Counter(); loss_sum = 0.; started = time.perf_counter()
            for ordinal, draw in enumerate(draws):
                ident = draw['window_id']; optimizer.zero_grad(set_to_none=True)
                _, result = score_loss(head, hidden_by_id[ident], evidence_by_id[ident], windows[ident], records[ident])
                weighted = result['loss']*draw['importance_weight_uniform_within_stratum']
                assert bool(torch.isfinite(weighted)); weighted.backward()
                evidence_grad = head.path_evidence.weight.grad
                assert evidence_grad is not None and bool(torch.isfinite(evidence_grad).all())
                evidence_norm = float(evidence_grad.double().norm())
                evidence_gradient_steps += int(evidence_norm > 0.)
                norm = torch.nn.utils.clip_grad_norm_(params, CONFIG['clip_norm'], error_if_nonfinite=True)
                assert not any(p.grad is not None for p in generator.parameters())
                optimizer.step(); step += 1
                row = dict(step=step, epoch=epoch+1, draw_ordinal=ordinal, window_id=ident,
                    species=windows[ident]['species'], importance_weight=draw['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()), log_partition=float(result['log_partition'].detach()),
                    gold_gain=float(result['gold_gain'].detach()), gradient_norm_before_clip=float(norm), phase_evidence_gradient_norm_before_clip=evidence_norm,
                    **records[ident]['counts'])
                log.write(json.dumps(row, allow_nan=False)+'\n'); counts.update(records[ident]['counts']); loss_sum += row['weighted_loss']
                if step % 64 == 0:
                    log.flush(); print(json.dumps(dict(phase='fit', step=step, total=4608)), flush=True)
                if step % 256 == 0:
                    with (FIT/('step_'+str(step).zfill(6)+'.pt')).open('xb') as handle:
                        torch.save(dict(arm='A26', configuration=CONFIG, protocol=PROTOCOL, step=step,
                            revision=R.C.I.F.REVISION, head=head.state_dict(), optimizer=optimizer.state_dict()), handle)
            epochs.append(dict(epoch=epoch+1, steps=1536, counts=dict(counts),
                               mean_weighted_loss=loss_sum/1536, seconds=time.perf_counter()-started))
    assert step == 4608 and evidence_gradient_steps > 0
    evidence_weight_norm = float(head.path_evidence.weight.detach().double().norm())
    assert evidence_weight_norm > 0. and math.isfinite(evidence_weight_norm)
    R.C.R.write_summary(FIT, dict(arm='A26', completed=True, configuration=CONFIG, protocol=PROTOCOL,
        feature_revision=R.C.I.F.REVISION, optimizer_steps=step,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(ROOT)), epochs=epochs,
        quality_trainable_parameters=CONFIG['quality_trainable_parameters'], generator_trainable_parameters=0,
        DEV_used_for_selection=False, test_setaria_used=False, new_candidate_search=False,
        new_backbone_forward=False, checkpoint_interval=256,
        immutable_hidden_cache_receipt=str((FIT/'cache_receipt.json').relative_to(ROOT)),
        phase_evidence_nonzero_gradient_steps=evidence_gradient_steps,
        phase_evidence_final_weight_norm=evidence_weight_norm, reference_used_in_evidence=False,
        wall_seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated()))
    (FIT/'STATUS').write_text('COMPLETED\n')
    readout(hidden_by_id, evidence_by_id, windows, records)


def check(fit_job):
    for folder in (FIT, READOUT):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == fit_job
    fitted = json.loads((FIT/'summary.json').read_text()); validate_fit(fitted)
    report = json.loads((READOUT/'summary.json').read_text())
    assert report['completed'] and report['configuration'] == CONFIG and report['protocol'] == PROTOCOL
    assert report['checkpoint'] == fitted['final_checkpoint'] and report['same_fit_immutable_hidden_used']
    assert report['same_fit_frozen_phase_evidence_used'] and not report['reference_used_in_evidence']
    assert not fitted['reference_used_in_evidence']
    assert report['optimizer_steps_this_run'] == 0 and not report['DEV_forward_used'] and not report['test_setaria_used']
    assert not report['new_candidate_search'] and not report['new_backbone_forward']
    assert report['threshold'] == 0. and not report['threshold_scan']
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True); validate_fit(state, True)
    head = new_head('cpu'); head.load_state_dict(state['head'], strict=True)
    assert all(bool(torch.isfinite(v).all()) for v in state['head'].values())
    for step in range(256, 4609, 256):
        assert (FIT/('step_'+str(step).zfill(6)+'.pt')).stat().st_size > 0
    cache = json.loads((FIT/'cache_receipt.json').read_text())
    assert cache['completed'] and cache['windows'] == 1491 and cache['hidden_bytes'] == 18761121792
    assert cache['optimizer_steps_this_check'] == 0 and cache['head_weights_unchanged'] and cache['generator_frozen']
    assert len(cache['witnesses']) == 8 and all(w['live_cached_equivalence'] and w['unknown_direct_gain_gradient_zero'] for w in cache['witnesses'])
    assert cache['evidence_bytes'] == 1827952 and cache['evidence_dimensions'] == 28
    assert cache['candidate_window_pairs'] == 16321 and not cache['reference_used_in_evidence']
    assert all(w['reference_used_in_evidence'] is False and w['evidence_finite'] for w in cache['witnesses'])
    draws, _, windows, originals, _ = R.population(); paired = paired_a21(originals)
    evidence_saved = R.C.R.rows(FIT/'phase_evidence.jsonl')
    # Materialization follows windows, while natural readout follows records.
    # Candidate vectors below remain checked and joined by the original ID.
    assert [r['window_id'] for r in evidence_saved] == list(windows)
    evidence_by_id = {}
    for rec in evidence_saved:
        ident = rec['window_id']
        assert rec['chains'] == [r['chain'] for r in originals[ident]['candidates']]
        assert len(rec['vectors']) == len(rec['chains'])
        assert all(len(v) == 28 and all(math.isfinite(x) for x in v) for v in rec['vectors'])
        evidence_by_id[ident] = rec['vectors']
    projection = state['head']['path_evidence.weight'][0].tolist()
    assert len(projection) == 28
    weight_norm = float(state['head']['path_evidence.weight'].double().norm())
    assert weight_norm > 0. and math.isclose(weight_norm, fitted['phase_evidence_final_weight_norm'], rel_tol=1e-12, abs_tol=1e-12)
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    assert ready['natural_counts'] == R.totals(originals, originals)
    assert ready['replayed_draw_counts'] == R.totals(originals, [d['window_id'] for d in draws])
    training = R.C.R.rows(FIT/'training.jsonl'); assert len(training) == 4608
    counts = [Counter() for _ in range(3)]; losses = [0.] * 3; evidence_gradient_steps = 0
    for i, row in enumerate(training):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; ident = draw['window_id']
        assert (row['step'], row['epoch'], row['draw_ordinal'], row['window_id'], row['species']) == (
            i+1, epoch+1, ordinal, ident, windows[ident]['species'])
        assert row['importance_weight'] == draw['importance_weight_uniform_within_stratum']
        assert all(row[k] == v for k, v in originals[ident]['counts'].items())
        assert all(math.isfinite(row[k]) for k in ('weighted_loss', 'log_partition', 'gold_gain', 'gradient_norm_before_clip', 'phase_evidence_gradient_norm_before_clip'))
        assert row['phase_evidence_gradient_norm_before_clip'] >= 0.
        evidence_gradient_steps += int(row['phase_evidence_gradient_norm_before_clip'] > 0.)
        counts[epoch].update(originals[ident]['counts']); losses[epoch] += row['weighted_loss']
    assert len(fitted['epochs']) == 3
    assert evidence_gradient_steps == fitted['phase_evidence_nonzero_gradient_steps'] > 0
    for i, count in enumerate(counts):
        assert dict(count) == fitted['epochs'][i]['counts'] == ready['replayed_draw_counts']
        assert fitted['epochs'][i]['steps'] == 1536 and losses[i]/1536 == fitted['epochs'][i]['mean_weighted_loss']
    saved = R.C.R.rows(READOUT/'windows.jsonl'); assert [r['window']['id'] for r in saved] == list(originals)
    selected = {s: {a: Counter() for a in ('A6', 'A21', 'A26')} for s in R.C.R.ALLOW}
    rows = {a: [] for a in ('A6', 'A21', 'A26')}
    excluded = set(COMPONENTS)|{'max_overlapping_known_nonmatch_gain'}
    for rec in saved:
        ident = rec['window']['id']; original = originals[ident]
        assert rec['window'] == original['window'] and rec['counts'] == original['counts']
        local = rec['candidates']; assert len(local) == len(original['candidates'])
        for a, b, evidence in zip(local, original['candidates'], evidence_by_id[ident]):
            assert {k: v for k, v in a.items() if k not in excluded} == {k: v for k, v in b.items() if k not in excluded}
            assert a['origin'] == 'free' and all(math.isfinite(a[k]) for k in COMPONENTS)
            assert abs(a['gain']-(a['additive']+a['nonadditive']-a['null'])) <= 1e-4
            assert abs(a['additive']-a['exon_additive']-a['junction_additive']) <= 1e-4
            assert abs(a['nonadditive']-a['recurrent_nonadditive']-a['phase_evidence_term']) <= 1e-4
            term = math.fsum(w*x for w, x in zip(projection, evidence))
            assert math.isclose(term, a['phase_evidence_term'], rel_tol=1e-5, abs_tol=1e-4)
            assert a['max_overlapping_known_nonmatch_gain'] == (R.D.overlapping_max(a['chain'], local) if a['label'] == 1 else None)
        candidates = dict(A6=original['candidates'], A21=paired[ident], A26=local)
        comparison = {a: R.S.natural_selection(v) for a, v in candidates.items()}
        assert comparison == rec['natural_selection']
        for arm, values in comparison.items():
            selected[rec['window']['species']][arm].update(values); rows[arm].extend(candidates[arm])
    actual = {s: {a: dict(c) for a, c in v.items()} for s, v in selected.items()}
    assert actual == report['natural_selection_per_species']
    assert {a: R.P.compact(v) for a, v in rows.items()} == report['score_profiles']
    assert len(saved) == report['windows'] == 1491 and len(rows['A26']) == report['candidate_window_pairs'] == 16321
    assert report['reference_injected_candidates'] == 0 and report['natural_candidate_geometry_labels_and_origins_unchanged']
    checks = gate(selected)
    assert checks == report['predeclared_TRAIN_preparation_gate'] and all(checks.values()) == report['TRAIN_preparation_supported']
    R.C.R.write_summary(CHECK, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL, source_fit_readout_job=fit_job,
        actual_4608_steps_draw_weights_and_original_natural_domain_verified=True,
        actual_last_head_shapes_finite_tensors_components_selection_profiles_and_gate_verified=True,
        cache_zero_update_equivalence_verified=True, original_unknown_semantics_unchanged=True,
        actual_frozen_evidence_geometry_order_and_projection_components_verified=True,
        actual_evidence_materialization_window_order_verified=True,
        actual_new_projection_gradient_and_final_weight_verified=True, reference_used_in_evidence=False,
        natural_selection_per_species=actual, original_frozen_TRAIN_gate=checks,
        TRAIN_preparation_supported=all(checks.values()), windows=1491, candidate_window_pairs=16321,
        optimizer_steps_this_run=0, new_head_or_backbone_inference=False,
        DEV_forward_used=False, test_setaria_used=False,
        limitation='TRAIN mechanism experiment; no independent gene F1, causal isolation or publication claim.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('fit-readout', 'check'), required=True)
    parser.add_argument('--fit-job'); args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'fit-readout':
        fit_readout()
    else:
        assert args.fit_job; check(args.fit_job)
