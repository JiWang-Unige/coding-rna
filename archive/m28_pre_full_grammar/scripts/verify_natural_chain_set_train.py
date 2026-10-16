#!/usr/bin/env python3
"""Read-only reconciliation of the closed A15 TRAIN fit and matched readout."""
import argparse
from collections import Counter
import json
import math
import torch
import run_natural_chain_set as R

OUT = R.ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-RESULT-CHECK'


def main(fit_job, readout_job):
    for folder, job in ((R.FIT, fit_job), (R.READOUT, readout_job)):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
    fitted = json.loads((R.FIT/'summary.json').read_text()); R.validate_fit(fitted)
    report = json.loads((R.READOUT/'summary.json').read_text())
    assert report['completed'] and report['configuration'] == R.CONFIG and report['protocol'] == R.PROTOCOL
    assert report['checkpoint'] == fitted['final_checkpoint']
    assert report['optimizer_steps_this_run'] == 0 and not report['DEV_forward_used'] and not report['test_setaria_used']
    assert not report['new_candidate_search'] and not report['new_backbone_forward']
    assert report['threshold'] == 0. and not report['threshold_scan']
    state = torch.load(R.ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    R.validate_fit(state, state=True)
    assert all(bool(torch.isfinite(v).all()) for v in state['head'].values())
    for step in (1536, 3072, 4608):
        assert (R.FIT/('step_'+str(step).zfill(6)+'.pt')).stat().st_size > 0
    draws, _, windows, originals, source = R.population()
    paired = R.a14_comparator(originals)
    ready = json.loads((R.PREFLIGHT/'summary.json').read_text())
    assert R.totals(source, source) == ready['source_counts']
    assert R.totals(originals, originals) == ready['natural_counts']
    assert R.totals(originals, [d['window_id'] for d in draws]) == ready['replayed_draw_counts']
    training = R.C.R.rows(R.FIT/'training.jsonl'); assert len(training) == 4608
    epochs = [Counter() for _ in range(3)]; losses = [0.] * 3
    for i, row in enumerate(training):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; ident = draw['window_id']
        assert (row['step'], row['epoch'], row['draw_ordinal'], row['window_id'], row['species']) == (
            i+1, epoch+1, ordinal, ident, windows[ident]['species'])
        assert row['importance_weight'] == draw['importance_weight_uniform_within_stratum']
        assert all(row[k] == v for k, v in originals[ident]['counts'].items())
        assert all(math.isfinite(row[k]) for k in ('weighted_loss', 'log_partition', 'gold_gain', 'gradient_norm_before_clip'))
        epochs[epoch].update(originals[ident]['counts']); losses[epoch] += row['weighted_loss']
    assert len(fitted['epochs']) == 3
    for i, counts in enumerate(epochs):
        assert dict(counts) == fitted['epochs'][i]['counts'] == ready['replayed_draw_counts']
        assert fitted['epochs'][i]['steps'] == 1536 and losses[i]/1536 == fitted['epochs'][i]['mean_weighted_loss']
    records = R.C.R.rows(R.READOUT/'windows.jsonl')
    assert [r['window']['id'] for r in records] == list(originals) and len(records) == 1491
    selected = {s: {a: Counter() for a in ('A6', 'A14', 'A15')} for s in R.C.R.ALLOW}
    rows = {a: [] for a in ('A6', 'A14', 'A15')}
    excluded = set(R.D.COMPONENTS)|{'max_overlapping_known_nonmatch_gain'}
    for rec in records:
        ident = rec['window']['id']; old = originals[ident]; species = rec['window']['species']
        assert rec['window'] == old['window'] and rec['counts'] == old['counts']
        local = rec['candidates']; assert len(local) == len(old['candidates'])
        for a, b in zip(local, old['candidates']):
            assert {k: v for k, v in a.items() if k not in excluded} == {k: v for k, v in b.items() if k not in excluded}
            assert a['origin'] == 'free' and all(math.isfinite(a[k]) for k in R.D.COMPONENTS)
            assert abs(a['gain']-(a['additive']+a['nonadditive']-a['null'])) <= 1e-4
            assert a['max_overlapping_known_nonmatch_gain'] == (R.D.overlapping_max(a['chain'], local) if a['label'] == 1 else None)
        comparison = {arm: R.S.natural_selection(values) for arm, values in
                      (('A6', old['candidates']), ('A14', paired[ident]), ('A15', local))}
        assert comparison == rec['natural_selection']
        for arm, values in comparison.items():
            selected[species][arm].update(values)
        rows['A6'].extend(old['candidates']); rows['A14'].extend(paired[ident]); rows['A15'].extend(local)
    actual = {s: {a: dict(c) for a, c in v.items()} for s, v in selected.items()}
    assert actual == report['natural_selection_per_species']
    assert {a: R.P.compact(v) for a, v in rows.items()} == report['score_profiles']
    assert len(rows['A15']) == report['candidate_window_pairs'] == 16321
    assert report['reference_injected_candidates'] == 0 and report['natural_candidate_geometry_labels_and_origins_unchanged']
    checks = {}
    for s, values in actual.items():
        a, b = values['A6'], values['A15']
        checks[s+'_natural_selected_correct_increased'] = b['selected_correct'] > a['selected_correct']
        checks[s+'_natural_selected_known_nonmatch_not_increased'] = b['selected_known_nonmatch'] <= a['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_increased'] = (
        sum(v['A15']['selected_five_or_more_correct'] for v in actual.values()) >
        sum(v['A6']['selected_five_or_more_correct'] for v in actual.values()))
    assert checks == report['predeclared_TRAIN_preparation_gate'] and all(checks.values()) == report['TRAIN_preparation_supported']
    result = dict(completed=True, source_jobs=dict(fit=fit_job, readout=readout_job),
        actual_natural_domain_and_original_labels_verified=True,
        actual_4608_steps_three_pass_draw_and_weight_identity_verified=True,
        actual_paired_selection_profiles_and_frozen_A6_gate_recomputed=True,
        all_saved_head_tensors_finite=True, windows=1491, candidate_window_pairs=16321,
        reference_injected_candidates=0, source_counts=ready['source_counts'], natural_counts=ready['natural_counts'],
        replayed_draw_counts=ready['replayed_draw_counts'], windows_losing_all_gold=ready['windows_losing_all_gold'],
        natural_selection_per_species=actual, original_frozen_TRAIN_gate=checks,
        TRAIN_preparation_supported=all(checks.values()), optimizer_steps_this_run=0,
        new_head_or_backbone_inference=False, DEV_forward_used=False, test_setaria_used=False,
        limitation='Matched TRAIN draw/optimizer exposure, reduced positive exposure; no independent gene F1 or causal proof.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-job', required=True); parser.add_argument('--readout-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.fit_job, args.readout_job)
