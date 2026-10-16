#!/usr/bin/env python3
"""Reconcile terminal A14 TRAIN output and diagnose injected-gold masking."""
from collections import Counter
import json
import math
import torch
import run_chain_set as R
from src.m28.chain_set import interval_log_partition
from src.m28.core import select_nonoverlapping

OUT = R.ROOT/'outputs/M28-A14-CHAIN-SET-RESULT-CHECK'


def selected_indices(rows, natural):
    indices = [i for i, r in enumerate(rows) if not natural or r['origin'] == 'free']
    chosen = select_nonoverlapping([rows[i]['chain'] for i in indices], [rows[i]['gain'] for i in indices])
    return {indices[i] for i in chosen}


def marginals(rows, indices):
    gain = torch.tensor([rows[i]['gain'] for i in indices], dtype=torch.float64, requires_grad=True)
    logz = interval_log_partition(gain, [rows[i]['chain'] for i in indices])
    values = torch.autograd.grad(logz, gain)[0].tolist()
    return dict(zip(indices, values))


def main():
    for folder, job in ((R.FIT, '13366600'), (R.READOUT, '13366602')):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
    fitted = json.loads((R.FIT/'summary.json').read_text()); R.validate_fit(fitted)
    report = json.loads((R.READOUT/'summary.json').read_text())
    assert report['completed'] and report['configuration'] == R.CONFIG and report['protocol'] == R.PROTOCOL
    assert report['checkpoint'] == fitted['final_checkpoint'] and not report['DEV_forward_used'] and not report['test_setaria_used']
    state = torch.load(R.ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    R.validate_fit(state, state=True)
    assert all(bool(torch.isfinite(v).all()) for v in state['head'].values())
    for step in (1536, 3072, 4608):
        assert (R.FIT/('step_'+str(step).zfill(6)+'.pt')).stat().st_size > 0
    draws, _, windows, originals, old_rows = R.population()
    training = R.C.R.rows(R.FIT/'training.jsonl')
    assert len(training) == 4608
    epoch_counts = [Counter() for _ in range(3)]
    for i, row in enumerate(training):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; ident = draw['window_id']
        assert (row['step'], row['epoch'], row['draw_ordinal'], row['window_id'], row['species']) == (
            i+1, epoch+1, ordinal, ident, windows[ident]['species'])
        assert row['importance_weight'] == draw['importance_weight_uniform_within_stratum']
        assert all(row[k] == v for k, v in originals[ident]['counts'].items())
        assert all(math.isfinite(row[k]) for k in ('weighted_loss', 'log_partition', 'gold_gain', 'gradient_norm_before_clip'))
        epoch_counts[epoch].update(originals[ident]['counts'])
    assert all(dict(counts) == fitted['epochs'][i]['counts'] for i, counts in enumerate(epoch_counts))
    records = R.C.R.rows(R.READOUT/'windows.jsonl')
    assert [r['window']['id'] for r in records] == list(originals) and len(records) == 1491
    selected = {s: {'A6': Counter(), 'A14': Counter()} for s in R.C.R.ALLOW}
    masking = {s: Counter() for s in R.C.R.ALLOW}
    probability = {s: [] for s in R.C.R.ALLOW}
    new_rows = []
    excluded = set(R.D.COMPONENTS)|{'max_overlapping_known_nonmatch_gain'}
    for rec in records:
        ident = rec['window']['id']; old = originals[ident]; species = rec['window']['species']
        assert rec['window'] == old['window'] and rec['counts'] == old['counts']
        local, prior = rec['candidates'], old['candidates']
        assert len(local) == len(prior)
        for a, b in zip(local, prior):
            assert {k: v for k, v in a.items() if k not in excluded} == {k: v for k, v in b.items() if k not in excluded}
            assert all(math.isfinite(a[k]) for k in R.D.COMPONENTS)
            assert abs(a['gain']-(a['additive']+a['nonadditive']-a['null'])) <= 1e-4
            assert a['max_overlapping_known_nonmatch_gain'] == (R.D.overlapping_max(a['chain'], local) if a['label'] == 1 else None)
        a, b = R.natural_selection(prior), R.natural_selection(local)
        assert rec['natural_selection'] == {'A6': a, 'A14': b}
        selected[species]['A6'].update(a); selected[species]['A14'].update(b)
        # Same saved last-head scores; adding TRAIN-injected chains is a
        # diagnostic counterfactual, never a real inference or accuracy result.
        free_selected, augmented_selected = selected_indices(local, True), selected_indices(local, False)
        injected = [i for i in augmented_selected if local[i]['origin'] == 'reference_injected']
        false_selected = [i for i in free_selected if local[i]['label'] == 0]
        masked = [i for i in false_selected if i not in augmented_selected and any(
            local[i]['chain'][0][0] < local[j]['chain'][-1][1] and
            local[i]['chain'][-1][1] > local[j]['chain'][0][0] for j in injected)]
        masking[species].update(windows=1, natural_selected_known_nonmatch=len(false_selected),
            natural_known_nonmatch_masked_by_selected_injected_gold=len(masked), windows_with_masking=int(bool(masked)),
            augmented_selected_injected_correct=len(injected),
            augmented_selected_known_nonmatch=sum(local[i]['label'] == 0 for i in augmented_selected))
        if masked:
            known = [i for i, r in enumerate(local) if r['label'] >= 0]
            natural_known = [i for i in known if local[i]['origin'] == 'free']
            full, natural = marginals(local, known), marginals(local, natural_known)
            assert all(-1e-12 <= p <= 1+1e-12 for p in list(full.values())+list(natural.values()))
            probability[species].extend((full[i], natural[i]) for i in masked)
        new_rows.extend(local)
    actual = {s: {name: dict(values) for name, values in pair.items()} for s, pair in selected.items()}
    assert actual == report['natural_selection_per_species']
    assert R.P.compact(old_rows) == report['score_profiles']['A6']
    assert R.P.compact(new_rows) == report['score_profiles']['A14']
    assert len(new_rows) == report['candidate_window_pairs'] == 16857
    checks = {}
    for s, values in actual.items():
        a, b = values['A6'], values['A14']
        checks[s+'_natural_selected_correct_increased'] = b['selected_correct'] > a['selected_correct']
        checks[s+'_natural_selected_known_nonmatch_not_increased'] = b['selected_known_nonmatch'] <= a['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_increased'] = (
        sum(v['A14']['selected_five_or_more_correct'] for v in actual.values()) >
        sum(v['A6']['selected_five_or_more_correct'] for v in actual.values()))
    assert checks == report['predeclared_TRAIN_preparation_gate'] and all(checks.values()) == report['TRAIN_preparation_supported']
    diagnostic = {}
    for s, values in masking.items():
        pairs = probability[s]
        diagnostic[s] = dict(**values, masked_pairs=len(pairs),
            mean_augmented_known_partition_marginal_on_masked_nonmatches=sum(a for a, _ in pairs)/len(pairs) if pairs else None,
            mean_natural_known_partition_marginal_on_masked_nonmatches=sum(b for _, b in pairs)/len(pairs) if pairs else None)
    result = dict(completed=True, actual_original_window_candidate_source_label_identity_verified=True,
        actual_4608_steps_three_pass_draw_and_weight_identity_verified=True,
        actual_natural_selection_score_profiles_and_frozen_gate_recomputed=True,
        all_saved_head_tensors_finite=True, windows=1491, candidate_window_pairs=16857,
        natural_selection_per_species=actual, original_frozen_TRAIN_gate=checks,
        TRAIN_preparation_supported=all(checks.values()), injected_gold_masking_diagnostic=diagnostic,
        optimizer_steps_this_run=0, new_head_or_backbone_inference=False, DEV_forward_used=False, test_setaria_used=False,
        limitation='Exploratory same-TRAIN counterfactual at unchanged saved scores. Injected gold is unavailable at inference; no new real annotation result, model fit, DEV, causal repair guarantee or biological false-positive claim.')
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
