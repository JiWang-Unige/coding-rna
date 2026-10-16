#!/usr/bin/env python3
"""A20 CPU preparation using closed A15 gains, never model/head training."""
from collections import Counter
import json
import math
import time
import torch
import recover_natural_chain_set as E
from src.m28.chain_set import chain_set_loss
from src.m28.labels import oriented_parts
from src.m28.partial_chain_set import partial_chain_set_loss, span_overlap

R = E.R
OUT = R.ROOT/'outputs/M28-A20-PARTIAL-CHAIN-SET-PREPARATION'
PROTOCOL = 'reports/M28-METHOD-RESTART/a20_partial_chain_set_preparation.md'


def aggregate(rows):
    fields = ('candidates', 'known_positive', 'known_nonmatch', 'unknown',
              'unknown_conflicting_with_gold', 'exact_protected_unknown_conflicting_with_gold')
    total = Counter()
    for row in rows:
        total.update({k: row[k] for k in fields})
    return dict(windows_or_draws=len(rows), **dict(total),
                windows_or_draws_with_unknown=sum(r['unknown'] > 0 for r in rows),
                windows_or_draws_with_unknown_gold_conflict=sum(r['unknown_conflicting_with_gold'] > 0 for r in rows),
                windows_or_draws_with_exact_protected_gold_conflict=sum(
                    r['exact_protected_unknown_conflicting_with_gold'] > 0 for r in rows))


def main():
    began = time.perf_counter()
    checked = json.loads((E.V.OUT/'summary.json').read_text())
    assert (E.V.OUT/'STATUS').read_text().strip() == 'COMPLETED'
    assert (E.V.OUT/'JOBID').read_text().strip() == '13387977'
    assert checked['completed'] and checked['source_jobs'] == dict(fit='13387908', readout='13387947')
    assert checked['actual_natural_domain_and_original_labels_verified']
    assert checked['actual_4608_steps_three_pass_draw_and_weight_identity_verified']
    assert checked['actual_paired_selection_profiles_and_frozen_A6_gate_recomputed']
    fitted = json.loads((R.FIT/'summary.json').read_text()); R.validate_fit(fitted)
    draws, _, windows, originals, _ = R.population()
    metadata_draws, _, metadata_windows, catalogs = R.C.R.Q.load_train_inputs()
    assert draws == metadata_draws and windows == metadata_windows
    records = R.C.R.rows(R.READOUT/'windows.jsonl')
    assert [r['window']['id'] for r in records] == list(originals) and len(records) == 1491
    identity = ('window_id', 'species', 'chain', 'origin', 'label', 'exon_class')
    results = []; gradients = []; null_gradients = []; witnesses = []
    for rec in records:
        w = rec['window']; ident = w['id']; local = rec['candidates']
        assert w == originals[ident]['window'] and rec['counts'] == originals[ident]['counts']
        assert [{k: r[k] for k in identity} for r in local] == [
            {k: r[k] for k in identity} for r in originals[ident]['candidates']]
        chains = [r['chain'] for r in local]; labels = [r['label'] for r in local]
        assert len({tuple(map(tuple, c)) for c in chains}) == len(chains)
        assert all(r['origin'] == 'free' and math.isfinite(r['gain']) for r in local)
        gain = torch.tensor([r['gain'] for r in local], dtype=torch.float64, requires_grad=True)
        score = dict(logits=gain, null=torch.tensor(0., dtype=torch.float64, requires_grad=True))
        baseline = chain_set_loss(score, chains, labels)
        partial = partial_chain_set_loss(score, chains, labels)
        old_gradient = torch.autograd.grad(baseline['loss'], gain, retain_graph=True)[0]
        grad, null_grad = torch.autograd.grad(partial['loss'], (gain, score['null']), allow_unused=True)
        assert bool(torch.isfinite(grad).all())
        assert all(bool(torch.isfinite(partial[k])) for k in ('loss', 'log_partition', 'observed_log_partition'))
        assert float(partial['loss']) >= -1e-10
        unknown = [i for i, y in enumerate(labels) if y == -1]
        positive = [i for i, y in enumerate(labels) if y == 1]
        conflict = [i for i in unknown if any(span_overlap(chains[i], chains[j]) for j in positive)]
        assert all(float(old_gradient[i]) == 0. for i in unknown)
        if not unknown:
            torch.testing.assert_close(partial['loss'], baseline['loss'], rtol=1e-12, atol=1e-12)
            torch.testing.assert_close(grad, old_gradient, rtol=1e-12, atol=1e-12)
        protected = {tuple(oriented_parts(catalogs[w['species']][tx], w))
                     for tx in w['overlapping_transcript_ids']}
        protected_conflict = [i for i in conflict if tuple(map(tuple, chains[i])) in protected]
        if conflict and len(witnesses) < 8:
            i = conflict[0]
            witnesses.append(dict(window_id=ident, species=w['species'], chain=chains[i], label=-1,
                saved_A15_gain=float(gain[i]), partial_gradient_wrt_saved_gain=float(grad[i]),
                original_A15_gradient_wrt_saved_gain=float(old_gradient[i]),
                exact_protected_geometry=i in protected_conflict))
        row = dict(window_id=ident, species=w['species'], candidates=len(local), known_positive=labels.count(1),
            known_nonmatch=labels.count(0), unknown=len(unknown), unknown_conflicting_with_gold=len(conflict),
            exact_protected_unknown_conflicting_with_gold=len(protected_conflict),
            original_known_only_loss=float(baseline['loss']), partial_observed_loss=float(partial['loss']),
            unknown_gradients_at_saved_gains=[float(grad[i]) for i in unknown],
            unknown_candidate_indices=unknown)
        results.append(row); gradients.extend(row['unknown_gradients_at_saved_gains'])
        if null_grad is not None:
            null_gradients.append(float(null_grad))
    assert sum(r['candidates'] for r in results) == 16321
    by_id = {r['window_id']: r for r in results}
    replay = [by_id[d['window_id']] for d in draws]
    unique = aggregate(results); per_draw = aggregate(replay)
    report = dict(completed=True, protocol=PROTOCOL, source_fit_job='13387908', source_readout_job='13387947',
        source_result_check_job='13387977', all_1491_natural_geometry_labels_origins_verified=True,
        actual_saved_A15_gains_used=True, arithmetic_precision='float64_recomposition_of_saved_FP32_gains',
        unique_windows=unique, original_one_pass_draws=per_draw,
        per_species_unique={s: aggregate([r for r in results if r['species'] == s]) for s in R.C.R.ALLOW},
        per_species_draws={s: aggregate([r for r in replay if r['species'] == s]) for s in R.C.R.ALLOW},
        all_original_unknown_gain_gradients_zero=True,
        partial_unknown_gain_gradient_range=[min(gradients), max(gradients)] if gradients else None,
        unknown_gold_span_conflicts_present=bool(unique['unknown_conflicting_with_gold']),
        unknown_conflict_witnesses=witnesses, optimizer_steps_this_run=0, model_or_backbone_forward=False,
        checkpoint_loaded=False, feature_tensor_or_DNA_loaded=False, labels_changed=False,
        new_candidate_search=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began,
        decision='No training: partial marginalization is not automatically the old zero-gradient protected policy.',
        limitation='Gain-variable derivatives on saved TRAIN candidates, not parameter gradients, causal model failure or performance.')
    assert all(math.isfinite(v) for v in gradients+null_gradients)
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in results:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
