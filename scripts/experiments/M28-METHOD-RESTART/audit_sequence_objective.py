#!/usr/bin/env python3
"""Saved-score TRAIN assignment counterfactuals; no fitting or model forward."""
from collections import Counter
import json
import math
import time
import numpy as np
from scipy.optimize import linear_sum_assignment
import analyze_query_beam_train as B
from src.m28.query_model import unmatched_query_policy

R = B.R
OUT = R.ROOT/'outputs/M28-A18-SEQUENCE-OBJECTIVE-AUDIT-RECOVERY-R1'
A5 = R.ROOT/'outputs/M28-A5-TRUE-PREFIX-DIAGNOSTIC'
A16 = R.ROOT/'outputs/M28-A16-LEGAL-SUPPORT-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a18_sequence_objective_audit.md'


def costs(boundary, termination, exons, logits, unmatched):
    n, q = boundary.shape
    if n == 0:
        return {arm: np.empty((0, q)) for arm in ('MEAN', 'EVENT_RATE', 'SEQUENCE_SUM')}
    delta = -logits if unmatched == 0 else np.logaddexp(0., -logits)
    existence = delta/(q if unmatched == 0 else n)
    return dict(MEAN=existence+(boundary+termination)/n,
                EVENT_RATE=existence+(2*boundary+termination)/n,
                SEQUENCE_SUM=existence+exons[:, None]*(2*boundary+termination)/n)


def assignment(matrix):
    rows, queries = linear_sum_assignment(matrix)
    assert list(rows) == list(range(len(matrix))) and len(set(queries)) == len(matrix)
    return queries.tolist()


def checks():
    b = np.array([[2., 2.]]); t = np.array([[.5, .5]])
    m = costs(b, t, np.array([3]), np.zeros(2), 0)
    assert np.array_equal(m['MEAN'], [[2.5, 2.5]])
    assert np.array_equal(m['EVENT_RATE'], [[4.5, 4.5]])
    assert np.array_equal(m['SEQUENCE_SUM'], [[13.5, 13.5]])
    b = np.array([[0., 1.], [.2, .9]]); t = np.zeros_like(b)
    m = costs(b, t, np.array([1, 10]), np.zeros(2), 0)
    assert assignment(m['MEAN']) == assignment(m['EVENT_RATE']) == [0, 1]
    assert assignment(m['SEQUENCE_SUM']) == [1, 0]
    logits = np.array([2., -2.]); b = np.zeros((1, 2))
    unknown = costs(b, b, np.array([1]), logits, -1)['MEAN']
    known = costs(b, b, np.array([1]), logits, 0)['MEAN']
    assert np.allclose(unknown, np.logaddexp(0., -logits)[None])
    assert np.array_equal(known, -logits[None]/2)
    empty = costs(np.empty((0, 16)), np.empty((0, 16)), np.array([]), np.zeros(16), 0)
    assert all(m.shape == (0, 16) and assignment(m) == [] for m in empty.values())
    print('three analytic objective/assignment/unknown checks and zero-positive regression PASS', flush=True)


def main():
    checks(); began = time.perf_counter()
    for folder, job in ((A5, '13359200'), (A16, '13367351'),
                        (R.ROOT/'outputs/M28-A16-LEGAL-SUPPORT-RESULT-CHECK', '13367450')):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
        assert json.loads((folder/'summary.json').read_text())['completed']
    source = json.loads((A16/'summary.json').read_text())
    assert source['generator_checkpoint'] == 'outputs/M28-A3-FIRST-FIT/step_004608.pt'
    assert source['windows'] == 8 and source['reference_window_pairs'] == 27 and source['query_profiles'] == 432
    assert source['optimizer_steps_this_run'] == 0 and not source['DEV_forward_used'] and not source['test_setaria_used']
    predictions = R.rows(B.OUT/'predictions.jsonl'); refs = R.rows(A5/'references.jsonl')
    old = R.rows(A5/'query_profiles.jsonl'); legal = R.rows(A16/'query_profiles.jsonl')
    key = lambda r: (r['window_id'], r['reference_id'], r['query_index'])
    original = {key(r): r for r in old}; profiles = {key(r): r for r in legal}
    expected = {(r['window_id'], r['reference_id'], q) for r in refs for q in range(16)}
    assert len(refs) == 27 and len(old) == len(legal) == len(expected) == 432
    assert set(original) == set(profiles) == expected
    assert len(predictions) == 8 and len({p['window']['id'] for p in predictions}) == 8
    assert {p['window']['id'] for p in predictions} == {r['window_id'] for r in R.rows(R.READY/'witnesses.jsonl')}
    windows, catalogs = B.load_selected_labels(predictions)
    draw_counts = Counter(r['window_id'] for r in R.rows(R.Q.DRAWS))
    reference_rows = []; window_rows = []; maximum_identity_error = 0.
    for pred in predictions:
        ident = pred['window']['id']; w = windows[ident]
        local = [r for r in refs if r['window_id'] == ident]; n = len(local)
        assert 0 <= n <= 16 and n == pred['expected_positive_chains'] and not pred['reference_used_in_decoding']
        targets = B.build_targets(w, catalogs[w['species']], pred['sequence'])
        assert {B.chain_key(r['chain']) for r in local} == {B.chain_key(c) for c in targets['positive_chains']}
        assert len(targets['positive_chains']) == n
        policy = unmatched_query_policy(targets, w['valid_bases'])
        boundary = np.empty((n, 16)); termination = np.empty((n, 16))
        exon_counts = np.array([r['exons'] for r in local])
        assert [r['query_index'] for r in pred['greedy_queries']] == list(range(16))
        logits = np.array([r['query_logit'] for r in pred['greedy_queries']])
        for i, ref in enumerate(local):
            assert ref['exons'] == len(ref['chain']) and ref['original_A3_draw_occurrences_per_pass'] == draw_counts[ident]
            for q in range(16):
                p = profiles[(ident, ref['reference_id'], q)]; oldp = original[(ident, ref['reference_id'], q)]
                assert p['reference_conditioned'] and oldp['grammar_supported'] and oldp['query_logit'] == logits[q]
                events = p['events']; assert len(events) == ref['exons']
                assert [[e['start'], e['end']] for e in events] == ref['chain']
                boundary[i, q] = -sum(e['teacher_entry_target_logp']+e['end_logp'] for e in events)/(2*ref['exons'])
                termination[i, q] = -sum(e['eos_logp'] for e in events)/ref['exons']
                raw = sum(e['entry_logp']+e['end_logp']+e['eos_logp'] for e in events)
                delta = sum(e['decoder_vs_teacher_lower_target_logp_change'] for e in events)
                error = max(abs(raw-p['raw_path_log_score']), abs(raw-oldp['path_log_score']),
                    abs(raw-(-ref['exons']*(2*boundary[i, q]+termination[i, q])+delta)))
                assert math.isfinite(error) and error <= 1e-5
                maximum_identity_error = max(maximum_identity_error, error)
        matrices = costs(boundary, termination, exon_counts, logits, policy['unmatched_label'])
        assert all(np.isfinite(m).all() for m in matrices.values())
        allocations = {arm: assignment(m) for arm, m in matrices.items()}
        window_rows.append(dict(window_id=ident, species=w['species'], reference_pairs=n,
            unmatched_policy=policy, original_draw_occurrences_per_pass=draw_counts[ident], assignments=allocations,
            matrix_costs={a: m.tolist() for a, m in matrices.items()},
            assigned_total_cost={a: float(m[np.arange(n), np.array(allocations[a], dtype=int)].sum())
                                 for a, m in matrices.items()}))
        for i, ref in enumerate(local):
            assigned = {a: qs[i] for a, qs in allocations.items()}
            reference_rows.append(dict(**ref, assigned_queries=assigned,
                event_rate_changed=assigned['EVENT_RATE'] != assigned['MEAN'],
                sequence_sum_changed=assigned['SEQUENCE_SUM'] != assigned['MEAN'],
                length_changed_from_event_rate=assigned['SEQUENCE_SUM'] != assigned['EVENT_RATE'],
                assigned_local_necessary_pass={a: original[(ident, ref['reference_id'], q)][
                    'local_width4_necessary_conditions_passed'] for a, q in assigned.items()}))
    assert len(reference_rows) == 27 and len(window_rows) == 8
    def summarize(values):
        return dict(reference_window_pairs=len(values),
            **{k: sum(r[k] for r in values) for k in ('event_rate_changed', 'sequence_sum_changed', 'length_changed_from_event_rate')},
            assigned_local_necessary_pass={a: sum(r['assigned_local_necessary_pass'][a] for r in values) for a in matrices})
    for name, values in (('windows.jsonl', window_rows), ('references.jsonl', reference_rows)):
        with (OUT/name).open('x') as handle:
            for row in values:
                handle.write(json.dumps(row, allow_nan=False)+'\n')
    report = dict(completed=True, protocol=PROTOCOL, windows=8, reference_window_pairs=27, query_profiles=432,
        positive_windows=sum(r['reference_pairs'] > 0 for r in window_rows),
        zero_positive_windows=sum(r['reference_pairs'] == 0 for r in window_rows),
        source_jobs=dict(A5='13359200', A16='13367351', A16_check='13367450'),
        source_checkpoint=source['generator_checkpoint'], all_original_identities_and_raw_scores_verified=True,
        score_decomposition_max_absolute_error=maximum_identity_error, pooled=summarize(reference_rows),
        per_species={s: summarize([r for r in reference_rows if r['species'] == s]) for s in R.ALLOW},
        by_exon_class={c: summarize([r for r in reference_rows if r['exon_class'] == c]) for c in ('one', 'two_to_four', 'five_or_more')},
        original_drawn_reference_pairs=sum(r['original_A3_draw_occurrences_per_pass'] > 0 for r in reference_rows),
        optimizer_steps_this_run=0, new_model_or_backbone_forward=False, new_candidate_search=False,
        checkpoint_written=False, DEV_forward_used=False, test_setaria_used=False, seconds=time.perf_counter()-began,
        limitation='Saved scalar FP32 probes recombined in float64, not native batched training Hungarian replay. Forced TRAIN pressure histories, not natural predictions, improved accuracy or causal attribution.')
    R.write_summary(OUT, report); print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
