"""A32 prospective gates and saved-pool joins, without real models or data."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_search_local_caller as S


def generator_reports():
    arab, rice = list(S.R.C.R.ALLOW)
    counts = {
        arab: dict(A3=dict(correct=1751, known_nonmatch=6305, five_or_more_correct=604),
                   A30=dict(correct=1752, known_nonmatch=6305, five_or_more_correct=605)),
        rice: dict(A3=dict(correct=629, known_nonmatch=6261, five_or_more_correct=184),
                   A30=dict(correct=630, known_nonmatch=6261, five_or_more_correct=185))}
    checks = S.P.R.gate(counts)
    report = dict(natural_counts_per_species=counts, predeclared_TRAIN_generator_gate=checks,
                  TRAIN_generator_supported=all(checks.values()))
    checked = dict(actual_saved_pool_identity_and_five_fixed_gate_checks_verified=True,
        actual_new_and_original_natural_chain_labels_recomputed=True,
        natural_counts_per_species=counts, original_frozen_TRAIN_generator_gate=checks,
        TRAIN_generator_supported=all(checks.values()))
    return report, checked


def value(chains, scores):
    return dict(chains=chains, scores=scores, selected=S.R.select_nonoverlapping(chains, scores))


class SearchLocalCallerTests(unittest.TestCase):
    def test_each_failed_original_generator_gate_closes_a32(self):
        report, checked = generator_reports(); S.require_generator_gate(report, checked)
        arab, rice = list(S.R.C.R.ALLOW)
        for species, field, number in ((arab, 'correct', 1751), (rice, 'correct', 629),
                (arab, 'five_or_more_correct', 603), (arab, 'known_nonmatch', 6306),
                (rice, 'known_nonmatch', 6262)):
            with self.subTest(species=species, field=field):
                report, checked = generator_reports()
                checked['natural_counts_per_species'][species]['A30'][field] = number
                checks = S.P.R.gate(checked['natural_counts_per_species'])
                report['predeclared_TRAIN_generator_gate'] = checks
                checked['original_frozen_TRAIN_generator_gate'] = checks
                report['TRAIN_generator_supported'] = checked['TRAIN_generator_supported'] = False
                self.assertEqual(sum(not flag for flag in checks.values()), 1)
                with self.assertRaisesRegex(ValueError, 'A32 not executed'):
                    S.require_generator_gate(report, checked)

    def test_source_gate_failure_precedes_weights_cache_and_population_access(self):
        with patch.object(S, 'sources', side_effect=ValueError('A32 not executed')), \
                patch.object(S, 'load_model') as model, patch.object(S.R.C.R, 'rows') as rows, \
                patch.object(S.torch, 'load') as cache:
            with self.assertRaisesRegex(ValueError, 'A32 not executed'):
                S.readout()
            model.assert_not_called(); rows.assert_not_called(); cache.assert_not_called()

    def test_unknown_competes_and_changed_pool_or_saved_baseline_is_rejected(self):
        w = {k: k for k in S.R.C.R.GEOMETRY}; w['species'] = next(iter(S.R.C.R.ALLOW))
        chains = [[[0, 9]], [[0, 12]], [[20, 29]]]
        scored = dict(window=w, reference_used_in_scoring=False, **value(chains, [1., 10., 1.]))
        raw = dict(window=w, chains=chains, reference_used_in_pool_generation=False)
        labelled = dict(window=w, candidates=[dict(chain=c, label=y) for c, y in zip(chains, [1, -1, 0])])
        old = value(chains[:1], [1.])
        baseline = dict(window=w, variants=dict(A31_A3=old))
        baseline_labelled = dict(window=w, labels=dict(A31_A3=[1]),
                                counts=dict(A31_A3=S.R.selected_counts(old, [1])))
        counts, labels = S.paired_counts(scored, raw, labelled, baseline, baseline_labelled)
        self.assertEqual(labels, [1, -1, 0])
        self.assertEqual(counts['A32_A30']['selected_correct'], 0)
        self.assertEqual(counts['A32_A30']['selected_unknown'], 1)
        self.assertEqual(counts['A32_A30']['selected_known_nonmatch'], 1)
        changed = dict(raw, chains=list(reversed(chains)))
        with self.assertRaisesRegex(ValueError, 'closed A30 pool'):
            S.paired_counts(scored, changed, labelled, baseline, baseline_labelled)
        changed = deepcopy(baseline_labelled); changed['counts']['A31_A3']['selected_correct'] = 0
        with self.assertRaisesRegex(ValueError, 'per-window counts'):
            S.paired_counts(scored, raw, labelled, baseline, changed)

    def test_five_caller_gates_require_strict_f1_and_long_gain_and_preserve_precision(self):
        refs = {s: dict(eligible=100) for s in S.R.C.R.ALLOW}
        old = dict(selected_correct=50, selected_known_nonmatch=50, selected_five_or_more_correct=10)
        counts = {s: dict(A31_A3=dict(old), A32_A30=dict(old)) for s in refs}
        checks = S.gate(counts, refs)
        self.assertEqual(len(checks), 5); self.assertEqual(sum(checks.values()), 2)
        for values in counts.values():
            values['A32_A30'] = dict(selected_correct=60, selected_known_nonmatch=40,
                                     selected_five_or_more_correct=11)
        self.assertTrue(all(S.gate(counts, refs).values()))
        rice = list(refs)[1]; counts[rice]['A32_A30']['selected_known_nonmatch'] = 61
        self.assertFalse(S.gate(counts, refs)[rice+'_proxy_precision_not_decreased'])


if __name__ == '__main__':
    unittest.main()
