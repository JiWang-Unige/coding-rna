"""Closed-pool oracle math and reference preservation, without real data."""
from collections import Counter
from copy import deepcopy
from itertools import combinations
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_local_caller_selector_bound as D


def row(ident, chain, strand='+', species='synthetic'):
    return dict(id=ident, chain=chain, species=species, seqid='chr', strand=strand)


class SelectorBoundTests(unittest.TestCase):
    def test_oracle_uses_spans_and_touching_is_compatible_and_maximal(self):
        values = [row('a', [[0, 2], [8, 10]]), row('b', [[4, 6]]), row('c', [[10, 12]]),
                  row('d', [[0, 12]]), row('e', [[12, 14]])]
        selected = D.oracle(values)
        expected = 0
        for n in range(len(values)+1):
            for subset in combinations(values, n):
                spans = sorted((r['chain'][0][0], r['chain'][-1][1]) for r in subset)
                if all(a[1] <= b[0] for a, b in zip(spans, spans[1:])):
                    expected = max(expected, n)
        self.assertEqual(expected, 3); self.assertEqual(len(selected), expected)
        self.assertFalse({'a', 'b'} <= {r['id'] for r in selected})

    def test_strands_and_species_do_not_compete(self):
        values = [row('plus', [[0, 10]]), row('minus', [[0, 10]], '-'),
                  row('other', [[0, 10]], species='other')]
        self.assertEqual(len(D.oracle(values)), 3)

    def test_strict_native_exceedance_excludes_equal_f1(self):
        baseline = dict(tp=2, predicted=4, reference=4)
        result = D.comparison(2, 4, baseline)
        self.assertTrue(result['ideal_F1_strictly_above_native_ANNEVO'])
        self.assertEqual(result['minimum_TP_at_precision_one_to_strictly_exceed_native_ANNEVO'], 2)
        self.assertEqual(result['maximum_nonexact_predictions_if_TP_equals_oracle_and_F1_strictly_exceeds_ANNEVO'], 1)
        self.assertGreater(4/7, .5); self.assertEqual(4/8, .5)
        self.assertFalse(D.comparison(1, 4, baseline)['ideal_F1_strictly_above_native_ANNEVO'])
        with self.assertRaisesRegex(ValueError, 'denominator'):
            D.comparison(2, 3, baseline)

    def test_full_unsupported_reference_is_retained_and_duplicate_key_is_rejected(self):
        rows = [dict(species='synthetic', id='a', seqid='chr', strand='+', CDS_intervals='[[0,2],[8,10]]',
                     unsupported_strata='long_chain', A6_state='not_free', A7_state='not_free',
                     A6_selected='False', A7_selected='False'),
                dict(species='synthetic', id='b', seqid='chr', strand='+', CDS_intervals='[[20,30]]',
                     unsupported_strata='', A6_state='selected', A7_state='selected',
                     A6_selected='True', A7_selected='True')]
        counts = Counter(not_free=1, selected=1)
        saved = dict(per_species=dict(synthetic={a: {k: counts[k] for k in D.STATES} for a in ('A6', 'A7')}))
        with patch.object(D, 'REFERENCES', dict(synthetic=2)):
            actual = D.prepare(rows, saved)
            self.assertEqual(len(actual), 2); self.assertEqual(actual[0]['unsupported_strata'], 'long_chain')
            changed = deepcopy(rows); changed[1]['CDS_intervals'] = changed[0]['CDS_intervals']
            with self.assertRaisesRegex(ValueError, 'exact-chain keys'):
                D.prepare(changed, saved)


if __name__ == '__main__':
    unittest.main()
