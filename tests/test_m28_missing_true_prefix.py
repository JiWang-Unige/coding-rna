import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from analyze_missing_true_prefix_train import summarize, aggregate


def profiles(passing=(), supported=True):
    values = []
    for q in range(16):
        if not supported:
            values.append(dict(query_index=q, grammar_supported=False, unsupported_reason='end_stop', events=[]))
            continue
        passed = q in passing
        exclusion = None if passed else dict(exon=1, stage='start', rank=5)
        values.append(dict(query_index=q, grammar_supported=True,
            events=[dict(entry_legal_rank=1 if passed else 5, end_same_kind_rank=2)],
            local_width4_necessary_conditions_passed=passed, first_local_exclusion=exclusion))
    return values


class MissingPrefixTests(unittest.TestCase):
    def test_pair_and_query_denominators_stay_separate(self):
        a = summarize({}, profiles(), 'missing')
        b = summarize({}, profiles((15,)), 'missing')
        self.assertEqual(a['mechanism'], 'all_queries_locally_excluded')
        self.assertEqual(b['mechanism'], 'some_query_passes_local_necessary_conditions')
        self.assertEqual(b['local_survivor_queries'], [15])
        self.assertEqual(b['best_query_max_necessary_local_rank'], 2)
        total = aggregate([a, b])
        self.assertEqual(total['reference_window_pairs'], 2)
        self.assertEqual(total['mechanisms']['all_queries_locally_excluded'], 1)
        self.assertEqual(total['query_first_exclusion_stages'], {'start': 31})

    def test_unsupported_is_not_rank_failure(self):
        row = summarize({}, profiles(supported=False), 'missing')
        self.assertEqual(row['mechanism'], 'grammar_unsupported')
        self.assertIsNone(row['best_query_max_necessary_local_rank'])
        self.assertEqual(row['query_first_exclusion_stages'], {})
        self.assertEqual(row['unsupported_reasons'], ['end_stop'])

    def test_natural_beam_control_must_have_a_necessary_local_survivor(self):
        with self.assertRaisesRegex(ValueError, 'Saved natural beam hit'):
            summarize({}, profiles(), 'natural_hit_control')
        value = summarize({}, profiles((15,)), 'natural_hit_control')
        self.assertEqual(value['local_survivor_queries'], [15])


if __name__ == '__main__':
    unittest.main()

