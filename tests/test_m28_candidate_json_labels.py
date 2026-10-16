"""Preserve the original label policy when natural chains come from JSON."""
import json
import unittest
import numpy as np
from src.m28.labels import candidate_targets


class CandidateJSONLabels(unittest.TestCase):
    def test_json_roundtrip_preserves_positive_protected_unknown_and_callable_policy(self):
        mask = np.ones(60, dtype=bool)
        mask[42] = False
        chains = [((1, 4), (9, 12)), ((14, 20),), ((23, 27),),
                  ((29, 36),), ((40, 46),), ((50, 55),)]
        targets = dict(positive_chains=[chains[0]], protected_chains=[chains[1]],
                       unknown_chain_spans=[(30, 34)], callable_bases=mask)
        expected = [1, -1, 0, -1, -1, 0]
        self.assertEqual(candidate_targets(chains, targets), expected)
        self.assertEqual(candidate_targets(json.loads(json.dumps(chains)), targets), expected)

    def test_exact_chain_not_only_outer_span_defines_positive(self):
        correct = ((1, 4), (9, 12))
        targets = dict(positive_chains=[correct], protected_chains=[],
                       unknown_chain_spans=[], callable_bases=np.ones(20, dtype=bool))
        candidates = [correct, ((1, 12),), ((1, 5), (9, 12))]
        self.assertEqual(candidate_targets(json.loads(json.dumps(candidates)), targets), [1, 0, 0])


if __name__ == '__main__':
    unittest.main()
