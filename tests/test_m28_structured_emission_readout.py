from pathlib import Path
import copy
import sys
import unittest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import readout_structured_emission_train as D
import verify_structured_emission_train as V


def value(chains, scores):
    return dict(chains=chains, scores=scores, selected=D.R.select_nonoverlapping(chains, scores))


def toy():
    species = next(iter(D.R.C.R.ALLOW))
    w = dict(id='synthetic', species=species)
    chains = [[[0, 9]], [[0, 18]], [[24, 30]]]; labels = [1, -1, 0]
    old = value(chains, [1., -1., -1.])
    prior = dict(window=w, variants=dict(A31_A3=old))
    lab = dict(window=w, labels=dict(A31_A3=labels),
        counts=dict(A31_A3=D.R.selected_counts(old, labels)), eligible_references=1, eligible_five_or_more=0)
    new = dict(window=w, reference_used_in_scoring=False, **value(chains, [1., 3., -1.]))
    return species, new, prior, lab


class StructuredReadoutTests(unittest.TestCase):
    def test_saved_logits_independent_full_window_replay_rejects_wrong_score(self):
        torch.manual_seed(2)
        logits = torch.randn(36, 15)
        chains = [[[0, 6], [12, 18]], [[24, 30]]]
        expected = D.local_path_scores(logits, chains).tolist()
        self.assertLessEqual(V.verify_direct_scores(logits, chains, expected), 1e-8)
        with self.assertRaises(AssertionError):
            V.verify_direct_scores(logits, chains, [expected[0]+1, expected[1]])

    def test_unknown_competes_and_is_not_reported_as_known_nonmatch(self):
        species, new, prior, lab = toy()
        rows, counts, refs = D.paired_readout([new], [prior], [lab])
        self.assertEqual(counts[species]['A36']['selected_correct'], 0)
        self.assertEqual(counts[species]['A36']['selected_unknown'], 1)
        self.assertEqual(counts[species]['A36']['selected_known_nonmatch'], 0)
        self.assertEqual(refs[species]['eligible'], 1)
        self.assertEqual(rows[0]['labels'], [1, -1, 0])
        bad = {**new, 'selected': [0]}
        with self.assertRaisesRegex(ValueError, 'reconcile'):
            D.paired_readout([bad], [prior], [lab])

    def test_reordered_or_injected_geometry_and_incomplete_join_are_rejected(self):
        _, new, prior, lab = toy()
        for chains in (new['chains'][::-1], new['chains']+[[[32, 38]]]):
            bad = {**new, 'chains': chains}
            with self.assertRaisesRegex(ValueError, 'geometry'):
                D.paired_readout([bad], [prior], [lab])
        with self.assertRaisesRegex(ValueError, 'paired'):
            D.paired_readout([new], [prior], [])

    def test_five_conditions_use_stronger_baseline_and_each_species(self):
        counts = {}
        for s in D.R.C.R.ALLOW:
            counts[s] = dict(A31_A3=dict(selected_correct=10, selected_known_nonmatch=4,
                selected_five_or_more_correct=3), A36=dict(selected_correct=11,
                selected_known_nonmatch=4, selected_five_or_more_correct=4),
                A6_auxiliary=dict(selected_correct=1, selected_known_nonmatch=10,
                    selected_five_or_more_correct=0))
        self.assertEqual(len(D.gate(counts)), 5)
        self.assertTrue(all(D.gate(counts).values()))
        s = next(iter(counts))
        for field, value in (('selected_correct', 10), ('selected_known_nonmatch', 5)):
            bad = copy.deepcopy(counts); bad[s]['A36'][field] = value
            self.assertFalse(all(D.gate(bad).values()))
        bad = copy.deepcopy(counts)
        for s in bad:
            bad[s]['A36']['selected_five_or_more_correct'] = 3
        self.assertFalse(all(D.gate(bad).values()))


if __name__ == '__main__':
    unittest.main()
