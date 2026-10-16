import copy
import math
import unittest

import torch

from src.m28.core import select_nonoverlapping
from src.m28.natural_competitor_pool import augment_training_pool
from src.m28.structured_emission import structured_emission_loss, trainable_path_scores
from census_natural_competitor_support import census


class NaturalCompetitorPoolTests(unittest.TestCase):
    def test_saved_membership_counts_not_training_gold_or_selected_identity(self):
        w = dict(id='artificial', species='arabidopsis_thaliana', split='train')
        old = dict(window=w, reference_used_in_scoring=False,
                   variants={'A31_A3': dict(chains=[[[0, 12]], [[24, 30]], [[32, 38]]])})
        new = dict(window=w, reference_used_in_scoring=False, original_ordinal=0,
                   chains=[[[0, 12]], [[0, 6]], [[32, 38]], [[40, 46]]], selected=[0, 1, 2, 3])
        lab = dict(window=w, labels={'A31_A3': [1, 0, -1], 'A42': [1, 0, -1, 1]},
                   eligible_references=2, eligible_five_or_more=0)
        result = census([old], [new], [lab])
        self.assertEqual((result['windows'], result['original_natural_candidates'],
                          result['A42_natural_selected_candidates']), (1, 3, 4))
        self.assertEqual(result['pooled'], dict(correct=dict(present=1, new=1),
            known_nonmatch=dict(present=0, new=1), unknown=dict(present=1, new=0)))
        self.assertEqual(result['original_natural_labels_pooled'], dict(correct=1, known_nonmatch=1, unknown=1))
        # A positive absent from the *natural* source may be an original injected
        # gold: a membership count must not assert an expanded training target.
        with self.assertRaises(ValueError): census([old], [], [lab])
        changed = copy.deepcopy(lab); changed['labels']['A42'][0] = 0
        with self.assertRaises(ValueError): census([old], [new], [changed])
        changed = copy.deepcopy(new); changed['original_ordinal'] = 1
        with self.assertRaises(ValueError): census([old], [changed], [lab])

    def test_whole_geometry_join_order_original_gold_and_unknown_preserved(self):
        original = [[[0, 6], [12, 18]], [[24, 30]], [[32, 38]]]
        natural = [[[24, 30]], [[0, 6]], [[0, 18]], [[0, 6], [12, 18]]]
        labels = [1, 0, -1]; natural_labels = [0, 0, -1, 1]
        before = copy.deepcopy((original, labels, natural, natural_labels))
        value = augment_training_pool(original, labels, natural, natural_labels)
        self.assertEqual(value, dict(chains=original+[natural[1], natural[2]],
                                    labels=[1, 0, -1, 0, -1], added_natural_indices=[1, 2]))
        self.assertEqual((original, labels, natural, natural_labels), before)
        self.assertEqual([c for c, y in zip(value['chains'], value['labels']) if y == 1], [original[0]])
        # The same outer span but different exon geometry is a new competitor.
        self.assertNotEqual(value['chains'][0], value['chains'][4])
        value['chains'][0][0][0] = 1
        self.assertEqual(original[0][0][0], 0)

    def test_new_known_competitor_enters_loss_unknown_stays_predictable(self):
        original = [[[0, 12]], [[24, 30]]]
        natural = [[[0, 6]], [[0, 18]], [[0, 12]]]
        value = augment_training_pool(original, [1, 0], natural, [0, -1, 1])
        logits = torch.zeros(40, 15, dtype=torch.float64, requires_grad=True)
        with torch.no_grad():
            logits[9:12, 14] = 2/3  # Original gold gain 2.
            logits[3:6, 14] = 1    # New known competitor gain 3.
            logits[15:18, 14] = 100/3  # Unknown natural competitor gain 100.
        before = logits.detach().clone()
        result = structured_emission_loss(logits, value['chains'], value['labels'])
        gains = trainable_path_scores(logits, value['chains'])
        torch.testing.assert_close(gains, torch.tensor([2., 0., 3., 100.], dtype=torch.float64),
                                   rtol=1e-12, atol=1e-10)
        # Two mutually exclusive known chains in the first region, independent
        # background candidate in the second; all compatible sets include null.
        expected = (math.log1p(math.exp(2)+math.exp(3)) + math.log(2) - 2)/3
        self.assertAlmostEqual(float(result['loss']), expected, places=12)
        self.assertEqual((result['known_candidates'], result['known_positives'],
                          result['known_nonmatches'], result['unknown_excluded']), (3, 1, 2, 1))
        gradient, = torch.autograd.grad(result['loss'], logits, retain_graph=True)
        original_result = structured_emission_loss(logits, original, [1, 0])
        original_gradient, = torch.autograd.grad(original_result['loss'], logits)
        self.assertEqual(float(original_gradient[3:6, 14].abs().sum()), 0.)
        self.assertGreater(float(gradient[3:6, 14].min()), 0.)
        self.assertLess(float(gradient[9:12, 14].max()), 0.)
        self.assertEqual(float(gradient[15:18, 14].abs().sum()), 0.)
        torch.testing.assert_close(logits.detach(), before, rtol=0, atol=0)
        # A43 support assembly is not a prediction-time label filter.
        selected = select_nonoverlapping(value['chains'], gains.detach().tolist())
        self.assertEqual([value['labels'][i] for i in selected], [-1])

    def test_real_support_errors_rejected_without_changing_gold(self):
        chain = [[0, 6]]; other = [[12, 18]]
        cases = [
            ([chain], [], [], []),
            ([chain], [1], [other], []),
            ([chain], [1], [chain], [0]),
            ([chain], [-1], [chain], [0]),
            ([chain], [0], [other], [1]),
            ([chain, chain], [1, 1], [], []),
            ([chain], [1], [other, other], [0, 0]),
            ([chain], [2], [], []),
        ]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                augment_training_pool(*args)

    def test_empty_original_draw_and_all_unknown_are_not_dropped(self):
        for old_chains, old_labels, new_chains, new_labels in (
                ([], [], [], []),
                ([], [], [[[0, 6]]], [-1]),
                ([[[0, 6]]], [-1], [[[6, 12]]], [-1])):
            value = augment_training_pool(old_chains, old_labels, new_chains, new_labels)
            logits = torch.zeros(18, 15, dtype=torch.float64, requires_grad=True)
            result = structured_emission_loss(logits, value['chains'], value['labels'])
            self.assertEqual(float(result['loss']), 0.)
            result['loss'].backward()
            self.assertEqual(float(logits.grad.abs().sum()), 0.)
        # A previously empty draw may acquire evaluable background competition;
        # the original per-known-candidate normalization then uses one, not zero.
        value = augment_training_pool([], [], [[[0, 6]]], [0])
        logits = torch.zeros(12, 15, dtype=torch.float64, requires_grad=True)
        result = structured_emission_loss(logits, value['chains'], value['labels'])
        self.assertAlmostEqual(float(result['loss']), math.log(2), places=12)
        self.assertEqual((result['known_candidates'], result['known_positives']), (1, 0))
        result['loss'].backward()
        self.assertGreater(float(logits.grad[0:3, 13].min()), 0.)


if __name__ == '__main__':
    unittest.main()
