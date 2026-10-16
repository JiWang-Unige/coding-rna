import unittest
import torch
from torch.nn import functional as F
from src.m28.chain_competition import chain_competition_loss as loss


def scores(values):
    return {'logits': torch.tensor(values, dtype=torch.float64, requires_grad=True),
            'null': torch.tensor(0., dtype=torch.float64, requires_grad=True)}


class ChainCompetitionTests(unittest.TestCase):
    def test_disjoint_positives_are_independent_not_one_correct_window_class(self):
        s = scores([2., 4.])
        r = loss(s, [[(0, 9)], [(12, 21)]], [1, 1])
        self.assertAlmostEqual(float(r['loss']), float(F.softplus(-s['logits']).mean()))
        self.assertEqual(r['conflict_pairs'], 0)

    def test_unknown_and_touching_negative_do_not_enter_positive_competition(self):
        s = scores([2., 100., 4.])
        r = loss(s, [[(0, 9)], [(3, 12)], [(9, 18)]], [1, -1, 0])
        expected = (F.softplus(-s['logits'][0])+F.softplus(s['logits'][2]))/2
        self.assertAlmostEqual(float(r['loss']), float(expected))
        r['loss'].backward()
        self.assertEqual(float(s['logits'].grad[1]), 0.)
        self.assertEqual(r['unknown_excluded'], 1)

    def test_actual_conflict_gradient_raises_correct_and_lowers_nonmatch(self):
        s = scores([-1., 2.])
        r = loss(s, [[(0, 9)], [(3, 12)]], [1, 0])
        expected = (torch.logsumexp(torch.tensor([0., -1., 2.], dtype=torch.float64), 0)+1+F.softplus(s['logits'][1]))/2
        self.assertAlmostEqual(float(r['loss']), float(expected))
        self.assertEqual(r['conflict_pairs'], 1)
        r['loss'].backward()
        self.assertLess(float(s['logits'].grad[0]), 0.)
        self.assertGreater(float(s['logits'].grad[1]), 0.)

    def test_null_anchor_prevents_pure_relative_ranking_shift_invariance(self):
        a = loss(scores([-2., -3.]), [[(0, 9)], [(3, 12)]], [1, 0])
        b = loss(scores([2., 1.]), [[(0, 9)], [(3, 12)]], [1, 0])
        self.assertNotEqual(float(a['positive_competition']), float(b['positive_competition']))
        self.assertGreater(float(b['negative_null_anchor']), float(a['negative_null_anchor']))

    def test_unknown_or_empty_pool_has_differentiable_zero(self):
        for values, chains, labels in [([100.], [[(0, 9)]], [-1]), ([], [], [])]:
            s = scores(values)
            r = loss(s, chains, labels)
            self.assertEqual(float(r['loss']), 0.)
            r['loss'].backward()
            self.assertTrue(torch.equal(s['logits'].grad, torch.zeros_like(s['logits'])))
            self.assertEqual(float(s['null'].grad), 0.)


if __name__ == '__main__':
    unittest.main()
