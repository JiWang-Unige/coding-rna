import itertools
import unittest
import torch
from torch.nn import functional as F
from src.m28.chain_set import chain_set_loss, interval_log_partition
from src.m28.core import select_nonoverlapping


def scores(values):
    return {'logits': torch.tensor(values, dtype=torch.float64, requires_grad=True),
            'null': torch.tensor(0., dtype=torch.float64, requires_grad=True)}


def compatible(chains, indices):
    selected = sorted((chains[i][0][0], chains[i][-1][1]) for i in indices)
    return all(a[1] <= b[0] for a, b in zip(selected, selected[1:]))


class ChainSetTests(unittest.TestCase):
    def test_partition_gradients_and_MAP_match_all_compatible_sets(self):
        chains = [[(0, 9)], [(3, 15)], [(9, 18)], [(18, 24)], [(0, 24)]]
        s = scores([-.5, 1., 2., -.25, .75])
        sets = [q for size in range(6) for q in itertools.combinations(range(5), size)
                if compatible(chains, q)]
        total = [s['logits'][list(q)].sum() for q in sets]
        expected = torch.logsumexp(torch.stack(total), 0)
        actual = interval_log_partition(s['logits'], chains)
        self.assertAlmostEqual(float(actual), float(expected), places=12)
        a = torch.autograd.grad(actual, s['logits'], retain_graph=True)[0]
        b = torch.autograd.grad(expected, s['logits'], retain_graph=True)[0]
        torch.testing.assert_close(a, b, rtol=1e-12, atol=1e-12)
        picked = select_nonoverlapping(chains, s['logits'].detach().tolist())
        self.assertAlmostEqual(sum(float(s['logits'][i]) for i in picked), max(map(float, total)))
        order = [4, 2, 0, 3, 1]
        self.assertAlmostEqual(float(interval_log_partition(s['logits'][order], [chains[i] for i in order])),
                               float(expected), places=12)

    def test_disjoint_and_touching_chains_reduce_exactly_to_A6_BCE(self):
        s = scores([2., -1., 3.])
        chains = [[(0, 9)], [(9, 18)], [(21, 30)]]
        result = chain_set_loss(s, chains, [1, 0, 1])
        expected = F.binary_cross_entropy_with_logits(s['logits'], torch.tensor([1., 0., 1.], dtype=torch.float64))
        self.assertAlmostEqual(float(result['loss']), float(expected), places=12)
        self.assertEqual(result['known_positives'], 2)

    def test_one_positive_competes_with_a_set_of_two_coexisting_nonmatches(self):
        s = scores([1.5, 1., 1.])
        chains = [[(0, 18)], [(0, 9)], [(9, 18)]]
        result = chain_set_loss(s, chains, [1, 0, 0])
        expected = (torch.logsumexp(torch.tensor([0., 1.5, 1., 1., 2.], dtype=torch.float64), 0)-1.5)/3
        self.assertAlmostEqual(float(result['loss']), float(expected), places=12)
        result['loss'].backward()
        self.assertLess(float(s['logits'].grad[0]), 0.)
        self.assertTrue(bool((s['logits'].grad[1:] > 0).all()))
        self.assertEqual(select_nonoverlapping(chains, s['logits'].detach().tolist()), [1, 2])

    def test_unknown_has_zero_gradient_and_empty_or_all_unknown_is_zero(self):
        s = scores([2., 100., -1.])
        result = chain_set_loss(s, [[(0, 9)], [(0, 30)], [(12, 21)]], [1, -1, 0])
        result['loss'].backward()
        self.assertEqual(float(s['logits'].grad[1]), 0.)
        for values, chains, labels in [([], [], []), ([100.], [[(0, 9)]], [-1])]:
            s = scores(values); result = chain_set_loss(s, chains, labels)
            result['loss'].backward()
            self.assertEqual(float(result['loss']), 0.)
            self.assertTrue(torch.equal(s['logits'].grad, torch.zeros_like(s['logits'])))
            self.assertEqual(float(s['null'].grad), 0.)
        self.assertEqual(select_nonoverlapping([[(0, 9)]], [0.]), [])

    def test_incompatible_gold_or_duplicate_known_geometry_is_rejected_not_filtered(self):
        for chains, labels in [([[(0, 9)], [(3, 12)]], [1, 1]),
                               ([[(0, 9)], [(0, 9)]], [1, 0])]:
            with self.assertRaises(ValueError):
                chain_set_loss(scores([1., 2.]), chains, labels)


if __name__ == '__main__':
    unittest.main()
