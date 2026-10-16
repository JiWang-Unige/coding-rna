import itertools
import unittest
import torch
from src.m28.chain_set import chain_set_loss
from src.m28.partial_chain_set import partial_chain_set_loss


def scores(values):
    return dict(logits=torch.tensor(values, dtype=torch.float64, requires_grad=True),
                null=torch.tensor(0., dtype=torch.float64, requires_grad=True))


def compatible(chains, indices):
    spans = sorted((chains[i][0][0], chains[i][-1][1]) for i in indices)
    return all(a[1] <= b[0] for a, b in zip(spans, spans[1:]))


class PartialChainSetTests(unittest.TestCase):
    def test_both_partitions_and_gradients_match_exhaustive_partial_event(self):
        chains = [[(0, 9)], [(3, 15)], [(9, 18)], [(18, 24)], [(0, 24)]]
        labels = [1, -1, 0, -1, 0]
        s = scores([-.5, 1., 2., -.25, .75])
        sets = [q for size in range(6) for q in itertools.combinations(range(5), size)
                if compatible(chains, q)]
        observed = [q for q in sets if all((i in q) == bool(y) for i, y in enumerate(labels) if y >= 0)]
        denominator = torch.logsumexp(torch.stack([s['logits'][list(q)].sum() for q in sets]), 0)
        numerator = torch.logsumexp(torch.stack([s['logits'][list(q)].sum() for q in observed]), 0)
        expected = (denominator-numerator)/3
        result = partial_chain_set_loss(s, chains, labels)
        torch.testing.assert_close(result['log_partition'], denominator, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(result['observed_log_partition'], numerator, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(result['loss'], expected, rtol=1e-12, atol=1e-12)
        a = torch.autograd.grad(result['loss'], s['logits'], retain_graph=True)[0]
        b = torch.autograd.grad(expected, s['logits'])[0]
        torch.testing.assert_close(a, b, rtol=1e-12, atol=1e-12)

    def test_without_unknown_matches_original_A15_loss_and_gradients(self):
        chains = [[(0, 18)], [(0, 9)], [(9, 18)]]
        s = scores([1.5, 1., 1.]); labels = [1, 0, 0]
        a = partial_chain_set_loss(s, chains, labels)['loss']
        b = chain_set_loss(s, chains, labels)['loss']
        torch.testing.assert_close(a, b, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(torch.autograd.grad(a, s['logits'], retain_graph=True)[0],
                                   torch.autograd.grad(b, s['logits'])[0], rtol=1e-12, atol=1e-12)

    def test_touching_disjoint_unknown_cancels_instead_of_becoming_negative(self):
        s = scores([2., 10., -1.]); chains = [[(0, 9)], [(9, 18)], [(18, 27)]]
        result = partial_chain_set_loss(s, chains, [1, -1, 0])
        result['loss'].backward()
        self.assertAlmostEqual(float(s['logits'].grad[1]), 0., places=12)
        self.assertEqual(result['unknown_conflicting_with_gold'], 0)

    def test_conflicting_unknown_has_gradient_despite_latent_label(self):
        s = scores([1., 2.]); chains = [[(0, 12)], [(3, 15)]]
        result = partial_chain_set_loss(s, chains, [1, -1])
        result['loss'].backward()
        # Gradient descent suppresses this unlabelled competitor. This is a
        # scientific semantics difference, not an implementation defect to hide.
        self.assertGreater(float(s['logits'].grad[1]), 0.)
        self.assertEqual(result['unknown_conflicting_with_gold'], 1)

    def test_all_unknown_empty_and_conflicting_gold_keep_correct_semantics(self):
        for values, chains, labels in [([], [], []), ([100., -100.], [[(0, 9)], [(3, 12)]], [-1, -1])]:
            s = scores(values); result = partial_chain_set_loss(s, chains, labels)
            result['loss'].backward()
            self.assertEqual(float(result['loss']), 0.)
            self.assertTrue(torch.equal(s['logits'].grad, torch.zeros_like(s['logits'])))
            self.assertEqual(float(s['null'].grad), 0.)
        with self.assertRaises(ValueError):
            partial_chain_set_loss(scores([1., 2.]), [[(0, 9)], [(3, 12)]], [1, 1])


if __name__ == '__main__':
    unittest.main()
