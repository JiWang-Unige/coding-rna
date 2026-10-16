"""A39 objective semantics only; no research data, model fit or optimizer."""
import unittest
import torch
from src.m28.core import chain_loss
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss


def score(values):
    return dict(logits=torch.tensor(values, dtype=torch.float64, requires_grad=True),
                null=torch.tensor(.2, dtype=torch.float64, requires_grad=True))


class NullAnchoredChainSetTests(unittest.TestCase):
    def test_disjoint_known_chains_keep_BCE_loss_and_gradient_scale(self):
        chains = [[[0, 10]], [[10, 20]], [[20, 30]], [[0, 30]]]
        value = score([.4, -.7, 1.2, 12.]); labels = [1, 1, 0, -1]
        result = null_anchored_chain_set_loss(value, chains, labels)
        reference = chain_loss(value, labels)
        actual = torch.autograd.grad(result['loss'], value['logits'], retain_graph=True)[0]
        expected = torch.autograd.grad(reference, value['logits'])[0]
        torch.testing.assert_close(result['loss'], reference, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(actual, expected, rtol=1e-12, atol=1e-12)
        self.assertEqual((result['known_candidates'], result['unknown_excluded']), (3, 1))
        self.assertEqual(float(actual[-1]), 0.)

    def test_conflict_gradient_is_exact_mean_of_absolute_and_set_derivatives(self):
        value = score([.4, .7, 12.]); labels = [1, 0, -1]
        result = null_anchored_chain_set_loss(value, [[[0, 20]], [[0, 19]], [[0, 20]]], labels)
        result['loss'].backward(); gain = value['logits'].detach()
        targets = gain.new_tensor([1., 0.])
        marginal = torch.cat([gain.new_zeros(1), gain[:2]]).softmax(0)[1:]
        expected = torch.cat([(gain[:2].sigmoid()+marginal-2*targets)/4, gain.new_zeros(1)])
        torch.testing.assert_close(value['logits'].grad, expected, rtol=1e-12, atol=1e-12)
        self.assertEqual(float(value['logits'].grad[2]), 0.)
        self.assertLess(float(value['logits'].grad[0]), 0.)
        self.assertGreater(float(value['logits'].grad[1]), 0.)

    def test_both_absolute_nonmatch_anchor_and_losing_positive_ranking_remain(self):
        value = score([14., 3., 3.]); chains = [[[0, 30]], [[0, 15]], [[15, 30]]]
        result = null_anchored_chain_set_loss(value, chains, [1, 0, 0])
        actual = torch.autograd.grad(result['loss'], value['logits'])[0]
        self.assertGreater(float(actual[1]), .15)
        self.assertGreater(float(actual[2]), .15)
        value = score([.5, 1.]); chains = [[[0, 20]], [[0, 19]]]
        result = null_anchored_chain_set_loss(value, chains, [1, 0])
        actual = torch.autograd.grad(result['loss'], value['logits'], retain_graph=True)[0]
        pointwise = torch.autograd.grad(result['mean_known_BCE'], value['logits'])[0]
        self.assertLess(float(actual[0]), float(pointwise[0]))

    def test_empty_unknown_and_background_keep_original_known_mean_semantics(self):
        for values, labels, chains in (([], [], []), ([7., -3.], [-1, -1], [[[0, 10]], [[0, 9]]]),
                                       ([.4, -.9], [0, 0], [[[0, 10]], [[10, 20]]])):
            value = score(values); result = null_anchored_chain_set_loss(value, chains, labels)
            result['loss'].backward()
            if any(y == 0 for y in labels):
                torch.testing.assert_close(value['logits'].grad, value['logits'].detach().sigmoid()/2,
                                           rtol=1e-12, atol=1e-12)
            else:
                self.assertEqual(float(result['loss']), 0.)
                self.assertTrue(bool((value['logits'].grad == 0).all()))
                self.assertEqual(float(value['null'].grad), 0.)


if __name__ == '__main__':
    unittest.main()
