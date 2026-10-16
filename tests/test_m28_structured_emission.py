import unittest
import torch
from src.m28.phase_evidence_chain import candidate_roles
from src.m28.structured_emission import trainable_path_scores, structured_emission_loss
from src.m28.local_path_score import local_path_scores


class StructuredEmissionTests(unittest.TestCase):
    def test_full_log_probability_score_and_gradient_identity(self):
        torch.manual_seed(0)
        logits = torch.randn(42, 15, dtype=torch.float64, requires_grad=True)
        chains = [[(0, 6), (12, 18)], [(24, 30)]]
        value = trainable_path_scores(logits, chains, 36)
        logp = logits[:36].log_softmax(-1)
        positions = torch.arange(36)
        direct = torch.stack([(logp[positions, candidate_roles(c, 36)]-logp[:, 0]).sum()
                              for c in chains])
        torch.testing.assert_close(value, direct, rtol=1e-12, atol=1e-10)
        torch.testing.assert_close(value, local_path_scores(logits, chains, 36), rtol=0, atol=0)
        actual, = torch.autograd.grad(value.sum(), logits, retain_graph=True)
        expected, = torch.autograd.grad(direct.sum(), logits)
        torch.testing.assert_close(actual, expected, rtol=0, atol=1e-14)
        self.assertEqual(float(actual[36:].abs().sum()), 0.)

    def test_two_coexisting_genes_have_independent_known_BCE(self):
        logits = torch.zeros(30, 15, dtype=torch.float64, requires_grad=True)
        chains = [[(0, 6)], [(12, 18)]]
        labels = [1, 0]
        result = structured_emission_loss(logits, chains, labels)
        gain = trainable_path_scores(logits, chains)
        direct = torch.nn.functional.binary_cross_entropy_with_logits(
            gain, torch.tensor(labels, dtype=torch.float64))
        torch.testing.assert_close(result['loss'], direct, rtol=0, atol=0)
        grad, = torch.autograd.grad(result['loss'], logits)
        self.assertLess(float(grad[0, 13]), 0.)
        self.assertGreater(float(grad[12, 13]), 0.)

    def test_unknown_competitor_changes_neither_loss_nor_gradient(self):
        torch.manual_seed(1)
        logits = torch.randn(30, 15, dtype=torch.float64, requires_grad=True)
        known = [[(0, 12)], [(0, 6)]]
        with_unknown = known+[[(6, 18)]]
        a = structured_emission_loss(logits, known, [1, 0])
        b = structured_emission_loss(logits, with_unknown, [1, 0, -1])
        torch.testing.assert_close(a['loss'], b['loss'], rtol=0, atol=0)
        ga, = torch.autograd.grad(a['loss'], logits, retain_graph=True)
        gb, = torch.autograd.grad(b['loss'], logits)
        torch.testing.assert_close(ga, gb, rtol=0, atol=0)
        self.assertEqual(b['unknown_excluded'], 1)

    def test_empty_and_all_unknown_keep_zero_gradient(self):
        for chains, labels in (([], []), ([[(0, 6)]], [-1])):
            logits = torch.randn(12, 15, requires_grad=True)
            result = structured_emission_loss(logits, chains, labels)
            self.assertEqual(float(result['loss']), 0.)
            result['loss'].backward()
            self.assertIsNotNone(logits.grad)
            self.assertEqual(float(logits.grad.abs().sum()), 0.)


if __name__ == '__main__':
    unittest.main()
