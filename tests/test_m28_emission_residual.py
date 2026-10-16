import unittest
import torch
from src.m28.core import select_nonoverlapping
from src.m28.emission_residual import EmissionResidual
from src.m28.local_path_score import local_path_scores
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss
from src.m28.phase_evidence_chain import candidate_roles


class EmissionResidualTest(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(0)
        self.hidden = torch.randn(24, 4)
        self.local = torch.randn(24, 15)
        self.chains = [[(1, 4), (6, 12)], [(13, 19)], [(2, 5), (8, 14)]]
        self.model = EmissionResidual(4)

    def test_zero_start_preserves_original_potential_and_selection(self):
        scores = self.model(self.hidden, self.local, self.chains, 20)
        original = local_path_scores(self.local, self.chains, 20)
        self.assertTrue(torch.equal(scores['logits'], original))
        self.assertTrue(torch.equal(scores['residual_path_gain'], torch.zeros_like(original)))
        self.assertEqual(scores['logits'].dtype, torch.float64)
        self.assertEqual(select_nonoverlapping(self.chains, scores['logits'].tolist()),
                         select_nonoverlapping(self.chains, original.tolist()))
        self.assertEqual(sum(p.numel() for p in EmissionResidual().parameters()), 1806)

    def test_role_correction_direct_arithmetic_and_padding(self):
        with torch.no_grad():
            self.model.margin.weight.copy_(torch.arange(56).reshape(14, 4)/113.)
            self.model.margin.bias.copy_(torch.arange(14)/37.)
        scores = self.model(self.hidden, self.local, self.chains, 20)
        direct = []
        margins = self.model.margin(self.hidden[:20]).double()
        for chain in self.chains:
            roles = candidate_roles(chain, 20)
            direct.append(sum((margins[p, int(roles[p])-1]
                               for p in range(chain[0][0], chain[-1][1])), margins.sum()*0))
        torch.testing.assert_close(scores['residual_path_gain'], torch.stack(direct), rtol=0, atol=1e-12)
        self.assertTrue(bool((scores['residual_local_logits'][:, 0] == 0).all()))
        edited_hidden, edited_local = self.hidden.clone(), self.local.clone()
        edited_hidden[20:] = 1e6; edited_local[20:] = -1e6
        padded = self.model(edited_hidden, edited_local, self.chains, 20)
        self.assertTrue(torch.equal(scores['logits'], padded['logits']))
        expected = torch.autograd.grad(torch.stack(direct).sum(), list(self.model.parameters()), retain_graph=True)
        actual = torch.autograd.grad(scores['logits'].sum(), list(self.model.parameters()))
        for a, b in zip(actual, expected):
            torch.testing.assert_close(a, b, rtol=1e-6, atol=1e-6)

    def test_zero_start_has_trainable_gradient_without_unknown_supervision(self):
        scores = self.model(self.hidden, self.local, self.chains, 20)
        scores['logits'].retain_grad()
        full = null_anchored_chain_set_loss(scores, self.chains, [1, 0, -1])
        known_scores = self.model(self.hidden, self.local, self.chains[:2], 20)
        known = null_anchored_chain_set_loss(known_scores, self.chains[:2], [1, 0])
        torch.testing.assert_close(full['loss'], known['loss'], rtol=0, atol=0)
        full['loss'].backward()
        self.assertEqual(float(scores['logits'].grad[2]), 0.)
        self.assertGreater(float(self.model.margin.weight.grad.norm()), 0.)
        self.assertTrue(all(p.grad is not None and bool(torch.isfinite(p.grad).all())
                            for p in self.model.parameters()))
        self.assertIsNone(self.hidden.grad); self.assertIsNone(self.local.grad)
        self.assertEqual(len(scores['logits']), 3)
        with self.assertRaisesRegex(ValueError, 'Detached'):
            self.model(self.hidden.clone().requires_grad_(), self.local, self.chains, 20)

    def test_empty_and_all_unknown_keep_differentiable_zero(self):
        for chains, labels in (([], []), (self.chains, [-1, -1, -1])):
            self.model.zero_grad(set_to_none=True)
            scores = self.model(self.hidden, self.local, chains, 20)
            loss = null_anchored_chain_set_loss(scores, chains, labels)['loss']
            self.assertEqual(float(loss), 0.)
            loss.backward()
            self.assertTrue(all(p.grad is not None and bool((p.grad == 0).all())
                                for p in self.model.parameters()))

    def test_saturated_correct_frozen_gains_can_have_zero_gradient(self):
        local = torch.zeros_like(self.local)
        for chain, sign in zip(self.chains[:2], (1., -1.)):
            roles = candidate_roles(chain, 20)
            for p in range(chain[0][0], chain[-1][1]):
                local[p, roles[p]] = sign*1000.
        scores = self.model(self.hidden, local, self.chains[:2], 20)
        scores['logits'].retain_grad()
        null_anchored_chain_set_loss(scores, self.chains[:2], [1, 0])['loss'].backward()
        self.assertTrue(bool((scores['logits'].grad == 0).all()))
        self.assertTrue(all(p.grad is not None and bool((p.grad == 0).all())
                            for p in self.model.parameters()))


if __name__ == '__main__':
    unittest.main()
