import math
import unittest
import torch
from src.m28.legal_support import conditional_support


class LegalSupportTests(unittest.TestCase):
    def test_exact_mass_decomposition_and_conditional_rank_invariance(self):
        raw = torch.tensor([.02, .08, .90], dtype=torch.float64).log()
        result = conditional_support(raw[:2], 0)
        self.assertAlmostEqual(float(result['log_mass']), math.log(.10), places=12)
        self.assertAlmostEqual(float(result['conditional_target_logp']), math.log(.2), places=12)
        normalized = raw[:2]-result['log_mass']
        self.assertEqual(int(raw[:2].argmax()), int(normalized.argmax()))
        self.assertAlmostEqual(float(normalized.exp().sum()), 1., places=12)

    def test_joint_boundary_EOS_mass_is_not_boundary_mask_mass(self):
        boundary = torch.tensor([.6, .1], dtype=torch.float64).log()
        branch = torch.tensor([.25, .75], dtype=torch.float64).log()
        result = conditional_support(boundary+branch, 1)
        self.assertAlmostEqual(float(result['log_mass'].exp()), .225, places=12)
        self.assertAlmostEqual(float(result['conditional_target_logp'].exp()), 1/3, places=12)
        self.assertNotAlmostEqual(float(result['log_mass']), float(torch.logsumexp(boundary, 0)))

    def test_conditioning_changes_gradients_not_just_target_score(self):
        logits = torch.tensor([1., 2., 9.], dtype=torch.float64, requires_grad=True)
        raw = torch.log_softmax(logits, 0)
        loss = -conditional_support(raw[:2], 0)['conditional_target_logp']
        actual = torch.autograd.grad(loss, logits)[0]
        expected = torch.tensor([-1., 0., 0.], dtype=torch.float64)
        expected[:2] += torch.softmax(logits[:2].detach(), 0)
        torch.testing.assert_close(actual, expected, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(float(actual[2]), 0., places=12)


if __name__ == '__main__':
    unittest.main()
