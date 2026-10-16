"""Only the new A44 loss unit, with independent artificial formulas."""
import unittest
import torch

from src.m28.labels import CDS_COLUMNS, IGNORE
from src.m28.local_structured_emission import position_normalized_local_structured_loss


class PositionNormalizedLossTests(unittest.TestCase):
    def test_valid_base_not_candidate_or_width_denominator_and_gradient(self):
        values = torch.sin(torch.arange(360, dtype=torch.float64)).reshape(24, 15)*.2
        logits = values.clone().requires_grad_(True)
        targets = dict(seg=torch.full((24,), IGNORE), region=torch.full((24,), IGNORE))
        targets['seg'][0] = 0; targets['seg'][3] = 13; targets['region'][6] = 1
        result = position_normalized_local_structured_loss(logits,
            [[(3, 12)], [(3, 15)], [(15, 18)]], [1, 0, -1], targets, 18,
            local_coefficient=1., structured_coefficient=1.)
        ref = values.clone().requires_grad_(True)
        logp = ref.log_softmax(-1)
        local = (-logp[0, 0]-logp[3, 13]
                 -torch.logsumexp(logp[6, list(CDS_COLUMNS)], 0))/3
        gains = []
        for end in (12, 15):
            p = torch.arange(3, end); roles = 1+(p-3)%3
            roles[:3] = 13; roles[-3:] = 14
            gains.append((ref[p, roles]-ref[p, 0]).sum())
        structured = (torch.logsumexp(torch.stack((gains[0]*0, *gains)), 0)-gains[0])/18
        importance = .37
        (result['loss']*importance).backward(); ((local+structured)*importance).backward()
        torch.testing.assert_close(result['local_loss'], local, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(result['structured_loss'], structured, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(result['structured_candidate_mean_loss']*2/18,
                                   structured, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(logits.grad, ref.grad, rtol=1e-12, atol=1e-12)
        self.assertEqual(result['local_supervised_positions'], 3)
        self.assertEqual(result['known_candidates'], 2)
        self.assertEqual(result['structured_normalizer_valid_bases'], 18)
        self.assertEqual(int(torch.count_nonzero(logits.grad[15:])), 0)

    def test_empty_and_all_unknown_partition_preserve_exact_local_control(self):
        targets = dict(seg=torch.full((24,), IGNORE), region=torch.full((24,), IGNORE))
        targets['region'][5] = 1
        for chains, labels, coefficient in (([], [], 0.), ([], [], 1.), ([[(3, 12)]], [-1], 1.)):
            with self.subTest(chains=chains, coefficient=coefficient):
                values = torch.cos(torch.arange(360, dtype=torch.float64)).reshape(24, 15)*.2
                logits = values.clone().requires_grad_(True)
                result = position_normalized_local_structured_loss(logits, chains, labels,
                    targets, 18, local_coefficient=1., structured_coefficient=coefficient)
                ref = values.clone().requires_grad_(True)
                local = -torch.logsumexp(ref.log_softmax(-1)[5, list(CDS_COLUMNS)], 0)
                result['loss'].backward(); local.backward()
                self.assertEqual(float(result['structured_loss']), 0.)
                torch.testing.assert_close(result['loss'], local, rtol=0, atol=0)
                torch.testing.assert_close(logits.grad, ref.grad, rtol=0, atol=0)
                self.assertEqual(int(torch.count_nonzero(logits.grad[18:])), 0)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
