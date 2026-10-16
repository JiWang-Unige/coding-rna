import math
import unittest
import torch
from src.m28.core import ChainHead
from src.m28.chain_set import chain_set_loss
from src.m28.junction_chain import JunctionChainHead


class JunctionChainTests(unittest.TestCase):
    def test_single_exon_preserves_original_fresh_seed_branches_and_scores(self):
        torch.manual_seed(0); old = ChainHead(8).double()
        torch.manual_seed(0); new = JunctionChainHead(8).double()
        for name in ('exon', 'sequence', 'additive', 'nonadditive', 'null'):
            self.assertEqual(list(getattr(old, name).state_dict()), list(getattr(new, name).state_dict()))
            for key, value in getattr(old, name).state_dict().items():
                self.assertTrue(torch.equal(value, getattr(new, name).state_dict()[key]))
        features = torch.randn(18, 8, dtype=torch.float64)
        a, b = old(features, [[[2, 11]]]), new(features, [[[2, 11]]])
        for key in a:
            torch.testing.assert_close(a[key], b[key], rtol=0, atol=0)
        self.assertEqual(float(b['junction_additive'][0]), 0.)

    def test_junction_uses_intronic_sides_and_half_open_interval_mean(self):
        head = JunctionChainHead(2).double()
        features = torch.arange(24, dtype=torch.float64).reshape(12, 2)
        observed = []
        handle = head.junction.register_forward_pre_hook(lambda module, args: observed.append(args[0].detach().clone()))
        head(features, [[[1, 4], [8, 11]]]); handle.remove()
        expected = torch.cat([features[4:8].mean(0), features[3], features[4],
                              features[7], features[8], features.new_tensor([math.log1p(4), 0.])])
        torch.testing.assert_close(observed[0], expected[None], rtol=0, atol=0)

    def test_intron_interior_reaches_chain_evidence_not_just_window_null(self):
        torch.manual_seed(0); old = ChainHead(8).double()
        torch.manual_seed(0); new = JunctionChainHead(8).double()
        features = torch.randn(16, 8, dtype=torch.float64, requires_grad=True)
        chain = [[[1, 4], [11, 14]]]
        a, b = old(features, chain), new(features, chain)
        ga = torch.autograd.grad((a['additive']+a['nonadditive']).sum(), features, retain_graph=True)[0]
        gb = torch.autograd.grad((b['additive']+b['nonadditive']).sum(), features)[0]
        self.assertEqual(float(ga[7].abs().sum()), 0.)
        self.assertGreater(float(gb[7].abs().sum()), 0.)
        self.assertTrue(bool(torch.isfinite(gb).all()))

    def test_padding_is_excluded_and_empty_pool_keeps_null_gradient(self):
        head = JunctionChainHead(8).double()
        features = torch.randn(12, 8, dtype=torch.float64)
        padded = torch.cat([features, torch.full((4, 8), 1e9, dtype=torch.float64)])
        chains = [[[1, 4], [8, 11]]]
        torch.testing.assert_close(head(features, chains), head(padded, chains, 12), rtol=0, atol=0)
        empty = head(padded, [], 12)
        self.assertEqual(len(empty['logits']), 0)
        loss = chain_set_loss(empty, [], [])['loss']; loss.backward()
        self.assertEqual(float(loss), 0.)
        self.assertTrue(torch.equal(head.null.weight.grad, torch.zeros_like(head.null.weight)))

    def test_touching_exons_have_zero_intronic_features_and_finite_score(self):
        head = JunctionChainHead(2).double(); features = torch.randn(8, 2, dtype=torch.float64)
        observed = []
        handle = head.junction.register_forward_pre_hook(lambda module, args: observed.append(args[0].detach().clone()))
        score = head(features, [[[1, 4], [4, 7]]]); handle.remove()
        for indices in (slice(0, 2), slice(4, 8)):
            self.assertTrue(torch.equal(observed[0][0, indices], torch.zeros_like(observed[0][0, indices])))
        self.assertTrue(bool(torch.isfinite(score['logits']).all()))
        with self.assertRaises(ValueError):
            head(features, [[[1, 5], [4, 7]]])

    def test_original_unknown_label_has_zero_direct_gain_gradient(self):
        head = JunctionChainHead(8).double(); features = torch.randn(16, 8, dtype=torch.float64)
        chains = [[[1, 4], [11, 14]], [[1, 5], [10, 14]]]
        scores = head(features, chains); scores['logits'].retain_grad()
        chain_set_loss(scores, chains, [1, -1])['loss'].backward()
        self.assertEqual(float(scores['logits'].grad[1]), 0.)
        self.assertNotEqual(float(scores['logits'].grad[0]), 0.)


if __name__ == '__main__':
    unittest.main()
