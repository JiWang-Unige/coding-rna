import unittest
import torch
from src.m28.context_chain import ContextChain
from src.m28.chain_set import chain_set_loss


class ContextChainTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(0)
        self.model = ContextChain(8, 16)
        self.features = torch.randn(1, 8, 8)
        self.dna = torch.randn(1, 48, 5)
        self.chains = [[[0,6], [12,18]], [[0,6], [24,30]], [[6,12], [18,24]]]

    def score(self, chains=None):
        return self.model(self.features, self.dna, [30], self.chains if chains is None else chains)

    def test_joint_forward_preserves_head_composition_order_and_valid_length(self):
        actual = self.score()
        hidden, _ = self.model.shared(self.features, self.dna, [30])
        expected = self.model.head(hidden[0], self.chains, 30)
        for key in actual:
            torch.testing.assert_close(actual[key], expected[key], rtol=0, atol=0)
        reordered = self.score(self.chains[::-1])
        for key in actual:
            torch.testing.assert_close(actual[key], reordered[key] if actual[key].ndim == 0
                                       else reordered[key].flip(0), rtol=0, atol=0)
        tail_features, tail_dna = self.features.clone(), self.dna.clone()
        tail_features[:,5:] = 1000.; tail_dna[:,30:] = 1000.
        padded = self.model(tail_features, tail_dna, [30], self.chains)
        for key in actual:
            torch.testing.assert_close(actual[key], padded[key], rtol=0, atol=0)

    def test_known_chain_loss_reaches_context_DNA_and_all_head_branches(self):
        score = self.score(); score['logits'].retain_grad()
        result = chain_set_loss(score, self.chains, [1,0,-1])
        result['loss'].backward()
        self.assertEqual(float(score['logits'].grad[2]), 0.)
        modules = list(self.model.shared.children())+list(self.model.head.children())
        for module in modules:
            params = list(module.parameters())
            self.assertTrue(all(p.grad is not None and bool(torch.isfinite(p.grad).all()) for p in params))
            self.assertGreater(sum(float(p.grad.double().square().sum()) for p in params), 0.)
        self.assertFalse(self.features.requires_grad)
        self.assertFalse(self.dna.requires_grad)

    def test_unknown_geometry_does_not_change_known_loss_or_joint_gradients(self):
        before = {k: v.detach().clone() for k, v in self.model.state_dict().items()}
        snapshots = []
        for chains, labels in ((self.chains[:2], [1,0]), (self.chains, [1,0,-1])):
            self.model.zero_grad(set_to_none=True)
            result = chain_set_loss(self.score(chains), chains, labels)
            result['loss'].backward()
            snapshots.append((result['loss'].detach().clone(),
                              {k: p.grad.detach().clone() for k, p in self.model.named_parameters()}))
        torch.testing.assert_close(snapshots[0], snapshots[1], rtol=1e-6, atol=1e-6)
        self.assertTrue(all(torch.equal(v, before[k]) for k, v in self.model.state_dict().items()))

    def test_empty_and_all_unknown_preserve_zero_loss_and_no_updates(self):
        for chains, labels in (([], []), (self.chains, [-1]*3)):
            self.model.zero_grad(set_to_none=True)
            score = self.score(chains)
            result = chain_set_loss(score, chains, labels)
            self.assertEqual(float(result['loss']), 0.)
            result['loss'].backward()
            self.assertTrue(all(p.grad is None or bool((p.grad == 0).all()) for p in self.model.parameters()))


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
