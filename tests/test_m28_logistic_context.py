import unittest
import torch
from src.m28.core import chain_loss
from src.m28.chain_set import chain_set_loss
from src.m28.context_chain import ContextChain


class LogisticContextTests(unittest.TestCase):
    def test_mean_known_BCE_derivative_and_conflicting_unknown_protection(self):
        gain = torch.tensor([-.4, .7, 12.], dtype=torch.float64, requires_grad=True)
        score = dict(logits=gain, null=gain.new_tensor(0., requires_grad=True))
        loss = chain_loss(score, [1,0,-1]); loss.backward()
        expected = torch.cat(((gain.detach()[:2].sigmoid()-gain.new_tensor([1,0]))/2,
                              gain.new_zeros(1)))
        torch.testing.assert_close(gain.grad, expected, rtol=1e-12, atol=1e-12)
        self.assertIsNone(score['null'].grad)  # Null already affects actual model logits.
        self.assertLess(float(gain.grad[0]), 0.); self.assertGreater(float(gain.grad[1]), 0.)

    def test_conflict_does_not_remove_pointwise_null_penalty(self):
        chains = [[[0,30]], [[0,15]], [[15,30]]]
        gain = torch.tensor([14., 3., 3.], dtype=torch.float64, requires_grad=True)
        score = dict(logits=gain, null=gain.new_tensor(0.))
        local = chain_loss(score, [1,0,0])
        structured = chain_set_loss(score, chains, [1,0,0])['loss']
        point = torch.autograd.grad(local, gain, retain_graph=True)[0]
        joint = torch.autograd.grad(structured, gain)[0]
        self.assertGreater(float(point[1]), .3)
        self.assertLess(float(joint[1]), .001)
        self.assertGreater(float(local), float(structured))

    def test_BCE_unknown_does_not_change_joint_context_gradients(self):
        torch.manual_seed(0); model = ContextChain(8,16)
        features, dna = torch.randn(1,8,8), torch.randn(1,48,5)
        chains = [[[0,6],[12,18]], [[0,6],[24,30]], [[6,12],[18,24]]]
        before = {k: v.detach().clone() for k,v in model.state_dict().items()}; snapshots = []
        for subset, labels in ((chains[:2],[1,0]), (chains,[1,0,-1])):
            model.zero_grad(set_to_none=True)
            score = model(features,dna,[30],subset); score['logits'].retain_grad()
            loss = chain_loss(score,labels); loss.backward()
            if len(labels) == 3:
                self.assertEqual(float(score['logits'].grad[2]), 0.)
            for module in list(model.shared.children())+list(model.head.children()):
                params = list(module.parameters())
                self.assertTrue(all(p.grad is not None and bool(torch.isfinite(p.grad).all()) for p in params))
                self.assertGreater(sum(float(p.grad.double().square().sum()) for p in params),0.)
            snapshots.append((loss.detach(), {k:p.grad.detach().clone() for k,p in model.named_parameters()}))
        torch.testing.assert_close(snapshots[0],snapshots[1],rtol=1e-6,atol=1e-6)
        self.assertTrue(all(torch.equal(v,before[k]) for k,v in model.state_dict().items()))

    def test_empty_all_unknown_and_no_gold_preserve_mean_semantics(self):
        for values, labels in (([],[]), ([7.,-3.],[-1,-1]), ([.4,-.9],[0,0])):
            gain = torch.tensor(values,requires_grad=True); null = torch.tensor(.2,requires_grad=True)
            loss = chain_loss(dict(logits=gain,null=null),labels); loss.backward()
            if any(y == 0 for y in labels):
                torch.testing.assert_close(gain.grad,gain.detach().sigmoid()/2)
            else:
                self.assertEqual(float(loss),0.); self.assertTrue(bool((gain.grad == 0).all()))
                self.assertEqual(float(null.grad),0.)


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
