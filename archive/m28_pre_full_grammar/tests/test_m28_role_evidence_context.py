"""Only new role-to-token wiring, frozen input, and empty/unknown semantics."""
import unittest
import torch
from src.m28.role_evidence_context import RoleEvidenceContext, RoleConditionedChainHead, role_margins
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss


class RoleTokens(unittest.TestCase):
    def test_margin_normalizer_shift_padding_and_candidate_span(self):
        torch.manual_seed(7)
        logits = torch.randn(30,15)
        chain = [[0,6],[12,18]]
        actual = role_margins(logits,chain,24)
        from src.m28.phase_evidence_chain import candidate_roles
        roles = candidate_roles(chain,24)
        logp = logits[:24].double().log_softmax(-1)
        direct = logp[torch.arange(24),roles]-logp[:,0]
        torch.testing.assert_close(actual.double(),direct,rtol=1e-6,atol=1e-6)
        shifted = logits+torch.arange(30)[:,None]
        torch.testing.assert_close(role_margins(shifted,chain,24),actual,rtol=1e-5,atol=2e-6)
        changed = logits.clone(); changed[18:24]=500.; changed[24:]=float('nan')
        torch.testing.assert_close(role_margins(changed,chain,24),actual,rtol=0,atol=0)

    def test_exon_and_junction_columns_follow_candidate_roles(self):
        head = RoleConditionedChainHead(8)
        logits = torch.arange(15,dtype=torch.float32)[None].expand(24,15).clone()
        captured = {}
        hooks = [head.exon[0].register_forward_pre_hook(lambda m,a:captured.update(exon=a[0].detach().clone())),
                 head.junction[0].register_forward_pre_hook(lambda m,a:captured.update(junction=a[0].detach().clone()))]
        try:
            head(torch.zeros(24,8), [[[0,6],[12,18]]], logits,24)
            torch.testing.assert_close(captured['exon'][:,-3:],torch.tensor([[8.5,13.,9.],[9.5,10.,14.]]),rtol=0,atol=0)
            torch.testing.assert_close(captured['junction'][:,-5:],torch.tensor([[6.,9.,6.,6.,10.]]),rtol=0,atol=0)
            head(torch.zeros(24,8), [[[0,6],[6,12]]], logits,24)
            torch.testing.assert_close(captured['junction'][:,-5:],torch.tensor([[0.,9.,0.,0.,10.]]),rtol=0,atol=0)
        finally:
            for hook in hooks: hook.remove()

    def test_joint_gradient_reaches_new_columns_but_not_local_source_or_unknown(self):
        torch.manual_seed(0)
        model = RoleEvidenceContext(8,8)
        state = {k:v.clone() for k,v in model.state_dict().items()}
        features, dna = torch.randn(1,8,8), torch.randn(1,48,5)
        logits = torch.randn(48,15,requires_grad=True)
        chains = [[[0,6],[12,18]],[[0,6]],[[24,30]]]
        score = model(features,dna,[48],chains,logits)
        score['logits'].retain_grad()
        null_anchored_chain_set_loss(score,chains,[1,0,-1])['loss'].backward()
        self.assertEqual(float(score['logits'].grad[-1]),0.)
        self.assertIsNone(logits.grad)
        self.assertGreater(float(model.head.exon[0].weight.grad[:,-3:].abs().sum()),0.)
        self.assertGreater(float(model.head.junction[0].weight.grad[:,-5:].abs().sum()),0.)
        self.assertTrue(all(p.grad is not None and bool(torch.isfinite(p.grad).all()) for p in model.parameters()))
        self.assertTrue(all(torch.equal(v,state[k]) for k,v in model.state_dict().items()))

    def test_parameter_count_and_empty_or_all_unknown_draw(self):
        model = RoleEvidenceContext()
        self.assertEqual(sum(p.numel() for p in model.parameters()),613348)
        self.assertEqual(sum(p.numel() for p in model.head.parameters()),249092)
        small = RoleEvidenceContext(8,8)
        for chains, labels in (([],[]),([[[0,6],[12,18]]],[-1])):
            small.zero_grad(set_to_none=True)
            score=small(torch.randn(1,4,8),torch.randn(1,24,5),[24],chains,torch.randn(24,15))
            loss=null_anchored_chain_set_loss(score,chains,labels)['loss']
            self.assertEqual(float(loss),0.)
            loss.backward()
            self.assertTrue(all(p.grad is None or bool((p.grad==0).all()) for p in small.parameters()))


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
