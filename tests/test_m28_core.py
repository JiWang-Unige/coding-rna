import sys
import unittest
from pathlib import Path
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.core import M28Core,ChainHead,SharedContext,chain_loss,select_nonoverlapping

class CoreTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(17)
        torch.set_num_threads(2)

    def test_cross_block_context_and_padding(self):
        shared=SharedContext(8,16).eval()
        tokens=torch.randn(1,16,8)
        dna=torch.nn.functional.one_hot(torch.zeros((1,96),dtype=torch.long),5).float()
        h,_=shared(tokens,dna,[72])
        changed=tokens.clone();changed[:,12:]=1000
        padded=dna.clone();padded[:,72:]=1000
        same,mask=shared(changed,padded,[72])
        torch.testing.assert_close(h,same)
        self.assertEqual(int(mask.sum()),72)
        changed=tokens.clone();changed[:,6]+=20
        distant,_=shared(changed,dna,[72])
        self.assertGreater(float((h[0,30]-distant[0,30]).abs().max()),1e-6)

    def test_shared_initialization_and_shapes(self):
        torch.manual_seed(3);c=M28Core("C0",8,16)
        torch.manual_seed(3);b=M28Core("B1",8,16)
        for x,y in zip(c.shared.parameters(),b.shared.parameters()): torch.testing.assert_close(x,y)
        tokens=torch.randn(2,16,8);dna=torch.randn(2,96,5)
        out=b(tokens,dna,[96,90])
        self.assertEqual(tuple(out["segmentation_logits"].shape),(2,96,15))
        self.assertEqual(tuple(out["endpoint_logits"].shape),(2,96,4))

    def test_chain_and_null_receive_gradients(self):
        head=ChainHead(16)
        h=torch.randn(96,16,requires_grad=True)
        chains=[[(3,12),(30,39)],[(3,12)],[(60,78)]]
        scores=head(h,chains)
        torch.testing.assert_close(scores["additive_only_logits"],scores["additive"]-scores["null"])
        loss=chain_loss(scores,[1,0,-1])
        loss.backward()
        self.assertGreater(float(head.null.weight.grad.abs().sum()),0.)
        self.assertGreater(float(head.sequence.weight_hh_l0.grad.abs().sum()),0.)
        self.assertGreater(float(h.grad.abs().sum()),0.)
        self.assertTrue(torch.isfinite(loss))

    def test_reference_negative_trains_null_competitor(self):
        head=ChainHead(16);h=torch.randn(96,16)
        out=head(h,[[(3,12)]])
        loss=chain_loss(out,[0]);loss.backward()
        self.assertLess(float(head.null.bias.grad),0.) # gradient descent raises null
        self.assertGreater(float(head.additive.bias.grad),0.)

    def test_nonadditive_changes_and_empty_candidates(self):
        head=ChainHead(16);h=torch.randn(96,16)
        out=head(h,[[(3,12),(30,39)],[(3,12),(30,45)]])
        self.assertNotEqual(float(out["nonadditive"][0]),float(out["nonadditive"][1]))
        empty=head(h,[])
        self.assertEqual(empty["logits"].numel(),0)
        self.assertEqual(float(chain_loss(empty,[])),0.)
        with self.assertRaises(ValueError): head(h,[[(30,40),(20,25)]])
        self.assertEqual(tuple(head.link_logits(h,[(12,30)]).shape),(1,))

    def test_multiple_genes_and_null_selector(self):
        chains=[[(0,9)],[(6,18)],[(20,29)]]
        self.assertEqual(select_nonoverlapping(chains,[2,1,3]),[0,2])
        self.assertEqual(select_nonoverlapping(chains,[-1,-2,-3]),[])
        self.assertEqual(select_nonoverlapping(chains,[0,0,0]),[])

if __name__=="__main__": unittest.main()
