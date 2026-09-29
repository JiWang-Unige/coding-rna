import copy
import sys
import unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.candidates import Budget,ORFIndex,generate,STOP
from src.m28.core import M28Core
from src.m28.labels import build_targets
from src.m28.training import balanced_binary_loss,connection_targets,inject_training,joint_losses,one_hot
from test_m28_c0_decode import example
from test_m28_labels import fixture

class CandidateTests(unittest.TestCase):
    def test_carry_matches_spliced_codon_scan(self):
        sequence="TAAACGTGACCCATGNTGATTAG"
        idx=ORFIndex(sequence)
        for a in range(len(sequence)):
            for b in range(a+1,len(sequence)+1):
                for carry in ("","A","T","TA","TG","CC"):
                    for terminal in (False,True):
                        if terminal and b-a<3: continue
                        joined=carry+sequence[a:b]
                        codons=[joined[i:i+3] for i in range(0,len(joined)-2,3)]
                        valid=set(joined)<=set("ACGT")
                        valid &= (len(joined)%3==0 and bool(codons) and codons[-1] in STOP and not any(x in STOP for x in codons[:-1])) if terminal else not any(x in STOP for x in codons)
                        expected=joined[-(len(joined)%3):] if len(joined)%3 else ""
                        self.assertEqual(idx.extend(a,b,carry,terminal),expected if valid else None)

    def test_reference_free_graph_all_splice_phases(self):
        for motif in ("GT","GC"):
            for phase in range(3):
                dna,_,expected=example(12+phase,motif)
                ep=np.full((len(dna),4),-8.)
                for chain in expected:
                    ep[chain[0][0],0]=8;ep[chain[-1][1]-1,1]=8
                    for (_,d),(a,_) in zip(chain,chain[1:]):
                        ep[d-1,2]=8;ep[a,3]=8
                result=generate(dna,ep,np.full(len(dna),.8))
                self.assertIn(expected[0],result["chains"])
                for chain in result["chains"]:
                    coding="".join(dna[a:b] for a,b in chain)
                    self.assertEqual(coding[:3],"ATG")
                    self.assertEqual(len(coding)%3,0)
                    self.assertIn(coding[-3:],STOP)
                    self.assertFalse(any(coding[i:i+3] in STOP for i in range(0,len(coding)-3,3)))

    def test_pruning_counts_and_empty_background(self):
        dna=("ATGAAATAA"+"GTCCCCAG")*60
        rng=np.random.default_rng(1)
        b=Budget(bin_bp=64,sites_per_bin=1,ends_per_start_kind=1,intron_links_per_donor=1,beam_per_carry=1,chains_per_kb=1)
        result=generate(dna,rng.normal(size=(len(dna),4)),np.full(len(dna),.6),budget=b)
        self.assertLess(result["counts"]["kept_sites"]["start"],result["counts"]["raw_sites"]["start"])
        self.assertLessEqual(len(result["chains"]),2)
        self.assertLessEqual(result["counts"]["intron_links_retained"],result["counts"]["intron_links_considered"])
        self.assertEqual(generate("C"*60,np.zeros((60,4)),np.full(60,.5))["chains"],[])

class SupervisionTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(0);torch.set_num_threads(2)

    def test_mask_and_class_normalization(self):
        logits=torch.tensor([0.,0.,0.],requires_grad=True)
        loss=balanced_binary_loss(logits,[1,0,-1]);loss.backward()
        self.assertLess(float(logits.grad[0]),0.)
        self.assertGreater(float(logits.grad[1]),0.)
        self.assertEqual(float(logits.grad[2]),0.)
        a=balanced_binary_loss(torch.tensor([.3,-.4]),[1,0])
        b=balanced_binary_loss(torch.tensor([.3]+[-.4]*9),[1]+[0]*9)
        torch.testing.assert_close(a,b)
        z=torch.randn(3,requires_grad=True);balanced_binary_loss(z,[-1]*3).backward()
        self.assertEqual(float(z.grad.abs().sum()),0.)

    def test_injection_is_separate_and_partial_links_are_known(self):
        w,c,dna=fixture(partial=True);t=build_targets(w,c,dna)
        self.assertEqual(connection_targets([(8,14),(8,17)],w,c,t),[1,-1])
        free=[[(2,9)]]
        original=copy.deepcopy(free)
        pool,sources=inject_training(free,[[(3,8),(14,21)]])
        self.assertEqual(free,original)
        self.assertEqual(sources,["free","reference_injected"])
        self.assertEqual(len(pool),2)

    def test_joint_head_gradients_and_masked_candidates(self):
        w,c,dna=fixture();dna=dna[:22]+"N"+dna[23:]
        t=build_targets(w,c,dna)
        model=M28Core("B1",8,16)
        feature=torch.randn(1,5,8)
        out=model(feature,one_hot(dna,30)[None],[24])
        free={"chains":[[(3,8)],[(22,24)]],"links":[(8,14),(8,17)]}
        original=copy.deepcopy(free)
        result=joint_losses(model,out,free,w,c,t)
        self.assertEqual(free,original)
        self.assertEqual(result["chain_labels"],[0,-1,1])
        groups={"local":model.segmentation,"endpoint":model.endpoint,
                "link":model.chain.link,"chain":model.chain.nonadditive}
        for name,layer in groups.items():
            params=list(layer.parameters())+list(model.shared.parameters())
            gradients=torch.autograd.grad(result["losses"][name],params,retain_graph=True,allow_unused=True)
            self.assertGreater(sum(float(g.abs().sum()) for g in gradients[:len(list(layer.parameters()))] if g is not None),0.)
            self.assertGreater(sum(float(g.abs().sum()) for g in gradients[len(list(layer.parameters())):] if g is not None),0.)
        scores=result["chain_scores"]["logits"]
        g=torch.autograd.grad(result["total"],scores,retain_graph=True)[0]
        self.assertGreater(float(g[0]),0.);self.assertEqual(float(g[1]),0.);self.assertLess(float(g[2]),0.)
        eg=torch.autograd.grad(result["total"],out["endpoint_logits"],retain_graph=True)[0][0]
        self.assertEqual(float(eg[~torch.as_tensor(t["endpoint_known"])].abs().sum()),0.)
        self.assertFalse(feature.requires_grad)

if __name__=="__main__": unittest.main()
