import sys
import unittest
from pathlib import Path
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.m28.candidates import generate,Budget
from src.m28.prefix_search import generate_prefix
from src.m28.prefix_value import PrefixValue
from src.m28.prefix_training import FirstCutObserver,prefix_targets,prefix_losses

def fixture():
    dna='ATGATGATGAAA'+'GT'+'C'*7+'AG'+'AAATAA'
    ep=np.zeros((len(dna),4));ep[3,0]=3;ep[6,0]=6
    chain=[(0,12),(23,29)]
    targets={'positive_chains':[chain],'protected_chains':[chain],
             'unknown_chain_spans':[],'callable_bases':np.ones(len(dna),dtype=bool)}
    return dna,ep,chain,targets

class PrefixTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(0);torch.set_num_threads(2)

    def test_paired_additive_control_is_literal_original(self):
        dna=('ATGAAATAA'+'GTCCCCAG')*16
        rng=np.random.default_rng(7)
        b=Budget(bin_bp=64,sites_per_bin=2,ends_per_start_kind=2,beam_per_carry=1,chains_per_kb=1)
        for _ in range(3):
            ep=rng.normal(size=(len(dna),4));cp=rng.uniform(.1,.9,len(dna))
            links=lambda pairs:np.array([.01*(a-d) for d,a in pairs])
            self.assertEqual(generate(dna,ep,cp,links,b),generate_prefix(dna,ep,cp,links,b))

    def test_prefix_priority_can_change_pruning_without_changing_budget(self):
        dna,ep,chain,_=fixture();b=Budget(beam_per_carry=1)
        old=generate(dna,ep,np.full(len(dna),.6),budget=b)
        scorer=lambda stage,pos,states:[-(s[1][0][0] if s[1] else pos) for s in states]
        new=generate_prefix(dna,ep,np.full(len(dna),.6),budget=b,rank_states_fn=scorer)
        self.assertNotIn(chain,old['chains']);self.assertIn(chain,new['chains'])
        self.assertEqual(old['budget'],new['budget']);self.assertFalse(new['reference_used'])

    def test_prefix_cache_and_gradients_use_frozen_features(self):
        model=PrefixValue(8)
        features=torch.randn(40,8,requires_grad=True)
        states=[(0.,((0,12),),''),(0.,((0,12),(23,29)),'')]
        cached=model.for_window(features)
        cached.score('donor',12,states[:1])
        a=cached.score('final',None,states)
        b=model.for_window(features).score('final',None,states)
        torch.testing.assert_close(a,b)
        a.sum().backward()
        self.assertIsNone(features.grad)
        for layer in (model.exon,model.sequence,model.stage,model.value):
            self.assertGreater(sum(float(p.grad.abs().sum()) for p in layer.parameters() if p.grad is not None),0.)

    def test_other_isoform_and_unknown_prefixes_are_protected(self):
        _,_,chain,t=fixture()
        alt=[(0,12),(26,29)];t['protected_chains'].append(alt)
        states=[(0.,((0,12),),''),(0.,((3,12),),'')]
        self.assertEqual(prefix_targets(states,'entry',23,t),[1,0])
        self.assertEqual(prefix_targets(states,'entry',26,t),[-1,0])
        t['unknown_chain_spans']=[(17,27)]
        self.assertEqual(prefix_targets(states,'entry',23,t),[1,-1])
        self.assertEqual(prefix_targets([(0.,tuple(alt),'')],'final',None,t),[-1])

    def test_first_cut_loss_does_not_inject_into_free_search(self):
        dna,ep,chain,t=fixture();b=Budget(beam_per_carry=1)
        observer=FirstCutObserver([chain])
        original=generate(dna,ep,np.full(len(dna),.6),budget=b)
        free=generate_prefix(dna,ep,np.full(len(dna),.6),budget=b,observer=observer)
        self.assertEqual(original,free);self.assertNotIn(chain,free['chains'])
        self.assertEqual(len(observer.cuts),1);self.assertEqual(observer.cuts[0]['stage'],'donor')
        model=PrefixValue(8);features=torch.randn(len(dna),8,requires_grad=True)
        result=prefix_losses(model,features,free,observer,t,window_weight=.7)
        self.assertGreater(result['known_rank_pairs'],0)
        self.assertIn('reference_injected',result['origins'])
        self.assertTrue(torch.isfinite(result['total']))
        result['total'].backward()
        self.assertIsNone(features.grad)
        self.assertGreater(sum(float(p.grad.abs().sum()) for p in model.parameters() if p.grad is not None),0.)

    def test_background_empty_search_has_zero_differentiable_loss(self):
        dna='C'*30;ep=np.zeros((30,4))
        t={'positive_chains':[],'protected_chains':[],'unknown_chain_spans':[],
           'callable_bases':np.ones(30,dtype=bool)}
        observer=FirstCutObserver([])
        free=generate_prefix(dna,ep,np.full(30,.2),observer=observer)
        model=PrefixValue(8);features=torch.randn(30,8,requires_grad=True)
        result=prefix_losses(model,features,free,observer,t)
        self.assertEqual(float(result['total']),0.)
        result['total'].backward()
        self.assertIsNone(features.grad)
        self.assertEqual(sum(float(p.grad.abs().sum()) for p in model.parameters() if p.grad is not None),0.)

if __name__=='__main__': unittest.main()
