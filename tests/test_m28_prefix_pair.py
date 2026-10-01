import copy
import unittest
import numpy as np
import torch
from src.m28.prefix_value import PrefixValue
from src.m28.prefix_inference import score_paired_pools,validate_additive_replay

class PairedTests(unittest.TestCase):
    def test_common_identity_has_exactly_same_gain_in_different_orders(self):
        torch.manual_seed(0);torch.set_num_threads(2)
        model=PrefixValue(8).eval();h=torch.randn(50,8)
        a=[(0,9)];b=[(12,21),(30,39)];c=[(21,30)]
        early=[a,b,c,a];late=[c,b]
        scores=score_paired_pools(model,h,early,late)
        self.assertEqual(scores['early'][1],scores['late'][1])
        self.assertEqual(scores['early'][2],scores['late'][0])
        self.assertEqual(scores['early'][0],scores['early'][3])
        self.assertTrue(np.isfinite(scores['early']+scores['late']).all())
        self.assertEqual(score_paired_pools(model,h,[],[]),{'early':[],'late':[]})
        self.assertEqual(len(score_paired_pools(model,h,[],[b])['late']),1)

    def test_literal_replay_accepts_json_lists_but_rejects_discrete_changes(self):
        actual={'chains':[[(0,9)],[(12,21),(30,39)]],'links':[(21,30)],
                'counts':{'paths':3},'budget':{'beam':2},'proposal_scores':[2.,1.],'reference_used':False}
        saved=copy.deepcopy(actual)
        saved['chains']=[[[0,9]],[[12,21],[30,39]]];saved['links']=[[21,30]]
        self.assertEqual(validate_additive_replay(actual,saved),0.)
        for key,value in [('chains',saved['chains'][::-1]),('links',[]),('counts',{'paths':4}),
                          ('budget',{'beam':3}),('reference_used',True),('proposal_scores',[2.,float('nan')]),
                          ('proposal_scores',[2.,1.01])]:
            bad=copy.deepcopy(saved);bad[key]=value
            with self.assertRaises(ValueError): validate_additive_replay(actual,bad)
        saved['proposal_scores'][0]+=1e-6
        self.assertLess(validate_additive_replay(actual,saved),1e-5)

if __name__=='__main__': unittest.main()
