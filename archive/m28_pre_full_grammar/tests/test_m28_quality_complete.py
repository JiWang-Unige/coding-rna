import unittest
import torch
from src.m28.quality_complete import decode_query_complete_pool
from src.m28.query_beam import decode_query_beam
from test_m28_query_decode import fixture


def run(model,features,sequence,quality,max_exons=128,limit=4,valid_bases=None):
    return decode_query_complete_pool(model,features,sequence,quality,max_exons,4,limit,valid_bases)


class QualityCompleteTests(unittest.TestCase):
    def test_future_quality_is_not_bounded_by_raw_prefix_and_path_control_is_paired(self):
        sequence,chain='ATGCGTTAGTTTAA',((0,4),(9,14))
        model,features=fixture(sequence,chain)
        features[3,1]=0;features[8,1]=6;features[8,2]=1
        original=decode_query_beam(model,features,sequence,128,4,1)
        self.assertTrue(all(q['counts']['score_bound_stopped_live'] for q in original))
        result=run(model,features,sequence,lambda cs:[float(len(c)) for c in cs],limit=1)
        for q,old in zip(result['queries'],original):
            self.assertEqual(q['quality_completed'][0]['chain'],chain)
            self.assertEqual(q['path_completed'][0]['chain'],old['completed'][0]['chain'])
            self.assertLess(q['quality_completed'][0]['path_log_score'],q['path_completed'][0]['path_log_score'])
        self.assertIn(chain,[r['chain'] for r in result['all_completed']])
        self.assertFalse(result['early_score_bound_used'])
        self.assertFalse(result['quality_used_for_live_history'])

    def test_quality_is_per_geometry_and_negative_scores_are_not_removed(self):
        sequence='ATGAAATAAATGCCCTAA'
        model,features=fixture(sequence,[(0,9),(9,18)])
        observed=[]
        def quality(chains):
            observed.extend(chains);return [-10. for c in chains]
        result=run(model,features,sequence,quality,max_exons=1)
        self.assertEqual(len(observed),len(set(observed)))
        self.assertEqual(observed,sorted(observed))
        for q in result['queries']:
            self.assertEqual(q['quality_completed'],sorted(q['quality_completed'],key=lambda r:r['chain']))
            self.assertTrue(all(r['quality_logit']==-10 for r in q['quality_completed']))

    def test_cap_and_real_carry_do_not_manufacture_completion(self):
        for sequence,chain,cap in [('ATGCGTTAGTTTAA',((0,4),(9,14)),1),
                                   ('ATGTGTTAGAAAAATAA',((0,4),(9,17)),128)]:
            model,features=fixture(sequence,chain)
            result=run(model,features,sequence,lambda cs:[len(c) for c in cs],max_exons=cap)
            self.assertNotIn(chain,[r['chain'] for r in result['all_completed']])

    def test_padding_softmask_and_existence_do_not_select_chains(self):
        sequence='ATGC'+'GTNAG'+'TTTAA'
        model,features=fixture(sequence,[(0,4),(9,14)])
        fn=lambda cs:[len(c) for c in cs]
        expected=run(model,features,sequence,fn)
        padded=torch.cat([features,torch.full((6,4),10000.)])
        self.assertEqual(expected,run(model,padded,sequence.lower(),fn,valid_bases=len(sequence)))
        with torch.no_grad():model.exists.bias.fill_(-100.)
        low=run(model,features,sequence,fn)
        self.assertEqual(expected['all_completed'],low['all_completed'])
        self.assertEqual([q['quality_completed'] for q in expected['queries']],[q['quality_completed'] for q in low['queries']])

    def test_empty_complete_pool_never_invokes_scorer(self):
        sequence='CCCAAATAA';model,features=fixture(sequence,[(0,9)])
        def unexpected(chains):raise AssertionError('No complete candidates to score')
        result=run(model,features,sequence,unexpected)
        self.assertEqual(result['all_completed'],[])

    def test_invalid_score_alignment_and_nonfinite_values_stop(self):
        sequence='ATGAAATAA';model,features=fixture(sequence,[(0,9)])
        for fn in (lambda cs:[],lambda cs:[float('nan') for c in cs]):
            with self.assertRaisesRegex(ValueError,'finite frozen quality'):
                run(model,features,sequence,fn)


if __name__=='__main__':unittest.main()
