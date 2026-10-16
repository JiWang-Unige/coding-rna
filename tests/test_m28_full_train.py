import copy
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import prepare_full_train as P
import cache_full_train as C


def window(start=0, strand='+'):
    species,seqid = 'arabidopsis_thaliana','NC_003070.9'
    return {'id':species+'|'+seqid+'|'+str(start)+'|'+strand,
            'species':species,'seqid':seqid,'strand':strand,'split':'train',
            'start':start,'end':start+24576,'width':24576,'valid_bases':24576,
            'pad_right_after_orientation':0,'label_policy':'CDS_assessable_primary_R4','positive_ids':[]}


def inputs():
    windows = [window(0,'+'),window(12288,'+'),window(0,'-'),window(12288,'-')]
    draws = [{'window_id':w['id'],'species':w['species'],'stratum':[w['seqid'],w['strand']],
              'importance_weight_uniform_within_stratum':1.}
             for w in (windows[0],windows[2],windows[2])]
    index = {w['id']:{'artifact':'synthetic/'+str(i)+'.pt'} for i,w in enumerate((windows[0],windows[2]))}
    return windows,index,draws


class FullTrainTests(unittest.TestCase):
    def test_all_windows_once_and_original_expected_stratum_risk(self):
        windows,index,draws = inputs()
        plan,order,sentinels = P.build_plan(windows,index,draws)
        self.assertEqual({r['window']['id'] for r in plan},{w['id'] for w in windows})
        self.assertEqual(sum(r['reused_existing'] for r in plan),2)
        self.assertEqual(len(sentinels),2)
        weights = {r['stratum'][1]:r['importance_weight_original_stratum_objective'] for r in order}
        self.assertAlmostEqual(weights['+'],2/3); self.assertAlmostEqual(weights['-'],4/3)
        self.assertAlmostEqual(sum(r['importance_weight_original_stratum_objective'] for r in order)/len(order),1.)
        self.assertEqual(P.build_plan(windows,index,draws),(plan,order,sentinels))
        self.assertTrue(all(set(r['window']) == set(P.GEOMETRY) for r in plan))

    def test_val_unallowed_duplicate_and_capacity_do_not_enter_full_training(self):
        windows,index,draws = inputs()
        for field,value in (('split','val'),('seqid','sealed'),('positive_ids',[str(i) for i in range(17)])):
            bad = copy.deepcopy(windows); bad[0][field] = value
            with self.assertRaises(ValueError): P.build_plan(bad,index,draws)
        with self.assertRaises(ValueError): P.build_plan(windows+[windows[0]],index,draws)

    def test_reuse_must_match_old_draw_ids_and_new_paths_do_not_collide(self):
        windows,index,draws = inputs()
        with self.assertRaises(ValueError): P.build_plan(windows,{**index,'unknown':{'artifact':'bad'}},draws)
        plan,_,_ = P.build_plan(windows,index,draws)
        self.assertEqual(len({r['artifact'] for r in plan}),len(plan))
        self.assertTrue(all(r['shard'] == r['ordinal'] % P.SHARDS for r in plan))
        self.assertTrue(all(r['artifact'] == index[r['window']['id']]['artifact']
                            for r in plan if r['reused_existing']))

    def test_native_feature_row_keeps_old_raw_oriented_DNA_but_normalizes_encoder_tokens(self):
        w = window(); w.update(end=4,valid_bases=4,pad_right_after_orientation=24572,strand='-')
        with patch.object(C,'encode_window',return_value=torch.zeros(1,4096,2048)) as encode:
            result = C.native_row(None,None,w,{w['species']:{w['seqid']:'ACGR'}})
        self.assertEqual(result['sequence'],'RCGT')
        self.assertEqual(encode.call_args.args[2],'NCGT')
        self.assertEqual(result['features'].dtype,torch.bfloat16)
        self.assertEqual(result['revision'],P.FEATURE_SPEC['revision'])
        self.assertNotIn('positive_ids',result['window'])

    def test_native_replay_requires_exact_same_BF16_features(self):
        old = torch.tensor([[1.,2.]],dtype=torch.bfloat16)
        self.assertEqual(C.exact_replay(old,old.clone()),0.)
        with self.assertRaisesRegex(ValueError,'differ'): C.exact_replay(old,old+1)
        with self.assertRaisesRegex(ValueError,'BF16'): C.exact_replay(old,old.float())


if __name__ == '__main__':
    unittest.main()
