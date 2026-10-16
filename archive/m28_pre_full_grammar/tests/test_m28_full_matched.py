import copy
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_full_train as R
from src.m28.training import one_hot
from src.m28.labels import build_targets


def dataset():
    windows = {}; catalogs = {s:{} for s in R.ALLOW}; entries = []; draws = []
    for species in sorted(R.ALLOW):
        for seqid in sorted(R.ALLOW[species]):
            for strand in ('+','-'):
                ident = species+'|'+seqid+'|0|'+strand; tid = ident+'tx'
                w = {'id':ident,'species':species,'seqid':seqid,'strand':strand,'split':'train',
                     'start':0,'end':24576,'width':24576,'valid_bases':24576,'pad_right_after_orientation':0,
                     'label_policy':'CDS_assessable_primary_R4','positive_ids':[tid],
                     'overlapping_transcript_ids':[tid],'chain_unknown_intervals':[],
                     'fine_label_unknown_intervals':[],'no_annotated_gene_either_strand':False}
                windows[ident] = w
                catalogs[species][tid] = {'id':tid,'split':'train','seqid':seqid,'strand':strand,
                                          'training_primary':True,'eligible_positive':True,
                                          'CDS':[[300,309,0]],'CDS_complete_supported':True,'gene_id':tid}
                entries.append({'ordinal':len(entries),'window':{k:w[k] for k in R.GEOMETRY},
                                'artifact':'synthetic/'+str(len(entries))+'.pt'})
                draws.extend([{'window_id':ident,'species':species,'stratum':[seqid,strand]}]*128)
    order = [{'window_id':w['id'],'species':w['species'],'stratum':[w['seqid'],w['strand']],
              'importance_weight_original_stratum_objective':1.} for w in windows.values()]
    random.Random(0).shuffle(order)
    return order,entries,windows,catalogs,draws


class FullMatchedTests(unittest.TestCase):
    def test_every_window_once_with_balanced_original_expected_stratum_risk(self):
        data = dataset(); facts = R.validate_dataset(*data)
        self.assertEqual(facts['windows'],12)
        self.assertTrue(all(v['effective_weight_sum']==6. for v in facts['per_species'].values()))
        self.assertEqual(R.STEPS,76536)

    def test_duplicates_missing_reordered_and_wrong_weights_are_rejected(self):
        for kind in ('duplicate','missing','reordered','weight'):
            order,entries,windows,catalogs,draws = dataset()
            if kind=='duplicate': order[0]=order[1]
            elif kind=='missing': order.pop()
            elif kind=='reordered': order.reverse()
            else: order[0]['importance_weight_original_stratum_objective']=1.1
            with self.assertRaises(ValueError): R.validate_dataset(order,entries,windows,catalogs,draws)

    def test_feature_and_R4_target_geometry_must_match(self):
        for kind in ('split','geometry','policy','eligibility'):
            order,entries,windows,catalogs,draws = dataset()
            w = next(iter(windows.values()))
            if kind=='split': w['split']='val'
            elif kind=='geometry': entries[0]['window']['start']=1
            elif kind=='policy': w['label_policy']='old'
            else: catalogs[w['species']][w['positive_ids'][0]]['eligible_positive']=False
            with self.assertRaises(ValueError): R.validate_dataset(order,entries,windows,catalogs,draws)

    def test_mixed_catalog_DEV_poison_is_not_parsed(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'catalog.jsonl'
            path.write_text('{"split": "val", broken DEV\n'+json.dumps({'split':'train','id':'allowed'})+'\n')
            self.assertEqual(list(R.train_records(path)),[{'split':'train','id':'allowed'}])

    def test_native_feature_summary_rejects_partial_or_DEV_cache(self):
        s = {'completed':True,'feature_spec':R.FEATURE_SPEC,'train_unique_windows':25512,
             'native_windows':24021,'reused_existing_windows':1491,'DEV_windows':0,
             'optimizer_steps_this_run':0,'test_setaria_used':False,'DEV_forward_used':False}
        R.validate_feature_summary(s)
        for k,v in (('train_unique_windows',25511),('DEV_windows',1),('optimizer_steps_this_run',1),
                    ('test_setaria_used',True)):
            bad = {**s,k:v}
            with self.assertRaises(ValueError): R.validate_feature_summary(bad)

    def test_shared_initialization_matches_between_fixed_arms(self):
        c0,q = R.new_model('C0','cpu'),R.new_model('Q','cpu')
        a,b = c0.state_dict(),q.state_dict()
        shared = [k for k in a if k.startswith(('shared.','segmentation.'))]
        self.assertTrue(shared); self.assertTrue(all(torch.equal(a[k],b[k]) for k in shared))
        self.assertEqual(sum(p.numel() for p in c0.parameters()),366191)
        self.assertEqual(sum(p.numel() for p in q.parameters()),649457)

    def test_actual_supervision_routes_C0_local_and_Q_known_mean_with_full_weight(self):
        sequence = 'CCC'+'ATGAAATAA'+'C'*24
        w = {'start':0,'end':36,'width':36,'valid_bases':36,'strand':'+',
             'overlapping_transcript_ids':['t'],'positive_ids':['t'],
             'fine_label_unknown_intervals':[],'chain_unknown_intervals':[],
             'no_annotated_gene_either_strand':False}
        tx = {'id':'t','gene_id':'g','CDS':[[3,12,0]],'training_primary':True,'CDS_complete_supported':True}
        targets = build_targets(w,{'t':tx},sequence)
        for arm in ('C0','Q'):
            model = R.new_model(arm,'cpu')
            outputs = model(torch.zeros(1,6,2048),one_hot(sequence,36)[None],[36])
            result = R.supervised_readout(arm,model,outputs,targets,36,2.)
            names = {'local'} if arm=='C0' else {'local','existence_bce','boundary_ce','termination_bce'}
            self.assertEqual(set(result['losses']),names)
            self.assertTrue(torch.allclose(result['total'],2*sum(result['losses'].values())))
            if arm=='Q': self.assertEqual(result['query']['query_labels'].count(1),1)
            result['total'].backward()
            grad = model.shared.project.weight.grad
            self.assertIsNotNone(grad); self.assertTrue(torch.isfinite(grad).all())
            if arm=='Q': self.assertIsNotNone(model.query.exists.weight.grad)


if __name__=='__main__':
    torch.set_num_threads(1)
    unittest.main()
