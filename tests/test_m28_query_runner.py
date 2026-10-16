import copy
import json
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_query as R
from src.m28.query_model import GeneQueryModel
from src.m28.query_output import query_window_candidates


def drawn_inputs():
    draws,windows,index = [],{},{}
    for species,unique in zip(R.ALLOW,(745,746)):
        seqid = sorted(R.ALLOW[species])[0]
        for j in range(768):
            ident = species+'|'+seqid+'|'+str((j%unique)*12288)+'|+'
            w = {'id':ident,'species':species,'seqid':seqid,'strand':'+','split':'train',
                 'label_policy':'CDS_assessable_primary_R4','positive_ids':[],
                 'overlapping_transcript_ids':[]}
            windows[ident] = w
            index[ident] = {'window_id':ident,'split':'train','species':species}
            draws.append({'window_id':ident,'species':species,'stratum':[seqid,'+'],
                          'importance_weight_uniform_within_stratum':1.})
    return draws,index,windows,{s:{} for s in R.ALLOW}


def dna_fixture():
    species = next(iter(R.ALLOW)); seqid = sorted(R.ALLOW[species])[0]
    sequence = 'C'*256+'ATGAAATAA'+'C'*(600-265)
    w = {'id':species+'|'+seqid+'|0|+','species':species,'seqid':seqid,'strand':'+','split':'train',
         'start':0,'end':600,'width':24576,'valid_bases':600,
         'label_policy':'CDS_assessable_primary_R4','overlapping_transcript_ids':['tx'],
         'positive_ids':['tx'],'fine_label_unknown_intervals':[],'chain_unknown_intervals':[],
         'no_annotated_gene_either_strand':False}
    tx = {'id':'tx','gene_id':'g','seqid':seqid,'strand':'+','split':'train',
          'CDS':[(256,265,'0')],'training_primary':True,
          'eligible_positive':True,'CDS_complete_supported':True}
    obsolete = dict(w,positive_ids=['nonexistent-old-target'],overlapping_transcript_ids=['obsolete'])
    cached = {'model':R.F.MODEL,'revision':R.F.REVISION,'window':obsolete,'sequence':sequence,
              'features':torch.zeros(4096,2048,dtype=torch.bfloat16),
              'labels':{'forbidden_old_label_source':True}}
    return w,{'tx':tx},cached


def timing_fixture():
    draws,windows,catalogs = [],{},{s:{} for s in R.ALLOW}
    for species in R.ALLOW:
        seqid = sorted(R.ALLOW[species])[0]
        definitions = [([],[],[]),(['single'],['single'],[(500,510)]),
                       (['a','b'],['a','b'],[]),(['long'],['long'],[]),
                       (['a','b'],['a','b'],[])]
        for suffix,parts in {'single':[(260,269)],'a':[(260,269)],'b':[(290,299)],
                             'long':[(260,264),(280,285),(300,306)]}.items():
            tid = species+suffix
            catalogs[species][tid] = {'id':tid,'gene_id':tid,'CDS':[(a,b,'0') for a,b in parts]}
        for j,(positives,overlap,unknown) in enumerate(definitions):
            ident = species+'-synthetic-'+str(j)
            windows[ident] = {'id':ident,'species':species,'seqid':seqid,'strand':'+','start':0,'end':600,
                             'width':24576,'valid_bases':600,'positive_ids':[species+p for p in positives],
                             'overlapping_transcript_ids':[species+p for p in overlap],
                             'chain_unknown_intervals':unknown}
            draws.append({'window_id':ident,'species':species})
    return draws,windows,catalogs


class QueryRunnerTests(unittest.TestCase):
    def test_train_join_rejects_split_stratum_weight_and_positive_mismatch(self):
        valid = drawn_inputs(); before = copy.deepcopy(valid[0])
        R.validate_train_inputs(*valid)
        self.assertEqual(valid[0],before)
        for case in ('window_split','index_split','sealed_chromosome','weight','stratum','eligibility','partial_draws'):
            draws,index,windows,catalogs = drawn_inputs()
            ident = draws[0]['window_id']; w = windows[ident]
            if case == 'window_split': w['split'] = 'val'
            if case == 'index_split': index[ident]['split'] = 'val'
            if case == 'sealed_chromosome': w['seqid'] = 'NC_003074.8'
            if case == 'weight': draws[0]['importance_weight_uniform_within_stratum'] = float('nan')
            if case == 'stratum': draws[0]['stratum'][1] = '-'
            if case == 'partial_draws': draws.pop()
            if case == 'eligibility':
                w['positive_ids'] = w['overlapping_transcript_ids'] = ['tx']
                catalogs[w['species']]['tx'] = {'id':'tx','split':'train','seqid':w['seqid'],
                    'strand':w['strand'],'training_primary':True,'eligible_positive':False}
            with self.assertRaises(ValueError):
                R.validate_train_inputs(draws,index,windows,catalogs)

    def test_feature_join_ignores_cached_labels_but_rejects_real_geometry_mismatch(self):
        w,_,cached = dna_fixture()
        R.validate_cached(cached,w)
        for key,value in (('split','val'),('strand','-'),('seqid','wrong'),('start',1),('valid_bases',599)):
            wrong = dict(w); wrong[key] = value
            with self.assertRaisesRegex(ValueError,'geometry'):
                R.validate_cached(cached,wrong)
        for key,value in (('revision','wrong'),('sequence',cached['sequence'][:-1])):
            wrong = dict(cached); wrong[key] = value
            with self.assertRaisesRegex(ValueError,'geometry'):
                R.validate_cached(wrong,w)

    def test_timing_categories_keep_original_ordinal_and_fail_if_missing(self):
        draws,windows,catalogs = timing_fixture()
        selected = R.select_timing_draws(draws,windows,catalogs)
        self.assertEqual(len(selected),8)
        self.assertEqual([p['ordinal'] for p in selected],[2,3,0,1,7,8,5,6])
        self.assertEqual([p['category'] for p in selected[:4]],
                         ['most_positive','most_positive_exons','zero_positive_unprotected','positive_protected'])
        windows[draws[0]['window_id']]['chain_unknown_intervals'] = [(400,410)]
        with self.assertRaisesRegex(ValueError,'Missing fixed TRAIN timing category'):
            R.select_timing_draws(draws,windows,catalogs)

    def test_preparation_identity_ignores_retired_cost_gate(self):
        prepared = {'completed':True,'configuration':copy.deepcopy(R.CONFIG),'optimizer_steps':0,
                    'feature_revision':R.F.REVISION,'windows':8,'budget_gate':{'passed':False}}
        R.validate_preparation(prepared)
        del prepared['budget_gate']
        R.validate_preparation(prepared)
        for key,value in (('completed',False),('windows',7),('optimizer_steps',1),
                          ('feature_revision','wrong')):
            wrong = copy.deepcopy(prepared); wrong[key] = value
            with self.assertRaisesRegex(ValueError,'preparation identity'):
                R.validate_preparation(wrong)
        wrong = copy.deepcopy(prepared); wrong['configuration']['num_queries'] = 11
        with self.assertRaisesRegex(ValueError,'preparation identity'):
            R.validate_preparation(wrong)

    def test_artifact_bytes_reports_without_retired_cost_cap(self):
        file = SimpleNamespace(is_file=lambda:True,
                               stat=lambda:SimpleNamespace(st_size=3*1024**3))
        directory = SimpleNamespace(rglob=lambda pattern:[file])
        with patch.object(R,'STAGE_DIRS',{'fit':directory}):
            self.assertEqual(R.artifact_bytes(),3*1024**3)

    def test_only_exact_last_A_checkpoint_identity_is_accepted(self):
        for state in (False,True):
            valid = {'configuration':copy.deepcopy(R.CONFIG),'arm':'A','completed':True,
                     'feature_revision':R.F.REVISION,'revision':R.F.REVISION,
                     'optimizer_steps':4608,'step':4608}
            R.check_last_identity(valid,state)
            for key,value in (('arm','B1'),('step' if state else 'optimizer_steps',1536),
                              ('revision' if state else 'feature_revision','wrong')):
                wrong = copy.deepcopy(valid); wrong[key] = value
                with self.assertRaisesRegex(ValueError,'last checkpoint'):
                    R.check_last_identity(wrong,state)
            wrong = copy.deepcopy(valid); wrong['configuration']['num_queries'] = 11
            with self.assertRaisesRegex(ValueError,'last checkpoint'):
                R.check_last_identity(wrong,state)

    def test_shared_inputs_use_R4_supervision_not_obsolete_cached_target(self):
        w,catalog,cached = dna_fixture()
        w['chain_unknown_intervals'] = [(400,410)]
        torch.manual_seed(0); model = GeneQueryModel(2048,16,128)
        features,dna = R.window_inputs(cached,w,'cpu')
        outputs = model(features,dna,[600])
        result = R.supervised_readout(model,outputs,w,catalog,cached['sequence'],.25)
        self.assertEqual(result['query']['components']['positive_queries'],1)
        self.assertEqual(result['query']['query_labels'].count(-1),15)
        self.assertEqual(list(result['query']['positive_chains'].values()),[((256,265),)])
        self.assertEqual(result['protection_reasons'],('unknown_chain_span',))
        self.assertTrue(torch.isfinite(result['total']))
        result['total'].backward()  # No optimizer is constructed or stepped.
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None))

    def test_inference_from_shared_context_has_no_label_input_and_preserves_all_queries(self):
        w,_,cached = dna_fixture()
        geometry = {k:w[k] for k in ('id','species','split','seqid','strand','start','end','width','valid_bases')}
        torch.manual_seed(0); model = GeneQueryModel(2048,16,128).eval()
        with torch.no_grad():
            model.query.exists.weight.zero_(); model.query.exists.bias.fill_(1.)
            features,dna = R.window_inputs(cached,geometry,'cpu')
            outputs = model(features,dna,[600])
            free = query_window_candidates(model.query,outputs['features'][0],cached['sequence'],0,w['id'],
                                           R.CONFIG['max_exons'],R.CONFIG['existence_logit_threshold'],600)
        self.assertFalse(free['reference_used'])
        self.assertEqual(free['counts']['queries'],16)
        self.assertEqual(free['counts']['completed_queries'],16)
        self.assertEqual(free['chains'],[((256,265),)]*16)
        self.assertEqual(free['score'],[1.]*16)
        self.assertEqual(len(json.loads(json.dumps(free,allow_nan=False))['raw_queries']),16)


if __name__ == '__main__':
    unittest.main()
