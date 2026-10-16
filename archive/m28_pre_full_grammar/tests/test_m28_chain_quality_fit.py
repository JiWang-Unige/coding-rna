import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_chain_quality as Q
import evaluate_chain_quality as E


class QualityFitTests(unittest.TestCase):
    def fit_receipt(self):
        return {'arm':'A6','completed':True,'configuration':Q.CONFIG,'protocol':Q.PROTOCOL,
                'feature_revision':Q.I.F.REVISION,'optimizer_steps':4608,
                'final_checkpoint':str((Q.FIT/'step_004608.pt').relative_to(Q.ROOT)),
                'DEV_used_for_selection':False,'test_setaria_used':False}

    def decoded(self):
        scopes = [{'species':s,'seqid':q,'strand':t,
                   'windows':[[a,min(a+Q.I.M.WIDTH,n)] for a in Q.I.M.starts(n)]}
                  for s,qs in Q.I.DEV_LENGTHS.items() for q,n in qs.items() for t in ('+','-')]
        return {**self.fit_receipt(),'inference_complete':True,'checkpoint':self.fit_receipt()['final_checkpoint'],
                'feature_spec':Q.R.FEATURE_SPEC,'input_mode':'fixed_DEV_raw_FASTA',
                'reference_used':False,'labels_or_feature_cache_used':False,
                'optimizer_steps_this_run':0,'windows':8690,'scopes':scopes}

    def test_only_last_matching_head_not_early_or_dev_selected(self):
        Q.validate_head(self.fit_receipt())
        for key,value in [('optimizer_steps',1536),('DEV_used_for_selection',True),
                          ('final_checkpoint','outputs/other/step_004608.pt')]:
            row = self.fit_receipt(); row[key] = value
            with self.assertRaises(ValueError): Q.validate_head(row)

    def test_all_full_dev_windows_required_and_no_cache(self):
        E.validate_complete(self.fit_receipt(),self.decoded())
        row = self.decoded(); row['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError): E.validate_complete(self.fit_receipt(),row)
        row = self.decoded(); row['labels_or_feature_cache_used'] = True
        with self.assertRaises(ValueError): E.validate_complete(self.fit_receipt(),row)

    def test_raw_record_uses_per_chain_head_logits_not_query_existence(self):
        class Head:
            def __call__(self,features,chains,length):
                return {'logits':torch.tensor([-2.,3.])}
        chains = [[(0,9)],[(15,24)]]; sources = [['greedy'],['beam4']]
        with patch.object(Q,'free_pool',return_value=(chains,sources)):
            row = Q.quality_record(None,Head(),torch.zeros(30,128),'A'*30,0,'window')
        self.assertEqual(row['chains'],chains)
        self.assertEqual(row['score'],[-2.,3.])
        self.assertEqual(row['counts'],{'free_candidates':2,'nonpositive_candidates':1})
        self.assertFalse(row['reference_used'])

    def test_actual_decoder_bridge_preserves_greedy_and_beam_without_labels(self):
        class Generator: query = object()
        greedy = [{'terminated':True,'chain':[(0,9)],'query_logit':-5.},
                  {'terminated':False,'chain':[(8,11)]}]
        beam = [{'completed':[{'chain':[(15,24)],'query_logit':10.}]}]
        with patch.object(Q,'decode_query_chains',return_value=greedy),\
             patch.object(Q,'decode_query_beam',return_value=beam):
            chains,sources = Q.free_pool(Generator(),torch.zeros(30,128),'A'*30)
        self.assertEqual(chains,[[(0,9)],[(15,24)]])
        self.assertEqual(sources,[['greedy'],['beam4']])

    def test_reconstructed_record_rejects_injected_or_duplicate_inference_chains(self):
        row = {'window_index':0,'window_id':'w','chains':[[[0,9]]],'score':[1.],
               'decoder_sources':[['greedy']],'reference_used':False,
               'score_kind':'A6_complete_chain_logit_minus_window_null','quality_logit_threshold':0.}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'scope.jsonl'; path.write_text(json.dumps(row)+'\n')
            scope = {'windows':[[0,30]],'species':'s','seqid':'q','strand':'+'}
            # Existing validator owns window identity; this test targets the new
            # source/geometry checks, not a fabricated full reference evaluation.
            with patch.object(E,'validated_scope_records',return_value=iter([row])):
                self.assertEqual(len(list(E.records(path,scope))),1)
            changed = copy.deepcopy(row); changed['decoder_sources']=[['reference_injected']]
            with patch.object(E,'validated_scope_records',return_value=iter([changed])):
                with self.assertRaises(ValueError): list(E.records(path,scope))
            changed = copy.deepcopy(row)
            changed['chains'] *= 2; changed['score'] *= 2; changed['decoder_sources'] *= 2
            with patch.object(E,'validated_scope_records',return_value=iter([changed])):
                with self.assertRaises(ValueError): list(E.records(path,scope))

    def test_existing_quality_head_only_trains_used_branches(self):
        head = Q.new_quality('cpu')
        self.assertEqual(sum(p.numel() for p in head.parameters() if p.requires_grad),165635)
        self.assertTrue(all(not p.requires_grad for p in head.link.parameters()))


if __name__=='__main__': unittest.main()
