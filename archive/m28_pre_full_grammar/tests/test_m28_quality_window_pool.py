import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import evaluate_quality_window_pool as W
from src.m28.export import as_transcript
from evaluate_pilot import aggregate


class QualityWindowTests(unittest.TestCase):
    def test_A6_aggregation_auxiliary_keeps_primary_ruler_separate(self):
        ref = as_transcript([(0,9)],'q','+','primary')
        old = as_transcript([(0,6)],'q','+','parent')
        annotation = {'genes':{'g':{'start':0,'end':9}}}
        raw = W.A.E.score([ref],[ref],[ref],30,annotation)
        with self.assertRaisesRegex(KeyError,'historical_parent_exact_chain'):
            aggregate({'species':raw})
        measured = W.F.score_species([ref],[ref],[ref],30,annotation,[old])
        result = aggregate({'species':measured})
        self.assertEqual(result['exact_chain']['tp'],1)
        self.assertEqual(result['exact_chain']['reference'],1)
        self.assertEqual(result['historical_parent_exact_chain']['tp'],0)
        self.assertEqual(result['historical_parent_exact_chain']['reference'],1)

    def test_parent_must_be_completed_exact_A6_evaluation_not_other_arm(self):
        fit = {'final_checkpoint':'fixed'}
        decoded = {'checkpoint':'fixed'}
        parent = {'completed':True,'reference_primary_chains':7728,'test_setaria_used':False,
                  'methods':{'A6':{}},'fit':fit,'protocol':W.Q.PROTOCOL,
                  'inference_receipt':{'checkpoint':'fixed'}}
        with patch.object(W.F,'validate_complete'):
            W.require_parent(fit,decoded,parent)
            for changes in ({'completed':False},{'reference_primary_chains':7727},
                            {'methods':{'A3':{}}},{'fit':{'final_checkpoint':'wrong'}}):
                with self.assertRaises(ValueError):
                    W.require_parent(fit,decoded,{**parent,**changes})

    def test_pool_uses_original_per_chain_quality_not_existence_and_preserves_free(self):
        scope = {'species':'plant','seqid':'q','strand':'+','windows':[[0,100],[40,140]],
                 'assembly_counts':{'free_records':2,'free_unique':1}}
        rows = [{'window_index':i,'window_id':'plant|q|'+str(scope['windows'][i][0])+'|+',
                 'chains':[[[70-40*i,79-40*i]]],'score':[2.5 if i==0 else -1.],
                 'reference_used':False,'decoder_sources':[['beam4']],
                 'score_kind':'A6_complete_chain_logit_minus_window_null','quality_logit_threshold':0.}
                for i in range(2)]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'scope.jsonl'
            p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
            assembled = W.assemble_scope(p,scope)
            self.assertEqual(assembled['predictions'],[{'chain':((70,79),),'score':2.5,'window_index':0}])
            self.assertEqual(assembled['counts']['free_unique'],1)
            with self.assertRaisesRegex(ValueError,'free pool'):
                W.assemble_scope(p,{**scope,'assembly_counts':{'free_records':2,'free_unique':2}})
            rows[0]['score_kind'] = 'query_existence_logit_minus_threshold'
            p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
            with self.assertRaisesRegex(ValueError,'score/pool/source'):
                W.assemble_scope(p,scope)

    def test_paired_reference_membership_retains_all_four_classes_and_unsupported(self):
        refs = [as_transcript([(i*12,i*12+9)],'q','+','r'+str(i)) for i in range(4)]
        refs[3]['M28_unsupported_strata_not_excluded'] = ['longer_than_window']
        before = {W.A.E.S.key(refs[i]) for i in (0,1)}
        after = {W.A.E.S.key(refs[i]) for i in (0,2)}
        rows = W.paired_references({'plant':refs},before,after)
        self.assertEqual(len(rows),4)
        self.assertEqual(W.membership(rows),{'both':1,'A7_only':1,'A6_only':1,'neither':1})
        self.assertEqual(rows[3]['unsupported_strata'],'longer_than_window')


if __name__=='__main__':
    unittest.main()
