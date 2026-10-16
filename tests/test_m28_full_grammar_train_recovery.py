import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import compare_full_grammar_train_recovery as R


def row(name,anchor,candidate,native=(False,False,False),seqid='TRAIN',strand='+'):
    value = dict(id=name,species='plant',seqid=seqid,strand=strand,CDS=[[10,19,'0']],
        CDS_segments=1,span_bp=9,exact_primary_hit=anchor,original_draw_R4_positive_chain=False,
        original_model_error_category='exact_primary_chain' if anchor else 'no_same_strand_CDS_overlap',
        native_exact_hits=dict(zip(R.D.METHODS,native)),model_hits={m:candidate for m in R.METHODS})
    value['model_hits'][R.ANCHOR]=anchor
    return value


class TrainRecoveryTests(unittest.TestCase):
    def test_gains_losses_and_net_reference_hits_preserve_denominator(self):
        rows = [row('a',True,True),row('b',True,False),row('c',False,True),
                row('d',False,True),row('e',False,False),row('f',False,False)]
        result = R.pair(rows,R.LOCAL)
        self.assertEqual(result,dict(reference_chains=6,control_hits=2,candidate_hits=3,
            retained=1,gained=2,lost=1,still_missed=2,net_reference_hits=1,reference_recall=.5))

    def test_native_groups_use_original_not_candidate_hits(self):
        rows = [row('a',True,False,(True,False,False)),row('b',True,False),
                row('c',False,True,(True,True,True)),row('d',False,True)]
        result = R.summarize(rows,R.LOCAL)
        groups = result['fixed_anchor_native_groups']
        self.assertEqual([groups[g]['reference_chains'] for g in R.GROUPS],[1,1,1,1])
        self.assertEqual(groups['anchor_missed_any_native']['candidate_hits'],1)
        self.assertEqual(result['all_three_native_exact_anchor_gap']['candidate_hits'],1)
        self.assertEqual(sum(g['gained'] for g in groups.values()),result['gained'])

    def test_matched_control_differs_from_shorter_update_anchor(self):
        a,b = row('a',False,True),row('b',True,False)
        a['model_hits'][R.CONTROL]=True; b['model_hits'][R.CONTROL]=True
        self.assertEqual(R.pair([a,b],R.LOCAL)['net_reference_hits'],0)
        self.assertEqual(R.pair([a,b],R.LOCAL,R.CONTROL)['net_reference_hits'],-1)

    def test_strata_partition_full_references_without_short_chain_filter(self):
        a,b = row('a',False,True),row('b',True,False)
        b['CDS']=[[10,19,'0'],[8000,8009,'2']]; b['CDS_segments']=2; b['span_bp']=7999
        result = R.describe([a,b],R.LOCAL)
        for groups in result['strata'].values():
            self.assertEqual(sum(v['reference_chains'] for v in groups.values()),2)
            self.assertEqual(sum(v['gained'] for v in groups.values()),1)
            self.assertEqual(sum(v['lost'] for v in groups.values()),1)
        self.assertEqual(result['strata']['span_bp']['over_6144']['lost'],1)

    def test_empty_group_rate_is_null_not_false_zero(self):
        value = R.summarize([],R.LOCAL)
        self.assertIsNone(value['reference_recall'])
        self.assertIsNone(value['all_three_native_exact_anchor_gap']['reference_recall'])
        self.assertEqual(value['net_reference_hits'],0)

    def test_phase_id_contig_and_strand_geometry_identity(self):
        a = row('a',False,False); b = copy.deepcopy(a)
        b['id']='other'; b['CDS'][0][2]='2'
        self.assertEqual(R.D.reference_key(a),R.D.reference_key(b))
        for field,value in (('species','rice'),('seqid','other'),('strand','-')):
            self.assertNotEqual(R.D.reference_key(a),R.D.reference_key({**a,field:value}))

    def test_exact_replay_deduplicates_predictions_and_rejects_count_drift(self):
        a,b = row('a',False,False),row('b',False,False,strand='-')
        pk = {R.D.reference_key(a)}
        official = dict(exact_chain=dict(tp=1,predicted=1,reference=2),prediction_records=2,duplicate_chain_records=1)
        R.assign_hits([a,b],R.LOCAL,pk,2,official)
        self.assertTrue(a['model_hits'][R.LOCAL]); self.assertFalse(b['model_hits'][R.LOCAL])
        for changed in ({**official,'prediction_records':1},
                        {**official,'exact_chain':dict(tp=2,predicted=1,reference=2)}):
            with self.assertRaises(ValueError): R.assign_hits([a,b],R.LOCAL,pk,2,changed)

    def test_all_eight_fixed_candidates_required_no_selected_subset(self):
        a = row('a',False,True)
        self.assertEqual(tuple(R.complete_family([a],R.METHODS)),R.METHODS)
        with self.assertRaises(ValueError): R.complete_family([a],R.METHODS[:-1])
        del a['model_hits'][R.COMBO]
        with self.assertRaises(ValueError): R.complete_family([a],R.METHODS)


if __name__ == '__main__': unittest.main()
