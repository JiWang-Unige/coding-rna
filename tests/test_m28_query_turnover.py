import copy
import sys
import unittest
from unittest.mock import patch
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_query_turnover as D


class TurnoverTests(unittest.TestCase):
    def test_complete_same_checkpoint_and_no_reference_inference_required(self):
        fit = {'arm':'A2','completed':True,'optimizer_steps':4608,'configuration':{'fixed':True},
               'final_checkpoint':'A2/step_004608.pt'}
        decoded = {**fit,'inference_complete':True,'windows':8690,'reference_used':False,
                   'test_setaria_used':False,'checkpoint':fit['final_checkpoint']}
        evaluated = {'fit':fit,'reference_primary_chains':7728,'test_setaria_used':False,'methods':{'A2':{}}}
        D.require_source('A2',fit,decoded,evaluated)
        for changes in ({'checkpoint':'wrong'},{'windows':8689},{'reference_used':True}):
            with self.assertRaises(ValueError): D.require_source('A2',fit,{**decoded,**changes},evaluated)
        with self.assertRaises(ValueError): D.require_source('A',fit,decoded,evaluated)

    def test_A_to_A2_intersections_and_state_transitions_are_not_net_only(self):
        records = []
        for old,new in ((True,True),(True,False),(False,True),(False,False)):
            row = {'A_state':'selected' if old else 'not_free','A2_state':'selected' if new else 'not_free'}
            row.update({'A_'+s:old for s in D.STAGES})
            row.update({'A2_'+s:new for s in D.STAGES})
            records.append(row)
        report = D.transition_summary(records)
        for c in report['stages'].values():
            self.assertEqual(c,{'both':1,'A2_only':1,'A_only':1,'neither':1})
        self.assertEqual(sum(r['reference_chains'] for r in report['A_to_A2_states']),4)
        self.assertEqual(report['states']['A'],{'selected':2,'not_free':2})

    def test_reference_mapping_retains_unsupported_and_rejects_non_nested_stages(self):
        ref = {'species':'plant','id':'r','gene_id':'g','seqid':'DEV','strand':'-',
               'CDS':[(20,29,'1'),(40,46,'2')],'M28_unsupported_strata_not_excluded':['longer_than_window']}
        chain = ((20,29),(40,46))
        sets = {'A':dict(zip(D.STAGES,[{chain},{chain},{chain}])),
                'A2':dict(zip(D.STAGES,[{chain},{chain},set()]))}
        row = D.reference_row(ref,sets,{'A':{chain:1.},'A2':{chain:0.}})
        self.assertEqual(row['A_state'],'selected')
        self.assertEqual(row['A2_state'],'owner_nonpositive')
        self.assertEqual(row['unsupported_strata'],'longer_than_window')
        self.assertEqual((row['exon_count'],row['span_bp'],row['coding_bp']),(2,26,15))
        bad = copy.deepcopy(sets); bad['A2']['free_unique'] = set()
        with self.assertRaises(ValueError): D.reference_row(ref,bad,{'A':{chain:1.},'A2':{chain:0.}})

    def test_A2_A3_turnover_keeps_both_arm_labels_and_all_four_memberships(self):
        ref = {'species':'plant','id':'r','gene_id':'g','seqid':'DEV','strand':'+','CDS':[(0,9,'0')]}
        chain = ((0,9),)
        with patch.object(D,'ARMS',('A2','A3')):
            records = []
            for old,new in ((True,True),(True,False),(False,True),(False,False)):
                sets = {arm:{stage:({chain} if present else set()) for stage in D.STAGES}
                        for arm,present in (('A2',old),('A3',new))}
                records.append(D.reference_row(ref,sets,{'A2':{chain:1.},'A3':{chain:1.}}))
            report = D.transition_summary(records)
            self.assertEqual(report['reference_chains'],4)
            for counts in report['stages'].values():
                self.assertEqual(counts,{'both':1,'A3_only':1,'A2_only':1,'neither':1})
            self.assertEqual(sum(row['reference_chains'] for row in report['A2_to_A3_states']),4)
            self.assertEqual(report['states']['A3'],{'selected':2,'not_free':2})


if __name__ == '__main__':
    unittest.main()
