import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from run_query_prefix_value import candidate_readout_report


class QueryPrefixReportTests(unittest.TestCase):
    def test_query_duplicates_are_not_genes_and_zero_is_not_a_positive_value(self):
        primary=((0,9),);long=((20,23),(30,33),(40,43),(50,53),(60,63))
        protected=((90,99),);negative=((100,109),)
        def row(q,c,v):return dict(query_index=q,chain=c,shared_complete_value=v)
        arms=dict(early=[row(0,primary,.2),row(1,primary,.4),row(2,long,0.),row(3,negative,.8),row(4,protected,-.5)],
                  late=[row(0,primary,.2),row(0,primary,.2),row(4,protected,-.5)])
        targets=dict(positive_chains=[primary,long],protected_chains=[primary,long,protected],
                     unknown_chain_spans=[],callable_bases=np.ones(120,dtype=bool))
        report,gains,equal=candidate_readout_report(arms,targets)
        self.assertTrue(equal);self.assertEqual(gains['early'][primary],.4)
        self.assertEqual(report['early']['unique_geometry'],4)
        self.assertEqual(report['early']['unique_query_geometry'],5)
        self.assertEqual(report['late']['completed_records'],3)
        self.assertEqual(report['late']['unique_geometry'],2)
        self.assertEqual(report['early']['positive_value_geometry'],2)
        self.assertEqual(report['early']['reference_relative_labels'],{'1':2,'0':1,'-1':1})
        self.assertEqual(report['early']['positive_value_reference_relative_labels'],{'1':1,'0':1,'-1':0})
        self.assertEqual(report['late']['positive_value_reference_relative_labels'],{'1':1,'0':0,'-1':0})

    def test_unknown_and_noncallable_spans_remain_unknown_after_slot_max(self):
        targets=dict(positive_chains=[],protected_chains=[],unknown_chain_spans=[(12,15)],
                     callable_bases=np.ones(30,dtype=bool))
        targets['callable_bases'][25]=False
        rows=[dict(query_index=0,chain=((10,19),),shared_complete_value=1.),
              dict(query_index=1,chain=((20,29),),shared_complete_value=1.)]
        report,_,equal=candidate_readout_report(dict(early=rows,late=rows),targets)
        self.assertTrue(equal)
        for arm in report.values():
            self.assertEqual(arm['reference_relative_labels'],{'1':0,'0':0,'-1':2})
            self.assertEqual(arm['positive_value_reference_relative_labels'],{'1':0,'0':0,'-1':2})

    def test_empty_pool_and_shared_query_value_mismatch_are_reported(self):
        targets=dict(positive_chains=[],protected_chains=[],unknown_chain_spans=[],callable_bases=np.ones(20,dtype=bool))
        report,gains,equal=candidate_readout_report(dict(early=[],late=[]),targets)
        self.assertTrue(equal);self.assertEqual(gains,dict(early={},late={}))
        self.assertEqual(report['early']['unique_geometry'],0)
        def row(v):return dict(query_index=0,chain=((0,9),),shared_complete_value=v)
        _,_,equal=candidate_readout_report(dict(early=[row(.5)],late=[row(.4)]),targets)
        self.assertFalse(equal)


if __name__=='__main__':unittest.main()
