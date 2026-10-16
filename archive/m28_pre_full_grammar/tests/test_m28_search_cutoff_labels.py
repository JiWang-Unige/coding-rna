import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from analyze_search_cutoff_labels_train import classify, node_at_loss
from src.m28.query_prefix_value import START, NEXT_START, DONOR, STOP


def targets():
    return dict(positive_chains=[[(0,6),(12,18)],[(24,33)]],
        protected_chains=[[(0,6),(12,18)],[(24,33)],[(40,46),(52,61)]],
        unknown_chain_spans=[(70,80)],callable_bases=np.ones(100,dtype=bool))


class SearchCutoffLabelTests(unittest.TestCase):
    def test_another_primary_is_set_positive_not_a_query_negative(self):
        r = classify(dict(chain=(),entry=24,event=START),targets(),[(0,6),(12,18)])
        self.assertEqual(r['label'],1)
        self.assertEqual(r['matching_primary_indices'],[1])
        self.assertFalse(r['target_primary_compatible'])

    def test_proper_donor_and_premature_stop_are_not_conflated(self):
        node = dict(chain=((0,6),),entry=None,event=DONOR)
        self.assertEqual(classify(node,targets(),[(0,6),(12,18)])['label'],1)
        self.assertEqual(classify(dict(node,event=STOP),targets(),[(0,6),(12,18)])['label'],0)

    def test_protected_unknown_and_noncallable_remain_unknown(self):
        for node in (dict(chain=((40,46),),entry=52,event=NEXT_START),
                     dict(chain=(),entry=75,event=START)):
            self.assertEqual(classify(node,targets(),[(0,6),(12,18)])['label'],-1)
        t = targets();t['callable_bases'][90]=False
        self.assertEqual(classify(dict(chain=((85,94),),entry=None,event=STOP),t,[(0,6),(12,18)])['label'],-1)

    def test_score_bound_cutoff_is_complete_but_target_is_live_donor(self):
        loss = dict(stage='retained_beam_score_bound',exon=1,
            target=dict(chain=[[0,6]],entry=None,path_log_score=-3.),
            completion_cutoff=dict(chain=[[24,33]],entry=None,path_log_score=-2.))
        self.assertEqual(node_at_loss(loss)['event'],STOP)
        self.assertEqual(node_at_loss(loss,True)['event'],DONOR)
        self.assertEqual(classify(node_at_loss(loss),targets(),[(0,6),(12,18)])['label'],1)


if __name__ == '__main__':
    unittest.main()
