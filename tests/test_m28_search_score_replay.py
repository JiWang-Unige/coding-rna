import sys
import unittest
from pathlib import Path
import torch
from src.m28.query_pointer import GeneQueryPointer
from src.m28.query_prefix_value import START, NEXT_START, DONOR, STOP
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from query_search_score_replay import replay_node_score, qualified_retention_margin
from trace_frozen_query_beam import original_query


def window(sequence):
    torch.manual_seed(0)
    model = GeneQueryPointer(num_queries=2, hidden=8)
    return model.for_window(torch.randn(len(sequence),8))


class SearchScoreReplayTests(unittest.TestCase):
    def test_actual_original_completions_have_the_same_gradient_bearing_score(self):
        for sequence in ('ATGAAATAA'*2, 'ATGCGTTAGTTTAA', 'ATGCGCTAGTTTAA'):
            w = window(sequence)
            for q in range(2):
                observed = original_query(w,sequence,q)
                self.assertGreater(len(observed['completed']),0)
                for row in observed['completed']:
                    score = replay_node_score(w,sequence,q,dict(chain=row['chain'],entry=None,event=STOP))
                    self.assertTrue(score.requires_grad)
                    self.assertAlmostEqual(float(score.detach()),row['path_log_score'],places=5)

    def test_donor_and_next_entry_preserve_original_carry_and_lower_bound(self):
        sequence = 'ATGCGTTAGTTTAA'
        w = window(sequence)
        prefix = dict(chain=((0,4),),entry=None,event=DONOR)
        donor = replay_node_score(w,sequence,0,prefix)
        entry = replay_node_score(w,sequence,0,dict(chain=((0,4),),entry=9,event=NEXT_START))
        state = w.advance(w.advance(w.initial[0],'start',0),'end',4)
        expected = torch.log_softmax(w.pointer_logits(state,'next_start',9),dim=-1)[9]
        torch.testing.assert_close(entry-donor,expected)
        with self.assertRaisesRegex(ValueError,'ORF grammar'):
            replay_node_score(w,sequence,0,dict(prefix,event=STOP))

    def test_known_nonmatch_margin_pushes_target_up_and_cutoff_down(self):
        target = torch.tensor(-5.,requires_grad=True)
        cutoff = torch.tensor(-4.,requires_grad=True)
        loss = qualified_retention_margin(target,cutoff,0,1.)
        self.assertEqual(float(loss),2.)
        loss.backward()
        self.assertEqual(float(target.grad),-1.)
        self.assertEqual(float(cutoff.grad),1.)

    def test_other_primary_and_unresolved_cutoffs_have_no_negative_gradient(self):
        for label in (1,-1):
            target = torch.tensor(-5.,requires_grad=True)
            cutoff = torch.tensor(-4.,requires_grad=True)
            loss = qualified_retention_margin(target,cutoff,label,1.)
            loss.backward()
            self.assertEqual(float(loss),0.)
            self.assertEqual(float(target.grad),0.)
            self.assertEqual(float(cutoff.grad),0.)

    def test_replayed_competing_paths_backpropagate_into_the_original_pointer(self):
        sequence = 'ATGAAATAA'*2
        w = window(sequence)
        target = replay_node_score(w,sequence,0,dict(chain=((0,9),),entry=None,event=STOP))
        cutoff = replay_node_score(w,sequence,0,dict(chain=((9,18),),entry=None,event=STOP))
        loss = qualified_retention_margin(target,cutoff,0,100.)
        loss.backward()
        self.assertIsNotNone(w.model.pointer_query.weight.grad)
        self.assertTrue(bool(torch.isfinite(w.model.pointer_query.weight.grad).all()))
        self.assertGreater(float(w.model.pointer_query.weight.grad.norm()),0.)


if __name__ == '__main__':
    unittest.main()
