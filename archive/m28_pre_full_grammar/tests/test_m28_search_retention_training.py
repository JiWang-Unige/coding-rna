import sys
import unittest
from pathlib import Path
import numpy as np
import torch
from src.m28.query_pointer import GeneQueryPointer
from src.m28.query_prefix_value import START,NEXT_START
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from query_search_retention_training import search_retention_loss
from trace_frozen_query_beam import original_query


def case():
    sequence = 'ATGAAATAA'*5
    torch.manual_seed(0)
    model = GeneQueryPointer(num_queries=2,hidden=8)
    w = model.for_window(torch.randn(len(sequence),8))
    natural = original_query(w,sequence,0)
    complete = {tuple(tuple(e) for e in r['chain']) for r in natural['completed']}
    reference = next(((a,a+9),) for a in range(0,45,9) if ((a,a+9),) not in complete)
    targets = dict(positive_chains=[reference],protected_chains=[],unknown_chain_spans=[],
                   callable_bases=np.ones(len(sequence),dtype=bool))
    return sequence,w,reference,targets,natural


def cutoff_reference(node):
    if node['event'] in (START,NEXT_START):
        return tuple(node['chain'])+((node['entry'],node['entry']+9),)
    return tuple(node['chain'])


class SearchRetentionTrainingTests(unittest.TestCase):
    def test_current_search_is_unchanged_and_qualified_loss_has_live_gradients(self):
        sequence,w,reference,targets,natural = case()
        out = search_retention_loss(w,sequence,{0:reference},targets,1.)
        self.assertEqual(original_query(w,sequence,0),natural)
        self.assertEqual(out['positive_queries'],1)
        self.assertEqual(out['reference_injected_candidates'],0)
        self.assertEqual(out['records'][0]['cutoff_annotation']['label'],0)
        self.assertGreater(float(out['loss'].detach()),0.)
        out['loss'].backward()
        self.assertTrue(bool(torch.isfinite(w.model.pointer_query.weight.grad).all()))
        self.assertGreater(float(w.model.pointer_query.weight.grad.norm()),0.)

    def test_other_primary_or_protected_current_competitor_is_not_trained_negative(self):
        for key,label in (('positive_chains',1),('protected_chains',-1)):
            sequence,w,reference,targets,_ = case()
            observed = search_retention_loss(w,sequence,{0:reference},targets,1.)
            other = cutoff_reference(observed['records'][0]['cutoff_node'])
            targets[key].append(other)
            safe = search_retention_loss(w,sequence,{0:reference},targets,1.)
            self.assertEqual(safe['records'][0]['cutoff_annotation']['label'],label)
            self.assertEqual(float(safe['loss'].detach()),0.)
            safe['loss'].backward()
            self.assertEqual(float(w.model.pointer_query.weight.grad.norm()),0.)

    def test_raw_unsupported_reference_remains_original_teacher_supervision(self):
        sequence = 'ATGCGTAGTTTAA'
        torch.manual_seed(0)
        w = GeneQueryPointer(num_queries=2,hidden=8).for_window(torch.randn(len(sequence),8))
        reference = ((0,4),(8,13))
        targets = dict(positive_chains=[reference],protected_chains=[],unknown_chain_spans=[],
                       callable_bases=np.ones(len(sequence),dtype=bool))
        out = search_retention_loss(w,sequence,{0:reference},targets,1.)
        self.assertEqual(out['positive_queries'],1)
        self.assertFalse(out['records'][0]['raw_grammar_supported'])
        self.assertEqual(out['records'][0]['unsupported_reason'],'next_start')
        self.assertEqual(out['records'][0]['supervision_disposition'],'original_teacher_only')
        self.assertEqual(float(out['loss'].detach()),0.)


if __name__ == '__main__':
    unittest.main()
