"""Only synthetic recovery records, no real inputs or CUDA."""
import copy
import io
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import resume_query_collaboration as R


def row(i=0):
    raw = [{'query_index':q,'query_logit':1.,'chain':[[0,9]],'terminated':True,
            'stop_reason':'canonical_stop'} for q in range(16)]
    return {'window_index':i,'window_id':'window'+str(i),'reference_used':False,
            'existence_logit_threshold':0.,'score_kind':'query_existence_logit_minus_threshold',
            'raw_queries':raw,'chains':[[[0,9]] for q in range(16)],
            'query_indices':list(range(16)),'score':[1.]*16,
            'counts':{'queries':16,'completed_queries':16,'unterminated_queries':0,'completed_nonpositive_gain':0}}


class PrefixTests(unittest.TestCase):
    def read(self,text):
        with patch.object(Path,'exists',return_value=True), patch.object(Path,'open',return_value=io.StringIO(text)):
            return R.read_prefix(Path('synthetic'),['window0','window1'])

    def test_complete_prefix_and_truncated_final_write(self):
        text = json.dumps(row())+'\n'
        saved,incomplete = self.read(text+json.dumps(row(1))[:50])
        self.assertEqual(saved,[row()]); self.assertTrue(incomplete)
        saved,incomplete = self.read(text+json.dumps(row(1))+'\n')
        self.assertEqual(saved,[row(),row(1)]); self.assertFalse(incomplete)

    def test_identity_duplicates_and_complete_malformed_line_stop(self):
        for value in (row(1),{**row(),'reference_used':True},{**row(),'score':[2.]*16}):
            with self.assertRaises(ValueError): self.read(json.dumps(value)+'\n')
        with self.assertRaises(ValueError): self.read(json.dumps(row())+'\n'+json.dumps(row())+'\n')
        with self.assertRaises(json.JSONDecodeError): self.read('broken\n')

    def test_replay_checks_discrete_history_logits_and_positive_operating_point(self):
        self.assertEqual(R.replay_equal(row(),row()),0)
        bad = copy.deepcopy(row()); bad['raw_queries'][0]['chain'] = [[0,6]]
        with self.assertRaises(ValueError): R.replay_equal(row(),bad)
        bad = copy.deepcopy(row()); bad['raw_queries'][0]['query_logit'] += 1e-3
        with self.assertRaises(ValueError): R.replay_equal(row(),bad)
        bad = copy.deepcopy(row()); bad['score'][0] = -1
        with self.assertRaises(ValueError): R.replay_equal(row(),bad)

    def test_reconstructed_counts_do_not_include_sentinel_replays(self):
        from collections import Counter
        c = Counter(); R.count_record(c,row()); R.count_record(c,row(1))
        self.assertEqual(c['queries'],32); self.assertEqual(c['completed_queries'],32)
        self.assertEqual(c['candidate_records'],32)
        self.assertEqual(c['all_queries_positive_completed'],2)
