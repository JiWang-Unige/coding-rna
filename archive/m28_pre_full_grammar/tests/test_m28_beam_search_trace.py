import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from trace_frozen_query_beam import trace_query, original_query, reconcile
from test_m28_query_decode import fixture


class BeamSearchTraceTests(unittest.TestCase):
    def test_exact_native_completion_and_carry_are_unchanged(self):
        for sequence, chain in [('ATGCGTTAGTTTAA', ((0,4),(9,14))),
                                ('ATGCGCTAGTTTAA', ((0,4),(9,14))),
                                ('ATGAAATAA', ((0,9),))]:
            model, features = fixture(sequence, chain)
            window = model.for_window(features)
            for q in range(2):
                observed, refs = trace_query(window, sequence, q, {0:chain})
                reconcile(original_query(window, sequence, q), observed)
                self.assertTrue(refs[0]['returned'])
                self.assertIsNone(refs[0]['first_loss'])

    def test_reference_fifth_start_is_observed_not_injected(self):
        sequence = 'ATGAAATAA'*5; chain = ((36,45),)
        model, features = fixture(sequence, chain)
        for a, score in zip(range(0,45,9), [6,5,4,3,2]):
            features[a,0] = score
        window = model.for_window(features)
        observed, refs = trace_query(window, sequence, 0, {0:chain})
        reconcile(original_query(window, sequence, 0), observed)
        self.assertEqual(refs[0]['first_loss']['stage'], 'local_start_topk')
        self.assertEqual(refs[0]['first_loss']['target_rank'], 5)
        self.assertFalse(refs[0]['returned'])

    def test_terminal_pool_and_donor_pool_losses_remain_distinct(self):
        sequence = 'ATGAAATAA'*2
        model, features = fixture(sequence, [(0,9),(9,18)])
        # Make the first terminal chain unambiguously better despite distinct
        # geometric end normalizers and the fixture's per-exon EOS signals.
        features[0,0] = 10; features[9,0] = 0
        features[8,1] = 5; features[17,1] = 1; features[8,2] = 1
        window = model.for_window(features)
        observed, refs = trace_query(window, sequence, 0, {0:((9,18),)}, completed_limit=1)
        reconcile(original_query(window, sequence, 0, completed_limit=1), observed)
        self.assertEqual(refs[0]['first_loss']['stage'], 'completion_pool_pruned')
        sequence = 'ATG'+'GTG'*10+'GTTAG'+'AAATAA'
        model, features = fixture(sequence, [(0,len(sequence))])
        # A real lower-ranked donor, not an invented grammar event.
        from src.m28.candidates import ORFIndex
        from src.m28.query_decode import _end_events, _legal_ends
        donors = [b for b, kind, _ in _legal_ends(ORFIndex(sequence), _end_events(sequence), 0, '', True)
                  if kind == 'donor']
        for b in donors:
            features[b-1,1] = 2; features[b-1,2] = -1
        chain = ((0,donors[4]), (38,len(sequence)))
        window = model.for_window(features)
        observed, refs = trace_query(window, sequence, 0, {0:chain})
        reconcile(original_query(window, sequence, 0), observed)
        self.assertEqual(refs[0]['first_loss']['stage'], 'continuation_pool_pruned')

    def test_locally_ranked_entry_can_lose_to_other_active_histories(self):
        sequence = 'ATG'+'GTG'*5+'GTTAG'+'AAG'*4+'AAATAA'
        chain = ((0,12),(23,41))
        model, features = fixture(sequence, chain)
        from src.m28.candidates import ORFIndex
        from src.m28.query_decode import _end_events, _legal_ends
        donors = [b for b,kind,_ in _legal_ends(ORFIndex(sequence), _end_events(sequence), 0, '', True)
                  if kind == 'donor']
        for j,b in enumerate(donors):
            features[b-1,1] = 6-j if j < 4 else 0
            features[b-1,2] = -1
        for a,score in zip((23,26,29,32),(3.,2.9,2.8,2.7)):
            features[a,0] = score
        window = model.for_window(features)
        observed, refs = trace_query(window, sequence, 0, {0:chain})
        reconcile(original_query(window, sequence, 0), observed)
        self.assertEqual(refs[0]['first_loss']['stage'], 'global_entry_pool_pruned')

    def test_score_bound_is_a_retained_beam_stop_not_completion_eviction(self):
        sequence = 'ATGCGTTAGTTTAA'+'ATGAAATAA'
        chain = ((0,4),(9,14))
        model, features = fixture(sequence, chain)
        features[14,0] = 6; features[22,1] = 6; features[22,2] = 1
        window = model.for_window(features)
        observed, refs = trace_query(window, sequence, 0, {0:chain}, completed_limit=1)
        reconcile(original_query(window, sequence, 0, completed_limit=1), observed)
        self.assertGreater(observed['counts']['score_bound_stopped_live'], 0)
        self.assertEqual(refs[0]['first_loss']['stage'], 'retained_beam_score_bound')


if __name__ == '__main__':
    unittest.main()
