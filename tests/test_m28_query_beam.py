import math
import unittest
import torch
from torch.nn import functional as F
from src.m28.query_beam import decode_query_beam
from src.m28.query_decode import decode_query_chains
from test_m28_query_decode import fixture


def run(model, features, sequence, max_exons=3, valid_bases=None):
    return decode_query_beam(model, features, sequence, max_exons, 4, 4, valid_bases)


class QueryBeamTests(unittest.TestCase):
    def test_second_start_survives_a_dead_greedy_start(self):
        sequence = 'ATGNAACCCATGAAATAA'
        model, features = fixture(sequence, [(9,18)])
        features[0,0], features[9,0] = 6,4
        greedy = decode_query_chains(model, features, sequence, 3)
        self.assertTrue(all(not r['terminated'] for r in greedy))
        rows = run(model, features, sequence)
        self.assertTrue(all(((9,18),) in [r['chain'] for r in q['completed']] for q in rows))
        self.assertTrue(all(((0,18),) not in [r['chain'] for r in q['completed']] for q in rows))

    def test_multiexon_gt_gc_and_actual_codon_carry(self):
        for intron in ('GTTAG','GCTAG','GTNAG'):
            sequence = 'ATGC' + intron + 'TTTAA'
            model, features = fixture(sequence, [(0,4),(9,14)])
            for query in run(model, features, sequence):
                self.assertIn(((0,4),(9,14)), [r['chain'] for r in query['completed']])
                self.assertTrue(all(r['terminated'] for r in query['completed']))

    def test_cap_ambiguity_and_short_intron_never_manufacture_completion(self):
        sequence = 'ATGCGTTAGTTTAA'
        model, features = fixture(sequence, [(0,4),(9,14)])
        # A genuine direct TAG stop at 9 is legal. A cap must not manufacture
        # the intended two-exon continuation or reinterpret donor 4 as STOP.
        for q in run(model, features, sequence, 1):
            chains = [r['chain'] for r in q['completed']]
            self.assertEqual(chains, [((0,9),)])
            self.assertNotIn(((0,4),(9,14)), chains)
        for sequence, chain, limit in [
                ('ATGNAATAA',[(0,9)],3),
                ('ATGCGTAGTTTAA',[(0,4),(8,13)],3),
                ('CCCAAATAA',[(0,9)],3)]:
            model, features = fixture(sequence, chain)
            self.assertTrue(all(not q['completed'] for q in run(model, features, sequence, limit)))

    def test_scores_are_exact_cumulative_boundary_and_eos_not_existence(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        model, features = fixture(sequence, [(0,9),(9,18)])
        rows = run(model, features, sequence, 1)
        window = model.for_window(features)
        expected = {}
        for a,b in [(0,9),(9,18)]:
            state = window.initial[0]
            start = F.log_softmax(window.pointer_logits(state,'start',0),dim=-1)[a]
            state = window.advance(state,'start',a)
            end = F.log_softmax(window.pointer_logits(state,'end',a),dim=-1)[b]
            state = window.advance(state,'end',b)
            eos = F.logsigmoid(model.terminate(state).squeeze(-1))
            expected[((a,b),)] = float(start+end+eos)
        actual = {r['chain']:r['path_log_score'] for r in rows[0]['completed']}
        self.assertEqual(set(actual),set(expected))
        for chain,score in actual.items():
            self.assertAlmostEqual(score,expected[chain],places=5)
            self.assertTrue(math.isfinite(score))
            self.assertLessEqual(score,0)
        self.assertEqual(rows[0]['completed'],sorted(
            rows[0]['completed'],key=lambda r:(-r['path_log_score'],r['chain'])))

    def test_padding_softmask_determinism_and_query_identity(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        model, features = fixture(sequence, [(0,9),(9,18)])
        expected = run(model,features,sequence)
        padded = torch.cat([features,torch.full((4,4),10000.)])
        self.assertEqual(expected,run(model,padded,sequence.lower(),valid_bases=len(sequence)))
        self.assertEqual(expected,run(model,features,sequence))
        self.assertEqual([q['query_index'] for q in expected],[0,1])
        for q in expected:
            self.assertLessEqual(len(q['completed']),4)
            self.assertEqual(len(q['completed']),len({r['chain'] for r in q['completed']}))

    def test_actual_carry_rejects_a_stop_created_at_a_splice(self):
        sequence = 'ATGTGTTAGAAAAATAA'
        model, features = fixture(sequence, [(0,4),(9,17)])
        rows = run(model,features,sequence)
        for q in rows:
            self.assertNotIn(((0,4),(9,17)),[r['chain'] for r in q['completed']])


if __name__ == '__main__':
    unittest.main()
