import unittest
import torch
from src.m28.query_component_beam import ARMS, completed_score, decode_query_component_pair
from src.m28.query_beam import decode_query_beam
from src.m28.query_training import teacher_chain_components
from test_m28_query_decode import fixture


def run(sequence, chain, cap=3, features=None, model=None):
    if model is None:
        model, features = fixture(sequence, chain)
    return decode_query_component_pair(model, features, sequence, cap, 4, 4)


class QueryComponentBeamTests(unittest.TestCase):
    def test_component_means_and_length_order_are_not_raw_sum(self):
        short = dict(chain=((0, 9),), path_log_score=-3., boundary_log_sum=-2., termination_log_sum=-1.)
        long = dict(chain=((0, 4), (9, 14)), path_log_score=-5., boundary_log_sum=-4., termination_log_sum=-1.)
        self.assertGreater(completed_score(short, 'RAW_SUM'), completed_score(long, 'RAW_SUM'))
        self.assertGreater(completed_score(short, 'COMPONENT_SUM'), completed_score(long, 'COMPONENT_SUM'))
        self.assertLess(completed_score(short, 'COMPONENT_MEAN'), completed_score(long, 'COMPONENT_MEAN'))
        self.assertEqual(completed_score(long, 'COMPONENT_MEAN'), -4/4-1/2)

    def test_single_exon_mean_equals_negative_teacher_components(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        model, features = fixture(sequence, [(0, 9), (9, 18)])
        output = run(sequence, (), cap=1, features=features, model=model)
        original = decode_query_beam(model, features, sequence, 1, 4, 4)
        for q, rows in enumerate(output['RAW_SUM']):
            self.assertEqual([r['chain'] for r in rows['completed']], [r['chain'] for r in original[q]['completed']])
            for a, b in zip(rows['completed'], original[q]['completed']):
                self.assertAlmostEqual(a['path_log_score'], b['path_log_score'], delta=1e-6)
        for rows in output['COMPONENT_MEAN']:
            for row in rows['completed']:
                teacher = teacher_chain_components(model.for_window(features), row['chain'])
                q = row['query_index']
                expected = -float(teacher['boundary_ce'][q]+teacher['termination_bce'][q])
                self.assertAlmostEqual(row['completed_score'], expected, delta=1e-6)

    def test_weighted_pools_share_search_and_no_existence_score(self):
        sequence = 'ATGCGTTAGTTTAA'
        model, features = fixture(sequence, [(0, 4), (9, 14)])
        output = run(sequence, (), features=features, model=model)
        for a, b in zip(output['COMPONENT_SUM'], output['COMPONENT_MEAN']):
            self.assertEqual(a['counts'], b['counts'])
            self.assertEqual(a['counts']['score_bound_stopped_live'], 0)
            for row in b['completed']:
                self.assertAlmostEqual(row['completed_score'],
                    row['boundary_log_sum']/(2*len(row['chain']))+row['termination_log_sum']/len(row['chain']))
        expected = output
        with torch.no_grad():
            model.exists.bias.add_(100.)
        actual = run(sequence, (), features=features, model=model)
        for arm in ARMS:
            for a, b in zip(expected[arm], actual[arm]):
                self.assertEqual([(r['chain'], r['completed_score']) for r in a['completed']],
                                 [(r['chain'], r['completed_score']) for r in b['completed']])

    def test_carry_gc_unknown_and_cap_use_unchanged_dna_grammar(self):
        for intron in ('GTTAG', 'GCTAG', 'GTNAG'):
            sequence = 'ATGC'+intron+'TTTAA'
            output = run(sequence, [(0, 4), (9, 14)])
            for arm in ARMS:
                self.assertTrue(all(((0, 4), (9, 14)) in [r['chain'] for r in q['completed']] for q in output[arm]))
                capped = run(sequence, [(0, 4), (9, 14)], cap=1)
                self.assertTrue(all(((0, 4), (9, 14)) not in [r['chain'] for r in q['completed']] for q in capped[arm]))
        for sequence, chain in [('ATGNAATAA', [(0, 9)]), ('CCCAAATAA', [(0, 9)]),
                                ('ATGTGTTAGAAAAATAA', [(0, 4), (9, 17)])]:
            for rows in run(sequence, chain).values():
                self.assertTrue(all(tuple(chain) not in [r['chain'] for r in q['completed']] for q in rows))

    def test_softmask_padding_ties_and_q_identity_are_deterministic(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        model, features = fixture(sequence, [(0, 9), (9, 18)])
        expected = run(sequence, (), features=features, model=model)
        padded = torch.cat([features, torch.full((4, 4), 10000.)])
        actual = decode_query_component_pair(model, padded, sequence.lower(), 3, 4, 4, len(sequence))
        self.assertEqual(expected, actual)
        self.assertEqual(expected, run(sequence, (), features=features, model=model))
        for rows in expected.values():
            self.assertEqual([q['query_index'] for q in rows], [0, 1])
            for q in rows:
                self.assertLessEqual(len(q['completed']), 4)
                self.assertEqual(len({r['chain'] for r in q['completed']}), len(q['completed']))


if __name__ == '__main__':
    unittest.main()
