"""A33 role precedence, admissible retained-history bound and raw grammar."""
import unittest
import torch
from src.m28.emission_query_beam import EmissionScores, decode_emission_query_beam, augment_pool
from src.m28.local_path_score import local_path_scores
from src.m28.phase_evidence_chain import candidate_roles
from test_m28_query_decode import fixture


def score_chain(scorer, chain):
    row = dict(chain=(), coding_bases=0, previous_end=-1, emission_gain=0.)
    for index, (a, b) in enumerate(chain):
        gain = scorer.extend(row, a, b, index == len(chain)-1)
        row = dict(chain=row['chain']+((a, b),), coding_bases=row['coding_bases']+b-a,
                   previous_end=b, emission_gain=gain)
    return row


class EmissionQueryBeamTests(unittest.TestCase):
    def test_complete_incremental_roles_match_direct_sum_all_phases_and_short_exons(self):
        torch.manual_seed(0)
        logits = torch.randn(80, 15, dtype=torch.float64)
        chains = [((2, 14),), ((2, 5), (10, 13)), ((2, 6), (11, 22)),
                  ((2, 7), (12, 13), (18, 19), (24, 28), (33, 36)),
                  ((2, 6), (11, 13), (18, 23))]
        scorer = EmissionScores(logits, 80)
        expected = local_path_scores(logits, chains)
        actual = torch.tensor([score_chain(scorer, c)['emission_gain'] for c in chains], dtype=torch.float64)
        torch.testing.assert_close(actual, expected, rtol=0, atol=1e-12)

    def test_positive_future_roles_invalidate_the_old_log_score_bound(self):
        logits = torch.zeros(24, 15, dtype=torch.float64)
        logits[12:, 1:] = 2.
        scorer = EmissionScores(logits, 24)
        row = dict(chain=((0, 4),), coding_bases=4, previous_end=4, emission_gain=-2.)
        self.assertGreater(scorer.upper_bound(row), 0.)
        extended = scorer.extend(row, 12, 17, True)
        self.assertGreater(extended, row['emission_gain'])
        self.assertGreaterEqual(scorer.upper_bound(row), extended)

    def test_raw_gt_gc_carry_cap_and_direct_final_gain(self):
        for intron in ('GTTAG', 'GCTAG', 'GTNAG'):
            sequence = 'ATGC'+intron+'TTTAA'; target = ((0, 4), (9, 14))
            model, features = fixture(sequence, target)
            logits = torch.zeros(len(sequence), 15)
            roles = candidate_roles(target, len(sequence)); logits[torch.arange(len(sequence)), roles] = 4.
            queries = decode_emission_query_beam(model, features, logits, sequence, 3, 4, 4)
            for query in queries:
                self.assertIn(target, [r['chain'] for r in query['completed']])
                for row in query['completed']:
                    self.assertAlmostEqual(row['emission_gain'], float(local_path_scores(logits, [row['chain']])[0]), places=10)
                capped = decode_emission_query_beam(model, features, logits, sequence, 1, 4, 4)[query['query_index']]
                self.assertNotIn(target, [r['chain'] for r in capped['completed']])
        sequence = 'ATGTGTTAGAAAAATAA'; target = ((0, 4), (9, 17))
        model, features = fixture(sequence, target)
        rows = decode_emission_query_beam(model, features, torch.zeros(17, 15), sequence, 3, 4, 4)
        self.assertTrue(all(target not in [r['chain'] for r in q['completed']] for q in rows))

    def test_padding_offset_and_softmask_do_not_change_scores_or_search(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        model, features = fixture(sequence, [(0, 9), (9, 18)])
        logits = torch.arange(18*15).reshape(18, 15).float()/100
        expected = decode_emission_query_beam(model, features, logits, sequence, 3, 4, 4)
        padded_features = torch.cat((features, torch.full((4, 4), 10000.)))
        shifted = logits.double()+torch.arange(18, dtype=torch.float64)[:, None]*8
        padded_logits = torch.cat((shifted, torch.full((4, 15), 1e9, dtype=torch.float64)))
        actual = decode_emission_query_beam(model, padded_features, padded_logits, sequence.lower(), 3, 4, 4, 18)
        for x, y in zip(expected, actual):
            self.assertEqual([r['chain'] for r in x['completed']], [r['chain'] for r in y['completed']])
            for a, b in zip(x['completed'], y['completed']):
                self.assertAlmostEqual(a['emission_gain'], b['emission_gain'], places=10)

    def test_original_pool_is_preserved_and_new_sources_are_not_old_beam(self):
        old = [((0, 9),), ((20, 29),)]
        queries = [dict(completed=[dict(chain=((0, 9),)), dict(chain=((40, 49),))])]
        chains, sources = augment_pool(old, [['greedy'], ['beam4']], queries)
        self.assertEqual(chains, [list(c) for c in (*old, ((40, 49),))])
        self.assertEqual(sources, [['emission_beam4', 'greedy'], ['beam4'], ['emission_beam4']])


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
