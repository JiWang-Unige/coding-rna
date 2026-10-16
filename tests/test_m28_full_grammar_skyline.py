"""Artificial proofs of the new exact max-only residual skyline."""
import unittest

import numpy as np
import torch

from src.m28.full_grammar_boundary_compiled import _state_table
from src.m28.full_grammar_skyline import _prune, _skyline_groups, decode_full_grammar_skyline as decode
from test_m28_emission_lattice import oracle_gain, direct_gain, prefer_chain


class FullGrammarSkylineTests(unittest.TestCase):
    def test_actual_graph_has_role_preserving_lower_budget_simulation(self):
        # Adjacent counts suffice by transitivity; every higher-count edge
        # must have a matching lower-count edge for each actual DNA symbol.
        for cap in (1, 2, 4, 128):
            states, destinations, roles, _ = _state_table(cap)
            groups, _ = _skyline_groups(cap)
            previous = {}
            comparisons = 0
            for high, (stage, count, carry) in enumerate(states):
                group = groups[high]
                if group in previous:
                    low = previous[group]
                    low_count = states[low][1]
                    self.assertLess(low_count, count)
                    for base in range(5):
                        lower = {(states[target], int(roles[base, low, edge]))
                            for edge, target in enumerate(destinations[base, low]) if target >= 0}
                        for edge, target in enumerate(destinations[base, high]):
                            if target < 0:
                                break
                            next_stage, next_count, next_carry = states[target]
                            expected_count = 0 if next_stage == 'B' else next_count - count + low_count
                            expected = ((next_stage, expected_count, next_carry), int(roles[base, high, edge]))
                            self.assertIn(expected, lower)
                            comparisons += 1
                previous[group] = high
            if cap > 1:
                self.assertGreater(comparisons, 0)

    def test_skyline_discards_equal_worse_counts_but_retains_better_and_actual_carry(self):
        states, _, _, background = _state_table(4)
        indices = {state: i for i, state in enumerate(states)}
        groups, size = _skyline_groups(4)
        entries = [(('E', 1, ''), 4.), (('E', 2, ''), 4.), (('E', 3, ''), 5.),
                   (('E', 4, ''), 3.), (('E', 2, 'A'), -9.)]
        following = np.full(len(states), -np.inf)
        following[background] = 0.
        for state, value in entries:
            following[indices[state]] = value
        candidates = np.array([indices[state] for state, _ in entries] + [background], dtype=np.int32)
        active = np.empty(len(states), dtype=np.int32)
        kept, removed = _prune(following, candidates, len(candidates), groups, np.empty(size), active)
        self.assertEqual(removed, 2)
        self.assertEqual(set(active[:kept]), {background, indices['E', 1, ''], indices['E', 3, ''], indices['E', 2, 'A']})
        self.assertEqual(following[indices['E', 2, '']], -np.inf)

    def test_complete_gene_set_score_and_legal_chains_match_independent_oracle(self):
        generator = torch.Generator().manual_seed(211)
        for sequence in ('ATGAAATAAATGCCCTAA', 'ATGCGTNAGTTTAA', 'ATGGCTAGAAATAA',
                'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA'):
            for cap in (1, 2, 4):
                for values in (torch.randint(-3, 4, (len(sequence), 15), generator=generator).float(),
                        torch.randn(len(sequence), 15, generator=generator).double()):
                    expected, legal = oracle_gain(sequence, values, cap)
                    actual = decode(sequence, values, cap, 3)
                    self.assertAlmostEqual(actual['gain'], expected, places=9)
                    last = 0
                    for chain, gain in zip(actual['chains'], actual['chain_gains']):
                        self.assertIn(chain, legal)
                        self.assertGreaterEqual(chain[0][0], last)
                        self.assertAlmostEqual(gain, direct_gain(values, chain), places=9)
                        last = chain[-1][1]

    def test_checkpoint_spans_preserve_new_ties_and_arrays(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        values = (torch.arange(len(sequence)*15).reshape(-1, 15).double()%11 - 5)/7
        reference = None
        for span in (1, 2, 3, 5, 7, 8192):
            actual = decode(sequence, values, 4, span)
            key = (actual['chains'], actual['gain'], actual['chain_gains'])
            if reference is None:
                reference = key
            self.assertEqual(key, reference)
            counts = actual['counts']
            rows = len([0] + list(range(span, len(sequence), span)) + [len(sequence)])
            self.assertEqual(counts['checkpoint_bytes'], rows*counts['state_count']*8)
            self.assertEqual(counts['maximum_block_predecessor_bytes'], min(span, len(sequence))*counts['state_count']*4)

    def test_zero_null_row_gauge_adjacent_genes_and_sparse_visits(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        values = torch.zeros(18, 15)
        null = decode(sequence, values, 128, 7)
        self.assertEqual((null['chains'], null['gain']), ((), 0.))
        self.assertLess(null['counts']['forward_active_state_visits'], len(sequence)*null['counts']['state_count'])
        values[:, 0] = -1.
        values[:, 13] = 2.
        values[:, 14] = 2.
        first = decode(sequence, values, 128, 7)
        shifted = decode(sequence, values + torch.arange(18)[:, None]*8, 128, 7)
        self.assertEqual((first['chains'], first['gain'], first['chain_gains']),
            (shifted['chains'], shifted['gain'], shifted['chain_gains']))
        self.assertEqual(first['chains'], (((0, 9),), ((9, 18),)))

    def test_128_exon_capacity_higher_budget_winner_and_cap127_exclusion(self):
        sequence = 'GTCAG'.join(['ATG'] + ['AAA']*126 + ['TAA'])
        chain = tuple((8*i, 8*i + 3) for i in range(128))
        values = prefer_chain(sequence, chain).float()
        actual = decode(sequence, values, 128, 31)
        self.assertEqual(actual['chains'], (chain,))
        self.assertEqual(actual['gain'], len(sequence)*10.)
        self.assertNotIn(chain, decode(sequence, values, 127, 31)['chains'])


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
