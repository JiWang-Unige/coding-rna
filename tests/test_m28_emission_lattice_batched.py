"""Artificial oracle/identity checks before a same-input full contig replay."""
import math
import unittest

import torch

import profile_emission_lattice_synthetic as S
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.emission_lattice import _Live, _Solution, decode_emission_lattice
from src.m28.emission_lattice_batched import _long_exit_winners, decode_batched_emission_lattice
from test_m28_emission_lattice import oracle_gain, prefer_chain
from test_m28_emission_lattice_owner import reverse_complement, sequence_with_genes


class BatchedEmissionLatticeTests(unittest.TestCase):
    def compare(self, sequence, potentials, cap, oracle=False):
        old = decode_emission_lattice(sequence, potentials, cap)
        actual = decode_batched_emission_lattice(sequence, potentials, cap)
        for key in ('chains', 'gain', 'chain_gains', 'reference_used'):
            self.assertEqual(actual[key], old[key])
        for new, original in (('logical_exon_extensions', 'exon_extensions'),
                ('logical_terminal_extensions', 'completed_extensions'),
                ('logical_cap_nonterminal_extensions', 'exon_cap_nonterminal_edges'),
                ('maximum_active_states', 'maximum_active_states')):
            self.assertEqual(actual['counts'][new], old['counts'][original])
        self.assertLessEqual(actual['counts']['state_update_proposals'], actual['counts']['logical_exon_extensions'])
        if oracle:
            expected, legal = oracle_gain(sequence, potentials, cap)
            self.assertEqual(actual['gain'], expected)
            self.assertTrue(all(chain in legal for chain in actual['chains']))
        return actual

    def test_every_exit_adds_FP64_before_first_arrival_and_matches_scalar(self):
        generator = torch.Generator().manual_seed(61)
        scores = [1., math.nextafter(1., math.inf), -1., math.nextafter(-1., -math.inf), 0., -0.]
        scores += torch.randn(15, generator=generator).double().tolist()
        deltas = [2., -2., 0., -0., 2. ** 53, -(2. ** 53), 2. ** -52]
        deltas += torch.randn(15, generator=generator).double().tolist()
        for values in (scores, list(reversed(scores)), scores[:2]):
            rows = [_Live(score, _Solution(), (), '') for score in values]
            winners, sums = _long_exit_winners(rows, deltas)
            for i, delta in enumerate(deltas):
                expected = max(range(len(rows)), key=lambda j: rows[j].score + delta)
                self.assertEqual(winners[i], expected)
                self.assertEqual(sums[i], rows[expected].score + delta)
        pair = [_Live(score, _Solution(), (), '') for score in scores[:2]]
        self.assertEqual(_long_exit_winners(pair, [2.]), ([0], [3.]))
        # Actual grammar histories A and C converge at acceptor18. The later
        # entry has the larger score; after final delta2 they tie at3. Earlier
        # acceptor9 paths are penalized without penalizing the later ATG.
        sequence = 'ATGAGTTAGATGCGTTAGCCTAA'
        potentials = torch.zeros(len(sequence), 15)
        potentials[3, 7] = 1.
        potentials[9, 13] = 1.
        potentials[12, 7] = 2. ** -52
        potentials[9, 11] = -4.
        potentials[20, 14] = 2.
        result = self.compare(sequence, potentials, 2)
        self.assertEqual(result['gain'], 3.)
        self.assertEqual(result['chains'], (((0, 4), (18, 23)),))

    def test_small_oracle_chain_identity_all_phases_carries_and_null(self):
        generator = torch.Generator().manual_seed(67)
        sequences = ['ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
            'ATGCGTNAGTTTAA', 'ATGCGTTAGAGTTAGAGCTAGTAA',
            'ATGTGTTAGATGCGTTAGAATAA', 'ATGNAATAA', 'AGATGCGTTAGTTTAA']
        for sequence in sequences:
            for cap in (1, 2, 4):
                for _ in range(3):
                    potentials = torch.randint(-3, 4, (len(sequence), 15), generator=generator).float()
                    result = self.compare(sequence, potentials, cap, oracle=True)
                    S.verify_selected(sequence, potentials, result, cap)
            self.assertEqual(self.compare(sequence, torch.zeros(len(sequence), 15), 4, oracle=True)['chains'], ())

    def test_known_short_actual_carry_GT_GC_N_negative_prefix_and_adjacent_genes(self):
        fixtures = [('ATGCGTTAGAGTTAGAGCTAGTAA', ((0, 4), (9, 10), (15, 16), (21, 24)), 4)]
        fixtures += [('ATGC' + intron + 'TTTAA', ((0, 4), (9, 14)), 2) for intron in ('GTTAG', 'GCTAG', 'GTNAG')]
        for sequence, chain, cap in fixtures:
            potentials = prefer_chain(sequence, chain).float()
            self.assertEqual(self.compare(sequence, potentials, cap, oracle=True)['chains'], (chain,))
        sequence, chain = 'ATGCGTTAGTTTAA', ((0, 4), (9, 14))
        potentials = prefer_chain(sequence, chain).float()
        potentials[:3, 13] = -10.
        self.assertEqual(self.compare(sequence, potentials, 2, oracle=True)['chains'], (chain,))
        sequence = 'ATGAAATAA' * 6
        potentials = torch.ones(len(sequence), 15); potentials[:, 0] = 0.
        self.assertEqual(len(self.compare(sequence, potentials, 1, oracle=True)['chains']), 6)

    def test_dense_small_random_FP32_integer_and_noninteger_identity(self):
        generator = torch.Generator().manual_seed(71)
        reduced = False
        for _ in range(12):
            sequence = 'ATG' + ''.join('ACGT'[i] for i in torch.randint(0, 4, (768,), generator=generator).tolist()) + 'TAA'
            for cap in (2, 4, 8):
                for potentials in (torch.randint(-3, 4, (len(sequence), 15), generator=generator).float(),
                                   torch.randn(len(sequence), 15, generator=generator)):
                    result = self.compare(sequence, potentials, cap)
                    reduced |= result['counts']['state_update_proposals'] < result['counts']['logical_exon_extensions']
        self.assertTrue(reduced)

    def test_128_exon_upper_bound_and_cap_are_not_phase_merged(self):
        sequence, potentials, chain, upper_bound = S.known_fixture()
        potentials = potentials.float()
        result = self.compare(sequence, potentials, 128)
        S.verify_selected(sequence, potentials, result, 128)
        self.assertEqual(result['chains'], (chain,))
        self.assertEqual(result['gain'], upper_bound)
        self.assertLess(self.compare(sequence, potentials, 127)['gain'], upper_bound)

    def test_30000bp_complete_exon_both_context_directions(self):
        chain, windows = ((3000, 33000),), [(0, 24576), (12288, 36864)]
        oriented = sequence_with_genes(36864, [chain])
        potentials = prefer_chain(oriented, chain).float()
        for strand in ('+', '-'):
            genomic = oriented if strand == '+' else reverse_complement(oriented)
            rows = [dict(window_index=i, sequence=genomic[a:b] if strand == '+' else reverse_complement(genomic[a:b]),
                local_potentials=potentials[a:b].clone() if strand == '+' else potentials[36864 - b:36864 - a].clone())
                for i, (a, b) in enumerate(windows)]
            joined = stitch_contig_emissions(genomic, windows, strand, rows)
            result = self.compare(joined['sequence'], joined['local_potentials'], 128)
            S.verify_selected(joined['sequence'], joined['local_potentials'], result, 128)
            self.assertEqual(result['chains'], (chain,))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
