"""Artificial joint-context integration, never actual raw inference or labels."""
import unittest

import torch

from src.m28.assembly import assemble_predictions
from src.m28.contig_emission_lattice import (context_segments, decode_contig_emission_lattice,
                                           stitch_contig_emissions)
from src.m28.emission_lattice import decode_emission_lattice
from test_m28_emission_lattice import direct_gain, oracle_gain, prefer_chain
from test_m28_emission_lattice_owner import WINDOWS, reverse_complement, sequence_with_genes


def records(sequence, windows, strand, matrices):
    return [dict(window_index=i, sequence=sequence[a:b] if strand == '+' else
                 reverse_complement(sequence[a:b]), local_potentials=matrix)
            for i, ((a, b), matrix) in enumerate(zip(windows, matrices))]


class ContigEmissionLatticeTests(unittest.TestCase):
    def check_mapping(self, length, windows):
        sequence = ('ACGTN' * length)[:length]
        for strand in ('+', '-'):
            matrices = [torch.arange((b - a) * 15).reshape(b - a, 15).float() + i * 1000000
                        for i, (a, b) in enumerate(windows)]
            saved = [matrix.clone() for matrix in matrices]
            joined = stitch_contig_emissions(sequence.lower(), windows, strand,
                                             iter(records(sequence, windows, strand, matrices)))
            expected = torch.empty(length, 15)
            for oriented_p in range(length):
                p = oriented_p if strand == '+' else length - 1 - oriented_p
                owners = [i for i, (a, b) in enumerate(windows) if a <= p < b]
                owner = min(owners, key=lambda i: (abs(2 * p + 1 - sum(windows[i])),
                                                  windows[i][0], i))
                a, b = windows[owner]
                local = p - a if strand == '+' else b - 1 - p
                expected[oriented_p] = matrices[owner][local]
            self.assertTrue(torch.equal(joined['local_potentials'], expected))
            self.assertEqual(joined['sequence'], sequence if strand == '+' else
                             sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1])
            for actual, old in zip(matrices, saved):
                self.assertTrue(torch.equal(actual, old))

    def test_base_owner_ties_negative_mapping_and_role_columns(self):
        windows = [(0, 9), (4, 13), (6, 15)]
        self.assertEqual(context_segments(15, windows), ((0, 0, 7), (1, 7, 10), (2, 10, 15)))
        self.check_mapping(15, windows)
        self.check_mapping(12, [(0, 8), (6, 12)])

    def test_original_width_shifted_last_window_and_single_short_contig(self):
        self.check_mapping(40000, [(0, 24576), (12288, 36864), (15424, 40000)])
        self.check_mapping(7, [(0, 7)])

    def test_incomplete_misordered_or_wrong_context_is_rejected(self):
        sequence, windows = 'ATGAAATAA', [(0, 6), (3, 9)]
        rows = records(sequence, windows, '+', [torch.zeros(6, 15), torch.zeros(6, 15)])
        malformed = [rows[:1], rows[::-1], rows + rows[:1],
                     [rows[0], dict(rows[1], sequence='CCCCCC')],
                     [rows[0], dict(rows[1], local_potentials=torch.zeros(7, 15))],
                     [rows[0], dict(rows[1], local_potentials=torch.zeros(6, 15).double())],
                     [rows[0], dict(rows[1], local_potentials=torch.full((6, 15), float('nan')))]]
        for bad in malformed:
            with self.assertRaises(ValueError):
                stitch_contig_emissions(sequence, windows, '+', bad)
        for bad in ([(0, 3), (4, 9)], [(1, 9)], [(0, 8)], [(0, 9), (3, 8)]):
            with self.assertRaises(ValueError):
                stitch_contig_emissions(sequence, bad, '+', rows)

    def test_prior_selected_only_owner_losses_under_new_context_objective(self):
        for strand, owned, nonowned, markers in [
                ('+', ((15000, 20001),), ((15010, 23002),), [(15005, 3, 5.), (21000, 3, 10.)]),
                ('-', ((10000, 16000),), ((0, 11001),), [(15000, 3, 5.), (100, 2, 10.)])]:
            local = sequence_with_genes(24576, [owned, nonowned])
            sequence = local + 'C' * 12288 if strand == '+' else reverse_complement(local) + 'G' * 12288
            matrices = [torch.zeros(24576, 15), torch.zeros(24576, 15)]
            for position, role, value in markers:
                matrices[0][position, role] = value
            rows = records(sequence, WINDOWS, strand, matrices)
            local_results = [decode_emission_lattice(row['sequence'], row['local_potentials'], 128)
                             for row in rows]
            old = assemble_predictions([dict(window_index=i, chains=r['chains'], score=r['chain_gains'])
                                        for i, r in enumerate(local_results)], WINDOWS, strand)
            self.assertEqual(old['predictions'], [])
            joint = decode_contig_emission_lattice(sequence, WINDOWS, strand, iter(rows), 128)
            expected = owned if strand == '+' else tuple((24576 - b, 24576 - a)
                                                          for a, b in reversed(owned))
            self.assertEqual(joint['genomic_chains'], (expected,))
            self.assertEqual(joint['chain_gains'], (5.,))
            self.assertEqual(joint['gain'], 5.)
            self.assertFalse(joint['reference_used'])
            # A different score definition, not recovery under old chain-owner scores:
            joined = stitch_contig_emissions(sequence, WINDOWS, strand, rows)
            shifted_nonowned = nonowned if strand == '+' else tuple((a + 12288, b + 12288)
                                                                    for a, b in nonowned)
            self.assertEqual(direct_gain(joined['local_potentials'], shifted_nonowned), 0.)

    def test_joint_selection_retains_cross_window_alternative(self):
        a, b, c = ((8000, 14000),), ((12000, 21000),), ((18001, 24001),)
        sequence = sequence_with_genes(36864, [a, b, c])
        matrices = [torch.zeros(24576, 15), torch.zeros(24576, 15)]
        for p, value, role in [(8100, 5., 2), (16000, 8., 2), (22000, 10., 1)]:
            matrices[0][p, role] = value
        matrices[1][22000 - 12288, 1] = 10.
        joint = decode_contig_emission_lattice(sequence, WINDOWS, '+',
                                               records(sequence, WINDOWS, '+', matrices), 128)
        self.assertEqual(joint['genomic_chains'], (a, c))
        self.assertEqual(joint['chain_gains'], (5., 10.))
        self.assertEqual(joint['gain'], 15.)

    def test_complete_gene_longer_than_window_across_context_seam_both_strands(self):
        chain = ((3000, 33000),)
        oriented = sequence_with_genes(36864, [chain])
        logits = prefer_chain(oriented, chain).float()
        for strand in ('+', '-'):
            sequence = oriented if strand == '+' else reverse_complement(oriented)
            matrices = [logits[a:b].clone() if strand == '+' else
                        logits[len(sequence) - b:len(sequence) - a].clone() for a, b in WINDOWS]
            rows = records(sequence, WINDOWS, strand, matrices)
            for row in rows:
                self.assertEqual(decode_emission_lattice(row['sequence'], row['local_potentials'], 128)['chains'], ())
            joint = decode_contig_emission_lattice(sequence, WINDOWS, strand, rows, 128)
            self.assertEqual(joint['chains'], (chain,))
            self.assertEqual(joint['genomic_chains'], (chain if strand == '+' else ((3864, 33864),),))
            self.assertEqual(joint['gain'], 300000.)
            self.assertEqual(joint['chain_gains'], (300000.,))

    def test_splice_codon_carry_and_cap_across_multiple_seams(self):
        sequence, chain = 'ATGCGTNAGTTTAA', ((0, 4), (9, 14))
        windows = [(0, 5), (3, 8), (6, 11), (9, 14)]
        logits = prefer_chain(sequence, chain).float()
        rows = records(sequence, windows, '+', [logits[a:b].clone() for a, b in windows])
        joined = stitch_contig_emissions(sequence, windows, '+', rows)
        self.assertTrue(torch.equal(joined['local_potentials'], logits))
        self.assertEqual(oracle_gain(sequence, logits, 2)[0], 140.)
        self.assertEqual(decode_contig_emission_lattice(sequence, windows, '+', rows, 2)['chains'], (chain,))
        self.assertEqual(decode_contig_emission_lattice(sequence, windows, '+', rows, 1)['chains'], ())

    def test_new_context_random_small_oracle_and_zero_null(self):
        generator = torch.Generator().manual_seed(23)
        for oriented in ('ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGCGTTAGAGTTAGAGCTAGTAA'):
            length = len(oriented)
            windows = [(0, length - 4), (4, length)]
            for strand in ('+', '-'):
                sequence = oriented if strand == '+' else reverse_complement(oriented)
                for _ in range(4):
                    matrices = [torch.randint(-3, 4, (b - a, 15), generator=generator).float()
                                for a, b in windows]
                    rows = records(sequence, windows, strand, matrices)
                    joined = stitch_contig_emissions(sequence, windows, strand, rows)
                    expected, legal = oracle_gain(joined['sequence'], joined['local_potentials'], 4)
                    joint = decode_contig_emission_lattice(sequence, windows, strand, rows, 4)
                    self.assertEqual(joint['gain'], expected)
                    last = 0
                    for chain, gain in zip(joint['chains'], joint['chain_gains']):
                        self.assertIn(chain, legal)
                        self.assertGreaterEqual(chain[0][0], last)
                        self.assertEqual(gain, direct_gain(joined['local_potentials'], chain))
                        last = chain[-1][1]
                    self.assertEqual(sum(joint['chain_gains']), joint['gain'])
                null_rows = records(sequence, windows, strand, [torch.zeros(b - a, 15) for a, b in windows])
                null = decode_contig_emission_lattice(sequence, windows, strand, null_rows, 4)
                self.assertEqual((null['chains'], null['gain']), ((), 0.))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
