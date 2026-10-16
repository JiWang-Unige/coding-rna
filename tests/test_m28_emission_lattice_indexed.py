"""Artificial SAME-target indexed DP equivalence; no scientific inputs."""
import unittest

import torch

import profile_emission_lattice_synthetic as S
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.emission_lattice import decode_emission_lattice
from src.m28.emission_lattice_indexed import decode_indexed_emission_lattice
from test_m28_emission_lattice import oracle_gain, prefer_chain
from test_m28_emission_lattice_owner import reverse_complement, sequence_with_genes


class IndexedEmissionLatticeTests(unittest.TestCase):
    def compare(self, sequence, potentials, cap, oracle=False):
        old = decode_emission_lattice(sequence, potentials, cap)
        indexed = decode_indexed_emission_lattice(sequence, potentials, cap)
        self.assertEqual(indexed, old)  # Includes every traversal count, chain and tie identity.
        S.verify_selected(sequence, potentials, indexed, cap)
        if oracle:
            expected, legal = oracle_gain(sequence, potentials, cap)
            self.assertEqual(indexed['gain'], expected)
            self.assertTrue(all(chain in legal for chain in indexed['chains']))
        return indexed

    def test_small_oracle_and_original_all_counts_ties_and_carries(self):
        generator = torch.Generator().manual_seed(31)
        sequences = ['ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
                     'ATGCGTNAGTTTAA', 'ATGCGTTAGAGTTAGAGCTAGTAA',
                     'ATGTGTTAGATGCGTTAGAATAA', 'ATGNAATAA']
        for sequence in sequences:
            for cap in (1, 2, 4):
                for _ in range(3):
                    potentials = torch.randint(-3, 4, (len(sequence), 15), generator=generator).float()
                    self.compare(sequence, potentials, cap, oracle=True)
            self.assertEqual(self.compare(sequence, torch.zeros(len(sequence), 15), 4, oracle=True)['chains'], ())

    def test_original_128_exon_known_upper_bound_and_cap(self):
        sequence, potentials, chain, upper_bound = S.known_fixture()
        result = self.compare(sequence, potentials, 128)
        self.assertEqual(result['chains'], (chain,))
        self.assertEqual(result['gain'], upper_bound)
        self.assertLess(self.compare(sequence, potentials, 127)['gain'], upper_bound)

    def test_joint_context_long_gene_no_width_truncation_both_strands(self):
        chain, windows = ((3000, 33000),), [(0, 24576), (12288, 36864)]
        oriented = sequence_with_genes(36864, [chain])
        potentials = prefer_chain(oriented, chain).float()
        for strand in ('+', '-'):
            genomic = oriented if strand == '+' else reverse_complement(oriented)
            rows = [dict(window_index=i,sequence=genomic[a:b] if strand == '+' else
                         reverse_complement(genomic[a:b]),local_potentials=potentials[a:b].clone()
                         if strand == '+' else potentials[36864 - b:36864 - a].clone())
                    for i, (a, b) in enumerate(windows)]
            joined = stitch_contig_emissions(genomic,windows,strand,rows)
            self.assertEqual(self.compare(joined['sequence'],joined['local_potentials'],128)['chains'], (chain,))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
