"""Metadata reuse preserves all prior phase-block answers and counters."""
import unittest

import torch

import test_m28_emission_lattice_phase_block as P
from src.m28.emission_lattice_cached_metadata import decode_cached_metadata_emission_lattice
from src.m28.long_exon_phase import long_exon_phase
from test_m28_emission_lattice import prefer_chain


class CachedMetadataEmissionLatticeTests(P.PhaseBlockEmissionLatticeTests):
    def compare(self, sequence, potentials, cap, oracle=False):
        previous = super().compare(sequence, potentials, cap, oracle=oracle)
        actual = decode_cached_metadata_emission_lattice(sequence, potentials, cap)
        for key in ('chains', 'gain', 'chain_gains', 'reference_used'):
            self.assertEqual(actual[key], previous[key])
        for key, value in previous['counts'].items():
            self.assertEqual(actual['counts'][key], value)
        acceptors = sum(sequence[i:i + 2] == 'AG' for i in range(len(sequence) - 2))
        self.assertLessEqual(actual['counts']['long_phase_metadata_evaluations'], 21 * acceptors)
        return actual

    def test_cached_None_long_phase_still_keeps_short_actual_carry_and_cap(self):
        sequence = 'ATGT' + 'GTTAG' + 'A' + 'GTTAG' + 'CAAATAA'
        chain = ((0, 4), (9, 10), (15, 22))
        self.assertIsNone(long_exon_phase(sequence, 9, 'T'))
        potentials = prefer_chain(sequence, chain).float()
        self.assertEqual(self.compare(sequence, potentials, 3, oracle=True)['chains'], (chain,))
        self.assertNotIn(chain, self.compare(sequence, potentials, 2, oracle=True)['chains'])


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
