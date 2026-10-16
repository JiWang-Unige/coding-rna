"""Same artificial fixtures, independent DNA oracle and previous batch counts."""
import math
import unittest

import torch

import test_m28_emission_lattice_batched as B
from src.m28.emission_lattice import _Live, _Solution
from src.m28.emission_lattice_phase_block import _phase_block_winners, decode_phase_block_emission_lattice


class PhaseBlockEmissionLatticeTests(B.BatchedEmissionLatticeTests):
    def compare(self, sequence, potentials, cap, oracle=False):
        previous = super().compare(sequence, potentials, cap, oracle=oracle)
        actual = decode_phase_block_emission_lattice(sequence, potentials, cap)
        for key in ('chains', 'gain', 'chain_gains', 'reference_used'):
            self.assertEqual(actual[key], previous[key])
        for key, value in previous['counts'].items():
            self.assertEqual(actual['counts'][key], value)
        self.assertLessEqual(actual['counts']['long_array_blocks'], actual['counts']['long_history_groups'])
        self.assertLessEqual(actual['counts']['long_array_blocks'], actual['counts']['long_nonempty_history_groups'])
        self.assertLessEqual(actual['counts']['long_nonempty_history_groups'], actual['counts']['long_history_groups'])
        return actual

    def test_every_exit_adds_FP64_before_first_arrival_and_matches_scalar(self):
        generator = torch.Generator().manual_seed(61)
        scores = [1., math.nextafter(1., math.inf), -1., math.nextafter(-1., -math.inf), 0., -0.]
        scores += torch.randn(15, generator=generator).double().tolist()
        deltas = [2., -2., 0., -0., 2. ** 53, -(2. ** 53), 2. ** -52]
        deltas += torch.randn(15, generator=generator).double().tolist()
        carries = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        block = []
        # Missing columns differ across count rows; even the first column is
        # absent in some rows. There is no synthetic zero-score fallback.
        for count in range(1, 6):
            rows = [_Live(scores[j], _Solution(), ((0, 3),) * count, carry)
                    for j, carry in enumerate(carries) if (j + count) % 3 != 0]
            block.append(rows)
        lookup, winners, values, width = _phase_block_winners(block, deltas)
        self.assertEqual(width, 16)
        for i, rows in enumerate(block):
            for k, delta in enumerate(deltas):
                expected = max(rows, key=lambda row: row.score + delta)
                self.assertIs(lookup[i][winners[i][k]], expected)
                self.assertEqual(values[i][k], expected.score + delta)
        pair = [_Live(score, _Solution(), ((0, 3),), carry) for score, carry in zip(scores[:2], ('A', 'C'))]
        lookup, winners, values, _ = _phase_block_winners([pair], [2.])
        self.assertIs(lookup[0][winners[0][0]], pair[0])
        self.assertEqual(values, [[3.]])
        sequence = 'ATGAGTTAGATGCGTTAGCCTAA'
        potentials = torch.zeros(len(sequence), 15)
        potentials[3, 7] = 1.; potentials[9, 13] = 1.
        potentials[12, 7] = 2. ** -52; potentials[9, 11] = -4.; potentials[20, 14] = 2.
        result = self.compare(sequence, potentials, 2)
        self.assertEqual(result['gain'], 3.)
        self.assertEqual(result['chains'], (((0, 4), (18, 23)),))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
