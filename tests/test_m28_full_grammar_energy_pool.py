"""Independent synthetic complete-set tests for the fixed E/M energy mean."""
import unittest

import torch

from src.m28.full_grammar_energy_pool import equal_role_energy_pool as pool
from src.m28.full_grammar_skyline import decode_full_grammar_skyline as decode
from test_m28_emission_lattice import exhaustive_chains, direct_gain, prefer_chain


def all_sets(sequence, cap):
    chains = exhaustive_chains(sequence, cap)
    result = [()]

    def visit(prefix, end):
        for chain in chains:
            if chain[0][0] >= end:
                selected = prefix + (chain,)
                result.append(selected)
                visit(selected, chain[-1][1])
    visit((), 0)
    return result


class FullGrammarEnergyPoolTests(unittest.TestCase):
    def test_all_complete_set_energies_MAP_and_normalized_log_pool(self):
        generator = torch.Generator().manual_seed(251)
        for sequence in ('ATGAAATAAATGCCCTAA', 'ATGCGTNAGTTTAA',
                         'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA'):
            for cap in (1, 2, 4):
                left = torch.randn(len(sequence), 15, generator=generator)
                right = torch.randn(len(sequence), 15, generator=generator)
                values = pool(left, right)
                sets = all_sets(sequence, cap)
                energies = [torch.tensor([sum(direct_gain(x, c) for c in selected)
                    for selected in sets], dtype=torch.float64) for x in (left, right, values)]
                self.assertTrue(torch.allclose(energies[2], (energies[0]+energies[1])/2,
                                              rtol=0, atol=1e-12))
                logp = [x-torch.logsumexp(x, 0) for x in energies]
                product = (logp[0]+logp[1])/2
                product -= torch.logsumexp(product, 0)
                self.assertTrue(torch.allclose(logp[2], product, rtol=0, atol=1e-12))
                actual = decode(sequence, values, cap, 3)
                self.assertIn(actual['chains'], sets)
                self.assertAlmostEqual(actual['gain'], float(energies[2].max()), places=9)
                self.assertAlmostEqual(sum(actual['chain_gains']), actual['gain'], places=9)

    def test_source_immutability_FP64_promotion_symmetry_and_row_gauge(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        left = (torch.arange(270).reshape(18, 15)%11-5).float()
        right = (torch.arange(270).reshape(18, 15)%7-3).float()
        originals = left.clone(), right.clone()
        values = pool(left, right)
        self.assertEqual(values.dtype, torch.float64)
        self.assertTrue(torch.equal(values, pool(right, left)))
        self.assertTrue(torch.equal(left, originals[0]))
        self.assertTrue(torch.equal(right, originals[1]))
        shifts = torch.arange(18)[:, None]*8
        first, shifted = (decode(sequence, v, 4, 3) for v in
                          (values, pool(left+shifts, right-shifts*2)))
        self.assertEqual((first['chains'], first['gain']), (shifted['chains'], shifted['gain']))
        huge = torch.full((1, 15), 3e38)
        self.assertTrue(torch.isfinite(pool(huge, huge)).all())
        self.assertTrue(torch.equal(pool(huge, huge), huge.double()))

    def test_pool_can_select_a_legal_chain_neither_parent_emitted(self):
        sequence = 'ATGATGATGAAATAA'
        left = torch.full((15, 15), -100.)
        left[:, 0] = 0.; left[:, 1:4] = 0.; left[12:, 14] = 1.
        right = left.clone()
        left[:3, 13], right[:3, 13] = 2., -2.
        left[3:6, 13], right[3:6, 13] = -2., 2.
        left[6:9, 13] = right[6:9, 13] = 1.
        self.assertEqual(decode(sequence, left, 1, 3)['chains'], (((0, 15),),))
        self.assertEqual(decode(sequence, right, 1, 3)['chains'], (((3, 15),),))
        actual = decode(sequence, pool(left, right), 1, 3)
        self.assertEqual(actual['chains'], (((6, 15),),))
        self.assertEqual(actual['gain'], 6.)

    def test_actual_codon_carry_and_exon_cap_remain_binding(self):
        sequence = 'ATGTGTTAGATGCGTTAGAATAA'
        chain = ((9, 13), (18, 23))
        values = prefer_chain(sequence, chain).float()
        values[:3, 13] = 1000.; values[3, 7] = 30.
        values[6:9, 14] = -10000.; values[9, 11] = -10000.
        self.assertNotIn(((0, 4), (18, 23)), exhaustive_chains(sequence, 3))
        self.assertEqual(decode(sequence, pool(values, values), 3, 3)['chains'], (chain,))
        sequence = 'ATGCGTTAGTTTAA'; chain = ((0, 4), (9, 14))
        values = prefer_chain(sequence, chain).float()
        self.assertEqual(decode(sequence, pool(values, values), 1, 3)['chains'], ())
        self.assertEqual(decode(sequence, pool(values, values), 2, 3)['chains'], (chain,))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
