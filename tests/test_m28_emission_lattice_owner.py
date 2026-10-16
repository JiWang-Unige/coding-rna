"""Synthetic original-width owner/assembly counterexamples, not a raw adapter."""
import unittest

import torch

from src.m28.assembly import assemble_predictions, nearest_owner
from src.m28.core import select_nonoverlapping
from src.m28.emission_lattice import decode_emission_lattice
from src.m28.phase_evidence_chain import candidate_roles


WINDOWS = [(0, 24576), (12288, 36864)]
STOP = {'TAA', 'TAG', 'TGA'}


def sequence_with_genes(length, chains):
    bases = list('C' * length)
    for chain in chains:
        assert len(chain) == 1
        a, b = chain[0]
        bases[a:a + 3] = 'ATG'
        bases[b - 3:b] = 'TAA'
    return ''.join(bases)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def direct_gain(sequence, potentials, chain):
    """Check actual concatenated CDS independently of decoder ORFIndex."""
    coding = ''.join(sequence[a:b] for a, b in chain)
    assert coding[:3] == 'ATG' and len(coding) % 3 == 0
    codons = [coding[i:i + 3] for i in range(0, len(coding), 3)]
    assert codons[-1] in STOP and not any(c in STOP for c in codons[:-1])
    roles = candidate_roles(chain, len(sequence))
    positions = torch.arange(chain[0][0], chain[-1][1])
    return float((potentials[positions, roles[positions]].double() -
                  potentials[positions, 0].double()).sum())


def predicted(result):
    return {row['chain']: row['score'] for row in result['predictions']}


class EmissionLatticeOwnerTests(unittest.TestCase):
    def check_selected_before_owner(self, strand, owned, nonowned, markers):
        sequence = sequence_with_genes(24576, [owned, nonowned])
        potentials = torch.zeros(24576, 15)
        for position, role, value in markers:
            potentials[position, role] = value
        self.assertEqual(direct_gain(sequence, potentials, owned), 5.)
        self.assertEqual(direct_gain(sequence, potentials, nonowned), 10.)
        natural = decode_emission_lattice(sequence, potentials, 128)
        self.assertEqual(natural['chains'], (nonowned,))
        self.assertEqual(natural['chain_gains'], (10.,))
        self.assertEqual(natural['gain'], 10.)
        self.assertFalse(natural['reference_used'])

        if strand == '+':
            global_owned, global_nonowned = owned, nonowned
            genomic = sequence + 'C' * 12288
            second_sequence = genomic[12288:36864]
            second_chain = tuple((a - 12288, b - 12288) for a, b in global_nonowned)
        else:
            global_owned = tuple((24576 - b, 24576 - a) for a, b in reversed(owned))
            global_nonowned = tuple((24576 - b, 24576 - a) for a, b in reversed(nonowned))
            genomic = reverse_complement(sequence) + 'G' * 12288
            second_sequence = reverse_complement(genomic[12288:36864])
            second_chain = tuple((36864 - b, 36864 - a)
                                 for a, b in reversed(global_nonowned))
        self.assertEqual(nearest_owner(global_owned, WINDOWS), 0)
        self.assertEqual(nearest_owner(global_nonowned, WINDOWS), 1)
        neutral = torch.zeros(24576, 15)
        self.assertEqual(direct_gain(second_sequence, neutral, second_chain), 0.)
        second = decode_emission_lattice(second_sequence, neutral, 128)
        self.assertEqual(second['chains'], ())
        full = assemble_predictions([
            dict(window_index=0, chains=[owned, nonowned], score=[5., 10.]),
            dict(window_index=1, chains=[second_chain], score=[0.])], WINDOWS, strand)
        selected_only = assemble_predictions([
            dict(window_index=0, chains=natural['chains'], score=natural['chain_gains']),
            dict(window_index=1, chains=second['chains'], score=second['chain_gains'])], WINDOWS, strand)
        self.assertEqual(predicted(full), {global_owned: 5.})
        self.assertEqual(predicted(selected_only), {})

    def test_positive_strand_selection_then_owner_loses_positive_owned_chain(self):
        self.check_selected_before_owner('+', ((15000, 20001),), ((15010, 23002),),
                                         [(15005, 3, 5.), (21000, 3, 10.)])

    def test_negative_strand_selection_then_owner_loses_positive_owned_chain(self):
        self.check_selected_before_owner('-', ((10000, 16000),), ((0, 11001),),
                                         [(15000, 3, 5.), (100, 2, 10.)])

    def test_owner_filtered_local_optimum_still_loses_cross_window_alternative(self):
        a, b, c = ((8000, 14000),), ((12000, 21000),), ((18001, 24001),)
        genomic = sequence_with_genes(36864, [a, b, c])
        potentials = torch.zeros(24576, 15)
        for position, value, role in [(8100, 5., 2), (16000, 8., 2), (22000, 10., 1)]:
            potentials[position, role] = value
        sequence = genomic[:24576]
        for chain, score, owner in [(a, 5., 0), (b, 8., 0), (c, 10., 1)]:
            self.assertEqual(direct_gain(sequence, potentials, chain), score)
            self.assertEqual(nearest_owner(chain, WINDOWS), owner)
        natural = decode_emission_lattice(sequence, potentials, 128)
        self.assertEqual(natural['chains'], (a, c))
        self.assertEqual(natural['gain'], 15.)

        second_chain = tuple((start - 12288, end - 12288) for start, end in c)
        second_sequence = genomic[12288:36864]
        second_potentials = torch.zeros(24576, 15)
        second_potentials[22000 - 12288, 1] = 10.
        self.assertEqual(direct_gain(second_sequence, second_potentials, second_chain), 10.)
        second = decode_emission_lattice(second_sequence, second_potentials, 128)
        self.assertEqual(second['chains'], (second_chain,))
        full = assemble_predictions([
            dict(window_index=0, chains=[a, b, c], score=[5., 8., 10.]),
            dict(window_index=1, chains=[second_chain], score=[10.])], WINDOWS, '+')

        # This is the proposed owner-before-local-DP shortcut, not an existing
        # A42 adapter: an owned local optimum need not be globally sufficient.
        local = select_nonoverlapping([a, b], [5., 8.])
        self.assertEqual([a, b][local[0]], b)
        self.assertEqual(len(local), 1)
        shortcut = assemble_predictions([
            dict(window_index=0, chains=[[a, b][i] for i in local],
                 score=[[5., 8.][i] for i in local]),
            dict(window_index=1, chains=second['chains'], score=second['chain_gains'])], WINDOWS, '+')
        self.assertEqual(predicted(full), {a: 5., c: 10.})
        self.assertEqual(predicted(shortcut), {c: 10.})


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
