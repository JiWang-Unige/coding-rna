"""Artificial equivalence and concrete exclusions before merging histories."""
import itertools
import math
import unittest

import torch

from src.m28.candidates import ORFIndex
from src.m28.emission_lattice import _Margins
from src.m28.long_exon_phase import long_exon_phase


STOP = {'TAA', 'TAG', 'TGA'}
CARRIES = [''] + list('ACGT') + [''.join(p) for p in itertools.product('ACGT', repeat=2)]


def direct_extension(sequence, a, b, carry, terminal):
    coding = carry + sequence[a:b]
    if set(coding) - set('ACGT'):
        return None
    codons = [coding[i:i + 3] for i in range(0, len(coding) - 2, 3)]
    if terminal:
        return '' if len(coding) % 3 == 0 and codons and codons[-1] in STOP and not any(c in STOP for c in codons[:-1]) else None
    if any(c in STOP for c in codons):
        return None
    n = len(coding) % 3
    return coding[-n:] if n else ''


class LongExonPhaseTests(unittest.TestCase):
    def compare(self, sequence, a):
        orf = ORFIndex(sequence)
        generator = torch.Generator().manual_seed(53)
        margins = _Margins(torch.randint(-3, 4, (len(sequence), 15), generator=generator).float())
        groups = {}
        for carry in CARRIES:
            phase = long_exon_phase(sequence, a, carry)
            actual, independent, role_deltas = [], [], []
            for b in range(a + 3, len(sequence) + 1):
                for terminal in (False, True):
                    expected = direct_extension(sequence, a, b, carry, terminal)
                    observed = orf.extend(a, b, carry, terminal)
                    self.assertEqual(observed, expected)
                    actual.append(observed)
                    independent.append(expected)
                    if observed is not None:
                        role_deltas.append((b, terminal, margins.exon_gain(a, b, len(carry), False, terminal)))
            if phase is None:
                self.assertTrue(all(value is None for value in independent))
            else:
                signature = (actual, role_deltas)
                if phase in groups:
                    self.assertEqual(signature, groups[phase])
                else:
                    groups[phase] = signature

    def test_all_first_codon_DNA_and_actual_carries(self):
        for bases in itertools.product('ACGT', repeat=3):
            self.compare('CCCCAG' + ''.join(bases) + 'CCGTGCAGCCCTAA', 6)

    def test_random_ambiguity_and_all_endpoints_against_independent_DNA(self):
        generator = torch.Generator().manual_seed(59)
        for _ in range(25):
            middle = ''.join('ACGTN'[i] for i in torch.randint(0, 5, (40,), generator=generator).tolist())
            self.compare('CCCCAG' + middle + 'CCCTAA', 6)

    def test_short_middle_exon_cannot_discard_ineligible_long_carry(self):
        sequence, a = 'CCCCAGAAGTCCCTAA', 6
        orf = ORFIndex(sequence)
        self.assertIsNone(long_exon_phase(sequence, a, 'T'))
        self.assertEqual(long_exon_phase(sequence, a, 'C'), 1)
        self.assertEqual(orf.extend(a, a + 1, 'T', False), 'TA')
        self.assertEqual(orf.extend(a, a + 1, 'C', False), 'CA')
        self.assertIsNone(orf.extend(a, a + 3, 'T', False))
        self.assertEqual(orf.extend(a, a + 3, 'C', False), 'G')
        # Both actual carries survive the short exon and differ at the next exon.
        self.assertIsNone(direct_extension('ACCC', 0, 4, 'TA', False))
        self.assertEqual(direct_extension('ACCC', 0, 4, 'CA', False), '')

    def test_FP64_entry_winner_is_not_first_arrival_after_common_delta(self):
        # A supported finite FP64 tie, not overflow or an adversarial input.
        earlier, later, delta = 1.0, math.nextafter(1.0, math.inf), 2.0
        self.assertEqual(later, earlier + float(torch.tensor(2.0 ** -52, dtype=torch.float32)))
        self.assertGreater(later, earlier)
        self.assertEqual(earlier + delta, later + delta)
        original_winner, best = None, -math.inf
        for name, score in (('earlier_actual_carry', earlier), ('later_actual_carry', later)):
            if score + delta > best:
                original_winner, best = name, score + delta
        preselected = max((('earlier_actual_carry', earlier), ('later_actual_carry', later)), key=lambda row: row[1])[0]
        self.assertEqual(original_winner, 'earlier_actual_carry')
        self.assertNotEqual(original_winner, preselected)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
