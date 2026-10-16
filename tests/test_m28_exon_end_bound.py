"""Artificial proof checks for skipping ONLY necessarily illegal end events."""
import itertools
import unittest

import torch

from src.m28.candidates import ORFIndex
from src.m28.exon_end_bound import exon_event_range, maximum_exon_end


STOP = {'TAA', 'TAG', 'TGA'}
CARRIES = [''] + list('ACGT') + [''.join(pair) for pair in itertools.product('ACGT', repeat=2)]


def direct_extension(sequence, a, b, carry, terminal):
    coding = carry + sequence[a:b]
    if set(coding) - set('ACGT'):
        return None
    codons = [coding[i:i + 3] for i in range(0, len(coding) - 2, 3)]
    if terminal:
        return '' if len(coding) % 3 == 0 and codons and codons[-1] in STOP and not any(c in STOP for c in codons[:-1]) else None
    if any(c in STOP for c in codons):
        return None
    leftover = len(coding) % 3
    return coding[-leftover:] if leftover else ''


def events_for(sequence):
    return sorted([(i, False) for i in range(1, len(sequence) - 1) if sequence[i:i + 2] in ('GT', 'GC')] +
                  [(i + 3, True) for i in range(len(sequence) - 2) if sequence[i:i + 3] in STOP])


class ExonEndBoundTests(unittest.TestCase):
    def compare(self, sequence):
        orf, events = ORFIndex(sequence), events_for(sequence)
        starts = {i for i in range(len(sequence) - 2) if sequence[i:i + 3] == 'ATG'}
        acceptors = {i + 2 for i in range(len(sequence) - 2) if sequence[i:i + 2] == 'AG'}
        for a in sorted(starts | acceptors):
            entries = [(True, '')] if a in starts else []
            if a in acceptors:
                entries += [(False, carry) for carry in CARRIES]
            for first, carry in entries:
                expected, original = [], []
                for b, terminal in events:
                    if b <= a or ((first or terminal) and b - a < 3):
                        continue
                    direct = direct_extension(sequence, a, b, carry, terminal)
                    actual = orf.extend(a, b, carry, terminal)
                    self.assertEqual(actual, direct)
                    if direct is not None:
                        expected.append((b, terminal, direct))
                    if actual is not None:
                        original.append((b, terminal, actual))
                lo, hi = exon_event_range(orf, events, a, carry)
                bounded = []
                for b, terminal in events[lo:hi]:
                    if (first or terminal) and b - a < 3:
                        continue
                    value = orf.extend(a, b, carry, terminal)
                    if value is not None:
                        bounded.append((b, terminal, value))
                self.assertEqual(bounded, original)
                self.assertEqual(bounded, expected)
                self.assertTrue(all(b <= maximum_exon_end(orf, a, carry) for b, _, _ in expected))

    def test_first_stop_retains_terminal_and_unfinished_prefix(self):
        sequence = 'ATGAAATAAGTCCCTAA'
        orf = ORFIndex(sequence)
        self.assertEqual(maximum_exon_end(orf, 0, ''), 9)
        self.assertEqual(orf.extend(0, 9, '', True), '')
        self.assertEqual(orf.extend(0, 8, '', False), 'TA')
        self.assertIsNone(orf.extend(0, 9, '', False))
        self.compare(sequence)

    def test_incoming_actual_carry_stop_keeps_short_middle_exon(self):
        sequence = 'ATGTGTTAGAACCCTAA'
        orf = ORFIndex(sequence)
        self.assertEqual(maximum_exon_end(orf, 9, 'T'), 10)
        self.assertEqual(orf.extend(9, 10, 'T', False), 'TA')
        self.assertIsNone(orf.extend(9, 11, 'T', False))
        self.assertNotEqual(maximum_exon_end(orf, 9, 'T'), maximum_exon_end(orf, 9, 'C'))
        self.compare(sequence)

    def test_ambiguity_is_outside_last_retained_endpoint(self):
        sequence = 'ATGCNCCCTAA'
        orf = ORFIndex(sequence)
        self.assertEqual(maximum_exon_end(orf, 0, ''), 4)
        self.assertEqual(maximum_exon_end(orf, 4, ''), 4)
        self.assertEqual(orf.extend(0, 4, '', False), 'C')
        self.assertIsNone(orf.extend(0, 5, '', False))
        self.compare('ATGCGTNAGTTTAAATGCNCCCTAA')

    def test_small_random_events_match_independent_DNA_and_original_order(self):
        generator = torch.Generator().manual_seed(29)
        for _ in range(25):
            middle = ''.join('ACGTN'[i] for i in torch.randint(0, 5, (72,), generator=generator).tolist())
            self.compare('ATG' + middle + 'AGCCCTAA')
        for sequence in ('ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA', 'ATGAAATAAATGCCCTAA'):
            self.compare(sequence)

    def test_long_complete_exon_has_no_window_or_length_cutoff(self):
        sequence = 'C' * 3000 + 'ATG' + 'C' * 29994 + 'TAA' + 'C' * 3864
        orf = ORFIndex(sequence)
        self.assertEqual(maximum_exon_end(orf, 3000, ''), 33000)
        lo, hi = exon_event_range(orf, events_for(sequence), 3000, '')
        self.assertIn((33000, True), events_for(sequence)[lo:hi])
        self.assertEqual(orf.extend(3000, 33000, '', True), '')


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
