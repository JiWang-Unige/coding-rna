"""Synthetic exhaustive-path oracle; no research checkpoint, DNA or reference."""
import unittest

import torch

from src.m28.emission_lattice import _Margins, decode_emission_lattice
from src.m28.phase_evidence_chain import candidate_roles


STOP = {'TAA', 'TAG', 'TGA'}


def exhaustive_chains(sequence, cap):
    """Naive DNA/CDS concatenation, deliberately independent of ORFIndex."""
    sequence = sequence.upper()
    starts = [i for i in range(len(sequence) - 2) if sequence[i:i + 3] == 'ATG']
    acceptors = [i + 2 for i in range(len(sequence) - 2) if sequence[i:i + 2] == 'AG']
    ends = sorted([(i, False) for i in range(1, len(sequence) - 1)
                   if sequence[i:i + 2] in ('GT', 'GC')] +
                  [(i + 3, True) for i in range(len(sequence) - 2)
                   if sequence[i:i + 3] in STOP])
    result = set()

    def visit(a, prefix, chain):
        for b, terminal in ends:
            if b <= a or ((not chain or terminal) and b - a < 3):
                continue
            coding = prefix + sequence[a:b]
            if any(base not in 'ACGT' for base in coding):
                continue
            codons = [coding[i:i + 3] for i in range(0, len(coding) - 2, 3)]
            new = chain + ((a, b),)
            if terminal:
                if len(coding) % 3 == 0 and codons[-1] in STOP and not any(c in STOP for c in codons[:-1]):
                    result.add(new)
            elif not any(c in STOP for c in codons) and len(new) < cap:
                for entry in acceptors:
                    if entry >= b + 5:
                        visit(entry, coding, new)
    for a in starts:
        visit(a, '', ())
    return sorted(result)


def direct_gain(logits, chain):
    roles = candidate_roles(chain, len(logits))
    positions = torch.arange(chain[0][0], chain[-1][1])
    return float((logits[positions, roles[positions]].double() - logits[positions, 0].double()).sum())


def oracle_gain(sequence, logits, cap):
    """Independent interval recurrence over ALL exhaustively enumerated chains."""
    chains = exhaustive_chains(sequence, cap)
    by_end = {}
    for chain in chains:
        by_end.setdefault(chain[-1][1], []).append((chain[0][0], direct_gain(logits, chain)))
    best = [0.] * (len(sequence) + 1)
    for end in range(1, len(sequence) + 1):
        best[end] = max([best[end - 1]] + [best[start] + score for start, score in by_end.get(end, [])])
    return best[-1], set(chains)


def prefer_chain(sequence, chain):
    logits = torch.full((len(sequence), 15), -100., dtype=torch.float64)
    logits[:, 0] = 0.
    roles = candidate_roles(chain, len(sequence))
    positions = torch.arange(chain[0][0], chain[-1][1])
    logits[positions, roles[positions]] = 10.
    return logits


class EmissionLatticeTests(unittest.TestCase):
    def check_oracle(self, sequence, logits, cap):
        expected, legal = oracle_gain(sequence, logits, cap)
        actual = decode_emission_lattice(sequence, logits, cap)
        self.assertAlmostEqual(actual['gain'], expected, places=9)
        self.assertEqual(len(actual['chains']), len(actual['chain_gains']))
        last = 0
        for chain, score in zip(actual['chains'], actual['chain_gains']):
            self.assertIn(chain, legal)
            self.assertGreaterEqual(chain[0][0], last)
            self.assertGreater(score, 0.)
            self.assertLessEqual(len(chain), cap)
            self.assertAlmostEqual(score, direct_gain(logits, chain), places=9)
            last = chain[-1][1]
        self.assertAlmostEqual(sum(actual['chain_gains']), actual['gain'], places=9)
        self.assertFalse(actual['reference_used'])
        return actual

    def test_role_precedence_matches_direct_sum_all_phases_and_short_exons(self):
        generator = torch.Generator().manual_seed(0)
        logits = torch.randn(80, 15, dtype=torch.float64, generator=generator)
        margins = _Margins(logits)
        chains = [((2, 14),), ((2, 5), (10, 13)), ((2, 6), (11, 22)),
                  ((2, 7), (12, 13), (18, 19), (24, 28), (33, 36)),
                  ((2, 6), (11, 13), (18, 23)), ((2, 5),)]
        for chain in chains:
            value, coding, previous = 0., 0, None
            for index, (a, b) in enumerate(chain):
                if previous is not None:
                    value += margins.interval(previous, a, 4 + (coding - 1) % 3)
                value += margins.exon_gain(a, b, coding % 3, index == 0, index == len(chain) - 1)
                coding += b - a
                previous = b
            self.assertAlmostEqual(value, direct_gain(logits, chain), places=11)

    def test_exhaustive_oracle_across_single_spliced_and_competing_genes(self):
        sequences = ['ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
                     'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA']
        generator = torch.Generator().manual_seed(12)
        for sequence in sequences:
            for cap in (1, 2, 4):
                logits = torch.randn(len(sequence), 15, dtype=torch.float64, generator=generator)
                self.check_oracle(sequence, logits, cap)

    def test_adjacent_genes_and_no_four_completion_limit(self):
        sequence = 'ATGAAATAA' * 6
        logits = torch.ones(len(sequence), 15, dtype=torch.float64)
        logits[:, 0] = 0.
        actual = self.check_oracle(sequence, logits, 1)
        self.assertEqual(actual['chains'], tuple(((i, i + 9),) for i in range(0, 54, 9)))
        self.assertEqual(actual['gain'], 54.)

    def test_zero_and_negative_objectives_return_null(self):
        sequence = 'ATGAAATAAATGCCCTAA'
        logits = torch.zeros(len(sequence), 15)
        self.assertEqual(self.check_oracle(sequence, logits, 4)['chains'], ())
        logits[:, 0] = 2.
        self.assertEqual(self.check_oracle(sequence, logits, 4)['chains'], ())

    def test_negative_prefix_with_positive_completion_is_retained(self):
        sequence = 'ATGCGTTAGTTTAA'
        chain = ((0, 4), (9, 14))
        logits = prefer_chain(sequence, chain)
        logits[:3, 13] = -10.
        self.assertLess(_Margins(logits).exon_gain(0, 4, 0, True, False), 0.)
        self.assertEqual(self.check_oracle(sequence, logits, 2)['chains'], (chain,))

    def test_one_base_middle_and_three_base_terminal_precedence(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        chain = ((0, 4), (9, 10), (15, 16), (21, 24))
        self.assertEqual(self.check_oracle(sequence, prefer_chain(sequence, chain), 4)['chains'], (chain,))
        sequence = 'ATGGTTAGAAATAA'
        chain = ((0, 3), (8, 14))
        self.assertEqual(self.check_oracle(sequence, prefer_chain(sequence, chain), 2)['chains'], (chain,))

    def test_actual_carry_not_only_phase_protects_legal_history(self):
        sequence = 'ATGTGTTAGATGCGTTAGAATAA'
        chain = ((9, 13), (18, 23))
        logits = prefer_chain(sequence, chain)
        logits[:3, 13] = 1000.
        logits[3, 7] = 30.
        logits[6:9, 14] = -10000.  # Prevent the earlier start ending at the first raw TAG.
        logits[9, 11] = -10000.  # Prevent a different legal earlier chain through this acceptor.
        legal = exhaustive_chains(sequence, 3)
        self.assertIn(chain, legal)
        self.assertNotIn(((0, 4), (18, 23)), legal)  # T + AA is internal STOP; C + AA is not.
        self.assertEqual(self.check_oracle(sequence, logits, 3)['chains'], (chain,))

    def test_gt_gc_and_intronic_ambiguity_but_not_ambiguous_coding(self):
        chain = ((0, 4), (9, 14))
        for intron in ('GTTAG', 'GCTAG', 'GTNAG'):
            sequence = 'ATGC' + intron + 'TTTAA'
            self.assertEqual(self.check_oracle(sequence, prefer_chain(sequence, chain), 2)['chains'], (chain,))
        sequence = 'ATGNAATAA'
        self.assertEqual(self.check_oracle(sequence, torch.ones(9, 15), 2)['chains'], ())

    def test_cap_does_not_repair_a_donor_into_stop(self):
        sequence = 'ATGCGTTAGTTTAA'
        chain = ((0, 4), (9, 14))
        logits = prefer_chain(sequence, chain)
        self.assertEqual(self.check_oracle(sequence, logits, 1)['chains'], ())
        self.assertEqual(self.check_oracle(sequence, logits, 2)['chains'], (chain,))

    def test_softmask_offset_and_deterministic_ties(self):
        sequence = 'ATGCGTTAGTTTAA'
        generator = torch.Generator().manual_seed(2)
        logits = torch.randint(-3, 4, (len(sequence), 15), generator=generator).double()
        expected = self.check_oracle(sequence, logits, 3)
        offset = torch.arange(len(sequence), dtype=torch.float64)[:, None] * 8
        actual = decode_emission_lattice(sequence.lower(), logits + offset, 3)
        self.assertEqual(actual, expected)
        self.assertEqual(decode_emission_lattice(sequence, logits, 3), expected)

    def test_integer_random_small_sequences_match_exhaustive_oracle(self):
        generator = torch.Generator().manual_seed(17)
        alphabet = 'ACGT'
        for _ in range(25):
            middle = ''.join(alphabet[i] for i in torch.randint(0, 4, (24,), generator=generator).tolist())
            sequence = 'ATG' + middle + 'TAA'
            logits = torch.randint(-3, 4, (len(sequence), 15), generator=generator).double()
            self.check_oracle(sequence, logits, 3)

    def test_misaligned_nonfinite_and_missing_cap_are_rejected(self):
        sequence = 'ATGAAATAA'
        for logits in (torch.zeros(8, 15), torch.zeros(9, 14), torch.full((9, 15), float('nan'))):
            with self.assertRaises(ValueError):
                decode_emission_lattice(sequence, logits, 2)
        with self.assertRaises(ValueError):
            decode_emission_lattice(sequence, torch.zeros(9, 15), 0)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
