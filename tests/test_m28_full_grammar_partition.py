"""Artificial exhaustive complete-chain sets, logZ and emission gradients."""
import itertools
import math
import unittest

import torch

from src.m28.emission_lattice import decode_emission_lattice
from src.m28.full_grammar_partition import full_grammar_log_partition
from src.m28.phase_evidence_chain import candidate_roles
from test_m28_emission_lattice import exhaustive_chains


def exhaustive_partition(sequence, logits, cap):
    # Naive concatenated-DNA oracle does not call ORFIndex or this new DP.
    chains = exhaustive_chains(sequence, cap)
    scores = []
    for chain in chains:
        roles = candidate_roles(chain, len(sequence))
        positions = torch.arange(chain[0][0], chain[-1][1])
        values = logits.double()
        scores.append((values[positions, roles[positions]] - values[positions, 0]).sum())
    zero = logits[0, 0].double()*0
    sets = []
    for selected in itertools.product((False, True), repeat=len(chains)):
        chosen = [chain for chain, yes in zip(chains, selected) if yes]
        if any(a[-1][1] > b[0][0] for a, b in zip(chosen, chosen[1:])):
            continue
        sets.append(sum((score for score, yes in zip(scores, selected) if yes), zero))
    return torch.logsumexp(torch.stack(sets), 0), chains, len(sets)


class FullGrammarPartitionTests(unittest.TestCase):
    def compare(self, sequence, logits, cap):
        logits = logits.detach().requires_grad_(True)
        result = full_grammar_log_partition(sequence, logits, cap)
        actual = result['log_partition']
        expected, chains, sets = exhaustive_partition(sequence, logits, cap)
        self.assertEqual(actual.dtype, torch.float64)
        self.assertAlmostEqual(float(actual.detach()), float(expected.detach()), places=10)
        grad, = torch.autograd.grad(actual, logits, retain_graph=True)
        oracle_grad, = torch.autograd.grad(expected, logits)
        tolerance = (1e-9, 1e-10) if logits.dtype == torch.float64 else (1e-6, 1e-6)
        torch.testing.assert_close(grad, oracle_grad, rtol=tolerance[0], atol=tolerance[1])
        self.assertFalse(result['reference_used'])
        self.assertFalse(result['enumerated_candidate_pool_used'])
        self.assertFalse(result['Viterbi_winner_pruning'])
        viterbi = decode_emission_lattice(sequence, logits.detach(), cap)
        self.assertGreaterEqual(float(actual.detach()) + 1e-10, viterbi['gain'])
        return result, chains, sets, grad

    def test_exhaustive_gene_sets_and_all_emission_gradients(self):
        sequences = ['ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
            'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA',
            'ATGCGTNAGTTTAA', 'ATGNAATAA', 'AGATGCGTTAGTTTAA']
        generator = torch.Generator().manual_seed(83)
        for sequence in sequences:
            for cap in (1, 2, 4):
                for logits in (torch.randint(-3, 4, (len(sequence), 15), generator=generator).double(),
                               torch.randn(len(sequence), 15, generator=generator, dtype=torch.float64),
                               torch.randn(len(sequence), 15, generator=generator)):
                    with self.subTest(sequence=sequence, cap=cap):
                        self.compare(sequence, logits, cap)

    def test_null_enters_once_and_six_adjacent_genes_have_64_sets(self):
        sequence = 'ATGAAATAA'*6
        result, chains, sets, _ = self.compare(sequence, torch.zeros(len(sequence), 15, dtype=torch.float64), 1)
        self.assertEqual(len(chains), 6)
        self.assertEqual(sets, 64)
        self.assertAlmostEqual(float(result['log_partition'].detach()), math.log(64), places=12)
        for sequence in ('CCCCCCCC', 'ATGNAATAA', 'ATGAAATA'):
            result, chains, sets, grad = self.compare(sequence, torch.zeros(len(sequence), 15, dtype=torch.float64), 4)
            self.assertEqual(chains, [])
            self.assertEqual(sets, 1)
            self.assertEqual(float(result['log_partition'].detach()), 0.)
            self.assertTrue(torch.equal(grad, torch.zeros_like(grad)))

    def test_two_donor_histories_same_state_are_summed_not_winner_pruned(self):
        sequence = 'ATGAAAGTAGTCCCCAGAAATAA'
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        result, chains, sets, _ = self.compare(sequence, logits, 2)
        self.assertEqual(chains, [((0, 6), (17, 23)), ((0, 9), (17, 23))])
        self.assertEqual(sets, 3)
        self.assertAlmostEqual(float(result['log_partition'].detach()), math.log(3), places=12)
        self.assertEqual(result['counts']['maximum_active_states'], 1)
        self.assertEqual(decode_emission_lattice(sequence, logits, 2)['gain'], 0.)

    def test_literal_spliced_roles_have_posterior_weighted_gradient(self):
        sequence = 'ATGCGTTAGTTTAA'
        roles = [13, 13, 13, 7, 4, 4, 4, 4, 4, 11, 3, 14, 14, 14]
        logits = torch.arange(14*15, dtype=torch.float64).reshape(14, 15)/100
        result, chains, sets, grad = self.compare(sequence, logits, 2)
        # The same DNA also has the unspliced ATG--TAG chain (0,9).
        # Keep that legal competitor; do not alter DNA or prune the oracle.
        self.assertEqual(chains, [((0, 4), (9, 14)), ((0, 9),)])
        self.assertEqual(sets, 3)
        score = sum(float(logits[i, role]-logits[i, 0]) for i, role in enumerate(roles))
        competitor = [13, 13, 13, 1, 2, 3, 14, 14, 14]
        competitor_score = sum(float(logits[i, role]-logits[i, 0]) for i, role in enumerate(competitor))
        logz = math.log(1+math.exp(score)+math.exp(competitor_score))
        self.assertAlmostEqual(float(result['log_partition'].detach()), logz, places=12)
        expected = torch.zeros_like(logits)
        probability = math.exp(score-logz)
        for i, role in enumerate(roles):
            expected[i, role] += probability
            expected[i, 0] -= probability
        competitor_probability = math.exp(competitor_score-logz)
        for i, role in enumerate(competitor):
            expected[i, role] += competitor_probability
            expected[i, 0] -= competitor_probability
        torch.testing.assert_close(grad, expected, rtol=1e-12, atol=1e-12)

    def test_short_actual_carry_internal_STOP_and_cap_remain_distinct(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        full = self.compare(sequence, torch.zeros(len(sequence), 15, dtype=torch.float64), 4)
        capped = self.compare(sequence, torch.zeros(len(sequence), 15, dtype=torch.float64), 3)
        self.assertIn(((0, 4), (9, 10), (15, 16), (21, 24)), full[1])
        self.assertNotIn(((0, 4), (9, 10), (15, 16), (21, 24)), capped[1])
        self.assertGreater(full[2], capped[2])
        sequence = 'ATGTGTTAGATGCGTTAGAATAA'
        _, chains, _, _ = self.compare(sequence, torch.zeros(len(sequence), 15, dtype=torch.float64), 3)
        self.assertIn(((9, 13), (18, 23)), chains)
        self.assertNotIn(((0, 4), (18, 23)), chains)

    def test_FP32_common_logit_shift_and_large_finite_gain(self):
        sequence = 'ATGCGTTAGTTTAA'
        logits = torch.arange(len(sequence)*15).reshape(len(sequence), 15).float()%7-3
        first = full_grammar_log_partition(sequence, logits.requires_grad_(True), 2)['log_partition']
        shifted = (logits.detach()+torch.arange(len(sequence))[:, None]*8).requires_grad_(True)
        second = full_grammar_log_partition(sequence.lower(), shifted, 2)['log_partition']
        self.assertEqual(float(first.detach()), float(second.detach()))
        torch.testing.assert_close(torch.autograd.grad(first, logits)[0],
                                   torch.autograd.grad(second, shifted)[0], rtol=0, atol=0)
        for scale in (-1000., 1000.):
            values = torch.full((len(sequence), 15), scale, dtype=torch.float64)
            values[:, 0] = 0.
            result, _, _, grad = self.compare(sequence, values, 2)
            self.assertTrue(bool(torch.isfinite(result['log_partition'])))
            self.assertTrue(bool(torch.isfinite(grad).all()))

    def test_shape_nonfinite_and_explicit_cap(self):
        sequence = 'ATGAAATAA'
        for values in (torch.zeros(8, 15), torch.zeros(9, 14),
                       torch.full((9, 15), float('nan')), torch.zeros(9, 15, dtype=torch.int64)):
            with self.assertRaises(ValueError): full_grammar_log_partition(sequence, values, 2)
        with self.assertRaises(ValueError): full_grammar_log_partition(sequence, torch.zeros(9, 15), 0)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
