"""Artificial partial-observation mass and gradients against complete-set enumeration."""
import itertools
import math
import unittest

import torch

from src.m28.full_grammar_observed_partition import (
    InfeasibleObservation, _Observation, allowed_roles_from_local_targets,
    observed_full_grammar_log_partition, partial_observation_nll)
from src.m28.full_grammar_partition_compiled_reverse import (
    _accumulate_reverse_nodes, compiled_reverse_full_grammar_log_partition)
from src.m28.labels import REGION_COLUMNS
from src.m28.phase_evidence_chain import candidate_roles
from test_m28_emission_lattice import exhaustive_chains


def exhaustive_observed(sequence, logits, cap, allowed):
    # Independent concatenated-DNA chain oracle, then all compatible subsets.
    chains = exhaustive_chains(sequence, cap)
    values, positions, scores = logits.double(), torch.arange(len(sequence)), []
    for selected in itertools.product((False, True), repeat=len(chains)):
        chosen = [chain for chain, yes in zip(chains, selected) if yes]
        if any(a[-1][1] > b[0][0] for a, b in zip(chosen, chosen[1:])):
            continue
        roles = torch.zeros(len(sequence), dtype=torch.long)
        for chain in chosen:
            a, b = chain[0][0], chain[-1][1]
            roles[a:b] = candidate_roles(chain, len(sequence))[a:b]
        if bool(allowed[positions, roles].all()):
            scores.append((values[positions, roles] - values[:, 0]).sum())
    return (torch.logsumexp(torch.stack(scores), 0) if scores else None), len(scores)


def exact_roles(roles):
    return torch.nn.functional.one_hot(roles, num_classes=15).bool()


class ObservedPartitionTests(unittest.TestCase):
    def compare(self, sequence, logits, cap, allowed, scale=1., check_nll=False):
        logits = logits.detach().requires_grad_(True)
        expected, sets = exhaustive_observed(sequence, logits, cap, allowed)
        if expected is None:
            with self.assertRaises(InfeasibleObservation):
                observed_full_grammar_log_partition(sequence, logits, cap, allowed)
            return None, sets, None
        result = observed_full_grammar_log_partition(sequence, logits, cap, allowed)
        value = result['log_partition']
        self.assertEqual(value.dtype, torch.float64)
        self.assertAlmostEqual(float(value.detach()), float(expected.detach()), places=10)
        actual_grad, = torch.autograd.grad(value * scale, logits, retain_graph=True)
        expected_grad, = torch.autograd.grad(expected * scale, logits, retain_graph=True)
        tolerance = (1e-9, 1e-10) if logits.dtype == torch.float64 else (1e-6, 1e-6)
        torch.testing.assert_close(actual_grad, expected_grad, rtol=tolerance[0], atol=tolerance[1])
        total = compiled_reverse_full_grammar_log_partition(sequence, logits, cap)
        self.assertLessEqual(float(value.detach()), float(total['log_partition'].detach()) + 1e-10)
        self.assertFalse(result['enumerated_candidate_pool_used'])
        self.assertFalse(result['Viterbi_winner_pruning'])
        if check_nll:
            nll = partial_observation_nll(sequence, logits, cap, allowed)
            all_expected, _ = exhaustive_observed(sequence, logits, cap, torch.ones_like(allowed))
            expected_loss = all_expected - expected
            self.assertAlmostEqual(float(nll['loss'].detach()), float(expected_loss.detach()), places=10)
            grad, = torch.autograd.grad(nll['loss'] * scale, logits)
            oracle_grad, = torch.autograd.grad(expected_loss * scale, logits)
            torch.testing.assert_close(grad, oracle_grad, rtol=tolerance[0], atol=tolerance[1])
        return result, sets, actual_grad

    def test_all_unknown_matches_original_partition_counts_and_zero_loss(self):
        generator = torch.Generator().manual_seed(101)
        for sequence in ('ATGCGTTAGTTTAA', 'ATGCGTTAGAGTTAGAGCTAGTAA', 'CCCCCCCC'):
            for dtype in (torch.float32, torch.float64):
                logits = torch.randn(len(sequence), 15, dtype=dtype, generator=generator).requires_grad_(True)
                allowed = torch.ones_like(logits, dtype=torch.bool)
                observed = observed_full_grammar_log_partition(sequence, logits, 128, allowed)
                original = compiled_reverse_full_grammar_log_partition(sequence, logits, 128)
                self.assertEqual(float(observed['log_partition']), float(original['log_partition']))
                for key, count in original['counts'].items():
                    self.assertEqual(observed['counts'][key], count)
                nll = partial_observation_nll(sequence, logits, 128, allowed)
                self.assertEqual(float(nll['loss']), 0.)
                self.assertTrue(torch.equal(torch.autograd.grad(nll['loss'], logits)[0], torch.zeros_like(logits)))
        self.assertGreater(len(_accumulate_reverse_nodes.nopython_signatures), 0)

    def test_unannotated_genes_remain_marginalized_not_negative(self):
        sequence = 'ATGAAATAA' * 6
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        allowed = torch.ones_like(logits, dtype=torch.bool)
        allowed[:3] = False; allowed[:3, 13] = True
        allowed[18:27] = False; allowed[18:27, 0] = True
        result, sets, gradient = self.compare(sequence, logits, 1, allowed, check_nll=True)
        self.assertEqual(sets, 16)
        self.assertAlmostEqual(float(result['log_partition']), math.log(16), places=12)
        self.assertGreater(float(gradient[9, 13]), 0.)
        self.assertEqual(float(gradient[18:27].abs().sum()), 0.)

    def test_fully_observed_spliced_roles_leave_one_set(self):
        sequence = 'ATGCGTTAGTTTAA'
        roles = torch.tensor([13, 13, 13, 7, 4, 4, 4, 4, 4, 11, 3, 14, 14, 14])
        logits = torch.arange(14 * 15, dtype=torch.float64).reshape(14, 15) / 100
        result, sets, _ = self.compare(sequence, logits, 2, exact_roles(roles), check_nll=True)
        self.assertEqual(sets, 1)
        expected = (logits[torch.arange(14), roles] - logits[:, 0]).sum()
        self.assertAlmostEqual(float(result['log_partition']), float(expected), places=12)

    def test_fine_coarse_unknown_and_padding_match_original_policy(self):
        targets = dict(seg=torch.tensor([13, -100, -100, -100, 4, -100, -100, -100, -100]),
                       region=torch.tensor([0, 1, 2, -100, 1, 0, -100, -100, -100]))
        allowed = allowed_roles_from_local_targets(targets, 7)
        self.assertEqual(tuple(allowed.shape), (7, 15))
        self.assertEqual(torch.where(allowed[0])[0].tolist(), [13])
        self.assertEqual(torch.where(allowed[4])[0].tolist(), [4])
        for row, label in ((1, 1), (2, 2), (5, 0)):
            self.assertEqual(torch.where(allowed[row])[0].tolist(), list(REGION_COLUMNS[label]))
        self.assertTrue(bool(allowed[[3, 6]].all()))
        targets['region'][-1] = 0
        with self.assertRaises(ValueError): allowed_roles_from_local_targets(targets, 7)

    def test_background_before_between_and_after_genes_is_constrained(self):
        sequence = 'CCCATGAAATAACCCATGCCCTAACCC'
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        for position in (0, 13, len(sequence) - 1):
            allowed = torch.ones_like(logits, dtype=torch.bool)
            allowed[position] = False; allowed[position, 1] = True
            _, sets, _ = self.compare(sequence, logits, 1, allowed)
            self.assertEqual(sets, 0)
        allowed = torch.ones_like(logits, dtype=torch.bool)
        allowed[3:6] = False; allowed[3:6, 13] = True
        allowed[15:18] = False; allowed[15:18, 13] = True
        _, sets, _ = self.compare(sequence, logits, 1, allowed)
        self.assertEqual(sets, 1)

    def test_intron_histories_are_filtered_at_maturity_and_after_activation(self):
        sequence = 'ATGAAAGTAGTCCCCAGAAATAA'
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        allowed = torch.ones_like(logits, dtype=torch.bool)
        role = int(candidate_roles(((0, 9), (17, 23)), len(sequence))[7])
        allowed[7] = False; allowed[7, role] = True
        _, sets, _ = self.compare(sequence, logits, 2, allowed)
        self.assertEqual(sets, 1)
        allowed = torch.ones_like(logits, dtype=torch.bool)
        allowed[11] = False; allowed[11, 0] = True
        result, sets, _ = self.compare(sequence, logits, 2, allowed)
        self.assertEqual(sets, 1)
        self.assertGreater(result['counts']['invalidated_observed_intron_states'], 0)
        self.assertGreater(result['counts']['rejected_observed_mature_states'], 0)

    def test_short_exon_actual_carry_and_cap_keep_observation_support(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        chain = ((0, 4), (9, 10), (15, 16), (21, 24))
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        allowed = exact_roles(candidate_roles(chain, len(sequence)))
        _, sets, _ = self.compare(sequence, logits, 4, allowed, check_nll=True)
        self.assertEqual(sets, 1)
        _, sets, _ = self.compare(sequence, logits, 3, allowed)
        self.assertEqual(sets, 0)

    def test_random_partial_masks_and_NLL_full_gradients(self):
        generator = torch.Generator().manual_seed(102)
        sequences = ['ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
            'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA', 'ATGCGTNAGTTTAA']
        for sequence in sequences:
            for cap in (1, 2, 4):
                chains = exhaustive_chains(sequence, cap)
                witness = candidate_roles(chains[-1], len(sequence)) if chains else torch.zeros(len(sequence), dtype=torch.long)
                allowed = torch.rand(len(sequence), 15, generator=generator) > .55
                allowed[torch.arange(len(sequence)), witness] = True
                for dtype, scale in ((torch.float32, -.7), (torch.float64, 2.3)):
                    logits = torch.randn(len(sequence), 15, dtype=dtype, generator=generator)
                    with self.subTest(sequence=sequence, cap=cap, dtype=dtype):
                        self.compare(sequence, logits, cap, allowed, scale, check_nll=True)

    def test_exon_role_precedence_matches_literal_chain_roles(self):
        chains = [((2, 14),), ((2, 5), (10, 13)), ((2, 6), (11, 22)),
            ((2, 7), (12, 13), (18, 19), (24, 28), (33, 36)),
            ((2, 6), (11, 13), (18, 23)), ((2, 5),)]
        for chain in chains:
            roles = candidate_roles(chain, 80)
            observation = _Observation(exact_roles(roles))
            coding = 0
            for ordinal, (a, b) in enumerate(chain):
                first, terminal = ordinal == 0, ordinal == len(chain) - 1
                self.assertTrue(observation.exon_allowed(a, b, coding % 3, first, terminal))
                for position in (a, b - 1):
                    allowed = exact_roles(roles); allowed[position] = False
                    self.assertFalse(_Observation(allowed).exon_allowed(a, b, coding % 3, first, terminal))
                coding += b - a

    def test_boundary_infeasibility_and_invalid_rows_are_explicit(self):
        sequence = 'ATGAAATA'
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        allowed = torch.ones_like(logits, dtype=torch.bool)
        allowed[0] = False; allowed[0, 13] = True
        self.compare(sequence, logits, 2, allowed)
        sequence = 'CCCCCCCC'
        logits = torch.zeros(8, 15, dtype=torch.float64)
        allowed = torch.ones(8, 15, dtype=torch.bool); allowed[3] = False
        self.compare(sequence, logits, 2, allowed)
        for mask in (torch.ones(7, 15, dtype=torch.bool), torch.ones(8, 15)):
            with self.assertRaises(ValueError): observed_full_grammar_log_partition(sequence, logits, 2, mask)
        with self.assertRaises(ValueError):
            observed_full_grammar_log_partition(sequence, torch.full((8, 15), float('nan')), 2, torch.ones(8, 15, dtype=torch.bool))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
