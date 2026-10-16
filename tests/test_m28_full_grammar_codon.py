"""Independent spliced-CDS enumeration of codon energies, mass and gradients."""
import itertools
import math
import unittest

import torch

from src.m28.full_grammar_boundary_compiled import (
    _state_table, compiled_boundary_full_grammar_log_partition as parent_partition)
from src.m28.full_grammar_codon import (
    CODONS, _codon_edge_indices, centered_codon_scores,
    codon_full_grammar_log_partition as partition, codon_coarse_channel_nll as nll)
from src.m28.full_grammar_codon_skyline import decode_full_grammar_codon as decode
from src.m28.full_grammar_skyline import decode_full_grammar_skyline as parent_decode
from src.m28.full_grammar_observed_partition import InfeasibleObservation
from src.m28.phase_evidence_chain import candidate_roles
from test_m28_emission_lattice import exhaustive_chains, prefer_chain

STOPS = {'TAA', 'TAG', 'TGA'}
REGIONS = torch.tensor([0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1])


def oracle_coefficients(raw):
    # Independent ordering and normalization; no transition table or carry.
    words = [''.join(x) for x in itertools.product('ACGT', repeat=3)]
    sense = [i for i, c in enumerate(words) if c not in STOPS]
    stop = [i for i, c in enumerate(words) if c in STOPS]
    means = {False: raw.double()[sense].mean(), True: raw.double()[stop].mean()}
    return {c: raw.double()[i] - means[c in STOPS] for i, c in enumerate(words)}


def chain_energy(sequence, logits, raw, chain):
    a, b = chain[0][0], chain[-1][1]
    roles = candidate_roles(chain, len(sequence))
    pos = torch.arange(a, b)
    energy = (logits[pos, roles[pos]].double() - logits[pos, 0].double()).sum()
    cds = ''.join(sequence[x:y] for x, y in chain)
    coefficients = oracle_coefficients(raw)
    return energy + sum(coefficients[cds[i:i+3]] for i in range(0, len(cds), 3))


def oracle_sets(sequence, logits, raw, cap, allowed=None, logq=None):
    chains = exhaustive_chains(sequence, cap)
    scores, selected_sets = [], []
    for flags in itertools.product((False, True), repeat=len(chains)):
        chosen = tuple(chain for chain, flag in zip(chains, flags) if flag)
        if any(a[-1][1] > b[0][0] for a, b in zip(chosen, chosen[1:])):
            continue
        roles = torch.zeros(len(sequence), dtype=torch.long)
        energy = logits.double().sum() * 0 + raw.double().sum() * 0
        for chain in chosen:
            a, b = chain[0][0], chain[-1][1]
            roles[a:b] = candidate_roles(chain, len(sequence))[a:b]
            energy = energy + chain_energy(sequence, logits, raw, chain)
        if allowed is not None and not bool(allowed[torch.arange(len(sequence)), roles].all()):
            continue
        if logq is not None:
            energy = energy + logq[torch.arange(len(sequence)), roles].sum()
        scores.append(energy)
        selected_sets.append(chosen)
    return scores, selected_sets


class CodonFactorTests(unittest.TestCase):
    def gradients_close(self, actual, expected, dtype):
        tol = (1e-9, 1e-10) if dtype == torch.float64 else (2e-6, 2e-6)
        for a, b in zip(actual, expected):
            torch.testing.assert_close(a, b, rtol=tol[0], atol=tol[1])

    def test_mass_MAP_and_both_gradients_against_independent_complete_sets(self):
        rng = torch.Generator().manual_seed(301)
        sequences = ('ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
            'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA',
            'ATGCGTNAGTTTAA', 'ATGNAATAA', 'ATGAAAGTAGTCCCCAGAAATAA')
        for seq in sequences:
            for cap in (1, 2, 4):
                for dtype in (torch.float32, torch.float64):
                    with self.subTest(sequence=seq, cap=cap, dtype=dtype):
                        logits = torch.randn(len(seq), 15, generator=rng, dtype=dtype).requires_grad_(True)
                        raw = torch.randn(64, generator=rng, dtype=dtype).requires_grad_(True)
                        for partial in (False, True):
                            mask = torch.ones(len(seq), 15, dtype=torch.bool)
                            if partial:
                                mask = torch.rand(len(seq), 15, generator=rng) > .55
                                mask[:, 0] = True
                            scores, _ = oracle_sets(seq, logits, raw, cap, mask)
                            expected = torch.logsumexp(torch.stack(scores), 0)
                            result = partition(seq, logits, cap, raw, mask, cuts=(2, len(seq)-2), checkpoint_bp=3)
                            self.assertAlmostEqual(float(result['log_partition']), float(expected), places=9)
                            actual = torch.autograd.grad(-.7 * result['log_partition'], (logits, raw), retain_graph=True)
                            oracle = torch.autograd.grad(-.7 * expected, (logits, raw), retain_graph=True)
                            self.gradients_close(actual, oracle, dtype)
                        scores, sets = oracle_sets(seq, logits, raw, cap)
                        map_result = decode(seq, logits, cap, raw, checkpoint_bp=3)
                        self.assertAlmostEqual(map_result['gain'], max(float(s) for s in scores), places=9)
                        self.assertIn(map_result['chains'], sets)
                        for chain, gain in zip(map_result['chains'], map_result['chain_gains']):
                            self.assertAlmostEqual(gain, float(chain_energy(seq, logits, raw, chain)), places=9)

    def test_zero_reproduces_parent_scores_logits_gradients_chains_and_state_count(self):
        rng = torch.Generator().manual_seed(302)
        for seq in ('ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGCGTTAGAGTTAGAGCTAGTAA'):
            for cap in (1, 2, 4, 128):
                for dtype in (torch.float32, torch.float64):
                    values = torch.randn(len(seq), 15, generator=rng, dtype=dtype).requires_grad_(True)
                    raw = torch.zeros(64, dtype=dtype, requires_grad=True)
                    for mask in (None, torch.ones(len(seq), 15, dtype=torch.bool)):
                        a = partition(seq, values, cap, raw, mask, checkpoint_bp=2)
                        b = parent_partition(seq, values, cap, mask, checkpoint_bp=2)
                        self.assertEqual(float(a['log_partition']), float(b['log_partition']))
                        torch.testing.assert_close(torch.autograd.grad(a['log_partition'], values, retain_graph=True)[0],
                            torch.autograd.grad(b['log_partition'], values, retain_graph=True)[0], rtol=0, atol=0)
                    a, b = decode(seq, values, cap, raw, 3), parent_decode(seq, values, cap, 3)
                    for key in ('chains', 'gain', 'chain_gains'):
                        self.assertEqual(a[key], b[key])
                    self.assertEqual(a['counts']['state_count'], len(_state_table(cap)[0]))
                    self.assertEqual(_codon_edge_indices(cap).shape, _state_table(cap)[1].shape)
        self.assertEqual(decode('ATGAAATAA', torch.zeros(9, 15), 1, torch.zeros(64))['chains'], ())

    def test_split_codons_short_exons_and_each_stop_are_scored_once(self):
        fixtures = (
            ('ATGCGTTAGTTTAA', ((0, 4), (9, 14)), ('ATG', 'CTT', 'TAA')),
            ('ATGCGTTAGAGTTAGAGCTAGTAA', ((0, 4), (9, 10), (15, 16), (21, 24)), ('ATG', 'CAA', 'TAA')),
            ('ATGGCTAGAAATAA', ((0, 3), (8, 14)), ('ATG', 'AAA', 'TAA')),
            ('ATGAAATAG', ((0, 9),), ('ATG', 'AAA', 'TAG')),
            ('ATGAAATGA', ((0, 9),), ('ATG', 'AAA', 'TGA')))
        raw = (torch.arange(64, dtype=torch.float64) / 7).requires_grad_(True)
        for seq, chain, words in fixtures:
            self.assertEqual(''.join(seq[a:b] for a, b in chain), ''.join(words))
            logits = prefer_chain(seq, chain).requires_grad_(True)
            mask = torch.nn.functional.one_hot(candidate_roles(chain, len(seq)), 15).bool()
            result = partition(seq, logits, len(chain), raw, mask, tuple(range(1, len(seq))), 2)
            expected = chain_energy(seq, logits, raw, chain)
            self.assertAlmostEqual(float(result['log_partition']), float(expected), places=10)
            self.gradients_close(torch.autograd.grad(result['log_partition'], (logits, raw), retain_graph=True),
                torch.autograd.grad(expected, (logits, raw), retain_graph=True), torch.float64)
            mapped = decode(seq, logits, len(chain), raw, 2)
            self.assertEqual(mapped['chains'], (chain,))
            self.assertEqual(mapped['counts']['completed_codon_count'], len(words))

    def test_all_chunk_cuts_and_checkpoint_spans_preserve_codon_gradient_history(self):
        seq = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        logits = (torch.arange(len(seq)*15).reshape(-1, 15).double()/43).requires_grad_(True)
        raw = (torch.arange(64).double().sin()/3).requires_grad_(True)
        whole = partition(seq, logits, 4, raw, checkpoint_bp=3)['log_partition']
        expected = torch.autograd.grad(whole, (logits, raw), retain_graph=True)
        reference_map = decode(seq, logits, 4, raw, 3)
        for cuts in [(i,) for i in range(1, len(seq))] + [tuple(range(1, len(seq)))]:
            result = partition(seq, logits, 4, raw, cuts=cuts, checkpoint_bp=2)
            self.assertEqual(float(result['log_partition']), float(whole))
            self.gradients_close(torch.autograd.grad(result['log_partition'], (logits, raw), retain_graph=True),
                                 expected, torch.float64)
        for span in (1, 2, 5, 256):
            result = partition(seq, logits, 4, raw, checkpoint_bp=span)
            self.assertEqual(float(result['log_partition']), float(whole))
            self.gradients_close(torch.autograd.grad(result['log_partition'], (logits, raw), retain_graph=True),
                                 expected, torch.float64)
            mapped = decode(seq, logits, 4, raw, span)
            for key in ('chains', 'gain', 'chain_gains'):
                self.assertEqual(mapped[key], reference_map[key])

    def test_original_fine_coarse_unknown_channel_absolute_mass_and_gradients(self):
        seq, a, b = 'ATGCGTTAGTTTAA', 3, 12
        logits = (torch.arange(14*15).reshape(14, 15).double()/47).requires_grad_(True)
        raw = torch.arange(64).double().cos().requires_grad_(True)
        targets = dict(seg=torch.tensor([7, -100, -100, -100, -100, -100, -100, -100, 14]),
                       region=torch.tensor([1, 2, -100, 2, 2, -100, 1, 1, 1]))
        for epsilon_value in (.1, 2/3):
            epsilon = torch.tensor(epsilon_value, dtype=torch.float64, requires_grad=True)
            fine, region = targets['seg'], targets['region']
            allowed = torch.ones(14, 15, dtype=torch.bool)
            logq = torch.zeros(14, 15, dtype=torch.float64)
            for j in range(b-a):
                if fine[j] >= 0:
                    allowed[a+j] = torch.arange(15) == fine[j]
                elif region[j] >= 0:
                    logq[a+j] = torch.where(REGIONS == region[j], torch.log1p(-epsilon), torch.log(epsilon/2))
            all_scores, _ = oracle_sets(seq, logits, raw, 2)
            observed_scores, _ = oracle_sets(seq, logits, raw, 2, allowed, logq)
            expected = torch.logsumexp(torch.stack(all_scores), 0) - torch.logsumexp(torch.stack(observed_scores), 0)
            result = nll(seq, logits, 2, raw, a, b, targets, epsilon, cuts=(a, b), checkpoint_bp=2)
            self.assertAlmostEqual(float(result['loss']), float(expected), places=10)
            self.gradients_close(torch.autograd.grad(result['loss'], (logits, raw, epsilon), retain_graph=True),
                torch.autograd.grad(expected, (logits, raw, epsilon), retain_graph=True), torch.float64)
            if epsilon_value == 2/3:
                fine_mass = partition(seq, logits, 2, raw, allowed)['log_partition']
                expected_uniform = torch.logsumexp(torch.stack(all_scores), 0) - fine_mass + 5*math.log(3)
                self.assertAlmostEqual(float(result['loss']), float(expected_uniform), places=10)
        unknown = dict(seg=torch.full((b-a,), -100), region=torch.full((b-a,), -100))
        loss = nll(seq, logits, 2, raw, a, b, unknown, .1)['loss']
        self.assertEqual(float(loss), 0.)
        for grad in torch.autograd.grad(loss, (logits, raw)):
            # Centering/chunk VJPs may accumulate machine-scale cancellation
            # residues even when the two log masses are exactly identical.
            torch.testing.assert_close(grad, torch.zeros_like(grad), rtol=0, atol=1e-12)

    def test_centering_row_shifts_and_nonzero_split_codon_can_change_MAP(self):
        seq, chain = 'ATGCGTTAGTTTAA', ((0, 4), (9, 14))
        logits, raw = torch.zeros(14, 15, dtype=torch.float64), torch.zeros(64, dtype=torch.float64)
        self.assertEqual(decode(seq, logits, 2, raw)['chains'], ())
        raw[CODONS.index('CTT')] = 12.
        self.assertEqual(decode(seq, logits, 2, raw)['chains'], (chain,))
        scores = centered_codon_scores(raw)
        self.assertAlmostEqual(float(scores[[c not in STOPS for c in CODONS]].sum()), 0., places=12)
        shifted = raw + torch.tensor([3. if c in STOPS else 7. for c in CODONS])
        torch.testing.assert_close(centered_codon_scores(shifted), scores, rtol=0, atol=2e-15)
        a = partition(seq, logits, 2, raw)['log_partition']
        b = partition(seq, logits + torch.arange(14)[:, None]*8, 2, shifted)['log_partition']
        self.assertAlmostEqual(float(a), float(b), places=12)

    def test_infeasible_partial_genes_carry_stops_and_invalid_input_remain_excluded(self):
        seq = 'ATGAAATA'
        mask = torch.ones(8, 15, dtype=torch.bool)
        mask[:3] = False; mask[:3, 13] = True
        with self.assertRaises(InfeasibleObservation):
            partition(seq, torch.zeros(8, 15), 2, torch.zeros(64), mask, cuts=(1, 2, 7))
        seq, chain = 'ATGTGTTAGATGCGTTAGAATAA', ((0, 4), (18, 23))
        mask = torch.nn.functional.one_hot(candidate_roles(chain, len(seq)), 15).bool()
        with self.assertRaises(InfeasibleObservation):
            partition(seq, torch.zeros(len(seq), 15), 3, torch.ones(64), mask)
        for raw in (torch.zeros(63), torch.full((64,), float('nan'))):
            with self.assertRaises(ValueError):
                partition('ATGAAATAA', torch.zeros(9, 15), 1, raw)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
