"""Independent complete-set enumeration of new read-only event posteriors."""
import itertools
import math
import unittest

import torch

from src.m28.full_grammar_gene_posterior import complete_gene_set_posteriors as posterior
from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition as partition
from test_m28_emission_lattice import exhaustive_chains, direct_gain, prefer_chain


def oracle(sequence, logits, cap):
    chains = exhaustive_chains(sequence, cap)
    sets, scores = [], []
    for keep in itertools.product((False, True), repeat=len(chains)):
        chosen = tuple(chain for chain, yes in zip(chains, keep) if yes)
        if any(a[-1][1] > b[0][0] for a, b in zip(chosen, chosen[1:])):
            continue
        sets.append(chosen); scores.append(sum(direct_gain(logits, g) for g in chosen))
    scores = torch.tensor(scores, dtype=torch.float64)
    log_z = float(torch.logsumexp(scores, 0))
    return sets, scores, log_z


class GenePosteriorTests(unittest.TestCase):
    def test_gene_event_and_exact_set_match_independent_enumeration(self):
        generator = torch.Generator().manual_seed(313)
        for sequence in ('ATGAAATAAATGCCCTAA', 'ATGATGTAA', 'ATGCGTNAGTTTAA',
                         'ATGGCTAGAAATAA', 'ATGCGTTAGAGTTAGAGCTAGTAA'):
            for cap in (1, 2, 4):
                for dtype in (torch.float32, torch.float64):
                    logits = torch.randn(len(sequence), 15, dtype=dtype, generator=generator)
                    sets, scores, log_z = oracle(sequence, logits, cap)
                    for chosen, score in zip(sets, scores.tolist()):
                        actual = posterior(sequence, logits, chosen, cap)
                        self.assertAlmostEqual(actual['log_partition'], log_z, places=9)
                        self.assertAlmostEqual(actual['queried_set_log_probability'], score-log_z, places=9)
                        for gene, mass in zip(chosen, actual['gene_log_probabilities']):
                            selected = torch.tensor([gene in s for s in sets], dtype=torch.bool)
                            expected = float(torch.logsumexp(scores[selected], 0))-log_z
                            self.assertAlmostEqual(mass, expected, places=9)
                        self.assertAlmostEqual(actual['empty_set_log_probability'], -log_z, places=9)

    def test_adjacent_and_overlapping_competition_not_sigmoid_or_set_product(self):
        sequence = 'ATGAAATAAATGCCCTAA'; values = torch.zeros(18, 15)
        actual = posterior(sequence, values, (((0, 9),), ((9, 18),)), 1)
        for probability in actual['gene_probabilities']:
            self.assertAlmostEqual(probability, .5, places=12)
        self.assertAlmostEqual(actual['queried_set_probability'], .25, places=12)
        sequence = 'ATGATGTAA'; values = torch.zeros(9, 15)
        for gene in (((0, 9),), ((3, 9),)):
            result = posterior(sequence, values, (gene,), 1)
            self.assertAlmostEqual(result['gene_probabilities'][0], 1/3, places=12)
            self.assertNotAlmostEqual(result['gene_probabilities'][0], .5)
        with self.assertRaises(ValueError):
            posterior(sequence, values, (((0, 9),), ((3, 9),)), 1)

    def test_true_carry_short_exons_GC_N_and_fixed_128_exon_cap(self):
        legal = (('ATGCGTNAGTTTAA', ((0, 4), (9, 14)), 2),
                 ('ATGGCTAGAAATAA', ((0, 3), (8, 14)), 2),
                 ('ATGCGTTAGAGTTAGAGCTAGTAA', ((0, 4), (9, 10), (15, 16), (21, 24)), 4))
        for sequence, chain, cap in legal:
            values = prefer_chain(sequence, chain)
            result = posterior(sequence, values, (chain,), cap)
            self.assertGreater(result['gene_probabilities'][0], .999999)
            self.assertAlmostEqual(result['log_partition'], float(partition(sequence, values, cap)['log_partition']), places=8)
        for sequence, chain, cap in (
                ('ATGTGTTAGATGCGTTAGAATAA', ((0, 4), (18, 23)), 3),
                ('ATGGTAGAAATAA', ((0, 3), (7, 13)), 2),
                ('ATGNAATAA', ((0, 9),), 1)):
            with self.assertRaises(ValueError):
                posterior(sequence, torch.zeros(len(sequence), 15), (chain,), cap)
        sequence = 'GTCAG'.join(['ATG']+['AAA']*126+['TAA'])
        chain = tuple((8*i, 8*i+3) for i in range(128))
        actual = posterior(sequence, prefer_chain(sequence, chain), (chain,), 128)
        self.assertGreater(actual['gene_probabilities'][0], .999999)
        with self.assertRaises(ValueError):
            posterior(sequence, prefer_chain(sequence, chain), (chain,), 127)

    def test_empty_set_row_shift_and_log_underflow_without_mutation(self):
        sequence = 'ATGAAATAA'
        values = torch.zeros(9, 15, dtype=torch.float64, requires_grad=True)
        first = posterior(sequence, values, (((0, 9),),), 128)
        shifted = posterior(sequence.lower(), values+torch.arange(9)[:, None]*8, (((0, 9),),), 128)
        for key in ('log_partition', 'gene_log_probabilities', 'queried_set_log_probability'):
            self.assertEqual(first[key], shifted[key])
        self.assertIsNone(values.grad)
        self.assertTrue(torch.equal(values, torch.zeros_like(values)))
        impossible = posterior('CCCCNNCC', torch.zeros(8, 15), (), 128)
        self.assertEqual((impossible['log_partition'], impossible['empty_set_probability']), (0., 1.))
        negative = torch.full((9, 15), -1000., dtype=torch.float64); negative[:, 0] = 0.
        rare = posterior(sequence, negative, (((0, 9),),), 1)
        self.assertEqual(rare['gene_probabilities'], (0.,))
        self.assertEqual(rare['gene_log_probabilities'], (-9000.,))
        self.assertFalse(rare['observation_model_used'])
        self.assertFalse(rare['posterior_threshold_applied'])
        self.assertEqual(rare['optimizer_steps'], 0)
