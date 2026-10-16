"""Independent combined FP64-pool and scaled-message event accounting."""
import math
import unittest
from unittest.mock import patch

import torch
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_energy_pool import equal_role_energy_pool
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries
from test_m28_full_grammar_gene_posterior import oracle


class PooledPosteriorTests(unittest.TestCase):
    def test_independent_complete_set_and_gene_event_masses(self):
        generator = torch.Generator().manual_seed(521)
        for sequence in ('ATGATGTAA', 'ATGAAATAAATGCCCTAA', 'ATGCGTTAGAGTTAGAGCTAGTAA'):
            left = torch.randn(len(sequence),15,generator=generator)
            right = torch.randn(len(sequence),15,generator=generator)
            originals = left.clone(), right.clone()
            pooled = equal_role_energy_pool(left, right)
            self.assertTrue(torch.equal(pooled, (left.double()+right.double())*.5))
            sets, scores, z = oracle(sequence, pooled, 4)
            genes = sorted({gene for selected in sets for gene in selected})
            with patch.multiple(G,_prefix_queries=_scaled_prefix_queries,_suffix_queries=_scaled_suffix_queries):
                for gene in genes:
                    result = G.complete_gene_set_posteriors(sequence, pooled, (gene,), 4)
                    wanted = float(torch.logsumexp(scores[torch.tensor([gene in selected for selected in sets])],0))-z
                    self.assertAlmostEqual(result['log_partition'], z, places=10)
                    self.assertAlmostEqual(result['reverse_log_partition'], z, places=10)
                    self.assertAlmostEqual(result['gene_log_probabilities'][0], wanted, places=10)
                    self.assertAlmostEqual(result['empty_set_log_probability'], -z, places=10)
            self.assertTrue(torch.equal(left, originals[0]) and torch.equal(right, originals[1]))

    def test_promote_before_mean_preserves_half_unit_complete_gene_energy(self):
        sequence = 'ATGAAATAA'
        left = torch.full((9,15),float(2**25),dtype=torch.float32); left[:,1:] += 4
        right = torch.ones(9,15,dtype=torch.float32); right[:,1:] = 2
        pooled = equal_role_energy_pool(left,right)
        wrong = (left+right)*.5
        self.assertEqual(float(pooled[0,1]-pooled[0,0]),2.5)
        self.assertEqual(float(wrong[0,1]-wrong[0,0]),4.)
        with patch.multiple(G,_prefix_queries=_scaled_prefix_queries,_suffix_queries=_scaled_suffix_queries):
            result = G.complete_gene_set_posteriors(sequence,pooled,(((0,9),),),128)
        self.assertEqual(result['chain_gains'][0],22.5)
        z = math.log1p(math.exp(22.5))
        self.assertAlmostEqual(result['log_partition'],z,places=12)
        self.assertAlmostEqual(result['gene_log_probabilities'][0],22.5-z,places=12)
        self.assertFalse(result['posterior_threshold_applied'])
        self.assertFalse(result['observation_model_used'])


if __name__ == '__main__':
    unittest.main()
