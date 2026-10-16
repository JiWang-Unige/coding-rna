"""Independent short-set and long factorized tests; no research data."""
import math
import unittest
from unittest.mock import patch

import numpy as np
import torch

import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_boundary_compiled import ALPHABET, _state_table
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries as prefix, _scaled_suffix_queries as suffix
from test_m28_full_grammar_gene_posterior import oracle
from test_m28_emission_lattice import direct_gain


def messages(sequence, values, starts, ends, cap, block):
    _,destinations,roles,background = _state_table(cap)
    codes = np.array([ALPHABET.find(x) for x in sequence],dtype=np.int8)
    gain = (values.double()-values[:,:1].double()).numpy()
    a = prefix(codes,gain,np.array(starts,dtype=np.int64),destinations,roles,background,block)
    b = suffix(codes,gain,np.array(ends,dtype=np.int64),destinations,roles,background,block)
    return a,b


class ScaledMassTests(unittest.TestCase):
    def test_independent_all_sets_and_gene_events_at_internal_cuts(self):
        generator = torch.Generator().manual_seed(503)
        for sequence in ('ATGATGTAA','ATGAAATAAATGCCCTAA','ATGCGTNAGTTTAA',
                         'ATGCGTTAGAGTTAGAGCTAGTAA'):
            for cap in (1,2,4):
                values = torch.randn(len(sequence),15,generator=generator).double()
                sets,scores,expected = oracle(sequence,values,cap)
                genes = sorted({g for selected in sets for g in selected})
                starts = sorted({0,*(g[0][0] for g in genes)})
                ends = sorted({len(sequence),*(g[-1][1] for g in genes)})
                for block in (1,3,7,256):
                    (z,a),(rz,b) = messages(sequence,values,starts,ends,cap,block)
                    self.assertAlmostEqual(z,expected,places=10)
                    self.assertAlmostEqual(rz,expected,places=10)
                    aa,bb = dict(zip(starts,a)),dict(zip(ends,b))
                    for gene in genes:
                        keep = torch.tensor([gene in selected for selected in sets])
                        wanted = float(torch.logsumexp(scores[keep],0))-expected
                        actual = math.fsum((aa[gene[0][0]],direct_gain(values,gene),bb[gene[-1][1]],-z))
                        self.assertAlmostEqual(actual,wanted,places=10)

    def test_query_and_scale_blocks_do_not_close_gene_paths(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAAATGAAATAA'
        values = ((torch.arange(len(sequence)*15).reshape(-1,15)%11)-5).double()/8
        reference = messages(sequence,values,[0,3,6,11,19],[3,6,11,19,len(sequence)],4,256)
        for block in (1,2,3,5,7,8192):
            actual = messages(sequence,values,[0,3,6,11,19],[3,6,11,19,len(sequence)],4,block)
            for (z,a),(rz,b) in zip(reference,actual):
                self.assertAlmostEqual(z,rz,places=10)
                self.assertTrue(np.allclose(a,b,rtol=0,atol=1e-10))
        shifted = messages(sequence,values+torch.arange(len(sequence))[:,None]*8,
            [0,3,6,11,19],[3,6,11,19,len(sequence)],4,3)
        self.assertAlmostEqual(reference[0][0],shifted[0][0],places=10)

    def test_long_analytic_product_of_independent_gene_choices(self):
        # N separates literal ORFs, so all legal sets are independent binary choices.
        count = 100000; sequence = 'ATGAAATAAN'*count; n = len(sequence)
        _,destinations,roles,background = _state_table(1)
        codes = np.tile(np.array([ALPHABET.find(x) for x in 'ATGAAATAAN'],dtype=np.int8),count)
        gain = np.zeros((n,15),dtype=np.float64)
        low = float(np.float32(.1)); gain[:n//2,1:] = low; gain[n//2:,1:] = 1.
        low_mass,high_mass = math.log1p(math.exp(9*low)),math.log1p(math.exp(9.))
        expected = math.fsum([low_mass]*(count//2)+[high_mass]*(count//2))
        starts = np.array([0,n//2,n-10],dtype=np.int64); ends = np.array([9,n//2+9,n-1,n],dtype=np.int64)
        z,a = prefix(codes,gain,starts,destinations,roles,background)
        rz,b = suffix(codes,gain,ends,destinations,roles,background)
        self.assertAlmostEqual(z,expected,delta=5e-8)
        self.assertAlmostEqual(rz,expected,delta=5e-8)
        self.assertTrue(math.isclose(z,rz,rel_tol=1e-12,abs_tol=1e-8))
        for x,y,energy,mass in ((0,9,9*low,low_mass),(n//2,n//2+9,9.,high_mass),(n-10,n-1,9.,high_mass)):
            event = math.fsum((a[np.searchsorted(starts,x)],energy,b[np.searchsorted(ends,y)],-z))
            self.assertAlmostEqual(event,energy-mass,delta=5e-8)

    def test_original_posterior_API_has_unchanged_input_energy_and_probability_rules(self):
        sequence = 'ATGATGTAA'; values = torch.zeros(9,15,dtype=torch.float64)
        original = values.clone()
        with patch.multiple(G,_prefix_queries=prefix,_suffix_queries=suffix):
            result = G.complete_gene_set_posteriors(sequence,values,(((0,9),),),128)
            self.assertAlmostEqual(result['gene_probabilities'][0],1/3,places=12)
            self.assertFalse(result['posterior_threshold_applied'])
            self.assertFalse(result['Viterbi_winner_pruning'])
            self.assertEqual(result['optimizer_steps'],0)
        self.assertTrue(torch.equal(values,original))


if __name__=='__main__':
    torch.set_num_threads(2)
    unittest.main()
