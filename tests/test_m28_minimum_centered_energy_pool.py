"""Independent fixed local-disagreement pooling and whole-set MAP checks."""
import unittest

import torch
from src.m28.full_grammar_energy_pool import minimum_centered_role_energy_pool as pool
from src.m28.full_grammar_skyline import decode_full_grammar_skyline as decode
from test_m28_full_grammar_energy_pool import all_sets
from test_m28_emission_lattice import direct_gain, prefer_chain


class MinimumCenteredEnergyTests(unittest.TestCase):
    def test_complete_sets_match_independent_local_minimum_and_parent_bound(self):
        rng=torch.Generator().manual_seed(251)
        for dna in ('ATGAAATAAATGCCCTAA','ATGCGTNAGTTTAA','ATGTGTTAGATGCGTTAGAATAA'):
            left=torch.randn(len(dna),15,generator=rng); right=torch.randn(len(dna),15,generator=rng)
            manual=torch.tensor([[min(float(left[i,j])-float(left[i,0]),
                float(right[i,j])-float(right[i,0])) for j in range(15)] for i in range(len(dna))],dtype=torch.float64)
            actual=pool(left,right)
            self.assertTrue(torch.equal(actual,manual))
            for cap in (1,2,4):
                legal=all_sets(dna,cap)
                gains=[sum(direct_gain(manual,c) for c in genes) for genes in legal]
                for genes,gain in zip(legal,gains):
                    for branch in (left,right):
                        self.assertLessEqual(gain,sum(direct_gain(branch,c) for c in genes)+1e-12)
                decoded=decode(dna,actual,cap,3)
                self.assertIn(decoded['chains'],legal)
                self.assertAlmostEqual(decoded['gain'],max(gains),places=9)

    def test_identity_symmetry_independent_row_gauge_and_immutability(self):
        left=(torch.arange(270).reshape(18,15)%11-5).float()
        right=(torch.arange(270).reshape(18,15)%7-3).float()
        a,b=left.clone(),right.clone(); result=pool(left,right)
        self.assertTrue(torch.equal(result,pool(right,left)))
        self.assertTrue(torch.equal(left,a)); self.assertTrue(torch.equal(right,b))
        shift=torch.arange(18)[:,None]*8
        self.assertTrue(torch.equal(result,pool(left+shift,right-2*shift)))
        self.assertTrue(torch.equal(result[:,0],torch.zeros(18,dtype=torch.float64)))
        centered=left.double()-left[:,0:1].double()
        self.assertTrue(torch.equal(pool(left,left),centered))
        # A role-local minimum is intentionally stricter than a whole-path min.
        x=torch.zeros(2,15); y=x.clone(); x[:,1]=torch.tensor([2.,-2.]); y[:,1]=-x[:,1]
        self.assertEqual(float(pool(x,y)[:,1].sum()),-4.)
        self.assertEqual(min(float(x[:,1].sum()),float(y[:,1].sum())),0.)

    def test_literal_carry_and_exon_budget_unchanged(self):
        dna='ATGTGTTAGATGCGTTAGAATAA'; chain=((9,13),(18,23))
        scores=prefer_chain(dna,chain).float()
        scores[:3,13]=1000.; scores[3,7]=30.; scores[6:9,14]=-10000.; scores[9,11]=-10000.
        self.assertEqual(decode(dna,pool(scores,scores),3,3)['chains'],(chain,))
        dna='ATGCGTTAGTTTAA'; chain=((0,4),(9,14)); scores=prefer_chain(dna,chain).float()
        self.assertEqual(decode(dna,pool(scores,scores),1,3)['chains'],())
        self.assertEqual(decode(dna,pool(scores,scores),2,3)['chains'],(chain,))


if __name__=='__main__':
    torch.set_num_threads(2)
    unittest.main()
