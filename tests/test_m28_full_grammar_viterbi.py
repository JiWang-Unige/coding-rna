"""New max/checkpoint recurrence vs independent complete gene-set oracle."""
import unittest
import torch
from src.m28.full_grammar_viterbi import decode_full_grammar_viterbi as decode
from test_m28_emission_lattice import oracle_gain,direct_gain,prefer_chain


class FullGrammarViterbiTests(unittest.TestCase):
    def test_complete_set_optimum_actual_carry_and_selected_scores_match_oracle(self):
        generator=torch.Generator().manual_seed(211)
        for sequence in ('ATGAAATAAATGCCCTAA','ATGCGTNAGTTTAA','ATGGCTAGAAATAA',
                'ATGCGTTAGAGTTAGAGCTAGTAA','ATGTGTTAGATGCGTTAGAATAA'):
            for cap in (1,2,4):
                for values in (torch.randint(-3,4,(len(sequence),15),generator=generator).float(),
                        torch.randn(len(sequence),15,generator=generator).double()):
                    expected,legal=oracle_gain(sequence,values,cap)
                    actual=decode(sequence,values,cap,checkpoint_bp=3)
                    self.assertAlmostEqual(actual['gain'],expected,places=9)
                    last=0
                    for chain,gain in zip(actual['chains'],actual['chain_gains']):
                        self.assertIn(chain,legal); self.assertGreaterEqual(chain[0][0],last)
                        self.assertAlmostEqual(gain,direct_gain(values,chain),places=9)
                        last=chain[-1][1]
                    self.assertFalse(actual['reference_used'])

    def test_checkpoint_spans_preserve_scores_roles_and_selected_gene_sets(self):
        sequence='ATGCGTTAGAGTTAGAGCTAGTAA'
        values=(torch.arange(len(sequence)*15).reshape(-1,15).double()%11-5)/7
        reference=None
        for span in (1,2,3,5,7,8192):
            actual=decode(sequence,values,4,span)
            key=(actual['chains'],actual['gain'],actual['chain_gains'])
            if reference is None: reference=key
            self.assertEqual(key,reference)
            s=actual['counts']; rows=len([0]+list(range(span,len(sequence),span))+[len(sequence)])
            self.assertEqual(s['checkpoint_bytes'],rows*s['state_count']*8)
            self.assertEqual(s['maximum_block_predecessor_bytes'],min(span,len(sequence))*s['state_count']*4)

    def test_zero_null_row_shift_and_adjacent_complete_genes(self):
        sequence='ATGAAATAAATGCCCTAA'; values=torch.zeros(18,15)
        result=decode(sequence,values,128,7)
        self.assertEqual((result['chains'],result['gain'],result['chain_gains']),((),0.,()))
        values[:,0]=-1.; values[:,13]=2.; values[:,14]=2.
        first=decode(sequence,values,128,7)
        shifted=decode(sequence,values+torch.arange(18)[:,None]*8,128,7)
        self.assertEqual((first['chains'],first['gain'],first['chain_gains']),
            (shifted['chains'],shifted['gain'],shifted['chain_gains']))
        self.assertEqual(first['chains'],(((0,9),),((9,18),)))

    def test_original_128_exon_capacity_and_cap127_exclusion(self):
        pieces=['ATG']+['AAA']*126+['TAA']; sequence='GTCAG'.join(pieces)
        chain=tuple((8*i,8*i+3) for i in range(128)); values=prefer_chain(sequence,chain).float()
        full=decode(sequence,values,128,31)
        self.assertEqual(full['chains'],(chain,)); self.assertEqual(full['gain'],len(sequence)*10.)
        self.assertNotIn(chain,decode(sequence,values,127,31)['chains'])


if __name__=='__main__':
    torch.set_num_threads(2); unittest.main()
