import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts/experiments/M28-METHOD-RESTART'))
from src.m28.candidates import Budget,generate
from r6_shadow import ReferenceTrace,generate_shadow

def anchor(chain):
    return {'reference_id':'fixed','species':'test','strand':'+','stratum':'selected_exact',
            'oriented_chain':chain}

def replay(dna,ep,cp,chain,budget=Budget()):
    observer=ReferenceTrace([anchor(chain)])
    original=generate(dna,ep,cp,budget=budget)
    shadow=generate_shadow(dna,ep,cp,budget=budget,observer=observer)
    if original!=shadow: raise AssertionError('Observer changed generator')
    return observer.results()[0]

class TraceTests(unittest.TestCase):
    def test_carry_and_motifs_replay(self):
        for motif in ('GT','GC'):
            for phase in range(3):
                cut=12+phase;coding='ATG'+'AAA'*12+'TAA';intron=motif+'C'*24+'AG'
                dna='CCC'+coding[:cut]+intron+coding[cut:]+'CCC'
                chain=[(3,3+cut),(3+cut+len(intron),3+len(coding)+len(intron))]
                ep=np.full((len(dna),4),-8.)
                ep[3,0]=8;ep[chain[-1][1]-1,1]=8;ep[chain[0][1]-1,2]=8;ep[chain[1][0],3]=8
                result=replay(dna,ep,np.full(len(dna),.8),chain)
                self.assertIsNone(result['first_actual_cut'])
                self.assertTrue(result['final']['survived'])
                self.assertTrue(all(x['retained'] for x in result['conditional_reference_edges']))

    def test_donor_first_cut_same_carry_and_score_components(self):
        dna='ATGATGATGAAA'+'GT'+'C'*7+'AG'+'AAATAA'
        ep=np.zeros((len(dna),4));ep[3,0]=3;ep[6,0]=6
        result=replay(dna,ep,np.full(len(dna),.6),[(0,12),(23,29)],Budget(beam_per_carry=1))
        self.assertEqual(result['first_actual_cut']['stage'],'donor_beam')
        detail=result['first_actual_cut']['detail']
        self.assertEqual(detail['rank'],3);self.assertAlmostEqual(detail['cutoff_gap'],6.)
        self.assertEqual(detail['surviving_competitors'][0]['prefix'],[[6,12]])
        self.assertTrue(result['mechanism_witness'])
        later=next(x for x in result['actual_beam_events'] if x['stage']=='entry' and x['exon_index']==1)
        self.assertFalse(later['present_before_beam'])
        self.assertFalse(result['conditional_reference_edges'][1]['actual_execution_claim'])

    def test_final_cut_competitors(self):
        dna=('ATGAAATAA'+'CCCC')*8
        ep=np.zeros((len(dna),4))
        for j in range(8): ep[j*13+8,1]=j
        result=replay(dna,ep,np.full(len(dna),.6),[(0,9)],Budget(sites_per_bin=64,ends_per_start_kind=64,chains_per_kb=1))
        self.assertEqual(result['first_actual_cut']['stage'],'final_chain_cut')
        detail=result['first_actual_cut']['detail']
        self.assertEqual(detail['rank'],8);self.assertAlmostEqual(detail['cutoff_gap'],7.)
        self.assertEqual(detail['surviving_competitors'][0]['prefix'],[[91,100]])
        self.assertTrue(result['mechanism_witness'])

    def test_single_exon_empty_links_does_not_prove_endpoint_support(self):
        dna='ATGAAATAA'+'CCCC'+'ATGAAATAA'
        ep=np.zeros((len(dna),4));ep[13,0]=4
        b=Budget(sites_per_bin=1)
        result=replay(dna,ep,np.full(len(dna),.6),[(0,9)],b)
        self.assertEqual(result['necessary_links'],[])
        self.assertEqual(result['first_actual_cut']['stage'],'start_endpoint')
        self.assertFalse(result['mechanism_witness'])

    def test_random_pruning_replay_and_ineligible_reference(self):
        dna=('ATGAAATAA'+'GTCCCCAG')*12
        rng=np.random.default_rng(19)
        for j in range(5):
            b=Budget(bin_bp=64,sites_per_bin=2,ends_per_start_kind=2,
                     intron_links_per_donor=2,beam_per_carry=1,chains_per_kb=1)
            replay(dna,rng.normal(size=(len(dna),4)),rng.uniform(.1,.9,len(dna)),[(0,9)],b)

if __name__=='__main__': unittest.main()
