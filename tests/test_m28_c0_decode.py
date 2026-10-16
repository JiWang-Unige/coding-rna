import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.c0_decode import native_probabilities,decode_c0,genomic_chain,owner_window
from src.m28.labels import build_targets

def example(cut=12,motif="GT"):
    coding="ATG"+"AAA"*12+"TAA"
    intron=motif+"C"*24+"AG"
    dna="CCC"+coding[:cut]+intron+coding[cut:]+"CCC"
    tx={"id":"p","gene_id":"g","CDS":[(3,3+cut,"0"),(3+cut+len(intron),3+len(coding)+len(intron),str((3-cut%3)%3))],
        "training_primary":True,"CDS_complete_supported":True,"strand":"+"}
    w={"start":0,"end":len(dna),"valid_bases":len(dna),"width":len(dna),"strand":"+",
       "overlapping_transcript_ids":["p"],"positive_ids":["p"],"fine_label_unknown_intervals":[],
       "chain_unknown_intervals":[],"no_annotated_gene_either_strand":False}
    labels=build_targets(w,{"p":tx},dna)["seg"]
    probabilities=np.full((len(dna),15),1e-6,dtype=np.float32)
    probabilities[np.arange(len(dna)),labels]=1-14e-6
    return dna,probabilities,[[(a,b) for a,b,_ in tx["CDS"]]]

class NativeBridgeTests(unittest.TestCase):
    def test_explicit_phase_permutation(self):
        p=np.arange(15,dtype=np.float32)[None]
        converted=native_probabilities(p)
        self.assertEqual(converted.tolist(),[[0,1,3,2,4,6,5,7,9,8,10,12,11,13,14]])
        np.testing.assert_array_equal(p,np.arange(15,dtype=np.float32)[None])

    def test_native_all_splice_phases_gt_gc(self):
        for motif in ("GT","GC"):
            for phase in (0,1,2):
                with self.subTest(motif=motif,phase=phase):
                    dna,p,expected=example(12+phase,motif)
                    before=p.copy()
                    result=decode_c0(dna,p)
                    self.assertEqual(result["chains"],expected)
                    self.assertEqual(result["partial_paths"],0)
                    np.testing.assert_array_equal(p,before)

    def test_background_and_multiple_genes(self):
        dna,p,expected=example()
        background=np.full((12,15),1e-6,dtype=np.float32);background[:,0]=1-14e-6
        self.assertEqual(decode_c0("C"*12,background)["chains"],[])
        combined=dna+"C"*12+dna
        probs=np.concatenate([p,background,p])
        shift=len(dna)+12
        second=[(a+shift,b+shift) for a,b in expected[0]]
        self.assertEqual(decode_c0(combined,probs)["chains"],expected+[second])

    def test_strand_and_tail_core(self):
        chain=[(3,15),(43,73)]
        self.assertEqual(genomic_chain(chain,100,180,"-"),[(107,137),(165,177)])
        self.assertEqual(genomic_chain(chain,100,180,"+"),[(103,115),(143,173)])
        windows=[(0,24),(12,36),(17,41)]
        self.assertEqual(owner_window([(20,40)],windows),2)
        # If any same-width window contains a chain, nearest-center ownership also contains it.
        for a in range(41):
            for b in range(a+1,42):
                if any(x<=a and b<=y for x,y in windows):
                    x,y=windows[owner_window([(a,b)],windows)]
                    self.assertTrue(x<=a and b<=y)

if __name__=="__main__": unittest.main()
