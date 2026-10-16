"""Tests that affect M28 sampling, coordinates and annotation policy."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/experiments/M28-METHOD-RESTART"))
import manifest as M

class ManifestTests(unittest.TestCase):
    def test_tail_and_short_padding_orientation(self):
        self.assertEqual(M.starts(31,12,6),[0,6,12,18,19])
        self.assertEqual(M.starts(5,12,6),[0])
        self.assertEqual(M.orient_interval(2,5,0,8,"-"),(3,6))
        self.assertEqual(M.rc("ATGCCN"),"NGGCAT")

    def test_gene_first_marginal_and_importance(self):
        rows=[{"id":"a","positive_ids":["x","y"]},{"id":"b","positive_ids":["y"]},
              {"id":"c","positive_ids":[]}]
        q=M.probabilities(rows)
        np.testing.assert_allclose(q,[.4375,.3125,.25])
        self.assertAlmostEqual(float(np.sum(q/(len(q)*q))),1.)
        draws=M.sample_stratum(rows,32,np.random.default_rng(0))
        self.assertTrue(all(d["stratum_probability"]>0 for d in draws))
        np.testing.assert_allclose(M.probabilities([{"positive_ids":[]},{"positive_ids":[]}]),[.5,.5])

    def test_length_allocation_and_no_prefix(self):
        self.assertEqual(M.allocate({"a":1,"b":3},8),{"a":2,"b":6})
        self.assertEqual(sum(M.allocate({"a":1,"b":1,"c":1}).values()),384)
        rows=[{"id":str(i),"positive_ids":[]} for i in range(100)]
        draws=M.sample_stratum(rows,32,np.random.default_rng(0))
        self.assertGreater(max(int(d["window_id"]) for d in draws),90)

    def test_spliced_phase_and_reverse(self):
        dna="ATGAAGTCCAGACCCTAA"
        tx={"CDS":[(0,5,"0"),(11,18,"1")],"strand":"+"}
        self.assertEqual(M.chain_checks(tx,dna),([],{"GT-AG":1}))
        reverse={"CDS":[(0,7,"1"),(13,18,"0")],"strand":"-"}
        self.assertEqual(M.chain_checks(reverse,M.rc(dna)),([],{"GT-AG":1}))
        tx["CDS"][1]=(11,18,"0")
        self.assertIn("GFF_phase_disagreement",M.chain_checks(tx,dna)[0])

    def test_sequence_completeness_not_parent_flag(self):
        tx={"CDS":[(0,9,"0")],"strand":"+","partial":True}
        self.assertEqual(M.chain_checks(tx,"ATGAAATAA")[0],[])
        self.assertIn("internal_stop",M.chain_checks({"CDS":[(0,12,"0")],"strand":"+"},"ATGTAACCCTAA")[0])
        self.assertEqual(M.merged([(1,4),(3,6),(8,9)]),[[1,6],[8,9]])

if __name__=="__main__": unittest.main()
