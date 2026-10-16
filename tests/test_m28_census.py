import importlib.util
import unittest
from pathlib import Path
P = Path(__file__).resolve().parents[1] / "scripts/experiments/M28-METHOD-RESTART/census.py"
spec = importlib.util.spec_from_file_location("m28census", P)
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

class CensusTests(unittest.TestCase):
    def test_half_open_exact_and_boundary(self):
        self.assertTrue(M.contained(0,10,10,10,20))
        self.assertFalse(M.contained(5,15,10,10,20))
        self.assertTrue(M.contained(5,15,10,5,20))
        self.assertFalse(M.contained(0,11,10,5,20))
        self.assertFalse(M.contained(19,20,10,10,19))
    def test_sorted_sample_then_prefix_bias(self):
        rows = [("A", i) for i in range(1000)] + [("B", i) for i in range(1000)]
        old = M.old_sample(rows,.12,40)
        new = M.uniform_cap(rows,40)
        self.assertEqual({r[0] for r in old}, {"A"})
        self.assertEqual({r[0] for r in new}, {"A","B"})
        self.assertEqual(len(new),40)
        self.assertEqual(new,M.uniform_cap(rows,40))
    def test_overlap_strand_and_touch(self):
        rows = [{"id":str(i),"strand":s,"CDS":[(a,b,"0")]} for i,(a,b,s) in enumerate([(0,10,"+"),(10,20,"+"),(8,12,"+"),(0,30,"-")])]
        self.assertEqual(M.overlapping_ids(rows),{"0","1","2"})
    def test_union(self):
        self.assertEqual(M.union_length([(0,10),(5,15),(20,21)]),16)
