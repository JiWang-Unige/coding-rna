import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from analyze_query_beam_train import compare_references,aggregate,label_counts


class TrainedBeamCountsTests(unittest.TestCase):
    def test_fixed_identity_denominator_and_four_states(self):
        refs=[{'reference_id':str(i),'chain':[(i*10,i*10+6)]} for i in range(4)]
        g={((0,6),),((20,26),)};b={((0,6),),((10,16),)}
        counts,rows=compare_references(refs,g,b)
        self.assertEqual(counts['reference_window_pairs'],4)
        self.assertEqual(counts['states'],dict(both=1,beam_only=1,greedy_only=1,neither=1))
        self.assertEqual(counts,aggregate(rows))
        self.assertEqual((counts['greedy_free_hits'],counts['beam_free_hits']),(2,2))

    def test_equal_geometry_identities_and_empty_windows_are_not_filtered(self):
        refs=[{'reference_id':k,'chain':[(0,6)]} for k in ('a','b')]
        counts,rows=compare_references(refs,{((0,6),)},set())
        self.assertEqual(counts['reference_window_pairs'],2)
        self.assertEqual(counts['states']['greedy_only'],2)
        self.assertEqual([r['reference_id'] for r in rows],['a','b'])
        counts,rows=compare_references([],set(),{((0,6),)})
        self.assertEqual(counts,aggregate(rows))
        self.assertEqual(counts['reference_window_pairs'],0)

    def test_protected_and_unknown_candidates_are_not_false_negative_labels(self):
        chains={((0,6),),((10,16),),((20,26),),((30,36),)}
        targets={'positive_chains':[[(0,6)]],'protected_chains':[[(10,16)]],
                 'unknown_chain_spans':[(20,26)],'callable_bases':np.ones(40,dtype=bool)}
        self.assertEqual(label_counts(chains,targets),{'1':1,'0':1,'-1':2})


if __name__=='__main__':
    unittest.main()
