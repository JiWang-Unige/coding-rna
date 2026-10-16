import copy
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import readout_logistic_context_train as D


class LogisticReadoutTests(unittest.TestCase):
    def fixture(self):
        chains = [[[0,3],[6,9],[12,15],[18,21],[24,27]], [[0,27]], [[30,36]]]
        labels = [1,0,-1]
        def arm(scores):
            return dict(chains=copy.deepcopy(chains),scores=scores,
                        selected=D.R.select_nonoverlapping(chains,scores))
        rows = [[],[],[],[]]
        for species in D.R.C.R.ALLOW:
            w = dict(id=species+'|toy',species=species,split='train')
            old, context, new = arm([4.,1.,2.]),arm([-1.,-2.,-3.]),arm([5.,.2,1.])
            oldcounts = D.R.selected_counts(old,labels)
            rows[0].append(dict(window=w,reference_used_in_scoring=False,**new))
            rows[1].append(dict(window=w,reference_used_in_scoring=False,variants=dict(A31_A3=old)))
            rows[2].append(dict(window=w,labels=dict(A31_A3=labels),counts=dict(A31_A3=oldcounts),
                                eligible_references=2,eligible_five_or_more=1))
            rows[3].append(dict(window=w,reference_used_in_scoring=False,**context))
        return rows

    def test_three_arms_unknown_and_full_reference_denominators(self):
        rows, counts, references = D.paired_readout(*self.fixture())
        self.assertEqual(len(rows),2)
        for species,value in counts.items():
            self.assertEqual(value['A38']['selected_correct'],1)
            self.assertEqual(value['A38']['selected_unknown'],1)
            self.assertEqual(value['A38']['selected_known_nonmatch'],0)
            self.assertEqual(value['A38']['selected_five_or_more_correct'],1)
            self.assertEqual(value['A37']['selected'],0)
            self.assertEqual(references[species],dict(eligible=2,five_or_more=1))
        checks = D.gate(counts)
        self.assertEqual(len(checks),5); self.assertEqual(sum(checks.values()),2)
        # A38 beating the weaker A37 is not sufficient against the original stronger arm.

    def test_missing_or_reordered_comparator_rejected(self):
        rows = self.fixture(); rows[3].pop()
        with self.assertRaises(ValueError): D.paired_readout(*rows)
        rows = self.fixture(); rows[3][0]['chains'] = rows[3][0]['chains'][::-1]
        with self.assertRaises(ValueError): D.paired_readout(*rows)


if __name__ == '__main__':
    unittest.main()
