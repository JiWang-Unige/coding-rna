import json
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from diagnose_pool_ceiling import ceiling


def reference(seqid, strand, chain, available=True):
    return {'species':'plant','seqid':seqid,'strand':strand,
            'CDS_intervals':json.dumps(chain),'A2_free_unique':str(available)}


class PoolCeilingTests(unittest.TestCase):
    def test_same_strand_conflict_does_not_merge_distinct_contigs_or_opposite_strands(self):
        records = [reference('chr1','+',[(0,12)]),reference('chr1','+',[(6,18)]),
                   reference('chr1','-',[(0,12)]),reference('chr2','+',[(0,12)])]
        result = ceiling(records,'A2','free_unique')
        self.assertEqual((result['reference_chains'],result['available_exact_reference_chains'],
                          result['maximum_nonoverlapping_TP'],result['overlap_conflict_loss']),(4,4,3,1))
        self.assertAlmostEqual(result['F1_ceiling_with_reference_known_perfect_reranking'],6/7)

    def test_missing_reference_stays_in_full_F1_denominator(self):
        records = [reference('chr1','+',[(0,12)]),reference('chr1','+',[(20,32)],False)]
        result = ceiling(records,'A2','free_unique')
        self.assertEqual(result['maximum_nonoverlapping_TP'],1)
        self.assertAlmostEqual(result['F1_ceiling_with_reference_known_perfect_reranking'],2/3)
        records[1]['A2_free_unique'] = 'unknown'
        with self.assertRaises(ValueError):
            ceiling(records,'A2','free_unique')


if __name__ == '__main__':
    unittest.main()
