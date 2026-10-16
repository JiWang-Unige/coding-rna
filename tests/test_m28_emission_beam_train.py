"""A33 complete partition, fixed integer gates and original raw grammar."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_emission_query_beam_train as A


class EmissionBeamTrainTests(unittest.TestCase):
    def test_four_shards_cover_original_order_once_and_missing_or_duplicate_fails(self):
        ids = ['w'+str(i) for i in range(1491)]
        parts = [[dict(ordinal=i, window=dict(id=ids[i])) for i in range(shard, 1491, 4)]
                 for shard in range(4)]
        self.assertEqual([len(p) for p in parts], [373, 373, 373, 372])
        flat = [r for p in parts for r in p]
        self.assertEqual([r['window']['id'] for r in A.complete_order(flat, ids)], ids)
        for broken in (flat[:-1], flat[:-1]+[flat[0]]):
            with self.assertRaises(ValueError):
                A.complete_order(broken, ids)

    def test_fixed_gate_does_not_confuse_recall_gain_with_precision_or_long_chain_gain(self):
        def row(tp, fp, long):
            return dict(selected_correct=tp, selected_known_nonmatch=fp,
                        selected_five_or_more_correct=long)
        counts = {s: dict(A31_A3=row(10, 5, 2), A33=row(12, 9, 2)) for s in ('a', 'b')}
        references = {s: dict(eligible=20) for s in counts}
        checks = A.gate(counts, references)
        self.assertEqual(len(checks), 5)
        self.assertTrue(checks['a_proxy_F1_strictly_increased'])
        self.assertFalse(checks['a_proxy_precision_not_decreased'])
        self.assertFalse(checks['pooled_five_or_more_selected_correct_strictly_increased'])
        for s in counts:
            counts[s]['A33'] = row(12, 5, 3)
        self.assertTrue(all(A.gate(counts, references).values()))
        counts['a']['A33'] = row(10, 5, 3)
        self.assertFalse(A.gate(counts, references)['a_proxy_F1_strictly_increased'])

    def test_raw_grammar_retains_gc_gt_n_intron_carry_but_rejects_split_stop_and_short_gap(self):
        for intron in ('GTTAG', 'GCTAG', 'GTNAG'):
            A.raw_grammar([[(0, 4), (9, 14)]], 'ATGC'+intron+'TTTAA')
        for chain, sequence in (([(0, 4), (9, 17)], 'ATGTGTTAGAAAAATAA'),
                                ([(0, 4), (8, 13)], 'ATGCGTAGTTTAA')):
            with self.assertRaises(ValueError):
                A.raw_grammar([chain], sequence)


if __name__ == '__main__':
    unittest.main()
