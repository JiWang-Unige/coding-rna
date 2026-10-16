import copy
import unittest

from build_natural_competitor_support import D, build_records


def fixture():
    species = 'arabidopsis_thaliana'
    values = dict(id='artificial', species=species, seqid=next(iter(D.R.C.R.ALLOW[species])),
                  split='train', strand='+', start=0, end=60,
                  tile_bases=24576, valid_bases=60, pad_bases=24516)
    window = {k: values.get(k, 0) for k in D.R.C.R.GEOMETRY}
    chains = [[[0, 12]], [[20, 26]], [[30, 36]], [[40, 46]]]
    labels = [1, 0, -1, 1]
    old = dict(window={**window, 'chain_unknown_intervals': []},
               candidates=[dict(window_id=window['id'], species=species, chain=c, label=y,
                                origin='free' if i < 3 else 'reference_injected', old_model_gain=99.)
                           for i, (c, y) in enumerate(zip(chains, labels))],
               counts=dict(windows=1, free_candidates=3, reference_injected_candidates=1,
                           zero_known_windows=0, label_1=2, label_0=1, **{'label_-1': 1}))
    prior = dict(window=window, reference_used_in_scoring=False,
                 variants={'A31_A3': dict(chains=chains[:3])})
    natural = dict(window=window, original_ordinal=0, reference_used_in_scoring=False,
                   chains=[[[0, 6]], [[12, 18]], [[40, 46]]], selected=[0, 1, 2])
    labelled = dict(window=window, labels={'A31_A3': labels[:3], 'A42': [0, -1, 1]},
                    eligible_references=2, eligible_five_or_more=0)
    draws = [dict(window_id=window['id'], species=species,
                  importance_weight_uniform_within_stratum=1.25) for _ in range(2)]
    receipts = [dict(step=i+1, epoch=1, draw_ordinal=i, window_id=window['id'],
                     species=species, importance_weight=1.25) for i in range(2)]
    return draws, receipts, [old], [prior], [natural], [labelled]


class NaturalCompetitorSupportTests(unittest.TestCase):
    def test_injected_gold_bridge_and_repeated_draws_preserve_original_supervision(self):
        args = fixture(); before = copy.deepcopy(args)
        records, stats = build_records(*args)
        self.assertEqual(args, before)
        self.assertEqual((stats['windows'], stats['original_draws'],
                          stats['original_candidate_window_pairs'], stats['candidate_window_pairs']), (1, 2, 4, 6))
        self.assertEqual(stats['label_counts_pooled'], {'1': 2, '0': 2, '-1': 2})
        self.assertEqual(stats['added_label_counts_pooled'], {'0': 1, '-1': 1})
        self.assertEqual(stats['A42_correct_matches_original_reference_injected_per_species']
                         ['arabidopsis_thaliana'], 1)
        self.assertEqual(stats['original_replayed_draw_counts'], dict(windows=2, free_candidates=6,
            reference_injected_candidates=2, zero_known_windows=0, label_1=4, label_0=2, **{'label_-1': 2}))
        self.assertEqual(stats['merged_replayed_draw_counts'], dict(windows=2, original_free_candidates=6,
            reference_injected_candidates=2, added_A42_natural_candidates=4, zero_known_windows=0,
            label_1=4, label_0=4, **{'label_-1': 4}))
        candidates = records[0]['candidates']
        self.assertEqual([r['chain'] for r in candidates if r['label'] == 1], [[[0, 12]], [[40, 46]]])
        self.assertEqual([r['origin'] for r in candidates], ['free']*3+['reference_injected']+['A42_natural']*2)
        self.assertEqual([r['natural_index'] for r in candidates[4:]], [0, 1])
        self.assertNotIn('chain_unknown_intervals', records[0]['window'])
        self.assertTrue(all('old_model_gain' not in r for r in candidates))

    def test_joint_source_and_draw_conflicts_rejected(self):
        cases = []
        args = fixture(); args[1][0]['importance_weight'] = 1.; cases.append(args)
        args = fixture(); args[1][1]['draw_ordinal'] = 0; cases.append(args)
        args = fixture(); args[3].clear(); cases.append(args)
        args = fixture(); args[2][0]['window']['start'] = 1; cases.append(args)
        args = fixture(); args[2][0]['candidates'][3]['label'] = 0; cases.append(args)
        args = fixture(); args[5][0]['labels']['A42'][2] = 0; cases.append(args)
        args = fixture(); args[4][0]['chains'][2] = [[48, 54]]; cases.append(args)
        args = fixture(); args[5][0]['eligible_references'] = 3; cases.append(args)
        for i, args in enumerate(cases):
            with self.subTest(case=i), self.assertRaises(ValueError):
                build_records(*args)

    def test_empty_original_draw_is_retained_when_natural_background_is_added(self):
        draws, receipts, originals, previous, natural, labelled = fixture()
        originals[0]['candidates'] = []
        originals[0]['counts'] = dict(windows=1, free_candidates=0, reference_injected_candidates=0,
                                      zero_known_windows=1, label_1=0, label_0=0, **{'label_-1': 0})
        previous[0]['variants']['A31_A3']['chains'] = []
        natural[0]['chains'] = [[[0, 6]]]; natural[0]['selected'] = [0]
        labelled[0]['labels'] = {'A31_A3': [], 'A42': [0]}; labelled[0]['eligible_references'] = 0
        records, stats = build_records(draws, receipts, originals, previous, natural, labelled)
        self.assertEqual(len(records), 1)
        self.assertEqual(stats['original_replayed_draw_counts']['zero_known_windows'], 2)
        self.assertEqual(stats['merged_replayed_draw_counts']['zero_known_windows'], 0)
        self.assertEqual(stats['merged_replayed_draw_counts']['added_A42_natural_candidates'], 2)
        self.assertEqual(stats['label_counts_pooled'], {'0': 1})


if __name__ == '__main__':
    unittest.main()
