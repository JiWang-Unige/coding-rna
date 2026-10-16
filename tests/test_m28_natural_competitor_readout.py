import copy
import unittest
from unittest.mock import patch
import numpy as np
import torch

import readout_natural_competitor_train as D
from verify_natural_competitor_train import verify_training_rows


class NaturalCompetitorReadoutTests(unittest.TestCase):
    def test_prediction_inputs_are_original_geometry_only_not_training_gold(self):
        species = sorted(D.R.C.R.ALLOW)[0]
        def window(ident, start, strand):
            return dict(id=ident, species=species, seqid=sorted(D.R.C.R.ALLOW[species])[0],
                split='train', strand=strand, start=start, end=start+30, width=24576,
                valid_bases=30, pad_right_after_orientation=24546)
        a, b = window('a', 0, '+'), window('b', 100, '-')
        draws = [dict(window_id=w['id'], species=species) for w in (b, a, b)]
        entries = [dict(window=w, artifact='artificial/'+w['id']) for w in (a, b)]
        with patch.object(D.A42.B, 'population', side_effect=AssertionError('Prediction must not open a candidate pool')):
            windows, index = D.ordered_prediction_inputs(draws, entries)
        self.assertEqual(windows, [dict(window=b), dict(window=a)])
        self.assertEqual(set(index), {'a', 'b'})
        for changed in ([entries[0]], entries+[entries[0]],
                        [dict(entries[0], window=dict(a, positive_ids=['gold'])), entries[1]]):
            with self.subTest(entries=changed), self.assertRaises(ValueError): D.ordered_prediction_inputs(draws, changed)

    def test_only_actual_A43_last_fit_header_accepted(self):
        value = dict(arm='A43', configuration=copy.deepcopy(D.F.CONFIG), protocol=D.F.PROTOCOL,
            optimizer_steps=4608, feature_revision=D.R.C.I.F.REVISION, completed=True,
            source_support_job='13424089', preflight_job='13424115', final_checkpoint=D.CONFIG['source_checkpoint'],
            trainable_parameters=366191, initialization='fresh_seed0_no_supervised_weight_inheritance',
            **{k: False for k in ('supervised_checkpoint_read', 'new_R4_labelling', 'generator_updated',
                'new_backbone_forward', 'DEV_used_for_selection', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened')})
        self.assertIsNone(D.validate_fit(value))
        for key, wrong in (('arm', 'A36'), ('optimizer_steps', 1536), ('source_support_job', 'wrong'),
                           ('final_checkpoint', 'early.pt'), ('DEV_used_for_selection', True)):
            changed = dict(value, **{key: wrong})
            with self.subTest(key=key), self.assertRaises(ValueError): D.validate_fit(changed)

    def test_natural_score_check_has_no_candidate_or_label_input_and_retains_null(self):
        local = torch.zeros(6, 15); local[:3, 13] = 1.
        value = dict(sequence='ATGTAA', local_logits=local)
        natural = dict(chains=[[[0, 6]]], chain_gains=[3.], gain=3., reference_used=False)
        with patch.object(D, 'candidate_targets', side_effect=AssertionError('No prediction labels')):
            self.assertIsNone(D.validate_prediction(value, natural))
            self.assertIsNone(D.validate_prediction(value, dict(chains=[], chain_gains=[], gain=0., reference_used=False)))
        with self.assertRaises(ValueError): D.validate_prediction(value, dict(natural, reference_used=True))

    def test_new_geometry_R4_labels_A43_unknown_and_full_denominators(self):
        w = dict(id='toy', species=sorted(D.R.C.R.ALLOW)[0], valid_bases=50)
        prior = dict(window=w, variants={'A31_A3': dict(chains=[[[0, 6]]], scores=[2.], selected=[0])})
        row = dict(window=w, chains=[[[10, 16]], [[30, 36]], [[40, 46]]], scores=[3., 4., 5.],
                   selected=[0, 1, 2], reference_used_in_scoring=False)
        targets = dict(positive_chains=[[(10, 16)], [(20, 26)]], protected_chains=[[(0, 6)]],
                       unknown_chain_spans=[(30, 36)], callable_bases=np.ones(50, dtype=bool))
        original = dict(window=w, labels={'A31_A3': [-1]},
            counts={'A31_A3': D.R.selected_counts(prior['variants']['A31_A3'], [-1])},
            eligible_references=2, eligible_five_or_more=0)
        actual = D.label_record(row, prior, original, targets)
        self.assertEqual(actual['labels'], {'A31_A3': [-1], 'A43': [1, -1, 0]})
        self.assertEqual(actual['counts']['A43']['selected_unknown'], 1)
        self.assertEqual(actual['eligible_references'], 2)
        with self.assertRaises(ValueError): D.label_record(row, prior, dict(original, eligible_references=1), targets)

    def test_original_strict_per_species_five_conditions_under_A43_identity(self):
        counts = {}
        for s, correct, nonmatch, long in zip(sorted(D.R.C.R.ALLOW), (1531, 478), (347, 251), (548, 152)):
            old = dict(selected_correct=correct, selected_known_nonmatch=nonmatch, selected_five_or_more_correct=long)
            counts[s] = dict(A31_A3=old, A43=dict(old))
        self.assertEqual(sum(D.gate(counts).values()), 2)
        for values in counts.values():
            values['A43']['selected_correct'] += 1; values['A43']['selected_five_or_more_correct'] += 1
        self.assertTrue(all(D.gate(counts).values()))
        counts[sorted(counts)[0]]['A43']['selected_known_nonmatch'] += 1
        self.assertEqual(sum(D.gate(counts).values()), 4)

    def test_actual_receipt_union_denominator_order_weight_and_empty_known_draw(self):
        draws = [dict(window_id='a', species='toy', importance_weight_uniform_within_stratum=2.),
                 dict(window_id='b', species='toy', importance_weight_uniform_within_stratum=.5)]
        records = {'a': dict(counts=dict(label_1=1, label_0=3, **{'label_-1': 2})),
                   'b': dict(counts=dict(label_1=0, label_0=0, **{'label_-1': 1}))}
        rows = []
        for i in range(6):
            epoch, ordinal = divmod(i, 2); d = draws[ordinal]
            rows.append(dict(step=i+1, epoch=epoch+1, draw_ordinal=ordinal, window_id=d['window_id'], species='toy',
                importance_weight=d['importance_weight_uniform_within_stratum'], weighted_loss=2. if ordinal == 0 else 0.,
                log_partition=7., gold_gain=3., gradient_norm_before_clip=1., **records[d['window_id']]['counts']))
        self.assertEqual(verify_training_rows(rows, draws, records), dict(label_1=1, label_0=3, **{'label_-1': 3}))
        for key, value in (('weighted_loss', 4.), ('importance_weight', 1.), ('label_0', 2), ('draw_ordinal', 1)):
            changed = copy.deepcopy(rows); changed[0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): verify_training_rows(changed, draws, records)


if __name__ == '__main__':
    unittest.main()
