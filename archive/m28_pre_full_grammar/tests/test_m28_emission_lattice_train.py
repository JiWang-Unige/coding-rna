"""New full TRAIN adapters: artificial geometry only, no research weights or DNA."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_emission_lattice_train as D
import verify_emission_lattice_train as V


def fixture():
    species = sorted(D.R.C.R.ALLOW)[0]; w = dict(id='toy', species=species, valid_bases=50)
    # Old geometry is protected; the newly generated, distinct geometry is primary.
    old = dict(chains=[[[0, 6]]], scores=[2.], selected=[0])
    new = dict(window=w, chains=[[[10, 16]], [[30, 36]], [[40, 46]]], scores=[3., 4., 5.], selected=[0, 1, 2],
               original_ordinal=0, shard=0, reference_used_in_scoring=False)
    targets = dict(positive_chains=[[(10, 16)], [(20, 26)]], protected_chains=[[(0, 6)]],
                   unknown_chain_spans=[(30, 36)], callable_bases=np.ones(50, dtype=bool))
    labels = dict(window=w, labels={'A31_A3': [-1]}, counts={'A31_A3': D.R.selected_counts(old, [-1])},
                  eligible_references=2, eligible_five_or_more=0)
    return new, dict(window=w, variants={'A31_A3': old}), labels, targets


class EmissionLatticeTrain(unittest.TestCase):
    def test_R4_label_metadata_is_not_natural_geometry_and_all_geometry_fields_match(self):
        natural = dict(id='toy', species=sorted(D.R.C.R.ALLOW)[0], split='train',
            seqid='toy_chromosome', strand='-', start=12, end=62, width=24576,
            valid_bases=50, pad_right_after_orientation=24526)
        target = dict(natural, label_policy='CDS_assessable_primary_R4',
            positive_ids=['primary'], overlapping_transcript_ids=['primary', 'alternative'],
            primary_complete_ids=['primary'], motifs={}, no_annotated_gene_either_strand=False,
            fine_label_unknown_intervals=[[30, 36]], chain_unknown_intervals=[[30, 36]])
        original = deepcopy(target)
        self.assertNotEqual(target, natural)  # The original whole-dictionary guard fails.
        self.assertIsNone(D.validate_target_geometry(target, natural))
        self.assertEqual(target, original)  # Full metadata remains available to build_targets.
        self.assertEqual(set(natural), set(D.R.C.R.GEOMETRY))
        for key in D.R.C.R.GEOMETRY:
            wrong = deepcopy(target)
            wrong[key] = wrong[key]+1 if isinstance(wrong[key], int) else wrong[key]+'_different'
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'geometry differ'):
                D.validate_target_geometry(wrong, natural)

    def test_new_geometry_labels_not_old_indices_unknown_kept_and_full_denominator(self):
        new, old, labels, targets = fixture()
        row = D.label_record(new, old, labels, targets)
        self.assertEqual(row['labels'], {'A31_A3': [-1], 'A42': [1, -1, 0]})
        self.assertEqual(row['counts']['A42'], V.counts(new, [1, -1, 0]))
        self.assertEqual(row['counts']['A42']['selected_unknown'], 1)
        self.assertEqual(row['eligible_references'], 2)  # Not reduced to the selected primary.
        self.assertNotIn('candidates', row['counts']['A42'])

    def test_old_policy_drift_and_stale_selection_rejected(self):
        new, old, labels, targets = fixture()
        wrong = deepcopy(labels); wrong['labels']['A31_A3'] = [0]
        with self.assertRaisesRegex(ValueError, 'frozen R4'): D.label_record(new, old, wrong, targets)
        wrong = deepcopy(labels); wrong['eligible_references'] = 1
        with self.assertRaisesRegex(ValueError, 'denominators'): D.label_record(new, old, wrong, targets)
        wrong = deepcopy(new); wrong['selected'] = [0, 2]
        with self.assertRaisesRegex(ValueError, 'reconcile'): D.label_record(wrong, old, labels, targets)

    def test_complete_shards_restore_original_order(self):
        new, _, _, _ = fixture(); other = deepcopy(new)
        other.update(original_ordinal=1, shard=1, window=dict(new['window'], id='toy2'))
        free = [dict(window=new['window']), dict(window=other['window'])]
        self.assertEqual(D.assemble_records([other, new], free), [new, other])
        for rows in ([new], [new, new], [new, dict(other, shard=0)],
                     [new, dict(other, reference_used_in_scoring=True)]):
            with self.subTest(rows=rows), self.assertRaises(ValueError): D.assemble_records(rows, free)

    def test_null_and_fp64_role_arithmetic_without_float32_sum_rounding(self):
        local = torch.zeros(6, 15); local[0, 13] = 2**24; local[1, 13] = 1.
        actual, error = V.verify_scores(local, [[[0, 6]]], [float(2**24+1)])
        self.assertEqual((actual, error), ([float(2**24+1)], 0.))
        with self.assertRaises(AssertionError): V.verify_scores(local, [[[0, 6]]], [float(2**24)])
        self.assertEqual(V.verify_scores(local, [], []), ([], 0.))
        self.assertEqual(V.counts(dict(chains=[], scores=[], selected=[]), [])['selected'], 0)
        sequence = 'ATGTAA'; value = dict(sequence=sequence, local_logits=torch.zeros(6, 15),
                                         original_chains=[], original_scores=[], original_selected=[])
        self.assertEqual(D.validate_prediction(value, dict(chains=[], chain_gains=[], gain=0., reference_used=False)), 0.)

    def test_original_five_conditions_strict_and_per_species(self):
        metrics = {}
        for s, correct, nonmatch, long in zip(sorted(D.R.C.R.ALLOW), (1531, 478), (347, 251), (548, 152)):
            old = dict(selected_correct=correct, selected_known_nonmatch=nonmatch, selected_five_or_more_correct=long)
            metrics[s] = {'A31_A3': old, 'A42': dict(old)}
        self.assertEqual((len(D.gate(metrics)), sum(D.gate(metrics).values())), (5, 2))
        for value in metrics.values():
            value['A42']['selected_correct'] += 1
            value['A42']['selected_five_or_more_correct'] += 1
        self.assertTrue(all(D.gate(metrics).values()))
        metrics[sorted(metrics)[0]]['A42']['selected_known_nonmatch'] += 1
        self.assertEqual(sum(D.gate(metrics).values()), 4)


if __name__ == '__main__':
    unittest.main()
