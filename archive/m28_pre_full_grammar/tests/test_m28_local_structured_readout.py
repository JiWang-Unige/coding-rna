"""A44 paired adapter and receipt checks; no research data or model weights."""
import copy
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
import torch

import readout_local_structured_emission_train as D
from verify_local_structured_emission_train import verify_training_rows


def fit_header(arm):
    return dict(arm=arm, configuration=copy.deepcopy(D.F.CONFIG), protocol=D.F.PROTOCOL,
        optimizer_steps=4608, feature_revision=D.R.C.I.F.REVISION, completed=True,
        source_support_job='13424089', source_preflight_job='13425772', position_checks_job='13425929',
        final_checkpoint=D.CONFIG['source_checkpoints'][arm], trainable_parameters=366191,
        initialization='fresh_seed0_no_supervised_weight_inheritance',
        same_capacity_fresh_seed0_paired_initialization_verified=True,
        loss_coefficients=D.F.CONFIG['loss_coefficients_by_arm'][arm],
        structured_candidate_supervision_used=arm == 'A44_JOINT',
        original_TRAIN_catalog_metadata_read=True, original_R4_local_targets_regenerated=True,
        **{k: False for k in ('supervised_checkpoint_read', 'old_cached_label_fields_used',
            'R4_policy_or_catalog_changed', 'generator_updated', 'new_backbone_forward',
            'DEV_used_for_selection', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened')})


class PairedReadoutTests(unittest.TestCase):
    def test_paired_source_identity_last_step_and_actual_R4_header(self):
        for arm in D.ARMS:
            value = fit_header(arm); self.assertIsNone(D.validate_fit(value, arm))
            for key, wrong in (('arm', D.ARMS[1-D.ARMS.index(arm)]), ('optimizer_steps', 1536),
                    ('loss_coefficients', [1., 9.]), ('original_TRAIN_catalog_metadata_read', False),
                    ('source_preflight_job', '13424115'), ('DEV_used_for_selection', True)):
                with self.subTest(arm=arm, key=key), self.assertRaises(ValueError):
                    D.validate_fit(dict(value, **{key: wrong}), arm)

    def test_64_task_mapping_exactly_covers_each_original_population(self):
        coverage = {a: [] for a in D.ARMS}
        for task in range(64):
            arm, shard = D.task_identity(task)
            coverage[arm].extend(range(shard, 1491, 32))
        for arm in D.ARMS: self.assertEqual(sorted(coverage[arm]), list(range(1491)))
        self.assertEqual(D.task_identity(31), ('A44_LOCAL', 31))
        self.assertEqual(D.task_identity(32), ('A44_JOINT', 0))
        for task in (-1, 64):
            with self.assertRaises(ValueError): D.task_identity(task)

    def test_paired_R4_join_keeps_unknown_and_full_reference_denominators(self):
        w = dict(id='toy', species=sorted(D.R.C.R.ALLOW)[0], valid_bases=50)
        baseline = dict(chains=[[[0, 6]]], scores=[2.], selected=[0])
        old = dict(window=w, variants={'A31_A3': baseline})
        rows = {'A44_LOCAL': dict(arm='A44_LOCAL', window=w, chains=[[[40, 46]]], scores=[5.], selected=[0]),
            'A44_JOINT': dict(arm='A44_JOINT', window=w, chains=[[[10, 16]], [[30, 36]], [[40, 46]]],
                             scores=[3., 4., 5.], selected=[0, 1, 2])}
        targets = dict(positive_chains=[[(10, 16)], [(20, 26)]], protected_chains=[[(0, 6)]],
                       unknown_chain_spans=[(30, 36)], callable_bases=np.ones(50, dtype=bool))
        original = dict(window=w, labels={'A31_A3': [-1]},
            counts={'A31_A3': D.R.selected_counts(baseline, [-1])}, eligible_references=2, eligible_five_or_more=0)
        actual = D.label_paired_record(rows, old, original, targets)
        self.assertEqual(actual['labels'], {'A31_A3': [-1], 'A44_LOCAL': [0], 'A44_JOINT': [1, -1, 0]})
        self.assertEqual(actual['eligible_references'], 2)
        self.assertEqual(actual['counts']['A44_JOINT']['selected_unknown'], 1)
        swapped = dict(A44_LOCAL=rows['A44_JOINT'], A44_JOINT=rows['A44_LOCAL'])
        with self.assertRaises(ValueError): D.label_paired_record(swapped, old, original, targets)
        with self.assertRaises(ValueError): D.label_paired_record(rows, old, dict(original, eligible_references=1), targets)

    def test_joint_must_pass_against_both_stronger_baseline_and_LOCAL(self):
        counts = {}
        for s in D.R.C.R.ALLOW:
            old = dict(selected_correct=10, selected_known_nonmatch=3, selected_five_or_more_correct=2)
            counts[s] = dict(A31_A3=old, A44_LOCAL=dict(old, selected_correct=12, selected_five_or_more_correct=3),
                            A44_JOINT=dict(old, selected_correct=11, selected_five_or_more_correct=3))
        checks, supported = D.gates(counts)
        self.assertTrue(all(checks['joint_vs_A31_A3'].values())); self.assertFalse(supported)
        for v in counts.values(): v['A44_JOINT'] = dict(v['A44_LOCAL'])
        self.assertFalse(D.gates(counts)[1])
        for v in counts.values(): v['A44_JOINT'].update(selected_correct=13, selected_five_or_more_correct=4)
        self.assertTrue(D.gates(counts)[1])
        counts[next(iter(counts))]['A44_JOINT']['selected_known_nonmatch'] = 4
        self.assertFalse(D.gates(counts)[1])

    def test_actual_paired_receipt_units_masks_importance_and_zero_known_draw(self):
        draws = [dict(window_id='a', species='toy', importance_weight_uniform_within_stratum=2.),
                 dict(window_id='b', species='toy', importance_weight_uniform_within_stratum=.5)]
        records = {'a': dict(counts=dict(label_1=1, label_0=3, **{'label_-1': 2})),
                   'b': dict(counts=dict(label_1=0, label_0=0, **{'label_-1': 1}))}
        windows = {'a': dict(valid_bases=18), 'b': dict(valid_bases=24)}; positions = {'a': 3, 'b': 0}
        for arm in D.ARMS:
            rows = []
            for i in range(6):
                epoch, ordinal = divmod(i, 2); draw = draws[ordinal]; ident = draw['window_id']
                numerator = 4. if arm == 'A44_JOINT' and ident == 'a' else 0.
                local = 1.5 if ident == 'a' else 0.; structured = numerator/windows[ident]['valid_bases']
                rows.append(dict(step=i+1, epoch=epoch+1, draw_ordinal=ordinal, arm=arm, window_id=ident,
                    species='toy', valid_bases=windows[ident]['valid_bases'], importance_weight=draw['importance_weight_uniform_within_stratum'],
                    loss_coefficients=D.F.CONFIG['loss_coefficients_by_arm'][arm], local_loss=local,
                    structured_loss=structured, structured_candidate_mean_loss=numerator/4 if numerator else 0.,
                    weighted_loss=(local+structured)*draw['importance_weight_uniform_within_stratum'],
                    log_partition=numerator, gold_gain=0., local_supervised_positions=positions[ident],
                    structured_known_candidates_used=4 if numerator else 0,
                    structured_unknown_chains_excluded=records[ident]['counts']['label_-1'] if arm == 'A44_JOINT' else 0,
                    gradient_norm_before_clip=1., **records[ident]['counts']))
            self.assertEqual(verify_training_rows(arm, rows, draws, records, windows, positions),
                             dict(label_1=1, label_0=3, **{'label_-1': 3}))
            for key, wrong in (('weighted_loss', 9.), ('valid_bases', 24), ('local_supervised_positions', 18),
                    ('importance_weight', 1.), ('loss_coefficients', [2., 1.]), ('draw_ordinal', 1)):
                changed = copy.deepcopy(rows); changed[0][key] = wrong
                with self.subTest(arm=arm, key=key), self.assertRaises(ValueError):
                    verify_training_rows(arm, changed, draws, records, windows, positions)

    def test_artificial_potential_cannot_cross_arm_or_last_model_identity(self):
        s = next(iter(D.R.C.R.ALLOW)); seqid = sorted(D.R.C.R.ALLOW[s])[0]
        w = dict(id='toy', species=s, split='train', seqid=seqid, strand='+', start=0, end=6,
                 width=24576, valid_bases=6, pad_right_after_orientation=24570)
        with tempfile.TemporaryDirectory() as folder, patch.object(D, 'ROOT', Path(folder)), \
                patch.object(D, 'POTENTIALS', Path(folder)/'potentials'):
            arm = 'A44_LOCAL'; path = D.potential_path(arm, 0); path.parent.mkdir(parents=True)
            value = dict(configuration=D.CONFIG, protocol=D.PROTOCOL, arm=arm, fit_job=D.FIT_JOBS[arm],
                original_ordinal=0, window=w, sequence='ATGTAA', local_logits=torch.zeros(6, 15),
                revision=D.R.C.I.F.REVISION, reference_used_in_scoring=False)
            receipt = dict(artifact=str(path.relative_to(D.ROOT)), original_ordinal=0, window=w)
            report = dict(feature_revision=D.R.C.I.F.REVISION)
            torch.save(value, path); self.assertEqual(D.load_potential(arm, 0, receipt, report)['arm'], arm)
            for key, wrong in (('arm', 'A44_JOINT'), ('fit_job', D.FIT_JOBS['A44_JOINT']), ('original_ordinal', 1)):
                torch.save(dict(value, **{key: wrong}), path)
                with self.subTest(key=key), self.assertRaises(ValueError): D.load_potential(arm, 0, receipt, report)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
