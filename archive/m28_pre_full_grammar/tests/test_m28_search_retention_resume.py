"""Continuation driver: synthetic CPU objectives/ledger; no source fit or real data."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
import torch
from src.m28.query_collaboration import CollaborativeGeneQueryModel
from src.m28.query_known_mean import joint_known_query_losses
from src.m28.training import one_hot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import resume_query_search_retention as S


class SearchRetentionResumeTests(unittest.TestCase):
    def test_live_source_cannot_load_weights_TRAIN_or_create_recovery(self):
        with patch.object(S.K, 'source_receipt', return_value=dict(job=S.K.SOURCE_JOB, state='RUNNING', exit_code='0:0')), \
                patch.object(S.K, 'restore_checkpoint') as restore, patch.object(S.R.Q, 'load_train_inputs') as inputs, \
                patch.object(S, 'bind_outputs') as bind:
            with self.assertRaisesRegex(ValueError, 'TIMEOUT'):
                S.fit()
            restore.assert_not_called(); inputs.assert_not_called(); bind.assert_not_called()

    def test_copied_and_continued_ledger_has_exact_original_epochs_no_duplicate(self):
        species = list(S.R.A.ALLOW); saved = S.accumulator()
        for i in range(4608):
            epoch, ordinal = divmod(i, 1536)
            row = dict(step=i+1, epoch=epoch+1, draw_ordinal=ordinal, species=species[ordinal % 2],
                losses=dict(local=float(epoch+1), existence_bce=1., boundary_ce=2., termination_bce=3.),
                search_retention_margin=.5, weighted_total=float(epoch+1)+6.5,
                query_label_counts={'-1': 0, '0': 15, '1': 1}, protection_reasons=[],
                load_forward_supervision_search_backward_optimizer_seconds=1.)
            S.accumulate(saved, row)
            if i == 3071:
                with self.assertRaisesRegex(ValueError, 'three original'):
                    S.epochs(saved)
                copied = copy.deepcopy(saved)
        results = S.epochs(saved)
        self.assertEqual(copied['steps'], [1536, 1536, 0])
        self.assertEqual(saved['total']-copied['total'], 1536)
        self.assertEqual([r['mean_losses']['local'] for r in results], [1., 2., 3.])
        self.assertEqual([r['steps'] for r in results], [1536]*3)
        self.assertEqual([r['seconds'] for r in results], [1536.]*3)
        for result in results:
            for s in species:
                self.assertEqual(result['supervision_by_species'][s], dict(query_0=11520, query_1=768, **{'query_-1': 0}))
        with self.assertRaisesRegex(ValueError, 'ordered'):
            S.accumulate(copied, row)  # Skipping 3073 directly to 4608 is not continuation.

    def test_update_preserves_four_loss_plus_weighted_extra_and_clip(self):
        torch.manual_seed(0); model = CollaborativeGeneQueryModel(4, 2, 8)
        oracle = copy.deepcopy(model)
        features = torch.randn(1, 8, 4); sequence = 'ATGAAATAA'*5
        dna = one_hot(sequence, len(sequence))[None]
        w = dict(id='synthetic', species='synthetic', valid_bases=len(sequence))
        seg = np.zeros(len(sequence), dtype=np.int64); seg[:3] = 13; seg[3:6] = [1, 2, 3]; seg[6:9] = 14
        region = np.zeros(len(sequence), dtype=np.int64); region[:9] = 1
        targets = dict(seg=seg, region=region,
            positive_chains=[[(0, 9)]], protected_chains=[[(0, 9)]], unknown_chain_spans=[],
            callable_bases=np.ones(len(sequence), dtype=bool))
        draw = dict(importance_weight_uniform_within_stratum=2.5)
        def base(m, outputs, window, catalog, seq, weight):
            return joint_known_query_losses(m, outputs, targets, len(sequence), weight, S.R.CONFIG['loss_coefficients'])
        def synthetic_extra(window, seq, assigned, target, margin):
            return dict(loss=window.initial.square().mean()+1., positive_queries=len(assigned),
                        records=[{'synthetic_math_only': True} for _ in assigned], reference_injected_candidates=0)
        optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
        expected_optimizer = torch.optim.AdamW(oracle.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
        outputs = oracle(features, dna, [len(sequence)])
        expected_base = base(oracle, outputs, w, {}, sequence, 2.5)
        expected_extra = synthetic_extra(oracle.query_window(outputs['features'][0], len(sequence)),
                                         sequence, expected_base['query']['positive_chains'], targets, 1.)
        expected_total = expected_base['total']+2.5*expected_extra['loss']
        expected_total.backward(); expected_norm = torch.nn.utils.clip_grad_norm_(oracle.parameters(), 1.)
        expected_optimizer.step()
        with patch.object(S.R.Q, 'window_inputs', return_value=(features, dna)), \
                patch.object(S.R.A, 'supervised_readout', side_effect=base), \
                patch.object(S.R.Q, 'build_targets', return_value=targets), \
                patch.object(S.R, 'search_retention_loss', side_effect=synthetic_extra), \
                patch.object(torch.cuda, 'synchronize'):
            row = S.current_update(model, optimizer, draw, w, dict(sequence=sequence), {})
        self.assertEqual(set(row['losses']), {'local', 'existence_bce', 'boundary_ce', 'termination_bce'})
        self.assertAlmostEqual(row['weighted_total'], float(expected_total.detach()), places=6)
        self.assertAlmostEqual(row['gradient_norm_before_clip'], float(expected_norm), places=6)
        self.assertTrue(all(torch.equal(v, oracle.state_dict()[k]) for k, v in model.state_dict().items()))
        self.assertEqual(row['gene_to_query'], expected_base['query']['gene_to_query'])
        self.assertGreater(row['search_retention_margin'], 0.)
        self.assertEqual(row['reference_injected_candidates'], 0)


if __name__ == '__main__':
    unittest.main()
