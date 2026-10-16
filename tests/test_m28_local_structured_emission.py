"""Artificial loss/gradient checks, not a model or original-data experiment."""
import math
import unittest
import torch

from src.m28.labels import IGNORE, CDS_COLUMNS, build_targets
from src.m28.local_structured_emission import local_structured_emission_loss


def unknown_local(length):
    return dict(seg=torch.full((length,), IGNORE, dtype=torch.long),
                region=torch.full((length,), IGNORE, dtype=torch.long))


def three_local_positions(length):
    targets = unknown_local(length)
    targets['seg'][0] = 0
    targets['seg'][3] = 13
    targets['region'][6] = 1
    return targets


def hand_single_exon_gain(logits, start, end):
    """Independent direct roles for the test's single-exon complete chains."""
    positions = torch.arange(start, end)
    roles = 1+(positions-start) % 3
    roles[:3] = 13
    roles[-3:] = 14
    values = logits.double()
    return (values[positions, roles]-values[positions, 0]).sum()


class LocalStructuredEmissionTests(unittest.TestCase):
    def test_separate_denominators_and_grouped_region_formula(self):
        logits = torch.zeros((24, 15), dtype=torch.float64, requires_grad=True)
        result = local_structured_emission_loss(logits,
            [[(3, 12)], [(3, 15)], [(12, 21)]], [1, 0, -1],
            three_local_positions(24), local_coefficient=2, structured_coefficient=3)
        local = math.log(15)-math.log(len(CDS_COLUMNS))/3
        structured = math.log(3)/2  # empty, gold, or overlapping nonmatch
        self.assertAlmostEqual(float(result['local_loss']), local, places=12)
        self.assertAlmostEqual(float(result['structured_loss']), structured, places=12)
        self.assertAlmostEqual(float(result['loss']), 2*local+3*structured, places=12)
        self.assertEqual(result['local_supervised_positions'], 3)
        self.assertEqual(result['known_candidates'], 2)
        self.assertEqual(result['unknown_excluded'], 1)

    def test_shared_FP32_gradient_matches_independent_formula(self):
        values = torch.sin(torch.arange(450, dtype=torch.float32)).reshape(30, 15)*.2
        logits = values.clone().requires_grad_(True)
        result = local_structured_emission_loss(logits,
            [[(3, 12)], [(3, 15)], [(18, 27)]], [1, 0, -1],
            three_local_positions(30), local_coefficient=2, structured_coefficient=3)
        result['loss'].backward()
        reference = values.clone().requires_grad_(True)
        logp = reference.log_softmax(-1)
        local = (-logp[0, 0]-logp[3, 13]
                 -torch.logsumexp(logp[6, list(CDS_COLUMNS)], dim=0))/3
        gold = hand_single_exon_gain(reference, 3, 12)
        wrong = hand_single_exon_gain(reference, 3, 15)
        structured = (torch.logsumexp(torch.stack((gold.new_zeros(()), gold, wrong)), 0)-gold)/2
        (2*local+3*structured).backward()
        self.assertEqual(result['loss'].dtype, torch.float64)
        self.assertEqual(logits.grad.dtype, torch.float32)
        self.assertTrue(bool(torch.isfinite(logits.grad).all()))
        torch.testing.assert_close(logits.grad, reference.grad, rtol=1e-6, atol=1e-7)
        self.assertEqual(int(torch.count_nonzero(logits.grad[18:27])), 0)

    def test_unknown_chain_does_not_erase_resolved_local_region(self):
        logits = torch.zeros((18, 15), dtype=torch.float64, requires_grad=True)
        targets = unknown_local(18)
        targets['region'][5] = 1
        result = local_structured_emission_loss(logits, [[(3, 12)]], [-1], targets,
            local_coefficient=1, structured_coefficient=1)
        self.assertEqual(result['known_candidates'], 0)
        self.assertEqual(float(result['structured_loss']), 0)
        self.assertEqual(result['local_supervised_positions'], 1)
        self.assertAlmostEqual(float(result['loss']), math.log(15/len(CDS_COLUMNS)), places=12)
        result['loss'].backward()
        self.assertTrue(bool((logits.grad[5, list(CDS_COLUMNS)] < 0).all()))
        self.assertEqual(int(torch.count_nonzero(logits.grad[:5])), 0)
        self.assertEqual(int(torch.count_nonzero(logits.grad[6:])), 0)

    def test_unknown_local_positions_do_not_erase_known_chain(self):
        logits = torch.zeros((24, 15), dtype=torch.float64, requires_grad=True)
        result = local_structured_emission_loss(logits, [[(3, 12)]], [1], unknown_local(24),
            valid_bases=18, local_coefficient=1, structured_coefficient=1)
        self.assertEqual(result['local_supervised_positions'], 0)
        self.assertEqual(float(result['local_loss']), 0)
        self.assertAlmostEqual(float(result['loss']), math.log(2), places=12)
        result['loss'].backward()
        self.assertGreater(int(torch.count_nonzero(logits.grad[3:12])), 0)
        self.assertEqual(int(torch.count_nonzero(logits.grad[18:])), 0)

    def test_empty_or_all_unknown_supervision_is_differentiable_zero(self):
        for chains, labels in (([], []), ([[(3, 12)]], [-1])):
            with self.subTest(chains=chains):
                logits = torch.zeros((18, 15), dtype=torch.float32, requires_grad=True)
                result = local_structured_emission_loss(logits, chains, labels, unknown_local(18),
                    local_coefficient=1, structured_coefficient=1)
                self.assertEqual(float(result['loss']), 0)
                result['loss'].backward()
                self.assertEqual(int(torch.count_nonzero(logits.grad)), 0)

    def test_original_target_builder_preserves_fine_unknown_N_and_padding(self):
        window = dict(start=0, end=18, strand='+', width=24, valid_bases=18,
            overlapping_transcript_ids=['toy'], positive_ids=['toy'],
            fine_label_unknown_intervals=[(6, 9)], chain_unknown_intervals=[],
            no_annotated_gene_either_strand=False)
        transcript = dict(id='toy', gene_id='toy_gene', training_primary=True,
            CDS=[(3, 12, 0)], CDS_complete_supported=True)
        sequence = 'AAAATGAAATAAAAANAA'
        self.assertEqual(len(sequence), 18)
        targets = build_targets(window, {'toy': transcript}, sequence)
        self.assertTrue(bool((targets['seg'][6:9] == IGNORE).all()))
        self.assertTrue(bool((targets['region'][6:9] == 1).all()))
        logits = torch.zeros((24, 15), dtype=torch.float32, requires_grad=True)
        result = local_structured_emission_loss(logits, [[(3, 12)]], [1], targets,
            valid_bases=18, local_coefficient=1, structured_coefficient=1)
        self.assertEqual(result['local_supervised_positions'], 17)
        result['loss'].backward()
        self.assertEqual(int(torch.count_nonzero(logits.grad[15])), 0)
        self.assertEqual(int(torch.count_nonzero(logits.grad[18:])), 0)
        self.assertGreater(int(torch.count_nonzero(logits.grad[16])), 0)

    def test_wrong_local_window_or_known_padding_is_rejected(self):
        logits = torch.zeros((24, 15), dtype=torch.float32, requires_grad=True)
        with self.assertRaisesRegex(ValueError, 'full oriented'):
            local_structured_emission_loss(logits, [], [], unknown_local(18),
                valid_bases=18, local_coefficient=1, structured_coefficient=1)
        for field in ('seg', 'region'):
            with self.subTest(field=field):
                targets = unknown_local(24)
                targets[field][18] = 0
                with self.assertRaisesRegex(ValueError, 'ignore padded'):
                    local_structured_emission_loss(logits, [], [], targets,
                        valid_bases=18, local_coefficient=1, structured_coefficient=1)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
