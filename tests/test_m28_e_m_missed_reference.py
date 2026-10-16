import copy
import unittest
import torch

import diagnose_full_grammar_e_m_missed_reference as F


def row(i, hit, legal=True):
    return dict(species='fixture', seqid='fixture', strand='+', reference_CDS=[(i*12, i*12+9)],
        CDS_segments=1, span_bp=9, exact_primary_hit=hit, single_chain_supported=legal,
        rows_by_role=[0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0] if legal else None,
        category='exact_primary_chain' if hit else ('miss' if legal else 'outside_fixed_single_chain_grammar'))


class MissedReferenceTests(unittest.TestCase):
    def test_paired_hit_groups_and_category_migrations(self):
        e = [row(i, i in (0,1)) for i in range(4)]
        m = [row(i, i in (0,2)) for i in range(4)]
        result = F.paired_rows(e, m)
        self.assertEqual(result['hit_groups'], dict(both_hit=1, E_only=1, M_only=1, both_missed=1))
        self.assertEqual(result['M_minus_E_exact_hits'], 0)
        self.assertEqual(sum(result['category_transitions'].values()), 4)

    def test_unsupported_and_long_span_retained_and_mismatched_geometry_rejected(self):
        rows = [row(0, False, False), row(1, True)]; rows[0]['span_bp'] = 30000
        measured = F.counts(rows)
        self.assertEqual(measured['reference'], 2)
        self.assertEqual(measured['individually_supported_reference'], 1)
        self.assertEqual(measured['strata']['span_longer_than_native_window_not_excluded']['reference'], 1)
        for wrong in (rows[:1], rows+[rows[0]], copy.deepcopy(rows)):
            if len(wrong) == 2:
                wrong[1]['rows_by_role'][1] += 1
            with self.assertRaises(ValueError):
                F.paired_rows(rows, wrong)

    def test_touching_nonoverlap_and_multiple_conflict_role_margins(self):
        gains, counts = [0.]*15, [0]*15; gains[1]=4.; counts[1]=4
        a = dict(parts=[(0,9)], id='a', gain=3., gain_by_role=[0.,3.]+[0.]*13, rows_by_role=counts)
        b = dict(parts=[(12,21)], id='b', gain=2., gain_by_role=[0.,2.]+[0.]*13, rows_by_role=counts)
        index = ([a,b], [9,21])
        self.assertEqual(F.D.conflicting([(9,12)], index), [])
        conflicts = F.D.conflicting([(3,18)], index)
        mass, terms, rows, margin, delta = F.R.margin_terms(4., gains, conflicts)
        self.assertEqual((mass, margin, delta[1]), (5., -1., -1.))
        self.assertEqual(sum(rows), 8)
        self.assertEqual(F.D.classify(False, True, 4., mass), 'supported_positive_reference_outscored_by_selected_conflicts')
        self.assertEqual(F.D.classify(False, True, 6., mass), 'positive_reference_improves_saved_selected_set')

    def test_actual_carry_literal_support_and_score_role_partition(self):
        sequence='ATGCGTNAGTTTAA'; parts=[(0,4),(9,14)]
        path = F.D.role_path(parts)
        self.assertTrue(F.D.supported(sequence, parts, path, 2))
        self.assertFalse(F.D.supported(sequence, parts, path, 1))
        self.assertFalse(F.D.supported('ATGCCCNAGTTTAA', parts, path, 2))
        values = torch.arange(len(sequence)*15, dtype=torch.float32).reshape(-1,15).numpy()
        gain, terms, rows = F.D.chain_score(values, parts, path)
        result = F.T.terms(dict(gain_by_role=terms, rows_by_role=rows, reference_gain_vs_background=gain))
        self.assertAlmostEqual(sum(g['gain'] for g in result.values()), gain)
        self.assertEqual(sum(g['rows'] for g in result.values()), len(sequence))
