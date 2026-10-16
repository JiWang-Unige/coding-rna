#!/usr/bin/env python3
"""Artificial-only four-cell identity, denominator and interaction arithmetic checks."""
import copy
import json
import os
import time
import unittest
import aggregate_full_grammar_context_prior_train as A
import check_full_grammar_continuation_train as C


def fixtures():
    original = C.fixtures()
    mapping = {'F00': 'M28_STRUCTURED_CONTINUE_final9216', 'F01': 'M28_NATIVE_CDS_PRIOR_final9216',
        'F10': 'M28_JOINT24576_final9216', 'F11': 'M28_JOINT24576_final9216'}
    reports = {cell: copy.deepcopy(original[label]) for cell, label in mapping.items()}
    for cell, increment in zip(A.SPECS, (0, 1, 2, 4)):
        report = reports[cell]
        report.update(job_id=A.SPECS[cell][2], checkpoint=str(A.SPECS[cell][1]))
        for seqid, row in report['chromosomes'].items():
            count = A.REFERENCE_COUNTS[seqid]
            row['exact_chain'] = A.E.C.S.prf(count//2+increment, count, count)
        C.totals(report)
    return reports


class ComparisonChecks(unittest.TestCase):
    def test_four_cells_micro_and_all_six_percentage_point_contrasts(self):
        reports = fixtures(); saved = copy.deepcopy(reports)
        result = A.build_comparison(reports)
        self.assertEqual(reports, saved)
        self.assertEqual(result['cells'], ['F00', 'F01', 'F10', 'F11'])
        self.assertEqual(result['pooled_reference_primary_chains'], 22688)
        self.assertEqual(result['reference_primary_chains'], {'arabidopsis_thaliana': 17739, 'oryza_sativa': 4949})
        for name, net_tp in (('context_without_prior', 2), ('prior_without_joint_context', 1),
                ('context_with_prior', 3), ('prior_with_joint_context', 2),
                ('combined_minus_control', 4), ('descriptive_interaction', 1)):
            row = result['exact_chain_F1_contrasts_percentage_points'][name]
            self.assertAlmostEqual(row['pooled_micro'], 600*net_tp/22688, delta=1e-12)
            for species, count in result['reference_primary_chains'].items():
                self.assertAlmostEqual(row['per_species'][species], 300*net_tp/count, delta=1e-12)
            self.assertAlmostEqual(row['macro_species'], 150*net_tp*(1/17739+1/4949), delta=1e-12)
        self.assertEqual(result['measurements']['F11']['chromosomes'], reports['F11']['chromosomes'])
        self.assertTrue(result['cumulative_optimizer_steps_matched'])

    def test_worst_species_switch_uses_cell_minima_not_minimum_contrast(self):
        reports = fixtures()
        for cell, multipliers in zip(A.SPECS, ((1, 3), (3, 1), (1, 2), (2, 3))):
            for species, seqid, _ in A.E.I.SCOPE:
                row = reports[cell]['chromosomes'][seqid]
                count = A.REFERENCE_COUNTS[seqid]
                numerator = multipliers[A.SPECIES.index(species)]
                row['exact_chain'] = A.E.C.S.prf(numerator*(count//4), count, count)
            C.totals(reports[cell])
        result = A.build_comparison(reports)
        bases = {s: sum(A.REFERENCE_COUNTS[q]//4 for t, q, _ in A.E.I.SCOPE if t == s)/
            sum(A.REFERENCE_COUNTS[q] for t, q, _ in A.E.I.SCOPE if t == s) for s in A.SPECIES}
        at, rice = (bases[s] for s in A.SPECIES)
        row = result['exact_chain_F1_contrasts_percentage_points']['descriptive_interaction']
        self.assertAlmostEqual(row['worst_species'], 100*(2*at-rice), delta=1e-12)
        self.assertAlmostEqual(row['per_species'][A.SPECIES[0]], -100*at, delta=1e-12)
        self.assertAlmostEqual(row['per_species'][A.SPECIES[1]], 300*rice, delta=1e-12)
        self.assertNotEqual(row['worst_species'], min(row['per_species'].values()))

    def test_missing_cell_or_chromosome_and_swapped_final_identity_rejected(self):
        reports = fixtures(); del reports['F11']
        with self.assertRaisesRegex(ValueError, 'All four'): A.build_comparison(reports)
        reports = fixtures(); del reports['F11']['chromosomes']['NC_089046.1']
        with self.assertRaisesRegex(ValueError, 'scope'): A.build_comparison(reports)
        for field in ('checkpoint', 'job_id'):
            reports = fixtures(); reports['F11'][field] = reports['F10'][field]
            with self.assertRaisesRegex(ValueError, 'identity'): A.build_comparison(reports)

    def test_reference_drift_rejected_even_if_reaggregated(self):
        for field in ('primary', 'CDS', 'background', 'strata', 'isoforms'):
            reports = fixtures(); report = reports['F11']; row = report['chromosomes']['NC_003070.9']
            if field == 'primary': row['reference_primary_records'] += 1
            elif field == 'CDS': row['CDS_base_unstranded_union'] = A.E.C.S.prf(40, 100, 101)
            elif field == 'background': row['background']['bp'] += 1
            elif field == 'strata': row['retained_unsupported_reference_strata'] = {}
            else: row['reference_assessable_isoform_records'] += 1
            C.totals(report)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'reference denominator'):
                A.build_comparison(reports)

    def test_count_aggregate_phase_and_data_scope_errors_rejected(self):
        for field in ('chromosome', 'species', 'pooled', 'macro', 'worst', 'phase', 'DEV', 'sealed', 'reference'):
            reports = fixtures(); report = reports['F11']
            if field == 'chromosome': report['chromosomes']['NC_003070.9']['exact_chain']['f1'] += .01
            elif field == 'species': report['per_species'][A.SPECIES[0]]['exact_chain']['tp'] += 1
            elif field == 'pooled': report['pooled']['exact_chain']['f1'] += .01
            elif field == 'macro': report['macro_species_exact_chain_F1'] += .01
            elif field == 'worst': report['worst_species_exact_chain_F1'] += .01
            elif field == 'phase': report['chromosomes']['NC_003070.9']['phase_checked_prediction_records'] -= 1
            else: report[{'DEV': 'DEV_used', 'sealed': 'sealed_data_used', 'reference': 'reference_used_in_inference'}[field]] = True
            with self.subTest(field=field), self.assertRaises(ValueError): A.build_comparison(reports)


if __name__ == '__main__':
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ComparisonChecks))
    with (A.CHECKS/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'], tests_run=result.testsRun,
            failures=len(result.failures), errors=len(result.errors), research_data_read=False, models_or_weights_read=False,
            native_forward_run=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
            seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False); handle.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)
