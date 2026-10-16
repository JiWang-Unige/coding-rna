"""Saved A42 competition joins and replacement arithmetic on artificial rows."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_emission_lattice_competition as M


def fixture():
    old = dict(chains=[[[0,6]], [[10,13],[18,21],[26,29],[34,37],[42,45]], [[60,66]], [[70,76]]],
               scores=[2.,2.,1.,1.], selected=[0,1,2,3])
    new = dict(chains=[[[0,3],[7,10]], [[10,25]], [[25,45]], [[50,56]], [[60,66]], [[70,76]]],
               scores=[3.,3.,3.,2.,1.,1.], selected=list(range(6)), gain=13., original_baseline_gain=6.,
               reference_used_in_scoring=False)
    labels = {'A31_A3':[1,1,0,-1], 'A42':[0,0,-1,1,0,-1]}
    saved, previous, labelled = [], [], []
    for ordinal, species in enumerate(M.D.R.C.R.ALLOW):
        w = dict(id=species+'|toy', species=species, split='train')
        row = deepcopy(new); row.update(window=w, original_ordinal=ordinal)
        saved.append(row)
        previous.append(dict(window=w, variants={'A31_A3':deepcopy(old)}, reference_used_in_scoring=False))
        labelled.append(dict(window=w, labels=deepcopy(labels),
            counts={a:M.V.counts(v, labels[a]) for a,v in (('A31_A3',old),('A42',new))},
            eligible_references=4, eligible_five_or_more=3))
    return saved, previous, labelled


def closed():
    counts = {}
    for s, correct, nonmatch, long in zip(sorted(M.D.R.C.R.ALLOW), (1531,478), (347,251), (548,152)):
        value = dict(selected_correct=correct, selected_known_nonmatch=nonmatch, selected_five_or_more_correct=long)
        counts[s] = {'A31_A3':value, 'A42':dict(value)}
    common = dict(completed=True, configuration=M.D.CONFIG, protocol=M.D.PROTOCOL, windows=1491,
        original_draws=1536, original_candidates=16321, optimizer_steps_this_run=0,
        new_geometry_labels_computed=True, model_or_backbone_forward=False, new_candidate_search=False,
        R4_policy_or_catalog_changed=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        selected_counts_per_species=counts, predeclared_TRAIN_gate=M.D.gate(counts), TRAIN_preparation_supported=False,
        eligible_reference_pairs=dict(arabidopsis_thaliana=dict(eligible=2172,five_or_more=854),
                                      oryza_sativa=dict(eligible=744,five_or_more=232)))
    report = dict(common, source_jobs={k:M.JOBS[k] for k in ('potentials','decode')},
        all_natural_search_and_score_files_closed_before_labels=True, original_R4_labels_exact=True,
        reference_injected_predictions=0)
    checked = dict(common, source_jobs={k:M.JOBS[k] for k in ('potentials','decode','readout')},
                   independent_neural_forward_replay=False, **{k:True for k in M.CPU_PROOFS})
    return report, checked


class EmissionLatticeCompetition(unittest.TestCase):
    def test_different_pools_chain_joins_unknowns_long_chains_and_exact_net_counts(self):
        saved, previous, labelled = fixture(); result = M.census(saved, previous, labelled)
        self.assertEqual((result['windows'],result['original_candidates'],result['new_selected_chains']),(2,8,12))
        for species in M.D.R.C.R.ALLOW:
            groups = result['per_species'][species]
            self.assertEqual(groups['1'],dict(retained=0,lost=2,gained=1))
            self.assertEqual(groups['correct_at_least_5_exons'],dict(retained=0,lost=1,gained=0))
            self.assertEqual(groups['0'],dict(retained=1,lost=0,gained=2))
            self.assertEqual(groups['-1'],dict(retained=1,lost=0,gained=1))
            self.assertEqual(result['lost_correct_selected_overlap_label_sets'][species]['correct'],
                             {'chains':2,'0':1,'-1,0':1})
            self.assertEqual(result['paired_window_pattern'][species],{'objective_increased_correct_decreased':1})
        self.assertEqual(result['pooled']['1'],dict(retained=0,lost=4,gained=2))
        self.assertEqual(result['minimum_selected_overlap_gain_minus_lost_chain_gain'],1.)

    def test_ties_are_not_forced_old_identity_and_touching_spans_are_not_blockers(self):
        saved, previous, labelled = fixture()
        for row in saved: row['scores'][0]=2.; row['gain']=12.
        result=M.census(saved,previous,labelled)
        self.assertEqual(result['minimum_selected_overlap_gain_minus_lost_chain_gain'],0.)
        self.assertEqual(result['pooled_lost_correct_selected_overlap_label_sets']['correct_at_least_5_exons'],
                         {'chains':2,'-1,0':2})
        self.assertEqual(result['lost_selected_positive_replacement_checks'],4)

    def test_identity_common_label_score_and_specific_replaceable_positive_rejected(self):
        saved, previous, labelled = fixture()
        with self.assertRaisesRegex(ValueError,'Complete'): M.census(saved[:-1],previous,labelled)
        wrong=deepcopy(saved); wrong[1]['original_ordinal']=0
        with self.assertRaisesRegex(ValueError,'identity'): M.census(wrong,previous,labelled)
        wrong=deepcopy(saved); wrong[0]['scores'][4]+=.1; wrong[0]['gain']+=.1
        with self.assertRaisesRegex(ValueError,'Common geometry'): M.census(wrong,previous,labelled)
        wrong=deepcopy(labelled); wrong[0]['labels']['A42'][4]=-1
        wrong[0]['counts']['A42']=M.V.counts(saved[0],wrong[0]['labels']['A42'])
        with self.assertRaisesRegex(ValueError,'Common geometry'): M.census(saved,previous,wrong)
        wrong=deepcopy(saved); wrong[0]['scores'][0]=1.; wrong[0]['gain']=11.
        # Total still beats the old set; local replacement independently finds the flaw.
        with self.assertRaisesRegex(ValueError,'would improve'): M.census(wrong,previous,labelled)

    def test_scientific_FAIL_analyzable_but_each_missing_CPU_proof_blocks_saved_rows(self):
        report, checked = closed()
        self.assertEqual(sum(M.validate_closed(report,checked).values()),2)
        for proof in M.CPU_PROOFS:
            wrong=dict(checked,**{proof:False})
            with self.subTest(proof=proof), patch.object(M.D.B.F.A,'require_jobs'), \
                    patch.object(M.D.R,'completed',side_effect=[report,wrong]), patch.object(M.D.R.C.R,'rows') as rows:
                with self.assertRaisesRegex(ValueError,'seven CPU'): M.main()
                rows.assert_not_called()


if __name__=='__main__':
    unittest.main()
