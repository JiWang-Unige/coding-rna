"""Two synthetic checks for the new A39 completion and four-arm adapter."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_null_anchored_set_attrition as M


class AttritionAdapter(unittest.TestCase):
    def test_four_arm_population_and_protected_transitions(self):
        chains = [[[0,2]], [[10,11],[14,15],[18,19],[22,23],[26,27]], [[10,28]],
                  [[40,41]], [[50,51]], [[60,61]], [[70,71]]]
        labels = [1,1,0,-1,1,0,1]
        score = {'A31_A3':[2,2,0,1,-1,1,0], 'A37':[-1,1,2,-1,-1,2,1],
                 'A38':[1,1,0,1,1,0,-1], 'A39':[-1,2,3,1,2,0,3]}
        selected = {'A31_A3':[0,1,3,5], 'A37':[2,5,6], 'A38':[0,1,3,4], 'A39':[2,3,4,6]}
        known = {'A31_A3':(4,2,1,1,1), 'A37':(3,1,2,0,0),
                 'A38':(4,3,0,1,1), 'A39':(4,2,1,1,0)}
        expected = {a: dict(candidates=7,pool_correct=4,pool_known_nonmatch=2,pool_unknown=1,
            **dict(zip(('selected','selected_correct','selected_known_nonmatch','selected_unknown',
                        'selected_five_or_more_correct'), v))) for a,v in known.items()}
        saved, previous, labelled, context, logistic = [], [], [], [], []
        for species in ('arabidopsis_thaliana','oryza_sativa'):
            window = dict(species=species, split='train', id=species+'|toy|0|+')
            rows = {a:dict(window=window, chains=chains, scores=score[a], selected=selected[a],
                          reference_used_in_scoring=False) for a in score}
            saved.append(rows['A39']); context.append(rows['A37']); logistic.append(rows['A38'])
            previous.append(dict(window=window, variants={'A31_A3':rows['A31_A3']}, reference_used_in_scoring=False))
            labelled.append(dict(window=window,labels=labels,counts=expected,eligible_references=4,eligible_five_or_more=1))
        actual = M.census(saved, previous, labelled, context, logistic)
        self.assertEqual((actual['windows'],actual['natural_candidates']),(2,14))
        for species in ('arabidopsis_thaliana','oryza_sativa'):
            self.assertEqual(actual['selected_counts_per_species'][species],expected)
            self.assertEqual(actual['eligible_reference_pairs'][species],dict(eligible=4,five_or_more=1))
            groups=actual['per_species'][species]
            self.assertEqual(groups['A31_A3_to_A39']['1']['lost_nonpositive'],1)
            self.assertEqual(groups['A31_A3_to_A39']['1']['lost_positive_not_selected'],1)
            self.assertEqual(groups['A31_A3_to_A39']['1']['gained'],2)
            self.assertEqual(groups['A31_A3_to_A39']['correct_at_least_5_exons']['lost_positive_not_selected'],1)
            self.assertEqual(groups['A37_to_A39']['-1']['gained'],1)
            self.assertEqual(groups['A38_to_A39']['-1']['kept'],1)
            for comparison,source in M.COMPARISONS.items():
                for key in ('1','0','-1','correct_at_least_5_exons'):
                    v=groups[comparison][key]
                    self.assertEqual(sum(v.get(k,0) for k in ('kept','gained','lost_nonpositive',
                        'lost_positive_not_selected','neither_nonpositive','neither_positive_not_selected')),v['pool'])
                    self.assertEqual(v['new_selected']-v['old_selected'],
                        v.get('gained',0)-v.get('lost_nonpositive',0)-v.get('lost_positive_not_selected',0))
        self.assertEqual(actual['pooled']['A31_A3_to_A39']['1']['lost_nonpositive'],2)

    def test_closed_FAIL_is_analysable_but_actual_composite_proof_is_required(self):
        fitted = dict(arm='A39', configuration=M.D.CONFIG, protocol=M.D.PROTOCOL, optimizer_steps=4608,
            feature_revision='b0483c23b6b63787b61a6d3a204a9b517d6ba345', completed=True,
            preflight_job=M.D.F.PREFLIGHT_JOB, trainable_parameters=612324,
            initialization=M.D.CONFIG['initialization'], final_checkpoint=str((M.D.F.OUT/'step_004608.pt').relative_to(M.I.ROOT)),
            reference_injected_training_candidates=0, **{k:False for k in ('supervised_checkpoint_read','generator_forward',
                'generator_updated','new_candidate_search','new_backbone_forward','new_R4_labelling',
                'DEV_used_for_selection','test_setaria_used','A33_outputs_opened')})
        counts = {}
        for species,correct,nonmatch,long in (('arabidopsis_thaliana',1531,347,548),('oryza_sativa',478,251,152)):
            v=dict(selected_correct=correct,selected_known_nonmatch=nonmatch,selected_five_or_more_correct=long)
            counts[species]={a:deepcopy(v) for a in ('A31_A3','A37','A38','A39')}
        gates=M.D.gate(counts)
        common=dict(completed=True,configuration=M.D.CONFIG,protocol=M.D.PROTOCOL,windows=1491,
            natural_candidates=16321,reference_injected_predictions=0,optimizer_steps_this_run=0,
            selected_counts_per_species=counts,eligible_reference_pairs={
                'arabidopsis_thaliana':dict(eligible=2172,five_or_more=854),'oryza_sativa':dict(eligible=744,five_or_more=232)},
            predeclared_TRAIN_gate=gates,TRAIN_preparation_supported=False,
            **{k:False for k in ('new_candidate_search','new_R4_labelling','DEV_used','test_setaria_used','A33_outputs_opened')})
        report=dict(common,source_fit_job=M.I.FIT_JOB,source_checks_job=M.I.CHECKS_JOB,
            source_preflight_job=M.D.F.PREFLIGHT_JOB,checkpoint=fitted['final_checkpoint'],
            source_context_comparator_jobs=dict(fit='13410061',readout='13410739',CPU='13410742'),
            source_logistic_comparator_jobs=dict(fit='13412630',readout='13412701',CPU='13412702'),
            source_baseline_jobs=dict(readout='13398521',CPU='13398524'),label_catalog_opened=False,new_backbone_forward=False,
            **{k:True for k in ('natural_score_file_closed_before_saved_label_arrays_open',
                'original_natural_geometry_order_sources_unchanged','original_saved_labels_reused',
                'model_weights_unchanged','all_parameters_frozen_and_gradients_absent')})
        checked=dict(common,source_jobs=dict(fit=M.I.FIT_JOB,readout=M.I.READOUT_JOB),
            model_or_backbone_forward=False,independent_neural_forward_replay=False,**{k:True for k in M.I.CPU_PROOFS})
        self.assertEqual(M.validate_completed(fitted,report,checked),gates)
        self.assertFalse(all(gates.values()))
        checked['actual_saved_fixed_composite_loss_and_importance_arithmetic_verified']=False
        with self.assertRaisesRegex(ValueError,'five CPU reconciliations'):
            M.validate_completed(fitted,report,checked)


if __name__ == '__main__':
    unittest.main()
