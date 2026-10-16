"""Two artificial A40 transition/closure checks; no research result access."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from test_m28_role_evidence_raw import receipts

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_role_evidence_attrition as M


class RoleEvidenceAttrition(unittest.TestCase):
    def test_three_arm_pairing_counts_unknown_and_long_chain_transitions(self):
        chains = [[[0,2]], [[10,11],[14,15],[18,19],[22,23],[26,27]], [[10,28]],
                  [[40,41]], [[50,51]], [[60,61]], [[70,71]]]
        labels = [1,1,0,-1,1,0,1]
        scores = dict(A31_A3=[2,2,0,1,-1,1,0], A39=[-1,2,3,1,2,0,3], A40=[1,1,0,-1,1,0,0])
        selected = dict(A31_A3=[0,1,3,5], A39=[2,3,4,6], A40=[0,1,4])
        known = dict(A31_A3=(4,2,1,1,1), A39=(4,2,1,1,0), A40=(3,3,0,0,1))
        expected = {a: dict(candidates=7,pool_correct=4,pool_known_nonmatch=2,pool_unknown=1,
            **dict(zip(('selected','selected_correct','selected_known_nonmatch','selected_unknown',
                        'selected_five_or_more_correct'), v))) for a,v in known.items()}
        saved, previous, labelled, prior = [], [], [], []
        for species in ('arabidopsis_thaliana','oryza_sativa'):
            window = dict(species=species,split='train',id=species+'|toy|0|+')
            rows = {a:dict(window=window,chains=chains,scores=scores[a],selected=selected[a],
                          reference_used_in_scoring=False) for a in scores}
            saved.append(rows['A40']); prior.append(rows['A39'])
            previous.append(dict(window=window,variants={'A31_A3':rows['A31_A3']},reference_used_in_scoring=False))
            labelled.append(dict(window=window,labels=labels,counts=expected,eligible_references=4,eligible_five_or_more=1))
        actual = M.census(saved,previous,labelled,prior)
        self.assertEqual((actual['windows'],actual['natural_candidates']),(2,14))
        for species in ('arabidopsis_thaliana','oryza_sativa'):
            self.assertEqual(actual['selected_counts_per_species'][species],expected)
            self.assertEqual(actual['eligible_reference_pairs'][species],dict(eligible=4,five_or_more=1))
            groups = actual['per_species'][species]
            self.assertEqual(groups['A39_to_A40']['1']['lost_nonpositive'],1)
            self.assertEqual(groups['A39_to_A40']['1']['gained'],2)
            self.assertEqual(groups['A39_to_A40']['correct_at_least_5_exons']['gained'],1)
            self.assertEqual(groups['A31_A3_to_A40']['-1']['lost_nonpositive'],1)
            self.assertEqual(groups['A31_A3_to_A40']['0']['lost_nonpositive'],1)
            for values in groups.values():
                for v in values.values():
                    self.assertEqual(sum(v.get(k,0) for k in ('kept','gained','lost_nonpositive',
                        'lost_positive_not_selected','neither_nonpositive','neither_positive_not_selected')),v['pool'])
                    self.assertEqual(v['new_selected']-v['old_selected'],
                        v.get('gained',0)-v.get('lost_nonpositive',0)-v.get('lost_positive_not_selected',0))
        broken = deepcopy(prior); broken[0]['chains'][0] = [[0,3]]
        with self.assertRaisesRegex(ValueError,'identity'):
            M.census(saved,previous,labelled,broken)

    def test_scientific_FAIL_can_be_analysed_but_missing_CPU_proof_blocks_score_access(self):
        fitted, report, checked = receipts()
        for value in (report, checked):
            for arms in value['selected_counts_per_species'].values():
                arms['A40'] = deepcopy(arms['A31_A3'])
            value['predeclared_TRAIN_gate'] = M.D.gate(value['selected_counts_per_species'])
            value['TRAIN_preparation_supported'] = False
        gates = M.validate_completed(fitted,report,checked)
        self.assertFalse(all(gates.values()))
        for proof in M.I.CPU_PROOFS:
            broken = {**checked,proof:False}
            with self.subTest(proof=proof), patch.object(M.D.F.A.A,'require_jobs'), \
                 patch.object(M.R,'completed',side_effect=[fitted,report,broken]), \
                 patch.object(M.R.C.R,'rows') as rows:
                with self.assertRaisesRegex(ValueError,'seven CPU reconciliations'):
                    M.main()
                rows.assert_not_called()


if __name__ == '__main__':
    unittest.main()
