import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_b2_recovery as D


class RecoveryDiagnosisTests(unittest.TestCase):
    def test_full_pair_required_and_checkpoint_mismatch_rejected(self):
        decoded = dict(inference_complete=True, windows=8690, reference_used=False,
                       checkpoint='B2/step_004608.pt', base_checkpoint='R5/B1/step_004608.pt',
                       feature_revision='fixed', control_replay=dict(discrete_exact=True),
                       paired_common_gain=dict(exact_equality=True, single_union_readout=True))
        evaluated = dict(experiment='M28-B2', methods={'B2_early': {}, 'B2_late': {}},
                         paired_inference_cost=copy.deepcopy(decoded))
        D.require_complete(decoded, evaluated)
        for field,value in [('inference_complete',False), ('windows',8689), ('reference_used',True)]:
            broken = dict(decoded, **{field:value})
            with self.assertRaises(ValueError):
                D.require_complete(broken, evaluated)
        broken = copy.deepcopy(evaluated)
        broken['paired_inference_cost']['checkpoint'] = 'different/step_004608.pt'
        with self.assertRaises(ValueError):
            D.require_complete(decoded, broken)

    def test_first_loss_and_zero_gain_are_distinct_states(self):
        cases = [(False,False,False,None,'not_free'),
                 (True,False,False,None,'lost_at_owner'),
                 (True,True,False,-1.,'owner_nonpositive'),
                 (True,True,False,0.,'owner_nonpositive'),
                 (True,True,False,1.,'owner_positive_not_selected'),
                 (True,True,True,1.,'selected')]
        for free,owner,selected,gain,expected in cases:
            self.assertEqual(D.recovery_state(free,owner,selected,gain), expected)
        with self.assertRaises(ValueError):
            D.recovery_state(True,True,False,None)

    def test_chain_identity_not_phase_or_id_and_unsupported_reference_retained(self):
        ref = dict(species='arabidopsis_thaliana', id='primary_minus', gene_id='gene1',
                   seqid='DEV', strand='-', CDS=[(20,29,'1'),(40,46,'2')],
                   M28_unsupported_strata_not_excluded=['longer_than_window'])
        chain = ((20,29),(40,46))
        sets = {'early': dict(zip(D.STAGES, [{chain},{chain},{chain}])),
                'late': dict(zip(D.STAGES, [{chain},set(),set()]))}
        row = D.reference_row(ref, sets, {'early':{chain:2.},'late':{}})
        self.assertEqual(row['early_state'], 'selected')
        self.assertEqual(row['late_state'], 'lost_at_owner')
        self.assertEqual(row['unsupported_strata'], 'longer_than_window')
        self.assertEqual((row['exon_count'],row['span_bp'],row['coding_bp']), (2,26,15))
        self.assertIsNone(row['late_owner_gain'])
        changed = dict(ref, id='different_ID', CDS=[(20,29,'0'),(40,46,'0')])
        self.assertEqual(D.reference_row(changed, sets, {'early':{chain:2.},'late':{}})['early_state'], 'selected')

    def test_pair_matrix_retains_both_gains_losses_and_neither(self):
        records = []
        for early,late in [(True,True),(True,False),(False,True),(False,False)]:
            row = dict(early_state='selected' if early else 'not_free',
                       late_state='selected' if late else 'not_free')
            row.update({'early_'+s:early for s in D.STAGES})
            row.update({'late_'+s:late for s in D.STAGES})
            records.append(row)
        report = D.transition_summary(records)
        self.assertEqual(report['reference_chains'], 4)
        for counts in report['stages'].values():
            self.assertEqual(counts, dict(both=1,early_only=1,late_only=1,neither=1))
        self.assertEqual(sum(r['reference_chains'] for r in report['late_to_early_states']), 4)


if __name__ == '__main__':
    unittest.main()
