"""A38 source routing and conditional execution only, not biological performance."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from test_m28_context_chain_raw import receipts as old_receipts, Generator, Scorer, manual_scores
from test_m28_context_chain_evaluation import decoded_receipt as old_decoded
from test_m28_query_decode import fixture

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_logistic_context as I
import evaluate_logistic_context as V


def receipts():
    fitted, prepared, checked = old_receipts()
    fitted.update(arm='A38',configuration=I.D.CONFIG,protocol=I.D.F.PROTOCOL,
        preflight_job=I.D.F.PREFLIGHT_JOB,
        final_checkpoint=str((I.D.F.OUT/'step_004608.pt').relative_to(I.ROOT)))
    for report in (prepared,checked):
        report.update(configuration=I.D.CONFIG,protocol=I.D.PROTOCOL)
        for values in report['selected_counts_per_species'].values():
            values['A38'] = values.pop('A37')
            values['A37'] = copy.deepcopy(values['A31_A3'])
        report['predeclared_TRAIN_gate'] = I.D.gate(report['selected_counts_per_species'])
    prepared.update(source_fit_job=I.FIT_JOB,source_checks_job=I.CHECKS_JOB,
        source_preflight_job=I.D.F.PREFLIGHT_JOB,checkpoint=fitted['final_checkpoint'])
    checked['source_jobs'] = dict(fit=I.FIT_JOB,readout=I.READOUT_JOB)
    return fitted,prepared,checked


def decoded_receipt():
    value = old_decoded()
    value.update(arm='A38',configuration=I.D.CONFIG,protocol=I.DEV_PROTOCOL,
        TRAIN_protocol=I.D.PROTOCOL,fit_protocol=I.D.F.PROTOCOL,
        source_TRAIN_jobs=dict(fit=I.FIT_JOB,readout=I.READOUT_JOB,CPU=I.CHECK_JOB),
        scorer_checkpoint=str((I.D.F.OUT/'step_004608.pt').relative_to(I.ROOT)),score_kind=I.SCORE_KIND)
    return value


class LogisticRawTests(unittest.TestCase):
    def test_only_fixed_A38_sources_and_objective_accepted(self):
        fitted,prepared,checked = receipts(); I.validate_preparation(fitted,prepared,checked)
        for key,value in (('arm','A37'),('configuration',old_receipts()[0]['configuration']),
                ('preflight_job','13409938'),('final_checkpoint','outputs/M28-A37-CONTEXT-CHAIN-FIT/step_004608.pt')):
            with self.subTest(key=key),self.assertRaises(ValueError):
                I.validate_preparation({**fitted,key:value},prepared,checked)
        for key,value in (('source_fit_job','13410061'),('source_checks_job','13410681'),
                ('source_preflight_job','13409938')):
            with self.subTest(key=key),self.assertRaises(ValueError):
                I.validate_preparation(fitted,{**prepared,key:value},checked)
        with self.assertRaises(ValueError):
            I.validate_preparation(fitted,prepared,{**checked,'source_jobs':dict(fit=I.FIT_JOB,readout='13410739')})

    def test_failed_original_gate_blocks_before_native_access_and_output(self):
        fitted,prepared,checked = receipts()
        for report in (prepared,checked):
            values = report['selected_counts_per_species']['arabidopsis_thaliana']
            values['A38']['selected_correct'] = values['A31_A3']['selected_correct']
            values['A37']['selected_correct'] = 0
        # Cached true flags cannot replace the original stronger-baseline arithmetic.
        with tempfile.TemporaryDirectory() as directory,patch.object(I.D.F.A,'require_jobs'),\
             patch.object(I.D.R,'completed',side_effect=[fitted,prepared,checked]),\
             patch.object(I,'load_models') as models,patch.object(I.C.I,'load_native') as native,\
             patch.object(I.C.I,'input_contigs') as inputs:
            out = Path(directory)/'unexecuted'
            with self.assertRaisesRegex(ValueError,'five-condition'):
                I.main(out)
            self.assertFalse(out.exists()); models.assert_not_called(); native.assert_not_called(); inputs.assert_not_called()

    def test_new_native_metadata_rejects_old_arm_and_source_before_reference(self):
        good = decoded_receipt(); V.validate_complete(good)
        for key,value in (('arm','A37'),('configuration',old_decoded()['configuration']),
                ('protocol',old_decoded()['protocol']),('score_kind',old_decoded()['score_kind']),
                ('source_TRAIN_jobs',old_decoded()['source_TRAIN_jobs']),
                ('scorer_checkpoint',old_decoded()['scorer_checkpoint'])):
            with self.subTest(key=key),patch.object(V.A,'reference_context') as reference,\
                 self.assertRaises(ValueError):
                V.validate_complete({**good,key:value})
            reference.assert_not_called()

    def test_actual_record_uses_A38_identity_with_unchanged_score_components(self):
        sequence = 'ATGAAATAA'; query,hidden = fixture(sequence,[(0,9)])
        generator = Generator(query); scorer = Scorer(lambda chains:manual_scores(chains,3.))
        row = I.local_record(generator,hidden,scorer,hidden[None],
            I.one_hot(sequence,I.M.WIDTH)[None],sequence,0,'synthetic')
        self.assertTrue(row['chains']); self.assertEqual(row['variant'],'A38')
        self.assertEqual(row['score_kind'],I.SCORE_KIND); self.assertEqual(row['score_threshold'],0.)
        self.assertFalse(row['reference_used']); self.assertTrue(all(v == 3. for v in row['score']))
        V.V38.verify_components({**row,'scores':row['score']})


if __name__ == '__main__':
    unittest.main()
