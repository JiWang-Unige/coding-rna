"""Synthetic A40 raw routing; no research weights, FASTA, labels or DEV inputs."""
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch
from test_m28_null_anchored_set_raw import receipts as old_receipts
from test_m28_context_chain_raw import Generator as OldGenerator, manual_scores
from test_m28_query_decode import fixture

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_role_evidence_context as I
from verify_context_chain_train import verify_components


def receipts():
    fitted,prepared,checked = old_receipts()
    fitted.update(arm='A40',configuration=I.D.CONFIG,protocol=I.D.PROTOCOL,
        preflight_job=I.D.F.PREFLIGHT_JOB,trainable_parameters=613348,initialization=I.D.CONFIG['initialization'],
        final_checkpoint=str((I.D.F.OUT/'step_004608.pt').relative_to(I.ROOT)),
        student_supervised_checkpoint_read=False,frozen_supervised_A3_checkpoint_read=True,
        frozen_local_source_checkpoint=I.D.CONFIG['frozen_local_source'],frozen_local_source_optimizer_steps=0,
        frozen_local_source_weights_unchanged=True,frozen_local_source_gradients_absent=True,generator_decoding=False)
    for report in (prepared,checked):
        report.update(configuration=I.D.CONFIG,protocol=I.D.PROTOCOL,local_logit_receipts=1491)
        for values in report['selected_counts_per_species'].values():
            values['A40'] = values.pop('A39')
            values['A39'] = deepcopy(values['A31_A3'])
        report['predeclared_TRAIN_gate'] = I.D.gate(report['selected_counts_per_species'])
    prepared.update(source_fit_job=I.FIT_JOB,source_preflight_job=I.D.F.PREFLIGHT_JOB,
        source_baseline_jobs=dict(readout='13398521',CPU='13398524'),
        source_previous_arm_jobs=dict(fit='13415276',readout='13415309',CPU='13415319'),
        checkpoint=fitted['final_checkpoint'],frozen_local_source_checkpoint=I.D.CONFIG['frozen_local_source'],
        local_logit_files_closed_before_saved_label_arrays_open=True,frozen_local_weights_unchanged=True)
    checked.update(source_jobs=dict(fit=I.FIT_JOB,readout=I.READOUT_JOB),**{k:True for k in I.CPU_PROOFS})
    return fitted,prepared,checked


class Generator(OldGenerator):
    def shared(self,features,dna,valid):
        self.inputs.append((features,dna,valid))
        return features,None

    def segmentation(self,hidden):
        logits = torch.arange(len(hidden)*15,dtype=torch.float32).reshape(len(hidden),15)/100
        self.local_inputs = hidden
        self.local_logits = logits
        return logits


class Scorer:
    def __init__(self,generator,function=manual_scores):
        self.generator,self.function,self.inputs = generator,function,[]

    def __call__(self,features,dna,valid,chains,local_logits):
        if valid != [features.shape[1]] or tuple(dna.shape) != (1,I.M.WIDTH,5):
            raise AssertionError('Original synthetic oriented DNA and valid-length contract changed')
        if local_logits is not self.generator.local_logits or local_logits.dtype != torch.float32:
            raise AssertionError('New score must receive the same actual frozen source logits')
        self.inputs.append((features,dna,valid,chains,local_logits))
        return self.function(chains)


class RoleEvidenceRawTests(unittest.TestCase):
    def test_only_actual_A40_identity_and_all_new_source_proofs_accepted(self):
        fitted,prepared,checked = receipts(); I.validate_preparation(fitted,prepared,checked)
        for key,value in (('arm','A39'),('configuration',old_receipts()[0]['configuration']),
                ('preflight_job','13415165'),('frozen_local_source_checkpoint','wrong'),
                ('frozen_local_source_optimizer_steps',1),('student_supervised_checkpoint_read',True)):
            with self.subTest(key=key),self.assertRaises(ValueError):
                I.validate_preparation({**fitted,key:value},prepared,checked)
        with self.assertRaises(ValueError):
            I.validate_preparation(fitted,{**prepared,'source_fit_job':'13415276'},checked)
        for key in I.CPU_PROOFS:
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'CPU reconciliation'):
                I.validate_preparation(fitted,prepared,{**checked,key:False})

    def test_failed_stronger_gate_or_missing_role_proof_block_before_native_access(self):
        with tempfile.TemporaryDirectory() as directory,patch.object(I.D.F.A.A,'require_jobs'),\
             patch.object(I.D.R,'completed') as completed,patch.object(I,'load_models') as models,\
             patch.object(I.C.I,'load_native') as native,patch.object(I.C.I,'input_contigs') as inputs:
            for failed_gate in (True,False):
                fitted,prepared,checked = receipts()
                if failed_gate:
                    for report in (prepared,checked):
                        values = report['selected_counts_per_species']['arabidopsis_thaliana']
                        values['A40']['selected_correct'] = values['A31_A3']['selected_correct']
                        values['A39']['selected_correct'] = 0
                    # Cached true flags and improvement over A39 do not change strict A31 arithmetic.
                else:
                    checked['actual_saved_frozen_A3_role_inputs_and_original_A31_scores_verified'] = False
                completed.side_effect = [fitted,prepared,checked]
                out = Path(directory)/'unexecuted'
                with self.assertRaisesRegex(ValueError,'five-condition' if failed_gate else 'CPU reconciliation'):
                    I.main(out)
                self.assertFalse(out.exists())
            models.assert_not_called(); native.assert_not_called(); inputs.assert_not_called()

    def test_missing_separate_protocol_does_not_load_models_or_inputs(self):
        fitted,_,_ = receipts()
        with tempfile.TemporaryDirectory() as directory,patch.object(I,'ROOT',Path(directory)),\
             patch.object(I,'readiness',return_value=fitted),patch.object(I,'load_models') as models,\
             patch.object(I.C.I,'load_native') as native,patch.object(I.C.I,'input_contigs') as inputs:
            out = Path(directory)/'unexecuted'
            with self.assertRaisesRegex(ValueError,'separate fixed'):
                I.main(out)
            self.assertFalse(out.exists()); models.assert_not_called(); native.assert_not_called(); inputs.assert_not_called()

    def test_loaded_local_source_must_equal_actual_fit_not_another_emission_model(self):
        fitted,_,_ = receipts()
        generator = torch.nn.Linear(1,1).eval().requires_grad_(False)
        source = {k:v.clone() for k,v in generator.state_dict().items()}
        state = dict(arm='A40',configuration=I.D.CONFIG,protocol=I.D.PROTOCOL,step=4608,
            revision=I.C.I.F.REVISION,model={},frozen_local_source_checkpoint=I.D.CONFIG['frozen_local_source'],
            frozen_local_source_model=source)
        source['weight'] = source['weight']+1
        with patch.object(I.C,'load_generator',return_value=generator),patch.object(torch,'load',return_value=state),\
             patch.object(I.D.F.A,'new_model') as student,self.assertRaisesRegex(ValueError,'original frozen fit source'):
            I.load_models(fitted)
        student.assert_not_called()

    def test_same_local_logits_natural_pool_components_and_nonpositive_candidates(self):
        sequence = 'ATGAAATAA'; query,hidden = fixture(sequence,[(0,9)])
        generator = Generator(query); features = hidden[None]; dna = I.one_hot(sequence,I.M.WIDTH)[None]
        chains,sources = I.C.free_pool(generator,hidden,sequence)
        for gain in (3.,0.,-1.):
            scorer = Scorer(generator,lambda c:manual_scores(c,gain))
            row = I.local_record(generator,hidden,scorer,features,dna,sequence,0,'synthetic')
            self.assertTrue(chains); self.assertEqual(row['chains'],chains); self.assertEqual(row['decoder_sources'],sources)
            self.assertEqual(row['variant'],'A40'); self.assertEqual(row['score_kind'],I.SCORE_KIND)
            self.assertEqual(row['score_threshold'],0.); self.assertFalse(row['reference_used'])
            self.assertEqual(row['frozen_local_source_checkpoint'],I.D.CONFIG['frozen_local_source'])
            self.assertTrue(all(v == gain for v in row['score'])); verify_components({**row,'scores':row['score']})
        empty = I.local_record(generator,hidden,Scorer(generator),features,dna,'C'*len(sequence),0,'synthetic')
        self.assertEqual(empty['chains'],[]); self.assertEqual(empty['score'],[])

    def test_shared_once_for_both_proposals_and_role_input_and_minus_coordinates(self):
        oriented = 'ATGCGTTAGTTTAA'; query,features = fixture(oriented,[(0,4),(9,14)])
        generator = Generator(query); scorer = Scorer(generator); raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone,tokenizer,sequence,device):
            seen.append(sequence)
            return (torch.cat([torch.zeros(1,4),features]) if sequence == 'N'+oriented else torch.zeros(15,4))[None]
        with tempfile.TemporaryDirectory() as directory,patch.object(I.C.I,'encode_window',side_effect=encode):
            transcripts,scopes = I.annotate_contig(generator,scorer,None,None,'synthetic','chr',raw,0,Path(directory),'cpu')
            self.assertEqual(seen,[raw[:-1].upper()+'N','N'+oriented])
            self.assertEqual(len(generator.inputs),2); self.assertEqual(len(scorer.inputs),2)
            for old,new in zip(generator.inputs,scorer.inputs):
                self.assertIs(old[0],new[0]); self.assertIs(old[1],new[1]); self.assertEqual(old[2],new[2])
            minus = [t for t in transcripts if t['strand'] == '-']
            self.assertEqual(len(minus),1); self.assertEqual(minus[0]['CDS'],[(0,5,'2'),(10,14,'0')])
            self.assertEqual(minus[0]['score'],28.)
            row = next(I.C.I.saved_records(Path(directory)/scopes[1]['file']))
            verify_components({**row,'scores':row['score']}); self.assertFalse(row['reference_used'])


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
