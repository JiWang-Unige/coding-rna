"""A41 raw routing with synthetic receipts and weights, no research data access."""
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch
from test_m28_query_decode import fixture
from test_m28_role_evidence_raw import Generator

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_emission_residual as I
from src.m28.emission_residual import EmissionResidual
from src.m28.local_path_score import local_path_scores


def receipts():
    fitted=dict(arm='A41',completed=True,configuration=I.D.CONFIG,protocol=I.D.PROTOCOL,
        optimizer_steps=4608,feature_revision=I.C.I.F.REVISION,preflight_job=I.D.F.PREFLIGHT_JOB,
        trainable_parameters=1806,initialization=I.D.CONFIG['initialization'],
        final_checkpoint=str((I.D.F.OUT/'step_004608.pt').relative_to(I.ROOT)),
        frozen_source_checkpoint=I.D.CONFIG['frozen_source'],frozen_source_optimizer_steps=0,
        frozen_source_weights_unchanged=True,frozen_source_gradients_absent=True,reference_injected_training_candidates=0,
        **{k:False for k in ('generator_decoding','generator_updated','new_candidate_search',
            'new_backbone_forward','new_R4_labelling','DEV_used_for_selection','test_setaria_used','A33_outputs_opened')})
    counts,refs={},{}
    for s,correct,wrong,five,eligible,long in (
            ('arabidopsis_thaliana',1531,347,548,2172,854),('oryza_sativa',478,251,152,744,232)):
        old=dict(selected_correct=correct,selected_known_nonmatch=wrong,selected_five_or_more_correct=five)
        counts[s]=dict(A31_A3=old,A41={**old,'selected_correct':correct+1,'selected_five_or_more_correct':five+1})
        refs[s]=dict(eligible=eligible,five_or_more=long)
    common=dict(completed=True,configuration=I.D.CONFIG,protocol=I.D.PROTOCOL,windows=1491,natural_candidates=16321,
        potential_receipts=1491,optimizer_steps_this_run=0,selected_counts_per_species=counts,eligible_reference_pairs=refs,
        predeclared_TRAIN_gate=I.D.gate(counts),TRAIN_preparation_supported=True,
        **{k:False for k in ('new_candidate_search','new_R4_labelling','DEV_used','test_setaria_used','A33_outputs_opened')})
    prepared=dict(deepcopy(common),source_fit_job=I.FIT_JOB,source_checks_job='13420720',
        source_preflight_job=I.D.F.PREFLIGHT_JOB,checkpoint=fitted['final_checkpoint'],original_draws=1536,
        reference_injected_predictions=0,label_catalog_opened=False,new_backbone_forward=False,
        **{k:True for k in ('all_scores_and_potential_files_closed_before_saved_labels_open',
            'original_natural_geometry_order_sources_unchanged','original_saved_labels_reused',
            'model_weights_unchanged','frozen_source_weights_unchanged','all_parameters_frozen_and_gradients_absent')})
    checked=dict(deepcopy(common),source_jobs=dict(fit=I.FIT_JOB,readout=I.READOUT_JOB),
        model_or_backbone_forward=False,independent_neural_forward_replay=False,**{k:True for k in I.CPU_PROOFS})
    return fitted,prepared,checked


class RecordingResidual(EmissionResidual):
    def __init__(self,hidden):
        super().__init__(hidden); self.inputs=[]

    def forward(self,hidden,local,chains,length):
        self.inputs.append((hidden,local,chains,length))
        return super().forward(hidden,local,chains,length)


class EmissionResidualRawTests(unittest.TestCase):
    def test_actual_fixed_residual_source_and_all_six_CPU_proofs_required(self):
        fitted,prepared,checked=receipts(); I.validate_preparation(fitted,prepared,checked)
        for key,value in (('arm','A40'),('optimizer_steps',1536),('frozen_source_checkpoint','wrong.pt'),
                ('frozen_source_optimizer_steps',1),('trainable_parameters',613348)):
            with self.subTest(key=key),self.assertRaises(ValueError):
                I.validate_preparation({**fitted,key:value},prepared,checked)
        for key in I.CPU_PROOFS:
            with self.subTest(proof=key),self.assertRaisesRegex(ValueError,'CPU reconciliation'):
                I.validate_preparation(fitted,prepared,{**checked,key:False})

    def test_failed_strict_gate_blocks_models_encoder_FASTA_and_output(self):
        with tempfile.TemporaryDirectory() as directory,patch.object(I.D.F.A.A,'require_jobs'),\
             patch.object(I.D.R,'completed') as completed,patch.object(I,'load_models') as models,\
             patch.object(I.C.I,'load_native') as native,patch.object(I.C.I,'input_contigs') as contigs:
            fitted,prepared,checked=receipts()
            for report in (prepared,checked):
                report['selected_counts_per_species']['arabidopsis_thaliana']['A41']['selected_correct']=1531
            completed.side_effect=[fitted,prepared,checked]
            out=Path(directory)/'unexecuted'
            with self.assertRaisesRegex(ValueError,'five-condition'): I.main(out)
            self.assertFalse(out.exists()); models.assert_not_called(); native.assert_not_called(); contigs.assert_not_called()

    def test_missing_separate_contract_blocks_before_models_or_inputs(self):
        fitted=receipts()[0]
        with tempfile.TemporaryDirectory() as directory,patch.object(I,'ROOT',Path(directory)),\
             patch.object(I,'readiness',return_value=fitted),patch.object(I,'load_models') as models,\
             patch.object(I.C.I,'load_native') as native,patch.object(I.C.I,'input_contigs') as contigs:
            out=Path(directory)/'unexecuted'
            with self.assertRaisesRegex(ValueError,'separate fixed'): I.main(out)
            self.assertFalse(out.exists()); models.assert_not_called(); native.assert_not_called(); contigs.assert_not_called()

    def test_original_baseline_and_species_denominators_cannot_be_rebased(self):
        fitted,prepared,checked=receipts()
        for change in ('baseline','denominator'):
            new,cpu=deepcopy(prepared),deepcopy(checked)
            for report in (new,cpu):
                if change=='baseline':
                    report['selected_counts_per_species']['arabidopsis_thaliana']['A31_A3']['selected_correct']-=1
                else:
                    report['eligible_reference_pairs']['arabidopsis_thaliana']['eligible']-=1
                    report['eligible_reference_pairs']['oryza_sativa']['eligible']+=1
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'baseline and full species denominators'):
                I.validate_preparation(fitted,new,cpu)

    def test_loaded_source_must_equal_frozen_fit_not_another_emission_model(self):
        fitted,_,_=receipts(); generator=torch.nn.Linear(1,1).eval().requires_grad_(False)
        state=dict(arm='A41',configuration=I.D.CONFIG,protocol=I.D.PROTOCOL,step=4608,
            revision=I.C.I.F.REVISION,frozen_source_checkpoint=I.D.CONFIG['frozen_source'],model={},
            frozen_source_model={k:v.clone() for k,v in generator.state_dict().items()})
        state['frozen_source_model']['weight']+=1
        with patch.object(I.C,'load_generator',return_value=generator),patch.object(torch,'load',return_value=state),\
             patch.object(I.D.F.A,'new_model') as residual,self.assertRaisesRegex(ValueError,'actual original frozen'):
            I.load_models(fitted)
        residual.assert_not_called()

    def test_same_natural_pool_zero_start_nonzero_residual_and_empty_pool(self):
        sequence='ATGAAATAA'; query,hidden=fixture(sequence,[(0,9)])
        generator=Generator(query); residual=RecordingResidual(hidden.shape[1]).eval()
        chains,sources=I.C.free_pool(generator,hidden,sequence)
        row=I.local_record(generator,hidden,residual,sequence,0,'synthetic')
        self.assertTrue(chains); self.assertEqual(row['chains'],chains); self.assertEqual(row['decoder_sources'],sources)
        self.assertEqual(row['score'],row['frozen_path_gain']); self.assertEqual(row['residual_path_gain'],[0.]*len(chains))
        expected=local_path_scores(generator.local_logits,chains,len(sequence)).tolist()
        self.assertEqual(row['score'],expected); self.assertEqual(row['variant'],'A41')
        with torch.no_grad(): residual.margin.bias.fill_(1.)
        changed=I.local_record(generator,hidden,residual,sequence,0,'synthetic')
        self.assertTrue(any(v!=0 for v in changed['residual_path_gain'])); I.verify_score_components(changed)
        empty=I.local_record(generator,hidden,residual,'C'*len(sequence),0,'synthetic')
        self.assertEqual(empty['chains'],[]); self.assertEqual(empty['score'],[]); I.verify_score_components(empty)
        good=dict(score=[float(2**24+1)],frozen_path_gain=[float(2**24)],residual_path_gain=[1.],null=0.)
        I.verify_score_components(good)
        with self.assertRaisesRegex(ValueError,'FP64'):
            I.verify_score_components({**good,'score':[float(torch.tensor(2**24+1,dtype=torch.float32))]})

    def test_encoder_and_shared_once_per_orientation_and_minus_CDS_phase(self):
        oriented='ATGCGTTAGTTTAA'; query,features=fixture(oriented,[(0,4),(9,14)])
        generator=Generator(query); residual=RecordingResidual(features.shape[1]).eval(); raw=I.M.rc(oriented).lower()+'r'
        seen=[]
        def encode(backbone,tokenizer,sequence,device):
            seen.append(sequence)
            return (torch.cat([torch.zeros(1,features.shape[1]),features]) if sequence=='N'+oriented
                    else torch.zeros(15,features.shape[1]))[None]
        with torch.no_grad(): residual.margin.bias.fill_(10.)
        with tempfile.TemporaryDirectory() as directory,patch.object(I.C.I,'encode_window',side_effect=encode):
            txs,scopes=I.annotate_contig(generator,residual,None,None,'synthetic','chr',raw,0,Path(directory),'cpu')
            self.assertEqual(seen,[raw[:-1].upper()+'N','N'+oriented])
            self.assertEqual(len(generator.inputs),2); self.assertEqual(len(residual.inputs),2)
            minus=[t for t in txs if t['strand']=='-']
            self.assertEqual(len(minus),1); self.assertEqual(minus[0]['CDS'],[(0,5,'2'),(10,14,'0')])
            row=next(I.C.I.saved_records(Path(directory)/scopes[1]['file'])); I.verify_score_components(row)
            self.assertFalse(row['reference_used']); self.assertEqual(row['score_kind'],I.SCORE_KIND)


if __name__=='__main__':
    torch.set_num_threads(2); unittest.main()
