"""A41 adapters and receipt arithmetic; artificial inputs, no research tensors."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
import torch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import readout_emission_residual_train as D
import verify_emission_residual_train as V


def fixture():
    chains=[[[0,6]],[[10,13],[14,17],[18,21],[22,25],[26,29]],[[40,46]],[[50,56]]]
    labels=[1,1,-1,0]; saved=[]; previous=[]; labelled=[]
    for species in D.R.C.R.ALLOW:
        window=dict(species=species,split='train',id=species+'|toy|0|+')
        old=dict(chains=chains,scores=[1.,2.,3.,4.],selected=[0,1,2,3])
        new=dict(window=window,chains=chains,scores=[1.,2.,-1.,0.],selected=[0,1],reference_used_in_scoring=False)
        saved.append(new); previous.append(dict(window=window,variants={'A31_A3':old},reference_used_in_scoring=False))
        labelled.append(dict(window=window,labels={'A31_A3':labels},counts={'A31_A3':D.R.selected_counts(old,labels)},
            eligible_references=5,eligible_five_or_more=2))
    return saved,previous,labelled


class EmissionResidualReadout(unittest.TestCase):
    def test_pairing_unknown_full_denominators_and_strict_five_conditions(self):
        saved,previous,labelled=fixture(); rows,counts,refs=D.paired_readout(saved,previous,labelled)
        self.assertEqual(len(rows),2)
        for species in D.R.C.R.ALLOW:
            self.assertEqual(refs[species],dict(eligible=5,five_or_more=2))
            self.assertEqual(counts[species]['A31_A3']['selected_unknown'],1)
            self.assertEqual(counts[species]['A41']['pool_unknown'],1)
            self.assertEqual(counts[species]['A41']['selected_unknown'],0)
        gates=D.gate(counts)
        self.assertEqual((len(gates),sum(gates.values())),(5,2))

    def test_missing_reordered_reference_input_or_stale_selection_rejected(self):
        saved,previous,labelled=fixture()
        with self.assertRaisesRegex(ValueError,'Complete'): D.paired_readout(saved[:-1],previous,labelled)
        with self.assertRaisesRegex(ValueError,'identity'): D.paired_readout(saved[::-1],previous,labelled)
        broken=deepcopy(saved); broken[0]['reference_used_in_scoring']=True
        with self.assertRaisesRegex(ValueError,'identity'): D.paired_readout(broken,previous,labelled)
        broken=deepcopy(saved); broken[0]['selected'].append(3)
        with self.assertRaisesRegex(ValueError,'reconcile'): D.paired_readout(broken,previous,labelled)

    def test_FP64_saved_potential_sum_not_FP32_rounding(self):
        window=dict(id='toy',valid_bases=6); frozen=torch.zeros(6,15); residual=torch.zeros(6,15)
        frozen[0,13]=2**24; residual[1,13]=1.
        row=dict(window=window,chains=[[[0,6]]],scores=[float(2**24+1)],
            frozen_path_gain=[float(2**24)],residual_path_gain=[1.],null=0.)
        receipt=dict(window=window,source_checkpoint=D.CONFIG['frozen_source'],revision=D.R.C.I.F.REVISION,
            frozen_local_logits=frozen,residual_local_logits=residual)
        self.assertEqual(V.verify_potentials(row,receipt,[float(2**24)]),0.)
        broken=deepcopy(row); broken['scores']=[float(2**24)]
        with self.assertRaises(AssertionError): V.verify_potentials(broken,receipt,[float(2**24)])
        broken=deepcopy(receipt); broken['residual_local_logits'][0,0]=1.
        with self.assertRaisesRegex(ValueError,'receipt'): V.verify_potentials(row,broken,[float(2**24)])

    def test_missing_zero_update_preflight_proof_rejected(self):
        flags=('original_natural_identity_labels_and_stronger_baseline_reconciled','initial_A31_gain_and_selection_exact',
            'residual_weights_unchanged','frozen_source_weights_unchanged','frozen_source_gradients_absent',
            'actual_both_species_residual_gradients_nonzero')
        report=dict(completed=True,configuration=D.F.A.CONFIG,protocol=D.PROTOCOL,windows=1491,original_draws=1536,
            natural_candidates=16321,eligible_reference_pairs={'eligible':2916,'five_or_more':1086},records=[{}]*8,
            optimizer_steps_this_run=0,checkpoints_written=False,**{k:True for k in flags},
            **{k:False for k in ('generator_decoding','generator_updated','new_candidate_search','new_backbone_forward',
                'new_R4_labelling','A40_outputs_opened','DEV_used','test_setaria_used','A33_outputs_opened')})
        D.F.validate_preflight(report)
        for flag in flags:
            with self.subTest(flag=flag),self.assertRaisesRegex(ValueError,'preparation'):
                D.F.validate_preflight({**report,flag:False})

    def test_optimizer_steps_and_residual_only_population(self):
        named=list(D.F.A.new_model('cpu').named_parameters())
        opt=dict(param_groups=[dict(params=[0,1],lr=D.CONFIG['lr'],weight_decay=D.CONFIG['weight_decay'],
            betas=tuple(D.CONFIG['betas']),eps=D.CONFIG['eps'])],state={})
        for ident,(_,param) in enumerate(named):
            opt['state'][ident]=dict(step=torch.tensor(1536.),exp_avg=torch.zeros_like(param),exp_avg_sq=torch.zeros_like(param))
        V.verify_optimizer(opt,named,1536)
        wrong=deepcopy(opt); wrong['state'][0]['step']=torch.tensor(1535.)
        with self.assertRaisesRegex(ValueError,'exposure'): V.verify_optimizer(wrong,named,1536)
        wrong=deepcopy(opt); wrong['param_groups'][0]['params'].append(2)
        with self.assertRaisesRegex(ValueError,'residual-only'): V.verify_optimizer(wrong,named,1536)


if __name__=='__main__':
    unittest.main()
