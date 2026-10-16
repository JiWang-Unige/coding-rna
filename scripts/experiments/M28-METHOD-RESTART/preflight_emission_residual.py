#!/usr/bin/env python3
"""A41 zero-update frozen-context emission correction; no A40 outputs or fit."""
import json
import math
import time
import torch
import preflight_context_chain as A
from src.m28.core import select_nonoverlapping
from src.m28.emission_residual import EmissionResidual
from src.m28.local_path_score import local_path_scores
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss

ROOT, R = A.ROOT, A.R
OUT = ROOT/'outputs/M28-A41-EMISSION-RESIDUAL-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a41_emission_residual_protocol.md'
CONFIG = dict(hidden=128, seed=0, precision='FP32', score_accumulation='FP64',
    model='zero_affine_14_role_margin_correction_on_frozen_A3_context',
    trainable_parameters=1806, frozen_source_parameters=649457,
    frozen_source='outputs/M28-A3-FIRST-FIT/step_004608.pt',
    initialization='all_residual_weights_and_biases_zero_intergenic_correction_fixed_zero',
    candidate_pool='original_16321_natural_A3_chains_no_reference_injection',
    loss='equal_mean_known_BCE_and_mean_known_compatible_set_NLL',
    loss_unknown='exclude_-1_from_both_losses_keep_in_natural_competition',
    score='original_FP64_A31_gain_plus_FP64_sum_of_FP32_role_corrections',
    null=0., threshold=0., threshold_scan=False, optimizer_steps_this_run=0,
    new_context_learning=False, local_cross_entropy=False,
    population='original_1536_draws_1491_unique_TRAIN_2380_12566_1375',
    comparison='frozen_A31_A3_same_natural_pool_and_same_zero_gain_operating_point')


def new_model(device):
    torch.manual_seed(CONFIG['seed'])
    return EmissionResidual(CONFIG['hidden']).to(device)


def main():
    draws,index,windows,records,_,baseline,references=A.population()
    old={row['window']['id']:row['variants']['A31_A3']
         for row in R.C.R.rows(A.B.R.READOUT/'natural_scores.jsonl')}
    model=new_model('cuda').train()
    source=R.C.load_generator().to('cuda').eval().requires_grad_(False)
    if (sum(p.numel() for p in model.parameters())!=CONFIG['trainable_parameters'] or
            sum(p.numel() for p in source.parameters())!=CONFIG['frozen_source_parameters']):
        raise ValueError('Declared correction and frozen source structures differ')
    before={name:{k:v.detach().clone() for k,v in obj.state_dict().items()}
            for name,obj in (('residual',model),('source',source))}
    weights={d['window_id']:d['importance_weight_uniform_within_stratum'] for d in draws}
    began=time.perf_counter(); torch.cuda.reset_peak_memory_stats(); receipts=[]
    for witness in R.P.select_gradient_windows(records):
        ident=witness['window_id']; w=windows[ident]; rec=records[ident]
        chains,labels=[r['chain'] for r in rec['candidates']],[r['label'] for r in rec['candidates']]
        started=time.perf_counter()
        cached=torch.load(ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
        features,dna=R.C.R.Q.window_inputs(cached,w,'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Original cached inputs must remain detached')
        with torch.no_grad():
            hidden,_=source.shared(features,dna,[w['valid_bases']])
            hidden=hidden[0]; local=source.segmentation(hidden)
        model.zero_grad(set_to_none=True)
        score=model(hidden,local,chains,w['valid_bases'])
        original=local_path_scores(local,chains,w['valid_bases'])
        expected=torch.tensor(old[ident]['scores'],dtype=torch.float64,device='cuda')
        torch.testing.assert_close(original,expected,rtol=1e-12,atol=1e-8)
        if (not torch.equal(score['logits'],original) or
                not bool((score['residual_local_logits']==0).all()) or
                not bool((score['residual_path_gain']==0).all()) or
                any(not bool(torch.isfinite(v).all()) for v in score.values()) or
                select_nonoverlapping(chains,score['logits'].tolist())!=old[ident]['selected']):
            raise ValueError('Finite zero residual must preserve original gain and selection exactly')
        result=null_anchored_chain_set_loss(score,chains,labels)
        weighted=result['loss']*weights[ident]
        if not bool(torch.isfinite(weighted)):
            raise ValueError('Nonfinite actual zero-start fixed composite objective')
        score['logits'].retain_grad(); weighted.backward()
        if any(float(score['logits'].grad[i])!=0. for i,y in enumerate(labels) if y==-1):
            raise ValueError('Unknown direct gain gradient must remain zero')
        if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in model.parameters()):
            raise ValueError('Missing/nonfinite affine correction gradient')
        norm=math.sqrt(sum(float(p.grad.double().square().sum()) for p in model.parameters()))
        if hidden.requires_grad or local.requires_grad or any(p.requires_grad or p.grad is not None
                for p in source.parameters()):
            raise ValueError('Original frozen source received gradient')
        torch.cuda.synchronize()
        receipts.append(dict(**witness,weighted_loss=float(weighted.detach()),
            unweighted_mean_known_BCE=float(result['mean_known_BCE'].detach()),
            unweighted_mean_known_set_NLL=float(result['mean_known_set_NLL'].detach()),
            unweighted_fixed_composite_loss=float(result['loss'].detach()),
            residual_gradient_norm=norm,
            nonzero_candidate_gain_derivatives=int((score['logits'].grad!=0).sum()),
            nonzero_gradient_role_rows=int((model.margin.weight.grad.double().norm(dim=1)>0).sum()),
            original_saved_A31_score_max_error=float((original-expected).abs().max()) if len(chains) else 0.,
            zero_residual_gain_and_original_selection_exact=True,
            unknown_direct_gain_gradient_zero=True,frozen_source_gradients_absent=True,
            load_forward_backward_seconds=time.perf_counter()-started))
        print(json.dumps(receipts[-1],allow_nan=False),flush=True)
        del cached,features,dna,hidden,local,score,original,expected,result,weighted
    if (len(receipts)!=8 or not all(any(r['species']==s and r['residual_gradient_norm']>0 for r in receipts)
            for s in R.C.R.ALLOW) or any(not torch.equal(v,before[name][k])
            for name,obj in (('residual',model),('source',source)) for k,v in obj.state_dict().items())):
        raise ValueError('Original eight metadata cases and unchanged states required')
    R.C.R.write_summary(OUT,dict(completed=True,arm='A41',protocol=PROTOCOL,configuration=CONFIG,
        windows=1491,original_draws=1536,natural_candidates=16321,
        natural_counts=R.totals(records,records),
        replayed_original_draw_counts=R.totals(records,[d['window_id'] for d in draws]),
        eligible_reference_pairs=references,stronger_baseline_counts=baseline['selected_counts_per_species'],
        original_natural_identity_labels_and_stronger_baseline_reconciled=True,
        initial_A31_gain_and_selection_exact=True,residual_weights_unchanged=True,
        frozen_source_weights_unchanged=True,frozen_source_gradients_absent=True,
        actual_both_species_residual_gradients_nonzero=True,
        zero_gradient_metadata_cases=sum(r['residual_gradient_norm']==0 for r in receipts),
        optimizer_steps_this_run=0,checkpoints_written=False,
        frozen_supervised_A3_checkpoint_read=True,original_TRAIN_catalog_metadata_read=True,
        reference_injected_training_candidates=0,generator_decoding=False,generator_updated=False,
        new_candidate_search=False,new_backbone_forward=False,new_R4_labelling=False,
        A40_outputs_opened=False,DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        records=receipts,seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Zero-update affine residual readiness only. Frozen supervised prior and extra '
            'trainable capacity are explicit; no learned improvement or independent generalization.'))


if __name__=='__main__':
    torch.set_num_threads(2); main()
