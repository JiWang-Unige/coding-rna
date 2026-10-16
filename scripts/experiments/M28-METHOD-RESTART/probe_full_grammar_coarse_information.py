#!/usr/bin/env python3
"""Read-only native TRAIN observation-information gradients at fixed9216."""
import argparse
import json
import math
import os
import time
import torch
import fit_full_grammar_fresh_joint24 as B
import fit_full_grammar_joint24 as J
import annotate_full_grammar_joint24 as A
from src.m28.full_grammar_coarse_channel import coarse_channel_nll

ROOT = B.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-INFORMATION-PROBE'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_coarse_information_probe.md'
CONFIG = dict(scope='two_original_TRAIN_views_first_coarse_bearing_fresh_draw_per_species',
    source_fit_job='13453777', source_cumulative_steps=9216, fresh_preparation_job='13453348',
    fresh_support_job='13453408', encoder='unchanged_joint24576_native_production',
    channel_examples=[0.1, 2/3], examples_not_fitted_or_selected=True,
    normalization='same_original_fine_plus_coarse_known_rows_in_all_three_losses',
    gradient='wrt_same_native_fifteen_role_logits_only_no_head_backward_or_update',
    entropy_shift_atol=1e-6, entropy_shift_rtol=1e-9, gradient_atol=1e-7,
    gradient_rtol=1e-7, max_exons=128, checkpoint_bp=256)
GROUPS = dict(background=[0], coding_body=[1,2,3], intron_body=[4,5,6],
    splice_coding_boundary=list(range(7,13)), start_stop=[13,14])


def summary(vector, denominator):
    value = vector.double()/denominator
    return dict(l2=float(torch.linalg.vector_norm(value)), max_abs=float(value.abs().max()),
        l2_by_role_group={k:float(torch.linalg.vector_norm(value[:,cols])) for k,cols in GROUPS.items()})


def main(selection='first'):
    began = time.perf_counter()
    out = OUT if selection == 'first' else ROOT/'outputs/M28-FULL-GRAMMAR-MAX-COARSE-INFORMATION-PROBE'
    config = CONFIG if selection == 'first' else dict(CONFIG,
        scope='two_original_TRAIN_views_max_coarse_positions_fresh_draw_per_species',
        selection='maximum_coarse_count_then_earliest_additional_step_not_native_loss_or_DEV')
    fitted = json.loads((J.FIT/'summary.json').read_text())
    if ((J.FIT/'STATUS').read_text().strip() != 'COMPLETED' or not fitted['completed']
            or fitted['job_id'] != '13453777' or fitted['configuration'] != J.CONFIG
            or fitted['cumulative_optimizer_steps'] != 9216
            or fitted['DEV_used'] or fitted['sealed_data_used']):
        raise ValueError('Completed fixed joint24576 final9216 source required')
    checkpoint = J.FIT/'step_009216.pt'
    state = torch.load(checkpoint, map_location='cpu', weights_only=True)
    A.validate_checkpoint(state)
    model = B.G.M28Core('C0',2048,128)
    model.load_state_dict(state['model'], strict=True)
    model.eval().requires_grad_(False).to('cuda')
    plan, windows, catalogs = B.plan_and_population()
    if selection == 'first':
        selected = [next(i for i,r in enumerate(plan) if r['species'] == s and r['coarse_positions'] > 0)
                    for s in B.Q.ALLOW]
    else:
        selected = [max((i for i,r in enumerate(plan) if r['species'] == s),
                        key=lambda i:(plan[i]['coarse_positions'],-i)) for s in B.Q.ALLOW]
    B.C.write(out/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'], configuration=config,
        protocol=PROTOCOL, selected_additional_steps=[i+1 for i in selected],
        checkpoint=str(checkpoint.relative_to(ROOT)), optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False))
    sequences, descriptors = B.C.P.read_train_DNA()
    backbone, tokenizer = B.G.load_native()
    records = []
    for ordinal in selected:
        expected = plan[ordinal]
        w = windows[expected['window_id']]
        chromosome, geometry, context, offset, targets = B.C.example(w,expected,sequences,catalogs)
        before = {k:torch.as_tensor(targets[k]).clone() for k in ('seg','region')}
        with torch.inference_mode():
            native = J.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
        logits = native.clone().requires_grad_()
        n = w['valid_bases']
        cuts = tuple(v for v in (offset,offset+n) if 0 < v < len(context))
        fine_only = dict(targets,region=torch.full_like(torch.as_tensor(targets['region']),-100))
        losses, gradients, counts = {},{},{}
        for name,labels,epsilon in (('informative',targets,.1),('uniform',targets,2/3),('fine_only',fine_only,.1)):
            result = coarse_channel_nll(context,logits,128,offset,offset+n,labels,epsilon,
                cuts=cuts,checkpoint_bp=256)
            value = result['loss']
            gradient = torch.autograd.grad(value,logits)[0]
            if not bool(torch.isfinite(value)) or float(value.detach()) < -1e-8 or not bool(torch.isfinite(gradient).all()):
                raise ValueError('Finite native loss and logit derivatives required')
            losses[name],gradients[name] = float(value.detach()),gradient
            counts[name] = {k:result[k] for k in ('fine_positions','coarse_positions','fully_unknown_positions')}
            del result,value
        denominator = expected['known_observation_denominator']
        if (counts['informative']['fine_positions'] != expected['fine_positions']
                or counts['informative']['coarse_positions'] != expected['coarse_positions']
                or counts['uniform'] != counts['informative']
                or denominator != max(1,expected['fine_positions']+expected['coarse_positions'])
                or counts['fine_only']['fine_positions'] != expected['fine_positions']
                or counts['fine_only']['coarse_positions'] != 0):
            raise ValueError('Original observations and fixed denominator differ')
        shift = expected['coarse_positions']*math.log(3.)
        error = losses['uniform']-losses['fine_only']-shift
        if (not math.isclose(losses['uniform']-losses['fine_only'],shift,
                abs_tol=config['entropy_shift_atol'],rel_tol=config['entropy_shift_rtol'])
                or not torch.allclose(gradients['uniform'],gradients['fine_only'],
                    atol=config['gradient_atol'],rtol=config['gradient_rtol'])):
            raise ValueError('Uniform coarse loss must retain entropy constant but fine-only gradient')
        if any(not torch.equal(torch.as_tensor(targets[k]),v) for k,v in before.items()):
            raise ValueError('Original targets changed during counterfactual calculation')
        delta = gradients['informative']-gradients['uniform']
        row = dict(window_id=w['id'],species=w['species'],seqid=w['seqid'],strand=w['strand'],
            additional_step=ordinal+1,core_start=w['start'],core_end=w['end'],
            native_context_bases=len(context),fine_positions=expected['fine_positions'],
            coarse_positions=expected['coarse_positions'],known_observation_denominator=denominator,
            raw_NLL=losses, normalized_NLL={k:v/denominator for k,v in losses.items()},
            uniform_minus_fine_entropy_shift=shift,entropy_shift_error=error,
            uniform_minus_fine_gradient_max_abs=float((gradients['uniform']-gradients['fine_only']).abs().max()),
            normalized_gradient={k:summary(v,denominator) for k,v in gradients.items()},
            informative_minus_uniform_normalized_gradient=summary(delta,denominator),
            original_targets_unchanged=True)
        records.append(row)
        print(json.dumps(row,allow_nan=False),flush=True)
        del logits,native,gradients,delta,targets,before,fine_only
    if (any(p.requires_grad or p.grad is not None for p in model.parameters())
            or any(p.requires_grad or p.grad is not None for p in backbone.parameters())
            or any(not torch.equal(p.cpu(),state['model'][k]) for k,p in model.state_dict().items())):
        raise ValueError('Fixed model and native encoder must remain unchanged without parameter gradients')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=config,protocol=PROTOCOL,
        checkpoint=str(checkpoint.relative_to(ROOT)), source_fit_job='13453777', source_steps=9216,
        records=records, original_TRAIN_DNA_sources=descriptors,
        native_forward=True,optimizer_steps_this_run=0,optimizer_instantiated=False,
        model_weights_unchanged=True,pretrained_encoder_frozen_no_grad=True,
        reference_used_by_encoder_or_head=False,original_targets_changed=False,
        predictions_modified=False,decoder_rerun=False,checkpoints_written=False,
        cached_features_used=False,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began,
        limitation='two_selected_coarse_bearing_TRAIN_views_at_one_fixed_head_not_performance_causality_or_population_necessity')
    B.C.write(out/'summary.json',report)
    if sum(p.stat().st_size for p in out.iterdir() if p.is_file()) > 256*1024:
        raise ValueError('Declared compact probe output ceiling exceeded')
    print(json.dumps(dict(completed=True,seconds=report['seconds']),allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--selection',choices=['first','max-coarse'],default='first')
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.selection)
