#!/usr/bin/env python3
"""Fixed fresh/joint continuation: uniform coarse, original denominator."""
import argparse
from collections import Counter
import copy
import json
import math
import os
import time
import torch
import fit_full_grammar_fresh_joint24 as B
import annotate_full_grammar_uniform_coarse as A
from src.m28.full_grammar_coarse_channel import coarse_channel_nll

ROOT, C, G, Q, P, J = B.ROOT, B.C, B.G, B.Q, B.P, B.J
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-FIT'
PREFLIGHT = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-PREFLIGHT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-VERIFICATION'
PROTOCOL = A.PROTOCOL
CONFIG = dict(B.CONFIG, **A.IDENTITY)


def objective(context, logits, offset, n, targets):
    cuts = tuple(x for x in (offset,offset+n) if 0 < x < len(context))
    return coarse_channel_nll(context,logits,128,offset,offset+n,targets,CONFIG['epsilon'],
        cuts=cuts,checkpoint_bp=256)['loss']


def native_identity_check():
    state = dict(configuration=CONFIG,step=9216,epoch=6,protocol=PROTOCOL,
                 encoder_revision=Q.F.REVISION)
    A.validate_checkpoint(state)
    for key,value in (('epsilon',.1),('sampling_seed',2),('unique_windows',1491),('encoder_blocks_bp',6144)):
        altered = copy.deepcopy(state)
        altered['configuration'][key] = value
        try:
            A.validate_checkpoint(altered)
        except ValueError:
            pass
        else:
            raise AssertionError('Raw adapter accepted changed '+key)
    return 'uniform_fresh_joint_identity_accepts_fixed_rejects_coarse_sampler_and_encoder_drift'


def preflight():
    began = time.perf_counter()
    fixture = native_identity_check()
    plan,windows,catalogs = B.plan_and_population()
    selected = [max((i for i,r in enumerate(plan) if r['species'] == s),
                    key=lambda i:(plan[i]['coarse_positions'],-i)) for s in Q.ALLOW]
    C.write(PREFLIGHT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=PROTOCOL,selected_additional_steps=[i+1 for i in selected],starting_fit_job='13434603'))
    state = C.starting_state()
    sequences,descriptors = C.P.read_train_DNA()
    backbone,tokenizer = G.load_native()
    records = []
    for ordinal in selected:
        expected = plan[ordinal]
        w = windows[expected['window_id']]
        chromosome,geometry,context,offset,targets = C.example(w,expected,sequences,catalogs)
        model,optimizer = C.load_model_optimizer(state)
        joined = J.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
        raw = objective(context,joined,offset,w['valid_bases'],targets)
        fine_only = dict(targets,region=torch.full_like(torch.as_tensor(targets['region']),-100))
        fine_loss = objective(context,joined,offset,w['valid_bases'],fine_only)
        uniform_gradient = torch.autograd.grad(raw,joined,retain_graph=True)[0]
        fine_gradient = torch.autograd.grad(fine_loss,joined)[0]
        shift = expected['coarse_positions']*math.log(3.)
        error = float(raw.detach()-fine_loss.detach())-shift
        maximum_difference = float((uniform_gradient-fine_gradient).abs().max())
        if (not math.isclose(float(raw.detach()-fine_loss.detach()),shift,abs_tol=1e-6,rel_tol=1e-9)
                or not torch.allclose(uniform_gradient,fine_gradient,atol=1e-7,rtol=1e-7)):
            raise ValueError('Actual source4608 uniform loss and fine-only derivatives differ')
        total,norm = C.update(model,optimizer,backbone,raw,expected['objective_multiplier'])
        if (any(int(s['step']) != 4609 for s in optimizer.state.values())
                or all(torch.equal(p.cpu(),state['model'][k]) for k,p in model.state_dict().items())):
            raise ValueError('Independent disposable original-AdamW update4609 required')
        fixture_state = dict(configuration=CONFIG,protocol=PROTOCOL,step=9216,epoch=6,
            encoder_revision=Q.F.REVISION,model=model.state_dict())
        inference_model = A.model_from_state(fixture_state,'cuda')
        model.eval()
        with torch.inference_mode():
            trained = J.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
            restored = J.raw_view(inference_model,backbone,tokenizer,chromosome,geometry,w['strand'])
        if not torch.equal(trained,restored):
            raise ValueError('Updated native logits differ after strict uniform raw-adapter load')
        row = dict(additional_step=ordinal+1,window_id=w['id'],species=w['species'],
            fine_positions=expected['fine_positions'],coarse_positions=expected['coarse_positions'],
            known_observation_denominator=expected['known_observation_denominator'],
            uniform_raw_NLL=float(raw.detach()),fine_only_raw_NLL=float(fine_loss.detach()),
            entropy_shift=shift,entropy_shift_error=error,uniform_minus_fine_gradient_max_abs=maximum_difference,
            weighted_total=total,gradient_norm_before_clip=norm,
            initial_optimizer_step=4608,final_disposable_optimizer_step=4609,
            actual_native_logits_raw_adapter_exact=True)
        records.append(row)
        print(json.dumps(row,allow_nan=False),flush=True)
        del model,optimizer,inference_model,joined,raw,fine_loss,fine_gradient,uniform_gradient,trained,restored,fixture_state
    C.write(PREFLIGHT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG,protocol=PROTOCOL,fixture=fixture,records=records,
        original_TRAIN_DNA_sources=descriptors,disposable_updates_only=True,optimizer_steps_this_run=2,
        formal_training_started=False,checkpoints_written=False,raw_adapter_fixture_metadata_not_formal9216=True,
        actual_native_logits_raw_adapter_exact=True,uniform_entropy_and_fine_only_gradient_match=True,
        saved_source_weights_and_optimizer_unchanged=True,pretrained_encoder_frozen_no_grad=True,
        cached_features_used=False,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began))


def fit(preflight_job):
    began = time.perf_counter()
    ready = json.loads((PREFLIGHT/'summary.json').read_text())
    if ((PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != preflight_job or ready['configuration'] != CONFIG
            or ready['optimizer_steps_this_run'] != 2 or not ready['actual_native_logits_raw_adapter_exact']
            or not ready['uniform_entropy_and_fine_only_gradient_match']
            or ready['checkpoints_written'] or ready['DEV_used'] or ready['sealed_data_used']):
        raise ValueError('Completed matching uniform-coarse native preflight required')
    plan,windows,catalogs = B.plan_and_population()
    C.write(FIT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=PROTOCOL,preflight_job=preflight_job,starting_fit_job='13434603',support_job='13453408'))
    state = C.starting_state()
    model,optimizer = C.load_model_optimizer(state)
    sequences,descriptors = C.P.read_train_DNA()
    backbone,tokenizer = G.load_native()
    epochs = []
    with (FIT/'training.jsonl').open('x',buffering=512*1024) as handle:
        for pass_index in range(3):
            pass_started = time.perf_counter()
            sums,species_counts = Counter(),Counter()
            for expected in plan[pass_index*1536:(pass_index+1)*1536]:
                started = time.perf_counter()
                w = windows[expected['window_id']]
                chromosome,geometry,context,offset,targets = C.example(w,expected,sequences,catalogs)
                optimizer.zero_grad(set_to_none=True)
                joined = J.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
                raw = objective(context,joined,offset,w['valid_bases'],targets)
                total,norm = C.update(model,optimizer,backbone,raw,expected['objective_multiplier'])
                step = 4608+expected['additional_step']
                row = dict(step=step,epoch=pass_index+4,
                    **{k:expected[k] for k in ('additional_step','pass_index','draw_ordinal','window_id','species',
                        'importance_weight','known_observation_denominator','fine_positions','coarse_positions','fully_unknown_positions')},
                    raw_NLL=float(raw.detach()),normalized_NLL=float(raw.detach())/expected['known_observation_denominator'],
                    weighted_total=total,gradient_norm_before_clip=norm,seconds=time.perf_counter()-started)
                handle.write(json.dumps(row,allow_nan=False)+'\n')
                sums.update({k:row[k] for k in ('raw_NLL','normalized_NLL','weighted_total')})
                species_counts[w['species']] += 1
                if step % 64 == 0:
                    handle.flush()
                    print(json.dumps(dict(step=step,final_step=9216,recent=row,
                        elapsed_seconds=time.perf_counter()-began)),flush=True)
                del raw,joined,targets
            checkpoint = FIT/('step_%06d.pt' % step)
            if checkpoint.exists():
                raise FileExistsError(checkpoint)
            torch.save(dict(configuration=CONFIG,protocol=PROTOCOL,epoch=pass_index+4,step=step,
                encoder_revision=Q.F.REVISION,model=model.state_dict(),optimizer=optimizer.state_dict()),checkpoint)
            epochs.append(dict(epoch=pass_index+4,steps=1536,by_species=dict(species_counts),
                mean_losses={k:v/1536 for k,v in sums.items()},seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 9216 or any(int(s['step']) != 9216 for s in optimizer.state.values()):
        raise ValueError('Complete fixed uniform-coarse4608 continuation required')
    if sum(p.stat().st_size for p in FIT.iterdir() if p.is_file()) > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared fit artifact ceiling exceeded')
    C.write(FIT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=PROTOCOL,preflight_job=preflight_job,starting_fit_job='13434603',preparation_job='13453348',
        support_job='13453408',cumulative_optimizer_steps=9216,optimizer_steps_this_run=4608,trained_parameters=366191,
        fresh_draws=4608,fresh_unique_windows=4219,draws_per_pass=1536,epochs=epochs,
        final_checkpoint=str(checkpoint.relative_to(ROOT)),original_TRAIN_DNA_sources=descriptors,
        fixed_fresh_draw_order_and_unit_importance_preserved=True,all_fixed_draws_retained=True,
        original_targets_changed=False,encoder_condition='JOINT24576',native_CDS_priors_used=False,
        context_or_inference_grammar_changed=False,only_intervention='uniform_coarse_information_in_additional_fit',
        pretrained_encoder_frozen_no_grad=True,saved_source_weights_and_optimizer_unchanged=True,
        cached_features_used=False,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began,
        interpretation='coarse_information_additional_learning_from_already_informative_trained_start_not_total_necessity_or_generalization'))


def verify(fit_job):
    began = time.perf_counter()
    report = json.loads((FIT/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
            or report['job_id'] != fit_job or report['configuration'] != CONFIG
            or report['optimizer_steps_this_run'] != 4608 or report['cumulative_optimizer_steps'] != 9216
            or report['DEV_used'] or report['sealed_data_used']):
        raise ValueError('Completed matching uniform-coarse final fit required')
    plan,_,_ = B.plan_and_population()
    rows = list(P.rows(FIT/'training.jsonl'))
    if len(rows) != 4608:
        raise ValueError('Complete fresh exposure log required')
    for i,(row,expected) in enumerate(zip(rows,plan)):
        if (row['step'] != 4609+i or row['epoch'] != 4+i//1536
                or any(row[k] != expected[k] for k in ('additional_step','pass_index','draw_ordinal','window_id',
                    'species','importance_weight','known_observation_denominator','fine_positions','coarse_positions','fully_unknown_positions'))):
            raise ValueError('Actual draw, observation or update order differs')
        values = [row[k] for k in ('raw_NLL','normalized_NLL','weighted_total','gradient_norm_before_clip','seconds')]
        raw = row['raw_NLL']
        if (not all(math.isfinite(x) and x >= -1e-8 for x in values)
                or not math.isclose(row['normalized_NLL'],raw/expected['known_observation_denominator'],rel_tol=1e-12,abs_tol=1e-12)
                or not math.isclose(row['weighted_total'],raw*expected['objective_multiplier'],rel_tol=1e-12,abs_tol=1e-12)):
            raise ValueError('Logged fixed-denominator uniform objective arithmetic differs')
        if expected['fine_positions']+expected['coarse_positions'] == 0:
            if abs(raw) > 1e-8 or row['gradient_norm_before_clip'] > 1e-6:
                raise ValueError('All-unknown core gained observation gradient')
    previous = C.starting_state()['model']
    for index in range(3):
        step = 6144+index*1536
        state = torch.load(FIT/('step_%06d.pt' % step),map_location='cpu',weights_only=True)
        if (state['configuration'] != CONFIG or state['protocol'] != PROTOCOL or state['step'] != step
                or state['epoch'] != index+4 or state['encoder_revision'] != Q.F.REVISION):
            raise ValueError('Checkpoint identity differs')
        model = G.M28Core('C0',2048,128)
        model.load_state_dict(state['model'],strict=True)
        if (not all(v.dtype == torch.float32 and bool(torch.isfinite(v).all()) for v in state['model'].values())
                or all(torch.equal(v,previous[k]) for k,v in state['model'].items())):
            raise ValueError('Nonfinite or unchanged checkpoint weights')
        optimizer = state['optimizer']
        ids = [p for group in optimizer['param_groups'] for p in group['params']]
        parameters = list(model.parameters())
        if len(ids) != len(parameters) or len(ids) != len(optimizer['state']):
            raise ValueError('Incomplete AdamW state')
        for group in optimizer['param_groups']:
            if any(group[k] != CONFIG[k] for k in ('lr','weight_decay','eps')) or list(group['betas']) != CONFIG['betas']:
                raise ValueError('Saved original AdamW hyperparameters differ')
        for ident,parameter in zip(ids,parameters):
            saved = optimizer['state'][ident]
            if (int(saved['step']) != step or any(saved[k].shape != parameter.shape
                    or not bool(torch.isfinite(saved[k]).all()) for k in ('exp_avg','exp_avg_sq'))):
                raise ValueError('Nonfinite or misregistered AdamW moments')
        part = rows[index*1536:(index+1)*1536]
        if (dict(Counter(x['species'] for x in part)) != report['epochs'][index]['by_species']
                or any(not math.isclose(sum(x[k] for x in part)/1536,report['epochs'][index]['mean_losses'][k],
                    rel_tol=1e-12,abs_tol=1e-12) for k in ('raw_NLL','normalized_NLL','weighted_total'))):
            raise ValueError('Complete pass summary differs from log')
        previous = state['model']
    A.validate_checkpoint(state)
    C.write(VERIFY/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],source_fit_job=fit_job,
        configuration=CONFIG,protocol=PROTOCOL,final_checkpoint=report['final_checkpoint'],
        verified_additional_optimizer_steps=4608,verified_cumulative_optimizer_steps=9216,
        fresh_draws=4608,fresh_unique_windows=4219,fixed_fresh_draw_order_observations_and_importance_match=True,
        logged_objective_arithmetic_match=True,encoder_condition='JOINT24576',native_CDS_priors_used=False,
        final_checkpoint_strict_load=True,all_checkpoint_weights_and_optimizer_moments_finite=True,
        native_forward=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,interpretation='record_integrity_not_chain_accuracy'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('preflight','fit','verify'))
    parser.add_argument('--source-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'preflight':
        preflight()
    elif args.source_job is None:
        parser.error('--source-job required for fit/verify')
    elif args.mode == 'fit':
        fit(args.source_job)
    else:
        verify(args.source_job)
