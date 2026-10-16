#!/usr/bin/env python3
"""Fixed importance-corrected TRAIN continuation and final9216 artifact verification."""
import argparse
from collections import Counter
import json
import math
import os
import time
import torch
import preflight_full_grammar_complete_core_draws as N

C,P,A = N.C,N.P,N.A
ROOT,CONFIG = N.ROOT,N.CONFIG
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-VERIFICATION'
LOG_FIELDS = ('additional_step','pass_index','draw_ordinal','window_id','species',
              'stratum_size','stratum_probability','importance_weight','known_observation_denominator',
              'fine_positions','coarse_positions','fully_unknown_positions')


def fit(preflight_job):
    began = time.perf_counter()
    ready = json.loads((N.PREFLIGHT/'summary.json').read_text())
    if ((N.PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != preflight_job or ready['configuration'] != CONFIG
            or ready['optimizer_steps_this_run'] != 2 or not ready['actual_native_logits_raw_adapter_exact']
            or ready['checkpoints_written'] or ready['formal_training_started']
            or ready['DEV_used'] or ready['sealed_data_used']):
        raise ValueError('Completed two-update weighted native preflight required')
    plan,windows,catalogs = N.plan_and_population()
    C.write(FIT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=A.PROTOCOL,preflight_job=preflight_job,starting_fit_job='13434603',support_job=A.IDENTITY['support_job']))
    state = C.starting_state()
    model,optimizer = C.load_model_optimizer(state)
    sequences,descriptors = C.P.read_train_DNA()
    backbone,tokenizer = C.G.load_native()
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
                joined = C.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
                raw,_ = C.objective('STRUCTURED_CONTINUE',context,joined,offset,w['valid_bases'],targets)
                total,norm = C.update(model,optimizer,backbone,raw,expected['objective_multiplier'])
                step = 4608+expected['additional_step']
                row = dict(step=step,epoch=pass_index+4,**{k:expected[k] for k in LOG_FIELDS},
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
            torch.save(dict(configuration=CONFIG,protocol=A.PROTOCOL,epoch=pass_index+4,step=step,
                encoder_revision=C.Q.F.REVISION,model=model.state_dict(),optimizer=optimizer.state_dict()),checkpoint)
            epochs.append(dict(epoch=pass_index+4,steps=1536,by_species=dict(species_counts),
                mean_losses={k:v/1536 for k,v in sums.items()},seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 9216 or any(int(s['step']) != 9216 for s in optimizer.state.values()):
        raise ValueError('Complete fixed 4608 additional weighted updates required')
    if sum(p.stat().st_size for p in FIT.iterdir() if p.is_file()) > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared fit artifact ceiling exceeded')
    C.write(FIT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=A.PROTOCOL,preflight_job=preflight_job,starting_fit_job='13434603',
        preparation_job=A.IDENTITY['preparation_job'],support_job=A.IDENTITY['support_job'],
        cumulative_optimizer_steps=9216,optimizer_steps_this_run=4608,trained_parameters=366191,
        fixed_draws=4608,selected_unique_windows=4107,draws_per_pass=1536,epochs=epochs,
        final_checkpoint=str(checkpoint.relative_to(ROOT)),original_TRAIN_DNA_sources=descriptors,
        fixed_draw_order_proposal_probabilities_and_importance_preserved=True,all_fixed_draws_retained=True,
        original_targets_changed=False,context_or_inference_grammar_changed=False,
        pretrained_encoder_frozen_no_grad=True,saved_source_weights_and_optimizer_unchanged=True,
        cached_features_used=False,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began,
        interpretation='TRAIN_importance_corrected_complete_core_proposal_not_pure_coverage_causality_or_generalization'))


def verify(fit_job):
    began = time.perf_counter()
    summary = json.loads((FIT/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed']
            or summary['job_id'] != fit_job or summary['configuration'] != CONFIG
            or summary['optimizer_steps_this_run'] != 4608 or summary['cumulative_optimizer_steps'] != 9216
            or summary['DEV_used'] or summary['sealed_data_used']):
        raise ValueError('Completed matching final complete-core-proposal fit required')
    plan,_,_ = N.plan_and_population()
    rows = list(P.P.rows(FIT/'training.jsonl'))
    if len(rows) != 4608:
        raise ValueError('All 4608 actual weighted update records required')
    for i,(row,expected) in enumerate(zip(rows,plan)):
        if (row['step'] != 4609+i or row['epoch'] != 4+i//1536
                or any(row[k] != expected[k] for k in LOG_FIELDS)):
            raise ValueError('Actual draw, probability, importance or update order differs')
        values = [row[k] for k in ('raw_NLL','normalized_NLL','weighted_total','gradient_norm_before_clip','seconds')]
        raw = row['raw_NLL']
        if (not all(math.isfinite(x) and x >= -1e-8 for x in values)
                or not math.isclose(row['normalized_NLL'],raw/expected['known_observation_denominator'],rel_tol=1e-12,abs_tol=1e-12)
                or not math.isclose(row['weighted_total'],raw*expected['objective_multiplier'],rel_tol=1e-12,abs_tol=1e-12)):
            raise ValueError('Logged importance-corrected objective arithmetic differs')
        if expected['fine_positions']+expected['coarse_positions'] == 0:
            if abs(raw) > 1e-8 or row['gradient_norm_before_clip'] > 1e-6:
                raise ValueError('All-unknown core gained an observation gradient')
    previous = C.starting_state()['model']
    for index in range(3):
        step = 6144+index*1536
        state = torch.load(FIT/('step_%06d.pt' % step),map_location='cpu',weights_only=True)
        if (state['configuration'] != CONFIG or state['protocol'] != A.PROTOCOL or state['step'] != step
                or state['epoch'] != index+4 or state['encoder_revision'] != C.Q.F.REVISION):
            raise ValueError('Checkpoint identity differs')
        model = C.G.M28Core('C0',2048,128)
        model.load_state_dict(state['model'],strict=True)
        if (not all(v.dtype == torch.float32 and bool(torch.isfinite(v).all()) for v in state['model'].values())
                or all(torch.equal(v,previous[k]) for k,v in state['model'].items())):
            raise ValueError('Nonfinite or unchanged checkpoint weights')
        optimizer = state['optimizer']
        ids = [p for group in optimizer['param_groups'] for p in group['params']]
        parameters = list(model.parameters())
        if len(ids) != len(parameters) or len(ids) != len(optimizer['state']):
            raise ValueError('Incomplete AdamW state')
        for ident,parameter in zip(ids,parameters):
            saved = optimizer['state'][ident]
            if (int(saved['step']) != step or any(saved[k].shape != parameter.shape
                    or not bool(torch.isfinite(saved[k]).all()) for k in ('exp_avg','exp_avg_sq'))):
                raise ValueError('Nonfinite or misregistered AdamW moments')
        part = rows[index*1536:(index+1)*1536]
        if (dict(Counter(x['species'] for x in part)) != summary['epochs'][index]['by_species']
                or any(not math.isclose(sum(x[k] for x in part)/1536,summary['epochs'][index]['mean_losses'][k],
                    rel_tol=1e-12,abs_tol=1e-12) for k in ('raw_NLL','normalized_NLL','weighted_total'))):
            raise ValueError('Complete pass summary differs from actual update log')
        previous = state['model']
    A.validate_checkpoint(state)
    C.write(VERIFY/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],source_fit_job=fit_job,
        configuration=CONFIG,protocol=A.PROTOCOL,final_checkpoint=summary['final_checkpoint'],
        verified_additional_optimizer_steps=4608,verified_cumulative_optimizer_steps=9216,
        fixed_draws=4608,selected_unique_windows=4107,
        fixed_draw_order_observations_probabilities_and_importance_match=True,logged_objective_arithmetic_match=True,
        final_checkpoint_strict_load=True,all_checkpoint_weights_and_optimizer_moments_finite=True,
        native_forward=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,interpretation='record_integrity_not_chain_accuracy'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('fit','verify'))
    parser.add_argument('--source-job',required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'fit':
        fit(args.source_job)
    else:
        verify(args.source_job)
