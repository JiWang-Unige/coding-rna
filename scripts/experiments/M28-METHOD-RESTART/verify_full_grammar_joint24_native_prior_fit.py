#!/usr/bin/env python3
"""Reconcile combined-cell original exposures and both AdamW step origins on CPU."""
import argparse
from collections import Counter
import json
import math
import os
import time
import torch
import fit_full_grammar_joint24_native_prior as F
import annotate_full_grammar_joint24_native_cds_prior as A

OUT = F.ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-VERIFICATION'


def require(value,message):
    if not value: raise ValueError(message)


def verify_rows(rows,plan):
    require(len(plan) == 1536 and len(rows) == 4608,'Complete original1536 draws and4608 new updates required')
    for i,row in enumerate(rows):
        expected = plan[i%1536]
        require((row['step'],row['additional_step'],row['epoch']) == (4609+i,1+i,4+i//1536),'Step discontinuity')
        for key in ('draw_ordinal','window_id','species','importance_weight','known_observation_denominator',
                    'fine_positions','coarse_positions','fully_unknown_positions'):
            require(row[key] == expected[key],'Original exposure differs: '+key)
        require(all(math.isfinite(row[k]) and row[k] >= -1e-8 for k in
            ('raw_NLL','normalized_NLL','weighted_total','gradient_norm_before_clip','prior_weight_l2','seconds')),
            'Invalid logged training numeric result')
        require(math.isclose(row['normalized_NLL'],row['raw_NLL']/expected['known_observation_denominator'],
            rel_tol=1e-12,abs_tol=1e-12) and math.isclose(row['weighted_total'],
            row['raw_NLL']*expected['objective_multiplier'],rel_tol=1e-12,abs_tol=1e-12),
            'Original likelihood normalization or importance differs')
        if expected['fine_positions']+expected['coarse_positions'] == 0:
            require(abs(row['raw_NLL']) <= 1e-8 and row['gradient_norm_before_clip'] <= 1e-6,
                'All-unknown core gained an observation gradient')


def main(fit_job):
    began = time.perf_counter()
    summary = json.loads((F.FIT/'summary.json').read_text())
    started = json.loads((F.FIT/'started.json').read_text())
    require((F.FIT/'STATUS').read_text().strip() == 'COMPLETED' and summary['completed']
        and summary['job_id'] == started['job_id'] == fit_job
        and summary['configuration'] == started['configuration'] == F.CONFIG
        and summary['protocol'] == F.PROTOCOL and summary['preflight_job'] == started['preflight_job'] == '13454690'
        and summary['bindings_job'] == started['bindings_job']
        and summary['optimizer_steps_this_run'] == summary['new_prior_optimizer_steps'] == 4608
        and summary['cumulative_optimizer_steps'] == summary['old_C0_optimizer_steps'] == 9216
        and not summary['DEV_used'] and not summary['sealed_data_used'],'Complete declared combined fit required')
    plan = [json.loads(x) for x in (F.J.F.READY/'training_plan.jsonl').read_text().splitlines()]
    rows = [json.loads(x) for x in (F.FIT/'training.jsonl').read_text().splitlines()]
    require(len({x['window_id'] for x in plan}) == 1491,'Original1491 unique exposure windows required')
    verify_rows(rows,plan)
    previous = F.C.starting_state()['model']; checkpoints = []
    for index in range(3):
        step = 6144+index*1536; checkpoint = F.FIT/('step_%06d.pt' % step)
        state = torch.load(checkpoint,map_location='cpu',weights_only=True)
        require(state['configuration'] == F.CONFIG and state['protocol'] == F.PROTOCOL
            and state['step'] == step and state['epoch'] == index+4 and state['encoder_revision'] == F.Q.F.REVISION,
            'Stored combined checkpoint identity differs')
        model = A.NativeCDSRolePrior(); model.load_state_dict(state['model'],strict=True)
        require(all(v.dtype == torch.float32 and bool(torch.isfinite(v).all()) for v in state['model'].values()),
            'Nonfinite or non-FP32 checkpoint weights')
        require(any(not torch.equal(state['model'][k],v) for k,v in previous.items()),'Pass did not change any weights')
        optimizer = state['optimizer']; groups = optimizer['param_groups']
        parameters = list(model.parameters()); ids = [p for group in groups for p in group['params']]
        require(len(groups) == 2 and len(groups[1]['params']) == 1
            and len(ids) == len(parameters) == len(optimizer['state']),'Separate complete old/new AdamW groups required')
        for group in groups:
            require(all(group[k] == F.CONFIG[k] for k in ('lr','weight_decay','eps'))
                and list(group['betas']) == F.CONFIG['betas'],'AdamW hyperparameters changed')
        for i,(ident,parameter) in enumerate(zip(ids,parameters)):
            saved = optimizer['state'][ident]
            require(int(saved['step']) == (step-4608 if i == len(ids)-1 else step),'Old/new optimizer step differs')
            require(all(saved[k].shape == parameter.shape and bool(torch.isfinite(saved[k]).all())
                for k in ('exp_avg','exp_avg_sq')),'Invalid AdamW moment tensor')
        part = rows[index*1536:(index+1)*1536]
        require(dict(Counter(x['species'] for x in part)) == summary['epochs'][index]['by_species'],'Species exposures differ')
        require(all(math.isclose(sum(x[k] for x in part)/1536,summary['epochs'][index]['mean_losses'][k],
            rel_tol=1e-12,abs_tol=1e-12) for k in ('raw_NLL','normalized_NLL','weighted_total')),'Summary loss differs')
        checkpoints.append(dict(step=step,epoch=index+4,file=str(checkpoint.relative_to(F.ROOT))))
        previous = state['model']
    A.validate_checkpoint(state)
    require(summary['final_checkpoint'] == str(checkpoint.relative_to(F.ROOT)),'Final checkpoint path differs')
    F.C.write(OUT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],source_fit_job=fit_job,
        configuration=F.CONFIG,protocol=F.PROTOCOL,final_checkpoint=summary['final_checkpoint'],checkpoints=checkpoints,
        verified_additional_optimizer_steps=4608,verified_old_C0_optimizer_steps=9216,
        verified_new_prior_optimizer_steps=4608,original_draw_order_observations_and_importance_match=True,
        logged_objective_arithmetic_match=True,final_checkpoint_strict_load=True,
        all_checkpoint_weights_and_optimizer_moments_finite=True,
        native_prior_weight_nonzero=bool(state['model']['native_cds_role_prior.weight'].count_nonzero()),
        optimizer_steps_this_run=0,native_forward=False,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,interpretation='combined_training_record_integrity_not_chain_accuracy'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--fit-job',required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.fit_job)
