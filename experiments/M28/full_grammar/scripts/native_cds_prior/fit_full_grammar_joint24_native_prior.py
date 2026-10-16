#!/usr/bin/env python3
"""The fourth fixed TRAIN cell: joint24576 plus both frozen native CDS priors."""
import argparse
from collections import Counter
import json
import os
import time
import torch
import fit_full_grammar_joint24 as J
import fit_full_grammar_native_cds_prior as N
import preflight_full_grammar_joint24_native_prior as V
import annotate_full_grammar_joint24_native_cds_prior as A

C,ROOT,Q = J.C,J.ROOT,J.Q
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-FIT'
BINDINGS = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-BINDINGS'
CONFIG,PROTOCOL = V.CONFIG,A.PROTOCOL


def require_preflight(value,status,job):
    if (status != 'COMPLETED' or not value['completed'] or value['job_id'] != job
            or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL
            or value['optimizer_steps_this_run'] != 2 or value['checkpoints_written']
            or value['formal_fit_started'] or not value['disposable_updates_only']
            or not value['pretrained_encoder_and_heads_frozen_no_grad']
            or value['DEV_used'] or value['sealed_data_used'] or len(value['records']) != 2
            or any(not r['zero_prior_original_joint_window_logits_exact']
                or not r['zero_prior_original_joint_loss_matches']
                or not r['strict_neural_weight_in_memory_roundtrip']
                or not r['updated_raw_adapter_potentials_exact'] for r in value['records'])):
        raise ValueError('Completed exact-zero and two-disposable-update combined native preflight required')


def fit(preflight_job,bindings_job):
    began = time.perf_counter()
    require_preflight(json.loads((V.OUT/'summary.json').read_text()),
        (V.OUT/'STATUS').read_text().strip(),preflight_job)
    checked = json.loads((BINDINGS/'summary.json').read_text())
    if ((BINDINGS/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['job_id'] != bindings_job or checked['configuration'] != CONFIG
            or checked['protocol'] != PROTOCOL):
        raise ValueError('Completed new fit, checkpoint and original-TRAIN evaluation bindings required')
    draws,windows,catalogs,plan = J.original_plan()
    C.write(FIT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=PROTOCOL,preflight_job=preflight_job,bindings_job=bindings_job,starting_fit_job='13434603'))
    state = C.starting_state()
    core,optimizer = N.load_model_optimizer(state)
    sequences,descriptors = C.P.read_train_DNA()
    full,tokenizer = A.N.load_native()
    adapter = A.N.NativeCDSHeadAdapter(core,full.score)
    step,epochs = 4608,[]
    torch.cuda.reset_peak_memory_stats()
    with (FIT/'training.jsonl').open('x',buffering=512*1024) as log:
        for pass_index in range(3):
            pass_started = time.perf_counter(); sums,species_counts = Counter(),Counter()
            for draw,expected in zip(draws,plan):
                started = time.perf_counter(); w = windows[draw['window_id']]
                chromosome,geometry,context,offset,targets = C.example(w,expected,sequences,catalogs)
                optimizer.zero_grad(set_to_none=True)
                joined = J.potentials(adapter,full.model,tokenizer,chromosome,geometry,w['strand'])
                raw,components = C.objective('STRUCTURED_CONTINUE',context,joined,offset,w['valid_bases'],targets)
                total,norm = C.update(core,optimizer,full,raw,expected['objective_multiplier'])
                step += 1
                row = dict(step=step,additional_step=step-4608,epoch=pass_index+4,
                    **{k:expected[k] for k in ('draw_ordinal','window_id','species','importance_weight',
                        'known_observation_denominator','fine_positions','coarse_positions','fully_unknown_positions')},
                    raw_NLL=float(raw.detach()),normalized_NLL=float(raw.detach())/expected['known_observation_denominator'],
                    weighted_total=total,gradient_norm_before_clip=norm,
                    prior_weight_l2=float(core.native_cds_role_prior.weight.detach().norm()),
                    seconds=time.perf_counter()-started)
                log.write(json.dumps(row,allow_nan=False)+'\n')
                sums.update({k:row[k] for k in ('raw_NLL','normalized_NLL','weighted_total')})
                species_counts[w['species']] += 1
                if step%64 == 0:
                    log.flush(); print(json.dumps(dict(step=step,final_step=9216,recent=row,
                        elapsed_seconds=time.perf_counter()-began)),flush=True)
                del raw,components,joined,targets
            checkpoint = FIT/('step_%06d.pt' % step)
            if checkpoint.exists(): raise FileExistsError(checkpoint)
            old = [p for name,p in core.named_parameters() if name != 'native_cds_role_prior.weight']
            if (any(int(optimizer.state[p]['step']) != step for p in old)
                    or int(optimizer.state[core.native_cds_role_prior.weight]['step']) != step-4608):
                raise ValueError('Every original/new AdamW group must preserve its own step origin')
            torch.save(dict(configuration=CONFIG,protocol=PROTOCOL,epoch=pass_index+4,step=step,
                encoder_revision=Q.F.REVISION,model=core.state_dict(),optimizer=optimizer.state_dict()),checkpoint)
            epochs.append(dict(epoch=pass_index+4,steps=1536,by_species=dict(species_counts),
                mean_losses={k:v/1536 for k,v in sums.items()},seconds=time.perf_counter()-pass_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step != 9216: raise ValueError('All4608 new original-draw updates required')
    if sum(p.stat().st_size for p in FIT.iterdir() if p.is_file()) > CONFIG['maximum_new_fit_bytes']:
        raise ValueError('Declared fit artifact ceiling exceeded')
    C.write(FIT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=PROTOCOL,preflight_job=preflight_job,bindings_job=bindings_job,starting_fit_job='13434603',
        cumulative_optimizer_steps=9216,optimizer_steps_this_run=4608,new_prior_optimizer_steps=4608,
        old_C0_optimizer_steps=9216,epochs=epochs,final_checkpoint=str(checkpoint.relative_to(ROOT)),
        trained_parameters=366219,original_draws_per_pass=len(draws),original_unique_windows=len(windows),
        original_TRAIN_DNA_sources=descriptors,original_draw_order_and_importance_preserved=True,
        original_targets_changed=False,context_or_inference_grammar_changed=False,
        encoder_only_change='JOINT24576_each_original_read_window',both_frozen_native_heads_used=True,
        pretrained_encoder_and_heads_frozen_no_grad=True,saved_source_weights_and_optimizer_unchanged=True,
        cached_features_used=False,DEV_used=False,sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),seconds=time.perf_counter()-began,
        interpretation='combined_original_TRAIN_cell_not_interaction_benefit_or_generalization'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight-job',required=True); parser.add_argument('--bindings-job',required=True)
    args = parser.parse_args(); torch.set_num_threads(2); fit(args.preflight_job,args.bindings_job)
