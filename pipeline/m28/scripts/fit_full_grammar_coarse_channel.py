#!/usr/bin/env python3
"""Reviewable original-cohort fit recipe; prepare never updates a model."""
import argparse
from collections import Counter
import json
import math
import os
import re
import time
import torch
import preflight_full_grammar_coarse_channel_glm as V
from src.m28.full_grammar_coarse_channel import coarse_channel_nll

G,N,Q,W,P,ROOT=V.G,V.N,V.Q,V.W,V.P,V.ROOT
READY=ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIT-PREPARATION'
FIT=ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT'
PROTOCOL='reports/M28-METHOD-RESTART/full_grammar_first_fit_contract.md'
CONFIG=dict(scope='original_two_species_six_TRAIN_chromosomes_1536_draws_1491_unique_windows',
    encoder_model=Q.F.MODEL,encoder_revision=Q.F.REVISION,encoder_frozen=True,
    encoder_blocks_bp=6144,encoder_token_bp=6,encoder_materialization='fresh_native_BF16_then_FP32_clone_each_draw',
    model='fresh_seed0_M28Core_C0_2048_128_366191_parameters_FP32',seed=0,
    draws_per_pass=1536,unique_windows=1491,passes=3,optimizer_steps=4608,
    training_order='original_paired_training_draws_in_original_order_each_pass',
    importance='original_importance_weight_uniform_within_stratum',
    context=G.CONFIG['view'],ownership=G.CONFIG['ownership'],reference_view_expansion=False,
    fine='unchanged_exact_role',coarse='normalized_three_class_correct_1_minus_epsilon_other_each_epsilon_over_two',
    unknown='likelihood_one',epsilon=0.1,epsilon_trainable=False,
    loss='importance_times_complete_gene_set_NLL_divided_by_max_1_fine_plus_coarse_core_rows',
    auxiliary_losses=False,max_exons=128,checkpoint_bp=256,role_weights=1,null=0,score_accumulation='FP64',
    interior_cuts='core_edges_as_in_actual_native_preflight',
    optimizer='AdamW',lr=3e-4,weight_decay=.01,betas=[.9,.999],eps=1e-8,clip_norm=1.,
    formal_checkpoint='last_step_004608_only_no_DEV_selection',automatic_resume=False,
    maximum_new_fit_bytes=64*1024*1024)


def source(root,job,config):
    value=json.loads((root/'summary.json').read_text())
    if ((root/'STATUS').read_text().strip()!='COMPLETED' or (root/'JOBID').read_text().strip()!=job
            or value['job_id']!=job or not value['completed'] or value['configuration']!=config
            or value['optimizer_steps_this_run']!=0 or value['DEV_used']):
        raise ValueError('Actual unchanged completed preparation source required: '+job)
    return value


def original_plan():
    support=source(N.S.OUT,'13429444',N.S.CONFIG)
    candidate=source(N.OUT,'13429625',N.CONFIG)
    native=source(V.OUT,'13429660',V.CONFIG)
    if (len(support['records'])!=1491 or not candidate['all_original_negative_candidates_feasible']
            or candidate['candidate_finite_events_and_gradients']!=24
            or native['actual_predetermined_negative_contexts']!=4
            or not native['new_model_weights_unchanged'] or not native['pretrained_encoder_frozen_no_grad']):
        raise ValueError('Original whole support and actual candidate/native evidence required')
    draws,index,windows,catalogs=Q.load_train_inputs()
    if [r['window_id'] for r in support['records']]!=list(windows):
        raise ValueError('Original unique TRAIN order changed')
    by_id={r['window_id']:r for r in support['records']}; plan=[]
    for ordinal,d in enumerate(draws):
        w=windows[d['window_id']]; r=by_id[w['id']]
        length=native['original_TRAIN_DNA_sources'][w['species']]['TRAIN_sequence_lengths'][w['seqid']]
        geometry=G.context_geometry(w['start'],w['end'],length); a,b=geometry['genomic_context']
        offset=w['start']-a if w['strand']=='+' else b-w['end']; n=w['valid_bases']
        fine,coarse,free=(r[k] for k in ('fine_positions','coarse_positions','fully_unknown_positions'))
        if (fine+coarse+free!=n or [a,b]!=r['genomic_context']
                or [offset,offset+n]!=r['core_oriented_interval']):
            raise ValueError('Original source observation or fixed view changed')
        weight=float(d['importance_weight_uniform_within_stratum']); denominator=max(1,fine+coarse)
        plan.append(dict(draw_ordinal=ordinal,window_id=w['id'],species=d['species'],stratum=d['stratum'],
            importance_weight=weight,fine_positions=fine,coarse_positions=coarse,fully_unknown_positions=free,
            known_observation_denominator=denominator,objective_multiplier=weight/denominator,
            genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
            original_hard_event_support_feasible=r['support_feasible']))
    return draws,index,windows,catalogs,plan,native


def prepare():
    began=time.perf_counter()
    N.S.write_json(READY/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,protocol=PROTOCOL))
    draws,_,windows,_,plan,native=original_plan()
    with (READY/'training_plan.jsonl').open('x') as handle:
        for row in plan: handle.write(json.dumps(row,allow_nan=False)+'\n')
    by_species={s:dict(draws=sum(d['species']==s for d in draws),unique_windows=sum(w['species']==s for w in windows.values()),
        original_importance_sum=sum(r['importance_weight'] for r in plan if r['species']==s),
        fine_exposures=sum(r['fine_positions'] for r in plan if r['species']==s),
        coarse_exposures=sum(r['coarse_positions'] for r in plan if r['species']==s),
        unknown_exposures=sum(r['fully_unknown_positions'] for r in plan if r['species']==s)) for s in Q.ALLOW}
    warm=max(r['seconds'] for r in native['records'][1:])
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,protocol=PROTOCOL,
        source_jobs=['13429444','13429625','13429660'],original_draws=len(draws),original_unique_windows=len(windows),
        proposed_optimizer_steps=4608,proposed_training_recipe_only=True,fit_submitted=False,
        original_draw_order_and_importance_preserved=True,original_R4_labels_and_fixed_views_preserved=True,
        by_species=by_species,known_denominator_range=[min(r['known_observation_denominator'] for r in plan),max(r['known_observation_denominator'] for r in plan)],
        objective_multiplier_range=[min(r['objective_multiplier'] for r in plan),max(r['objective_multiplier'] for r in plan)],
        no_known_observation_draws=sum(r['fine_positions']+r['coarse_positions']==0 for r in plan),
        hard_negative_draws=sum(not r['original_hard_event_support_feasible'] for r in plan),
        unique_hard_negative_windows=len({r['window_id'] for r in plan if not r['original_hard_event_support_feasible']}),
        observed_steady_native_max_window_seconds=warm,steady_projection_seconds=warm*4608,
        projection_is_four_context_estimate_not_full_cohort_runtime_measurement=True,
        intended_fit_partition='private-teodoro-gpu',intended_fit_GPU='one_3090_at_least_24GiB',
        intended_fit_CPU=2,intended_fit_RAM_GiB=16,intended_fit_wall_hours=30,
        actual_raw_FASTA_read=False,feature_tensors_read=False,models_or_weights_read=False,
        optimizer_constructed=False,optimizer_steps_this_run=0,checkpoints_written=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,A44_execution_chain_reopened=False,
        first_fit_authority_not_inferred_from_this_record=True,seconds=time.perf_counter()-began,
        interpretation='audited_proposed_original_cohort_recipe_not_training_or_model_performance')
    N.S.write_json(READY/'summary.json',summary); print(json.dumps(summary,allow_nan=False),flush=True)


def fit():
    # The caller must already have scientific/first-fit authority. Artifact
    # equality below verifies inputs, not user consent; no approval-file scheme.
    began=time.perf_counter(); prepared=json.loads((READY/'summary.json').read_text())
    if ((READY/'STATUS').read_text().strip()!='COMPLETED' or not prepared['completed']
            or prepared['configuration']!=CONFIG or prepared['protocol']!=PROTOCOL
            or prepared['optimizer_steps_this_run']!=0):
        raise ValueError('Matching complete proposed recipe preparation required')
    draws,index,windows,catalogs,plan,_=original_plan()
    stored=[json.loads(x) for x in (READY/'training_plan.jsonl').read_text().splitlines()]
    if stored!=plan: raise ValueError('Original draw plan or loss multipliers changed')
    N.S.write_json(FIT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,protocol=PROTOCOL,
        preparation_job=prepared['job_id'],source_jobs=['13429444','13429625','13429660']))
    sequences,descriptors=P.read_train_DNA(); backbone,tokenizer=G.load_native()
    torch.manual_seed(0); model=G.M28Core('C0',2048,128).to('cuda').train()
    if sum(p.numel() for p in model.parameters())!=366191: raise ValueError('Fresh C0 capacity changed')
    before={k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
    parameters=list(model.parameters())
    optimizer=torch.optim.AdamW(parameters,lr=CONFIG['lr'],weight_decay=CONFIG['weight_decay'],
        betas=tuple(CONFIG['betas']),eps=CONFIG['eps'])
    checked=set(); opened=set(); step=0; epochs=[]; torch.cuda.reset_peak_memory_stats()
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(CONFIG['passes']):
            epoch_started=time.perf_counter(); sums=Counter(); species_counts=Counter()
            for d,expected in zip(draws,plan):
                t=time.perf_counter(); w=windows[d['window_id']]; n=w['valid_bases']
                chromosome=sequences[w['species']][w['seqid']]
                geometry=G.context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
                raw=chromosome[w['start']:w['end']]; core=raw if w['strand']=='+' else Q.M.rc(raw)
                if w['id'] not in checked:
                    cached=W.load_cache(w,index,opened)
                    if cached['sequence']!=core: raise ValueError('Actual original core DNA changed')
                    del cached; checked.add(w['id'])
                targets=P.P.build_targets(w,catalogs[w['species']],core)
                fine=torch.as_tensor(targets['seg'])[:n]; region=torch.as_tensor(targets['region'])[:n]
                counts=(int((fine>=0).sum()),int(((fine<0)&(region>=0)).sum()),int(((fine<0)&(region<0)).sum()))
                offset=w['start']-a if w['strand']=='+' else b-w['end']
                if (counts!=tuple(expected[k] for k in ('fine_positions','coarse_positions','fully_unknown_positions'))
                        or [a,b]!=expected['genomic_context'] or [offset,offset+n]!=expected['core_oriented_interval']):
                    raise ValueError('Actual original observations or proposed fixed view changed')
                context=re.sub('[^ACGT]','N',chromosome[a:b])
                if w['strand']=='-': context=Q.M.rc(context)
                optimizer.zero_grad(set_to_none=True); tensors=[]
                for left,right in geometry['genomic_windows']:
                    dna_string=re.sub('[^ACGT]','N',chromosome[left:right])
                    if w['strand']=='-': dna_string=Q.M.rc(dna_string)
                    features=G.encode_window(backbone,tokenizer,dna_string,'cuda').clone()
                    dna=G.one_hot(dna_string,Q.M.WIDTH)[None].to('cuda')
                    tensors.append(model(features,dna,[right-left])['segmentation_logits'][0,:right-left])
                joined=G.owned_context_logits(tensors,geometry,w['strand']).to('cpu')
                cuts=tuple(x for x in (offset,offset+n) if 0<x<len(context))
                result=coarse_channel_nll(context,joined,128,offset,offset+n,targets,CONFIG['epsilon'],cuts=cuts,checkpoint_bp=256)
                raw_NLL=result['loss']; total=raw_NLL*expected['objective_multiplier']
                if not bool(torch.isfinite(total)) or float(raw_NLL.detach()) < -1e-8:
                    raise ValueError('Nonfinite or negative structured observation likelihood')
                total.backward()
                if any(p.grad is None for p in parameters): raise ValueError('Missing trainable C0 gradient')
                norm=torch.nn.utils.clip_grad_norm_(parameters,CONFIG['clip_norm'],error_if_nonfinite=True)
                if any(p.requires_grad or p.grad is not None for p in backbone.parameters()):
                    raise ValueError('Released GLM must remain frozen')
                optimizer.step(); step+=1; torch.cuda.synchronize()
                record=dict(step=step,epoch=epoch+1,draw_ordinal=expected['draw_ordinal'],window_id=w['id'],species=w['species'],
                    importance_weight=expected['importance_weight'],known_observation_denominator=expected['known_observation_denominator'],
                    fine_positions=counts[0],coarse_positions=counts[1],fully_unknown_positions=counts[2],
                    raw_NLL=float(raw_NLL.detach()),normalized_NLL=float(raw_NLL.detach())/expected['known_observation_denominator'],
                    weighted_total=float(total.detach()),gradient_norm_before_clip=float(norm),seconds=time.perf_counter()-t)
                log.write(json.dumps(record,allow_nan=False)+'\n'); sums.update({k:record[k] for k in ('raw_NLL','normalized_NLL','weighted_total')})
                species_counts[w['species']]+=1
                if step%64==0:
                    log.flush(); print(json.dumps(dict(step=step,total_steps=4608,recent=record,elapsed_seconds=time.perf_counter()-began)),flush=True)
                del result,raw_NLL,total,joined,tensors,features,dna,targets,fine,region
            checkpoint=FIT/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists(): raise FileExistsError(checkpoint)
            torch.save(dict(configuration=CONFIG,protocol=PROTOCOL,epoch=epoch+1,step=step,
                encoder_revision=Q.F.REVISION,model=model.state_dict(),optimizer=optimizer.state_dict()),checkpoint)
            epochs.append(dict(epoch=epoch+1,steps=1536,by_species=dict(species_counts),
                mean_losses={k:v/1536 for k,v in sums.items()},seconds=time.perf_counter()-epoch_started,
                checkpoint=str(checkpoint.relative_to(ROOT))))
    if step!=4608 or len(checked)!=1491 or opened!={index[ident]['artifact'] for ident in windows}:
        raise ValueError('Complete original three-pass exposure and all original core DNA checks required')
    if all(torch.equal(v.cpu(),before[k]) for k,v in model.state_dict().items()):
        raise ValueError('Optimizer updates did not change the fresh model')
    if any(p.requires_grad or p.grad is not None for p in backbone.parameters()): raise ValueError('Released GLM changed freeze state')
    size=sum(p.stat().st_size for p in FIT.iterdir() if p.is_file())
    if size>CONFIG['maximum_new_fit_bytes']: raise ValueError('Proposed fit artifact ceiling exceeded')
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,protocol=PROTOCOL,
        preparation_job=prepared['job_id'],optimizer_steps=step,optimizer_steps_this_run=step,epochs=epochs,
        final_checkpoint=str(checkpoint.relative_to(ROOT)),trained_parameters=366191,new_model_weights_changed=True,
        pretrained_encoder_frozen_no_grad=True,original_R4_policy_changed=False,observation_rows_deleted=False,
        candidate_observation_model_changed=True,epsilon_learned=False,original_TRAIN_DNA_sources=descriptors,
        original_draw_order_and_importance_preserved=True,original_unique_windows=1491,original_draws_per_pass=1536,
        original_core_cache_DNA_checks=len(checked),cached_features_used=False,old_cached_labels_used=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,A44_execution_chain_reopened=False,
        seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),fit_file_bytes_before_summary=size,
        interpretation='fixed_original_TRAIN_discovery_fit_not_performance_generalization_or_publication_success')
    N.S.write_json(FIT/'summary.json',summary); print(json.dumps(summary,allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('mode',choices=('prepare','fit'))
    args=parser.parse_args(); torch.set_num_threads(2)
    prepare() if args.mode=='prepare' else fit()
