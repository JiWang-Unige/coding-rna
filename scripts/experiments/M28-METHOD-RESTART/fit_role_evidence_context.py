#!/usr/bin/env python3
"""One fixed A40 fit; frozen A3 role inputs, fresh role-conditioned student."""
import argparse
from collections import Counter
import json
import time
import torch
import preflight_role_evidence_context as A
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss

ROOT, CONFIG, PROTOCOL = A.ROOT, A.CONFIG, A.PROTOCOL
OUT = ROOT/'outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-FIT'
PREFLIGHT_JOB = '13418097'


def ready():
    A.A.require_jobs([PREFLIGHT_JOB])
    report=A.M.R.completed(A.OUT,PREFLIGHT_JOB)
    if (report['configuration']!=CONFIG or report['protocol']!=PROTOCOL or report['windows']!=1491 or
            report['original_draws']!=1536 or report['natural_candidates']!=16321 or len(report['records'])!=8 or
            report['eligible_reference_pairs']!={'eligible':2916,'five_or_more':1086} or
            report['optimizer_steps_this_run']!=0 or report['reference_injected_training_candidates']!=0 or
            report['frozen_supervised_A3_checkpoint_read'] is not True or
            any(report[k] is not True for k in ('original_natural_identity_labels_and_stronger_baseline_reconciled',
                'fresh_weights_unchanged','frozen_local_weights_unchanged','frozen_local_gradients_absent',
                'original_A31_local_source_scores_reconciled','actual_role_column_gradients_verified',
                'both_species_known_junction_gradient_coverage')) or
            any(report[k] is not False for k in ('student_supervised_checkpoint_read','checkpoints_written',
                'generator_decoding','generator_updated','new_candidate_search','new_backbone_forward',
                'new_R4_labelling','DEV_used','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual fixed zero-update role-conditioned preparation required')


def main():
    ready()
    draws,index,windows,records,_,_,_=A.A.population()
    student=A.new_model('cuda').train(); source=A.local_source()
    params=list(student.parameters())
    if sum(p.numel() for p in params)!=CONFIG['trainable_parameters']:
        raise ValueError('Fixed fresh student parameter count differs')
    source_before={k:v.detach().clone() for k,v in source.state_dict().items()}
    optimizer=torch.optim.AdamW(params,lr=CONFIG['lr'],weight_decay=CONFIG['weight_decay'],
        betas=tuple(CONFIG['betas']),eps=CONFIG['eps'])
    began=time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step=0; epochs=[]
    with (OUT/'training.jsonl').open('x') as handle:
        for epoch in range(CONFIG['proposed_passes']):
            started=time.perf_counter(); counts=Counter(); total=0.
            for ordinal,draw in enumerate(draws):
                ident=draw['window_id']; w=windows[ident]; rec=records[ident]
                cached=torch.load(ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
                features,dna=A.R.C.R.Q.window_inputs(cached,w,'cuda')
                if features.requires_grad or dna.requires_grad:
                    raise ValueError('Original encoder features and DNA must stay detached')
                logits=A.local_logits(source,features,dna,w['valid_bases'])
                optimizer.zero_grad(set_to_none=True)
                chains,labels=[r['chain'] for r in rec['candidates']],[r['label'] for r in rec['candidates']]
                score=student(features,dna,[w['valid_bases']],chains,logits)
                if any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
                    raise ValueError('Finite fixed FP32 role-conditioned components required')
                result=null_anchored_chain_set_loss(score,chains,labels)
                weighted=result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(weighted)):
                    raise ValueError('Nonfinite original fixed composite objective')
                weighted.backward()
                if logits.requires_grad or any(p.requires_grad or p.grad is not None for p in source.parameters()):
                    raise ValueError('Frozen local source received gradient')
                norm=torch.nn.utils.clip_grad_norm_(params,CONFIG['clip_norm'],error_if_nonfinite=True)
                optimizer.step(); step+=1
                row=dict(step=step,epoch=epoch+1,draw_ordinal=ordinal,window_id=ident,species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'],
                    weighted_loss=float(weighted.detach()),
                    unweighted_mean_known_BCE=float(result['mean_known_BCE'].detach()),
                    unweighted_mean_known_set_NLL=float(result['mean_known_set_NLL'].detach()),
                    unweighted_fixed_composite_loss=float(result['loss'].detach()),
                    log_partition=float(result['log_partition'].detach()),gold_gain=float(result['gold_gain'].detach()),
                    gradient_norm_before_clip=float(norm),frozen_local_source_gradients_absent=True,**rec['counts'])
                handle.write(json.dumps(row,allow_nan=False)+'\n'); total+=row['weighted_loss']; counts.update(rec['counts'])
                if step%64==0:
                    handle.flush(); print(json.dumps(dict(step=step,total=4608,seconds=time.perf_counter()-began)),flush=True)
                del cached,features,dna,logits,score,result,weighted
            state=student.state_dict()
            if (any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in state.values()) or
                    any(not torch.equal(v,source_before[k]) for k,v in source.state_dict().items())):
                raise ValueError('Finite student state and unchanged frozen A3 source required')
            checkpoint=OUT/('step_'+str(step).zfill(6)+'.pt')
            with checkpoint.open('xb') as file:
                torch.save(dict(arm='A40',configuration=CONFIG,protocol=PROTOCOL,step=step,
                    revision=A.R.C.I.F.REVISION,model=state,optimizer=optimizer.state_dict(),
                    frozen_local_source_model=source.state_dict(),frozen_local_source_checkpoint=CONFIG['frozen_local_source']),file)
            epochs.append(dict(epoch=epoch+1,steps=1536,counts=dict(counts),mean_weighted_loss=total/1536,
                seconds=time.perf_counter()-started,checkpoint=str(checkpoint.relative_to(ROOT))))
    if step!=4608 or step!=CONFIG['proposed_optimizer_steps']:
        raise ValueError('Incomplete fixed A40 exposure')
    exposure={name:int(optimizer.state[p]['step']) if p in optimizer.state else 0
              for name,p in student.named_parameters()}
    A.R.C.R.write_summary(OUT,dict(completed=True,arm='A40',configuration=CONFIG,protocol=PROTOCOL,
        preflight_job=PREFLIGHT_JOB,feature_revision=A.R.C.I.F.REVISION,feature_spec=A.R.C.R.FEATURE_SPEC,
        optimizer_steps=step,optimizer_steps_per_parameter=exposure,epochs=epochs,
        trainable_parameters=CONFIG['trainable_parameters'],
        final_checkpoint=str((OUT/'step_004608.pt').relative_to(ROOT)),
        initialization=CONFIG['initialization'],main_evaluation='last_step_004608_only_same_natural_TRAIN_pool_stronger_A31_A3',
        student_supervised_checkpoint_read=False,frozen_supervised_A3_checkpoint_read=True,
        frozen_local_source_checkpoint=CONFIG['frozen_local_source'],frozen_local_source_weights_unchanged=True,
        frozen_local_source_gradients_absent=True,frozen_local_source_optimizer_steps=0,
        original_TRAIN_catalog_metadata_read=True,pretrained_encoder_frozen=True,
        reference_injected_training_candidates=0,generator_decoding=False,generator_updated=False,
        new_candidate_search=False,new_backbone_forward=False,new_R4_labelling=False,
        DEV_used_for_selection=False,test_setaria_used=False,A33_outputs_opened=False,
        wall_seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One fixed fresh role-conditioned TRAIN student using frozen supervised A3 local inputs; '
            'not complete readout, matched-capacity causal evidence, gene F1 or independent generalization.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight-job',required=True); args=parser.parse_args()
    if args.preflight_job!=PREFLIGHT_JOB: parser.error('Only actual recorded A40 preparation allowed')
    torch.set_num_threads(2); main()
