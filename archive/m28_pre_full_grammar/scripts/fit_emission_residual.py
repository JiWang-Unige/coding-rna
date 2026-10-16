#!/usr/bin/env python3
"""One fixed A41 fit after actual zero-start readiness and closed A40 results."""
from collections import Counter
import json
import time
import torch
import preflight_emission_residual as A
import diagnose_role_evidence_attrition as H
import run_local_path_caller as R
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss

ROOT, PROTOCOL = A.ROOT, A.PROTOCOL
OUT = ROOT/'outputs/M28-A41-EMISSION-RESIDUAL-FIT'
PREFLIGHT_JOB = '13420270'
CONFIG = {**{k:v for k,v in A.CONFIG.items() if k!='optimizer_steps_this_run'},
    'proposed_passes':3,'proposed_optimizer_steps':4608,'optimizer':'AdamW',
    'lr':3e-4,'weight_decay':.01,'betas':[.9,.999],'eps':1e-8,'clip_norm':1.,
    'proposed_main_checkpoint':'last_step_004608_only'}


def validate_preflight(report):
    if (report['completed'] is not True or report['configuration']!=A.CONFIG or report['protocol']!=PROTOCOL or
            report['windows']!=1491 or report['original_draws']!=1536 or report['natural_candidates']!=16321 or
            report['eligible_reference_pairs']!={'eligible':2916,'five_or_more':1086} or len(report['records'])!=8 or
            report['optimizer_steps_this_run']!=0 or report['checkpoints_written'] is not False or
            any(report[k] is not True for k in ('original_natural_identity_labels_and_stronger_baseline_reconciled',
                'initial_A31_gain_and_selection_exact','residual_weights_unchanged','frozen_source_weights_unchanged',
                'frozen_source_gradients_absent','actual_both_species_residual_gradients_nonzero')) or
            any(report[k] is not False for k in ('generator_decoding','generator_updated','new_candidate_search',
                'new_backbone_forward','new_R4_labelling','A40_outputs_opened','DEV_used','test_setaria_used',
                'A33_outputs_opened'))):
        raise ValueError('Actual unchanged A41 zero-update preparation required')


def ready():
    A.A.require_jobs([PREFLIGHT_JOB,H.I.FIT_JOB,H.I.READOUT_JOB,H.I.CHECK_JOB])
    validate_preflight(R.completed(A.OUT,PREFLIGHT_JOB))
    # This reconciles a completed prior arm; it neither grants authority nor tunes A41.
    return H.validate_completed(H.R.completed(H.D.F.OUT,H.I.FIT_JOB),
        H.R.completed(H.D.OUT,H.I.READOUT_JOB),H.R.completed(H.D.CHECK,H.I.CHECK_JOB))


def source_model():
    source=R.C.load_generator().to('cuda').eval().requires_grad_(False)
    if sum(p.numel() for p in source.parameters())!=CONFIG['frozen_source_parameters']:
        raise ValueError('Original frozen A3 structure differs')
    return source


@torch.no_grad()
def local_evidence(source,features,dna,length):
    hidden,_=source.shared(features,dna,[length])
    hidden=hidden[0]
    return hidden,source.segmentation(hidden)


def main():
    previous_gates=ready()
    draws,index,windows,records,_,_,_=A.A.population()
    model=A.new_model('cuda').train(); source=source_model(); params=list(model.parameters())
    before={k:v.detach().clone() for k,v in source.state_dict().items()}
    optimizer=torch.optim.AdamW(params,lr=CONFIG['lr'],weight_decay=CONFIG['weight_decay'],
        betas=tuple(CONFIG['betas']),eps=CONFIG['eps'])
    began=time.perf_counter(); torch.cuda.reset_peak_memory_stats(); step=0; epochs=[]
    with (OUT/'training.jsonl').open('x') as handle:
        for epoch in range(CONFIG['proposed_passes']):
            started=time.perf_counter(); counts=Counter(); total=0.
            for ordinal,draw in enumerate(draws):
                ident=draw['window_id']; w=windows[ident]; rec=records[ident]
                cached=torch.load(ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
                features,dna=R.C.R.Q.window_inputs(cached,w,'cuda')
                if features.requires_grad or dna.requires_grad:
                    raise ValueError('Original encoder inputs must remain detached')
                hidden,local=local_evidence(source,features,dna,w['valid_bases'])
                optimizer.zero_grad(set_to_none=True)
                chains,labels=[r['chain'] for r in rec['candidates']],[r['label'] for r in rec['candidates']]
                scores=model(hidden,local,chains,w['valid_bases'])
                result=null_anchored_chain_set_loss(scores,chains,labels)
                weighted=result['loss']*draw['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(weighted)) or any(not bool(torch.isfinite(v).all()) for v in scores.values()):
                    raise ValueError('Finite unchanged A41 composite objective and potentials required')
                weighted.backward()
                if any(p.requires_grad or p.grad is not None for p in source.parameters()):
                    raise ValueError('Frozen A3 source received gradient')
                if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                    raise ValueError('Actual affine parameter gradient missing or nonfinite')
                norm=torch.nn.utils.clip_grad_norm_(params,CONFIG['clip_norm'],error_if_nonfinite=True)
                optimizer.step(); step+=1
                row=dict(step=step,epoch=epoch+1,draw_ordinal=ordinal,window_id=ident,species=w['species'],
                    importance_weight=draw['importance_weight_uniform_within_stratum'],weighted_loss=float(weighted.detach()),
                    unweighted_mean_known_BCE=float(result['mean_known_BCE'].detach()),
                    unweighted_mean_known_set_NLL=float(result['mean_known_set_NLL'].detach()),
                    unweighted_fixed_composite_loss=float(result['loss'].detach()),
                    log_partition=float(result['log_partition'].detach()),gold_gain=float(result['gold_gain'].detach()),
                    gradient_norm_before_clip=float(norm),frozen_source_gradients_absent=True,**rec['counts'])
                handle.write(json.dumps(row,allow_nan=False)+'\n'); counts.update(rec['counts']); total+=row['weighted_loss']
                if step%64==0:
                    handle.flush(); print(json.dumps(dict(step=step,total=4608,seconds=time.perf_counter()-began)),flush=True)
                del cached,features,dna,hidden,local,scores,result,weighted
            state=model.state_dict()
            if (any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in state.values()) or
                    any(not torch.equal(v,before[k]) for k,v in source.state_dict().items())):
                raise ValueError('Finite residual and unchanged frozen source required')
            path=OUT/('step_'+str(step).zfill(6)+'.pt')
            with path.open('xb') as file:
                torch.save(dict(arm='A41',configuration=CONFIG,protocol=PROTOCOL,step=step,
                    revision=R.C.I.F.REVISION,model=state,optimizer=optimizer.state_dict(),
                    frozen_source_model=source.state_dict(),frozen_source_checkpoint=CONFIG['frozen_source']),file)
            epochs.append(dict(epoch=epoch+1,steps=1536,counts=dict(counts),mean_weighted_loss=total/1536,
                seconds=time.perf_counter()-started,checkpoint=str(path.relative_to(ROOT))))
    if step!=4608: raise ValueError('Complete fixed three-pass exposure required')
    exposure={k:int(optimizer.state[p]['step']) for k,p in model.named_parameters()}
    R.C.R.write_summary(OUT,dict(completed=True,arm='A41',configuration=CONFIG,protocol=PROTOCOL,
        preflight_job=PREFLIGHT_JOB,feature_revision=R.C.I.F.REVISION,feature_spec=R.C.R.FEATURE_SPEC,
        optimizer_steps=step,optimizer_steps_per_parameter=exposure,epochs=epochs,trainable_parameters=1806,
        final_checkpoint=str((OUT/'step_004608.pt').relative_to(ROOT)),initialization=CONFIG['initialization'],
        previous_closed_jobs=dict(fit=H.I.FIT_JOB,readout=H.I.READOUT_JOB,CPU=H.I.CHECK_JOB),
        previous_closed_five_conditions=previous_gates,frozen_source_checkpoint=CONFIG['frozen_source'],
        frozen_source_weights_unchanged=True,frozen_source_gradients_absent=True,frozen_source_optimizer_steps=0,
        frozen_supervised_A3_checkpoint_read=True,original_TRAIN_catalog_metadata_read=True,
        reference_injected_training_candidates=0,generator_decoding=False,generator_updated=False,
        new_candidate_search=False,new_backbone_forward=False,new_R4_labelling=False,
        DEV_used_for_selection=False,test_setaria_used=False,A33_outputs_opened=False,
        wall_seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One fixed affine correction fit on original TRAIN and a frozen supervised prior; '
            'not increased nonlinear expressivity, matched exposure, complete annotation or generalization.'))


if __name__=='__main__':
    torch.set_num_threads(2); main()
