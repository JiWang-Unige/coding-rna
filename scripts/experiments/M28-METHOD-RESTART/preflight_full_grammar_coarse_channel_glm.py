#!/usr/bin/env python3
"""Four predetermined original negative contexts; native GLM, zero updates."""
import json
import math
import os
import re
import time
import torch
import preflight_full_grammar_glm_context as G
import probe_full_grammar_coarse_channel_train as N
import trace_full_grammar_train_support_failures as T
from src.m28.full_grammar_coarse_channel import coarse_channel_nll

Q,W,P,ROOT=G.Q,G.W,G.P,G.ROOT
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-GLM-PREFLIGHT'
CONFIG=dict(native_wiring_configuration=G.CONFIG,candidate_probe_configuration=N.CONFIG,
    selection='same_four_first_original_ordinal_hard_negatives_per_species_and_strand_from_13429503',
    loss='unnormalized_candidate_coarse_channel_NLL',epsilon=0.1,
    epsilon_diagnostic_only_not_fit_selection=True,optimizer_steps=0)


def main():
    began=time.perf_counter(); sources=[]
    for root,job,config in ((N.OUT,'13429625',N.CONFIG),(T.OUT,'13429503',T.CONFIG),(G.OUT,'13429053',G.CONFIG)):
        s=json.loads((root/'summary.json').read_text())
        if ((root/'STATUS').read_text().strip()!='COMPLETED' or (root/'JOBID').read_text().strip()!=job
                or not s['completed'] or s['job_id']!=job or s['configuration']!=config
                or s['optimizer_steps_this_run']!=0 or s['DEV_used']):
            raise ValueError('Actual complete unchanged candidate, selection and native sources required')
        sources.append(s)
    candidate,selection,native=sources
    if (not candidate['all_original_negative_candidates_feasible'] or candidate['candidate_finite_events_and_gradients']!=24
            or not native['all_eight_NLL_role_and_shared_parameter_gradients_completed']
            or not native['new_model_weights_unchanged'] or not native['pretrained_encoder_frozen_no_grad']):
        raise ValueError('Actual all-24 candidate support and unchanged prior native wiring required')
    rows={r['window_id']:r for r in candidate['records']}; ids=[r['window_id'] for r in selection['records']]
    if len(ids)!=4 or not set(ids)<=set(rows): raise ValueError('The same four predetermined original negatives required')
    N.S.write_json(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        source_jobs=['13429625','13429503','13429053'],selected_window_ids=ids))
    _,index,windows,catalogs=Q.load_train_inputs(); sequences,descriptors=P.read_train_DNA()
    backbone,tokenizer=G.load_native(); torch.manual_seed(0)
    model=G.M28Core('C0',2048,128).to('cuda').train()
    if sum(p.numel() for p in model.parameters())!=366191: raise ValueError('Original fresh C0 capacity changed')
    before={k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
    opened=set(); records=[]; torch.cuda.reset_peak_memory_stats()
    with (OUT/'windows.jsonl').open('x') as handle:
        for ident in ids:
            began_window=time.perf_counter(); original=rows[ident]; w=windows[ident]; n=w['valid_bases']
            chromosome=sequences[w['species']][w['seqid']]
            geometry=G.context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
            raw=chromosome[w['start']:w['end']]; core=raw if w['strand']=='+' else Q.M.rc(raw)
            cached=W.load_cache(w,index,opened)
            if cached['sequence']!=core: raise ValueError('Original core DNA changed')
            del cached
            targets=P.P.build_targets(w,catalogs[w['species']],core)
            fine=torch.as_tensor(targets['seg']); region=torch.as_tensor(targets['region'])
            counts=(int((fine[:n]>=0).sum()),int(((fine[:n]<0)&(region[:n]>=0)).sum()),int(((fine[:n]<0)&(region[:n]<0)).sum()))
            offset=w['start']-a if w['strand']=='+' else b-w['end']
            if ([a,b]!=original['genomic_context'] or [offset,offset+n]!=original['core_oriented_interval']
                    or counts!=(original['fine_positions'],original['coarse_positions'],original['fully_unknown_positions'])):
                raise ValueError('Original fixed view or all observation rows changed')
            context=re.sub('[^ACGT]','N',chromosome[a:b])
            if w['strand']=='-': context=Q.M.rc(context)
            model.zero_grad(set_to_none=True); tensors=[]; encoding=[]
            for left,right in geometry['genomic_windows']:
                dna_string=re.sub('[^ACGT]','N',chromosome[left:right])
                if w['strand']=='-': dna_string=Q.M.rc(dna_string)
                features=G.encode_window(backbone,tokenizer,dna_string,'cuda').clone()
                dna=G.one_hot(dna_string,Q.M.WIDTH)[None].to('cuda')
                logits=model(features,dna,[right-left])['segmentation_logits'][0]
                logits.retain_grad(); tensors.append(logits)
                encoding.append(dict(genomic_window=[left,right],feature_shape=list(features.shape),
                    features_require_grad=features.requires_grad,feature_is_inference_tensor=torch.is_inference(features)))
            values=[v[:right-left] for v,(left,right) in zip(tensors,geometry['genomic_windows'])]
            joined=G.owned_context_logits(values,geometry,w['strand']).to('cpu'); joined.retain_grad()
            if tuple(joined.shape)!=(len(context),15) or not bool(torch.isfinite(joined).all()):
                raise ValueError('Actual finite neural context potentials required')
            epsilon=torch.tensor(0.1,dtype=torch.float64,requires_grad=True)
            cuts=tuple(x for x in (offset,offset+n) if 0<x<len(context))
            result=coarse_channel_nll(context,joined,128,offset,offset+n,targets,epsilon,cuts=cuts,checkpoint_bp=256)
            masses={k:float(result[k].detach()) for k in ('loss','log_partition_all','log_partition_observed','background_correction')}
            if not all(math.isfinite(v) for v in masses.values()) or masses['loss'] < -1e-8:
                raise ValueError('Finite neural candidate likelihood required')
            result['loss'].backward(); torch.cuda.synchronize(); gradient=joined.grad
            if gradient is None or not bool(torch.isfinite(gradient).all()) or epsilon.grad is None or not bool(torch.isfinite(epsilon.grad)):
                raise ValueError('Actual neural context and epsilon first gradients must be finite')
            gauge=float(gradient.sum(1).abs().max()); bound=float(gradient.abs().max())
            if gauge>1e-5 or bound>1+1e-6: raise ValueError('Actual neural probability gradient bounds violated')
            for value,(left,right),(_,lo,hi) in zip(tensors,geometry['local_windows'],geometry['base_owned_segments']):
                expected=torch.zeros(value.shape,dtype=torch.float32)
                x,y=(lo-left,hi-left) if w['strand']=='+' else (right-hi,right-lo)
                g,h=(lo,hi) if w['strand']=='+' else (len(context)-hi,len(context)-lo)
                expected[x:y]=gradient[g:h]
                if value.grad is None or not torch.equal(value.grad.cpu(),expected):
                    raise ValueError('Actual selected owner-row VJP differs')
            groups={}
            for name,module in (('project',model.shared.project),('context',model.shared.context),
                    ('dna',model.shared.dna),('fuse',model.shared.fuse),('segmentation',model.segmentation)):
                parameters=list(module.parameters())
                if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in parameters):
                    raise ValueError('Missing or nonfinite native parameter gradient: '+name)
                groups[name]=math.sqrt(sum(float(p.grad.double().square().sum()) for p in parameters))
                if groups[name]<=0: raise ValueError('Zero native parameter gradient: '+name)
            row=dict(original_ordinal=original['original_ordinal'],window_id=ident,status='FINITE_NATIVE_CANDIDATE_NLL_AND_GRADIENTS',
                geometry=geometry,core_oriented_interval=[offset,offset+n],fine_positions=counts[0],coarse_positions=counts[1],
                fully_unknown_positions=counts[2],encoding=encoding,original_hard_event_support_feasible=False,**masses,
                all_context_role_gradients_finite=True,shared_parameter_gradient_l2=groups,
                all_shared_parameter_gradients_finite=True,owned_window_row_VJP_exact=True,
                row_gauge_max_abs=gauge,maximum_abs_role_gradient=bound,epsilon_gradient=float(epsilon.grad),
                all_core_observation_rows_preserved=True,original_labels_changed=False,
                coarse_observation_rows_dropped=False,inference_grammar_changed=False,seconds=time.perf_counter()-began_window)
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
            del result,joined,gradient,tensors,values,features,dna,logits,targets,epsilon
    if any(not torch.equal(v.cpu(),before[k]) for k,v in model.state_dict().items()):
        raise ValueError('Zero-update diagnostic changed fresh model state')
    if any(p.requires_grad or p.grad is not None for p in backbone.parameters()):
        raise ValueError('Released GLM must remain frozen with no gradients')
    if opened!={index[ident]['artifact'] for ident in ids}: raise ValueError('Only the same four original negative caches required')
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_jobs=['13429625','13429503','13429053'],
        records=records,actual_predetermined_negative_contexts=len(records),
        actual_released_GLM_forward=True,actual_fresh_C0_forward=True,pretrained_encoder_frozen_no_grad=True,
        new_model_weights_unchanged=True,optimizer_constructed=False,optimizer_steps_this_run=0,checkpoints_written=False,
        original_R4_policy_changed=False,observation_rows_deleted=False,inference_grammar_changed=False,
        candidate_observation_model_changed=True,fit_error_probability_selected=False,
        approved_fit_context_or_raw_chromosome_protocol_selected=False,loss_normalizer_or_fit_coefficients_selected=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,A44_execution_chain_reopened=False,
        cached_features_used=False,old_cached_labels_used=False,original_TRAIN_DNA_sources=descriptors,
        cached_artifacts_opened=sorted(opened),peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),seconds=time.perf_counter()-began,
        interpretation='candidate_native_neural_gradient_wiring_not_fit_accuracy_or_generalization')
    N.S.write_json(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(2); main()
