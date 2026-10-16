#!/usr/bin/env python3
"""Candidate likelihood and first gradients on all original hard negatives."""
from collections import Counter
import json
import math
import os
import re
import time
import torch
import probe_full_grammar_train_core_support as S
import check_full_grammar_coarse_channel as C
from src.m28.full_grammar_context_view import context_geometry
from src.m28.full_grammar_coarse_channel import coarse_channel_nll
from src.m28.full_grammar_observed_partition import InfeasibleObservation

Q,W,P,ROOT=S.Q,S.W,S.P,S.ROOT
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-TRAIN-PROBE'
CONFIG=dict(source_configuration=S.CONFIG,math_configuration=C.CONFIG,
    scope='all_24_original_hard_event_negative_TRAIN_windows_no_resampling',
    potentials='zero_CPU_FP64',max_exons=128,checkpoint_bp=256,
    epsilon=0.1,epsilon_diagnostic_only_not_fit_selection=True,
    fine='unchanged_exact_role',coarse='normalized_three_class_observation_channel',
    targets='all_original_R4_fine_coarse_unknown_rows_unchanged',
    context='unchanged_raw_coordinate_neighbor_view_no_annotation_expansion',
    gradient_gauge_atol=1e-8,posterior_difference_bound_atol=1e-9,optimizer_steps=0)


def main():
    began=time.perf_counter(); source=json.loads((S.OUT/'summary.json').read_text())
    math_source=json.loads((C.OUT/'summary.json').read_text())
    for root,s,job,config in ((S.OUT,source,'13429444',S.CONFIG),(C.OUT,math_source,'13429599',C.CONFIG)):
        if ((root/'STATUS').read_text().strip()!='COMPLETED' or (root/'JOBID').read_text().strip()!=job
                or not s['completed'] or s['job_id']!=job or s['configuration']!=config
                or s['optimizer_steps_this_run']!=0):
            raise ValueError('Matching actual completed support and mathematical sources required')
    if math_source['tests_run']!=5 or math_source['failures'] or math_source['errors']:
        raise ValueError('All five original artificial mathematical checks required')
    failures=[r for r in source['records'] if not r['support_feasible']]
    if len(source['records'])!=1491 or len(failures)!=24: raise ValueError('Original whole TRAIN support population required')
    S.write_json(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        source_jobs=['13429444','13429599'],selected_window_ids=[r['window_id'] for r in failures]))
    _,index,windows,catalogs=Q.load_train_inputs(); sequences,descriptors=P.read_train_DNA()
    if [r['window_id'] for r in source['records']]!=list(windows): raise ValueError('Original TRAIN order changed')
    records=[]; opened=set()
    with (OUT/'windows.jsonl').open('x') as handle:
        for original in failures:
            began_window=time.perf_counter(); w=windows[original['window_id']]; n=w['valid_bases']
            chromosome=sequences[w['species']][w['seqid']]
            geometry=context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
            raw=chromosome[w['start']:w['end']]; core=raw if w['strand']=='+' else Q.M.rc(raw)
            cached=W.load_cache(w,index,opened)
            if cached['sequence']!=core: raise ValueError('Original core DNA changed')
            del cached
            targets=P.P.build_targets(w,catalogs[w['species']],core)
            fine=torch.as_tensor(targets['seg']); region=torch.as_tensor(targets['region'])
            counts=(int((fine[:n]>=0).sum()),int(((fine[:n]<0)&(region[:n]>=0)).sum()),
                int(((fine[:n]<0)&(region[:n]<0)).sum()))
            context=re.sub('[^ACGT]','N',chromosome[a:b])
            if w['strand']=='-': context=Q.M.rc(context)
            offset=w['start']-a if w['strand']=='+' else b-w['end']
            if ([a,b]!=original['genomic_context'] or [offset,offset+n]!=original['core_oriented_interval']
                    or counts!=(original['fine_positions'],original['coarse_positions'],original['fully_unknown_positions'])
                    or context[offset:offset+n]!=re.sub('[^ACGT]','N',core)):
                raise ValueError('Original view or observation rows changed')
            before_fine,before_region=fine.clone(),region.clone()
            logits=torch.zeros((len(context),15),dtype=torch.float64,requires_grad=True)
            epsilon=torch.tensor(CONFIG['epsilon'],dtype=torch.float64,requires_grad=True)
            row=dict(original_ordinal=original['ordinal'],window_id=w['id'],species=w['species'],strand=w['strand'],
                genomic_context=[a,b],core_oriented_interval=[offset,offset+n],context_length=len(context),
                fine_positions=counts[0],coarse_positions=counts[1],fully_unknown_positions=counts[2],
                original_hard_event_support_feasible=False,raw_core_matches_original_cache=True)
            try:
                result=coarse_channel_nll(context,logits,128,offset,offset+n,targets,epsilon,checkpoint_bp=256)
                gx,ge=torch.autograd.grad(result['loss'],(logits,epsilon))
                values={k:float(result[k].detach()) for k in ('loss','log_partition_all','log_partition_observed',
                    'weighted_relative_log_partition','background_correction')}
                gauge=float(gx.sum(1).abs().max()); maximum=float(gx.abs().max()); derivative=float(ge)
                if (not all(math.isfinite(v) for v in values.values()) or values['loss'] < -1e-8
                        or not bool(torch.isfinite(gx).all()) or not math.isfinite(derivative)
                        or gauge>CONFIG['gradient_gauge_atol']
                        or maximum>1+CONFIG['posterior_difference_bound_atol']
                        or not -counts[1]/CONFIG['epsilon']-1e-6 <= derivative <= counts[1]/(1-CONFIG['epsilon'])+1e-6):
                    raise ValueError('Candidate likelihood or actual first gradient violates probability identity')
                if (result['fine_positions'],result['coarse_positions'],result['fully_unknown_positions'])!=counts:
                    raise ValueError('Candidate must use every original observation row')
                row.update(status='FINITE_CANDIDATE_LIKELIHOOD_AND_GRADIENT',candidate_event_feasible=True,
                    candidate_NLL=values.pop('loss'),**values,all_context_role_gradients_finite=True,
                    context_role_gradient_max_abs=maximum,context_role_gradient_nonzero=bool(gx.abs().max()>0),
                    row_gauge_max_abs=gauge,epsilon_gradient=derivative,
                    epsilon_gradient_finite=True,epsilon_gradient_nonzero=derivative!=0,
                    partition_statistics_all=result['statistics_all'],partition_statistics_observed=result['statistics_observed'])
                del result,gx,ge
            except InfeasibleObservation as exc:
                row.update(status='STILL_INFEASIBLE_EXACT_FINE_EVENT',candidate_event_feasible=False,
                    infeasible_reason=str(exc),candidate_NLL=None)
            if not torch.equal(fine,before_fine) or not torch.equal(region,before_region):
                raise ValueError('Original targets were modified')
            row.update(all_core_observation_rows_preserved=True,original_labels_changed=False,
                coarse_observation_rows_dropped=False,inference_grammar_changed=False,
                seconds=time.perf_counter()-began_window)
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
            del logits,epsilon,targets
    if opened!={index[r['window_id']]['artifact'] for r in failures}: raise ValueError('Only the 24 actual negative core caches required')
    feasible=sum(r['candidate_event_feasible'] for r in records)
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_jobs=['13429444','13429599'],
        records=records,actual_original_hard_negatives=24,candidate_finite_events_and_gradients=feasible,
        candidate_still_infeasible_events=24-feasible,all_original_negative_candidates_feasible=feasible==24,
        statuses=dict(Counter(r['status'] for r in records)),
        candidate_finite_by_species=dict(Counter(r['species'] for r in records if r['candidate_event_feasible'])),
        original_hard_event_failures_retained=24,candidate_observation_model_changed=True,
        fit_error_probability_selected=False,loss_normalizer_selected=False,
        whole_TRAIN_candidate_numerical_NLL_evaluated=False,full_chromosome_protocol_proven=False,
        original_TRAIN_DNA_sources=descriptors,cached_artifacts_opened=sorted(opened),
        seconds=time.perf_counter()-began,
        interpretation='candidate_fixed_view_likelihood_first_gradients_not_fit_or_performance',**S.FLAGS)
    S.write_json(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(2); main()
