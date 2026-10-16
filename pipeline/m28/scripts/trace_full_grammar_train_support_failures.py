#!/usr/bin/env python3
"""Locate fixed TRAIN counterexamples without changing observations or grammar."""
from collections import Counter
import json
import os
import re
import time
import numpy as np
import torch

import probe_full_grammar_train_core_support as S
from src.m28.full_grammar_context_view import context_geometry
from src.m28.full_grammar_boundary_compiled import ALPHABET,_state_table,_forward_messages,_advance
from src.m28.labels import oriented_parts

Q,W,P,ROOT=S.Q,S.W,S.P,S.ROOT
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-TRACE'
SOURCE_JOB='13429444'
CONFIG=dict(source_configuration=S.CONFIG,
    selection='first_original_ordinal_infeasible_event_per_species_and_strand_up_to_four',
    diagnostic='first_empty_prefix_or_missing_final_background_same_exact_automaton',
    primary_checks='original_primary_at_first_blocked_base_within_the_same_fixed_DNA_view',
    max_exons=128,checkpoint_bp=256,role_weights=1,null=0,score_accumulation='FP64',
    original_observations='all_R4_fine_coarse_unknown_rows_unchanged',optimizer_steps=0)


def trace(context,mask):
    states,destinations,roles,bg=_state_table(128)
    codes=np.array([ALPHABET.find(base) for base in context],dtype=np.int8)
    gain=np.zeros((len(context),15),dtype=np.float64)
    initial=np.full(len(states),-np.inf,dtype=np.float64); initial[bg]=0.
    positions=np.array([0]+list(range(256,len(context),256))+[len(context)],dtype=np.int64)
    final,checkpoints=_forward_messages(codes,gain,mask,initial,positions,destinations,roles)
    if np.isnan(final).any() or np.isposinf(final).any() or np.isfinite(final[bg]):
        raise ValueError('The original zero-support failure must reproduce exactly')
    empty=np.where(~np.isfinite(checkpoints).any(axis=1))[0]
    if not len(empty):
        active=np.where(np.isfinite(final))[0]
        return dict(mechanism='ONLY_UNFINISHED_STATES_AT_CLOSED_OUTER_RIGHT_BOUNDARY',
            first_blocked_oriented_base=None,first_empty_prefix_bp=None,
            final_supported_states=len(active),final_states_by_stage=dict(Counter(states[i][0] for i in active)),
            final_state_examples=[list(states[i]) for i in active[:8]],same_zero_support_failure_reproduced=True)
    block=int(empty[0]); alpha=checkpoints[block-1].copy()
    for pos in range(int(positions[block-1]),int(positions[block])):
        following=_advance(alpha,codes[pos],gain[pos],mask[pos],destinations,roles)
        if not np.isfinite(following).any():
            active=np.where(np.isfinite(alpha))[0]
            possible=set()
            for i in active:
                for edge in range(3):
                    if destinations[codes[pos],i,edge]<0:
                        break
                    possible.add(int(roles[codes[pos],i,edge]))
            required=set(np.where(mask[pos])[0].tolist())
            if not len(active) or possible & required:
                raise ValueError('Exact first-dead-base transition evidence contradicts the kernel')
            return dict(mechanism='OBSERVATION_EXCLUDES_ALL_DNA_LEGAL_TRANSITIONS' if possible else
                'NO_DNA_LEGAL_TRANSITION_FROM_REQUIRED_PRIOR_STATES',
                first_blocked_oriented_base=pos,first_empty_prefix_bp=pos+1,
                actual_base=context[pos],DNA_legal_role_options=sorted(possible),allowed_observed_roles=sorted(required),
                prior_supported_states=len(active),prior_states_by_stage=dict(Counter(states[i][0] for i in active)),
                prior_state_examples=[list(states[i]) for i in active[:8]],
                nearby_oriented_DNA=context[max(0,pos-8):min(len(context),pos+9)],
                nearby_oriented_interval=[max(0,pos-8),min(len(context),pos+9)],
                same_zero_support_failure_reproduced=True)
        alpha=following
    raise ValueError('The first empty checkpoint must contain the first empty prefix')


def main():
    began=time.perf_counter(); source=json.loads((S.OUT/'summary.json').read_text())
    if ((S.OUT/'JOBID').read_text().strip()!=SOURCE_JOB or (S.OUT/'STATUS').read_text().strip()!='COMPLETED'
            or source['job_id']!=SOURCE_JOB or not source['completed'] or source['configuration']!=S.CONFIG
            or not source['full_cohort_support_evaluated'] or source['optimizer_steps_this_run']!=0
            or source['original_train_unique_windows']!=1491 or source['DEV_used']):
        raise ValueError('Actual completed unchanged full TRAIN census required')
    failures=[r for r in source['records'] if not r['support_feasible']]
    selected={}
    for r in failures:
        if r['support_evidence']!='NEW_ZERO_SUPPORT':
            raise ValueError('A positive source support cannot become a failed census event')
        selected.setdefault((r['species'],r['strand']),r)
    selected=list(selected.values())
    S.write_json(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        source_job=SOURCE_JOB,full_population_infeasible=len(failures),
        selected_window_ids=[r['window_id'] for r in selected]))
    draws,index,windows,catalogs=Q.load_train_inputs()
    if [r['window_id'] for r in source['records']]!=list(windows):
        raise ValueError('Complete original TRAIN identity or order changed')
    sequences,descriptors=P.read_train_DNA(); opened=set(); records=[]
    with torch.inference_mode(), (OUT/'windows.jsonl').open('x') as handle:
        for original in selected:
            t=time.perf_counter(); w=windows[original['window_id']]; n=w['valid_bases']
            chromosome=sequences[w['species']][w['seqid']]
            geometry=context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
            raw=chromosome[w['start']:w['end']]; core=raw if w['strand']=='+' else Q.M.rc(raw)
            cached=W.load_cache(w,index,opened)
            if cached['sequence']!=core:
                raise ValueError('Actual original core DNA changed')
            del cached
            targets=P.P.build_targets(w,catalogs[w['species']],core)
            allowed=P.P.allowed_roles_from_local_targets(targets,n)
            fine=int((torch.as_tensor(targets['seg'][:n])>=0).sum())
            free=int(allowed.all(1).sum()); coarse=n-fine-free
            genomic=re.sub('[^ACGT]','N',chromosome[a:b])
            context=genomic if w['strand']=='+' else Q.M.rc(genomic)
            offset=w['start']-a if w['strand']=='+' else b-w['end']
            if ([a,b]!=original['genomic_context'] or [offset,offset+n]!=original['core_oriented_interval']
                    or (fine,coarse,free)!=(original['fine_positions'],original['coarse_positions'],original['fully_unknown_positions'])
                    or context[offset:offset+n]!=re.sub('[^ACGT]','N',core)):
                raise ValueError('Original fixed view or observations changed')
            mask=np.ones((len(context),15),dtype=np.bool_); mask[offset:offset+n]=allowed.numpy()
            result=trace(context,mask); pos=result['first_blocked_oriented_base']; primary_checks=[]
            if pos is not None:
                x=a+pos if w['strand']=='+' else b-1-pos
                in_core=offset<=pos<offset+n
                result.update(genomic_base_zero_based=x,first_blocked_base_in_original_core=in_core,
                    original_fine_label=int(targets['seg'][pos-offset]) if in_core else None,
                    original_coarse_label=int(targets['region'][pos-offset]) if in_core else None,
                    primary_transcripts_at_base=[tid for tid in original['original_primary_ids']
                        if catalogs[w['species']][tid]['CDS'][0][0]<=x<catalogs[w['species']][tid]['CDS'][-1][1]])
                for tid in result['primary_transcripts_at_base']:
                    tx=catalogs[w['species']][tid]
                    parts=oriented_parts(tx,dict(start=a,end=b,strand=w['strand']))
                    if not 0<=parts[0][0]<parts[-1][1]<=len(context):
                        primary_checks.append(dict(transcript_id=tid,DNA_checked=False,reason='outside_fixed_view'))
                        continue
                    geometry_failures=W.geometry_failures(parts,len(context))
                    DNA_failures,coding_legal=W.DNA_checks(context,parts)
                    primary_checks.append(dict(transcript_id=tid,gene_id=tx['gene_id'],
                        CDS_complete_supported=tx['CDS_complete_supported'],eligible_positive=tx['eligible_positive'],
                        oriented_CDS_parts=[list(p) for p in parts],DNA_checked=True,
                        geometry_failures=geometry_failures,direct_DNA_failures=DNA_failures,
                        actual_carry_coding_legal=coding_legal,
                        complete_original_chain_legal=not geometry_failures and not DNA_failures))
            row=dict(original_ordinal=original['ordinal'],window_id=w['id'],species=w['species'],strand=w['strand'],
                genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
                original_primary_ids=original['original_primary_ids'],outside_primary_ids=original['outside_primary_ids'],
                all_core_observation_rows_preserved=True,primary_DNA_checks_at_blocked_base=primary_checks,
                **result,seconds=time.perf_counter()-t)
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
    if opened!={index[r['window_id']]['artifact'] for r in selected}:
        raise ValueError('Only the predeclared failed witnesses may open their original cache')
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_job=SOURCE_JOB,
        full_population_infeasible=len(failures),selected_failed_witnesses=len(selected),records=records,
        mechanism_counts=dict(Counter(r['mechanism'] for r in records)),
        all_failed_windows_have_complete_primary_spans_inside_view=all(not r['outside_primary_ids'] for r in failures),
        diagnostic_covers_all_failed_mechanisms=False,original_TRAIN_DNA_sources=descriptors,
        cached_artifacts_opened=sorted(opened),seconds_after_module_import=time.perf_counter()-began,
        interpretation='localized_original_counterexamples_not_label_repair_or_full_failure_causal_census',**S.FLAGS)
    S.write_json(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
