#!/usr/bin/env python3
"""All original TRAIN core events in the fixed raw view; no model or fit."""
import argparse
from collections import Counter
import json
import math
import os
import re
import time

import torch
import preflight_full_grammar_glm_context as G
import probe_full_grammar_context_extent_warnings as E
from src.m28.full_grammar_context_view import context_geometry
from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition
from src.m28.full_grammar_observed_partition import InfeasibleObservation

Q,W,P,ROOT = G.Q,G.W,G.P,G.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-CORE-SUPPORT'
SHARDS = 8
CONFIG = dict(source_view_configuration=G.CONFIG,
    scope='all_1491_original_TRAIN_unique_core_events_no_resampling',
    partition='original_window_insertion_order_ordinal_modulo_eight',shards=SHARDS,
    reuse='eight_finite_neural_and_six_finite_zero_source_supports_no_partition_rerun',
    new_potentials='zero_CPU_FP32',targets='unchanged_R4_all_fine_coarse_unknown_core_rows',
    max_exons=128,checkpoint_bp=256,role_weights=1,null=0,score_accumulation='FP64',
    context='unchanged_raw_coordinate_neighbor_view_no_annotation_expansion',optimizer_steps=0)
FLAGS = dict(original_R4_policy_changed=False,observation_rows_deleted=False,
    primary_transcripts_dropped=False,context_expanded_from_annotations=False,
    model_or_backbone_forward=False,weights_read=False,optimizer_steps_this_run=0,
    checkpoints_written=False,cached_features_used=False,old_cached_labels_used=False,
    DEV_used=False,independent_species_used=False,A33_outputs_opened=False,
    A44_execution_chain_reopened=False,complete_original_primary_coverage_proven=False)


def write_json(path,value):
    with path.open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False)


def source_supports():
    reused = {}
    for root,job,config in ((G.OUT,'13429053',G.CONFIG),(E.OUT,'13429243',E.CONFIG)):
        s=json.loads((root/'summary.json').read_text())
        if ((root/'JOBID').read_text().strip()!=job or (root/'STATUS').read_text().strip()!='COMPLETED'
                or not s['completed'] or s['job_id']!=job or s['configuration']!=config
                or s['optimizer_steps_this_run']!=0 or s['DEV_used']):
            raise ValueError('Matching actual completed fixed-view source required')
        rows=s['records']
        if len(rows)!=(8 if root==G.OUT else 6):
            raise ValueError('All fourteen original source events required')
        for row in rows:
            ident=row['window_id']
            if ident in reused:
                raise ValueError('Source support identities overlap')
            if root==G.OUT:
                mass=row['log_partition_observed']; geometry=row['geometry']['genomic_context']
                if not row['all_context_role_gradients_finite'] or not row['raw_core_matches_original_cache']:
                    raise ValueError('Finite actual neural source evidence required')
                kind='REUSED_FINITE_NEURAL_SUPPORT'
            else:
                mass=row['zero_potential_observed_log_count']; geometry=row['diagnostic_genomic_context']
                if not row['core_observation_feasible'] or not row['all_core_observation_rows_preserved']:
                    raise ValueError('Finite original zero-potential source evidence required')
                kind='REUSED_ZERO_SUPPORT'
            if not math.isfinite(mass):
                raise ValueError('Source observed log partition is not finite')
            reused[ident]=dict(row=row,source_job=job,kind=kind,mass=mass,geometry=geometry)
    return reused


def run_shard(shard):
    began=time.perf_counter(); out=OUT/f'shard_{shard}'
    reused=source_supports()
    draws,index,windows,catalogs=Q.load_train_inputs()
    selected=[(i,w) for i,w in enumerate(windows.values()) if i%SHARDS==shard]
    if len(draws)!=1536 or len(windows)!=1491 or not set(reused)<=set(windows):
        raise ValueError('Original complete TRAIN population required')
    write_json(out/'started.json',dict(configuration=CONFIG,shard=shard,
        array_job_id=os.environ['SLURM_ARRAY_JOB_ID'],slurm_job_id=os.environ['SLURM_JOB_ID'],
        selected_window_ids=[w['id'] for _,w in selected],source_jobs=['13429053','13429243'],
        TRAIN_allowlists={s:sorted(v) for s,v in Q.ALLOW.items()}))
    sequences,descriptors=P.read_train_DNA(); records=[]; opened=set()
    with torch.inference_mode(), (out/'windows.jsonl').open('x') as handle:
        for ordinal,w in selected:
            t=time.perf_counter(); ident=w['id']; n=w['valid_bases']
            chromosome=sequences[w['species']][w['seqid']]
            geometry=context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
            raw=chromosome[w['start']:w['end']]; core=raw if w['strand']=='+' else Q.M.rc(raw)
            cached=W.load_cache(w,index,opened)
            if cached['sequence']!=core:
                raise ValueError('Original raw core differs from original cache: '+ident)
            del cached
            targets=P.P.build_targets(w,catalogs[w['species']],core)
            allowed=P.P.allowed_roles_from_local_targets(targets,n)
            fine=int((torch.as_tensor(targets['seg'][:n])>=0).sum())
            free=int(allowed.all(1).sum()); coarse=n-fine-free
            genomic=re.sub('[^ACGT]','N',chromosome[a:b])
            context=genomic if w['strand']=='+' else Q.M.rc(genomic)
            offset=w['start']-a if w['strand']=='+' else b-w['end']
            if context[offset:offset+n]!=re.sub('[^ACGT]','N',core):
                raise ValueError('Oriented raw core or coordinates changed')
            primary=[tid for tid in w['overlapping_transcript_ids'] if catalogs[w['species']][tid]['training_primary']]
            outside=[tid for tid in primary if catalogs[w['species']][tid]['CDS'][0][0]<a
                     or catalogs[w['species']][tid]['CDS'][-1][1]>b]
            reason=None
            if ident in reused:
                previous=reused[ident]; old=previous['row']
                if (previous['geometry']!=[a,b] or old['core_oriented_interval']!=[offset,offset+n]
                        or (old['fine_positions'],old['coarse_positions'],old['fully_unknown_positions'])!=(fine,coarse,free)):
                    raise ValueError('Original source view or observation counts changed: '+ident)
                feasible=True; mass=previous['mass']; kind=previous['kind']; source=previous['source_job']
            else:
                mask=torch.ones((len(context),15),dtype=torch.bool); mask[offset:offset+n]=allowed
                potentials=torch.zeros((len(context),15),dtype=torch.float32)
                cuts=tuple(x for x in (offset,offset+n) if 0<x<len(context))
                print(json.dumps(dict(stage='new_core_support_started',ordinal=ordinal,window_id=ident)),flush=True)
                try:
                    result=compiled_boundary_full_grammar_log_partition(context,potentials,128,
                        allowed=mask,cuts=cuts,checkpoint_bp=256)
                    mass=float(result['log_partition'])
                    if not math.isfinite(mass) or mass < -1e-8:
                        raise ValueError('Finite nonnegative zero-potential observed log count required')
                    feasible=True
                except InfeasibleObservation as error:
                    feasible=False; mass=None; reason=str(error)
                del mask,potentials
                kind='NEW_ZERO_SUPPORT'; source=f"{os.environ['SLURM_ARRAY_JOB_ID']}_{shard}"
            row=dict(ordinal=ordinal,window_id=ident,species=w['species'],seqid=w['seqid'],strand=w['strand'],
                valid_core_bases=n,genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
                context_length=len(context),fine_positions=fine,coarse_positions=coarse,fully_unknown_positions=free,
                original_primary_ids=primary,outside_primary_ids=outside,raw_core_matches_original_cache=True,
                all_core_observation_rows_preserved=True,support_feasible=feasible,infeasible_reason=reason,
                observed_log_partition=mass,support_evidence=kind,support_source_job=source,
                status='FEASIBLE_CORE_EVENT' if feasible else 'INFEASIBLE_CORE_EVENT',seconds=time.perf_counter()-t)
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
    if opened!={index[w['id']]['artifact'] for _,w in selected}:
        raise ValueError('Exactly this shard original TRAIN caches required')
    summary=dict(completed=True,configuration=CONFIG,shard=shard,array_job_id=os.environ['SLURM_ARRAY_JOB_ID'],
        slurm_job_id=os.environ['SLURM_JOB_ID'],original_draws=len(draws),original_train_unique_windows=len(windows),
        records=records,feasible=sum(r['support_feasible'] for r in records),
        infeasible=sum(not r['support_feasible'] for r in records),
        evidence_counts=dict(Counter(r['support_evidence'] for r in records)),
        cached_artifacts_opened=sorted(opened),original_TRAIN_DNA_sources=descriptors,
        full_cohort_observation_support_proven=False,seconds_after_module_import=time.perf_counter()-began,**FLAGS)
    write_json(out/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


def aggregate(array_job):
    began=time.perf_counter(); source_supports()
    draws,index,windows,catalogs=Q.load_train_inputs(); ids=list(windows)
    rows=[]; shards=[]
    for shard in range(SHARDS):
        out=OUT/f'shard_{shard}'; s=json.loads((out/'summary.json').read_text())
        part=[json.loads(line) for line in (out/'windows.jsonl').read_text().splitlines()]
        started=json.loads((out/'started.json').read_text())
        expected=ids[shard::SHARDS]
        if ((out/'STATUS').read_text().strip()!='COMPLETED' or (out/'JOBID').read_text().strip()!=f'{array_job}_{shard}'
                or not s['completed'] or s['array_job_id']!=array_job or s['shard']!=shard
                or s['configuration']!=CONFIG or started['configuration']!=CONFIG
                or started['selected_window_ids']!=expected or [r['window_id'] for r in part]!=expected
                or part!=s['records'] or any(s[k]!=v for k,v in FLAGS.items())
                or set(s['cached_artifacts_opened'])!={index[i]['artifact'] for i in expected}):
            raise ValueError('Complete matching unchanged shard and every artifact required')
        rows.extend(part); shards.append(dict(shard=shard,slurm_job_id=s['slurm_job_id'],
            windows=len(part),feasible=s['feasible'],infeasible=s['infeasible'],
            seconds_after_module_import=s['seconds_after_module_import']))
    rows.sort(key=lambda r:r['ordinal'])
    if [r['window_id'] for r in rows]!=ids or [r['ordinal'] for r in rows]!=list(range(1491)):
        raise ValueError('All original 1491 windows exactly once required')
    evidence=Counter(r['support_evidence'] for r in rows)
    if evidence!=Counter(NEW_ZERO_SUPPORT=1477,REUSED_FINITE_NEURAL_SUPPORT=8,REUSED_ZERO_SUPPORT=6):
        raise ValueError('Exactly 1477 new and fourteen unchanged source supports required')
    for r in rows:
        if (r['fine_positions']+r['coarse_positions']+r['fully_unknown_positions']!=r['valid_core_bases']
                or not r['raw_core_matches_original_cache'] or not r['all_core_observation_rows_preserved']
                or r['support_feasible']!=(r['status']=='FEASIBLE_CORE_EVENT')
                or (r['support_feasible'] and (not math.isfinite(r['observed_log_partition']) or r['infeasible_reason'] is not None))
                or (not r['support_feasible'] and (r['observed_log_partition'] is not None or not r['infeasible_reason']))):
            raise ValueError('Per-core evidence or unchanged observation coverage differs')
    counts={s:dict(windows=sum(r['species']==s for r in rows),
        feasible=sum(r['species']==s and r['support_feasible'] for r in rows),
        infeasible=sum(r['species']==s and not r['support_feasible'] for r in rows)) for s in Q.ALLOW}
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,array_job_id=array_job,
        original_draws=len(draws),original_train_unique_windows=len(rows),source_jobs=['13429053','13429243'],
        records=rows,shards=shards,species_counts=counts,evidence_counts=dict(evidence),
        feasible=sum(r['support_feasible'] for r in rows),infeasible=sum(not r['support_feasible'] for r in rows),
        full_cohort_support_evaluated=True,full_cohort_observation_support_proven=all(r['support_feasible'] for r in rows),
        original_primary_extent_warning_windows=sum(bool(r['outside_primary_ids']) for r in rows),
        seconds_after_module_import=time.perf_counter()-began,
        interpretation='fixed_view_all_original_TRAIN_core_support_not_fit_or_complete_primary_or_model_performance',**FLAGS)
    write_json(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--shard',type=int,choices=range(SHARDS)); group.add_argument('--aggregate')
    args=parser.parse_args(); torch.set_num_threads(2)
    if args.shard is not None:
        run_shard(args.shard)
    else:
        aggregate(args.aggregate)
