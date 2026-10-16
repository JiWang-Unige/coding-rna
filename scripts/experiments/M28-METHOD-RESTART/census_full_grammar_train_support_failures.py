#!/usr/bin/env python3
"""All 24 original counterexamples; reuse four traces, localize the other 20."""
from collections import Counter
import json
import os
import re
import time
import numpy as np
import torch

import trace_full_grammar_train_support_failures as T
from src.m28.full_grammar_context_view import context_geometry
from src.m28.labels import oriented_parts

S,Q,W,P,ROOT=T.S,T.Q,T.W,T.P,T.ROOT
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-CENSUS'
CONFIG=dict(source_configuration=S.CONFIG,source_trace_configuration=T.CONFIG,
    scope='all_24_actual_original_TRAIN_counterexamples_in_original_ordinal_order',
    reuse='four_actual_localizations_rechecked_against_unchanged_DNA_and_observations',
    new_localizations=20,optimizer_steps=0)


def primary_checks(context,a,b,w,ids,catalog):
    checks=[]
    for tid in ids:
        tx=catalog[tid]; parts=oriented_parts(tx,dict(start=a,end=b,strand=w['strand']))
        if not 0<=parts[0][0]<parts[-1][1]<=len(context):
            checks.append(dict(transcript_id=tid,DNA_checked=False,reason='outside_fixed_view'))
            continue
        geometry=W.geometry_failures(parts,len(context)); failures,coding=W.DNA_checks(context,parts)
        checks.append(dict(transcript_id=tid,gene_id=tx['gene_id'],CDS_complete_supported=tx['CDS_complete_supported'],
            eligible_positive=tx['eligible_positive'],oriented_CDS_parts=[list(p) for p in parts],DNA_checked=True,
            geometry_failures=geometry,direct_DNA_failures=failures,actual_carry_coding_legal=coding,
            complete_original_chain_legal=not geometry and not failures))
    return checks


def main():
    began=time.perf_counter()
    source=json.loads((S.OUT/'summary.json').read_text()); old=json.loads((T.OUT/'summary.json').read_text())
    for root,summary,job,config in ((S.OUT,source,'13429444',S.CONFIG),(T.OUT,old,'13429503',T.CONFIG)):
        if ((root/'STATUS').read_text().strip()!='COMPLETED' or (root/'JOBID').read_text().strip()!=job
                or not summary['completed'] or summary['job_id']!=job or summary['configuration']!=config
                or summary['optimizer_steps_this_run']!=0 or summary['DEV_used']):
            raise ValueError('Actual complete unchanged original sources required')
    failures=[r for r in source['records'] if not r['support_feasible']]
    reused={r['window_id']:r for r in old['records']}
    if len(failures)!=24 or len(reused)!=4 or not set(reused)<={r['window_id'] for r in failures}:
        raise ValueError('All original 24 failures and exactly four existing traces required')
    S.write_json(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        source_jobs=['13429444','13429503'],selected_window_ids=[r['window_id'] for r in failures],
        reused_window_ids=list(reused)))
    draws,index,windows,catalogs=Q.load_train_inputs(); sequences,descriptors=P.read_train_DNA()
    if [r['window_id'] for r in source['records']]!=list(windows):
        raise ValueError('Original full TRAIN identities or order changed')
    opened=set(); records=[]
    with torch.inference_mode(), (OUT/'windows.jsonl').open('x') as handle:
        for original in failures:
            t=time.perf_counter(); w=windows[original['window_id']]; n=w['valid_bases']
            chromosome=sequences[w['species']][w['seqid']]
            geometry=context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
            raw=chromosome[w['start']:w['end']]; core=raw if w['strand']=='+' else Q.M.rc(raw)
            cached=W.load_cache(w,index,opened)
            if cached['sequence']!=core:
                raise ValueError('Original core DNA changed')
            del cached
            targets=P.P.build_targets(w,catalogs[w['species']],core)
            allowed=P.P.allowed_roles_from_local_targets(targets,n)
            fine=int((torch.as_tensor(targets['seg'][:n])>=0).sum()); free=int(allowed.all(1).sum())
            context=re.sub('[^ACGT]','N',chromosome[a:b])
            if w['strand']=='-': context=Q.M.rc(context)
            offset=w['start']-a if w['strand']=='+' else b-w['end']
            if ([a,b]!=original['genomic_context'] or [offset,offset+n]!=original['core_oriented_interval']
                    or (fine,n-fine-free,free)!=(original['fine_positions'],original['coarse_positions'],original['fully_unknown_positions'])
                    or context[offset:offset+n]!=re.sub('[^ACGT]','N',core)):
                raise ValueError('Original view or observations changed')
            mask=np.ones((len(context),15),dtype=np.bool_); mask[offset:offset+n]=allowed.numpy()
            if w['id'] in reused:
                row=dict(reused[w['id']]); pos=row['first_blocked_oriented_base']
                if row['genomic_context']!=[a,b] or row['core_oriented_interval']!=[offset,offset+n]:
                    raise ValueError('Reused exact trace geometry changed')
                if pos is not None:
                    x,y=row['nearby_oriented_interval']
                    if context[x:y]!=row['nearby_oriented_DNA'] or np.where(mask[pos])[0].tolist()!=row['allowed_observed_roles']:
                        raise ValueError('Reused localized DNA or original observations changed')
                    checks=primary_checks(context,a,b,w,row['primary_transcripts_at_base'],catalogs[w['species']])
                    if checks!=row['primary_DNA_checks_at_blocked_base']:
                        raise ValueError('Reused actual primary evidence changed')
                row.update(localization_evidence='REUSED_TRACE',source_localization_job='13429503')
            else:
                result=T.trace(context,mask); pos=result['first_blocked_oriented_base']; checks=[]
                if pos is not None:
                    x=a+pos if w['strand']=='+' else b-1-pos; inside=offset<=pos<offset+n
                    tids=[tid for tid in original['original_primary_ids']
                        if catalogs[w['species']][tid]['CDS'][0][0]<=x<catalogs[w['species']][tid]['CDS'][-1][1]]
                    result.update(genomic_base_zero_based=x,first_blocked_base_in_original_core=inside,
                        original_fine_label=int(targets['seg'][pos-offset]) if inside else None,
                        original_coarse_label=int(targets['region'][pos-offset]) if inside else None,
                        primary_transcripts_at_base=tids)
                    checks=primary_checks(context,a,b,w,tids,catalogs[w['species']])
                row=dict(original_ordinal=original['ordinal'],window_id=w['id'],species=w['species'],strand=w['strand'],
                    genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
                    original_primary_ids=original['original_primary_ids'],outside_primary_ids=original['outside_primary_ids'],
                    all_core_observation_rows_preserved=True,primary_DNA_checks_at_blocked_base=checks,**result,
                    seconds=time.perf_counter()-t,localization_evidence='NEW_TRACE',source_localization_job=os.environ['SLURM_JOB_ID'])
            row['this_run_seconds']=time.perf_counter()-t
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
    if opened!={index[r['window_id']]['artifact'] for r in failures}:
        raise ValueError('Only all 24 original counterexample core caches required')
    unique={}
    for r in records:
        for c in r['primary_DNA_checks_at_blocked_base']:
            stable={k:v for k,v in c.items() if k!='oriented_CDS_parts'}
            if c['transcript_id'] in unique and unique[c['transcript_id']]!=stable:
                raise ValueError('The same original transcript has inconsistent direct DNA evidence')
            unique[c['transcript_id']]=stable
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        source_jobs=['13429444','13429503'],records=records,all_original_failures_localized=True,
        actual_counterexamples=24,new_localizations=20,reused_localizations=4,
        mechanism_counts=dict(Counter(r['mechanism'] for r in records)),
        first_blocked_fine_labels=dict(Counter(str(r.get('original_fine_label')) for r in records)),
        first_blocked_coarse_labels=dict(Counter(str(r.get('original_coarse_label')) for r in records)),
        first_blocked_in_core=sum(r.get('first_blocked_base_in_original_core',False) for r in records),
        unique_checked_primary_transcripts=list(unique.values()),
        unique_primary_DNA_failure_counts=dict(Counter(f for c in unique.values() for f in c.get('direct_DNA_failures',[]))),
        unique_primary_geometry_failure_counts=dict(Counter(f for c in unique.values() for f in c.get('geometry_failures',[]))),
        original_TRAIN_DNA_sources=descriptors,cached_artifacts_opened=sorted(opened),
        seconds_after_module_import=time.perf_counter()-began,
        interpretation='all_original_counterexample_localizations_not_annotation_error_or_model_performance',**S.FLAGS)
    S.write_json(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
