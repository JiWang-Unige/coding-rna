#!/usr/bin/env python3
"""Inspect the actual interrupted TRAIN window; no grammar or label changes."""
import json
import os
import re
import time
import torch
import trace_full_grammar_train_support_failures as T
from src.m28.candidates import ORFIndex, STOP
from src.m28.full_grammar_context_view import context_geometry
from src.m28.labels import oriented_parts

S,Q,W,P,ROOT=T.S,T.Q,T.W,T.P,T.ROOT
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-CDS-CHECK-DISAGREEMENT'
FAILED=ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-CENSUS'
CONFIG=dict(scope='first_actual_unwritten_counterexample_of_failed_13429578',
    diagnostic='independent_concatenated_CDS_versus_actual_ORFIndex_each_exon',
    grammar_changed=False,labels_changed=False,optimizer_steps=0)


def main():
    began=time.perf_counter()
    if ((FAILED/'JOBID').read_text().strip()!='13429578'
            or (FAILED/'STATUS').read_text().strip()!='FAILED exit=1'):
        raise ValueError('Actual failed source required')
    source=json.loads((S.OUT/'summary.json').read_text())
    partial=[json.loads(x) for x in (FAILED/'windows.jsonl').read_text().splitlines()]
    bad=[r for r in source['records'] if not r['support_feasible']]
    if len(partial)!=21 or [r['window_id'] for r in partial]!=[r['window_id'] for r in bad[:21]]:
        raise ValueError('Actual unchanged 21-record prefix required')
    original=bad[21]
    S.write_json(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        window_id=original['window_id'],source_jobs=['13429444','13429578']))
    _,_,windows,catalogs=Q.load_train_inputs(); sequences,descriptors=P.read_train_DNA()
    w=windows[original['window_id']]; chromosome=sequences[w['species']][w['seqid']]
    geometry=context_geometry(w['start'],w['end'],len(chromosome)); a,b=geometry['genomic_context']
    if [a,b]!=original['genomic_context']: raise ValueError('Original fixed view changed')
    context=re.sub('[^ACGT]','N',chromosome[a:b])
    if w['strand']=='-': context=Q.M.rc(context)
    records=[]
    for tid in original['original_primary_ids']:
        tx=catalogs[w['species']][tid]; parts=oriented_parts(tx,dict(start=a,end=b,strand=w['strand']))
        coding=''.join(context[x:y] for x,y in parts)
        direct=[]
        if len(coding)%3: direct.append('length_not_mod3')
        if coding[:3]!='ATG': direct.append('non_ATG_start')
        if coding[-3:] not in STOP: direct.append('no_terminal_stop')
        if any(coding[i:i+3] in STOP for i in range(0,max(0,len(coding)-3),3)): direct.append('internal_stop')
        if set(coding)-set('ACGT'): direct.append('ambiguous_CDS')
        orf=ORFIndex(context); carry=''; actual=coding[:3]=='ATG'; steps=[]
        for i,(x,y) in enumerate(parts):
            next_carry=orf.extend(x,y,carry,i==len(parts)-1)
            steps.append(dict(exon=i,interval=[x,y],length=y-x,incoming_carry=carry,
                outgoing_carry=next_carry,head=context[x:min(x+3,y)],tail=context[max(x,y-3):y],
                terminal=i==len(parts)-1))
            if next_carry is None: actual=False; break
            carry=next_carry
        error=None
        try: old_failures,old_actual=W.DNA_checks(context,parts)
        except ValueError as exc:
            if str(exc)!='Actual-carry and direct concatenated-CDS grammar verdicts disagree': raise
            error=str(exc); old_failures=None; old_actual=None
        records.append(dict(transcript_id=tid,CDS_complete_supported=tx['CDS_complete_supported'],
            eligible_positive=tx['eligible_positive'],oriented_CDS_parts=[list(p) for p in parts],
            geometry_failures=W.geometry_failures(parts,len(context)),
            concatenated_CDS_length=len(coding),concatenated_CDS_start=coding[:3],
            concatenated_CDS_stop=coding[-3:],direct_coding_failures=direct,
            actual_carry_coding_legal=actual,carry_steps=steps,
            original_helper_error=error,original_helper_failures=old_failures,original_helper_actual=old_actual))
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        source_jobs=['13429444','13429578'],original_ordinal=original['ordinal'],
        window_id=original['window_id'],records=records,original_TRAIN_DNA_sources=descriptors,
        seconds=time.perf_counter()-began,interpretation='diagnostic_of_actual_helper_disagreement_not_repair',**S.FLAGS)
    S.write_json(OUT/'summary.json',summary); print(json.dumps(summary,allow_nan=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(2); main()
