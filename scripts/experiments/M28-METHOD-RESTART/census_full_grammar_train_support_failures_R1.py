#!/usr/bin/env python3
"""One same-condition census repair: report split terminal STOP explicitly."""
import json
import torch
import census_full_grammar_train_support_failures as C
import check_full_grammar_CDS_disagreement as D
from src.m28.candidates import ORFIndex, STOP


def DNA_checks(sequence,parts):
    coding=''.join(sequence[x:y] for x,y in parts); failures=[]
    if len(coding)%3: failures.append('length_not_mod3')
    if coding[:3]!='ATG': failures.append('non_ATG_start')
    if coding[-3:] not in STOP: failures.append('no_terminal_stop')
    if any(coding[i:i+3] in STOP for i in range(0,max(0,len(coding)-3),3)): failures.append('internal_stop')
    if set(coding)-set('ACGT'): failures.append('ambiguous_CDS')
    # ORFIndex and the unchanged full grammar require the terminal STOP
    # within the last exon. Concatenation alone missed the observed 2bp exon.
    if coding[-3:] in STOP and parts[-1][1]-parts[-1][0]<3:
        failures.append('split_terminal_stop')
    orf=ORFIndex(sequence); carry=''; actual=coding[:3]=='ATG'
    for i,(x,y) in enumerate(parts):
        carry=orf.extend(x,y,carry,i==len(parts)-1)
        if carry is None: actual=False; break
    if actual!=(not failures):
        raise ValueError('Actual-carry and direct concatenated-CDS grammar verdicts disagree')
    for (_,y),(z,_) in zip(parts,parts[1:]):
        if sequence[y:y+2] not in ('GT','GC') or sequence[z-2:z]!='AG': failures.append('nonstandard_splice')
    return failures,actual


def primary_checks(context,a,b,w,ids,catalog):
    checks=[]
    for tid in ids:
        tx=catalog[tid]; parts=C.oriented_parts(tx,dict(start=a,end=b,strand=w['strand']))
        if not 0<=parts[0][0]<parts[-1][1]<=len(context):
            checks.append(dict(transcript_id=tid,DNA_checked=False,reason='outside_fixed_view')); continue
        geometry=C.W.geometry_failures(parts,len(context)); failures,coding=DNA_checks(context,parts)
        checks.append(dict(transcript_id=tid,gene_id=tx['gene_id'],CDS_complete_supported=tx['CDS_complete_supported'],
            eligible_positive=tx['eligible_positive'],oriented_CDS_parts=[list(p) for p in parts],DNA_checked=True,
            geometry_failures=geometry,direct_DNA_failures=failures,actual_carry_coding_legal=coding,
            complete_original_chain_legal=not geometry and not failures))
    return checks


if __name__=='__main__':
    evidence=json.loads((D.OUT/'summary.json').read_text())
    if ((D.OUT/'STATUS').read_text().strip()!='COMPLETED' or (D.OUT/'JOBID').read_text().strip()!='13429614'
            or evidence['job_id']!='13429614' or not evidence['completed'] or evidence['configuration']!=D.CONFIG):
        raise ValueError('Actual diagnostic source required for the one bounded repair')
    failure=[r for r in evidence['records'] if r['original_helper_error'] is not None]
    if (len(failure)!=1 or failure[0]['transcript_id']!='rna-XM_015764666.3'
            or failure[0]['direct_coding_failures'] or failure[0]['actual_carry_coding_legal']
            or failure[0]['geometry_failures']!=['split_terminal_codon']
            or failure[0]['carry_steps'][-1]['length']!=2
            or failure[0]['carry_steps'][-1]['incoming_carry']!='T'
            or failure[0]['carry_steps'][-1]['head']!='GA'):
        raise ValueError('Exact actual split-STOP implementation counterexample required')
    C.OUT=C.ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-CENSUS-R1'
    C.CONFIG=dict(C.CONFIG,implementation_repair='include_existing_within_last_exon_STOP_condition_in_direct_verdict',
        failed_job_retained='13429578',actual_disagreement_job='13429614',scientific_configuration_changed=False)
    C.primary_checks=primary_checks
    torch.set_num_threads(2); C.main()
