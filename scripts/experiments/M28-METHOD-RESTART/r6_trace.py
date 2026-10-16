#!/usr/bin/env python3
"""R6 one fixed DEV subset; no fitting, re-ranking, GFF or feature generation."""
import argparse
import json
import sys
import time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import torch
import feature_smoke as F
import manifest as M
from infer_pilot import DEV_LENGTHS,read_rows
from src.m28.core import M28Core
from src.m28.candidates import generate,Budget
from src.m28.labels import CDS_COLUMNS
from src.m28.training import one_hot
from r6_shadow import ReferenceTrace,generate_shadow

def normalized(value):
    return json.loads(json.dumps(value))

def compare(left,right,tolerance):
    for key in ('chains','links','counts','budget','reference_used'):
        if normalized(left[key])!=normalized(right[key]):
            raise ValueError('Replay mismatch in '+key)
    a=np.asarray(left['proposal_scores']);b=np.asarray(right['proposal_scores'])
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Nonfinite or misaligned proposal replay')
    error=float(np.max(np.abs(a-b))) if len(a) else 0.
    if error>tolerance: raise ValueError('Proposal replay error '+str(error))
    return error

def main(sample_path,out):
    sample=json.loads(sample_path.read_text())
    if sample['anchor_count']!=20 or sample['window_count']!=23 or sample['missing_strata']:
        raise ValueError('Frozen sample identity mismatch')
    anchors=sample['anchors'];windows=sample['windows']
    if len(anchors)!=20 or len(windows)!=23 or len({w['window_id'] for w in windows})!=23:
        raise ValueError('Frozen sample multiplicity mismatch')
    if out.exists(): raise FileExistsError(out)
    out.mkdir(parents=True)
    began=time.perf_counter()
    wanted={w['window_id'] for w in windows}
    pilot=M.C.ROOT/'outputs/M28-PILOT-R5'
    index={r['window_id']:r for r in read_rows(pilot/'features/index.jsonl')
           if r['split']=='val' and r['window_id'] in wanted}
    if set(index)!=wanted: raise ValueError('Sample cache index mismatch')
    saved={}
    for filename in sorted({w['scope_file'] for w in windows}):
        with (pilot/'infer_B1'/filename).open() as f:
            for line in f:
                record=json.loads(line)
                if record['window_id'] in wanted:
                    if record['window_id'] in saved: raise ValueError('Duplicate R5 replay row')
                    saved[record['window_id']]=record
    if set(saved)!=wanted: raise ValueError('Missing frozen R5 replay row')
    fit=json.loads((pilot/'B1/summary.json').read_text())
    if fit['optimizer_steps']!=4608 or fit['feature_revision']!=F.REVISION:
        raise ValueError('Wrong completed fit')
    checkpoint=M.C.ROOT/fit['final_checkpoint']
    if str(checkpoint.relative_to(M.C.ROOT))!='outputs/M28-PILOT-R5/B1/step_004608.pt':
        raise ValueError('Checkpoint path changed')
    state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    if state['arm']!='B1' or state['step']!=4608 or state['revision']!=F.REVISION:
        raise ValueError('Wrong frozen checkpoint')
    torch.set_num_threads(2)
    model=M28Core('B1',2048,128)
    model.load_state_dict(state['model'],strict=True);del state
    model.eval().requires_grad_(False).to('cuda')
    torch.cuda.reset_peak_memory_stats()
    cases=[];replay=[]
    with torch.inference_mode(),(out/'window_traces.jsonl').open('x') as handle:
        for w in windows:
            ident=w['window_id'];r=index[ident]
            cached=torch.load(M.C.ROOT/r['artifact'],map_location='cpu',weights_only=True)
            cw=cached['window'];n=cw['valid_bases'];sequence=cached['sequence']
            expected=(w['species'],w['seqid'],w['strand'],w['start'],w['end'])
            observed=(cw['species'],cw['seqid'],cw['strand'],cw['start'],cw['end'])
            if (observed!=expected or cw['id']!=ident or cached['revision']!=F.REVISION or
                len(sequence)!=n or tuple(cached['features'].shape)!=(4096,2048) or
                w['seqid'] not in DEV_LENGTHS[w['species']]):
                raise ValueError('DEV feature geometry mismatch')
            watched=[]
            for anchor in anchors:
                if ident not in anchor['window_ids']: continue
                chain=anchor['chain']
                oriented=[(a-w['start'],b-w['start']) for a,b in chain] if w['strand']=='+' else [
                    (w['end']-b,w['end']-a) for a,b in reversed(chain)]
                if not all(0<=a<b<=n for a,b in oriented): raise ValueError('Anchor outside window')
                watched.append(dict(anchor,oriented_chain=oriented))
            if not watched: raise ValueError('Unassigned sample window')
            outputs=model(cached['features'].float()[None].to('cuda'),
                          one_hot(sequence,cw['width'])[None].to('cuda'),[n])
            h=outputs['features'][0,:n]
            p=outputs['segmentation_logits'][0,:n].softmax(-1)
            ep=outputs['endpoint_logits'][0,:n].cpu().numpy()
            cp=p[:,list(CDS_COLUMNS)].sum(-1).cpu().numpy()
            score_links=lambda pairs:model.chain.link_logits(h,pairs).cpu().numpy()
            original=generate(sequence,ep,cp,link_score_fn=score_links,budget=Budget())
            observer=ReferenceTrace(watched)
            traced=generate_shadow(sequence,ep,cp,link_score_fn=score_links,budget=Budget(),observer=observer)
            shadow_error=compare(original,traced,0.)
            saved_error=compare(original,saved[ident],1e-5)
            traced_cases=observer.results()
            for case in traced_cases:
                case.update(window_id=ident,window_index=w['window_index'],
                            is_owner=w['window_index']==next(a['owner_index'] for a in watched if a['reference_id']==case['reference_id']))
                local_key=tuple(tuple(x) for x in case['oriented_chain'])
                seen=local_key in {tuple(tuple(x) for x in c) for c in saved[ident]['chains']}
                if seen!=case['final']['survived']: raise ValueError('Trace final presence mismatch')
                if seen:
                    k=next(i for i,c in enumerate(saved[ident]['chains']) if tuple(tuple(x) for x in c)==local_key)
                    case['frozen_R5_chain_gain']=saved[ident]['score'][k]
                cases.append(case)
            item={'window_id':ident,'shadow_exact':True,'saved_proposal_max_abs':saved_error,
                  'case_traces':traced_cases}
            handle.write(json.dumps(item,allow_nan=False)+'\n');handle.flush()
            replay.append({'window_id':ident,'shadow_max_abs':shadow_error,'saved_max_abs':saved_error})
            if (out/'window_traces.jsonl').stat().st_size>32*1024**2: raise RuntimeError('Trace output cap')
            print(json.dumps({'window':len(replay),'total':23,'seconds':time.perf_counter()-began}),flush=True)
            del cached,outputs,h,p,original,traced,observer
    positive_seen=defaultdict(bool)
    for c in cases:
        if c['stratum']=='selected_exact' and c['is_owner'] and not c['final']['survived']:
            raise ValueError('Positive owner control lost')
        if c['stratum']=='nonowner_positive' and not c['is_owner'] and c.get('frozen_R5_chain_gain',-1)>0:
            positive_seen[c['reference_id']]=True
    nonowners={a['reference_id'] for a in anchors if a['stratum']=='nonowner_positive'}
    if set(positive_seen)!=nonowners: raise ValueError('Nonowner positive witness missing')
    species_witnesses={s:sorted({c['reference_id'] for c in cases if c['species']==s and c['mechanism_witness']})
                       for s in DEV_LENGTHS}
    decision='prepare_B2_prefix_search_contract' if all(species_witnesses.values()) else 'do_not_expand_trace_no_two_species_witness'
    report={'experiment':'M28-R6-CANDIDATE-TRACE','completed':True,'fitting':False,'new_forward_pass':True,
            'windows':23,'anchors':20,'case_window_pairs':len(cases),'checkpoint':str(checkpoint.relative_to(M.C.ROOT)),
            'feature_revision':F.REVISION,'reference_used_for_observation_only':True,
            'reference_used_by_candidate_ranking':False,'test_setaria_used':False,
            'sample_manifest':str(sample_path.relative_to(M.C.ROOT)),'replay':replay,
            'species_mechanism_witnesses':species_witnesses,'decision':decision,
            'first_cut_counts_case_window_pairs':dict(Counter(c['first_actual_cut']['stage'] if c['first_actual_cut'] else 'survived_final' for c in cases)),
            'wall_seconds':time.perf_counter()-began,
            'peak_gpu_allocated_bytes':torch.cuda.max_memory_allocated(),
            'peak_gpu_reserved_bytes':torch.cuda.max_memory_reserved(),
            'small_sample_no_population_attribution':True,'no_new_performance_ranking':True}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    if sum(p.stat().st_size for p in out.iterdir() if p.is_file())>32*1024**2: raise RuntimeError('Total output cap')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sample',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.sample.resolve(),a.output_dir.resolve())
