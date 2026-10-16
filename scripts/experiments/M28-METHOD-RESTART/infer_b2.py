#!/usr/bin/env python3
"""Full DEV paired B2 inference. No reference annotation or TRAIN target read."""
import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
from infer_pilot import DEV_LENGTHS,dev_groups,read_rows
from fit_b2 import BASE,CONTRACT,check_output_cap
from src.m28.core import M28Core
from src.m28.candidates import generate,Budget
from src.m28.labels import CDS_COLUMNS
from src.m28.training import one_hot
from src.m28.prefix_search import generate_prefix
from src.m28.prefix_value import PrefixValue
from src.m28.prefix_inference import chain_key,score_paired_pools,validate_additive_replay

def main(fit,out):
    fit=fit.resolve();out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    fitted=json.loads((fit/'summary.json').read_text())
    if not fitted['completed'] or fitted['optimizer_steps']!=4608 or fitted['feature_revision']!=F.REVISION:
        raise ValueError('Only completed B2 last checkpoint may be evaluated')
    checkpoint=M.C.ROOT/fitted['final_checkpoint']
    state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    if (state['arm']!='B2' or state['step']!=4608 or state['seed']!=0 or state['revision']!=F.REVISION or
        state['base_checkpoint']!=BASE):
        raise ValueError('B2 checkpoint identity mismatch')
    torch.set_num_threads(2);model=PrefixValue(128)
    model.load_state_dict(state['prefix'],strict=True);del state
    model.eval().requires_grad_(False).to('cuda')
    state=torch.load(M.C.ROOT/BASE,map_location='cpu',weights_only=True)
    if state['arm']!='B1' or state['step']!=4608 or state['revision']!=F.REVISION:
        raise ValueError('Frozen B1 source mismatch')
    base=M28Core('B1',2048,128);base.load_state_dict(state['model'],strict=True);del state
    base.eval().requires_grad_(False).to('cuda')
    pilot=M.C.ROOT/'outputs/M28-PILOT-R5'
    groups=dev_groups(read_rows(pilot/'features/index.jsonl'))
    old=json.loads((pilot/'infer_B1/summary.json').read_text())
    saved_scopes={(s['species'],s['seqid'],s['strand']):s for s in old['scopes']}
    if (not old['inference_complete'] or old['windows']!=8690 or old['checkpoint']!=BASE or
        set(saved_scopes)!=set(groups) or len(old['scopes'])!=4):
        raise ValueError('Incomplete frozen R5 replay source')
    torch.cuda.reset_peak_memory_stats()
    began=time.perf_counter();done=0;common=0;max_error=0.;arms={'early':{'scopes':[]},'late':{'scopes':[]}}
    with torch.inference_mode():
        for scope_number,(key,rows) in enumerate(groups.items()):
            species,seqid,strand=key
            filenames={arm:arm+'_scope_'+str(scope_number).zfill(2)+'.jsonl' for arm in arms}
            windows=[];candidates=defaultdict(int);timing=defaultdict(float)
            source=pilot/'infer_B1'/saved_scopes[key]['file']
            with source.open() as replay, (out/filenames['early']).open('x') as early_file, (out/filenames['late']).open('x') as late_file:
                handles={'early':early_file,'late':late_file}
                for i,r in enumerate(rows):
                    t=time.perf_counter()
                    saved_line=replay.readline()
                    if not saved_line: raise ValueError('Partial R5 replay file')
                    saved=json.loads(saved_line)
                    if saved['window_id']!=r['window_id'] or saved['window_index']!=i:
                        raise ValueError('R5 replay window identity mismatch')
                    cached=torch.load(M.C.ROOT/r['artifact'],map_location='cpu',weights_only=True)
                    w=cached['window'];n=w['valid_bases'];sequence=cached['sequence']
                    start=int(r['window_id'].split('|')[2])
                    if (cached['revision']!=F.REVISION or w['id']!=r['window_id'] or
                        (w['species'],w['seqid'],w['strand'],w['start'],w['end'])!=
                        (species,seqid,strand,start,min(start+24576,DEV_LENGTHS[species][seqid])) or
                        len(sequence)!=n or tuple(cached['features'].shape)!=(4096,2048)):
                        raise ValueError('DEV feature/geometry mismatch')
                    windows.append([w['start'],w['end']])
                    features=cached['features'].float()[None].to('cuda')
                    dna=one_hot(sequence,w['width'])[None].to('cuda')
                    torch.cuda.synchronize();timing['load_and_saved_control_seconds']+=time.perf_counter()-t
                    t=time.perf_counter();outputs=base(features,dna,[n]);h=outputs['features'][0,:n]
                    p=outputs['segmentation_logits'][0,:n].softmax(-1)
                    ep=outputs['endpoint_logits'][0,:n].cpu().numpy()
                    cp=p[:,list(CDS_COLUMNS)].sum(-1).cpu().numpy()
                    torch.cuda.synchronize();timing['frozen_forward_seconds']+=time.perf_counter()-t
                    links=lambda pairs:base.chain.link_logits(h,pairs).cpu().numpy()
                    t=time.perf_counter()
                    late=generate(sequence,ep,cp,links,Budget())
                    error=validate_additive_replay(late,saved);max_error=max(max_error,error)
                    torch.cuda.synchronize();timing['late_control_and_replay_seconds']+=time.perf_counter()-t
                    t=time.perf_counter();early_scorer=model.for_window(h)
                    early=generate_prefix(sequence,ep,cp,links,Budget(),rank_states_fn=early_scorer.priorities)
                    del early_scorer
                    torch.cuda.synchronize();timing['early_search_seconds']+=time.perf_counter()-t
                    t=time.perf_counter();scores=score_paired_pools(model,h,early['chains'],late['chains'])
                    shared=set(map(chain_key,early['chains']))&set(map(chain_key,late['chains']))
                    a=dict(zip(map(chain_key,early['chains']),scores['early']))
                    b=dict(zip(map(chain_key,late['chains']),scores['late']))
                    if any(a[c]!=b[c] for c in shared): raise ValueError('Paired final common gain mismatch')
                    common+=len(shared)
                    torch.cuda.synchronize();timing['shared_union_readout_seconds']+=time.perf_counter()-t
                    t=time.perf_counter()
                    for arm,free in (('early',early),('late',late)):
                        record=dict(free,score=scores[arm],window_index=i,window_id=w['id'])
                        handles[arm].write(json.dumps(record,allow_nan=False)+'\n');candidates[arm]+=len(free['chains'])
                    done+=1;timing['serialization_seconds']+=time.perf_counter()-t
                    if done%32==0:
                        for f in handles.values(): f.flush()
                        check_output_cap()
                        print(json.dumps({'windows':done,'total':8690,'seconds':time.perf_counter()-began,
                                          'R5_proposal_max_abs':max_error,'common_chain_records':common}),flush=True)
                    del cached,features,dna,outputs,h,p,early,late,saved,scores
                if replay.readline(): raise ValueError('Extra R5 replay records')
            if windows!=saved_scopes[key]['windows']: raise ValueError('Frozen scope grid changed')
            for arm in arms:
                arms[arm]['scopes'].append({'species':species,'seqid':seqid,'strand':strand,'windows':windows,
                                            'file':filenames[arm],'candidate_records':candidates[arm]})
            print(json.dumps({'completed_scope':list(key),'candidates':dict(candidates),'seconds':dict(timing)}),flush=True)
            for arm in arms: arms[arm]['scopes'][-1]['shared_pair_seconds']=dict(timing)
    if done!=8690: raise ValueError('Incomplete paired DEV')
    report={'experiment':'M28-B2-PAIRED-INFER','inference_complete':True,'windows':done,'optimizer_steps':4608,
            'checkpoint':str(checkpoint.relative_to(M.C.ROOT)),'base_checkpoint':BASE,'feature_revision':F.REVISION,
            'contract':CONTRACT,'reference_used':False,'test_setaria_used':False,'arms':arms,
            'control_replay':{'windows':done,'discrete_exact':True,'proposal_max_abs':max_error,'tolerance':1e-5},
            'paired_common_gain':{'chain_records':common,'exact_equality':True,'single_union_readout':True},
            'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
            'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),'artifact_bytes_before_summary':check_output_cap()}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    check_output_cap();print(json.dumps({k:v for k,v in report.items() if k!='arms'},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fit-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.fit_dir,a.output_dir)
