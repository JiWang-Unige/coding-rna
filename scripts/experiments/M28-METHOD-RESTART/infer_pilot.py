#!/usr/bin/env python3
"""Full frozen DEV inference. No reference annotations, labels or injected candidates."""
import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
from src.m28.core import M28Core
from src.m28.candidates import generate,Budget
from src.m28.c0_decode import decode_c0
from src.m28.labels import CDS_COLUMNS
from src.m28.training import one_hot

DEV_LENGTHS={'arabidopsis_thaliana':{'NC_003074.8':23459830},
             'oryza_sativa':{'NC_089041.1':29936421}}

def read_rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def dev_groups(index):
    actual=[r for r in index if r['split']=='val']
    by_id={r['window_id']:r for r in actual}
    expected={}
    for species,lengths in DEV_LENGTHS.items():
        for seqid,length in lengths.items():
            for strand in ('+','-'):
                expected[species,seqid,strand]=[
                    species+'|'+seqid+'|'+str(a)+'|'+strand for a in M.starts(length)]
    ids={ident for group in expected.values() for ident in group}
    if len(actual)!=8690 or len(actual)!=len(by_id) or set(by_id)!=ids:
        raise ValueError('Incomplete or mismatched frozen DEV window grid')
    return {key:[by_id[ident] for ident in group] for key,group in expected.items()}

def main(arm,pilot,out):
    pilot=pilot.resolve();out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists(): raise FileExistsError('Completed inference exists')
    fitted=json.loads((pilot/arm/'summary.json').read_text())
    if fitted['optimizer_steps']!=4608 or fitted['feature_revision']!=F.REVISION:
        raise ValueError('Only the completed last checkpoint may be evaluated')
    checkpoint=M.C.ROOT/fitted['final_checkpoint']
    state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    if state['arm']!=arm or state['step']!=4608 or state['revision']!=F.REVISION:
        raise ValueError('Checkpoint identity mismatch')
    torch.set_num_threads(2)
    model=M28Core(arm,2048,128)
    model.load_state_dict(state['model'],strict=True);del state
    model.eval().requires_grad_(False).to('cuda')
    groups=dev_groups(read_rows(pilot/'features/index.jsonl'))
    torch.cuda.reset_peak_memory_stats()
    began=time.perf_counter();done=0;scopes=[]
    with torch.inference_mode():
        for scope_number,((species,seqid,strand),rows) in enumerate(groups.items()):
            file=out/('scope_'+str(scope_number).zfill(2)+'.jsonl')
            windows=[];candidate_count=0;timing=defaultdict(float)
            with file.open('x') as handle:
                for i,r in enumerate(rows):
                    t=time.perf_counter()
                    cached=torch.load(M.C.ROOT/r['artifact'],map_location='cpu',weights_only=True)
                    w=cached['window'];n=w['valid_bases'];sequence=cached['sequence']
                    expected_start=int(r['window_id'].split('|')[2])
                    if (cached['revision']!=F.REVISION or w['id']!=r['window_id'] or
                        (w['species'],w['seqid'],w['strand'],w['start'],w['end']) !=
                        (species,seqid,strand,expected_start,min(expected_start+24576,DEV_LENGTHS[species][seqid])) or
                        len(sequence)!=n or tuple(cached['features'].shape)!=(4096,2048)):
                        raise ValueError('DEV feature/geometry mismatch')
                    windows.append([w['start'],w['end']])
                    feat=cached['features'].float()[None].to('cuda')
                    dna=one_hot(sequence,w['width'])[None].to('cuda')
                    torch.cuda.synchronize();timing['load_seconds']+=time.perf_counter()-t
                    t=time.perf_counter()
                    outputs=model(feat,dna,[n])
                    h=outputs['features'][0,:n]
                    p=outputs['segmentation_logits'][0,:n].softmax(-1)
                    torch.cuda.synchronize();timing['shared_forward_seconds']+=time.perf_counter()-t
                    t=time.perf_counter()
                    if arm=='C0':
                        decoded=decode_c0(sequence,p.cpu().numpy())
                        free={'chains':decoded['chains'],'score':decoded['gains'],
                              'counts':{'partial_paths':decoded['partial_paths']},'reference_used':False}
                    else:
                        free=generate(sequence,outputs['endpoint_logits'][0,:n].cpu().numpy(),
                                      p[:,list(CDS_COLUMNS)].sum(-1).cpu().numpy(),
                                      link_score_fn=lambda pairs:model.chain.link_logits(h,pairs).cpu().numpy(),budget=Budget())
                        scores=model.chain(h,free['chains'])
                        free['score']=scores['logits'].cpu().tolist()
                        free['additive_only_score']=scores['additive_only_logits'].cpu().tolist()
                        del scores
                    if any(not torch.isfinite(torch.as_tensor(free[k])).all()
                           for k in ('score','additive_only_score') if k in free):
                        raise ValueError('Nonfinite inference score')
                    torch.cuda.synchronize();timing['candidate_decode_score_seconds']+=time.perf_counter()-t
                    t=time.perf_counter()
                    handle.write(json.dumps(dict(free,window_index=i,window_id=w['id']))+'\n')
                    candidate_count+=len(free['chains']);done+=1
                    timing['serialization_seconds']+=time.perf_counter()-t
                    if done%64==0:
                        handle.flush();print(json.dumps({'arm':arm,'windows':done,'total':8690,
                                                         'seconds':time.perf_counter()-began}),flush=True)
                    del cached,outputs,h,p,free,feat,dna
            scopes.append({'species':species,'seqid':seqid,'strand':strand,'windows':windows,
                           'file':file.name,'candidate_records':candidate_count,'seconds':dict(timing)})
    assert done==8690
    report={'arm':arm,'inference_complete':True,'windows':done,'checkpoint':str(checkpoint.relative_to(M.C.ROOT)),
            'optimizer_steps':4608,'feature_revision':F.REVISION,'reference_used':False,'test_setaria_used':False,
            'scopes':scopes,'wall_loop_seconds':time.perf_counter()-began,
            'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='scopes'},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--arm',choices=('C0','B1'),required=True)
    p.add_argument('--pilot-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.arm,a.pilot_dir,a.output_dir)
