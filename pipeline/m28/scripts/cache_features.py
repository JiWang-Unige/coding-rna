#!/usr/bin/env python3
"""One fixed R4 feature cache for both arms; train draws plus the complete DEV grid."""
import argparse
import json
import re
import time
from pathlib import Path
import torch
from transformers import AutoModelForTokenClassification,AutoTokenizer
import feature_smoke as F
import manifest as M

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def main(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'index.jsonl').exists(): raise FileExistsError('Do not overwrite cache index')
    base=M.C.ROOT/'outputs/M28-LABEL-POLICY-R4'
    train={};val=[]
    for species in ('arabidopsis_thaliana','oryza_sativa'):
        train.update({w['id']:w for w in rows(base/species/'windows_train.jsonl')})
        val+=rows(base/species/'windows_val.jsonl')
    draws=rows(M.C.ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl')
    assert len(draws)==1536 and len(val)==8690
    selected={d['window_id']:train[d['window_id']] for d in draws}
    ntrain=len(selected);selected.update({w['id']:w for w in val})
    meta=M.C.metadata(M.C.yaml.safe_load(M.C.CONFIG.read_text()))
    sequences={s:M.read_allowed(r['path']/'genome.fa',r['lengths']) for s,r in meta.items()}
    oldroot=M.C.ROOT/'outputs/M28-GLM-FEATURE-SMOKE-R3/repair1'
    old=json.loads((oldroot/'summary.json').read_text())
    assert old['revision']==F.REVISION
    reuse={r['window_id']:r for r in old['records']}
    tokenizer=AutoTokenizer.from_pretrained(F.MODEL,revision=F.REVISION,trust_remote_code=True,local_files_only=True)
    full=AutoModelForTokenClassification.from_pretrained(F.MODEL,revision=F.REVISION,trust_remote_code=True,local_files_only=True,
                                                       attn_implementation='sdpa',torch_dtype=torch.bfloat16)
    backbone=full.model;del full
    backbone.eval().requires_grad_(False).to('cuda')
    torch.cuda.reset_peak_memory_stats()
    began=time.perf_counter();stored=0;new_bytes=0;reused=0;forward_seconds=0.
    with (out/'index.jsonl').open('x') as index:
        for ordinal,w in enumerate(selected.values()):
            if w['id'] in reuse:
                artifact=oldroot/reuse[w['id']]['artifact'];reused+=1
            else:
                raw=sequences[w['species']][w['seqid']][w['start']:w['end']]
                sequence=raw if w['strand']=='+' else M.rc(raw)
                padded=re.sub('[^ACGT]','N',sequence)+'N'*w['pad_right_after_orientation']
                whole=tokenizer(padded,add_special_tokens=False,return_tensors='pt',truncation=False)['input_ids']
                blocks=[tokenizer(padded[a:a+F.BLOCK],add_special_tokens=False,return_tensors='pt',truncation=False)['input_ids']
                        for a in range(0,len(padded),F.BLOCK)]
                if any(tuple(x.shape)!=(1,1024) for x in blocks) or not torch.equal(torch.cat(blocks,dim=1),whole):
                    raise ValueError('Tokenizer alignment changed')
                features=[]
                for ids in blocks:
                    torch.cuda.synchronize();start=time.perf_counter()
                    with torch.inference_mode():
                        h=backbone(input_ids=ids.to('cuda'),attention_mask=torch.ones_like(ids,device='cuda')).last_hidden_state
                    torch.cuda.synchronize();forward_seconds+=time.perf_counter()-start
                    if h.shape!=(1,1024,2048) or not torch.isfinite(h).all(): raise ValueError('Invalid frozen features')
                    features.append(h[0].detach().to('cpu',dtype=torch.bfloat16))
                artifact=out/('window_'+str(ordinal).zfill(5)+'.pt')
                geometry={k:w[k] for k in ('id','species','split','seqid','start','end','width','valid_bases','strand','pad_right_after_orientation')}
                torch.save({'features':torch.cat(features,dim=0).contiguous(),'sequence':sequence,
                            'window':geometry,'model':F.MODEL,'revision':F.REVISION},artifact)
                new_bytes+=artifact.stat().st_size
            stored+=1
            if new_bytes>180*1024**3: raise RuntimeError('R4 feature cache exceeds frozen 180 GiB cap')
            index.write(json.dumps({'window_id':w['id'],'split':w['split'],'species':w['species'],
                                    'artifact':str(artifact.relative_to(M.C.ROOT)),'bytes':artifact.stat().st_size,
                                    'reused_R3':w['id'] in reuse})+'\n')
            if stored%64==0:
                index.flush();print(json.dumps({'completed':stored,'total':len(selected),'new_bytes':new_bytes,
                                               'seconds':time.perf_counter()-began}),flush=True)
    report={'revision':F.REVISION,'model':F.MODEL,'train_unique_windows':ntrain,'DEV_windows':len(val),
            'total_windows':stored,'reused_R3_windows':reused,'new_cache_bytes':new_bytes,
            'forward_seconds':forward_seconds,'extraction_loop_seconds':time.perf_counter()-began,
            'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
            'optimizer_steps':0,'backbone_trainable_parameters':sum(p.numel() for p in backbone.parameters() if p.requires_grad),
            'test_setaria_used':False,'feature_only_no_reference_inputs':True}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();main(a.output_dir)
