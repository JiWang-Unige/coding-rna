#!/usr/bin/env python3
"""B2 first fit: new prefix head only, 4608 updates on original TRAIN draws."""
import argparse
import json
import time
from collections import Counter
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
from b2_preparation import train_records
from src.m28.core import M28Core
from src.m28.candidates import Budget
from src.m28.labels import build_targets,CDS_COLUMNS
from src.m28.training import one_hot
from src.m28.prefix_search import generate_prefix
from src.m28.prefix_value import PrefixValue
from src.m28.prefix_training import FirstCutObserver,prefix_losses

CONTRACT='reports/M28-METHOD-RESTART/b2_first_fit_contract.md'
BASE='outputs/M28-PILOT-R5/B1/step_004608.pt'
OUTPUTS=('M28-B2-FIRST-FIT','M28-B2-PAIRED-INFER','M28-B2-EVALUATION')

def check_output_cap():
    size=sum(p.stat().st_size for name in OUTPUTS
             for p in (M.C.ROOT/'outputs'/name).rglob('*') if p.is_file())
    if size>2*1024**3: raise RuntimeError('B2 combined artifact cap exceeded')
    return size

def main(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'training.jsonl').exists(): raise FileExistsError('B2 fit cannot resume or overwrite')
    with (M.C.ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl').open() as f:
        draws=[json.loads(line) for line in f if line.strip()]
    ids={d['window_id'] for d in draws}
    if (len(draws)!=1536 or len(ids)!=1491 or
        Counter(d['species'] for d in draws)!=Counter({'arabidopsis_thaliana':768,'oryza_sativa':768})):
        raise ValueError('Original TRAIN draw contract mismatch')
    index={r['window_id']:r for r in train_records(M.C.ROOT/'outputs/M28-PILOT-R5/features/index.jsonl')
           if r['window_id'] in ids}
    label_root=M.C.ROOT/'outputs/M28-LABEL-POLICY-R4'
    windows={};catalogs={}
    for species in ('arabidopsis_thaliana','oryza_sativa'):
        windows.update({w['id']:w for w in train_records(label_root/species/'windows_train.jsonl') if w['id'] in ids})
        needed={ident for w in windows.values() if w['species']==species for ident in w['overlapping_transcript_ids']}
        catalogs[species]={t['id']:t for t in train_records(label_root/species/'transcripts.jsonl') if t['id'] in needed}
        if set(catalogs[species])!=needed: raise ValueError('TRAIN target identity mismatch')
    if set(index)!=ids or set(windows)!=ids: raise ValueError('TRAIN window/index mismatch')
    for d in draws:
        w=windows[d['window_id']]
        if ((w['species'],w['seqid'],w['strand'])!=(d['species'],*d['stratum']) or
            not 0<float(d['importance_weight_uniform_within_stratum'])<float('inf')):
            raise ValueError('TRAIN stratum/importance mismatch')
        if any(t['seqid']!=w['seqid'] or t['strand']!=w['strand']
               for t in (catalogs[w['species']][i] for i in w['overlapping_transcript_ids'])):
            raise ValueError('TRAIN transcript/window mismatch')
    state=torch.load(M.C.ROOT/BASE,map_location='cpu',weights_only=True)
    if state['arm']!='B1' or state['step']!=4608 or state['revision']!=F.REVISION:
        raise ValueError('Frozen B1 base mismatch')
    torch.set_num_threads(2)
    base=M28Core('B1',2048,128);base.load_state_dict(state['model'],strict=True);del state
    base.eval().requires_grad_(False).to('cuda')
    torch.manual_seed(0);model=PrefixValue(128).to('cuda').train()
    optimizer=torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.01,betas=(.9,.999),eps=1e-8)
    torch.cuda.reset_peak_memory_stats()
    began=time.perf_counter();step=0;cache={};epochs=[]
    with (out/'training.jsonl').open('x') as log, (out/'fixed_free_and_first_cuts.jsonl').open('x') as fixed:
        for epoch in range(3):
            sums=Counter();counts=Counter();epoch_began=time.perf_counter()
            for ordinal,d in enumerate(draws):
                start=time.perf_counter();ident=d['window_id'];w=windows[ident];n=w['valid_bases']
                cached=torch.load(M.C.ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
                if (cached['revision']!=F.REVISION or len(cached['sequence'])!=n or
                    tuple(cached['features'].shape)!=(4096,2048) or
                    any(cached['window'][k]!=w[k] for k in ('id','species','seqid','strand','start','end','valid_bases'))):
                    raise ValueError('TRAIN feature/geometry mismatch')
                sequence=cached['sequence']
                features=cached['features'].float()[None].to('cuda')
                dna=one_hot(sequence,w['width'])[None].to('cuda')
                torch.cuda.synchronize();load_seconds=time.perf_counter()-start
                t=time.perf_counter()
                with torch.no_grad():
                    outputs=base(features,dna,[n]);h=outputs['features'][0,:n]
                torch.cuda.synchronize();forward_seconds=time.perf_counter()-t
                targets=build_targets(w,catalogs[w['species']],sequence)
                t=time.perf_counter();reused=ident in cache
                if not reused:
                    observer=FirstCutObserver(targets['positive_chains'])
                    with torch.no_grad():
                        p=outputs['segmentation_logits'][0,:n].softmax(-1)
                        free=generate_prefix(sequence,outputs['endpoint_logits'][0,:n].cpu().numpy(),
                                             p[:,list(CDS_COLUMNS)].sum(-1).cpu().numpy(),
                                             lambda pairs:base.chain.link_logits(h,pairs).cpu().numpy(),
                                             Budget(),observer=observer)
                    cache[ident]=(free,observer)
                    fixed.write(json.dumps({'window_id':ident,'free':free,'first_cuts':observer.cuts,
                                            'reference_only_observation':True},allow_nan=False)+'\n')
                    del p
                free,observer=cache[ident]
                torch.cuda.synchronize();proposal_seconds=time.perf_counter()-t
                t=time.perf_counter();optimizer.zero_grad(set_to_none=True)
                result=prefix_losses(model,h,free,observer,targets,d['importance_weight_uniform_within_stratum'])
                if not bool(torch.isfinite(result['total'])): raise ValueError('Nonfinite B2 loss')
                result['total'].backward()
                norm=torch.nn.utils.clip_grad_norm_(model.parameters(),1.,error_if_nonfinite=True)
                if any(p.grad is not None for p in base.parameters()): raise ValueError('Frozen B1 received gradient')
                optimizer.step();step+=1
                torch.cuda.synchronize();backward_seconds=time.perf_counter()-t
                losses={k:float(result[k].detach()) for k in ('final','rank','total')}
                detail={k:result[k] for k in ('first_cut_events','known_rank_pairs','unknown_competitors')}
                record={'step':step,'epoch':epoch+1,'draw_ordinal':ordinal,'window_id':ident,'species':w['species'],
                        'importance_weight':d['importance_weight_uniform_within_stratum'],'losses':losses,
                        'gradient_norm_before_clip':float(norm),'base_parameters_with_gradient':0,
                        'fixed_control_reused':reused,'free_candidates':len(free['chains']),
                        'reference_injected':result['origins'].count('reference_injected'),
                        'final_label_counts':{str(y):result['final_labels'].count(y) for y in (-1,0,1)},
                        'seconds':{'load':load_seconds,'frozen_forward':forward_seconds,
                                   'control_generation':proposal_seconds,'loss_backward_optimizer':backward_seconds},**detail}
                log.write(json.dumps(record,allow_nan=False)+'\n');sums.update(losses);counts.update(detail)
                if step%64==0:
                    log.flush();fixed.flush();check_output_cap()
                    print(json.dumps({'step':step,'epoch':epoch+1,'seconds':time.perf_counter()-began,
                                      'recent_losses':losses,'fixed_windows':len(cache)}),flush=True)
                del cached,features,dna,outputs,h,targets,result
            checkpoint=out/('step_'+str(step).zfill(6)+'.pt')
            torch.save({'arm':'B2','seed':0,'step':step,'hidden':128,'revision':F.REVISION,
                        'base_checkpoint':BASE,'prefix':model.state_dict(),'optimizer':optimizer.state_dict(),
                        'contract':CONTRACT},checkpoint)
            epochs.append({'epoch':epoch+1,'steps':1536,'mean_losses':{k:v/1536 for k,v in sums.items()},
                           'training_exposure_counts':dict(counts),'seconds':time.perf_counter()-epoch_began,
                           'checkpoint':str(checkpoint.relative_to(M.C.ROOT))})
            print(json.dumps(epochs[-1],allow_nan=False),flush=True)
    if step!=4608 or len(cache)!=1491: raise ValueError('Incomplete B2 first fit')
    report={'experiment':'M28-B2-FIRST-FIT','completed':True,'arm':'B2','optimizer_steps':step,
            'seed':0,'feature_revision':F.REVISION,'base_checkpoint':BASE,'contract':CONTRACT,
            'trainable_parameters':sum(p.numel() for p in model.parameters()),'base_trainable_parameters':0,
            'inherited_base_optimizer_steps':4608,
            'inherited_base_cost_source':'reports/M28-METHOD-RESTART/r5_result.json',
            'fixed_control_unique_windows':len(cache),'fixed_control_observation_policy':'off_policy_additive',
            'DEV_used_for_selection':False,'DEV_forward_used':False,'test_setaria_used':False,
            'epochs':epochs,'final_checkpoint':str(checkpoint.relative_to(M.C.ROOT)),
            'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
            'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),'artifact_bytes_before_summary':check_output_cap()}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    check_output_cap();print(json.dumps(report,indent=2,allow_nan=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.output_dir)
