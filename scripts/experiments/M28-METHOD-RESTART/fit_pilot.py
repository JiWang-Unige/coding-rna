#!/usr/bin/env python3
"""Frozen-contract C0/B1 first fit: 3 passes of the same 1536 TRAIN draws."""
import argparse
import json
import time
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
from src.m28.core import M28Core
from src.m28.candidates import generate,Budget
from src.m28.labels import build_targets,local_loss,CDS_COLUMNS
from src.m28.training import joint_losses,one_hot

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def main(arm,feature_dir,out,contract='reports/M28-METHOD-RESTART/r4_pilot_contract.md'):
    out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'training.jsonl').exists(): raise FileExistsError('Training history already exists')
    source=json.loads((feature_dir/'summary.json').read_text())
    if source['revision']!=F.REVISION or source['DEV_windows']!=8690:
        raise ValueError('Unexpected shared feature cache')
    index={r['window_id']:r for r in rows(feature_dir/'index.jsonl')}
    base=M.C.ROOT/'outputs/M28-LABEL-POLICY-R4'
    catalogs={};windows={}
    for s in ('arabidopsis_thaliana','oryza_sativa'):
        catalogs[s]={t['id']:t for t in rows(base/s/'transcripts.jsonl')}
        windows.update({w['id']:w for w in rows(base/s/'windows_train.jsonl')})
    draws=rows(M.C.ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl')
    if len(draws)!=1536 or any(index[d['window_id']]['split']!='train' for d in draws):
        raise ValueError('TRAIN-only 1536 draw contract violated')
    torch.manual_seed(0);torch.set_num_threads(2)
    model=M28Core(arm,2048,128).to('cuda').train()
    optimizer=torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.01,betas=(.9,.999),eps=1e-8)
    torch.cuda.reset_peak_memory_stats()
    begin=time.perf_counter();step=0;epoch_summaries=[]
    free_handle=(out/'free_before_reference.jsonl').open('x') if arm=='B1' else None
    with (out/'training.jsonl').open('x') as log:
        for epoch in range(3):
            sums={};epoch_begin=time.perf_counter()
            for ordinal,d in enumerate(draws):
                start=time.perf_counter()
                ident=d['window_id'];w=windows[ident]
                cached=torch.load(M.C.ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
                if cached['revision']!=F.REVISION or any(cached['window'][k]!=w[k] for k in ('id','start','end','strand','valid_bases')):
                    raise ValueError('Feature/window alignment mismatch')
                n=w['valid_bases'];sequence=cached['sequence']
                if len(sequence)!=n or tuple(cached['features'].shape)!=(4096,2048):
                    raise ValueError('Feature geometry mismatch')
                features=cached['features'].float()[None].to('cuda')
                dna=one_hot(sequence,w['width'])[None].to('cuda')
                torch.cuda.synchronize();load_seconds=time.perf_counter()-start
                optimizer.zero_grad(set_to_none=True)
                start=time.perf_counter()
                outputs=model(features,dna,[n])
                torch.cuda.synchronize();forward_seconds=time.perf_counter()-start
                proposal_seconds=0.;detail={}
                if arm=='B1':
                    start=time.perf_counter()
                    with torch.no_grad():
                        h=outputs['features'][0,:n].detach()
                        p=outputs['segmentation_logits'][0,:n].softmax(-1)
                        free=generate(sequence,outputs['endpoint_logits'][0,:n].detach().cpu().numpy(),
                                      p[:,list(CDS_COLUMNS)].sum(-1).cpu().numpy(),
                                      link_score_fn=lambda pairs:model.chain.link_logits(h,pairs).cpu().numpy(),
                                      budget=Budget())
                    torch.cuda.synchronize();proposal_seconds=time.perf_counter()-start
                    free_handle.write(json.dumps({'step':step+1,'window_id':ident,'free':free})+'\n')
                    free_handle.flush()
                start=time.perf_counter()
                targets=build_targets(w,catalogs[w['species']],sequence)
                weight=d['importance_weight_uniform_within_stratum']
                if arm=='B1':
                    joint=joint_losses(model,outputs,free,w,catalogs[w['species']],targets,window_weight=weight)
                    losses=joint['losses'];loss=joint['total']
                    detail={'free_candidates':len(free['chains']),'training_candidates':len(joint['training_chains']),
                            'reference_injected':joint['origins'].count('reference_injected'),
                            'chain_positive':joint['chain_labels'].count(1),'chain_negative':joint['chain_labels'].count(0),
                            'chain_unknown':joint['chain_labels'].count(-1),'origins':joint['origins']}
                else:
                    losses={'local':local_loss(outputs['segmentation_logits'][0],targets)}
                    loss=losses['local']*weight
                if not bool(torch.isfinite(loss)): raise ValueError('Nonfinite training loss')
                loss.backward()
                norm=torch.nn.utils.clip_grad_norm_(model.parameters(),1.,error_if_nonfinite=True)
                optimizer.step();step+=1
                torch.cuda.synchronize();backward_step_seconds=time.perf_counter()-start
                values={k:float(v.detach()) for k,v in losses.items()}
                record={'step':step,'epoch':epoch+1,'draw_ordinal':ordinal,'window_id':ident,'species':w['species'],
                        'importance_weight':weight,'losses':values,'weighted_total':float(loss.detach()),
                        'gradient_norm_before_clip':float(norm),'load_seconds':load_seconds,
                        'shared_forward_seconds':forward_seconds,'proposal_seconds':proposal_seconds,
                        'supervision_backward_optimizer_seconds':backward_step_seconds,**detail}
                log.write(json.dumps(record)+'\n')
                for k,v in dict(values,weighted_total=record['weighted_total']).items(): sums[k]=sums.get(k,0.)+v
                if step%64==0:
                    log.flush();print(json.dumps({'arm':arm,'step':step,'epoch':epoch+1,
                                                 'seconds':time.perf_counter()-begin,'recent':values}),flush=True)
                if arm=='B1': del joint,free,h,p
                del outputs,features,dna,loss,losses,cached,targets
            checkpoint=out/('step_'+str(step).zfill(6)+'.pt')
            torch.save({'arm':arm,'seed':0,'feature_dim':2048,'hidden':128,'step':step,'revision':F.REVISION,
                        'model':model.state_dict(),'optimizer':optimizer.state_dict(),
                        'contract':contract},checkpoint)
            epoch_summaries.append({'epoch':epoch+1,'steps':1536,'mean_losses':{k:v/1536 for k,v in sums.items()},
                                    'seconds':time.perf_counter()-epoch_begin,
                                    'checkpoint':str(checkpoint.relative_to(M.C.ROOT))})
            print(json.dumps(epoch_summaries[-1]),flush=True)
    if free_handle: free_handle.close()
    assert step==4608
    report={'arm':arm,'contract':contract,'optimizer_steps':step,'seed':0,'feature_revision':F.REVISION,
            'trainable_parameters':sum(p.numel() for p in model.parameters()),'epochs':epoch_summaries,
            'training_loop_seconds':time.perf_counter()-begin,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
            'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),'DEV_used_for_selection':False,
            'final_checkpoint':str(checkpoint.relative_to(M.C.ROOT)),'test_setaria_used':False,
            'main_evaluation':'last_step_004608_only'}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--arm',choices=('C0','B1'),required=True)
    p.add_argument('--feature-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--contract',default='reports/M28-METHOD-RESTART/r4_pilot_contract.md')
    a=p.parse_args();main(a.arm,a.feature_dir,a.output_dir,a.contract)
