#!/usr/bin/env python3
"""Finite GPU timing using existing twelve TRAIN feature caches; zero optimizer steps."""
import argparse
import json
import time
from collections import Counter
from pathlib import Path
import torch
import manifest as M
from src.m28.core import M28Core
from src.m28.candidates import generate, Budget
from src.m28.labels import build_targets, local_loss, CDS_COLUMNS
from src.m28.training import one_hot, joint_losses

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def main(feature_dir,out):
    out.mkdir(parents=True,exist_ok=True)
    source=json.loads((feature_dir/'summary.json').read_text())
    catalogs={s:{r['id']:r for r in rows(M.C.ROOT/'outputs/M28-DATA-MANIFEST-R2'/s/'transcripts.jsonl')}
              for s in ('arabidopsis_thaliana','oryza_sativa')}
    torch.set_num_threads(2)
    report={'optimizer_steps':0,'feature_revision':source['revision'],'GPU':torch.cuda.get_device_name(),
            'feature_extraction_repeated':False,'val_test_setaria_used':False,'arms':{}}
    for arm in ('C0','B1'):
        torch.manual_seed(0)
        model=M28Core(arm,source['hidden_size']).to('cuda').train()
        torch.cuda.reset_peak_memory_stats()
        records=[]
        for r in source['records']:
            load_start=time.perf_counter()
            cached=torch.load(feature_dir/r['artifact'],map_location='cpu',weights_only=True)
            w=cached['window'];seq=cached['sequence'];n=w['valid_bases']
            feat=cached['features'].float()[None].to('cuda')
            dna=one_hot(seq,w['width'])[None].to('cuda')
            torch.cuda.synchronize()
            load_seconds=time.perf_counter()-load_start
            model.zero_grad(set_to_none=True)
            start=time.perf_counter()
            outputs=model(feat,dna,[n])
            torch.cuda.synchronize()
            forward_seconds=time.perf_counter()-start
            detail={};proposal_seconds=0.
            if arm=='B1':
                start=time.perf_counter()
                with torch.no_grad():
                    h=outputs['features'][0,:n].detach()
                    p=outputs['segmentation_logits'][0,:n].softmax(-1)
                    free=generate(seq,outputs['endpoint_logits'][0,:n].detach().cpu().numpy(),
                                  p[:,list(CDS_COLUMNS)].sum(-1).cpu().numpy(),
                                  link_score_fn=lambda pairs:model.chain.link_logits(h,pairs).cpu().numpy(),
                                  budget=Budget())
                torch.cuda.synchronize()
                proposal_seconds=time.perf_counter()-start
                # Free candidate evidence precedes reference supervision.
                with (out/('free_'+str(r['ordinal']).zfill(2)+'.json')).open('x') as f:
                    json.dump(dict(free,window_id=w['id']),f)
            targets=build_targets(w,catalogs[w['species']],seq)
            start=time.perf_counter()
            weight=cached['draw']['importance_weight_uniform_within_stratum']
            if arm=='C0':
                loss=local_loss(outputs['segmentation_logits'][0],targets)*weight
            else:
                joint=joint_losses(model,outputs,free,w,catalogs[w['species']],targets,window_weight=weight)
                loss=joint['total']
                # Check actual positive/negative/unknown connection directions, including rice multi-exon windows.
                lg=torch.autograd.grad(loss,joint['link_logits'],retain_graph=True)[0]
                labels=torch.as_tensor(joint['link_labels'],device=lg.device)
                assert bool((lg[labels==1]<0).all() and (lg[labels==0]>0).all() and (lg[labels<0]==0).all())
                detail={'free_chains':len(free['chains']),'training_chains':len(joint['training_chains']),
                        'link_labels':dict(Counter(map(str,joint['link_labels']))),
                        'chain_labels':dict(Counter(map(str,joint['chain_labels']))),
                        'reference_injected':joint['origins'].count('reference_injected'),
                        'link_gradient_direction_pass':True}
            if not torch.isfinite(loss): raise ValueError('Nonfinite joint loss')
            loss.backward()
            torch.cuda.synchronize()
            backward_seconds=time.perf_counter()-start
            gradients=[p.grad for p in model.parameters() if p.grad is not None]
            if not gradients or not all(bool(torch.isfinite(g).all()) for g in gradients):
                raise ValueError('Missing/nonfinite model gradient')
            record={'window_id':w['id'],'species':w['species'],'load_seconds':load_seconds,
                    'shared_forward_seconds':forward_seconds,'proposal_seconds':proposal_seconds,
                    'supervision_and_backward_seconds':backward_seconds,'weighted_loss':float(loss.detach()),**detail}
            records.append(record);print(json.dumps(dict(arm=arm,**record)),flush=True)
            if arm=='B1': del joint,lg,labels,h,p,free
            del outputs,loss,feat,dna,gradients
        report['arms'][arm]={'records':records,'trainable_parameters':sum(p.numel() for p in model.parameters()),
                             'peak_allocated_bytes':torch.cuda.max_memory_allocated(),
                             'peak_reserved_bytes':torch.cuda.max_memory_reserved()}
        del model
        torch.cuda.empty_cache()
    for species in catalogs:
        items=[r for r in report['arms']['B1']['records'] if r['species']==species]
        assert any(r['link_labels'].get('1',0)>0 and r['link_labels'].get('0',0)>0 for r in items), species
    report['actual_positive_and_negative_link_gradients_both_species']=True
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--feature-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.feature_dir,a.output_dir)
