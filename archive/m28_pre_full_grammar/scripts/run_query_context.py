#!/usr/bin/env python3
"""Matched A3 head with whole-window features: zero-update checks or one fixed fit."""
import argparse
import json
import time
from collections import Counter
import torch
import run_query_known_mean as A
import cache_draw_context as C
import cache_full_train as N

ROOT = N.ROOT
CONFIG = {**A.CONFIG,'encoder_context':C.FULL_SPEC}
PROPOSAL = 'reports/M28-METHOD-RESTART/a10_matched_context_protocol.md'
ENGINEERING = ROOT/'outputs/M28-A10-CONTEXT-ENGINEERING'
FIT = ROOT/'outputs/M28-A10-CONTEXT-FIT'


def load_inputs():
    if (C.OUT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Whole-window TRAIN cache must be terminal complete')
    prepared = json.loads((C.OUT/'summary.json').read_text())
    if (not prepared['completed'] or prepared['feature_spec'] != C.FULL_SPEC or
        prepared['unique_TRAIN_windows'] != 1491 or prepared['draws'] != 1536 or
        prepared['draw_order_weights_unchanged'] is not True or
        prepared['full_window_replay_exact'] is not True or
        prepared['optimizer_steps_this_run'] != 0 or prepared['test_setaria_used'] is not False):
        raise ValueError('Complete matching TRAIN representation required')
    draws, _, windows, catalogs = A.Q.load_train_inputs()
    if draws != N.rows(C.OUT/'draws.jsonl'):
        raise ValueError('Original draw list/order/weights changed')
    rows = N.rows(C.OUT/'index.jsonl')
    index = {r['window_id']:r for r in rows}
    if len(rows) != 1491 or len(index) != 1491:
        raise ValueError('All unique original draws required once in the new index')
    A.Q.validate_train_inputs(draws,index,windows,catalogs)
    for ordinal, row in enumerate(rows):
        expected = windows[row['window_id']]
        if (row['ordinal'] != ordinal or
            any(row['window'][k] != expected[k] for k in N.GEOMETRY) or
            (ROOT/row['artifact']).stat().st_size != row['bytes'] or row['bytes'] <= 0):
            raise ValueError('New feature geometry/artifact is incomplete or changed')
    return draws,index,windows,catalogs


def window_inputs(cached, window):
    if cached.get('feature_spec') != C.FULL_SPEC:
        raise ValueError('Four-block features cannot enter a whole-window head')
    return A.Q.window_inputs(cached,window,'cuda')


def check_last_identity(record, state=False):
    if (record['configuration'] != CONFIG or record['arm'] != 'A10_CONTEXT' or
        record['revision' if state else 'feature_revision'] != A.F.REVISION or
        record['step' if state else 'optimizer_steps'] != 4608 or
        (not state and not record['completed'])):
        raise ValueError('Only the matching completed last A10 context checkpoint is valid')


def engineering():
    draws,index,windows,catalogs = load_inputs()
    chosen = {}
    for d in draws:
        chosen.setdefault((d['species'],*d['stratum']),d)
    if len(chosen) != 12:
        raise ValueError('One existing draw per original twelve TRAIN strata required')
    model = A.new_model().train()
    if sum(p.numel() for p in model.parameters()) != 649457:
        raise ValueError('Matched A3 model parameter count changed')
    before = {k:v.detach().clone() for k,v in model.state_dict().items()}
    records = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    old_context_rejected = False
    for d in chosen.values():
        w = windows[d['window_id']]
        cached = torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
        if not old_context_rejected:
            # Same valid tensors/window/revision, but wrong encoding metadata.
            try:
                window_inputs({**cached,'feature_spec':N.FEATURE_SPEC},w)
            except ValueError:
                old_context_rejected = True
            else:
                raise ValueError('Old context with the same shape/revision was accepted')
        features,dna = window_inputs(cached,w)
        model.zero_grad(set_to_none=True)
        outputs = model(features,dna,[w['valid_bases']])
        result = A.supervised_readout(model,outputs,w,catalogs[w['species']],
                                     cached['sequence'],d['importance_weight_uniform_within_stratum'])
        result['total'].backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(),float('inf'),error_if_nonfinite=True)
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Frozen feature/DNA inputs required')
        labels = result['query']['query_labels']
        records.append({'window_id':w['id'],'species':w['species'],'stratum':d['stratum'],
                        'valid_bases':w['valid_bases'],'weighted_total':float(result['total'].detach()),
                        'gradient_norm':float(norm),'query_label_counts':{str(k):labels.count(k) for k in (-1,0,1)}})
        del cached,features,dna,outputs,result
    if not any(r['gradient_norm'] > 0 for r in records):
        raise ValueError('No nonzero head gradient observed in the fixed TRAIN checks')
    unchanged = all(torch.equal(before[k],v) for k,v in model.state_dict().items())
    if not unchanged:
        raise ValueError('Zero-update engineering altered the head state')
    N.write_summary(ENGINEERING,{
        'experiment':'M28-A10-CONTEXT-ENGINEERING','completed':True,'configuration':CONFIG,
        'matched_head_parameters':649457,'initialization':'same A3 seed0 constructor; no supervised inheritance',
        'TRAIN_windows':len(records),'records':records,'head_state_unchanged':True,
        'old_context_rejected':old_context_rejected,'original_draw_order_weights_unchanged':True,
        'zero_gradient_windows':sum(r['gradient_norm'] == 0 for r in records),
        'optimizer_steps_this_run':0,'head_checkpoint_read':False,'checkpoints_written':False,
        'new_backbone_forward':False,'DEV_forward_used':False,'test_setaria_used':False,
        'seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
        'limitations':'Real TRAIN wiring/gradient checks, not fitting, inference accuracy or scientific PASS.'})


def fit():
    if (ENGINEERING/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Matching zero-update head engineering must complete before fitting')
    checked = json.loads((ENGINEERING/'summary.json').read_text())
    if (not checked['completed'] or checked['configuration'] != CONFIG or
        not checked['head_state_unchanged'] or checked['optimizer_steps_this_run'] != 0):
        raise ValueError('Matched original head wiring required')
    draws,index,windows,catalogs = load_inputs()
    model = A.new_model().train()
    optimizer = torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.01,
                                 betas=(.9,.999),eps=1e-8)
    step = 0; epochs = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            sums = Counter(); counts = {s:Counter() for s in A.ALLOW}
            epoch_began = time.perf_counter()
            for ordinal,d in enumerate(draws):
                w = windows[d['window_id']]
                cached = torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
                features,dna = window_inputs(cached,w)
                optimizer.zero_grad(set_to_none=True)
                outputs = model(features,dna,[w['valid_bases']])
                result = A.supervised_readout(model,outputs,w,catalogs[w['species']],
                                             cached['sequence'],d['importance_weight_uniform_within_stratum'])
                result['total'].backward()
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(),1.,error_if_nonfinite=True)
                optimizer.step(); step += 1
                values = {k:float(v.detach()) for k,v in result['losses'].items()}
                labels = result['query']['query_labels']
                label_counts = {str(k):labels.count(k) for k in (-1,0,1)}
                total = float(result['total'].detach())
                log.write(json.dumps({'step':step,'epoch':epoch+1,'draw_ordinal':ordinal,
                                     'window_id':w['id'],'species':w['species'],
                                     'importance_weight':d['importance_weight_uniform_within_stratum'],
                                     'losses':values,'weighted_total':total,'gradient_norm_before_clip':float(norm),
                                     'query_label_counts':label_counts,'protection_reasons':result['protection_reasons']},
                                    allow_nan=False)+'\n')
                sums.update(dict(values,weighted_total=total))
                counts[w['species']].update({'query_'+k:v for k,v in label_counts.items()})
                counts[w['species']].update(result['protection_reasons'])
                if step % 64 == 0:
                    log.flush(); print(json.dumps({'step':step,'total':4608,
                                                  'seconds':time.perf_counter()-began}),flush=True)
                del cached,features,dna,outputs,result
            checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
            with checkpoint.open('xb') as handle:
                torch.save({'arm':'A10_CONTEXT','configuration':CONFIG,'step':step,
                            'revision':A.F.REVISION,'model':model.state_dict(),
                            'optimizer':optimizer.state_dict(),'proposal':PROPOSAL},handle)
            epochs.append({'epoch':epoch+1,'steps':1536,'mean_losses':{k:v/1536 for k,v in sums.items()},
                           'supervision_by_species':{s:dict(c) for s,c in counts.items()},
                           'seconds':time.perf_counter()-epoch_began,
                           'checkpoint':str(checkpoint.relative_to(ROOT))})
    if step != 4608:
        raise ValueError('Incomplete fixed matched context fit')
    N.write_summary(FIT,{
        'arm':'A10_CONTEXT','completed':True,'configuration':CONFIG,'feature_revision':A.F.REVISION,
        'proposal':PROPOSAL,'optimizer_steps':step,'epochs':epochs,
        'final_checkpoint':str(checkpoint.relative_to(ROOT)),'main_evaluation':'last_step_004608_only',
        'trainable_parameters':sum(p.numel() for p in model.parameters()),
        'initialization':'same A3 seed0 constructor; no supervised checkpoint inheritance',
        'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
        'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
        'DEV_used_for_selection':False,'test_setaria_used':False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('engineering','fit'),required=True)
    torch.set_num_threads(2)
    {'engineering':engineering,'fit':fit}[parser.parse_args().mode]()
