#!/usr/bin/env python3
"""Fixed A single-fit pipeline. Approval scope is recorded in the proposal.

Only the explicit prepare/fit/infer CLI runs touch real inputs or CUDA. Imports
and synthetic helper tests do neither. This cached-feature DEV pilot is not a
standalone raw-FASTA caller, independent species test, or publication result.
"""
import argparse
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
from infer_pilot import DEV_LENGTHS, dev_groups, read_rows
from query_training_census import ALLOW, train_records
from src.m28.labels import build_targets, oriented_parts
from src.m28.training import one_hot
from src.m28.query_model import GeneQueryModel, joint_query_losses
from src.m28.query_output import query_window_candidates

ROOT = M.C.ROOT
FEATURES = ROOT/'outputs/M28-PILOT-R4/features'
LABELS = ROOT/'outputs/M28-LABEL-POLICY-R4'
DRAWS = ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl'
PROPOSAL = 'reports/M28-METHOD-RESTART/a_first_fit_proposal.md'
CONFIG = {'feature_dim':2048, 'hidden':128, 'num_queries':16, 'seed':0,
          'loss_coefficients':[1.,1.,1.,1.], 'max_exons':128,
          'existence_logit_threshold':0.}
STEPS = 4608
STAGE_DIRS = {mode:ROOT/'outputs'/name for mode,name in (
    ('prepare','M28-A-QUERY-PREPARATION'), ('fit','M28-A-FIRST-FIT'),
    ('infer','M28-A-DEV-INFER'), ('evaluate','M28-A-EVALUATION'))}


def artifact_bytes():
    size = sum(p.stat().st_size for out in STAGE_DIRS.values()
               for p in out.rglob('*') if p.is_file())
    return size


def validate_train_inputs(draws, index, windows, catalogs):
    ids = {d['window_id'] for d in draws}
    if (len(draws) != 1536 or len(ids) != 1491 or
        Counter(d['species'] for d in draws) != Counter({s:768 for s in ALLOW}) or
        set(index) != ids or set(windows) != ids):
        raise ValueError('Original TRAIN draw/window/index identity mismatch')
    for d in draws:
        ident = d['window_id']; w = windows[ident]
        if (w['split'] != 'train' or index[ident]['split'] != 'train' or
            w['seqid'] not in ALLOW[w['species']] or
            w['label_policy'] != 'CDS_assessable_primary_R4' or
            (w['species'],w['seqid'],w['strand']) != (d['species'],*d['stratum']) or
            index[ident]['species'] != w['species']):
            raise ValueError('TRAIN split/stratum/R4 label mismatch')
        weight = float(d['importance_weight_uniform_within_stratum'])
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError('Original TRAIN importance mismatch')
        if not set(w['positive_ids']) <= set(w['overlapping_transcript_ids']):
            raise ValueError('Positive TRAIN target join mismatch')
        for tid in w['overlapping_transcript_ids']:
            tx = catalogs[w['species']][tid]
            if (tx['split'] != 'train' or
                (tx['seqid'],tx['strand']) != (w['seqid'],w['strand'])):
                raise ValueError('TRAIN transcript/window identity mismatch')
            if tid in w['positive_ids'] and not (tx['training_primary'] and tx['eligible_positive']):
                raise ValueError('TRAIN positive eligibility mismatch')


def load_train_inputs():
    summary = json.loads((FEATURES/'summary.json').read_text())
    if (summary['revision'] != F.REVISION or summary['train_unique_windows'] != 1491 or
        summary['DEV_windows'] != 8690):
        raise ValueError('Shared feature cache identity mismatch')
    draws = read_rows(DRAWS)
    ids = {d['window_id'] for d in draws}
    index = {r['window_id']:r for r in train_records(FEATURES/'index.jsonl') if r['window_id'] in ids}
    windows, catalogs = {}, {}
    for species in ALLOW:
        local = {w['id']:w for w in train_records(LABELS/species/'windows_train.jsonl') if w['id'] in ids}
        needed = {tid for w in local.values() for tid in w['overlapping_transcript_ids']}
        # Existing mixed catalog: pre-parse TRAIN filter, not read_rows.
        catalogs[species] = {t['id']:t for t in train_records(LABELS/species/'transcripts.jsonl') if t['id'] in needed}
        if set(catalogs[species]) != needed:
            raise ValueError('TRAIN transcript join coverage mismatch')
        windows.update(local)
    validate_train_inputs(draws,index,windows,catalogs)
    return draws,index,windows,catalogs


def validate_cached(cached, expected):
    keys = ('id','species','split','seqid','strand','start','end','width','valid_bases')
    n = expected['valid_bases']
    if (cached['revision'] != F.REVISION or cached['model'] != F.MODEL or
        any(cached['window'][k] != expected[k] for k in keys) or
        expected['width'] != 24576 or not 0 < n <= 24576 or
        expected['end']-expected['start'] != n or len(cached['sequence']) != n or
        tuple(cached['features'].shape) != (4096,2048)):
        raise ValueError('A feature/window/sequence geometry mismatch')
    # Cached R2/R3 label fields are deliberately not read or used as targets.


def window_inputs(cached, expected, device):
    validate_cached(cached,expected)
    features = cached['features'].float()[None].to(device)
    dna = one_hot(cached['sequence'],expected['width'])[None].to(device)
    return features,dna


def supervised_readout(model, outputs, window, catalog, sequence, weight):
    targets = build_targets(window,catalog,sequence)
    result = joint_query_losses(model,outputs,targets,window['valid_bases'],
                                weight,CONFIG['loss_coefficients'])
    if not bool(torch.isfinite(result['total'])):
        raise ValueError('Nonfinite A training loss')
    return result


def select_timing_draws(draws, windows, catalogs):
    profiles = []
    for ordinal,d in enumerate(draws):
        w = windows[d['window_id']]; catalog = catalogs[w['species']]
        positive = [tuple(oriented_parts(catalog[t],w)) for t in w['positive_ids']]
        protected = {tuple(oriented_parts(catalog[t],w)) for t in w['overlapping_transcript_ids']}
        protected = {p for p in protected if p and p[0][0] >= 0 and p[-1][1] <= w['valid_bases']}
        uncertain = bool(w['chain_unknown_intervals']) or bool(protected-set(positive))
        profiles.append({'ordinal':ordinal,'draw':d,'positive':len(positive),
                         'exons':sum(map(len,positive)),'metadata_protected':uncertain})
    selected = []
    for species in ALLOW:
        local = [p for p in profiles if p['draw']['species'] == species]
        categories = (
            ('most_positive',local,lambda p:p['positive']),
            ('most_positive_exons',local,lambda p:p['exons']),
            ('zero_positive_unprotected',[p for p in local if not p['positive'] and not p['metadata_protected']],lambda p:0),
            ('positive_protected',[p for p in local if p['positive'] and p['metadata_protected']],lambda p:0))
        for category,candidates,key in categories:
            if not candidates:
                raise ValueError('Missing fixed TRAIN timing category: '+species+' '+category)
            chosen = min(candidates,key=lambda p:(-key(p),p['ordinal']))
            selected.append({'species':species,'category':category,**chosen})
    return selected


def budget_gate(records):
    train = max(r['load_to_backward_seconds'] for r in records)
    decode = max(r['seconds']['decode'] for r in records)
    if not all(math.isfinite(t) and t >= 0 for t in (train,decode)):
        raise ValueError('Nonfinite timing budget input')
    fit_seconds, infer_decode_seconds = train*STEPS+1800, decode*8690+1800
    return {'fit_projected_seconds':fit_seconds,'infer_decode_projected_seconds':infer_decode_seconds,
            'passed':fit_seconds <= 4*3600 and infer_decode_seconds <= 6*3600,
            'inference_projection_scope':'decode component only, not a complete runtime guarantee'}


def check_last_identity(record, state=False):
    if (record['configuration'] != CONFIG or record['arm'] != 'A' or
        record['revision' if state else 'feature_revision'] != F.REVISION or
        record['step' if state else 'optimizer_steps'] != STEPS or
        (not state and not record['completed'])):
        raise ValueError('Only the matching completed A last checkpoint may be evaluated')


def validate_preparation(prepared):
    if (not prepared['completed'] or prepared['configuration'] != CONFIG or prepared['optimizer_steps'] != 0 or
        prepared['feature_revision'] != F.REVISION or prepared['windows'] != 8):
        raise ValueError('A requires completed fixed TRAIN preparation identity')


def new_model():
    torch.manual_seed(CONFIG['seed'])
    return GeneQueryModel(CONFIG['feature_dim'],CONFIG['num_queries'],CONFIG['hidden']).to('cuda')


def write_summary(out, report):
    report['artifact_bytes_before_summary'] = artifact_bytes()
    with (out/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    artifact_bytes()
    print(json.dumps({k:v for k,v in report.items() if k not in ('scopes','records')},allow_nan=False),flush=True)


def prepare():
    out = STAGE_DIRS['prepare']; out.mkdir(parents=True,exist_ok=True)
    draws,index,windows,catalogs = load_train_inputs()
    sample = select_timing_draws(draws,windows,catalogs)
    model = new_model().train(); torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter(); records = []
    with (out/'windows.jsonl').open('x') as log:
        for selected in sample:
            d = selected['draw']; w = windows[d['window_id']]; seconds = {}
            t0 = time.perf_counter()
            cached = torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
            features,dna = window_inputs(cached,w,'cuda')
            torch.cuda.synchronize(); seconds['load'] = time.perf_counter()-t0
            model.zero_grad(set_to_none=True); t = time.perf_counter()
            outputs = model(features,dna,[w['valid_bases']])
            torch.cuda.synchronize(); seconds['forward'] = time.perf_counter()-t
            t = time.perf_counter()
            result = supervised_readout(model,outputs,w,catalogs[w['species']],
                                        cached['sequence'],d['importance_weight_uniform_within_stratum'])
            torch.cuda.synchronize(); seconds['matching_and_loss_readout'] = time.perf_counter()-t
            t = time.perf_counter(); result['total'].backward()
            norm = torch.nn.utils.clip_grad_norm_(model.parameters(),1.,error_if_nonfinite=True)
            torch.cuda.synchronize(); seconds['backward'] = time.perf_counter()-t
            full_step = time.perf_counter()-t0
            t = time.perf_counter()
            free = query_window_candidates(model.query,outputs['features'][0],cached['sequence'],0,w['id'],
                                           CONFIG['max_exons'],CONFIG['existence_logit_threshold'],w['valid_bases'])
            torch.cuda.synchronize(); seconds['decode'] = time.perf_counter()-t
            record = {'species':w['species'],'category':selected['category'],'draw_ordinal':selected['ordinal'],
                      'window_id':w['id'],'optimizer_steps':0,'load_to_backward_seconds':full_step,
                      'query_label_counts':{str(y):result['query']['query_labels'].count(y) for y in (-1,0,1)},
                      'protection_reasons':result['protection_reasons'],'gradient_norm_before_clip':float(norm),
                      'losses':{k:float(v.detach()) for k,v in result['losses'].items()},
                      'raw_queries':free['raw_queries'],'query_counts':free['counts'],'seconds':seconds}
            t = time.perf_counter(); log.write(json.dumps(record,allow_nan=False)+'\n'); log.flush()
            seconds['serialization'] = time.perf_counter()-t
            records.append(record); artifact_bytes()
            del cached,features,dna,outputs,result,free
    report = {'arm':'A','completed':True,'configuration':CONFIG,'feature_revision':F.REVISION,
              'proposal':PROPOSAL,'optimizer_steps':0,'fitting':False,'windows':len(records),
              'unique_windows':len({r['window_id'] for r in records}),'records':records,
              'budget_gate':budget_gate(records),'wall_seconds':time.perf_counter()-began,
              'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
              'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
              'DEV_forward_used':False,'test_setaria_used':False}
    write_summary(out,report)


def fit():
    prepared = json.loads((STAGE_DIRS['prepare']/'summary.json').read_text())
    validate_preparation(prepared)
    out = STAGE_DIRS['fit']; out.mkdir(parents=True,exist_ok=True)
    draws,index,windows,catalogs = load_train_inputs()
    model = new_model().train()
    optimizer = torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.01,betas=(.9,.999),eps=1e-8)
    torch.cuda.reset_peak_memory_stats(); began = time.perf_counter(); step = 0; epochs = []
    with (out/'training.jsonl').open('x') as log:
        for epoch in range(3):
            sums = Counter(); counts = {s:Counter() for s in ALLOW}; epoch_began = time.perf_counter()
            for ordinal,d in enumerate(draws):
                t0 = time.perf_counter(); w = windows[d['window_id']]
                cached = torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
                features,dna = window_inputs(cached,w,'cuda')
                torch.cuda.synchronize(); load_seconds = time.perf_counter()-t0
                optimizer.zero_grad(set_to_none=True); t = time.perf_counter()
                outputs = model(features,dna,[w['valid_bases']])
                result = supervised_readout(model,outputs,w,catalogs[w['species']],
                                            cached['sequence'],d['importance_weight_uniform_within_stratum'])
                result['total'].backward()
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(),1.,error_if_nonfinite=True)
                optimizer.step(); step += 1
                torch.cuda.synchronize(); step_seconds = time.perf_counter()-t
                labels = result['query']['query_labels']
                label_counts = {str(y):labels.count(y) for y in (-1,0,1)}
                values = {k:float(v.detach()) for k,v in result['losses'].items()}
                record = {'step':step,'epoch':epoch+1,'draw_ordinal':ordinal,'window_id':w['id'],'species':w['species'],
                          'importance_weight':d['importance_weight_uniform_within_stratum'],'losses':values,
                          'weighted_total':float(result['total'].detach()),'gradient_norm_before_clip':float(norm),
                          'query_label_counts':label_counts,'protection_reasons':result['protection_reasons'],
                          'gene_to_query':result['query']['gene_to_query'],'load_seconds':load_seconds,
                          'forward_supervision_backward_optimizer_seconds':step_seconds}
                log.write(json.dumps(record,allow_nan=False)+'\n')
                sums.update(dict(values,weighted_total=record['weighted_total']))
                counts[w['species']].update({'query_'+k:v for k,v in label_counts.items()})
                counts[w['species']].update(result['protection_reasons'])
                if step % 64 == 0:
                    log.flush(); artifact_bytes()
                    print(json.dumps({'step':step,'seconds':time.perf_counter()-began,'recent':values}),flush=True)
                del cached,features,dna,outputs,result
            checkpoint = out/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists(): raise FileExistsError(checkpoint)
            torch.save({'arm':'A','configuration':CONFIG,'step':step,'revision':F.REVISION,
                        'model':model.state_dict(),'optimizer':optimizer.state_dict(),'proposal':PROPOSAL},checkpoint)
            artifact_bytes()
            epochs.append({'epoch':epoch+1,'steps':1536,'mean_losses':{k:v/1536 for k,v in sums.items()},
                           'supervision_by_species':{s:dict(c) for s,c in counts.items()},
                           'seconds':time.perf_counter()-epoch_began,'checkpoint':str(checkpoint.relative_to(ROOT))})
    if step != STEPS: raise ValueError('Incomplete A fit')
    write_summary(out,{'arm':'A','completed':True,'configuration':CONFIG,'feature_revision':F.REVISION,
                       'proposal':PROPOSAL,'optimizer_steps':step,'epochs':epochs,
                       'final_checkpoint':str(checkpoint.relative_to(ROOT)),'main_evaluation':'last_step_004608_only',
                       'trainable_parameters':sum(p.numel() for p in model.parameters()),
                       'initialization':'fresh common C0 seed0, no supervised checkpoint inheritance',
                       'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                       'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
                       'DEV_used_for_selection':False,'test_setaria_used':False})


def infer():
    # Geometry/index only. No TRAIN loader, label catalog or evaluator calls.
    fitted = json.loads((STAGE_DIRS['fit']/'summary.json').read_text()); check_last_identity(fitted)
    checkpoint = ROOT/fitted['final_checkpoint']
    state = torch.load(checkpoint,map_location='cpu',weights_only=True); check_last_identity(state,state=True)
    model = new_model(); model.load_state_dict(state['model'],strict=True); del state
    model.eval().requires_grad_(False)
    groups = dev_groups(read_rows(FEATURES/'index.jsonl'))
    out = STAGE_DIRS['infer']; out.mkdir(parents=True,exist_ok=True)
    began = time.perf_counter(); done = 0; scopes = []; torch.cuda.reset_peak_memory_stats()
    with torch.inference_mode():
        for scope_number,((species,seqid,strand),rows) in enumerate(groups.items()):
            file = out/('scope_'+str(scope_number).zfill(2)+'.jsonl')
            windows = []; counts = Counter(); seconds = defaultdict(float)
            with file.open('x') as handle:
                for i,row in enumerate(rows):
                    t = time.perf_counter()
                    cached = torch.load(ROOT/row['artifact'],map_location='cpu',weights_only=True)
                    a = int(row['window_id'].split('|')[2]); b = min(a+24576,DEV_LENGTHS[species][seqid])
                    w = {'id':row['window_id'],'species':species,'split':'val','seqid':seqid,'strand':strand,
                         'start':a,'end':b,'width':24576,'valid_bases':b-a}
                    features,dna = window_inputs(cached,w,'cuda'); windows.append([a,b])
                    torch.cuda.synchronize(); seconds['load'] += time.perf_counter()-t
                    t = time.perf_counter(); outputs = model(features,dna,[b-a])
                    torch.cuda.synchronize(); seconds['shared_forward'] += time.perf_counter()-t
                    t = time.perf_counter()
                    free = query_window_candidates(model.query,outputs['features'][0],cached['sequence'],i,w['id'],
                                                   CONFIG['max_exons'],CONFIG['existence_logit_threshold'],b-a)
                    torch.cuda.synchronize(); seconds['decode'] += time.perf_counter()-t
                    t = time.perf_counter(); handle.write(json.dumps(free,allow_nan=False)+'\n')
                    seconds['serialization'] += time.perf_counter()-t
                    counts.update(free['counts']); counts['candidate_records'] += len(free['chains'])
                    counts['all_queries_positive_completed'] += int(
                        sum(score > 0 for score in free['score']) == CONFIG['num_queries'])
                    done += 1
                    if done % 32 == 0:
                        handle.flush(); artifact_bytes()
                        print(json.dumps({'windows':done,'total':8690,'seconds':time.perf_counter()-began}),flush=True)
                    del cached,features,dna,outputs,free
            scopes.append({'species':species,'seqid':seqid,'strand':strand,'windows':windows,
                           'file':file.name,'candidate_records':counts['candidate_records'],
                           'query_counts':dict(counts),'seconds':dict(seconds)})
    if done != 8690: raise ValueError('Incomplete A full DEV inference')
    write_summary(out,{'arm':'A','completed':True,'inference_complete':True,'windows':done,
                       'optimizer_steps':STEPS,'configuration':CONFIG,'feature_revision':F.REVISION,
                       'proposal':PROPOSAL,'checkpoint':str(checkpoint.relative_to(ROOT)),
                       'reference_used':False,'test_setaria_used':False,'scopes':scopes,
                       'wall_seconds':time.perf_counter()-began,
                       'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                       'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('prepare','fit','infer'),required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    {'prepare':prepare,'fit':fit,'infer':infer}[args.mode]()
