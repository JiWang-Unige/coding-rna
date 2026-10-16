#!/usr/bin/env python3
"""Fixed A3 known-query-mean fit and full reference-free DEV inference; no parameter search."""
import argparse
import json
import time
from collections import Counter, defaultdict
import torch
import run_query as Q
from run_query import (ROOT, FEATURES, F, M, ALLOW, STEPS, load_train_inputs, window_inputs,
                       DEV_LENGTHS, dev_groups, read_rows,
                       query_window_candidates)
from src.m28.query_collaboration import CollaborativeGeneQueryModel
from src.m28.query_known_mean import joint_known_query_losses

PROPOSAL = 'reports/M28-METHOD-RESTART/a3_known_query_mean_protocol.md'
CONFIG = {**Q.CONFIG, 'query_collaboration':{
    'heads':4, 'dropout':0., 'pre_norm':True, 'residual':True, 'FFN':False},
    'existence_reduction':'mean_over_known_queries'}
STAGE_DIRS = {mode:ROOT/'outputs'/name for mode,name in (
    ('engineering','M28-A3-QUERY-ENGINEERING'), ('fit','M28-A3-FIRST-FIT'),
    ('infer','M28-A3-DEV-INFER'), ('evaluate','M28-A3-EVALUATION'))}


def artifact_bytes():
    return sum(p.stat().st_size for out in STAGE_DIRS.values()
               for p in out.rglob('*') if p.is_file())


def check_last_identity(record, state=False):
    if (record['configuration'] != CONFIG or record['arm'] != 'A3' or
        record['revision' if state else 'feature_revision'] != F.REVISION or
        record['step' if state else 'optimizer_steps'] != STEPS or
        (not state and not record['completed'])):
        raise ValueError('Only the matching completed A3 last checkpoint may be evaluated')


def new_model():
    torch.manual_seed(CONFIG['seed'])
    return CollaborativeGeneQueryModel(CONFIG['feature_dim'],CONFIG['num_queries'],CONFIG['hidden']).to('cuda')


def write_summary(out, report):
    report['artifact_bytes_before_summary'] = artifact_bytes()
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('scopes','records')},allow_nan=False),flush=True)


def supervised_readout(model, outputs, window, catalog, sequence, weight):
    targets = Q.build_targets(window, catalog, sequence)
    result = joint_known_query_losses(model, outputs, targets, window['valid_bases'],
                                     weight, CONFIG['loss_coefficients'])
    if not bool(torch.isfinite(result['total'])):
        raise ValueError('Nonfinite A3 training loss')
    return result


def fit():
    if (STAGE_DIRS['engineering']/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('A3 synthetic collaboration checks must complete before fitting')
    prepared = json.loads((Q.STAGE_DIRS['prepare']/'summary.json').read_text())
    Q.validate_preparation(prepared)
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
            torch.save({'arm':'A3','configuration':CONFIG,'step':step,'revision':F.REVISION,
                        'model':model.state_dict(),'optimizer':optimizer.state_dict(),'proposal':PROPOSAL},checkpoint)
            artifact_bytes()
            epochs.append({'epoch':epoch+1,'steps':1536,'mean_losses':{k:v/1536 for k,v in sums.items()},
                           'supervision_by_species':{s:dict(c) for s,c in counts.items()},
                           'seconds':time.perf_counter()-epoch_began,'checkpoint':str(checkpoint.relative_to(ROOT))})
    if step != STEPS: raise ValueError('Incomplete A3 fit')
    write_summary(out,{'arm':'A3','completed':True,'configuration':CONFIG,'feature_revision':F.REVISION,
                       'proposal':PROPOSAL,'optimizer_steps':step,'epochs':epochs,
                       'final_checkpoint':str(checkpoint.relative_to(ROOT)),'main_evaluation':'last_step_004608_only',
                       'trainable_parameters':sum(p.numel() for p in model.parameters()),
                       'initialization':'same complete A2 seed0 initialization; no supervised checkpoint inheritance',
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
    if done != 8690: raise ValueError('Incomplete A3 full DEV inference')
    write_summary(out,{'arm':'A3','completed':True,'inference_complete':True,'windows':done,
                       'optimizer_steps':STEPS,'configuration':CONFIG,'feature_revision':F.REVISION,
                       'proposal':PROPOSAL,'checkpoint':str(checkpoint.relative_to(ROOT)),
                       'reference_used':False,'test_setaria_used':False,'scopes':scopes,
                       'wall_seconds':time.perf_counter()-began,
                       'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                       'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('fit','infer'),required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    {'fit':fit,'infer':infer}[args.mode]()
