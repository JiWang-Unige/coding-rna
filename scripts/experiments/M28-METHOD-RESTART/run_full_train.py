#!/usr/bin/env python3
"""Matched full-TRAIN C0/query fit. No DEV reference, model selection or search."""
import argparse
import json
import math
import random
import time
from collections import Counter
import torch
import run_query as Q
from prepare_full_train import OUT as PLAN, FEATURES, FEATURE_SPEC, GEOMETRY
from query_training_census import ALLOW, train_records
from src.m28.core import M28Core
from src.m28.query_collaboration import CollaborativeGeneQueryModel
from src.m28.query_known_mean import joint_known_query_losses
from src.m28.labels import build_targets, local_loss

ROOT = Q.ROOT
LABELS = ROOT/'outputs/M28-LABEL-POLICY-R4'
READY = ROOT/'outputs/M28-FULL-TRAIN-READINESS'
SMOKE = ROOT/'outputs/M28-FULL-TRAIN-GPU-SMOKE'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_train_matched_protocol.md'
N = 25512
PASSES = 3
STEPS = N*PASSES
COMMON = {'feature_dim':2048,'hidden':128,'seed':0,'passes':PASSES,
          'windows_per_pass':N,'batch_size':1,'optimizer':'AdamW',
          'lr':3e-4,'weight_decay':.01,'betas':[.9,.999],'eps':1e-8,
          'clip_norm':1.,'precision':'FP32','training_order_seed':0,
          'importance_rule':'N_total * original_stratum_draw_count / (1536 * natural_stratum_windows)'}
CONFIGS = {
    'C0':{**COMMON,'objective':'local_only'},
    'Q':{**COMMON,'num_queries':16,'max_exons':128,'existence_logit_threshold':0.,
         'loss_coefficients':[1.,1.,1.,1.],'existence_reduction':'mean_over_known_queries',
         'query_collaboration':{'heads':4,'dropout':0.,'pre_norm':True,'residual':True,'FFN':False}}}


def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]


def validate_dataset(order, entries, windows, catalogs, old_draws):
    ids = [d['window_id'] for d in order]
    if (len(ids) != len(set(ids)) or set(ids) != set(windows) or
        len(entries) != len(windows) or
        len(old_draws) != 1536 or
        Counter(d['species'] for d in old_draws) != Counter({s:768 for s in ALLOW})):
        raise ValueError('Full natural TRAIN population once per pass and original strata required')
    strata = Counter((w['species'],w['seqid'],w['strand']) for w in windows.values())
    old_counts = Counter((d['species'],*d['stratum']) for d in old_draws)
    expected_strata = {(s,q,t) for s,qs in ALLOW.items() for q in qs for t in ('+','-')}
    if set(strata) != expected_strata or set(old_counts) != expected_strata:
        raise ValueError('Six TRAIN chromosomes on both strands required')
    ordered_windows = sorted(windows.values(),key=lambda w:(w['species'],w['seqid'],w['strand'],w['start']))
    for ordinal,(entry,w) in enumerate(zip(entries,ordered_windows)):
        if (entry['ordinal'] != ordinal or entry['window'] != {k:w[k] for k in GEOMETRY} or
            w['split'] != 'train' or w['seqid'] not in ALLOW[w['species']] or
            w['label_policy'] != 'CDS_assessable_primary_R4' or
            w['width'] != 24576 or w['end']-w['start'] != w['valid_bases'] or
            len(w['positive_ids']) > 16 or len(set(w['positive_ids'])) != len(w['positive_ids']) or
            not set(w['positive_ids']) <= set(w['overlapping_transcript_ids'])):
            raise ValueError('Full TRAIN feature/label geometry or capacity mismatch')
        for ident in w['overlapping_transcript_ids']:
            tx = catalogs[w['species']][ident]
            if (tx['split'] != 'train' or
                (tx['seqid'],tx['strand']) != (w['seqid'],w['strand']) or
                (ident in w['positive_ids'] and not (tx['training_primary'] and tx['eligible_positive']))):
                raise ValueError('TRAIN transcript identity or positive eligibility mismatch')
    expected_order = [w['id'] for w in ordered_windows]
    random.Random(0).shuffle(expected_order)
    if ids != expected_order:
        raise ValueError('Frozen full-window seed0 order changed')
    for d in order:
        w = windows[d['window_id']]; key = w['species'],w['seqid'],w['strand']
        expected = len(windows)*old_counts[key]/(1536*strata[key])
        weight = d['importance_weight_original_stratum_objective']
        if ((d['species'],*d['stratum']) != key or not math.isfinite(weight) or
            weight <= 0 or weight != expected):
            raise ValueError('Original expected stratum risk must be preserved')
    return {'windows':len(windows),
            'per_species':{s:{'windows':sum(w['species']==s for w in windows.values()),
                              'unique_positive_chains':len({tid for w in windows.values() if w['species']==s for tid in w['positive_ids']}),
                              'effective_weight_sum':sum(d['importance_weight_original_stratum_objective'] for d in order if d['species']==s)}
                           for s in ALLOW}}


def load_metadata():
    if (PLAN/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Complete full-TRAIN plan required')
    planned = json.loads((PLAN/'summary.json').read_text())
    if (not planned['completed'] or planned['feature_spec'] != FEATURE_SPEC or
        planned['natural_TRAIN_windows'] != N or planned['test_setaria_used'] is not False):
        raise ValueError('Matching fixed TRAIN-only plan identity required')
    entries = rows(PLAN/'feature_plan.jsonl'); order = rows(PLAN/'training_order.jsonl')
    windows = {}; catalogs = {}
    for species in ALLOW:
        local = list(train_records(LABELS/species/'windows_train.jsonl'))
        for w in local:
            if w['id'] in windows or w['species'] != species:
                raise ValueError('Duplicate or mismatched TRAIN window')
            windows[w['id']] = w
        needed = {tid for w in local for tid in w['overlapping_transcript_ids']}
        records = list(train_records(LABELS/species/'transcripts.jsonl'))
        selected = [t for t in records if t['id'] in needed]
        catalogs[species] = {t['id']:t for t in selected}
        if len(selected) != len(catalogs[species]) or set(catalogs[species]) != needed:
            raise ValueError('Full TRAIN transcript join incomplete or duplicated')
    old_draws = rows(Q.DRAWS)
    facts = validate_dataset(order,entries,windows,catalogs,old_draws)
    if facts['windows'] != N:
        raise ValueError('Full TRAIN window count changed')
    return order,entries,windows,catalogs,facts


def validate_feature_summary(s):
    if (not s['completed'] or s['feature_spec'] != FEATURE_SPEC or
        s['train_unique_windows'] != N or s['native_windows'] != 24021 or
        s['reused_existing_windows'] != 1491 or s['DEV_windows'] != 0 or
        s['optimizer_steps_this_run'] != 0 or s['test_setaria_used'] is not False or
        s['DEV_forward_used'] is not False):
        raise ValueError('Complete frozen full-TRAIN-only native features required')


def feature_index(entries):
    if (FEATURES/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Full native TRAIN materialization must finish before model execution')
    validate_feature_summary(json.loads((FEATURES/'summary.json').read_text()))
    records = rows(FEATURES/'index.jsonl'); index = {r['window_id']:r for r in records}
    if len(records) != N or len(index) != N:
        raise ValueError('Do not omit or duplicate natural TRAIN features')
    for r in entries:
        cached = index[r['window']['id']]
        if (cached['ordinal'] != r['ordinal'] or cached['split'] != 'train' or
            cached['artifact'] != r['artifact'] or cached['species'] != r['window']['species']):
            raise ValueError('Materialized features differ from the frozen TRAIN plan')
    return index


def new_model(arm, device):
    torch.manual_seed(COMMON['seed'])
    if arm == 'C0': model = M28Core('C0',COMMON['feature_dim'],COMMON['hidden'])
    elif arm == 'Q': model = CollaborativeGeneQueryModel(COMMON['feature_dim'],16,COMMON['hidden'])
    else: raise ValueError('Fixed matched arm required')
    return model.to(device)


def supervised_readout(arm, model, outputs, targets, n, weight):
    if arm == 'Q':
        result = joint_known_query_losses(model,outputs,targets,n,weight,CONFIGS['Q']['loss_coefficients'])
    elif arm == 'C0':
        loss = local_loss(outputs['segmentation_logits'][0],targets)
        result = {'total':loss*weight,'losses':{'local':loss},'protection_reasons':{}}
    else: raise ValueError('Fixed matched arm required')
    if not bool(torch.isfinite(result['total'])): raise ValueError('Nonfinite full-TRAIN loss')
    return result


def write_summary(out, report):
    with (out/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(report,allow_nan=False),flush=True)


def readiness():
    order,entries,windows,catalogs,facts = load_metadata()
    timing = Q.select_timing_draws(order,windows,catalogs)
    chosen = [{'window_id':r['draw']['window_id'],'species':r['species'],'category':r['category'],
               'positive_chains':r['positive'],'positive_exons':r['exons']} for r in timing]
    with (READY/'witnesses.jsonl').open('x') as f:
        for r in chosen: f.write(json.dumps(r)+'\n')
    write_summary(READY,{'completed':True,'experiment':'M28-FULL-TRAIN-READINESS',
                         'configurations':CONFIGS,'feature_spec':FEATURE_SPEC,'dataset':facts,
                         'fixed_worst_case_witnesses':chosen,'timing_prefix_windows':64,
                         'optimizer_steps':0,'DNA_features_checkpoints_read':False,
                         'DEV_used_for_selection':False,'test_setaria_used':False,
                         'planned_optimizer_steps_per_arm':STEPS,'protocol':PROTOCOL})


def require_readiness():
    if (READY/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Full TRAIN metadata and targeted tests must complete first')
    s = json.loads((READY/'summary.json').read_text())
    if not s['completed'] or s['configurations'] != CONFIGS or s['dataset']['windows'] != N:
        raise ValueError('Matching fixed full-TRAIN readiness required')


def smoke():
    require_readiness()
    order,entries,windows,catalogs,facts = load_metadata(); index = feature_index(entries)
    witnesses = rows(READY/'witnesses.jsonl')
    selected = [{'window_id':d['window_id'],'category':'fixed_order_prefix'} for d in order[:64]]+witnesses
    weights = {d['window_id']:d['importance_weight_original_stratum_objective'] for d in order}
    records = []; arms = []
    for arm in ('C0','Q'):
        model = new_model(arm,'cuda').train(); torch.cuda.reset_peak_memory_stats(); began = time.perf_counter()
        for r in selected:
            ident = r['window_id']; w = windows[ident]; started = time.perf_counter()
            cached = torch.load(ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
            features,dna = Q.window_inputs(cached,w,'cuda')
            targets = build_targets(w,catalogs[w['species']],cached['sequence'])
            model.zero_grad(set_to_none=True)
            outputs = model(features,dna,[w['valid_bases']])
            result = supervised_readout(arm,model,outputs,targets,w['valid_bases'],weights[ident])
            result['total'].backward()
            norm = torch.nn.utils.clip_grad_norm_(model.parameters(),1.,error_if_nonfinite=True)
            torch.cuda.synchronize()
            records.append({'arm':arm,**r,'losses':{k:float(v.detach()) for k,v in result['losses'].items()},
                            'positive_chains':len(targets['positive_chains']),
                            'positive_exons':sum(map(len,targets['positive_chains'])),
                            'gradient_norm_before_clip':float(norm),'seconds':time.perf_counter()-started})
            del cached,features,dna,targets,outputs,result
        prefix = [r['seconds'] for r in records if r['arm']==arm and r['category']=='fixed_order_prefix']
        arms.append({'arm':arm,'backward_windows':len(selected),'optimizer_steps':0,
                     'trainable_parameters':sum(p.numel() for p in model.parameters()),
                     'fixed_prefix_mean_seconds':sum(prefix)/len(prefix),
                     'projected_fit_seconds_excluding_AdamW_not_runtime_guarantee':STEPS*sum(prefix)/len(prefix),
                     'wall_seconds':time.perf_counter()-began,
                     'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                     'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()})
        del model; torch.cuda.empty_cache()
    write_summary(SMOKE,{'completed':True,'experiment':'M28-FULL-TRAIN-GPU-SMOKE',
                        'configurations':CONFIGS,'feature_spec':FEATURE_SPEC,'dataset':facts,
                        'arms':arms,'records':records,'optimizer_steps':0,
                        'DEV_used_for_selection':False,'test_setaria_used':False,
                        'scope':'readiness/throughput only; no fit, accuracy or scientific PASS'})


def fit(arm):
    require_readiness()
    if (SMOKE/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Native full-TRAIN forward/backward smoke must complete before fit')
    tested = json.loads((SMOKE/'summary.json').read_text())
    if not tested['completed'] or tested['configurations'] != CONFIGS or tested['optimizer_steps'] != 0:
        raise ValueError('Matching zero-update native smoke required')
    order,entries,windows,catalogs,facts = load_metadata(); index = feature_index(entries)
    out = ROOT/'outputs'/('M28-FULL-TRAIN-'+arm+'-FIT')
    model = new_model(arm,'cuda').train()
    optimizer = torch.optim.AdamW(model.parameters(),lr=COMMON['lr'],weight_decay=COMMON['weight_decay'],
                                 betas=tuple(COMMON['betas']),eps=COMMON['eps'])
    torch.cuda.reset_peak_memory_stats(); began = time.perf_counter(); step = 0; epochs = []
    with (out/'training.jsonl').open('x') as log:
        for epoch in range(PASSES):
            sums = Counter(); counts = {s:Counter() for s in ALLOW}; epoch_began = time.perf_counter()
            for ordinal,d in enumerate(order):
                ident = d['window_id']; w = windows[ident]; started = time.perf_counter()
                cached = torch.load(ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
                features,dna = Q.window_inputs(cached,w,'cuda')
                targets = build_targets(w,catalogs[w['species']],cached['sequence'])
                weight = d['importance_weight_original_stratum_objective']
                optimizer.zero_grad(set_to_none=True)
                outputs = model(features,dna,[w['valid_bases']])
                result = supervised_readout(arm,model,outputs,targets,w['valid_bases'],weight)
                result['total'].backward()
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(),COMMON['clip_norm'],error_if_nonfinite=True)
                optimizer.step(); step += 1; torch.cuda.synchronize()
                values = {k:float(v.detach()) for k,v in result['losses'].items()}
                label_counts = (Counter(result['query']['query_labels']) if arm=='Q' else {})
                record = {'step':step,'epoch':epoch+1,'window_ordinal':ordinal,'window_id':ident,
                          'species':w['species'],'importance_weight':weight,'losses':values,
                          'weighted_total':float(result['total'].detach()),'gradient_norm_before_clip':float(norm),
                          'positive_chains':len(targets['positive_chains']),
                          'query_label_counts':{str(k):v for k,v in label_counts.items()},
                          'protection_reasons':result['protection_reasons'],'seconds':time.perf_counter()-started}
                log.write(json.dumps(record,allow_nan=False)+'\n')
                sums.update(dict(values,weighted_total=record['weighted_total']))
                counts[w['species']].update({'windows':1,'positive_chains':record['positive_chains'],
                                             **{'query_'+str(k):v for k,v in label_counts.items()}})
                if step % 64 == 0:
                    log.flush(); print(json.dumps({'arm':arm,'step':step,'total':STEPS,'seconds':time.perf_counter()-began}),flush=True)
                del cached,features,dna,targets,outputs,result
            checkpoint = out/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists(): raise FileExistsError(checkpoint)
            torch.save({'arm':arm,'configuration':CONFIGS[arm],'step':step,'revision':Q.F.REVISION,
                        'model':model.state_dict(),'optimizer':optimizer.state_dict(),'protocol':PROTOCOL},checkpoint)
            epochs.append({'epoch':epoch+1,'steps':N,'mean_losses':{k:v/N for k,v in sums.items()},
                           'supervision_by_species':{s:dict(c) for s,c in counts.items()},
                           'seconds':time.perf_counter()-epoch_began,'checkpoint':str(checkpoint.relative_to(ROOT))})
    if step != STEPS: raise ValueError('Incomplete fixed full-TRAIN fit')
    write_summary(out,{'completed':True,'arm':arm,'configuration':CONFIGS[arm],'feature_spec':FEATURE_SPEC,
                       'feature_revision':Q.F.REVISION,'optimizer_steps':step,'dataset':facts,'epochs':epochs,
                       'trainable_parameters':sum(p.numel() for p in model.parameters()),'protocol':PROTOCOL,
                       'final_checkpoint':str(checkpoint.relative_to(ROOT)),'main_evaluation':'last_step_076536_only',
                       'initialization':'seed0 from scratch; identical shared/segmentation initial weights across arms',
                       'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                       'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
                       'DEV_used_for_selection':False,'test_setaria_used':False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('readiness','smoke','fit'),required=True)
    parser.add_argument('--arm',choices=('C0','Q'))
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'fit':
        if args.arm is None: parser.error('--arm is required for fit')
        fit(args.arm)
    else: {'readiness':readiness,'smoke':smoke}[args.mode]()
