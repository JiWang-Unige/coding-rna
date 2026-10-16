#!/usr/bin/env python3
"""Prepare every existing natural TRAIN window; no DNA/model execution."""
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
import feature_smoke as F
import manifest as M
from query_training_census import ALLOW, LABELS, DRAWS, train_records

ROOT = M.C.ROOT
OUT = ROOT/'outputs/M28-FULL-TRAIN-PLAN-R2'
FEATURES = ROOT/'outputs/M28-FULL-TRAIN-FEATURES'
OLD_FEATURES = ROOT/'outputs/M28-PILOT-R4/features'
CENSUS = ROOT/'outputs/M28-A-QUERY-TRAIN-CENSUS'
SHARDS = 16
GEOMETRY = ('id','species','split','seqid','strand','start','end','width','valid_bases','pad_right_after_orientation')
FEATURE_SPEC = {'model':F.MODEL,'revision':F.REVISION,'window_bp':24576,'block_bp':6144,
                'feature_shape':[4096,2048],'dtype':'torch.bfloat16','attn_implementation':'sdpa'}


def build_plan(windows, cached_index, draws):
    by_id = {}
    for w in windows:
        n = w['valid_bases']
        if (w['split'] != 'train' or w['species'] not in ALLOW or
            w['seqid'] not in ALLOW[w['species']] or w['strand'] not in ('+','-') or
            w['label_policy'] != 'CDS_assessable_primary_R4' or
            w['width'] != M.WIDTH or not 0 < n <= M.WIDTH or w['end']-w['start'] != n or
            w['pad_right_after_orientation'] != M.WIDTH-n or
            w['id'] != w['species']+'|'+w['seqid']+'|'+str(w['start'])+'|'+w['strand'] or
            w['id'] in by_id):
            raise ValueError('Unique allowed TRAIN R4 geometry required')
        if len(w['positive_ids']) > 16 or len(set(w['positive_ids'])) != len(w['positive_ids']):
            raise ValueError('Do not drop positive chains to fit Q16')
        by_id[w['id']] = w
    strata = Counter((w['species'],w['seqid'],w['strand']) for w in windows)
    draw_counts = Counter()
    for d in draws:
        w = by_id[d['window_id']]
        if (w['species'],w['seqid'],w['strand']) != (d['species'],*d['stratum']):
            raise ValueError('Original TRAIN stratum identity mismatch')
        old_weight = float(d['importance_weight_uniform_within_stratum'])
        if not math.isfinite(old_weight) or old_weight <= 0:
            raise ValueError('Original finite positive importance required')
        draw_counts[w['species'],w['seqid'],w['strand']] += 1
    if (set(cached_index) != {d['window_id'] for d in draws} or
        not set(cached_index) <= set(by_id) or set(draw_counts) != set(strata)):
        raise ValueError('Original cached draws and all twelve TRAIN strata required')
    entries = []; weights = {}
    for ordinal,w in enumerate(sorted(windows,key=lambda r:(r['species'],r['seqid'],r['strand'],r['start']))):
        key = w['species'],w['seqid'],w['strand']
        weight = len(windows)*draw_counts[key]/(len(draws)*strata[key])
        weights[w['id']] = weight
        old = cached_index.get(w['id']); shard = ordinal % SHARDS
        artifact = (old['artifact'] if old else
                    str((FEATURES/('shard_'+str(shard).zfill(2))/('window_'+str(ordinal).zfill(5)+'.pt')).relative_to(ROOT)))
        entries.append({'ordinal':ordinal,'shard':shard,
                        'window':{k:w[k] for k in GEOMETRY},
                        'artifact':artifact,'reused_existing':old is not None})
    order = list(entries); random.Random(0).shuffle(order)
    training = [{'window_id':r['window']['id'],'species':r['window']['species'],
                 'stratum':[r['window']['seqid'],r['window']['strand']],
                 'importance_weight_original_stratum_objective':weights[r['window']['id']]} for r in order]
    sentinels = {}
    for r in entries:
        w = r['window']; key = w['species'],w['seqid'],w['strand']
        if r['reused_existing']: sentinels.setdefault(key,r)
    if set(sentinels) != set(strata):
        raise ValueError('Each TRAIN stratum needs an existing feature replay sentinel')
    return entries,training,[sentinels[k] for k in sorted(sentinels)]


def main():
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (OLD_FEATURES,CENSUS)):
        raise ValueError('Existing cache and TRAIN census must be complete')
    source = json.loads((OLD_FEATURES/'summary.json').read_text())
    census = json.loads((CENSUS/'summary.json').read_text())
    if (source['model'] != F.MODEL or source['revision'] != F.REVISION or
        source['train_unique_windows'] != 1491 or source['DEV_windows'] != 8690):
        raise ValueError('Original frozen released feature identity required')
    windows = [w for species in ALLOW for w in train_records(LABELS/species/'windows_train.jsonl')]
    draws = [json.loads(line) for line in DRAWS.read_text().splitlines() if line.strip()]
    old_rows = list(train_records(OLD_FEATURES/'index.jsonl'))
    old_index = {r['window_id']:r for r in old_rows}
    if (len(draws) != 1536 or len({d['window_id'] for d in draws}) != 1491 or
        Counter(d['species'] for d in draws) != Counter({s:768 for s in ALLOW}) or
        len(old_rows) != len(old_index) or len(old_index) != 1491 or len(windows) != 25512):
        raise ValueError('Unchanged original draws and full natural TRAIN population required')
    entries,training,sentinels = build_plan(windows,old_index,draws)
    if len(sentinels) != 12:
        raise ValueError('Six TRAIN chromosomes on both strands required')
    lengths = {s:{} for s in ALLOW}; per_species = {}
    ids = {d['window_id'] for d in draws}
    for species in ALLOW:
        ws = [w for w in windows if w['species'] == species]
        if {w['seqid'] for w in ws} != ALLOW[species]:
            raise ValueError('Frozen TRAIN chromosome set changed')
        for seqid in sorted(ALLOW[species]):
            rows = [w for w in ws if w['seqid'] == seqid]
            length = max(w['end'] for w in rows); lengths[species][seqid] = length
            expected = [(a,min(a+M.WIDTH,length)) for a in M.starts(length)]
            for strand in ('+','-'):
                actual = sorted((w['start'],w['end']) for w in rows if w['strand'] == strand)
                if actual != expected:
                    raise ValueError('Full natural TRAIN window grid changed')
        seen = {tid for w in ws for tid in w['positive_ids']}
        drawn = {tid for w in ws if w['id'] in ids for tid in w['positive_ids']}
        previous = census['views'][species]
        if (len(ws) != previous['all_natural_TRAIN_windows']['counts']['windows'] or
            len(seen) != previous['all_natural_TRAIN_windows']['unique_positive_transcript_ids_exposed'] or
            len(drawn) != previous['original_768_draws_with_replacement']['unique_positive_transcript_ids_exposed']):
            raise ValueError('Full TRAIN exposure differs from completed R4 census')
        per_species[species] = {'natural_windows':len(ws),'original_unique_cached_windows':sum(w['id'] in ids for w in ws),
                               'natural_unique_positive_chains':len(seen),'original_unique_positive_chains':len(drawn),
                               'maximum_positive_chains_per_window':max(len(w['positive_ids']) for w in ws),
                               'effective_weight_sum_one_full_pass':sum(
                                   d['importance_weight_original_stratum_objective'] for d in training if d['species'] == species)}
    for r in entries:
        if r['reused_existing'] and not (ROOT/r['artifact']).is_file():
            raise FileNotFoundError(r['artifact'])
    report = {'experiment':'M28-FULL-TRAIN-PLAN','completed':True,'feature_spec':FEATURE_SPEC,
              'label_source':str(LABELS.resolve().relative_to(ROOT)),'old_feature_source':str(OLD_FEATURES.relative_to(ROOT)),
              'TRAIN_allowlists':{s:sorted(qs) for s,qs in ALLOW.items()},'TRAIN_lengths':lengths,
              'natural_TRAIN_windows':len(entries),'reused_existing_windows':len(old_index),
              'missing_native_feature_windows':sum(not r['reused_existing'] for r in entries),
              'shards':SHARDS,'per_species':per_species,'training_order_seed':0,
              'importance_rule':'N_total * original_stratum_draw_count / (1536 * natural_stratum_windows)',
              'preserves_original_expected_stratum_risk_not_identical_optimizer_updates':True,
              'planned_passes_for_new_matched_exposure_study':3,'planned_new_optimizer_steps':3*len(entries),
              'optimizer_steps_this_run':0,'Q16_capacity_sufficient_for_full_TRAIN_metadata':True,
              'DEV_reference_parsed':False,'DNA_features_checkpoints_read':False,
              'DEV_forward_used':False,'test_setaria_used':False,
              'limitations':'metadata coverage only; not performance, complete biological labels, or independent generalization'}
    OUT.mkdir(parents=True,exist_ok=True)
    for name,records in (('feature_plan.jsonl',entries),('training_order.jsonl',training),('replay_sentinels.jsonl',sentinels)):
        with (OUT/name).open('x') as f:
            for row in records: f.write(json.dumps(row,allow_nan=False)+'\n')
    with (OUT/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(report,indent=2,allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
