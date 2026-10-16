#!/usr/bin/env python3
"""Whole-window features for exactly the original 1491 unique TRAIN draws.

No reference catalog, supervised head, optimizer, DEV or test reads.
"""
import json
import math
import re
import time
from collections import Counter
import torch
import cache_full_train as C
from probe_full_window_encoder import FULL_SPEC, OUT as PREFLIGHT

OUT = C.ROOT/'outputs/M28-A10-DRAW-CONTEXT-CACHE'


def load_draw_plan():
    if (PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Whole-window native engineering preflight must be complete')
    checked = json.loads((PREFLIGHT/'summary.json').read_text())
    if (not checked['completed'] or checked['windows'] != 8 or
        not checked['four_block_replay_exact'] or
        checked['whole_window_feature_spec'] != FULL_SPEC):
        raise ValueError('Matching eight-window native preflight required')
    draws = C.rows(C.Q.DRAWS)
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    if (len(draws) != 1536 or len(ids) != 1491 or
        Counter(d['species'] for d in draws) != Counter({s:768 for s in C.ALLOW})):
        raise ValueError('Original TRAIN draws, not expanded TRAIN or DEV, required')
    _, plan = C.load_plan()
    needed = set(ids)
    matches = [r for r in plan if r['window']['id'] in needed]
    if len(matches) != 1491 or {r['window']['id'] for r in matches} != needed:
        raise ValueError('Each original unique draw requires one TRAIN plan entry')
    by_id = {r['window']['id']:r for r in matches}
    for d in draws:
        row = by_id[d['window_id']]; w = row['window']
        weight = float(d['importance_weight_uniform_within_stratum'])
        if (w['split'] != 'train' or w['seqid'] not in C.ALLOW[w['species']] or
            (w['species'],w['seqid'],w['strand']) != (d['species'],*d['stratum']) or
            not row['reused_existing'] or not math.isfinite(weight) or weight <= 0):
            raise ValueError('Original TRAIN source, stratum and importance required')
    shared = {r['window']['id']:r for r in checked['records'] if r['window']['id'] in needed}
    if not shared:
        raise ValueError('At least one existing native whole-window replay witness required')
    return draws, [by_id[ident] for ident in ids], shared


@torch.inference_mode()
def encode_full_window(backbone, tokenizer, sequence):
    if not 0 < len(sequence) <= 24576:
        raise ValueError('Valid oriented TRAIN DNA required')
    padded = re.sub('[^ACGT]','N',sequence)+'N'*(24576-len(sequence))
    whole = tokenizer(padded,add_special_tokens=False,
                      return_tensors='pt',truncation=False)['input_ids']
    blocks = [tokenizer(padded[a:a+C.F.BLOCK],add_special_tokens=False,
                        return_tensors='pt',truncation=False)['input_ids']
              for a in range(0,len(padded),C.F.BLOCK)]
    if (tuple(whole.shape) != (1,4096) or
        any(tuple(ids.shape) != (1,1024) for ids in blocks) or
        not torch.equal(torch.cat(blocks,dim=1),whole)):
        raise ValueError('Original four-block and whole-window six-mer IDs must agree')
    ids = whole.to('cuda')
    hidden = backbone(input_ids=ids,attention_mask=torch.ones_like(ids)).last_hidden_state
    if tuple(hidden.shape) != (1,4096,2048) or not bool(torch.isfinite(hidden).all()):
        raise ValueError('Finite untruncated whole-window representation required')
    return hidden[0].to('cpu',dtype=torch.bfloat16).contiguous()


@torch.inference_mode()
def main():
    draws, selected, witnesses = load_draw_plan()
    # Preserve the original list/order/weights for a later matched fit.
    with (OUT/'draws.jsonl').open('x') as handle:
        for d in draws:
            handle.write(json.dumps(d,allow_nan=False)+'\n')
    backbone, tokenizer = C.load_native()
    if (int(backbone.config.max_position_embeddings) < 4096 or
        any(p.requires_grad for p in backbone.parameters())):
        raise ValueError('Same frozen context-capable native encoder required')
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    total_bytes = 0; replayed = []; seconds = []; per_species = Counter()
    with (OUT/'index.jsonl').open('x') as index:
        for ordinal, entry in enumerate(selected):
            w = entry['window']
            # Only DNA/geometry/shape are needed, not the old feature payload.
            old = torch.load(C.ROOT/entry['artifact'],map_location='cpu',
                             weights_only=True,mmap=True)
            C.Q.validate_cached(old,w)
            sequence = old['sequence']; del old
            torch.cuda.synchronize(); started = time.perf_counter()
            features = encode_full_window(backbone,tokenizer,sequence)
            torch.cuda.synchronize(); elapsed = time.perf_counter()-started
            if w['id'] in witnesses:
                record = witnesses[w['id']]
                prior = torch.load(C.ROOT/record['artifact'],map_location='cpu',weights_only=True)
                if (prior['sequence'] != sequence or prior['window'] != w or
                    prior['feature_spec'] != FULL_SPEC):
                    raise ValueError('Existing whole-window replay identity differs')
                C.exact_replay(prior['features'],features)
                replayed.append(w['id']); del prior
            path = OUT/('window_'+str(ordinal).zfill(4)+'.pt')
            with path.open('xb') as handle:
                torch.save({'features':features,'sequence':sequence,'window':w,
                            'model':C.F.MODEL,'revision':C.F.REVISION,
                            'feature_spec':FULL_SPEC},handle)
            size = path.stat().st_size; total_bytes += size
            row = {'ordinal':ordinal,'window_id':w['id'],'species':w['species'],
                   'split':'train','window':w,'artifact':str(path.relative_to(C.ROOT)),
                   'source_feature_artifact':entry['artifact'],'bytes':size,
                   'encoder_seconds':elapsed}
            index.write(json.dumps(row,allow_nan=False)+'\n')
            seconds.append(elapsed); per_species[w['species']] += 1
            if (ordinal+1) % 64 == 0:
                index.flush()
                print(json.dumps({'windows':ordinal+1,'total':1491,
                                  'seconds':time.perf_counter()-began}),flush=True)
            del features
    if set(replayed) != set(witnesses):
        raise ValueError('Every shared original native witness must replay exactly')
    C.write_summary(OUT,{
        'experiment':'M28-A10-DRAW-CONTEXT-CACHE','completed':True,
        'feature_spec':FULL_SPEC,'unique_TRAIN_windows':len(selected),
        'source_draws':str(C.Q.DRAWS.relative_to(C.ROOT)), 'draws':len(draws),
        'per_species_unique_windows':dict(per_species),
        'draw_order_weights_unchanged':True,'features_finite':True,
        'full_window_replay_exact':True,'replayed_window_ids':replayed,
        'stored_feature_bytes':total_bytes,'wall_loop_seconds':time.perf_counter()-began,
        'encoder_seconds_min':min(seconds),'encoder_seconds_max':max(seconds),
        'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
        'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
        'optimizer_steps_this_run':0,'backbone_trainable_parameters':0,
        'head_checkpoint_read':False,'train_reference_catalog_read':False,
        'DEV_forward_used':False,'test_setaria_used':False,
        'limitations':'TRAIN feature preparation only, not fitted-head performance. Never mix with the four-block cache or reuse its trained heads as a clean comparator.'})


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
