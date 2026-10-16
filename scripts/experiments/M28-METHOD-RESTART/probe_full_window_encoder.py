#!/usr/bin/env python3
"""Eight fixed TRAIN windows: native four-block replay and whole-window encoding.

Engineering only. No head, optimizer, reference catalog, DEV or test execution.
"""
import json
import re
import time
from collections import Counter
import torch
import cache_full_train as C

OUT = C.ROOT/'outputs/M28-A10-ENCODER-CONTEXT-PREFLIGHT'
SELECTED = C.ROOT/'outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC'
FULL_SPEC = {**C.FEATURE_SPEC, 'block_bp':24576,
             'encoding_mode':'one_4096_token_forward'}


def selected_inputs():
    if (SELECTED/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Original eight-window TRAIN diagnostic must be closed')
    predictions = C.rows(SELECTED/'predictions.jsonl')
    ids = {r['window']['id'] for r in predictions}
    if (len(predictions) != 8 or len(ids) != 8 or
        Counter(r['window']['species'] for r in predictions) !=
        Counter({s:4 for s in C.ALLOW})):
        raise ValueError('Reuse exactly the original eight TRAIN windows')
    _, entries = C.load_plan()
    matches = [r for r in entries if r['window']['id'] in ids]
    if len(matches) != 8 or {r['window']['id'] for r in matches} != ids:
        raise ValueError('Each selected TRAIN window needs exactly one cache entry')
    plan = {r['window']['id']:r for r in matches}
    selected = []
    for pred in predictions:
        w = pred['window']; entry = plan[w['id']]
        if (any(w[k] != entry['window'][k] for k in C.GEOMETRY) or
            w['split'] != 'train' or w['seqid'] not in C.ALLOW[w['species']] or
            w['valid_bases'] != 24576):
            raise ValueError('Unchanged full-width allowed TRAIN geometry required')
        # Only geometry and oriented DNA are retained; no query or label fields.
        selected.append((entry, pred['sequence']))
    return selected


def memory_receipt(seconds):
    return {'seconds':seconds,
            'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
            'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()}


@torch.inference_mode()
def main():
    selected = selected_inputs()
    began = time.perf_counter()
    backbone, tokenizer = C.load_native()
    if int(backbone.config.max_position_embeddings) < 4096:
        raise ValueError('Pinned native encoder cannot support 4096 positions')
    if any(p.requires_grad for p in backbone.parameters()):
        raise ValueError('Frozen encoder required')
    torch.cuda.synchronize()
    load_seconds = time.perf_counter()-began
    records = []
    with (OUT/'records.jsonl').open('x') as log:
        for ordinal, (entry, source_dna) in enumerate(selected):
            w = entry['window']
            old = torch.load(C.ROOT/entry['artifact'],map_location='cpu',weights_only=True)
            C.Q.validate_cached(old,w)
            if old['sequence'] != source_dna:
                raise ValueError('Original selected oriented DNA differs from feature cache')
            # Match the existing TRAIN-cache normalizer exactly, not a new caller.
            sequence = re.sub('[^ACGT]','N',old['sequence'])
            whole = tokenizer(sequence,add_special_tokens=False,
                              return_tensors='pt',truncation=False)['input_ids']
            blocks = [tokenizer(sequence[a:a+C.F.BLOCK],add_special_tokens=False,
                                return_tensors='pt',truncation=False)['input_ids']
                      for a in range(0,len(sequence),C.F.BLOCK)]
            if (tuple(whole.shape) != (1,4096) or
                any(tuple(ids.shape) != (1,1024) for ids in blocks) or
                not torch.equal(torch.cat(blocks,dim=1),whole)):
                raise ValueError('Same 4096 six-mer IDs required in both encoding modes')

            torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
            started = time.perf_counter()
            replay = C.encode_window(backbone,tokenizer,sequence,'cuda')[0]
            replay = replay.to('cpu',dtype=torch.bfloat16).contiguous()
            torch.cuda.synchronize()
            control = memory_receipt(time.perf_counter()-started)
            C.exact_replay(old['features'],replay)
            del old

            ids = whole.to('cuda')
            torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
            started = time.perf_counter()
            hidden = backbone(input_ids=ids,attention_mask=torch.ones_like(ids)).last_hidden_state
            if tuple(hidden.shape) != (1,4096,2048) or not bool(torch.isfinite(hidden).all()):
                raise ValueError('Full-window hidden output must be finite and untruncated')
            features = hidden[0].to('cpu',dtype=torch.bfloat16).contiguous()
            del hidden, ids
            torch.cuda.synchronize()
            full = memory_receipt(time.perf_counter()-started)
            difference = (features.float()-replay.float()).abs()
            drift = {'mean_abs_difference':float(difference.mean()),
                     'max_abs_difference':float(difference.max())}
            del replay, difference
            path = OUT/('window_'+str(ordinal).zfill(2)+'.pt')
            with path.open('xb') as handle:
                torch.save({'features':features,'sequence':source_dna,'window':w,
                            'model':C.F.MODEL,'revision':C.F.REVISION,
                            'feature_spec':FULL_SPEC},handle)
            row = {'window':w,'source_feature_artifact':entry['artifact'],
                   'artifact':str(path.relative_to(C.ROOT)), 'bytes':path.stat().st_size,
                   'token_shape':[1,4096], 'feature_shape':list(features.shape),
                   'four_block_replay_exact':True, 'four_block':control,
                   'whole_window':full, 'feature_drift_not_accuracy':drift}
            log.write(json.dumps(row,allow_nan=False)+'\n'); log.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
            records.append(row)
            del features
    C.write_summary(OUT,{
        'experiment':'M28-A10-ENCODER-CONTEXT-PREFLIGHT', 'completed':True,
        'scope':'eight inherited metadata-selected TRAIN windows; engineering only',
        'source_selection':str((SELECTED/'predictions.jsonl').relative_to(C.ROOT)),
        'control_feature_spec':C.FEATURE_SPEC, 'whole_window_feature_spec':FULL_SPEC,
        'windows':len(records), 'per_species':dict(Counter(r['window']['species'] for r in records)),
        'four_block_replay_exact':True, 'records':records,
        'stored_feature_bytes':sum(r['bytes'] for r in records),
        'GPU_name':torch.cuda.get_device_name(),
        'encoder_load_seconds':load_seconds, 'wall_seconds':time.perf_counter()-began,
        'optimizer_steps_this_run':0, 'backbone_trainable_parameters':0,
        'head_checkpoint_read':False, 'head_forward_used':False,
        'train_reference_catalog_read':False, 'reference_labels_used_in_encoding':False,
        'DEV_forward_used':False, 'test_setaria_used':False,
        'limitations':'Both context attention and position treatment change. Features are not interchangeable with the four-block cache. Control-first single execution includes checks and materialization; no speedup, accuracy, causality or scientific PASS claim.'})


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
