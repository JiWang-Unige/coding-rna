#!/usr/bin/env python3
"""TRAIN-only native feature replay, independent missing-feature shards and final index."""
import argparse
import json
import re
import time
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
import run_query as Q
from annotate_query import encode_window
from prepare_full_train import OUT as PLAN, FEATURES, SHARDS, GEOMETRY, FEATURE_SPEC, ALLOW

ROOT = M.C.ROOT
REPLAY = ROOT/'outputs/M28-FULL-TRAIN-REPLAY'


def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]


def load_plan():
    if (PLAN/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Complete TRAIN-only preparation required')
    summary = json.loads((PLAN/'summary.json').read_text())
    entries = rows(PLAN/'feature_plan.jsonl')
    if (not summary['completed'] or summary['feature_spec'] != FEATURE_SPEC or
        summary['natural_TRAIN_windows'] != 25512 or len(entries) != 25512 or
        summary['missing_native_feature_windows'] != 24021 or summary['shards'] != SHARDS or
        summary['test_setaria_used'] is not False):
        raise ValueError('Complete fixed full-TRAIN feature plan required')
    for i,r in enumerate(entries):
        w = r['window']
        if (r['ordinal'] != i or r['shard'] != i % SHARDS or w['split'] != 'train' or
            w['seqid'] not in ALLOW[w['species']] or set(w) != set(GEOMETRY)):
            raise ValueError('TRAIN-only feature geometry and fixed shard assignment required')
    return summary,entries


def read_sequences(summary, entries):
    needed = {}
    for r in entries:
        w = r['window']; needed.setdefault(w['species'],{})[w['seqid']] = summary['TRAIN_lengths'][w['species']][w['seqid']]
    return {species:M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',lengths)
            for species,lengths in needed.items()}


def load_native():
    from transformers import AutoModelForTokenClassification, AutoTokenizer
    tokenizer = AutoTokenizer.from_pretrained(F.MODEL,revision=F.REVISION,trust_remote_code=True,local_files_only=True)
    full = AutoModelForTokenClassification.from_pretrained(
        F.MODEL,revision=F.REVISION,trust_remote_code=True,local_files_only=True,
        attn_implementation='sdpa',torch_dtype=torch.bfloat16)
    if int(getattr(full,'k',getattr(full.config,'k',0))) != 6:
        raise ValueError('Released six-mer encoder identity changed')
    backbone = full.model; del full
    backbone.eval().requires_grad_(False).to('cuda')
    return backbone,tokenizer


@torch.inference_mode()
def native_row(backbone, tokenizer, w, sequences):
    raw = sequences[w['species']][w['seqid']][w['start']:w['end']]
    sequence = raw if w['strand'] == '+' else M.rc(raw)
    # Match the old cache: normalized encoder tokens, original oriented DNA stored.
    encoded = encode_window(backbone,tokenizer,re.sub('[^ACGT]','N',sequence),'cuda')
    features = encoded[0].to('cpu',dtype=torch.bfloat16).contiguous()
    return {'features':features,'sequence':sequence,'window':{k:w[k] for k in GEOMETRY},
            'model':F.MODEL,'revision':F.REVISION}


def exact_replay(old, new):
    if old.dtype != torch.bfloat16 or new.dtype != torch.bfloat16 or old.shape != new.shape:
        raise ValueError('Identical BF16 feature geometry required')
    if not torch.equal(old,new):
        error = float((old.float()-new.float()).abs().max())
        raise ValueError('Native TRAIN features differ from old cache; max_abs_error='+str(error))
    return 0.


def write_summary(out, report):
    with (out/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(report,indent=2,allow_nan=False),flush=True)


def replay():
    summary,entries = load_plan()
    sentinels = rows(PLAN/'replay_sentinels.jsonl')
    if len(sentinels) != 12:
        raise ValueError('Twelve fixed TRAIN replay sentinels required')
    sequences = read_sequences(summary,sentinels); backbone,tokenizer = load_native()
    began = time.perf_counter(); records = []
    torch.cuda.reset_peak_memory_stats()
    for r in sentinels:
        w = r['window']; old = torch.load(ROOT/r['artifact'],map_location='cpu',weights_only=True)
        Q.validate_cached(old,w)
        new = native_row(backbone,tokenizer,w,sequences)
        if old['sequence'] != new['sequence']:
            raise ValueError('Native TRAIN oriented DNA differs from old cache')
        error = exact_replay(old['features'],new['features'])
        records.append({'window_id':w['id'],'artifact':r['artifact'],'features_exact':True,'max_abs_error':error})
        del old,new
    write_summary(REPLAY,{'experiment':'M28-FULL-TRAIN-REPLAY','completed':True,'feature_spec':FEATURE_SPEC,
                          'windows':12,'features_exact':True,'maximum_abs_error':0.,'records':records,
                          'optimizer_steps_this_run':0,'train_annotation_parsed_by_encoder':False,
                          'reference_labels_used_in_encoding':False,'DEV_forward_used':False,'test_setaria_used':False,
                          'wall_seconds':time.perf_counter()-began,
                          'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                          'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()})


def shard(number):
    if not 0 <= number < SHARDS:
        raise ValueError('Fixed full-TRAIN shard index required')
    summary,entries = load_plan()
    if (REPLAY/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Native replay must complete before mixing old/new features')
    replayed = json.loads((REPLAY/'summary.json').read_text())
    if not replayed['features_exact'] or replayed['feature_spec'] != FEATURE_SPEC:
        raise ValueError('Exact released encoder replay identity required')
    selected = [r for r in entries if r['shard'] == number]
    sequences = read_sequences(summary,selected); backbone,tokenizer = load_native()
    out = FEATURES/('shard_'+str(number).zfill(2))
    began = time.perf_counter(); native = reused = stored_bytes = 0
    torch.cuda.reset_peak_memory_stats()
    with (out/'index.jsonl').open('x') as index:
        for ordinal,r in enumerate(selected):
            w = r['window']; artifact = ROOT/r['artifact']
            if r['reused_existing']:
                reused += 1
            else:
                if artifact.exists(): raise FileExistsError(artifact)
                cached = native_row(backbone,tokenizer,w,sequences)
                Q.validate_cached(cached,w)
                torch.save(cached,artifact); del cached
                native += 1
            size = artifact.stat().st_size; stored_bytes += size
            index.write(json.dumps({'ordinal':r['ordinal'],'window_id':w['id'],'species':w['species'],'split':'train',
                                    'artifact':r['artifact'],'bytes':size,'reused_existing':r['reused_existing']})+'\n')
            if (ordinal+1) % 64 == 0:
                index.flush(); print(json.dumps({'shard':number,'windows':ordinal+1,'total':len(selected),
                                                'seconds':time.perf_counter()-began}),flush=True)
    write_summary(out,{'experiment':'M28-FULL-TRAIN-FEATURE-SHARD','completed':True,'feature_spec':FEATURE_SPEC,
                       'shard':number,'windows':len(selected),'native_windows':native,'reused_existing_windows':reused,
                       'referenced_bytes':stored_bytes,'optimizer_steps_this_run':0,
                       'backbone_trainable_parameters':sum(p.numel() for p in backbone.parameters() if p.requires_grad),
                       'train_annotation_parsed_by_encoder':False,'reference_labels_used_in_encoding':False,
                       'DEV_forward_used':False,'test_setaria_used':False,'wall_seconds':time.perf_counter()-began,
                       'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                       'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()})


def finalize():
    summary,plan = load_plan(); merged = {}; summaries = []
    for number in range(SHARDS):
        out = FEATURES/('shard_'+str(number).zfill(2))
        if (out/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('All sixteen feature shards must be terminal complete')
        s = json.loads((out/'summary.json').read_text())
        if (not s['completed'] or s['feature_spec'] != FEATURE_SPEC or s['shard'] != number or
            s['optimizer_steps_this_run'] != 0 or s['backbone_trainable_parameters'] != 0 or
            s['test_setaria_used'] is not False):
            raise ValueError('Fixed TRAIN-only frozen feature shard identity mismatch')
        local = rows(out/'index.jsonl')
        if len(local) != s['windows']: raise ValueError('Partial feature shard index')
        if (sum(not r['reused_existing'] for r in local) != s['native_windows'] or
            sum(r['reused_existing'] for r in local) != s['reused_existing_windows']):
            raise ValueError('Feature shard native/reuse counts differ from materialized index')
        for r in local:
            i = r['ordinal']; expected = plan[i]
            if (i in merged or expected['shard'] != number or
                r['window_id'] != expected['window']['id'] or r['split'] != 'train' or
                r['artifact'] != expected['artifact'] or r['reused_existing'] != expected['reused_existing'] or
                not (ROOT/r['artifact']).is_file() or (ROOT/r['artifact']).stat().st_size != r['bytes'] or r['bytes'] <= 0):
                raise ValueError('Incomplete or changed planned feature artifact')
            merged[i] = r
        summaries.append(s)
    if set(merged) != set(range(len(plan))):
        raise ValueError('Do not omit natural TRAIN windows')
    native_sentinels = {}
    for r in plan:
        w = r['window']
        if not r['reused_existing']: native_sentinels.setdefault((w['species'],w['seqid'],w['strand']),r)
    if len(native_sentinels) != 12:
        raise ValueError('Every TRAIN stratum must have new feature materialization')
    for r in native_sentinels.values():
        cached = torch.load(ROOT/r['artifact'],map_location='cpu',weights_only=True)
        Q.validate_cached(cached,r['window'])
        if cached['features'].dtype != torch.bfloat16 or not torch.isfinite(cached['features']).all():
            raise ValueError('New materialized BF16 TRAIN features invalid')
        del cached
    ordered = [merged[i] for i in sorted(merged)]
    with (FEATURES/'index.jsonl').open('x') as f:
        for r in ordered: f.write(json.dumps(r)+'\n')
    write_summary(FEATURES,{'experiment':'M28-FULL-TRAIN-FEATURES','completed':True,'feature_spec':FEATURE_SPEC,
                           'model':F.MODEL,'revision':F.REVISION,'train_unique_windows':len(ordered),'DEV_windows':0,
                           'native_windows':sum(s['native_windows'] for s in summaries),
                           'reused_existing_windows':sum(s['reused_existing_windows'] for s in summaries),
                           'referenced_cache_bytes':sum(r['bytes'] for r in ordered),
                           'native_materialization_checked_strata':12,'native_replay_source':str(REPLAY.relative_to(ROOT)),
                           'optimizer_steps_this_run':0,'reference_labels_used_in_encoding':False,
                           'test_setaria_used':False,'DEV_forward_used':False,'shards':summaries,
                           'scope':'all existing six-chromosome TRAIN features only; not new model training or accuracy'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('replay','shard','finalize'),required=True)
    parser.add_argument('--shard',type=int)
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'shard': shard(args.shard)
    else: {'replay':replay,'finalize':finalize}[args.mode]()
