#!/usr/bin/env python3
"""Real TRAIN-only raw-FASTA encoding/head/decode/export check, never model accuracy."""
import json
import math
import time
from pathlib import Path
import torch
import infer_full_train as I
from cache_full_train import PLAN, read_sequences
from evaluate_query import E, validate_gff_roundtrip

OUT = I.ROOT/'outputs/M28-FULL-TRAIN-NATIVE-INFERENCE-SMOKE'


@torch.inference_mode()
def main():
    planned = json.loads((PLAN/'summary.json').read_text())
    selected = {}
    with (PLAN/'replay_sentinels.jsonl').open() as f:
        for line in f:
            row = json.loads(line)
            selected.setdefault(row['window']['species'],row)
    if set(selected) != set(I.R.ALLOW):
        raise ValueError('One pre-existing replay sentinel per TRAIN species required')
    for species,row in selected.items():
        w = row['window']
        if (w['split'] != 'train' or w['seqid'] not in I.R.ALLOW[species] or
            w['end']-w['start'] != w['valid_bases'] or w['valid_bases'] != I.M.WIDTH):
            raise ValueError('Full-width fixed TRAIN sentinel geometry required')
    sequences = read_sequences(planned,list(selected.values()))
    samples = []
    for species,row in selected.items():
        w = row['window']
        sequence = sequences[species][w['seqid']][w['start']:w['end']]
        samples.append({'species':species,'source_window':w,'sequence':sequence,
                        'seqid':'TRAIN_SENTINEL_'+w['seqid']+'_'+str(w['start'])})
    del sequences
    backbone,tokenizer = I.load_native()
    if any(p.requires_grad for p in backbone.parameters()):
        raise ValueError('Native encoder must stay frozen')
    began = time.perf_counter(); results = []
    for arm in ('C0','Q'):
        out = OUT/arm; out.mkdir(exist_ok=False)
        model = I.R.new_model(arm,'cuda').eval().requires_grad_(False)
        torch.cuda.reset_peak_memory_stats(); started = time.perf_counter()
        transcripts,scopes,lengths = [],[],{}
        for ordinal,sample in enumerate(samples):
            txs,ss = I.annotate_contig(arm,model,backbone,tokenizer,sample['species'],sample['seqid'],
                                       sample['sequence'],ordinal,out,'cuda')
            transcripts.extend(txs); scopes.extend(ss); lengths[sample['seqid']] = len(sample['sequence'])
        records = []
        for scope in scopes:
            rows = list(I.saved_records(out/scope['file']))
            if (len(rows) != 1 or scope['windows'] != [(0,I.M.WIDTH)] or
                rows[0]['window_index'] != 0 or rows[0]['reference_used'] is not False or
                len(rows[0]['chains']) != len(rows[0]['score']) or
                any(not math.isfinite(x) for x in rows[0]['score'])):
                raise ValueError('Native free-record completion, geometry or provenance mismatch')
            if arm=='Q' and (rows[0]['counts']['queries'] != 16 or len(rows[0]['raw_queries']) != 16 or
                             rows[0]['existence_logit_threshold'] != 0.):
                raise ValueError('Native Q count/operating point changed')
            records.append({'species':scope['species'],'strand':scope['strand'],
                            'free_candidates':len(rows[0]['chains']),'assembly_counts':scope['assembly_counts']})
        gff = out/(arm+'.gff3'); I.write_gff3(gff,transcripts,lengths,'M28_NATIVE_INFERENCE_ENGINEERING')
        parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(gff,lengths),complete_only=False)
        validate_gff_roundtrip(parsed,transcripts)
        torch.cuda.synchronize()
        results.append({'arm':arm,'windows':sum(len(s['windows']) for s in scopes),
                        'scopes':scopes,'records':records,'GFF3_roundtrip_exact':True,
                        'primary_CDS_predictions':len(transcripts),
                        'head_initialization':'fixed seed0, untrained; no final-checkpoint impersonation',
                        'optimizer_steps':0,'trainable_parameters':sum(p.numel() for p in model.parameters() if p.requires_grad),
                        'seconds':time.perf_counter()-started,
                        'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
                        'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved()})
        del model; torch.cuda.empty_cache()
    report = {'completed':True,'experiment':'M28-FULL-TRAIN-NATIVE-INFERENCE-SMOKE','engineering_only':True,
              'feature_spec':I.R.FEATURE_SPEC,'configurations':I.R.CONFIGS,
              'source_windows':[s['source_window'] for s in samples],
              'coordinate_system':'each cropped TRAIN sentinel starts at zero; not whole-chromosome annotation',
              'arms':results,'total_direction_windows':sum(a['windows'] for a in results),
              'raw_DNA_encoding_used':True,'reference_or_labels_used':False,'feature_cache_used':False,
              'checkpoints_read_or_written':False,'optimizer_steps':0,'DEV_forward_used':False,
              'test_setaria_used':False,'backbone_trainable_parameters':0,
              'seconds':time.perf_counter()-began,
              'limitations':'tests native inference integration only; trained last-checkpoint main and full DEV accuracy pending'}
    with (OUT/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
