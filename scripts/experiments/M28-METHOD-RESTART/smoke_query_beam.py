#!/usr/bin/env python3
"""Eight fixed TRAIN engineering windows; no fitting, reference accuracy or DEV."""
import argparse
import json
import math
import time
from collections import Counter
from pathlib import Path
import torch
import run_full_train as R
from src.m28.candidates import STOP
from src.m28.query_decode import decode_query_chains
from src.m28.query_beam import decode_query_beam

OUT = R.ROOT/'outputs/M28-A5-QUERY-BEAM-ENGINEERING'
LIMITS = {'beam_width':4,'completed_limit_per_query':4,'max_exons':128}


def check_chain(chain, sequence):
    previous = -1
    for j,(a,b) in enumerate(chain):
        if not 0 <= a < b <= len(sequence) or (j and a-previous < 5):
            raise ValueError('Beam output has invalid exon/intron geometry')
        if j and (sequence[previous:previous+2] not in ('GT','GC') or sequence[a-2:a] != 'AG'):
            raise ValueError('Beam output has noncanonical splice')
        previous = b
    cds = ''.join(sequence[a:b] for a,b in chain)
    if (not chain or len(cds)%3 or not cds.startswith('ATG') or cds[-3:] not in STOP or
        any(c not in 'ACGT' for c in cds) or
        any(cds[i:i+3] in STOP for i in range(0,len(cds)-3,3))):
        raise ValueError('Beam output violates complete CDS/codon-carry grammar')


@torch.inference_mode()
def main(out):
    R.require_readiness()
    feature_summary = json.loads((R.FEATURES/'summary.json').read_text())
    R.validate_feature_summary(feature_summary)
    witnesses = R.rows(R.READY/'witnesses.jsonl')
    plan = {r['window']['id']:r for r in R.rows(R.PLAN/'feature_plan.jsonl')}
    if (len(witnesses) != 8 or len({r['window_id'] for r in witnesses}) != 8 or
        Counter(r['species'] for r in witnesses) != Counter({s:4 for s in R.ALLOW})):
        raise ValueError('Exactly the pre-existing eight TRAIN sentinels required')
    model = R.new_model('Q','cuda').eval().requires_grad_(False)
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter(); records = []
    with (out/'records.jsonl').open('x') as log:
        for witness in witnesses:
            entry = plan[witness['window_id']]; w = entry['window']
            if (w['split'] != 'train' or w['seqid'] not in R.ALLOW[w['species']] or
                w['valid_bases'] != 24576):
                raise ValueError('Full-width allowed TRAIN window required')
            cached = torch.load(R.ROOT/entry['artifact'],map_location='cpu',weights_only=True)
            features,dna = R.Q.window_inputs(cached,w,'cuda')
            outputs = model(features,dna,[w['valid_bases']])
            hidden = outputs['features'][0]
            sequence = cached['sequence'].upper()
            torch.cuda.synchronize(); started = time.perf_counter()
            greedy = decode_query_chains(model.query,hidden,sequence,128)
            torch.cuda.synchronize(); greedy_seconds = time.perf_counter()-started
            started = time.perf_counter()
            branched = decode_query_beam(model.query,hidden,sequence,128,4,4)
            torch.cuda.synchronize(); beam_seconds = time.perf_counter()-started
            if len(greedy) != 16 or len(branched) != 16:
                raise ValueError('Q16 identity changed')
            gc = {tuple(r['chain']) for r in greedy if r['terminated']}
            bc = {tuple(r['chain']) for q in branched for r in q['completed']}
            for chain in gc|bc:
                check_chain(chain,sequence)
            for q in branched:
                if len(q['completed']) > 4 or not math.isfinite(q['query_logit']):
                    raise ValueError('Completion budget or score mismatch')
                if any(not math.isfinite(r['path_log_score']) or r['path_log_score'] > 1e-6
                       for r in q['completed']):
                    raise ValueError('Path score must be finite and nonpositive')
            record = {'window':w,'category':witness['category'],'greedy_completed_records':sum(
                         r['terminated'] for r in greedy),
                      'beam_completed_records':sum(len(q['completed']) for q in branched),
                      'greedy_unique_chains':len(gc),'beam_unique_chains':len(bc),
                      'unique_both':len(gc&bc),'unique_only_beam':len(bc-gc),
                      'unique_only_greedy':len(gc-bc),'greedy_seconds':greedy_seconds,
                      'beam_seconds':beam_seconds,'beam_queries':branched,
                      'greedy_queries':greedy,'all_returned_chains_canonical':True}
            log.write(json.dumps(record,allow_nan=False)+'\n'); log.flush()
            records.append({k:v for k,v in record.items() if k not in ('greedy_queries','beam_queries')})
            print(json.dumps(records[-1],allow_nan=False),flush=True)
            del cached,features,dna,outputs,hidden,greedy,branched
    report = {'experiment':'M28-A5-QUERY-BEAM-ENGINEERING','completed':True,
              'engineering_only':True,'limits':LIMITS,'model_configuration':R.CONFIGS['Q'],
              'head_initialization':'seed0, untrained, eval, requires_grad false',
              'trainable_head_parameters':sum(p.numel() for p in model.parameters() if p.requires_grad),
              'optimizer_steps':0,'checkpoints_read_or_written':False,
              'feature_spec':R.FEATURE_SPEC,'feature_cache_used':True,
              'raw_DNA_encoding_this_run':False,'new_backbone_forward':False,
              'reference_targets_used_in_decoding':False,
              'sentinel_selection':'pre-existing TRAIN readiness metadata; not a random accuracy cohort',
              'DEV_forward_used':False,'test_setaria_used':False,
              'windows':len(records),'records':records,
              'seconds':time.perf_counter()-began,
              'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
              'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
              'score_kind':'cumulative geometric boundary and EOS log score, not gene correctness probability',
              'limitations':'untrained candidate/grammar/runtime check only; no accuracy, independent generalization, new fit or gate'}
    R.write_summary(out,report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.output_dir)
