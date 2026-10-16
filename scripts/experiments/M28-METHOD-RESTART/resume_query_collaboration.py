#!/usr/bin/env python3
"""One fixed A2 timeout recovery: validate/copy saved prefix and decode missing DEV windows."""
import json
import math
import time
from collections import Counter
from pathlib import Path
import torch
import run_query_collaboration as Q

SOURCE = Q.STAGE_DIRS['infer']
OUT = Q.ROOT/'outputs/M28-A2-DEV-INFER-RECOVERY'
SOURCE_JOB = '13348644'


def validate_record(row, index, window_id):
    raw = row['raw_queries']
    completed = [r for r in raw if r['terminated']]
    if (row['window_index'] != index or row['window_id'] != window_id or
        row.get('reference_used') is not False or
        row['existence_logit_threshold'] != Q.CONFIG['existence_logit_threshold'] or
        row['score_kind'] != 'query_existence_logit_minus_threshold' or
        len(raw) != Q.CONFIG['num_queries'] or
        [r['query_index'] for r in raw] != list(range(Q.CONFIG['num_queries'])) or
        any(not math.isfinite(r['query_logit']) for r in raw) or
        row['chains'] != [r['chain'] for r in completed] or
        row['query_indices'] != [r['query_index'] for r in completed] or
        row['score'] != [r['query_logit']-Q.CONFIG['existence_logit_threshold'] for r in completed]):
        raise ValueError('Saved prefix identity/provenance/query outputs changed')
    expected = {'queries':len(raw),'completed_queries':len(completed),
                'unterminated_queries':len(raw)-len(completed),
                'completed_nonpositive_gain':sum(score <= 0 for score in row['score'])}
    if row['counts'] != expected:
        raise ValueError('Saved prefix counts disagree with raw queries')


def read_prefix(path, expected_ids):
    saved, incomplete = [], False
    if not path.exists(): return saved,incomplete
    with path.open() as handle:
        for line in handle:
            # The producer appends a newline per record. An interrupted final
            # write is never promoted to a completed window, even if JSON parses.
            if not line.endswith('\n'):
                incomplete = True
                break
            row = json.loads(line)
            i = len(saved)
            if i >= len(expected_ids): raise ValueError('Extra saved prefix window')
            validate_record(row,i,expected_ids[i]); saved.append(row)
    return saved,incomplete


def replay_equal(old, new):
    new = json.loads(json.dumps(new,allow_nan=False))
    a,b = old['raw_queries'],new['raw_queries']
    if len(a) != len(b): raise ValueError('Replay query count mismatch')
    error = max(abs(x['query_logit']-y['query_logit']) for x,y in zip(a,b))
    aa = [{k:v for k,v in r.items() if k != 'query_logit'} for r in a]
    bb = [{k:v for k,v in r.items() if k != 'query_logit'} for r in b]
    if (aa != bb or old['counts'] != new['counts'] or error > 1e-5 or
        [x > 0 for x in old['score']] != [x > 0 for x in new['score']]):
        raise ValueError('Same-checkpoint saved-prefix replay mismatch; do not reuse')
    return error


def count_record(counts, free):
    counts.update(free['counts'])
    counts['candidate_records'] += len(free['chains'])
    counts['all_queries_positive_completed'] += int(
        sum(score > 0 for score in free['score']) == Q.CONFIG['num_queries'])


def main():
    fitted = json.loads((Q.STAGE_DIRS['fit']/'summary.json').read_text())
    Q.check_last_identity(fitted)
    if (SOURCE/'JOBID').read_text().strip() != SOURCE_JOB:
        raise ValueError('Only the observed A2 timeout source is authorized for recovery')
    checkpoint = Q.ROOT/fitted['final_checkpoint']
    state = torch.load(checkpoint,map_location='cpu',weights_only=True)
    Q.check_last_identity(state,state=True)
    model = Q.new_model(); model.load_state_dict(state['model'],strict=True); del state
    model.eval().requires_grad_(False)
    groups = Q.dev_groups(Q.read_rows(Q.FEATURES/'index.jsonl'))
    began = time.perf_counter(); done = reused = fresh = 0
    scopes, replays, incomplete_files = [],[],[]
    torch.cuda.reset_peak_memory_stats()
    with torch.inference_mode():
        for number,((species,seqid,strand),rows) in enumerate(groups.items()):
            name = 'scope_'+str(number).zfill(2)+'.jsonl'
            saved,incomplete = read_prefix(SOURCE/name,[r['window_id'] for r in rows])
            if incomplete: incomplete_files.append(name)
            counts = Counter(); windows = []
            print(json.dumps({'scope':name,'reusable_prefix':len(saved),'missing':len(rows)-len(saved)}),flush=True)
            with (OUT/name).open('x') as handle:
                for i,row in enumerate(rows):
                    a = int(row['window_id'].split('|')[2])
                    b = min(a+24576,Q.DEV_LENGTHS[species][seqid]); windows.append([a,b])
                    # One sentinel per existing oriented scope detects wrong
                    # weights/decoder before its stored rows are reused.
                    if i >= len(saved) or i == 0:
                        cached = torch.load(Q.ROOT/row['artifact'],map_location='cpu',weights_only=True)
                        w = {'id':row['window_id'],'species':species,'split':'val','seqid':seqid,
                             'strand':strand,'start':a,'end':b,'width':24576,'valid_bases':b-a}
                        features,dna = Q.window_inputs(cached,w,'cuda')
                        outputs = model(features,dna,[b-a])
                        generated = Q.query_window_candidates(model.query,outputs['features'][0],
                            cached['sequence'],i,w['id'],Q.CONFIG['max_exons'],
                            Q.CONFIG['existence_logit_threshold'],b-a)
                        if i < len(saved):
                            error = replay_equal(saved[i],generated)
                            replays.append({'window_id':row['window_id'],'discrete_equal':True,
                                            'max_abs_query_logit_difference':error})
                        del cached,features,dna,outputs
                    if i < len(saved):
                        free = saved[i]; reused += 1
                    else:
                        free = generated; fresh += 1
                    handle.write(json.dumps(free,allow_nan=False)+'\n')
                    count_record(counts,free); done += 1
                    if done % 32 == 0:
                        handle.flush()
                        print(json.dumps({'windows':done,'reused':reused,'new':fresh,
                                          'total':8690,'recovery_seconds':time.perf_counter()-began}),flush=True)
                del saved
            scopes.append({'species':species,'seqid':seqid,'strand':strand,'windows':windows,
                           'file':name,'candidate_records':counts['candidate_records'],
                           'query_counts':dict(counts)})
    if done != 8690: raise ValueError('Incomplete recovered A2 full DEV grid')
    report = {'arm':'A2','completed':True,'inference_complete':True,'windows':done,
              'optimizer_steps':Q.STEPS,'configuration':Q.CONFIG,'feature_revision':Q.F.REVISION,
              'proposal':Q.PROPOSAL,'checkpoint':fitted['final_checkpoint'],
              'reference_used':False,'test_setaria_used':False,'scopes':scopes,
              'wall_seconds_recovery_only':time.perf_counter()-began,
              'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
              'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
              'recovery':{'source':str(SOURCE.relative_to(Q.ROOT)),'source_job':SOURCE_JOB,
                          'source_root_state':'TIMEOUT','source_elapsed_seconds':7227,
                          'reused_windows':reused,'new_windows':fresh,'sentinel_replays':replays,
                          'discarded_incomplete_final_records':incomplete_files,
                          'source_files_unchanged':True,'model_or_operating_point_changed':False}}
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(report['recovery']),flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
