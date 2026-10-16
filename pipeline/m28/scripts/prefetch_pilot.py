#!/usr/bin/env python3
"""R4 read-only cache read-ahead, inside an existing Slurm allocation.

No tensors/models/labels are loaded and no cache or scientific output is changed.
The helper ends on parent STATUS completion, the fixed end count, or its timeout.
"""
import argparse
import json
import os
import time
from pathlib import Path

ROOT = Path('/home/users/j/jwang/coding-rna')
PILOT = ROOT / 'outputs/M28-PILOT-R4'

def rows(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]

def ordered_paths(stage):
    index = {r['window_id']: r for r in rows(PILOT / 'features/index.jsonl')}
    if stage == 'fit':
        draws = rows(ROOT / 'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl')
        assert len(draws) == 1536
        selected = [index[d['window_id']] for d in draws] * 3
        assert all(r['split'] == 'train' for r in selected)
    else:
        selected = [r for r in index.values() if r['split'] == 'val']
        def key(r):
            species, seqid, start, strand = r['window_id'].split('|')
            return species, seqid, strand, int(start)
        selected.sort(key=key)
        assert len(selected) == 8690
        assert {(r['species'], r['window_id'].split('|')[1]) for r in selected} == {
            ('arabidopsis_thaliana', 'NC_003074.8'), ('oryza_sativa', 'NC_089041.1')}
    return [(r['window_id'], ROOT / r['artifact']) for r in selected]

def progress(log, arm, field):
    with log.open('rb') as handle:
        handle.seek(max(0, log.stat().st_size - 65536))
        lines = handle.read().splitlines()
    for line in reversed(lines):
        try:
            record = json.loads(line)
        except (ValueError, UnicodeDecodeError):
            continue
        if isinstance(record, dict) and record.get('arm') == arm and field in record:
            return int(record[field])
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=('fit', 'infer'), required=True)
    parser.add_argument('--arm', choices=('C0', 'B1'), required=True)
    parser.add_argument('--job-id', required=True)
    parser.add_argument('--seconds', type=int, default=3300)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    paths = ordered_paths(args.stage)
    out = PILOT / (args.arm if args.stage == 'fit' else 'infer_' + args.arm)
    log = PILOT / 'logs' / ('M28' + args.stage.upper() + '-' + args.arm + '_' + args.job_id + '.out')
    field = 'step' if args.stage == 'fit' else 'windows'
    done = progress(log, args.arm, field)
    print(json.dumps({'stage': args.stage, 'arm': args.arm, 'parent_job': args.job_id,
                      'completed_at_start': done, 'total': len(paths), 'lookahead': 128,
                      'next_window': paths[done][0] if done < len(paths) else None}), flush=True)
    if args.inspect:
        return
    if os.environ.get('SLURM_JOB_ID') != args.job_id:
        raise RuntimeError('Read-ahead must run in the declared existing allocation')
    began = time.monotonic()
    fetched = done
    count = 0
    buffer = bytearray(1024 * 1024)
    while time.monotonic() - began < args.seconds:
        if (out / 'STATUS').read_text().strip() != 'RUNNING':
            break
        done = progress(log, args.arm, field)
        if done >= len(paths):
            break
        fetched = max(fetched, done)
        end = min(done + 128, len(paths))
        while fetched < end and time.monotonic() - began < args.seconds:
            with paths[fetched][1].open('rb', buffering=0) as handle:
                while handle.readinto(buffer):
                    pass
            fetched += 1
            count += 1
            if count % 128 == 0:
                print(json.dumps({'read_files': count, 'prefetched_through': fetched,
                                  'seconds': time.monotonic() - began}), flush=True)
        time.sleep(1)
    print(json.dumps({'finished': True, 'read_files': count,
                      'seconds': time.monotonic() - began}), flush=True)

if __name__ == '__main__':
    main()
