#!/usr/bin/env python3
"""A17 label-only recovery of all closed free pools; never regenerates a chain."""
import argparse
from collections import Counter
import json
import shutil
import time
import torch
import recover_query_local_context as E

R = E.R
FAILED = R.READOUT
R.READOUT = R.ROOT/'outputs/M28-A17-LOCAL-CONTEXT-TRAIN-READOUT-LABEL-RECOVERY-R2'
E.V.OUT = R.ROOT/'outputs/M28-A17-LOCAL-CONTEXT-RESULT-CHECK-LABEL-RECOVERY-R2'


def labels():
    if ((FAILED/'JOBID').read_text().strip() != '13387946' or
        not (FAILED/'STATUS').read_text().strip().startswith('FAILED') or
        (FAILED/'summary.json').exists() or (FAILED/'windows.jsonl').stat().st_size):
        raise ValueError('Only the original failed-before-labels A17 readout may be recovered')
    if (R.FIT/'STATUS').read_text().strip() != 'COMPLETED' or (R.FIT/'JOBID').read_text().strip() != '13387910':
        raise ValueError('The same completed R1 fit is required')
    fitted = json.loads((R.FIT/'summary.json').read_text())
    R.validate_fit(fitted)
    raw = R.C.R.rows(FAILED/'free_pool.jsonl')
    draws, index, windows, catalogs = R.Q.load_train_inputs()
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    if len(draws) != 1536 or len(ids) != 1491 or [r['window']['id'] for r in raw] != ids:
        raise ValueError('All original ordered 1491 saved natural pools are required')
    original_free = {r['window']['id']: r for r in R.C.R.rows(R.C.FIT/'free_pool.jsonl')}
    prior = json.loads((R.D.OUT/'summary.json').read_text())
    if ((R.D.OUT/'STATUS').read_text().strip() != 'COMPLETED' or
        (R.D.OUT/'JOBID').read_text().strip() != '13360650' or
        not prior['completed'] or prior['source_generator'] != R.C.GENERATOR or prior['windows'] != 1491):
        raise ValueError('The original closed A12 baseline is required')
    old = {r['window']['id']: r for r in R.C.R.rows(R.D.OUT/'windows.jsonl')}
    for free in raw:
        ident = free['window']['id']
        w = windows[ident]
        if free['reference_used_in_pool_generation'] or w['split'] != 'train' or w['seqid'] not in R.C.R.ALLOW[w['species']]:
            raise ValueError('Original reference-free allowed TRAIN pools required')
        if any(free['window'][k] != w[k] or free['window'][k] != old[ident]['window'][k]
               for k in R.C.R.GEOMETRY):
            raise ValueError('Saved original geometry changed')
        if len({tuple(map(tuple, c)) for c in free['chains']}) != len(free['chains']):
            raise ValueError('Saved natural pool is not geometry-deduplicated')
    shutil.copyfile(FAILED/'free_pool.jsonl', R.READOUT/'free_pool.jsonl')
    if R.C.R.rows(R.READOUT/'free_pool.jsonl') != raw:
        raise ValueError('Recovered pool differs from the original saved pool')
    totals = {s: {a: Counter() for a in ('A3', 'A17')} for s in R.A.ALLOW}
    began = time.perf_counter()
    with (R.READOUT/'windows.jsonl').open('x') as log:
        for ordinal, free in enumerate(raw):
            ident = free['window']['id']
            w = windows[ident]
            cached = torch.load(R.ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            R.Q.validate_cached(cached, w)
            targets = R.build_targets(w, catalogs[w['species']], cached['sequence'])
            baseline = [r for r in old[ident]['candidates'] if r['origin'] == 'free']
            if ([r['chain'] for r in baseline] != original_free[ident]['chains'] or
                R.candidate_targets([r['chain'] for r in baseline], targets) != [r['label'] for r in baseline]):
                raise ValueError('Actual original natural candidates or labels changed')
            values = R.candidate_targets(free['chains'], targets)
            local = [dict(chain=c, label=y) for c, y in zip(free['chains'], values)]
            counts = dict(A3=R.counts(baseline), A17=R.counts(local))
            log.write(json.dumps(dict(window=free['window'], candidates=local, natural_counts=counts), allow_nan=False)+'\n')
            for arm in ('A3', 'A17'):
                totals[w['species']][arm].update(counts[arm])
            if (ordinal+1) % 128 == 0:
                log.flush()
                print(json.dumps(dict(labelled_saved_TRAIN_windows=ordinal+1, total=1491)), flush=True)
            del cached, targets
    actual = {s: {a: dict(c) for a, c in arms.items()} for s, arms in totals.items()}
    checks = R.gate(actual)
    R.write(R.READOUT, dict(completed=True, arm='A17', checkpoint=fitted['final_checkpoint'], windows=1491,
        natural_counts_per_species=actual, original_A3_actual_natural_pool_labels_recomputed=True,
        predeclared_TRAIN_generator_gate=checks, TRAIN_generator_supported=all(checks.values()),
        reference_used_in_pool_generation=False, reference_injected_candidates=0,
        new_candidate_search=False, optimizer_steps_this_run=0, new_backbone_forward=False,
        DEV_forward_used=False, test_setaria_used=False, seconds=time.perf_counter()-began,
        source_free_pool_job='13387946', source_failed_readout=str(FAILED.relative_to(R.ROOT)),
        original_saved_natural_pool_identity_verified=True, engineering_recovery='labels_only_JSON_nested_tuple_fix',
        limitation='Same original natural pools and last model; no new generation, fitting or independent gene F1.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('labels', 'check'), required=True)
    parser.add_argument('--readout-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'labels':
        labels()
    else:
        if not args.readout_job:
            raise ValueError('Actual recovered readout job ID required')
        E.V.main('13387910', args.readout_job)
