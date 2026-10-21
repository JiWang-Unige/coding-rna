"""Summarise lr pilots / training logs: paired per-window comparison against the lr=0 (offset-only) baseline, same draws."""
import json, sys, glob, os
import numpy as np
base = os.path.expanduser('~/coding-rna/outputs/M29-A1-RESIDUAL')

def load(d):
    rows = {}
    for l in open(f'{base}/{d}/training.jsonl'):
        r = json.loads(l); rows[r['step']] = r          # last record per step wins (resume duplicates)
    return rows

def stats(rows, lo, hi):
    nll, fail = [], 0
    n = 0
    for s in range(lo, hi + 1):
        if s not in rows: continue
        for w in rows[s]['windows']:
            v = w.get('raw_NLL', 0.); k = max(1, w.get('known', 1))
            nll.append(v / k if not w.get('skipped_all_unknown') else 0.)
            fail += v > 0.5; n += 1
    nll = np.array(nll)
    return dict(windows=n, mean_normNLL=float(nll.mean()), median=float(np.median(nll)), windows_with_NLL_gt_0p5=int(fail))

if __name__ == '__main__':
    names = sys.argv[1:]
    b = load(names[0])
    for d in names:
        r = load(d)
        steps = sorted(r)
        print(d, 'steps', len(steps), 'grad_norm med', np.median([r[s]['grad_norm'] for s in steps]),
              'all', stats(r, 1, 100), 'last50', stats(r, 51, 100), 'baseline_last50', stats(b, 51, 100) if d != names[0] else '')
