"""Paired comparison on repeated windows: NLL of window w at training step s vs the same window in the offset-only baseline (pilot_lr0, 100 draws-batches)."""
import json, os, sys, numpy as np
B = os.path.expanduser('~/coding-rna/outputs/M29-A1-RESIDUAL')
base = {}
for l in open(f'{B}/pilot_lr0/training.jsonl'):
    r = json.loads(l)
    for w in r['windows']:
        if w.get('known', 0) > 0:
            base.setdefault(w['window'], w['raw_NLL'])
rows = {}
for l in open(f'{B}/{sys.argv[1]}/training.jsonl'):
    r = json.loads(l); rows[r['step']] = r
lo_hi = [(101, 500), (501, 1000), (1001, 1500), (1501, 2000), (2001, 2500), (2501, 3000), (3001, 3500), (3501, 4000), (4001, 4608)]
for lo, hi in lo_hi:
    p = [(base[w['window']], w['raw_NLL'], w['known']) for s in range(lo, hi + 1) if s in rows for w in rows[s]['windows'] if w.get('known', 0) > 0 and w['window'] in base]
    if not p: continue
    b = np.array([x[0] for x in p]); a = np.array([x[1] for x in p]); k = np.array([x[2] for x in p])
    print(f'steps {lo}-{hi}: pairs {len(p)}  mean normNLL base {np.mean(b/k):.5f} -> now {np.mean(a/k):.5f}  sum NLL {b.sum():.0f}->{a.sum():.0f}  #NLL>0.5 {int((b>0.5).sum())}->{int((a>0.5).sum())}  better {int((a<b-0.5).sum())} worse {int((a>b+0.5).sum())}')
