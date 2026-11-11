"""M41 CPU dry run (no training, no GPU): decode STORED potentials with the frozen skyline and with the M41 grammar at its initial prior, and compare
the decodes (exact-chain F1 vs the reference inside the slice, by exon-count stratum; predicted intron-length quantiles).

Potentials: outputs/M36-BOOTSTRAP-ADAPT/zeroshot_A/pot/<contig>_+.npy = T10 step 3072 (pc2) on the Dryas octopetala M30 set-A contigs (M36 adaptation
contigs: not DEV, not evaluation set B).  The DEV potentials of T2 were never kept (inference deletes them), so a DEV dry run needs the GPU encoder and is
part of the GPU smoke (protocol 33 section 5).  Only '+' strand rows and '+' strand references are used (oriented = genomic coordinates).
This is an engineering check of the decoder on real potentials (runs, time, memory, how far the untrained prior moves the decode), NOT an accuracy result.

Usage (ont):  nice -n 19 taskset -c 28-31 ~/envs/generanno/bin/python m41_dryrun.py --pot .../NC_141366.1_+.npy --contig NC_141366.1 --start 0 --end 4000000
              --prior ../results/intron_prior_fit.json --out ../results/dryrun_NC_141366.1_0-4M.json
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts'))
import m30_common as M30           # noqa: E402  (a1_common -> pipeline/m28 + M29 helpers on sys.path)
sys.path.insert(0, str(HERE))
import m41_grammar as GR           # noqa: E402
import m41_prior as P              # noqa: E402


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'


def f1(tp, p, r):
    return 100 * 2 * tp / (p + r) if p + r else 0.


def score(chains, refs, a):
    pred = {}
    for c in chains:
        ex = tuple((x + a, y + a) for x, y in c)
        if sum(y - x for x, y in ex) >= 96:
            pred.setdefault(ex, len(ex))
    st = {s: [0, 0, 0] for s in ('1', '2-5', '>=6')}
    for k, n in pred.items():
        st[strat(n)][1] += 1
        st[strat(n)][0] += k in refs
    for k, n in refs.items():
        st[strat(n)][2] += 1
    tp, p, r = (sum(v[i] for v in st.values()) for i in range(3))
    return dict(F1=round(f1(tp, p, r), 3), tp=tp, pred=p, ref=r, strata={s: dict(F1=round(f1(*v), 3), tp=v[0], pred=v[1], ref=v[2]) for s, v in st.items()})


def intron_q(chains):
    L = np.array([c - b for ch in chains for (_, b), (c, _) in zip(ch, ch[1:])])
    if not len(L):
        return {}
    return dict(n=int(len(L)), q=[int(x) for x in np.percentile(L, [1, 5, 25, 50, 75, 95, 99])], over_5kb=int((L > 5000).sum()), max=int(L.max()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pot', required=True)
    ap.add_argument('--contig', required=True)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--end', type=int, required=True)
    ap.add_argument('--prior', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    torch.set_num_threads(1)
    MC = M30.register()
    a, b = args.start, args.end
    seq = re.sub('[^ACGT]', 'N', MC.read_sequence(args.contig).upper())[a:b]
    pot = np.load(args.pot, mmap_mode='r')
    lg = torch.from_numpy(np.ascontiguousarray(pot[a:b], dtype=np.float32))
    assert lg.shape == (len(seq), 15)
    ref = MC.build_reference(args.contig)
    refs = {}
    for t in ref['refs']:
        k = MC.key(t)
        if k[1] == '+' and all(x >= a and y <= b for x, y, _ in t['CDS']):
            refs[tuple(tuple(e) for e in k[2])] = len(t['CDS'])
    fits = {f['short_stages']: f for f in json.load(open(args.prior))['fits']}
    rows = {}
    from m29_skyline_reuse import decode_full_grammar_skyline
    t0 = time.perf_counter()
    r0 = decode_full_grammar_skyline(seq, lg, 128)
    rows['frozen'] = dict(seconds=round(time.perf_counter() - t0, 1), **score(r0['chains'], refs, a), introns=intron_q(r0['chains']), chains=len(r0['chains']))
    variants = [('m41_off', GR.OFF, None)]
    for r, gauge in ((1, 'shape'), (1, 'absolute'), (4, 'shape')):
        f = fits.get(r)
        if f:
            variants.append((f'r{r}_{gauge}', GR.IntronLength(True, r), P.bonus(torch.tensor(f['theta'], dtype=torch.float64), r, f['p0'], gauge)))
    for name, il, bon in variants:
        t0 = time.perf_counter()
        res = GR.decode_skyline(seq, lg, 128, bon, il)
        row = dict(seconds=round(time.perf_counter() - t0, 1), **score(res['chains'], refs, a), introns=intron_q(res['chains']), chains=len(res['chains']),
                   state_count=res['counts']['state_count'], checkpoint_MB=round(res['counts']['checkpoint_bytes'] / 2**20, 1))
        if il.enabled:
            row['branches'] = dict(S=res['intron_branches'].count('S'), L=res['intron_branches'].count('L'))
            same = len(set(res['chains']) & set(r0['chains']))
            row['chains_identical_to_frozen'] = same
        else:
            row['identical_to_frozen'] = res['chains'] == r0['chains'] and res['gain'] == r0['gain']
        rows[name] = row
        print(json.dumps({name: row}), flush=True)
    out = dict(pot=args.pot, model='M34 pc2_T10 step 3072 (M36 zeroshot_A potentials)', contig=args.contig, strand='+', start=a, end=b,
               reference_chains_plus_strand=len(refs), reference_intron_q=intron_q(list(refs)), rows=rows,
               note='engineering dry run of the decoder on stored potentials of a head trained under the FROZEN grammar; not an accuracy result')
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(dict(frozen=rows['frozen'])))


if __name__ == '__main__':
    main()
