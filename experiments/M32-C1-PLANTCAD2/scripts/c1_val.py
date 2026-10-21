"""Internal-validation score for a c1_infer run on region(s) of the held-out TRAIN chromosomes (never DEV).
Reference = M26 CDS-assessable primary chains whose every CDS exon lies inside the region; prediction = unique exact chains (CDS >= 96 bp,
the rule used by the gate/Carbon-A scoring and by M29-A1/A2) from the region GFF.  One compact jsonl line is appended per evaluation.
Usage: c1_val.py <run_dir> <step> <regions> [--metrics FILE]
"""
import argparse, json, sys, time
from collections import defaultdict

import m32_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m29_common as MC            # noqa: E402
import c1_infer as I               # noqa: E402


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'


def f1(tp, p, r):
    return 100 * 2 * tp / (p + r) if p + r else 0.


def score_region(run, seq, a, b):
    L = MC.LENGTH[seq]
    b, name = I.resolve(seq, a, b, L)
    ref = MC.build_reference(seq)
    keys, nex = set(), {}
    for t in ref['refs']:
        if all(x >= a and y <= b for x, y, _ in t['CDS']):
            k = MC.key(t); keys.add(k); nex[k] = len(t['CDS'])
    pred = {}
    for tx in MC.load_txs(M.Path(run) / name / 'prediction' / 'full_grammar.gff3', seq):
        ex = [(x, y) for x, y, _ in tx['CDS']]
        if sum(y - x for x, y in ex) < 96:
            continue
        pred.setdefault(MC.key(tx), len(ex))
    tp = sum(1 for k in pred if k in keys)
    st = defaultdict(lambda: [0, 0, 0])
    for k, n in pred.items():
        st[strat(n)][1] += 1
        st[strat(n)][0] += k in keys
    for k in keys:
        st[strat(nex[k])][2] += 1
    return dict(seq=seq, start=a, end=b, tp=tp, pred=len(pred), ref=len(keys), F1=f1(tp, len(pred), len(keys)),
                strata={s: dict(tp=v[0], pred=v[1], ref=v[2], F1=f1(*v)) for s, v in st.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('step', type=int)
    ap.add_argument('regions')
    ap.add_argument('--metrics', default=None)
    ap.add_argument('--tag', default='internal_validation_region')
    args = ap.parse_args()
    t0 = time.time()
    rows = [score_region(args.run, *I.parse_region(r)) for r in args.regions.split(',')]
    tp, p, r = (sum(x[k] for x in rows) for k in ('tp', 'pred', 'ref'))
    out = dict(step=args.step, tag=args.tag, pooled_F1=f1(tp, p, r), pooled_tp=tp, pooled_pred=p, pooled_ref=r,
               pooled_precision=100 * tp / max(p, 1), pooled_recall=100 * tp / max(r, 1),
               per_region={f"{x['seq']}:{x['start']}-{x['end']}": dict(F1=round(x['F1'], 3), tp=x['tp'], pred=x['pred'], ref=x['ref']) for x in rows},
               time=time.strftime('%Y-%m-%dT%H:%M:%S'), score_seconds=round(time.time() - t0, 1))
    (M.Path(args.run) / 'valscore.json').write_text(json.dumps(dict(out, strata=[x['strata'] for x in rows]), indent=1))
    line = json.dumps(out)
    if args.metrics:
        with open(args.metrics, 'a') as h:
            h.write(line + '\n')
    print(line)


if __name__ == '__main__':
    main()
