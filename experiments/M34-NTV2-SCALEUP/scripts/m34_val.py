"""M34 internal-validation scoring (pooled exact-chain F1, CDS >= 96 bp rule, M26 CDS-assessable primary references wholly inside the region) of an
m34_infer run, for legacy held-out chromosomes and new-species validation pseudo-contigs.  Wraps the M31 `c1_val.score_region` (unchanged).
One compact line per call is appended to --metrics; full detail incl. strata goes to <run>/valscore_<tag>.json.
Usage: m34_val.py <run_dir> <step> <regions> --tag select|own_species|transfer [--arm ARM] [--metrics FILE]
"""
import argparse, json, time

import m34_train_common as T
T.register()
import c1_val as V                      # M31 (unchanged)
import c1_infer as I
M = T.M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('step', type=int)
    ap.add_argument('regions')
    ap.add_argument('--tag', required=True)
    ap.add_argument('--arm', default=None)
    ap.add_argument('--metrics', default=None)
    args = ap.parse_args()
    t0 = time.time()
    rows = [V.score_region(args.run, *I.parse_region(r)) for r in args.regions.split(',')]
    tp, p, r = (sum(x[k] for x in rows) for k in ('tp', 'pred', 'ref'))
    out = dict(arm=args.arm, step=args.step, tag=args.tag, pooled_F1=V.f1(tp, p, r), pooled_tp=tp, pooled_pred=p, pooled_ref=r,
               pooled_precision=100 * tp / max(p, 1), pooled_recall=100 * tp / max(r, 1),
               per_region={f"{x['seq']}:{x['start']}-{x['end']}": dict(F1=round(x['F1'], 3), tp=x['tp'], pred=x['pred'], ref=x['ref']) for x in rows},
               time=time.strftime('%Y-%m-%dT%H:%M:%S'), score_seconds=round(time.time() - t0, 1))
    (M.Path(args.run) / f'valscore_{args.tag}.json').write_text(json.dumps(dict(out, strata=[x['strata'] for x in rows]), indent=1))
    line = json.dumps(out)
    if args.metrics:
        with open(args.metrics, 'a') as h:
            h.write(line + '\n')
    print(line)


if __name__ == '__main__':
    main()
