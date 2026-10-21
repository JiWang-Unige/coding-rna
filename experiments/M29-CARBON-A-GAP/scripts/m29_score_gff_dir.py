"""Score a directory tree  <root>/<variant>/<seq>/prediction.gff  with the frozen ruler (both Carbon-A operating points)."""
import json
import sys
from pathlib import Path
import m29_common as M
from m29_r3_posterior_filter import carbon_conf


def score_variant(root, variant, seqs=None, fname='prediction.gff'):
    seqs = seqs or [s for _, s, _ in M.SCOPE]
    tot = {'unfiltered': [0, 0, 0], 'conf>=0.1|NA': [0, 0, 0]}
    per = {}
    for sp, seq, L in M.SCOPE:
        if seq not in seqs:
            continue
        gff = Path(root) / variant / seq / fname
        ref = M.build_reference(seq)
        rk = {M.key(t) for t in ref['refs']}
        conf = carbon_conf(gff)
        txs = M.load_txs(gff, seq)
        uniq, filt = {}, {}
        for tx in txs:
            k = M.key(tx)
            uniq.setdefault(k, None)
            c = conf[tx['gene_id']]
            if c is None or c >= 0.1:
                filt.setdefault(k, None)
        row = {}
        for name, d in (('unfiltered', uniq), ('conf>=0.1|NA', filt)):
            tp = sum(1 for k in d if k in rk)
            row[name] = dict(tp=tp, pred=len(d), ref=len(rk), F1=100 * 2 * tp / (len(d) + len(rk)))
            tot[name][0] += tp; tot[name][1] += len(d); tot[name][2] += len(rk)
        per[seq] = row
    pooled = {k: dict(tp=v[0], pred=v[1], ref=v[2], F1=100 * 2 * v[0] / (v[1] + v[2])) for k, v in tot.items()}
    return dict(variant=variant, pooled=pooled, per_contig=per)


if __name__ == '__main__':
    root = Path(sys.argv[1]); out = []
    for v in sys.argv[2].split(','):
        r = score_variant(root, v, sys.argv[3].split(',') if len(sys.argv) > 3 else None)
        out.append(r)
        print(v, {k: (d['tp'], d['pred'], round(d['F1'], 3)) for k, d in r['pooled'].items()}, flush=True)
    Path(root / f'scores_{sys.argv[2].replace(",", "+")}.json').write_text(json.dumps(out, indent=1))
