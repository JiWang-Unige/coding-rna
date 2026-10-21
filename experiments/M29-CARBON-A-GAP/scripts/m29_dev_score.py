"""Score DEV predictions (GFF3) with the frozen ruler: python m29_dev_score.py name=path_ara.gff[:fmt] name=path_rice.gff ...
Each argument  label:ara_gff:rice_gff[:carbon]  ('carbon' = also report conf>=0.1|NA operating point).  Writes outputs/M29-CARBON-A-GAP/G0/dev_scores.json"""
import json, sys
from pathlib import Path
import m29_common as M
from m29_r3_posterior_filter import carbon_conf

res = {}
for arg in sys.argv[1:]:
    label, ara, rice, *rest = arg.split(':')
    carbon = bool(rest)
    rows = {}
    for (sp, seq, L), path in zip(M.DEV_SCOPE, (ara, rice)):
        ref = M.build_reference(seq)
        rk = {M.key(t) for t in ref['refs']}
        txs = M.load_txs(path, seq)
        conf = carbon_conf(path) if carbon else None
        pts = {'all': {}}
        if carbon:
            pts['conf>=0.1|NA'] = {}
        for tx in txs:
            k = M.key(tx)
            pts['all'].setdefault(k, 1)
            if carbon:
                c = conf[tx['gene_id']]
                if c is None or c >= 0.1:
                    pts['conf>=0.1|NA'].setdefault(k, 1)
        rows[seq] = {n: dict(tp=sum(1 for k in d if k in rk), pred=len(d), ref=len(rk)) for n, d in pts.items()}
    out = {}
    for n in rows[M.DEV_SCOPE[0][1]]:
        tp = sum(rows[s][n]['tp'] for s in rows); pr = sum(rows[s][n]['pred'] for s in rows); rf = sum(rows[s][n]['ref'] for s in rows)
        out[n] = dict(pooled=dict(tp=tp, pred=pr, ref=rf, F1=100 * 2 * tp / (pr + rf)),
                      per_seq={s: dict(rows[s][n], F1=100 * 2 * rows[s][n]['tp'] / (rows[s][n]['pred'] + rows[s][n]['ref'])) for s in rows})
        print(label, n, out[n]['pooled'], {s: round(v['F1'], 3) for s, v in out[n]['per_seq'].items()}, flush=True)
    res[label] = out
(M.M29 / 'G0' / 'dev_scores.json').write_text(json.dumps(res, indent=1))
