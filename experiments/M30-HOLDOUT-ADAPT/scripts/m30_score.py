"""M30 scoring on the evaluation chromosomes (frozen ruler: M26 CDS-assessable primaries, exact CDS chain), CDS-base F1, background burden,
exon-count strata, and gene-unit paired bootstrap (protocol §6)."""
import json, sys
from collections import defaultdict
import numpy as np
import m30_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
M29 = M.register()
import m29_r4_error_types as R4          # noqa: E402  (Index, ov_bp)
from m29_r3_posterior_filter import carbon_conf   # noqa: E402


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'


def load_gff_chains(gff, contig, cds96=True, conf_filter=None):
    """unique exact chains from a GFF; conf_filter: None | 'carbon' (keep conf NA or >=0.1)"""
    conf = carbon_conf(gff) if conf_filter else None
    out = {}
    for tx in M29.load_txs(gff, contig):
        ex = [(a, b) for a, b, _ in tx['CDS']]
        if cds96 and sum(b - a for a, b in ex) < 96:
            continue
        if conf_filter:
            c = conf[tx['gene_id']]
            if not (c is None or c >= 0.1):
                continue
        out.setdefault(M29.key(tx), dict(nex=len(ex), exons=ex, strand=tx['strand']))
    return out


class Evaluator:
    def __init__(self, contigs):
        self.contigs = list(contigs)
        self.ref = {c: M29.build_reference(c) for c in self.contigs}
        self.refkeys, self.nex, self.unit_of_key, self.units = {}, {}, {}, []
        self.index = {}
        for c in self.contigs:
            keys = set()
            idx = []
            for t in self.ref[c]['refs']:
                k = M29.key(t)
                keys.add(k); self.nex[k] = len(t['CDS'])
                u = len(self.units); self.units.append(('gene', c, k)); self.unit_of_key[k] = u
                idx.append(dict(strand=t['strand'], exons=[(a, b) for a, b, _ in t['CDS']], id=u))
            self.refkeys[c] = keys
            self.index[c] = R4.Index(idx)
        self.n_gene_units = len(self.units)
        self.orphan_bin = {}
        self.methods = {}
        self._masks = {}

    def _unit_for_chain(self, c, k, v):
        a, b = v['exons'][0][0], v['exons'][-1][1]
        best, bu = 0, None
        for r in self.index[c].query(v['strand'], a, b):
            o = R4.ov_bp(r['exons'], v['exons'])
            if o > best:
                best, bu = o, r['id']
        if bu is not None:
            return bu
        key = (c, a // 100_000)
        if key not in self.orphan_bin:
            self.orphan_bin[key] = len(self.units) + len(self.orphan_bin)
        return self.orphan_bin[key]

    def add(self, name, chains_by_contig):
        pred = defaultdict(int); tp = defaultdict(int)
        strata = {s: [0, 0, 0] for s in ('1', '2-5', '>=6')}
        per = {}
        for c in self.contigs:
            ch = chains_by_contig[c]
            t = 0
            for k, v in ch.items():
                u = self._unit_for_chain(c, k, v)
                pred[u] += 1
                strata[strat(v['nex'])][1] += 1
                if k in self.refkeys[c]:
                    tp[self.unit_of_key[k]] += 1; t += 1; strata[strat(v['nex'])][0] += 1
            per[c] = dict(tp=t, pred=len(ch), ref=len(self.refkeys[c]))
        for k, n in self.nex.items():
            strata[strat(n)][2] += 1
        self.methods[name] = dict(pred=pred, tp=tp, strata=strata, per=per, chains=chains_by_contig)
        T = sum(p['tp'] for p in per.values()); P = sum(p['pred'] for p in per.values()); R = sum(p['ref'] for p in per.values())
        return dict(tp=T, pred=P, ref=R, F1=100 * 2 * T / (P + R))

    def arrays(self, name):
        m = self.methods[name]
        U = len(self.units) + len(self.orphan_bin)
        tp = np.zeros(U); pr = np.zeros(U)
        for u, v in m['pred'].items(): pr[u] = v
        for u, v in m['tp'].items(): tp[u] = v
        return tp, pr

    def bootstrap(self, pairs, B=2000, seed=0):
        """pairs: list of (a, b) method names -> CI of F1(a)-F1(b); also CI of each F1."""
        U = len(self.units) + len(self.orphan_bin)
        ref = np.zeros(U); ref[:self.n_gene_units] = 1
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, U, size=(B, U))
        F1 = {}
        for name in self.methods:
            tp, pr = self.arrays(name)
            T = tp[idx].sum(1); P = pr[idx].sum(1); R = ref[idx].sum(1)
            F1[name] = 100 * 2 * T / (P + R)
        out = {n: dict(F1=float(100 * 2 * self.arrays(n)[0].sum() / (self.arrays(n)[1].sum() + ref.sum())), lo=float(np.percentile(v, 2.5)), hi=float(np.percentile(v, 97.5))) for n, v in F1.items()}
        diffs = {}
        for a, b in pairs:
            d = F1[a] - F1[b]
            point = out[a]['F1'] - out[b]['F1']
            diffs[f'{a} - {b}'] = dict(point=point, lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5)), frac_positive=float((d > 0).mean()))
        return out, diffs

    def truth_masks(self, c):
        if c in self._masks:
            return self._masks[c]
        L = M29.LENGTH[c]; ref = self.ref[c]
        cds = {'+': np.zeros(L, bool), '-': np.zeros(L, bool)}
        for t in list(ref['refs']) + list(ref['isoforms']):
            for a, b, _ in t['CDS']:
                cds[t['strand']][a:b] = True
        gene = np.zeros(L, bool)
        for a, b in ref['gene_spans']:
            gene[max(0, a - 1):b] = True
        self._masks[c] = (cds, gene)
        return self._masks[c]

    def base_metrics(self, name):
        tp = pp = rp = bg = bgsize = 0
        for c in self.contigs:
            cds, gene = self.truth_masks(c)
            L = M29.LENGTH[c]
            pred = {'+': np.zeros(L, bool), '-': np.zeros(L, bool)}
            for v in self.methods[name]['chains'][c].values():
                for a, b in v['exons']:
                    pred[v['strand']][a:b] = True
            for st in '+-':
                tp += int((pred[st] & cds[st]).sum()); pp += int(pred[st].sum()); rp += int(cds[st].sum())
            anyp = pred['+'] | pred['-']
            bg += int((anyp & ~gene).sum()); bgsize += int((~gene).sum())
        return dict(base_P=100 * tp / max(pp, 1), base_R=100 * tp / max(rp, 1), base_F1=100 * 2 * tp / max(pp + rp, 1), bg_pred_kb=bg / 1000, bg_burden_pct=100 * bg / max(bgsize, 1))

    def table(self, name):
        m = self.methods[name]
        T = sum(p['tp'] for p in m['per'].values()); P = sum(p['pred'] for p in m['per'].values()); R = sum(p['ref'] for p in m['per'].values())
        st = {s: dict(tp=v[0], pred=v[1], ref=v[2], F1=100 * 2 * v[0] / max(v[1] + v[2], 1), recall=100 * v[0] / max(v[2], 1), precision=100 * v[0] / max(v[1], 1)) for s, v in m['strata'].items()}
        return dict(pooled=dict(tp=T, pred=P, ref=R, F1=100 * 2 * T / (P + R)), per_contig={c: dict(p, F1=100 * 2 * p['tp'] / (p['pred'] + p['ref'])) for c, p in m['per'].items()},
                    strata=st, base=self.base_metrics(name))
