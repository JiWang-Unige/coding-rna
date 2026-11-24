"""M42 gate-1 scoring on the evaluation chromosomes E (protocol 35 §4.3-§4.4).  The ONLY place that reads the NEW release on E (and the OLD release on E, to class the
reference genes).  Ruler: `m30_score.Evaluator` (M26 CDS-assessable primaries of the NEW release, exact CDS chain) with `load_gff_chains` (unique exact chains, CDS >= 96 bp,
no posterior filter) -- the M36 / P2b counting.
Reference-gene classes (gate-0 rules of protocol 35 §3.4, applied on E to the ruler's primaries):
  novel     : no CDS base overlaps any OLD-release CDS on the contig (either strand)   -> missing-locus recall; split by miniprot evidence (coverage >= .5, identity >= .4)
  identical : the exact CDS chain is a CDS chain of the OLD release                     -> did adaptation break correct old annotation
  revised   : everything else
False-positive burden: predicted chains whose CDS span overlaps no NEW-release gene (any biotype, either strand).
Paired 1-Mb block bootstrap (block = (contig, first CDS start // 1 Mb); 2,000 draws, default_rng(0)) of U minus every other run, for F1 and missing-locus recall.
usage: m42_score.py <species_dir> --run NAME=DIR [--run ...] --out RESULT.json
       (species_dir: m42_config.json, new / old GFF, miniprot GFF, gate0/split.json; the scoring tree must have data/m1_screen/<species>/genome.fa + lengths.json)
"""
import argparse
import bisect
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M30-HOLDOUT-ADAPT' / 'scripts'))
sys.path.append(str(HERE))
import m42_common as Z                                   # noqa: E402

BLK, DRAWS, SEED = 1_000_000, 2000, 0
EV_COV, EV_ID = 0.50, 0.40


def opn(p):
    p = str(p)
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def attrs(f):
    return dict(kv.split('=', 1) for kv in f.strip().split(';') if '=' in kv)


def write_reference(new_gff, contigs, path):
    """NEW release restricted to the evaluation contigs -> data/m1_screen/<species>/reference.gff3 (what the frozen ruler reads)"""
    keep = set(contigs)
    with opn(new_gff) as f, open(path, 'w') as o:
        o.write('##gff-version 3\n')
        for line in f:
            if line.startswith('#'):
                continue
            if line.split('\t', 1)[0] in keep:
                o.write(line)


def old_cds(old_gff, contigs):
    """OLD release on E: merged CDS intervals per contig (0-based half-open) and the set of exact CDS chains (contig, strand, ((a, b), ...))."""
    keep = set(contigs)
    tx2gene, gene_pc, strand, cds = {}, set(), {}, defaultdict(list)
    for line in opn(old_gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] not in keep:
            continue
        a = attrs(f[8])
        if f[2] == 'CDS' and 'Parent' in a:
            cds[(f[0], a['Parent'])].append((int(f[3]) - 1, int(f[4])))
            strand[(f[0], a['Parent'])] = f[6]
    merged, chains = defaultdict(list), set()
    for (c, tx), ex in cds.items():
        ex = tuple(sorted(ex))
        chains.add((c, strand[(c, tx)], ex))
        merged[c] += list(ex)
    out = {}
    for c, iv in merged.items():
        m = []
        for a, b in sorted(iv):
            if m and a <= m[-1][1]:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        out[c] = ([x for x, _ in m], [y for _, y in m])
    return out, chains


def overlaps(iv, a, b):
    if iv is None:
        return False
    s, e = iv
    i = bisect.bisect_right(s, b - 1) - 1
    return i >= 0 and e[i] > a


def evidence(mp_gff, plen, contigs):
    keep = set(contigs)
    ev = defaultdict(list)
    for line in opn(mp_gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'mRNA' or f[0] not in keep:
            continue
        a = attrs(f[8])
        t = a.get('Target', '').split()
        if len(t) < 3 or t[0] not in plen:
            continue
        if (int(t[2]) - int(t[1]) + 1) / plen[t[0]] >= EV_COV and float(a.get('Identity', 0)) >= EV_ID:
            ev[f[0]].append((int(f[3]) - 1, int(f[4])))
    out = {}
    for c, iv in ev.items():
        iv.sort()
        out[c] = ([x for x, _ in iv], [y for _, y in iv], np.maximum.accumulate([y for _, y in iv]).tolist())
    return out


def ev_hit(ev, c, a, b):
    if c not in ev:
        return False
    s, _, emax = ev[c]
    i = bisect.bisect_right(s, b - 1) - 1
    return i >= 0 and emax[i] > a


def protein_lengths(q):
    plen, name = {}, None
    for line in opn(q):
        if line.startswith('>'):
            name = line[1:].split()[0]
            plen[name] = 0
        else:
            plen[name] += len(line.strip())
    return plen


def block(c, k):
    return (c, k[2][0][0] // BLK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('species_dir')
    ap.add_argument('--run', action='append', required=True, help='NAME=DIR (DIR/<contig>/prediction/full_grammar.gff3)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    D = Path(args.species_dir)
    cfg = json.load(open(D / 'm42_config.json'))
    E = json.load(open(D / 'gate0' / 'split.json'))['E']
    ref_path = Z.M1 / 'reference.gff3'
    if not ref_path.exists():
        write_reference(D / cfg['new_gff'], E, ref_path)
    import m30_score as SC                               # registers nothing for this species: register it here
    Z.register(SC.M29)
    ev_ = SC.Evaluator(E)
    runs = dict(x.split('=', 1) for x in args.run)
    chains = {n: {c: SC.load_gff_chains(str(Path(d) / c / 'prediction' / 'full_grammar.gff3'), c) for c in E} for n, d in runs.items()}
    oc, ochains = old_cds(D / cfg['old_gff'], E)
    ev = evidence(D / cfg['miniprot_gff'], protein_lengths(cfg['queries']), E)
    cls = {}
    for c in E:
        for k in ev_.refkeys[c]:
            ex = k[2]
            if not any(overlaps(oc.get(c), a, b) for a, b in ex):
                cls[k] = 'novel_ev' if ev_hit(ev, c, ex[0][0], ex[-1][1]) else 'novel_noev'
            elif (c, k[1], tuple(ex)) in ochains:
                cls[k] = 'identical'
            else:
                cls[k] = 'revised'
    genes = {}
    for c in E:
        sp = sorted((int(a), int(b)) for a, b in ev_.ref[c]['gene_spans'])
        genes[c] = ([a for a, _ in sp], np.maximum.accumulate([b for _, b in sp]).tolist())
    counts = defaultdict(int)
    for v in cls.values():
        counts[v] += 1
    res = dict(contigs=E, runs=runs, ref_classes=dict(counts), novel=counts['novel_ev'] + counts['novel_noev'], per_run={})
    # per-block vectors: overall (tp, pred, ref) and novel (tp, ref)
    bl = set()
    for c in E:
        bl.update(block(c, k) for k in ev_.refkeys[c])
        for ch in chains.values():
            bl.update(block(c, k) for k in ch[c])
    bl = sorted(bl)
    idx = {b: i for i, b in enumerate(bl)}
    nb = len(bl)
    rf, nrf = np.zeros(nb), np.zeros(nb)
    for c in E:
        for k in ev_.refkeys[c]:
            rf[idx[block(c, k)]] += 1
            nrf[idx[block(c, k)]] += cls[k].startswith('novel')
    vec = {}
    for n, ch in chains.items():
        t = ev_.add(n, ch)
        tp, pr, ntp = np.zeros(nb), np.zeros(nb), np.zeros(nb)
        rec = defaultdict(int)
        fp_burden = 0
        for c in E:
            for k, v in ch[c].items():
                i = idx[block(c, k)]
                pr[i] += 1
                a, b = v['exons'][0][0], v['exons'][-1][1]
                s, emax = genes[c]
                j = bisect.bisect_right(s, b - 1) - 1
                if not (j >= 0 and emax[j] > a):
                    fp_burden += 1
                if k in ev_.refkeys[c]:
                    tp[i] += 1
                    rec[cls[k]] += 1
                    if cls[k].startswith('novel'):
                        ntp[i] += 1
        vec[n] = (tp, pr, ntp)
        nov = rec['novel_ev'] + rec['novel_noev']
        res['per_run'][n] = dict(F1=t['F1'], precision=100 * t['tp'] / max(t['pred'], 1), recall=100 * t['tp'] / max(t['ref'], 1), tp=t['tp'], pred=t['pred'], ref=t['ref'],
                                 novel_recall=100 * nov / max(res['novel'], 1), novel_tp=nov,
                                 recall_by_class={k: 100 * rec[k] / max(counts[k], 1) for k in counts}, fp_no_gene_overlap=fp_burden,
                                 strata=ev_.table(n)['strata'])
    rng = np.random.default_rng(SEED)
    f1 = lambda v, ii: 200 * v[0][ii].sum() / max(v[1][ii].sum() + rf[ii].sum(), 1)
    nr = lambda v, ii: 100 * v[2][ii].sum() / max(nrf[ii].sum(), 1)
    dr = {n: (np.empty(DRAWS), np.empty(DRAWS)) for n in vec}
    for d in range(DRAWS):
        ii = rng.integers(0, nb, nb)
        for n, v in vec.items():
            dr[n][0][d], dr[n][1][d] = f1(v, ii), nr(v, ii)
    full = np.arange(nb)
    if 'U' in vec:
        res['paired_U_minus'] = {}
        for n in vec:
            if n == 'U':
                continue
            dF, dN = dr['U'][0] - dr[n][0], dr['U'][1] - dr[n][1]
            res['paired_U_minus'][n] = dict(F1=f1(vec['U'], full) - f1(vec[n], full), F1_CI95=[float(np.percentile(dF, 2.5)), float(np.percentile(dF, 97.5))],
                                            novel_recall=nr(vec['U'], full) - nr(vec[n], full), novel_recall_CI95=[float(np.percentile(dN, 2.5)), float(np.percentile(dN, 97.5))])
    if all(x in res['per_run'] for x in ('U', 'Z', 'F')):           # protocol 35 §9 (2026-11-24, user decision m42-L-arm): L dropped, rule = U vs Z and F
        p = res['per_run']
        rule = dict(U_gt_Z=p['U']['F1'] > p['Z']['F1'], U_gt_F=p['U']['F1'] > p['F']['F1'],
                    novel_U_gt_Z=p['U']['novel_recall'] > p['Z']['novel_recall'],
                    precision_not_below_Z_minus_2=p['U']['precision'] >= p['Z']['precision'] - 2.0)
        res['gate1_rule'] = dict(rule, continue_to_gate2=all(rule.values()))
    res['blocks'] = nb
    json.dump(res, open(args.out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'per_run'}, indent=1))
    for n, r in res['per_run'].items():
        print(n, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if k not in ('strata', 'recall_by_class')}, {k: round(v, 1) for k, v in r['recall_by_class'].items()})


if __name__ == '__main__':
    main()
