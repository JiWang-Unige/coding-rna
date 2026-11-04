"""M38-VIPSANIA scorer: Vipsania GFF3 -> our chain representation -> the project's existing metric code (nothing of the metric is re-implemented).

Protocol: docs/review_20261016/28_animal_C2_and_Vipsania_proposal.md section 6 (C1 criterion; run order DEV + set B first).  Vipsania annotates the WHOLE nuclear genome
(vip_devb.sh); only the DEV chromosomes / set-B contigs are scored here, exactly like our own runs:

  DEV  (pooled exact-chain F1 of Arabidopsis NC_003074.8 + rice NC_089041.1, CDS >= 96 bp rule, M26 CDS-assessable primary references wholly in the region):
        per-region numbers from a1_eval.summarize (the M31/M34 DEV ruler), paired 1-Mb block bootstrap by running `m34_boot_pair.py` itself (2000 reps, seed 0) with the
        Vipsania run as the reference run against T2 / T10 / k50U / k50F (/ T2_s2 when present), plus its two built-in external comparators (GENERanno fresh-joint, Carbon-A filtered).
  set B (Dryas octopetala NC_141367.1/69/71/73 + Eschscholzia californica NW_028238309.1/461.1/483.1):
        `m34_p2b_score.run_scoring` (m30_score.Evaluator + load_gff_chains, same bootstrap) with the Vipsania run as `new` against T2 and T10 (and T2_s2 / the other Vipsania mode).

NOTE (protocol fact): the DEV chromosomes are NC_003074.8 and NC_089041.1 (m29_common.DEV_SCOPE, m31_common.DEV_SEQS; the numbers 70.107 / 77.406 are over them).  NC_003071.7 and
NC_089044.1 are the INTERNAL-VALIDATION chromosomes of M31/M34 (never DEV).

Two Vipsania modes are scored when given: `ft` (--finetune, the C1-criterion mode) and `pre` (pretrained only, descriptive).

Sub-commands (run with the generanno env for everything except `convert`, which is stdlib only):
  convert   --gff F --out DIR --contigs A,B [--slice slice.json]          normalise one GFF3 into DIR/<contig>/prediction/full_grammar.gff3 (the m34_infer layout all loaders read)
  dev       --out DIR [--ft KEY=GFF ...] [--pre KEY=GFF ...]               KEY in arabidopsis_thaliana, oryza_sativa
  setb      --out DIR [--ft KEY=GFF ...] [--pre KEY=GFF ...]               KEY in dryas_octopetala, eschscholzia_californica
  validate  --out DIR                                                      re-score already-scored external GFFs through the same converter + code path and compare
  smoke     --gff GFF --slice slice.json [--out JSON]                      PARSER CHECK on a TRAINING slice (not an evaluation)
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]                                           # .../experiments
OUTS = Path(os.path.expanduser('~/coding-rna/outputs'))
M34 = OUTS / 'M34-NTV2-SCALEUP'
DEV_CONTIGS = {'arabidopsis_thaliana': ['NC_003074.8'], 'oryza_sativa': ['NC_089041.1']}
B_CONTIGS = {'dryas_octopetala': ['NC_141367.1', 'NC_141369.1', 'NC_141371.1', 'NC_141373.1'],
             'eschscholzia_californica': ['NW_028238309.1', 'NW_028238461.1', 'NW_028238483.1']}
DEV_SEQS = tuple(c for cs in DEV_CONTIGS.values() for c in cs)
B_SEQS = tuple(c for cs in B_CONTIGS.values() for c in cs)
# existing arms (name -> run directory holding <contig>/prediction/full_grammar.gff3); only the complete ones are used
DEV_REF_RUNS = {'T2': M34 / 'pc2_T2/devM34_T2', 'T10': M34 / 'pc2_T10/devM34_T10', 'k50U': M34 / 'pc2_T10_k50U/devM34_T10_k50U',
                'k50F': M34 / 'pc2_T10_k50F/devM34_T10_k50F', 'T2_s2': M34 / 'pc2_T2_s2/devM34_T2_s2'}
B_REF_RUNS = {'T2': M34 / 'p2b/T2', 'T10': M34 / 'p2b/T10', 'T2_s2': M34 / 'pc2_T2_s2/setB'}
SOURCE = 'Vipsania'


# ======================================================================================================================================
# part 1 -- converter (stdlib only)
# ======================================================================================================================================
def parse_attrs(text):
    """GFF3 `k=v;k=v` (and GTF `k "v";`) attributes -> dict (same reading as the frozen ruler's parse_attrs)"""
    attrs = {}
    for item in text.strip().rstrip(';').split(';'):
        item = item.strip()
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
        else:
            k, v = item.split(None, 1)
        attrs[k] = v.strip().strip('"')
    return attrs


def read_gff(path, keep=None, slice_info=None):
    """Stream a GFF3 (gene / mRNA|transcript / CDS|exon with Parent links).  keep: set of output seqids to read (None = all).
    slice_info: dict(name, contig, start1): rows of seqid `name` are renamed to `contig` and shifted by start1 - 1 (smoke slices).
    -> {seqid: dict(genes={gid: dict(start, end, strand)}, tx={tid: dict(gene, strand, cds=[(start1, end1, phase)])})}   coordinates 1-based inclusive as in the file."""
    shift = (slice_info['start1'] - 1) if slice_info else 0
    pending = []                                                # CDS rows, assigned after every parent is known (children may precede parents)
    data = defaultdict(lambda: dict(genes={}, tx={}))
    with open(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                raise ValueError(f'{path}:{n}: expected 9 columns, got {len(f)}')
            seqid, _src, feature, s, e, _score, strand, phase, attrs_s = f
            if slice_info:
                if seqid != slice_info['name']:
                    continue
                seqid = slice_info['contig']
            if keep is not None and seqid not in keep:
                continue
            if feature not in ('gene', 'mRNA', 'transcript', 'CDS'):
                continue
            a = parse_attrs(attrs_s)
            s, e = int(s) + shift, int(e) + shift
            d = data[seqid]
            if feature == 'gene':
                gid = a.get('ID') or a.get('gene_id')
                if not gid:
                    raise ValueError(f'{path}:{n}: gene without ID')
                d['genes'][gid] = dict(start=s, end=e, strand=strand)
            elif feature in ('mRNA', 'transcript'):
                tid = a.get('ID') or a.get('transcript_id')
                gid = a.get('gene_id') or a.get('Parent') or tid
                if not tid or ',' in gid:
                    raise ValueError(f'{path}:{n}: ambiguous transcript id / parent')
                d['tx'][tid] = dict(gene=gid, strand=strand, cds=[])
            else:
                parents = (a.get('Parent') or a.get('transcript_id') or a.get('gene_id') or '').split(',')
                if not parents or not parents[0]:
                    raise ValueError(f'{path}:{n}: CDS without Parent')
                for p in parents:
                    pending.append((seqid, p, s, e, strand, phase))
    for seqid, p, s, e, strand, phase in pending:
        d = data[seqid]
        if p not in d['tx']:                                    # CDS parented directly on a gene (or an mRNA row missing): one transcript per parent id
            d['tx'][p] = dict(gene=p, strand=strand, cds=[])
        if d['tx'][p]['strand'] != strand:
            raise ValueError(f'{path}: transcript {p} has CDS rows on both strands')
        d['tx'][p]['cds'].append((s, e, phase))
    for d in data.values():
        for t in d['tx'].values():
            t['cds'].sort()
    return data


def chains_of(seq_data, seqid, cds96=True):
    """The frozen ruler's unique exact chains of one contig, from the parsed structure (independent re-derivation, used as a cross-check of the written files):
    one transcript per gene (longest total CDS, then smallest id), key = (seqid, strand, ((start0, end), ...))."""
    by_gene = defaultdict(list)
    for tid, t in seq_data['tx'].items():
        if t['cds']:
            by_gene[t['gene']].append((tid, t))
    out = {}
    for gid, txs in by_gene.items():
        tid, t = min(txs, key=lambda x: (-sum(e - s + 1 for s, e, _ in x[1]['cds']), x[0]))
        exons = tuple((s - 1, e) for s, e, _ in t['cds'])
        if cds96 and sum(b - a for a, b in exons) < 96:
            continue
        out.setdefault((seqid, t['strand'], exons), dict(nex=len(exons), exons=list(exons), strand=t['strand']))
    return out


def write_contig(seq_data, seqid, out_dir, gene_filter=None, source=SOURCE):
    """-> path of <out_dir>/<seqid>/prediction/full_grammar.gff3 (gene / mRNA / exon / CDS, Parent links, gene_biotype=protein_coding); returns (path, n_genes, n_tx)"""
    path = Path(out_dir) / seqid / 'prediction' / 'full_grammar.gff3'
    path.parent.mkdir(parents=True, exist_ok=True)
    by_gene = defaultdict(list)
    for tid, t in seq_data['tx'].items():
        if t['cds'] and (gene_filter is None or gene_filter(t['gene'])):
            by_gene[t['gene']].append((tid, t))
    rows, n_tx = [], 0
    for gid in sorted(by_gene, key=lambda g: (min(s for _, t in by_gene[g] for s, _, _ in t['cds']), g)):
        txs = sorted(by_gene[gid])
        strand = txs[0][1]['strand']
        gs = min(s for _, t in txs for s, _, _ in t['cds'])
        ge = max(e for _, t in txs for _, e, _ in t['cds'])
        rows.append(f'{seqid}\t{source}\tgene\t{gs}\t{ge}\t.\t{strand}\t.\tID={gid};gene_biotype=protein_coding')
        for tid, t in txs:
            n_tx += 1
            ts, te = t['cds'][0][0], t['cds'][-1][1]
            rows.append(f'{seqid}\t{source}\tmRNA\t{ts}\t{te}\t.\t{strand}\t.\tID={tid};Parent={gid}')
            for i, (s, e, ph) in enumerate(t['cds'], 1):
                rows.append(f'{seqid}\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\tID={tid}.exon{i};Parent={tid}')
                rows.append(f'{seqid}\t{source}\tCDS\t{s}\t{e}\t.\t{strand}\t{ph if ph in ("0", "1", "2") else "."}\tID={tid}.cds{i};Parent={tid}')
    with open(path, 'w') as h:
        h.write('##gff-version 3\n')
        h.write('\n'.join(rows) + ('\n' if rows else ''))
    return path, len(by_gene), n_tx


def convert(gff, out_dir, contigs, slice_info=None, gene_filter=None):
    """Normalise one GFF3 into the per-contig layout.  -> {contig: dict(genes, transcripts, cds_chains)}"""
    data = read_gff(gff, keep=set(contigs), slice_info=slice_info)
    report = {}
    for c in contigs:
        d = data.get(c, dict(genes={}, tx={}))
        p, ng, nt = write_contig(d, c, out_dir, gene_filter=gene_filter)
        report[c] = dict(file=str(p), genes=ng, transcripts=nt, chains=len(chains_of(d, c, cds96=False)), chains_cds96=len(chains_of(d, c, cds96=True)))
    return report


def sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_meta(gff):
    p = Path(gff).parent / 'run.json'
    return json.load(open(p)) if p.exists() else None


# ======================================================================================================================================
# part 2 -- the project's code (generanno env)
# ======================================================================================================================================
_P = {}


def project():
    """Import the frozen/project modules once (sys.path set up the way m34_boot_pair / m34_p2b_score do it)."""
    if _P:
        return _P
    sys.path.insert(0, str(EXP / 'M29-A1-RESIDUAL' / 'scripts'))
    import a1_common as K
    sys.path.insert(0, str(K.CA_SRC))
    import a1_eval as E
    import a1_compare as AC
    sys.path.insert(0, str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts'))
    sys.path.insert(0, str(EXP / 'M34-NTV2-SCALEUP' / 'scripts'))
    import m34_p2b_score as P2B
    import m29_common as M29
    assert tuple(E.DEV) == DEV_SEQS, (E.DEV, DEV_SEQS)           # the DEV chromosomes of the frozen ruler
    _P.update(K=K, E=E, AC=AC, P2B=P2B, M29=M29)
    return _P


def run_dir_ok(path, contigs):
    return all((Path(path) / c / 'prediction' / 'full_grammar.gff3').exists() for c in contigs)


def dev_summary(run_dir, cds96=True):
    """pooled / per-region / per-species / strata numbers of a run directory, by the DEV ruler (a1_eval.load_chains + a1_eval.summarize)"""
    P = project()
    E, AC = P['E'], P['AC']
    chains = {s: E.load_chains(str(run_dir), s, cds96=cds96) for s in DEV_SEQS}
    r = E.summarize(chains, AC.refs, DEV_SEQS)
    r['pooled']['precision'] = 100 * r['pooled']['tp'] / max(r['pooled']['pred'], 1)
    r['pooled']['recall'] = 100 * r['pooled']['tp'] / max(r['pooled']['ref'], 1)
    return r


def boot_pair(out_json, specs):
    """the paired 1-Mb block bootstrap of m34_boot_pair.py, run as is (first spec = reference run) -> its JSON"""
    env = dict(os.environ, M34_BACKEND=os.environ.get('M34_BACKEND', 'ntv2'))
    cmd = [sys.executable, str(EXP / 'M34-NTV2-SCALEUP' / 'scripts' / 'm34_boot_pair.py'), str(out_json)] + [f'{n}={d}' for n, d in specs]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=str(EXP / 'M34-NTV2-SCALEUP' / 'scripts'), env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'm34_boot_pair failed rc={r.returncode}\n{r.stderr[-3000:]}')
    res = json.load(open(out_json))
    res['_seconds'] = round(time.time() - t0, 1)
    return res


def parse_kv(items, allowed):
    out = {}
    for it in items or []:
        k, _, v = it.partition('=')
        assert k in allowed and v, f'expected KEY=GFF with KEY in {sorted(allowed)}, got {it!r}'
        out[k] = v
    return out


def refuse_overwrite(path, overwrite):
    if Path(path).exists() and not overwrite:
        raise SystemExit(f'{path} exists: DEV / set B are evaluated once per run (pre-registered); pass --overwrite only to redo a failed or deliberate re-scoring')


def cmd_dev(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result_path = out / 'vip_dev_result.json'
    refuse_overwrite(result_path, args.overwrite)
    variants = {'VIP_ft': parse_kv(args.ft, DEV_CONTIGS), 'VIP_pre': parse_kv(args.pre, DEV_CONTIGS)}
    variants = {n: g for n, g in variants.items() if g}
    assert variants, 'give --ft and/or --pre'
    for n, g in variants.items():
        assert set(g) == set(DEV_CONTIGS), f'{n}: need a GFF for every DEV genome {sorted(DEV_CONTIGS)}'
    P = project()
    runs, conv, meta = {}, {}, {}
    for name, g in variants.items():
        runs[name] = out / 'runs' / name
        conv[name] = {}
        meta[name] = {}
        for key, gff in g.items():
            conv[name].update(convert(gff, runs[name], DEV_CONTIGS[key]))
            meta[name][key] = dict(gff=str(gff), sha256=sha256(gff), run_json=read_meta(gff))
    refs = {n: d for n, d in DEV_REF_RUNS.items() if run_dir_ok(d, DEV_SEQS)}
    res = dict(label='M38-VIPSANIA DEV', protocol='docs/review_20261016/28_animal_C2_and_Vipsania_proposal.md section 6',
               ruler='a1_eval.load_chains (unique exact chains, CDS>=96 bp) + a1_eval.summarize; refs = M26 CDS-assessable primaries (DEV 7,728)',
               contigs=DEV_CONTIGS, converted=conv, inputs=meta, reference_runs={n: str(d) for n, d in refs.items()}, variants={})
    for name in variants:
        s = dev_summary(runs[name])
        s_nocds96 = dev_summary(runs[name], cds96=False)
        specs = [(name, runs[name])] + [(n, runs[m]) for n, m in ((o, o) for o in variants if o != name)] + [(n, d) for n, d in refs.items()]
        bp = boot_pair(out / f'boot_pair_{name}.json', specs)
        res['variants'][name] = dict(pooled=s['pooled'], per_region=s['per_contig'], by_species=s['by_species'], by_strata=s['by_strata'],
                                     pooled_no_cds96_rule=s_nocds96['pooled'], bootstrap=bp,
                                     paired_vs=dict({n: dict(point=bp[n]['diff_ref_minus_this_point'], CI95=bp[n]['diff_ref_minus_this_CI95'], F1_other=bp[n]['F1'])
                                                     for n in bp if isinstance(bp[n], dict) and 'diff_ref_minus_this_point' in bp[n]}))
    res['time'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    json.dump(res, open(result_path, 'w'), indent=1)
    for name, v in res['variants'].items():
        p = v['pooled']
        print(f"DEV {name}: pooled F1 {p['F1']:.3f} (tp {p['tp']} pred {p['pred']} ref {p['ref']}; P {p['precision']:.2f} R {p['recall']:.2f}) | " +
              ' | '.join(f"{c} {x['F1']:.2f}" for c, x in v['per_region'].items()) + f" | CI {v['bootstrap'][name]['CI95']}")
        print('   paired (this - other): ' + ' | '.join(f"{n} {d['point']:+.2f} [{d['CI95'][0]:+.2f},{d['CI95'][1]:+.2f}]" for n, d in v['paired_vs'].items()))
    return 0


def cmd_setb(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result_path = out / 'vip_setb_result.json'
    refuse_overwrite(result_path, args.overwrite)
    variants = {'VIP_ft': parse_kv(args.ft, B_CONTIGS), 'VIP_pre': parse_kv(args.pre, B_CONTIGS)}
    variants = {n: g for n, g in variants.items() if g}
    assert variants, 'give --ft and/or --pre'
    for n, g in variants.items():
        assert set(g) == set(B_CONTIGS), f'{n}: need a GFF for every set-B genome {sorted(B_CONTIGS)}'
    P = project()
    P2B = P['P2B']
    runs, conv, meta = {}, {}, {}
    for name, g in variants.items():
        runs[name] = out / 'runs' / name
        conv[name], meta[name] = {}, {}
        for key, gff in g.items():
            conv[name].update(convert(gff, runs[name], B_CONTIGS[key]))
            meta[name][key] = dict(gff=str(gff), sha256=sha256(gff), run_json=read_meta(gff))
    refs = {n: d for n, d in B_REF_RUNS.items() if run_dir_ok(d, B_SEQS)}
    assert 'T2' in refs and 'T10' in refs, f'T2 / T10 set-B runs missing: {refs}'
    res = dict(label='M38-VIPSANIA set B', protocol='docs/review_20261016/19_M34_scaleup_plan.md section 3.1 P2b scoring; 13_M30_protocol.md section 2 (set B)',
               contigs=B_CONTIGS, converted=conv, inputs=meta, reference_runs={n: str(d) for n, d in refs.items()}, variants={})
    for name in variants:
        extra = [(o, [str(runs[o])]) for o in variants if o != name] + [(n, [str(d)]) for n, d in refs.items() if n not in ('T2',)]
        r = P2B.run_scoring((name, [str(runs[name])]), ('T2', [str(refs['T2'])]), extra, label=f'M38 {name} vs T2/T10 (set B)')
        r.pop('verdict', None)                                  # the P2b criterion (T10 - T2 >= +2.0 ...) does not apply to Vipsania
        r.pop('m36_zero_shot_arm_Z', None)
        res['variants'][name] = dict(table=r['runs'][name]['table'], bootstrap_F1=r['runs'][name]['bootstrap_F1'], paired_new_minus=r['paired_new_minus'],
                                     others={n: dict(table=x['table'], bootstrap_F1=x['bootstrap_F1']) for n, x in r['runs'].items() if n != name},
                                     bootstrap=r['bootstrap'], ruler=r['ruler'], score_seconds=r['score_seconds'])
    res['time'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    json.dump(res, open(result_path, 'w'), indent=1)
    for name, v in res['variants'].items():
        t = v['table']
        print(f'set B {name}: ' + ' | '.join(f"{sc} F1 {t[sc]['F1']:.2f} (tp {t[sc]['tp']} pred {t[sc]['pred']} ref {t[sc]['ref']})" for sc in P2B.SCOPES) +
              ' | per contig: ' + ', '.join(f"{c} {x['F1']:.2f}" for c, x in t['per_contig'].items()))
        print('   paired (this - other, pooled): ' + ' | '.join(f"{n} {d['pooled']['point']:+.2f} [{d['pooled']['CI95'][0]:+.2f},{d['pooled']['CI95'][1]:+.2f}]" for n, d in v['paired_new_minus'].items()))
    return 0


# ======================================================================================================================================
# validation: re-score external GFFs that were scored before, through the converter and the same code paths
# ======================================================================================================================================
def close(a, b, tol=1e-9):
    return abs(a - b) <= tol


def cmd_validate(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    P = project()
    E, AC, P2B, M29 = P['E'], P['AC'], P['P2B'], P['M29']
    rep = dict(checks={}, notes=[])
    ok_all = True

    def check(name, got, want, tol=1e-9):
        nonlocal ok_all
        good = close(got, want, tol) if isinstance(want, float) else got == want
        ok_all &= bool(good)
        rep['checks'][name] = dict(got=got, expected=want, ok=bool(good))
        print(f"[{'OK ' if good else 'FAIL'}] {name}: got {got} expected {want}", flush=True)

    pub = json.load(open(M34 / 'boot_pair_T10_T2_DEV.json'))
    pub_k50u = json.load(open(M34 / 'boot_pair_k50U_DEV.json'))
    # ---- 1. GENERanno fresh-joint DEV GFFs (published 70.107): GFF3 of our own M28 writer -> converter -> DEV ruler
    fj = {s: AC.DEVGFF['B_freshjoint'][s] for s in DEV_SEQS}
    conv_fj = {}
    for s in DEV_SEQS:
        conv_fj.update(convert(fj[s], out / 'runs' / 'FJ', [s]))
    r_fj = dev_summary(out / 'runs' / 'FJ')                                         # CDS >= 96 bp rule on: the rule of our own arms (E.load_chains default)
    r_fj_all = dev_summary(out / 'runs' / 'FJ', cds96=False)                        # rule off: the path of the published external rows (AC.gff_chains)
    orig = {s: AC.gff_chains(fj[s], s) for s in DEV_SEQS}
    r_orig = E.summarize(orig, AC.refs, DEV_SEQS)
    check('fresh-joint via converter, rule off: pooled F1 == published GENERanno fresh-joint', r_fj_all['pooled']['F1'], pub['GENERanno_fresh_joint']['F1'])
    check('fresh-joint via converter, rule off: pooled F1 rounded (published 70.107)', round(r_fj_all['pooled']['F1'], 3), 70.107)
    check('fresh-joint original path (AC.gff_chains): pooled F1 == published', r_orig['pooled']['F1'], pub['GENERanno_fresh_joint']['F1'])
    orig96 = {s: {k: v for k, v in orig[s].items() if sum(b - a for a, b in v['exons']) >= 96} for s in DEV_SEQS}
    r_orig96 = E.summarize(orig96, AC.refs, DEV_SEQS)
    check('fresh-joint via converter, rule on == original chains filtered by the rule (rule effect +%.3f)' % (r_fj['pooled']['F1'] - r_orig['pooled']['F1']),
          r_fj['pooled']['F1'], r_orig96['pooled']['F1'])
    for s in DEV_SEQS:                                                              # converter keeps the chain sets identical
        a = E.load_chains(str(out / 'runs' / 'FJ'), s, cds96=False)
        check(f'fresh-joint {s}: chain set equal to original GFF', sorted(a) == sorted(orig[s]), True)
        mine = chains_of(read_gff(fj[s], keep={s})[s], s, cds96=False)             # independent stdlib derivation
        check(f'fresh-joint {s}: stdlib chain derivation equal to the ruler loader', sorted(mine) == sorted(orig[s]), True)
    # ---- 2. Carbon-A DEV GFFs (different dialect: one shared ID for all CDS rows; confidence= on the gene): unfiltered 75.24 and filtered 77.406
    from m29_r3_posterior_filter import carbon_conf                                 # noqa: E402  (K.CA_SRC is on the path)
    cu_dirs, cf_dirs = out / 'runs' / 'CarbonA_unfiltered', out / 'runs' / 'CarbonA_filtered'
    for s in DEV_SEQS:
        gff = AC.DEVGFF['CarbonA'][s]
        conf = carbon_conf(gff)
        convert(gff, cu_dirs, [s])
        convert(gff, cf_dirs, [s], gene_filter=lambda g, conf=conf: (conf[g] is None or conf[g] >= 0.1))
    ru = dev_summary(cu_dirs, cds96=False)
    rf = dev_summary(cf_dirs, cds96=False)
    cu0 = E.summarize({s: AC.carbon_dev(s, False) for s in DEV_SEQS}, AC.refs, DEV_SEQS)
    cf0 = E.summarize({s: AC.carbon_dev(s, True) for s in DEV_SEQS}, AC.refs, DEV_SEQS)
    check('Carbon-A unfiltered via converter == AC.carbon_dev(unfiltered)', ru['pooled']['F1'], cu0['pooled']['F1'])
    check('Carbon-A unfiltered pooled F1 rounded (published 75.24)', round(ru['pooled']['F1'], 2), 75.24)
    check('Carbon-A filtered via converter == published', rf['pooled']['F1'], pub['CarbonA_filtered']['F1'])
    check('Carbon-A filtered pooled F1 rounded (published 77.41)', round(rf['pooled']['F1'], 2), 77.41)
    check('Carbon-A filtered via converter == AC.carbon_dev(filtered)', rf['pooled']['F1'], cf0['pooled']['F1'])
    # ---- 3. the paired bootstrap through the same boot_pair code: converted runs against the arms of the published JSON
    specs = [('T10', DEV_REF_RUNS['T10']), ('T2', DEV_REF_RUNS['T2']), ('FJ', out / 'runs' / 'FJ'), ('CarbonA_filtered_conv', cf_dirs)]
    bp = boot_pair(out / 'boot_pair_validate.json', specs)
    check('bootstrap: blocks', bp['blocks'], pub['blocks'])
    check('bootstrap: T10 F1', bp['T10']['F1'], pub['T10']['F1'])
    check('bootstrap: T10 CI95', bp['T10']['CI95'], pub['T10']['CI95'])
    check('bootstrap: T2 F1', bp['T2']['F1'], pub['T2']['F1'])
    check('bootstrap: T10 - T2 point', bp['T2']['diff_ref_minus_this_point'], pub['T2']['diff_ref_minus_this_point'])
    check('bootstrap: T10 - T2 CI95', bp['T2']['diff_ref_minus_this_CI95'], pub['T2']['diff_ref_minus_this_CI95'])
    check('bootstrap: FJ (converted, rule on) point F1 == in-process E.summarize of the same run', bp['FJ']['F1'], r_fj['pooled']['F1'])
    check('bootstrap: Carbon-A filtered (converted) F1 == built-in', bp['CarbonA_filtered_conv']['F1'], bp['CarbonA_filtered']['F1'])
    check('bootstrap: T10 - Carbon-A filtered (converted) CI95 == published', bp['CarbonA_filtered_conv']['diff_ref_minus_this_CI95'], pub['CarbonA_filtered']['diff_ref_minus_this_CI95'])
    # a published run with the pooled number of k50U (so that a non-T10 arm also goes through the converter): k50U GFF -> converter -> same pooled F1
    conv_k = {}
    for s in DEV_SEQS:
        conv_k.update(convert(Path(DEV_REF_RUNS['k50U']) / s / 'prediction' / 'full_grammar.gff3', out / 'runs' / 'k50U_conv', [s]))
    check('k50U GFF via converter: pooled F1 == published', dev_summary(out / 'runs' / 'k50U_conv')['pooled']['F1'], pub_k50u['k50U']['F1'])
    # ---- 4. set B: the stored T10 / T2 B-set GFFs through the converter -> run_scoring -> the stored p2b_result.json
    p2b = json.load(open(M34 / 'p2b' / 'p2b_result.json'))
    for name in ('T10', 'T2'):
        for contig in B_SEQS:
            convert(B_REF_RUNS[name] / contig / 'prediction' / 'full_grammar.gff3', out / 'runs' / f'B_{name}_conv', [contig])
    new = ('T10c', [str(out / 'runs' / 'B_T10_conv')])
    old = ('T2c', [str(out / 'runs' / 'B_T2_conv')])
    r = P2B.run_scoring(new, old, label='validate')
    for sc in P2B.SCOPES:
        check(f'set B T10 via converter pooled/species F1 [{sc}]', r['runs']['T10c']['table'][sc]['F1'], p2b['runs']['T10']['table'][sc]['F1'])
        check(f'set B T2 via converter F1 [{sc}]', r['runs']['T2c']['table'][sc]['F1'], p2b['runs']['T2']['table'][sc]['F1'])
    check('set B T10 - T2 pooled point', r['paired_new_minus']['T2c']['pooled']['point'], p2b['paired_new_minus']['T2']['pooled']['point'])
    check('set B T10 - T2 pooled CI95', r['paired_new_minus']['T2c']['pooled']['CI95'], p2b['paired_new_minus']['T2']['pooled']['CI95'])
    rep.update(ok=bool(ok_all), dev_fresh_joint_pooled=r_fj['pooled'], dev_carbon_unfiltered_pooled=ru['pooled'], dev_carbon_filtered_pooled=rf['pooled'],
               setb_T10_pooled=r['runs']['T10c']['table']['pooled'], setb_T2_pooled=r['runs']['T2c']['table']['pooled'],
               time=time.strftime('%Y-%m-%dT%H:%M:%S'))
    json.dump(rep, open(out / 'validation.json', 'w'), indent=1, default=str)
    print('VALIDATION', 'PASSED' if ok_all else 'FAILED', f'({sum(c["ok"] for c in rep["checks"].values())}/{len(rep["checks"])} checks)')
    return 0 if ok_all else 1


# ======================================================================================================================================
# smoke: parser check on a TRAINING slice (never an evaluation)
# ======================================================================================================================================
def cmd_smoke(args):
    sl = json.load(open(args.slice))
    contig, a, b = sl['contig'], sl['start1'] - 1, sl['end1']                     # [a, b) 0-based half-open window in contig coordinates
    out = Path(args.out_dir)
    rep = convert(args.gff, out, [contig], slice_info=sl)
    P = project()
    E, AC = P['E'], P['AC']
    pred_all = E.load_chains(str(out), contig)
    inside = lambda k: k[2][0][0] >= a and k[2][-1][1] <= b
    pred = {k: v for k, v in pred_all.items() if inside(k)}
    ref = {k for k in AC.refs.keys[contig] if inside(k)}
    tp = sum(1 for k in pred if k in ref)
    res = dict(label='PARSER CHECK ONLY (training slice, not an evaluation)', slice=sl, converted=rep, chains_in_slice=len(pred), chains_total=len(pred_all), reference_in_slice=len(ref),
               tp=tp, precision=100 * tp / max(len(pred), 1), recall=100 * tp / max(len(ref), 1), F1=200 * tp / max(len(pred) + len(ref), 1),
               stdlib_chain_derivation_equal=sorted(chains_of(read_gff(args.gff, slice_info=sl)[contig], contig)) == sorted(pred_all))
    res['gff_features'] = feature_census(args.gff)
    if args.out:
        json.dump(res, open(args.out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'gff_features'}))
    return 0


def feature_census(gff):
    from collections import Counter
    c, parents, phases = Counter(), Counter(), Counter()
    with open(gff) as h:
        for line in h:
            if line[0] == '#' or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            c[f[2]] += 1
            if f[2] == 'CDS':
                phases[f[7]] += 1
                parents['multi' if ',' in f[8] else 'single'] += 1
    return dict(features=dict(c), cds_phase=dict(phases), cds_parent=dict(parents))


def cmd_convert(args):
    sl = json.load(open(args.slice)) if args.slice else None
    rep = convert(args.gff, args.out, args.contigs.split(','), slice_info=sl)
    print(json.dumps(rep, indent=1))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('convert'); p.add_argument('--gff', required=True); p.add_argument('--out', required=True); p.add_argument('--contigs', required=True); p.add_argument('--slice', default=None)
    p.set_defaults(fn=cmd_convert)
    for name, fn in (('dev', cmd_dev), ('setb', cmd_setb)):
        p = sub.add_parser(name)
        p.add_argument('--out', required=True)
        p.add_argument('--ft', nargs='+', default=None, metavar='KEY=GFF', help='--finetune annotation (the C1-criterion mode)')
        p.add_argument('--pre', nargs='+', default=None, metavar='KEY=GFF', help='pretrained-only annotation (descriptive)')
        p.add_argument('--overwrite', action='store_true')
        p.set_defaults(fn=fn)
    p = sub.add_parser('validate'); p.add_argument('--out', required=True); p.set_defaults(fn=cmd_validate)
    p = sub.add_parser('smoke'); p.add_argument('--gff', required=True); p.add_argument('--slice', required=True); p.add_argument('--out-dir', required=True); p.add_argument('--out', default=None)
    p.set_defaults(fn=cmd_smoke)
    args = ap.parse_args()
    sys.exit(args.fn(args) or 0)


if __name__ == '__main__':
    main()
