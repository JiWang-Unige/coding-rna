"""Parity gate scoring + difference attribution.

Scores outputs/M29-CARBON-A-GAP/R6/<variant>/<seq>/full_grammar.gff3 with the frozen ruler, optionally applying Carbon-A's CDS>=96 bp
output rule, and explains every Carbon-A-unfiltered chain that the gate decoder did NOT reproduce (grammar differences recorded in 06 §4.1):
  incomplete (record-edge truncated) | N_in_CDS | intron<5 bp | AT-AC intron | start/stop codon split by an intron | >128 exons | other
and the gate-only chains (not in Carbon-A).
Usage: m29_gate_diff.py <variant> [seq,seq,...]
"""
import json
import sys
from collections import Counter
import m29_common as M
from m29_load_preds import tool_chains

COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def cds_len(c):
    return sum(b - a for a, b in c['exons'])


def reasons(c, seqstr, status):
    r = []
    if status.get(c['id']) == 'incomplete':
        r.append('incomplete_record_edge_or_short_intron')
    if len(c['exons']) > 128:
        r.append('over_128_exons')
    ex = c['exons'] if c['strand'] == '+' else c['exons'][::-1]
    lens = [b - a for a, b in ex]
    if lens[0] < 3:
        r.append('start_codon_split_by_intron')
    if lens[-1] < 3:
        r.append('stop_codon_split_by_intron')
    for (a, b), (c2, d) in zip(c['exons'], c['exons'][1:]):
        s = seqstr[b:c2]
        s = s if c['strand'] == '+' else rc(s)
        if len(s) < 5:
            r.append('intron<5bp')
        if s[:2] == 'AT' and s[-2:] == 'AC':
            r.append('AT-AC_intron')
        elif s[:2] not in ('GT', 'GC') or s[-2:] != 'AG':
            r.append('other_noncanonical_intron')
    cds = ''.join(seqstr[a:b] for a, b in c['exons'])
    if 'N' in cds:
        r.append('N_in_CDS')
    return sorted(set(r)) or ['other']


def status_map(seq):
    st = {}
    for line in open(M.gff_path('CarbonA_unfiltered', seq)):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'gene':
            a = dict(x.split('=', 1) for x in f[8].split(';'))
            st[a['ID']] = a['status']
    return st


def main():
    variant = sys.argv[1]
    seqs = sys.argv[2].split(',') if len(sys.argv) > 2 else [s for _, s, _ in M.SCOPE]
    tot = {k: [0, 0, 0] for k in ('gate_raw', 'gate_cds96', 'carbonA_unfiltered')}
    perc, why_only_ca, only_gate_cnt, only_gate_len = {}, Counter(), 0, Counter()
    ident = ca_only = gate_only = 0
    for sp, seq, L in M.SCOPE:
        if seq not in seqs:
            continue
        rk = {M.key(t) for t in M.build_reference(seq)['refs']}
        gff = M.OUTS / 'M29-CARBON-A-GAP/R6' / variant / seq / 'full_grammar.gff3'
        gate = {}
        for tx in M.load_txs(gff, seq):
            gate.setdefault(M.key(tx), dict(key=M.key(tx), strand=tx['strand'], exons=[(a, b) for a, b, _ in tx['CDS']], id=tx['id']))
        ca = tool_chains('CarbonA_unfiltered', seq)
        gstat = status_map(seq)
        seqstr = M.read_sequence(seq)
        stat2 = {v['id']: gstat.get('gene' + v['id'][4:10], 'complete') for v in ca.values()}
        for name, d in (('gate_raw', gate), ('gate_cds96', {k: v for k, v in gate.items() if cds_len(v) >= 96}),
                        ('carbonA_unfiltered', ca)):
            tp = sum(1 for k in d if k in rk)
            tot[name][0] += tp; tot[name][1] += len(d); tot[name][2] += len(rk)
        g96 = {k for k, v in gate.items() if cds_len(v) >= 96}
        ident += len(set(ca) & g96)
        for k in set(ca) - g96:
            for rs in [reasons(ca[k], seqstr, stat2)]:
                why_only_ca['+'.join(rs)] += 1
        for k in g96 - set(ca):
            only_gate_cnt += 1
        gate_only += len(g96 - set(ca)); ca_only += len(set(ca) - g96)
        hits_ca_only = len((set(ca) - g96) & rk); hits_gate_only = len((g96 - set(ca)) & rk)
        perc[seq] = dict(gate_chains=len(gate), gate_cds96=len(g96), carbonA=len(ca), identical=len(set(ca) & g96),
                         only_carbon=len(set(ca) - g96), only_gate=len(g96 - set(ca)), tp_only_carbon=hits_ca_only, tp_only_gate=hits_gate_only,
                         gate_tp=sum(1 for k in g96 if k in rk), carbonA_tp=sum(1 for k in ca if k in rk))
        print(seq, perc[seq], flush=True)
    out = {name: dict(tp=v[0], pred=v[1], ref=v[2], F1=100 * 2 * v[0] / (v[1] + v[2])) for name, v in tot.items()}
    res = dict(variant=variant, pooled=out, per_contig=perc, why_carbon_chain_not_reproduced=dict(why_only_ca), identical=ident,
               only_carbon=ca_only, only_gate=gate_only)
    (M.M29 / 'R6' / f'diff_{variant}.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
