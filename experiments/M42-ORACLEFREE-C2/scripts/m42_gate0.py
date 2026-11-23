"""M42 gate 0 (protocol 35 §2-§3): train/eval split, positives from the OLD annotation, frozen mask variants V1/V2/V3, matched random masks,
oracle mask O, and the mask audit on training chromosomes A only.  The NEW annotation is read only for contigs in A (never for E).

usage: m42_gate0.py <species_dir>     (expects chrom12.fa / genome fasta, old.gff(.gz), new.gff(.gz), miniprot.gff, queries.faa(.gz), contigs.txt)
       config is read from <species_dir>/m42_config.json: {"genome": ..., "old_gff": ..., "new_gff": ..., "miniprot_gff": ..., "queries": ..., "contigs": [...]}
Output: <species_dir>/gate0/{split.json, positives.bed, mask_V1.bed, mask_V2.bed, mask_V3.bed, audit.json}
"""
import gzip, json, random, re, sys
from collections import defaultdict
from pathlib import Path

D = Path(sys.argv[1]); C = json.load(open(D / 'm42_config.json')); OUT = D / 'gate0'; OUT.mkdir(exist_ok=True)
SEED = 20261124
FLANK_POS = 1000                    # positives' span +-1 kb never unknown (§3.3)
EV_COV, EV_ID = 0.50, 0.40          # evidence hit (§3.2)
EXPAND = {'V1': 2000, 'V2': 10000}  # §3.3
STOPS = {'TAA', 'TAG', 'TGA'}
CODON = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
                                                 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}


def opn(p):
    p = str(p); return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def read_fasta(p, keep):
    seqs, cur, buf = {}, None, []
    for l in opn(p):
        if l.startswith('>'):
            if cur in keep: seqs[cur] = ''.join(buf).upper()
            cur, buf = l[1:].split()[0], []
        elif cur in keep: buf.append(l.strip())
    if cur in keep: seqs[cur] = ''.join(buf).upper()
    return seqs


def rc(s):
    return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))


def attrs(f):
    return dict(kv.split('=', 1) for kv in f.strip().split(';') if '=' in kv)


def read_gff(p, contigs):
    """-> {gene_id: dict(contig, strand, chains=[(tx_id, tuple(sorted CDS (s,e)) , partial)])} for protein-coding genes on contigs"""
    genes, tx2gene, cds = {}, {}, defaultdict(list); partial = set()
    for l in opn(p):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] not in contigs: continue
        a = attrs(f[8])
        if f[2] == 'gene' and a.get('gene_biotype') == 'protein_coding':
            genes[a['ID']] = dict(contig=f[0], strand=f[6], txs=[])
        elif f[2] in ('mRNA', 'transcript') and 'Parent' in a:
            tx2gene[a['ID']] = a['Parent']
        elif f[2] == 'CDS' and 'Parent' in a:
            cds[a['Parent']].append((int(f[3]), int(f[4])))
            if a.get('partial') == 'true' or 'partial' in a.get('Note', '').lower(): partial.add(a['Parent'])
    for tx, ex in cds.items():
        g = tx2gene.get(tx)
        if g in genes: genes[g]['txs'].append((tx, tuple(sorted(ex)), tx in partial))
    return {g: v for g, v in genes.items() if v['txs']}


def rep(g):
    """representative transcript = longest CDS (ties: first)"""
    return max(g['txs'], key=lambda t: sum(e - s + 1 for s, e in t[1]))


def cds_seq(seq, ex, strand):
    s = ''.join(seq[a - 1:b] for a, b in ex)
    return rc(s) if strand == '-' else s


def valid_positive(seq, g, tx):
    _, ex, part = tx
    if part: return 'partial'
    s = cds_seq(seq, ex, g['strand'])
    if len(s) < 96: return 'short'
    if len(s) % 3: return 'frame'
    if 'N' in s: return 'gap'
    if not s.startswith('ATG') or s[-3:] not in STOPS: return 'ends'
    if any(s[i:i + 3] in STOPS for i in range(0, len(s) - 3, 3)): return 'internal_stop'
    for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
        intr = seq[b1:a2 - 1]
        if g['strand'] == '-': intr = rc(intr)
        if (intr[:2], intr[-2:]) not in (('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')): return 'splice'
    return None


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return [tuple(x) for x in out]


def subtract(iv, cut):
    out, cut = [], merge(cut); j = 0
    for s, e in merge(iv):
        cur = s
        while j < len(cut) and cut[j][1] < cur: j += 1
        k = j
        while k < len(cut) and cut[k][0] <= e:
            if cut[k][0] > cur: out.append((cur, cut[k][0] - 1))
            cur = max(cur, cut[k][1] + 1); k += 1
        if cur <= e: out.append((cur, e))
    return out


def blen(iv):
    return sum(e - s + 1 for s, e in iv)


def covered(iv_sorted, s, e):
    """bases of [s,e] inside merged sorted intervals"""
    n = 0
    for a, b in iv_sorted:
        if b < s: continue
        if a > e: break
        n += min(b, e) - max(a, s) + 1
    return n


# ---- split (§2): primary chromosomes by length desc; odd ranks = A, even = E
contigs = C['contigs']
seqs = read_fasta(D / C['genome'], set(contigs))
order = sorted(contigs, key=lambda c: -len(seqs[c]))
A = [c for i, c in enumerate(order) if i % 2 == 0]; E = [c for i, c in enumerate(order) if i % 2 == 1]
json.dump(dict(order=order, lengths={c: len(seqs[c]) for c in order}, A=A, E=E, rule='protocol 35 §2'), open(OUT / 'split.json', 'w'), indent=1)
print('A', A, '\nE', E, flush=True)

# ---- positives from OLD annotation on A (§3.1)
old = read_gff(D / C['old_gff'], set(A))
pos, why = [], defaultdict(int)
for gid, g in old.items():
    tx = rep(g); r = valid_positive(seqs[g['contig']], g, tx)
    if r: why[r] += 1; continue
    pos.append((g['contig'], g['strand'], gid, tx[1]))
print('old protein-coding genes on A', len(old), 'positives', len(pos), dict(why), flush=True)
with open(OUT / 'positives.bed', 'w') as o:
    for c, st, gid, ex in pos:
        for a, b in ex: o.write(f'{c}\t{a - 1}\t{b}\t{gid}\t0\t{st}\n')
protect = defaultdict(list)   # positive gene spans +-1 kb (never unknown)
for c, st, gid, ex in pos: protect[c].append((max(1, ex[0][0] - FLANK_POS), min(len(seqs[c]), ex[-1][1] + FLANK_POS)))
protect = {c: merge(v) for c, v in protect.items()}
nonpos = {c: subtract([(1, len(seqs[c]))], protect.get(c, [])) for c in A}
nonpos_bp = sum(blen(v) for v in nonpos.values())

# ---- evidence (§3.2): miniprot mRNA records, coverage >= 0.5 and identity >= 0.4
plen = {}
name = None
for l in opn(C['queries']):
    if l.startswith('>'): name = l[1:].split()[0]; plen[name] = 0
    else: plen[name] += len(l.strip())
ev = defaultdict(list); n_ev = 0
for l in opn(D / C['miniprot_gff']):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'mRNA' or f[0] not in A: continue
    a = attrs(f[8]); t = a.get('Target', '').split()
    if len(t) < 3 or t[0] not in plen: continue
    cov = (int(t[2]) - int(t[1]) + 1) / plen[t[0]]; ident = float(a.get('Identity', 0))
    if cov >= EV_COV and ident >= EV_ID: ev[f[0]].append((int(f[3]), int(f[4]))); n_ev += 1
ev = {c: merge(v) for c, v in ev.items()}
print('evidence hits on A', n_ev, 'merged spans', sum(len(v) for v in ev.values()), flush=True)

# ---- frozen masks (§3.3)
masks = {}
for k, w in EXPAND.items():
    masks[k] = {c: merge(subtract([(max(1, s - w), min(len(seqs[c]), e + w)) for s, e in ev.get(c, [])], protect.get(c, []))) for c in A}
masks['V3'] = nonpos
for k in ('V1', 'V2', 'V3'):
    with open(OUT / f'mask_{k}.bed', 'w') as o:
        for c in A:
            for s, e in masks[k][c]: o.write(f'{c}\t{s - 1}\t{e}\tunknown_{k}\n')


def random_like(mask, seed):
    """same total bp and same segment-length multiset, placed uniformly in non-positive space (no overlap check beyond merging)"""
    rnd = random.Random(seed); segs = [e - s + 1 for c in A for s, e in mask[c]]
    space = [(c, s, e) for c in A for s, e in nonpos[c]]; wts = [e - s + 1 for _, s, e in space]; out = defaultdict(list)
    for L in segs:
        for _ in range(100):
            c, s, e = rnd.choices(space, wts)[0]
            if e - s + 1 >= L: x = rnd.randint(s, e - L + 1); out[c].append((x, x + L - 1)); break
    return {c: merge(out[c]) for c in A}


rand = {k: random_like(masks[k], SEED + i) for i, k in enumerate(('V1', 'V2'))}
json.dump(dict(frozen_rules=dict(flank_pos=FLANK_POS, ev_cov=EV_COV, ev_id=EV_ID, expand=EXPAND, seed=SEED),
               n_positives=len(pos), positive_rejects=dict(why), evidence_hits=n_ev,
               unknown_bp={k: sum(blen(masks[k][c]) for c in A) for k in masks}, nonpos_bp=nonpos_bp),
          open(OUT / 'masks_frozen.json', 'w'), indent=1)
print('masks frozen; now reading NEW annotation on A only', flush=True)

# ---- audit (§3.4): NEW annotation, A contigs only
new = read_gff(D / C['new_gff'], set(A))
old_cds = defaultdict(list); old_chains = defaultdict(set)
for g in old.values():
    for tx in g['txs']:
        old_cds[g['contig']] += list(tx[1]); old_chains[(g['contig'], g['strand'])].add(tx[1])
old_cds = {c: merge(v) for c, v in old_cds.items()}
pos_chains = {(c, st, ex) for c, st, _, ex in pos}
cls = defaultdict(list)
for gid, g in new.items():
    ex = rep(g)[1]; c = g['contig']
    if not any(covered(old_cds.get(c, []), a, b) for a, b in ex): cls['novel'].append((c, g['strand'], ex))
    elif any(tx[1] in old_chains[(c, g['strand'])] for tx in g['txs']): cls['identical'].append((c, g['strand'], ex))
    else: cls['revised'].append((c, g['strand'], ex))
oracle = {c: merge([(ex[0][0], ex[-1][1]) for cc, _, ex in cls['novel'] if cc == c]) for c in A}


def audit(mask):
    tot = inb = full = 0
    for c, st, ex in cls['novel']:
        m = mask.get(c, []); n = sum(covered(m, a, b) for a, b in ex); L = sum(b - a + 1 for a, b in ex)
        tot += L; inb += n; full += covered(m, ex[0][0], ex[-1][1]) == ex[-1][1] - ex[0][0] + 1
    nn = max(1, len(cls['novel']))
    return dict(novel_cds_bp_in_unknown=round(inb / max(1, tot), 4), novel_loci_fully_unknown=round(full / nn, 4), n_full=full,
                unknown_frac_of_nonpos=round(sum(blen(mask.get(c, [])) for c in A) / nonpos_bp, 4))


ident_pos = sum((c, st, ex) in pos_chains for c, st, ex in cls['identical'])
ev_hit = sum(any(covered(ev.get(c, []), a, b) for a, b in ex) for c, st, ex in cls['novel'])
res = dict(n_new_genes_A=len(new), n_novel=len(cls['novel']), n_revised=len(cls['revised']), n_identical=len(cls['identical']),
           positive_retention_of_identical=round(ident_pos / max(1, len(cls['identical'])), 4),
           novel_with_evidence=ev_hit, novel_without_evidence=len(cls['novel']) - ev_hit,
           variants={k: audit(masks[k]) for k in masks}, random={k: audit(rand[k]) for k in rand}, oracle=audit(oracle))
# §3.5 mechanical gate
cand = []
for k in ('V1', 'V2'):
    v, r = res['variants'][k], res['random'][k]
    ok = v['novel_loci_fully_unknown'] >= 2 * r['novel_loci_fully_unknown'] and v['unknown_frac_of_nonpos'] <= 0.30
    res['variants'][k]['passes_3_5'] = ok
    if ok: cand.append(k)
if len(cand) == 2: chosen = 'V1' if res['variants']['V1']['novel_loci_fully_unknown'] >= res['variants']['V2']['novel_loci_fully_unknown'] else 'V2'
elif cand: chosen = cand[0]
else: chosen = 'V3'
res['gate_3_5'] = dict(chosen=chosen, evidence_better_than_random=bool(cand), enough_novel=len(cls['novel']) >= 100,
                       note='V3 always enters gate 1 as the no-negative reference')
json.dump(res, open(OUT / 'audit.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
