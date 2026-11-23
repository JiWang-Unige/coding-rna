"""M42 gate 1 data preparation (protocol 35 §4; CPU, pure Python, no model).  Reads ONLY the frozen gate-0 outputs and the genome:
  <species_dir>/gate0/{split.json, positives.bed, mask_V1.bed, mask_V3.bed}, <species_dir>/<genome>
and writes the files the GPU side reads (no annotation of the target species beyond the frozen OLD-release positives; the NEW release never leaves the Mac):
  <out>/m1_screen/<species>/{genome.fa, lengths.json}                  sequences of the 12 chromosomes (A and E; E is only used for inference)
  <out>/m42/<species>/positives.json                                   seed rows (seqid, strand, CDS [[a0, b0, phase]] 0-based half-open incl. stop, id)
  <out>/m42/<species>/background_<V>.json                              {contig: [[a, b], ...]} = A contig minus mask_<V> (0-based half-open), V in V1, V3
  <out>/m42/<species>/split.json, eval_regions.txt, prep.json          (copy of the frozen split; comma list of E contigs; md5 of everything)
usage: m42_prep.py <species_dir> <species_key> <out_root>
"""
import hashlib
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

D, SP, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
G0 = D / 'gate0'
cfg = json.load(open(D / 'm42_config.json'))
split = json.load(open(G0 / 'split.json'))
A, E, L = split['A'], split['E'], split['lengths']


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def phases(ex, strand):
    order = ex if strand == '+' else ex[::-1]
    out, off = {}, 0
    for a, b in order:
        out[(a, b)] = (3 - off % 3) % 3
        off += b - a
    return [[a, b, str(out[(a, b)])] for a, b in ex]


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def complement(iv, n):
    out, cur = [], 0
    for a, b in merge(iv):
        if a > cur:
            out.append([cur, a])
        cur = max(cur, b)
    if cur < n:
        out.append([cur, n])
    return out


m1 = OUT / 'm1_screen' / SP
m42 = OUT / 'm42' / SP
m1.mkdir(parents=True, exist_ok=True)
m42.mkdir(parents=True, exist_ok=True)

# genome (12 chromosomes, exactly the gate-0 sequences) + lengths
src = D / cfg['genome']
seen = {}
with open(src) as f, open(m1 / 'genome.fa', 'w') as o:
    cur = None
    for line in f:
        if line.startswith('>'):
            cur = line[1:].split()[0]
            seen[cur] = 0
            o.write(f'>{cur}\n')
        else:
            seen[cur] += len(line.strip())
            o.write(line)
assert seen == {c: L[c] for c in seen} and set(seen) == set(A) | set(E), 'genome does not match split.json lengths'
json.dump({c: L[c] for c in split['order']}, open(m1 / 'lengths.json', 'w'), indent=1)

# positives: one bed line per CDS exon (0-based half-open, stop codon included), name = gene id
ex = defaultdict(list)
strand, contig = {}, {}
for line in open(G0 / 'positives.bed'):
    c, a, b, gid, _, st = line.rstrip('\n').split('\t')
    assert c in A, c
    ex[gid].append((int(a), int(b)))
    strand[gid], contig[gid] = st, c
rows = []
for gid in ex:
    e = sorted(ex[gid])
    rows.append(dict(seqid=contig[gid], strand=strand[gid], CDS=phases(e, strand[gid]), id=gid, gene_id=gid, source='old_release_positive'))
rows.sort(key=lambda r: (r['seqid'], r['CDS'][0][0], r['strand'], r['id']))
json.dump(rows, open(m42 / 'positives.json', 'w'))

# known background = A contig minus the frozen unknown mask
bg_stats = {}
for v in ('V1', 'V3'):
    mask = defaultdict(list)
    for line in open(G0 / f'mask_{v}.bed'):
        c, a, b, _ = line.rstrip('\n').split('\t')
        mask[c].append((int(a), int(b)))
    bg = {c: complement(mask.get(c, []), L[c]) for c in A}
    json.dump(bg, open(m42 / f'background_{v}.json', 'w'))
    bg_stats[v] = dict(unknown_bp=sum(b - a for c in A for a, b in merge(mask.get(c, []))), background_bp=sum(b - a for c in A for a, b in bg[c]))

shutil.copy(G0 / 'split.json', m42 / 'split.json')
(m42 / 'eval_regions.txt').write_text(','.join(E) + '\n')
files = {str(p.relative_to(OUT)): md5(p) for p in sorted(list(m1.iterdir()) + list(m42.iterdir())) if p.name != 'prep.json'}
info = dict(species=SP, A=A, E=E, A_bp=sum(L[c] for c in A), E_bp=sum(L[c] for c in E), positives=len(rows), positive_exons=sum(len(r['CDS']) for r in rows),
            mask=bg_stats, gate0_inputs={p: md5(G0 / p) for p in ('split.json', 'positives.bed', 'mask_V1.bed', 'mask_V3.bed', 'masks_frozen.json')}, files=files)
json.dump(info, open(m42 / 'prep.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in info.items() if k != 'files'}, indent=1))
