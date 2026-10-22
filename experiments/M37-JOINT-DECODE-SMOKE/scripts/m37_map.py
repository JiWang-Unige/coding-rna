"""Base-level map between assembly A (minimap2 target, NC_089044.1) and assembly B (query, CM000135.2) from the PAF `cs` tags (protocol §3).
Kept alignments: primary (tp:A:P), identity >= 0.90, block length >= 1 kb.  Identity = 1 - de (minimap2 `de` tag: gap-compressed per-base
divergence, the standard minimap2 identity; decided 2026-10-21 23:06 before any §4 result, see protocol §8): the raw matches / block-length ratio
(PAF columns 10/11) counts every gap base and falls to 0.75-0.90 for the long asm10 chains between indica and japonica although their aligned
bases are 98.5 % identical; both ratios are stored per block and the coverage under the raw ratio is reported in the stats as well.
Where kept blocks overlap, a base belongs to the block with more matching bases (blocks are applied in descending order of matches; an
already-mapped base is never reassigned), so the map is a partial bijection: b2a[a2b[a]] == a wherever a2b[a] >= 0.
Substituted bases (cs `*`) are mapped (homologous SNP positions); bases inside insertions / deletions are unmapped.
Arrays (npz):  a2b int64[LA] (-1 unmapped)  a2b_strand int8[LA] (+1 same strand, -1 inverted, 0 unmapped)  b2a / b2a_strand likewise
               a_state uint8[LA] (0 unmapped, 1 match, 2 mismatch, 3 A base absent from B inside a kept block)  a_ins uint8[LA] (# B-insertion
               events anchored at the preceding A base)   blocks: structured table of every PAF record with the kept flag.
Usage: m37_map.py <paf> <A.fa> <B.fa> <out.npz> [--min-identity 0.9 --min-block 1000]
"""
import argparse, json, re, sys
from pathlib import Path
import numpy as np

CS_RE = re.compile(r'([:*+\-~])([A-Za-z0-9]+)')
COMP = str.maketrans('ACGTN', 'TGCAN')


def read_fa_one(path):
    name, buf = None, []
    for line in open(path):
        if line[0] == '>':
            if name is not None:
                break
            name = line[1:].split()[0]
        else:
            buf.append(line.strip())
    return name, ''.join(buf).upper()


def parse_paf(path):
    rows = []
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        tags = {t[:2]: t[5:] for t in f[12:]}
        rows.append(dict(qname=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], tname=f[5], tlen=int(f[6]), ts=int(f[7]), te=int(f[8]),
                         matches=int(f[9]), blen=int(f[10]), mapq=int(f[11]), tp=tags.get('tp', ''), cs=tags.get('cs', ''),
                         de=float(tags['de']) if 'de' in tags else None))
    return rows


def walk_cs(cs, ts, qs, qe, strand):
    """yield (op, t_start, q_start, n) with per-op coordinates: ':' and '*' -> n aligned pairs t_start+i <-> q_start + dq*i;
    '-' (deletion from query: target-only bases) -> n target bases; '+' (insertion in query) -> n query bases anchored at target t_start-1."""
    t = ts
    q, dq = (qs, 1) if strand == '+' else (qe - 1, -1)
    for op, arg in CS_RE.findall(cs):
        if op == ':':
            n = int(arg)
            yield op, t, q, n
            t += n; q += dq * n
        elif op == '*':
            yield op, t, q, 1
            t += 1; q += dq
        elif op == '-':                       # bases present in target (A) only
            n = len(arg)
            yield op, t, q, n
            t += n
        elif op == '+':                       # bases present in query (B) only
            n = len(arg)
            yield op, t, q, n
            q += dq * n
        else:
            raise ValueError(f'unsupported cs op {op} (introns are not expected in asm alignments)')
    return t, q


def build_map(paf_rows, LA, LB, min_identity=0.9, min_block=1000, A=None, B=None):
    a2b = np.full(LA, -1, np.int64); b2a = np.full(LB, -1, np.int64)
    a2b_s = np.zeros(LA, np.int8); b2a_s = np.zeros(LB, np.int8)
    a_state = np.zeros(LA, np.uint8); a_ins = np.zeros(LA, np.uint8)
    blocks = []
    kept = []
    for i, r in enumerate(paf_rows):
        raw = r['matches'] / r['blen'] if r['blen'] else 0.
        ident = 1. - r['de'] if r['de'] is not None else raw
        ok = r['tp'] == 'P' and ident >= min_identity and r['blen'] >= min_block and bool(r['cs'])
        blocks.append((i, r['qs'], r['qe'], r['strand'], r['ts'], r['te'], r['matches'], r['blen'], ident, raw, ok))
        if ok:
            kept.append((r['matches'], i, r))
    kept.sort(key=lambda x: (-x[0], x[1]))
    n_check = n_bad = 0
    for _, i, r in kept:
        sgn = 1 if r['strand'] == '+' else -1
        for op, t, q, n in walk_cs(r['cs'], r['ts'], r['qs'], r['qe'], r['strand']):
            if op in ':*':
                ta = np.arange(t, t + n); qa = q + sgn * np.arange(n)
                free = (a2b[ta] < 0) & (b2a[qa] < 0)
                ta, qa = ta[free], qa[free]
                a2b[ta] = qa; b2a[qa] = ta; a2b_s[ta] = sgn; b2a_s[qa] = sgn
                a_state[ta] = 1 if op == ':' else 2
                if A is not None and op == ':' and n:        # verify on the real sequences: matched runs must be identical (complement if inverted)
                    sa = A[t:t + n]
                    sb = B[q:q + n] if sgn > 0 else B[q - n + 1:q + 1][::-1].translate(COMP)
                    n_check += n; n_bad += sum(x != y for x, y in zip(sa, sb)) if sa != sb else 0
            elif op == '-':
                ta = np.arange(t, t + n)
                sel = a2b[ta] < 0
                a_state[ta[sel]] = np.where(a_state[ta[sel]] == 0, 3, a_state[ta[sel]])
            else:
                if 0 < t <= LA:
                    a_ins[t - 1] = min(255, int(a_ins[t - 1]) + 1)
    assert np.all(b2a[a2b[a2b >= 0]] == np.nonzero(a2b >= 0)[0]), 'map is not a partial bijection'
    r_ok = lambda b: paf_rows[b[0]]['tp'] == 'P' and b[9] >= min_identity and b[7] >= min_block
    btab = np.array(blocks, dtype=[('i', 'i8'), ('qs', 'i8'), ('qe', 'i8'), ('strand', 'U1'), ('ts', 'i8'), ('te', 'i8'), ('matches', 'i8'), ('blen', 'i8'), ('identity', 'f8'), ('identity_raw_matches_over_blen', 'f8'), ('kept', '?')])
    stats = dict(LA=LA, LB=LB, paf_records=len(paf_rows), kept_blocks=int(btab['kept'].sum()), kept_inverted=int(sum(1 for b in blocks if b[10] and b[3] == '-')),
                 kept_target_span=int(sum(b[5] - b[4] for b in blocks if b[10])),
                 blocks_kept_under_raw_ratio=int(sum(1 for b in blocks if r_ok(b))), target_span_kept_under_raw_ratio=int(sum(b[5] - b[4] for b in blocks if r_ok(b))),
                 identity_definition='1 - de (gap-compressed)', min_identity=min_identity, min_block=min_block,
                 A_mapped=int((a2b >= 0).sum()), B_mapped=int((b2a >= 0).sum()), A_mapped_frac=float((a2b >= 0).mean()), B_mapped_frac=float((b2a >= 0).mean()),
                 A_inverted=int((a2b_s < 0).sum()), A_match=int((a_state == 1).sum()), A_mismatch=int((a_state == 2).sum()), A_deleted_in_B=int((a_state == 3).sum()),
                 B_insertion_events=int(a_ins.sum()), sequence_check_bases=n_check, sequence_check_mismatches=n_bad)
    return dict(a2b=a2b, a2b_strand=a2b_s, b2a=b2a, b2a_strand=b2a_s, a_state=a_state, a_ins=a_ins, blocks=btab), stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('paf'); ap.add_argument('A_fa'); ap.add_argument('B_fa'); ap.add_argument('out')
    ap.add_argument('--min-identity', type=float, default=0.9); ap.add_argument('--min-block', type=int, default=1000)
    a = ap.parse_args()
    an, A = read_fa_one(a.A_fa); bn, B = read_fa_one(a.B_fa)
    rows = parse_paf(a.paf)
    assert all(r['tname'] == an and r['qname'] == bn for r in rows), 'PAF target must be A, query must be B'
    m, st = build_map(rows, len(A), len(B), a.min_identity, a.min_block, A, B)
    st['B_nonN'] = B.count('A') + B.count('C') + B.count('G') + B.count('T')
    st['B_mapped_frac_of_nonN'] = st['B_mapped'] / st['B_nonN']
    assert st['sequence_check_mismatches'] == 0, st
    np.savez(a.out, **m)
    json.dump(st, open(Path(a.out).with_suffix('.json'), 'w'), indent=1)
    print(json.dumps(st))


if __name__ == '__main__':
    main()
