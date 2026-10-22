"""Mapping round trip on synthetic assemblies through the real minimap2 (-cx asm10 --cs) and m37_map.build_map.
B is built from A with known edits: SNPs + short indels (segment 1), an inverted segment (2), a deleted A segment and a B-only insertion (3).
Checks: partial bijection, strand signs, exact agreement with the known truth on >= 99.5 % of mapped bases, inverted / deleted / inserted regions."""
import os, random, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import m37_map as MM
MINIMAP = os.path.expanduser('~/envs/m37/bin/minimap2')
COMP = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(COMP)[::-1]


def build():
    rng = random.Random(0)
    LA = 90_000
    A = ''.join(rng.choice('ACGT') for _ in range(LA))
    truth = []                                             # per B base: (a, sign) or (-1, 0)
    Bp = []
    # segment 1: A[0:30000] with SNPs every ~150 bp, 1-3 bp indels every ~2 kb
    i = 0
    while i < 30_000:
        if rng.random() < 1 / 150:
            Bp.append(rng.choice([c for c in 'ACGT' if c != A[i]])); truth.append((i, 1)); i += 1
        elif rng.random() < 1 / 4000:
            n = rng.randint(1, 3); i += n                                  # deletion of A bases
        elif rng.random() < 1 / 4000:
            for _ in range(rng.randint(1, 3)):
                Bp.append(rng.choice('ACGT')); truth.append((-1, 0))     # insertion of B-only bases
        else:
            Bp.append(A[i]); truth.append((i, 1)); i += 1
    # segment 2: inverted A[30000:45000]
    seg = A[30_000:45_000]
    Bp.append(rc(seg)); truth += [(45_000 - 1 - j, -1) for j in range(len(seg))]
    # segment 3: A[45000:50000] deleted; B-only 800-bp insertion; then A[50000:90000]
    ins = ''.join(rng.choice('ACGT') for _ in range(800))
    Bp.append(ins); truth += [(-1, 0)] * 800
    Bp.append(A[50_000:]); truth += [(j, 1) for j in range(50_000, LA)]
    B = ''.join(Bp)
    assert len(B) == len(truth)
    return A, B, truth


def test_roundtrip():
    if not os.path.exists(MINIMAP):
        print('SKIP: minimap2 not installed'); return
    A, B, truth = build()
    with tempfile.TemporaryDirectory() as d:
        fa = Path(d) / 'A.fa'; fb = Path(d) / 'B.fa'; paf = Path(d) / 'ab.paf'
        fa.write_text('>A\n' + A + '\n'); fb.write_text('>B\n' + B + '\n')
        with open(paf, 'w') as f:
            subprocess.run([MINIMAP, '-cx', 'asm10', '--cs', '-t', '1', str(fa), str(fb)], stdout=f, stderr=subprocess.DEVNULL, check=True)
        rows = MM.parse_paf(str(paf))
    m, st = MM.build_map(rows, len(A), len(B), 0.9, 1000, A, B)
    assert st['sequence_check_mismatches'] == 0
    a2b, b2a, b2a_s = m['a2b'], m['b2a'], m['b2a_strand']
    ok = np.nonzero(a2b >= 0)[0]
    assert np.all(b2a[a2b[ok]] == ok)                                   # partial bijection
    tr_a = np.array([t[0] for t in truth]); tr_s = np.array([t[1] for t in truth])
    mapped = b2a >= 0
    agree = (b2a[mapped] == tr_a[mapped]) & (b2a_s[mapped] == tr_s[mapped])
    assert agree.mean() >= 0.995, agree.mean()
    assert mapped.mean() > 0.97, mapped.mean()
    inv = (tr_s == -1)
    assert (b2a_s[inv & mapped] == -1).all() and (b2a_s[(tr_s == 1) & mapped] == 1).all()
    assert (a2b[45_000:50_000] >= 0).mean() < 0.01                       # deleted A segment unmapped
    ins_pos = np.nonzero(tr_s == 0)[0]
    assert (b2a[ins_pos] >= 0).mean() < 0.05                              # B-only bases unmapped
    assert st['kept_inverted'] >= 1 and st['A_mismatch'] > 100
    print('test_roundtrip OK', {k: st[k] for k in ('kept_blocks', 'kept_inverted', 'A_mapped_frac', 'A_mismatch', 'B_insertion_events')}, 'agreement', float(agree.mean()))


if __name__ == '__main__':
    test_roundtrip()
