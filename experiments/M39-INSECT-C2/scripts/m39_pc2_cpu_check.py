"""CPU check of the PlantCAD2 (pc2) backend on the M39 insect data (decision pc2-insect-transfer, 2026-11-11).  No GPU, no model weights; reads the
prepared genomes and the frozen M39 schedules (the pc2 arms reuse them unchanged: same draws, same labels, same head geometry).

Checks
  1 tokenizer : the adapter's lookup table (`pc2_adapter.token_ids`, ACGT -> 3..6, everything else -> [UNK] 2) equals the PlantCAD2 HF tokenizer on
                insect strings (incl. N and IUPAC codes, both strands).
  2 geometry  : every draw of the three schedules -> frozen context geometry -> pc2 8,192-bp window grid (stride 6,144); window lengths in (0, 8192] <= head
                width 8,196; base-owned segments tile the context exactly; number of encoder windows per draw (vs NT-v2 12,288-bp windows).
  3 sequence  : per species N fraction and GC of the drawn contexts; windows that are mostly / entirely N (PlantCAD2 sees only [UNK] there).
  4 regions   : the evaluation regions (S, D, validation crops) on the pc2 inference grid (window count, shortest window).
Output: JSON on stdout (and --out FILE).
usage: M39_BACKEND=pc2 ~/envs/plantcad2/bin/python m39_pc2_cpu_check.py [--out FILE] [--tok-samples 200]
"""
import argparse
import json
import os
import random
import sys
from collections import Counter, defaultdict

os.environ.setdefault('M39_BACKEND', 'pc2')
assert os.environ['M39_BACKEND'] == 'pc2', 'run with M39_BACKEND=pc2'
import numpy as np                                            # noqa: E402

import m39_train_common as T                                  # noqa: E402  (pc2: imports M32 pc2_adapter as T.A)
import m39_backend as BK                                      # noqa: E402

A, G, MF, C = T.A, T.G, T.MF, T.C
NT_WIDTH, NT_STRIDE = 12288, 6144


def tokenizer_check(seqs, n, rng):
    from transformers import AutoTokenizer
    tok = AutoTokenizer.from_pretrained(A.MODEL_DIR, trust_remote_code=True)
    strings = []
    for s in seqs:
        for _ in range(n // len(seqs) + 1):
            a = rng.randrange(0, max(1, len(s) - 2048))
            strings.append(s[a:a + 2048])
    strings += ['ACGTNRYKMSWBDHV' * 8, 'N' * 100, 'acgt' * 10]
    bad, n_unk_input = 0, 0
    for x in strings:
        ids = tok(x, add_special_tokens=False)['input_ids']
        lut = A.LUT[np.frombuffer(x.encode('ascii'), dtype=np.uint8)].tolist()
        n_unk_input += sum(1 for c in x if c not in 'ACGTacgt')
        if ids != lut:
            bad += 1
    return dict(strings=len(strings), mismatches=bad, non_acgt_chars=n_unk_input, unk_id=tok.unk_token_id,
                vocab={c: tok.convert_tokens_to_ids(c) for c in 'acgtACGTN'})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out')
    ap.add_argument('--tok-samples', type=int, default=200)
    args = ap.parse_args()
    rng = random.Random(0)
    tab = C.species_table()
    rows = defaultdict(list)
    for arm in T.S.ARMS:
        for r in T.read_schedule(T.SCHEDULES / f'schedule_{arm}.jsonl'):
            rows[r['species']].append((arm, r))
    out = dict(backend=T.BACKEND, pc_width=A.PC_WIDTH, pc_stride=A.PC_STRIDE, head_width=A.HEAD_W, init=str(T.INIT), species={}, problems=[])
    tok_seqs = []
    win_pc, win_nt, minlen = Counter(), Counter(), A.PC_WIDTH
    for key in sorted(rows):
        sp = tab[key]
        contigs = sorted({r['seqid'] for _, r in rows[key]})
        seqs = MF.read_allowed(sp['dir'] / 'genome.fa', {c: sp['lengths'][c] for c in contigs})
        tok_seqs.append(seqs[contigs[0]])
        n_tot = gc = acgt = 0
        n_win = mostly_n = all_n = 0
        for arm, r in rows[key]:
            chrom = seqs[r['seqid']]
            geo = G.context_geometry(r['start'], r['end'], len(chrom))
            a, b = geo['genomic_context']
            pg = BK.geometry(a, b)
            lens = [y - x for x, y in pg['local_windows']]
            if min(lens) <= 0 or max(lens) > A.PC_WIDTH or A.PC_WIDTH > A.HEAD_W:
                out['problems'].append(dict(draw=r['id'], arm=arm, lens=lens))
            cover = [(x, y) for _, x, y in pg['base_owned_segments']]          # (window index, start, end), local coordinates
            if cover[0][0] != 0 or cover[-1][1] != b - a or any(p[1] != q[0] or p[0] >= p[1] for p, q in zip(cover, cover[1:])):
                out['problems'].append(dict(draw=r['id'], arm=arm, segments='do not tile'))
            minlen = min(minlen, min(lens))
            win_pc[len(lens)] += 1
            win_nt[len(T.MF.starts(b - a, NT_WIDTH, NT_STRIDE))] += 1
            if arm == 'I10_T':                                   # sequence statistics once per species (T schedule)
                for w in BK.item_windows(chrom, pg, '+')['gen']:
                    k = sum(w.count(c) for c in 'ACGT')
                    g = w.count('G') + w.count('C')
                    n_tot += len(w); acgt += k; gc += g; n_win += 1
                    if k < 0.5 * len(w):
                        mostly_n += 1
                    if k == 0:
                        all_n += 1
        out['species'][key] = dict(draws=len(rows[key]), windows_T=n_win, N_frac=round(1 - acgt / n_tot, 5), GC=round(gc / max(acgt, 1), 4),
                                   windows_mostly_N=mostly_n, windows_all_N=all_n)
        print(key, json.dumps(out['species'][key]), file=sys.stderr, flush=True)
    out['windows_per_draw_pc2'] = dict(sorted(win_pc.items()))
    out['windows_per_draw_ntv2'] = dict(sorted(win_nt.items()))
    out['shortest_training_window'] = minlen
    reg = {}
    for name, regions in (('S', T.select_regions()), ('D', T.d_regions()), ('val', [T.val_contig(k) for k in T.new_keys()])):
        for c in regions:
            sp = tab[[k for k in tab if c in tab[k]['val'] + tab[k]['heldout_d']][0]]
            L = int(sp['lengths'][c])
            st = MF.starts(L, A.PC_WIDTH, A.PC_STRIDE)
            reg[c] = dict(set=name, length=L, pc2_windows=len(st), shortest=min(min(s + A.PC_WIDTH, L) - s for s in st))
    out['regions'] = reg
    out['tokenizer'] = tokenizer_check(tok_seqs, args.tok_samples, rng)
    if out['tokenizer']['mismatches']:
        out['problems'].append('tokenizer mismatch')
    out['ok'] = not out['problems']
    s = json.dumps(out, indent=1)
    if args.out:
        open(args.out, 'w').write(s)
    print(s)


if __name__ == '__main__':
    main()
