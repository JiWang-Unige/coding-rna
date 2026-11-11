"""Repeat (soft-mask) index of the TRAINING contigs of the 10 M34 species, for the MLM repeat down-weighting (CPU only).

The M34 prepared genomes of the 8 new species were upper-cased by m34_prep.py, so the soft-masking is read back from the md5-verified NCBI downloads
(`data/m34_train/raw/<acc>_<asm>/<stem>_genomic.fna.gz`, NCBI RefSeq soft-masking); the legacy genomes (`data/m1_screen/<sp>/genome.fa`) are still soft-masked.
For every training contig the upper-cased source sequence is checked to be identical (md5) to the sequence the trainer reads.
Output: ~/coding-rna/outputs/M40-SEMISUP-MLM/repeats/<species>.npz  ({contig: int64 [2k] flattened [start, end) lowercase runs}) + summary.json
usage: m40_repeats.py build [species ...]
"""
import gzip
import hashlib
import json
import sys
import time

import numpy as np

import m40_common as W

C, T = W.C, W.T


def runs_lower(seq):
    """[start, end) runs of lowercase letters -> flattened int64 array."""
    b = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    low = (b >= 97) & (b <= 122)
    d = np.diff(np.concatenate([[0], low.astype(np.int8), [0]]))
    s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    return np.stack([s, e], 1).reshape(-1).astype(np.int64)


def read_fasta(path, wanted):
    op = gzip.open if str(path).endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with op(path, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if name in wanted:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = ''.join(buf)
    return seqs


def source(key):
    tab = C.species_table()
    if key in C.LEGACY:
        return tab[key]['dir'] / 'genome.fa'
    row = C.new_species()[key]
    d = C.raw_dir(row)
    return d / f"{row['accession']}_{row['assembly_name']}_genomic.fna.gz"


def build(key):
    t0 = time.time()
    sp = C.species_table()[key]
    contigs = list(sp['train'])
    soft = read_fasta(source(key), set(contigs))
    assert set(soft) == set(contigs), set(contigs) - set(soft)
    ref = T.MF.read_allowed(sp['dir'] / 'genome.fa', {c: sp['lengths'][c] for c in contigs})      # what the trainer reads (upper case)
    out, stats = {}, {}
    for c in contigs:
        up = soft[c].upper()
        assert len(up) == len(ref[c]) and hashlib.md5(up.encode()).hexdigest() == hashlib.md5(ref[c].encode()).hexdigest(), f'{key} {c}: source differs from trainer sequence'
        r = runs_lower(soft[c])
        out[c] = r
        stats[c] = dict(length=len(up), repeat_bp=int((r[1::2] - r[0::2]).sum()), runs=len(r) // 2)
    W.REPEATS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(W.REPEATS / f'{key}.npz', **out)
    tot = sum(v['length'] for v in stats.values())
    rep = sum(v['repeat_bp'] for v in stats.values())
    return dict(species=key, source=str(source(key)), contigs=stats, repeat_fraction=rep / tot, seconds=round(time.time() - t0, 1))


class RepeatIndex:
    """Per-species soft-mask runs; `mask(species, contig, a, b)` -> bool [b-a] (True = repeat)."""

    def __init__(self, species):
        self.runs = {}
        for key in species:
            z = np.load(W.REPEATS / f'{key}.npz')
            self.runs[key] = {c: z[c] for c in z.files}

    def mask(self, species, contig, a, b):
        r = self.runs[species][contig]
        s, e = r[0::2], r[1::2]
        out = np.zeros(b - a, dtype=bool)
        i = np.searchsorted(e, a, side='right')
        j = np.searchsorted(s, b, side='left')
        for x, y in zip(s[i:j], e[i:j]):
            out[max(x, a) - a:min(y, b) - a] = True
        return out


def main():
    assert sys.argv[1] == 'build'
    keys = sys.argv[2:] or list(C.species_table())
    p = W.REPEATS / 'summary.json'
    summ = json.load(open(p)) if p.exists() else {}
    for key in keys:
        summ[key] = build(key)
        print(json.dumps({k: summ[key][k] for k in ('species', 'repeat_fraction', 'seconds')}), flush=True)
        json.dump(summ, open(p, 'w'), indent=1)


if __name__ == '__main__':
    main()
