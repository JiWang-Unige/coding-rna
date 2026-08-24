"""FLOOR baseline (no training): a naive 3-forward-frame ORF finder. Marks ATG..in-frame-stop
ORFs of length >= min_len as CDS. Deliberately simple/weak -> a non-trivial LOWER bracket so we
can assert floor < screen_anchor < pretrained_ceiling. Forward strand only (underpredicts vs a
both-strand reference, which is fine for a floor). Pure stdlib + numpy.

Usage: python -m src.screen_anchor.floor_baseline --fasta <genome.fa> --out <cds.gff> [--min-len 300]
"""
import argparse

import numpy as np

from .data import read_fasta, CLASS_CDS
from .gff_io import labels_to_cds_gff

STOPS = {"TAA", "TAG", "TGA"}


def orf_labels(seq, min_len=300):
    """Return int8 array: 1 where inside a >=min_len forward-frame ORF, else 0."""
    L = len(seq)
    lab = np.zeros(L, dtype=np.int8)
    su = seq.upper()
    for frame in range(3):
        i = frame
        while i + 3 <= L:
            if su[i:i + 3] == "ATG":
                j = i
                while j + 3 <= L:
                    if su[j:j + 3] in STOPS:
                        end = j + 3
                        if end - i >= min_len:
                            lab[i:end] = CLASS_CDS
                        i = end
                        break
                    j += 3
                else:
                    break
                continue
            i += 3
    return lab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=300)
    args = ap.parse_args()
    seqs = read_fasta(args.fasta)
    pred = {sid: orf_labels(s, args.min_len) for sid, s in seqs.items()}
    n = labels_to_cds_gff(pred, args.out, source="floor_orf")
    print(f"FLOOR ORF baseline: {len(seqs)} seqs -> {n} predicted gene-bodies -> {args.out}")


if __name__ == "__main__":
    main()
