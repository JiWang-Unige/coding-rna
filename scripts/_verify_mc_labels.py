"""M8-CK1 self-review: verify build_labels_multiclass + collapse_mc_to_3class.
(1) MC labels: all 8 classes present with sane proportions (phase 0/1/2 ~balanced; donor/acceptor
    sparse but nonzero; intron/UTR present for UTR-rich species).
(2) CORRECTNESS: collapse_mc_to_3class(MC) must match the original 3-class build_labels (CDS regions
    identical; genebody = intron∪UTR∪splice = span−CDS). Report mismatch fraction (should be ~0).
Run: python scripts/_verify_mc_labels.py
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from src.screen_anchor import data as D

MC_NAMES = ["intergenic", "CDS_p0", "CDS_p1", "CDS_p2", "intron", "UTR", "donor", "acceptor"]
SPECIES = ["arabidopsis_thaliana", "oryza_sativa", "saccharomyces_cerevisiae"]

for sp in SPECIES:
    g = f"data/m1_screen/{sp}/genome.fa"; gff = f"data/m1_screen/{sp}/reference.gff3"
    try:
        seqs = D.read_fasta(g)
    except FileNotFoundError:
        print(f"{sp}: genome.fa MISSING, skip"); continue
    seqlen = {k: len(v) for k, v in seqs.items()}
    mc = D.build_labels_multiclass(gff, seqlen)
    tri = D.build_labels(gff, seqlen)
    mc_all = np.concatenate([mc[s] for s in seqs])
    tri_all = np.concatenate([tri[s] for s in seqs])
    coll = D.collapse_mc_to_3class(mc_all)
    tot = len(mc_all)
    cnt = np.bincount(mc_all.astype(int), minlength=8)
    print(f"\n=== {sp} ({tot/1e6:.1f}Mb) MC class fractions ===")
    for i, nm in enumerate(MC_NAMES):
        print(f"  {i} {nm:10s} {cnt[i]/tot:.5f} ({cnt[i]})")
    # correctness: collapse(MC) vs original 3-class
    mism = int((coll != tri_all).sum())
    print(f"  collapse(MC) vs 3-class build_labels: mismatch {mism}/{tot} = {mism/tot:.6f}")
    # per-3class agreement
    for c, nm in [(0, "intergenic"), (1, "CDS"), (2, "genebody")]:
        tri_c = (tri_all == c); coll_c = (coll == c)
        inter = int((tri_c & coll_c).sum()); union = int((tri_c | coll_c).sum())
        print(f"    class {c}({nm}) IoU {inter/union if union else 1.0:.5f}")
