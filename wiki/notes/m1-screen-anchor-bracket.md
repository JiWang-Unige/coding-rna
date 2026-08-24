# Note (tried-once): M1 same-budget screen_anchor bracket (CDS gene-body F1)

- slug: `m1-screen-anchor-bracket` · added: 2026-06-10
- refs: tiberius-2024 helixer-2025 annevo-2026

## What I ran
Trained 2 random-init reference archs (Tiberius-like CNN+biLSTM, Helixer-like dilated-conv) under one frozen small-sample screen protocol (yeast+fly, chrom-split, window 2048, sample 0.3, 8 epochs, 3 seeds, class-weighted CE), eval --span-mode cds; plus an ORF FLOOR and the pretrained-inference ceiling.

## Quick result
base-weighted CDS gene_body_F1_unconstrained: FLOOR(ORF)=0.3735 < screen_anchor=0.5579 (tiberius seed-mean 0.5576 / helixer 0.5579) < pretrained_ceiling=0.9213.

## Takeaway
The same-budget bar (0.56) is FAR below the pretrained-inference ceiling (0.92). Track A from-scratch candidates must strictly exceed 0.5579 (same protocol), NOT 0.92 — using the pretrained value would be the unfair small-sample-vs-large-sample-SOTA comparison.

## Next direction
Run Track A architecture portfolio vs 0.5579 (primary track foundation_probe->semi-CRF). Re-derive on frozen typical-intergenic species before heavy reliance (yeast/fly are gene-dense outliers).
