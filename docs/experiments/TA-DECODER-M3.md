# TA-DECODER-M3 — Track A structured-decoder focused batch

- approach_family: structured-decoder (on per-base backbone)
- date: 2026-06-11 | track: A-screen (M3) | profile: screen (non-claim)
- parent: M1-SAMEBUDGET-SCREEN-ANCHOR (screen_anchor 0.5576)

## Hypothesis
A structured decoder on the fixed tiberius_like backbone beats the same-budget per-base anchor on BOTH base-F1 and gene-level coherence (per-base models fragment: ratio 1.8-153).

## Architecture / mechanism
Backbone FIXED (tiberius_like CNN+biLSTM emissions). Decoder varies (major_axis=decoder):
- CONSTR: constrained-Viterbi post-processing (drop sub-min CDS, fill small gaps) — WON.
- CRF: linear-chain CRF (learned transitions + Viterbi) — correct (tests 5/5) but dropped, W=2048 too slow.
- semi-CRF: semi-Markov segment DP — dropped, pure-python intractable.

## Data / protocol
Frozen same-budget screen protocol (= anchor refs): yeast+fly, chrom-split, window 2048, sample 0.3, 8 epochs, patience 3, 3 seeds, class-weighted CE; eval --span-mode cds.

## Result
CONSTR seed-mean base-w gene_body_F1_unconstrained=0.5791 (>gate 0.5676, >anchor 0.5576); gene_count_ratio 2.74->1.12 (<1.25 claim guardrail). primary_progress_gate MET; R5 satisfied.

## Findings
Structured decoding validated (coherence + F1). But CONSTR is post-processing, not learned structure; CRF/semi-CRF need vectorization (tractability) to test the core learned-structure bet. High seed variance (0.5319/0.5779/0.6275).

## Lineage / next
Pivot: 1-1 reviewer split; user gate on promote-CONSTR-Track-B-now vs vectorize-learned-decoders-first. Both schedule vectorized CRF/semi-CRF batch.
