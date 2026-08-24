# TA-DECODER-VEC-M3 — vectorized LEARNED structured decoder (Track A)

- approach_family: structured-decoder (learned, on per-base backbone)
- date: 2026-06-11 | track: A-screen (M3) | profile: screen (non-claim)
- parent: TA-DECODER-M3 (CONSTR post-processing won; learned decoders untested for tractability)

## Hypothesis
A LEARNED structured decoder (CRF), made tractable via vectorization, beats both the per-base anchor AND the cheap CONSTR post-processing — i.e. learned structure is worth more than post-processing.

## Architecture / mechanism
Backbone FIXED (tiberius_like). CRF-vec = vectorized linear-chain CRF: log-space ASSOCIATIVE-SCAN partition (O(log W)), per-token NLL normalization (+aux class-weighted CE), batched Viterbi predict. semi-CRF deferred (vectorize segment DP).

## Result (same frozen screen protocol; base-w gene_body_F1_unconstrained, CDS span)
LADDER: anchor 0.5576 (ratio 2.74) < CONSTR 0.5791 (1.12) < **CRF-vec 0.6186 (0.88)**. CRF-vec per-seed 0.6605/0.6153/0.5799. Beats gate 0.5676, anchor, AND CONSTR (mean + 2/3 seeds). primary_progress_gate MET; R5 satisfied.

## Findings
LEARNED structured decoding > post-processing > per-base — core architecture bet validated. Vectorization (parallel-scan + per-token-NLL + batched-predict) turned an intractable CRF (epoch1 >28min) into ~20min/run. CAVEAT: high seed variance (s2 loses to CONSTR); robustness is a Track B question.

## Lineage / next
Pivot (2/3 consensus): promote CRF-vec to Track B (seeds≥5-8 + CI + scale-data scalability test; keep CONSTR baseline). LAUNCH pending user. Next branches: vectorize semi-CRF; foundation-probe path.
