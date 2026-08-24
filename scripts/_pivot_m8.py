"""tri-review (docs/07) + pivot (docs/08) for TB-GBF1-MULTICLASS-M8. python3 scripts/_pivot_m8.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D07 = """

## Tri-Review: TB-GBF1-MULTICLASS-M8
### Review mode
- independent_parallel_cli, one identical full-scope prompt (/tmp/tri_review_TB-GBF1-MULTICLASS-M8/prompt_full_scope.md)
- A Claude CLI · success | B Codex CLI (gpt-5.5) · success | C Antigravity (agy) · FAILED (timeout 420s)
- **Quorum: 2/3 DEGRADED_REVIEW; confidence ceiling Medium.** (pivot = route-direction to a BOUNDED screen-profile next step, not a claim/abandon; 2 independent reviewers sufficient; the full >24h unfreeze scale-up remains separately human-gated.)

### Reviewer A · Claude — unfreeze-finetune-backbone (BOUNDED screen first), Confidence Medium
- "frozen features cap gbF1" is currently INFERRED not MEASURED -> first step MUST be a bounded screen-profile partial/staged unfreeze (top N layers + low LR), doubling as sanity + route entry; NOT a direct >24h full run.
- Lever rank: (1) staged unfreeze — ANNEVO ceiling 0.8976 is END-TO-END trained vs frozen-head 0.74; frozen-vs-end-to-end is the natural explanation of the 0.16 gap; partial unfreeze directly tests this. (2) backbone-only self-train = clean control but if frozen CAPACITY is the bottleneck it stays limited.
- M8 negative + 3c clean-positive both sound.

### Reviewer B · Codex (gpt-5.5) — unfreeze-finetune-backbone (staged preflight), Confidence Medium-High
- Lever rank: (1) staged unfreeze; (2) backbone-only domain-adaptation (masked/self-sup or pseudo-label, then 3c head) — good if labels scarce / overfit worry; (3) different foundation model (higher cost, only if unfreeze fails or SegmentNT channels lack plant resolution); (4) evidence/multi-task (auxiliary splice/phase/ORF OK later; RNA/protein evidence breaks ab-initio purity — not the mainline); (5) accept-frozen-ceiling NOT recommended (gap too big -> route can't reach north star).
- Biggest risk: unfreeze comparability/leakage — clean species/chrom split, no test labels in early-stop/decode-tuning, stay raw-DNA ab-initio, same 3-class collapse ruler, no test-truth gene_count calibration.
- Next: staged-unfreeze Track-B PREFLIGHT — unfreeze last N SegmentNT layers + existing 3c FP-aware constrained head, SHORT budget on clean plant split; success = gbF1 directionally > frozen 3c 0.7392, spec not collapsing, gene_count sane (avoid mc-style under-call). If no directional gbF1 gain -> backbone-only domain-adapt or different foundation model.

### Cross-reviewer agreement (2/2)
- M8 multi-class bet REFUTED (mc gbF1 <= 3c, worse coherence); multi-class NOT scaled. 3c-candidate clean Pareto-over-anchor on plants is the validated leakage-free lead.
- Next axis = STAGED UNFREEZE / fine-tune SegmentNT, but the FIRST step is a BOUNDED screen-profile preflight (measure "frozen caps gbF1" before any >24h spend), success = gbF1 directionally > 0.7392 while spec held + gene_count sane.
- accept-frozen-ceiling rejected; evidence/RNA breaks ab-initio purity (not mainline).

### Aggregated recommendation
- [x] unfreeze-finetune-backbone — via a BOUNDED screen-profile staged-unfreeze PREFLIGHT first (NOT direct >24h). multi-class dropped. 3c clean-positive = current honest lead.

### Confidence
Medium (2/3 DEGRADED; both independent reviewers align; main uncertainty = mc under-prediction confound + whether unfreeze sacrifices specificity).

### Raw outputs
- /tmp/tri_review_TB-GBF1-MULTICLASS-M8/output_a_claude.md
- /tmp/tri_review_TB-GBF1-MULTICLASS-M8/output_b_codex.md
- /tmp/tri_review_TB-GBF1-MULTICLASS-M8/output_c_antigravity.md (timeout, empty)
"""

D08 = """

## Pivot Decision: TB-GBF1-MULTICLASS-M8 (2026-06-12)
### Inputs: docs/07#tri-review-tb-gbf1-multiclass-m8 (2/3 DEGRADED, 2-0 consensus), docs/06 result. Track B, NON-CLAIM.
### DECISION (autonomous, 2/2 consensus): DROP multi-class (M8 bet refuted); NEXT AXIS = STAGED UNFREEZE / fine-tune SegmentNT, entered via a BOUNDED screen-profile PREFLIGHT (not a direct >24h scale-up).
- M8 result: multi-class structured output did NOT recover gbF1 on CLEAN held-out plants {arabidopsis,rice}: mc gbF1 0.7189 NOT > 3c 0.7392 (worse + gcount 0.66 under-prediction). The gbF1->ANNEVO-ceiling gap (~0.16) is NOT closed by richer decoder labels -> structural; the most likely cause (both reviewers) is the FROZEN features (ANNEVO ceiling 0.8976 is end-to-end-trained; frozen-head caps ~0.74).
- Clean POSITIVE (new honest headline): 3c-candidate (frozen SegmentNT + FP-aware + constrained) PARETO-beats the raw-DNA same-budget anchor on CLEAN plants on BOTH co-primary axes (spec 0.9663 vs 0.9045 +0.062; gbF1 0.7392 vs 0.6960 +0.043) — leakage-free (SegmentNT backbone excludes plants), replacing M7's chicken-contaminated +0.155 with an honest clean dual-axis win. The foundation-feature route is VALIDATED clean; the open problem is gbF1 headroom to SOTA.
### NEXT = ④ STAGED-UNFREEZE PREFLIGHT (bounded screen; the bounded screen compute is within autonomy, but the IMPLEMENTATION is non-trivial -> user go-ahead recommended; the full >24h unfreeze scale-up AFTER is a hard user gate):
  - Mechanism: unfreeze the TOP N layers of SegmentNT + the existing 3c FP-aware constrained head; low LR; backprop into the backbone. Tests directly whether frozen features are the gbF1 cap.
  - **IMPLEMENTATION REALITY (key)**: SegmentNT is JAX/Haiku; the current head is torch with SEPARATE frozen jax feature extraction. Unfreezing needs an in-process trainable path: either (a) a JAX/Haiku head + jax fine-tune of SegmentNT, or (b) a torch port of SegmentNT (e.g. HF AutoModel if available). This is a substantial new implementation (the M7/M8 jax-extract / torch-head split was deliberate because they don't coexist). So ④ is a real new goal, not a config flag.
  - Bounded-screen success (both reviewers): gbF1 directionally > frozen 3c 0.7392 on clean plants, intergenic_specificity not collapsing, gene_count sane (avoid mc-style under-call). NOT set near ANNEVO. If no directional gbF1 gain -> backbone-only domain-adaptation (masked/self-sup or pseudo-label) or a different foundation model.
  - Leakage discipline (codex): clean species/chrom split; no test labels in early-stop/decode-tuning; stay raw-DNA ab-initio; same 3-class collapse ruler; no test-truth gene_count calibration. Pre-claim: keep evaluating on segmentation-clean species (plants) — chicken/fly stay contaminated.
### NOT pursued: multi-class scaling (refuted); accept-frozen-ceiling (gap too big -> can't reach north star); evidence/RNA (breaks ab-initio purity). chicken-macrochromosome stratum DEFERRED (the binding finding is the clean-species gbF1 negative; chicken is contaminated so its robustness doesn't change the route pivot).
### Anti-tuning: gbF1 gap 0.16 >> 0.05 -> structural; unfreeze is an architecture axis (not lr/batch) -> compliant.
### Optional /revise-goal (human-gated): record the clean 3c-candidate dual-axis result (spec 0.966 / gbF1 0.739 on clean plants) + the frozen-feature gbF1 ceiling finding.
### Docs updated: 04,06,07,08,10,00.
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/07_tri_review.md", D07)
app("docs/08_pivot_decisions.md", D08)
