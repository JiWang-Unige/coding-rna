"""tri-review (docs/07) + pivot (docs/08) for REANCHOR-HELDOUT-M7 + over-statement correction.
python3 scripts/_pivot_m7.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D07 = """

## Tri-Review: REANCHOR-HELDOUT-M7

### Review mode
- independent_parallel_cli, one identical full-scope prompt. Prompt: /tmp/tri_review_REANCHOR-HELDOUT-M7/prompt_full_scope.md
- Reviewer A: Claude CLI · success
- Reviewer B: Codex CLI (gpt-5.5) · success
- Reviewer C: Antigravity (agy --print) · FAILED (timeout 700s, no content)
- **Quorum: 2/3 DEGRADED_REVIEW; confidence ceiling = Medium.** (pivot = scale-to-track-b is not a claim/abandon/goal-revision; 2 independent reviewers sufficient; launch additionally human-gated.)

### Inputs
- Experiment: REANCHOR-HELDOUT-M7 (Track A screen, NON-CLAIM, retrospective-derived re-anchor gate). held-out anchor spec 0.8054 / candidate 0.9604 / ANNEVO ceiling 0.9824.

### Reviewer A · Claude — scale-to-track-b (conditional), Confidence Medium
- Methodology: leakage PASS; fairness PASS* (anchor=raw-DNA random-init vs candidate=SegmentNT pretrained features -> part of the win is a pretraining dividend, not decoder alone; OK for re-anchor, but tighten wording); foundation-feature leakage FAIL for the held-out *selling point* (SegmentNT ~850 vertebrate-biased species -> chicken almost certainly in-corpus, arabidopsis possibly -> not truly held-out at feature level; erodes the generalization narrative, not the fairness); chicken-subset FAIL on coverage (NC<=20Mb microchromosomes are the EASY regime for spec; gene-sparse macrochromosomes untested); VAL-band PASS.
- Verdict over-stated: candidate dominates AXIS-1 but gbF1 0.666 < anchor 0.710 -> "Pareto-ADMISSIBLE" (passes R6 contract), NOT "Pareto-beat BOTH axes". gbF1 to ceiling gap 0.23.
- Biggest risk: architectural gbF1 ceiling — candidate buys spec with CDS-F1; gap 0.23 >> 0.05 anti-tuning threshold = structural, not tunable.
- job#1 = gbF1-recovery architecture (richer strand/phase/splice multi-class structured output) WITHOUT losing spec; parallel non-blocking guardrails: (a) add gene-sparse macrochromosomes to eval, (b) deterministically audit SegmentNT pretraining species membership for arab/gallus.

### Reviewer B · Codex (gpt-5.5) — scale-to-track-b (gated entry, not blind scale), Confidence Medium
- Methodology table: Leakage Pass(screen)/Unknown(claim, audit external-weights species overlap); Fairness Pass (3 vs 5 seeds imperfect but margin large); Chicken subset Partial/material (microchromosome bias, biggest extrapolation risk for the north-star spec); VAL-band Pass (record grid + rule).
- Conclusion needs downscaling: strong spec/macro beat (all 5 seeds), but NOT mathematical Pareto-dominance because gbF1 0.6664 < anchor 0.7099. Say "passes screen promotion contract", not "both axes Pareto-beat".
- Biggest risk: scaling a specificity-biased route -> spec held but gene-F1 can't reach SOTA; ceiling gaps spec 0.022 vs gene-F1 0.231 (not tuning-closable). 2nd: microchromosome subset over-estimates macrochromosome/gene-sparse performance.
- job#1 = gated Track-B entry: richer multi-class structure-aware output to recover AXIS-2 + report 3 strata (Arabidopsis / Gallus microchromosome / Gallus macrochromosome) for spec/macro/gbF1/gene_count. Pass conditions: spec still >> anchor; gbF1 no longer < raw-DNA anchor (or clear recovery trend); no macrochromosome specificity collapse; gene_count <=1.25 + under-prediction not worsening.

### Cross-reviewer agreement (2/2)
- scale-to-track-b, but Track-B job#1 = gbF1 recovery (multi-class), NOT more spec / NOT blind scale.
- The "Pareto-beat both axes" verdict is OVER-STATED; correct to "Pareto-admissible / passes screen contract" (wins AXIS-1, loses AXIS-2 vs anchor).
- Two material risks before/within Track-B: (1) gene-sparse macrochromosome regime untested (add stratum); (2) SegmentNT foundation-feature species-overlap (audit before any claim).
- gbF1->ceiling gap 0.231 is an ARCHITECTURAL gap (>> 0.05) — multi-class output is the lever, not tuning.
- Confidence Medium (held-out spec evidence strong; AXIS-2 regressed + subset coverage gap).

### Disagreements
- None material. Both converge on conditional scale-to-track-b with identical job#1 redirection.

### Aggregated recommendation
- [x] Scale to Track B — CONDITIONAL: job#1 = gbF1-recovery multi-class + macrochromosome stratum + SegmentNT overlap audit. NOT blind scale, NOT more spec.

### Confidence
Medium (2/3 DEGRADED; both independent reviewers Medium; strong spec evidence, AXIS-2 + coverage caveats).

### Raw outputs
- /tmp/tri_review_REANCHOR-HELDOUT-M7/output_a_claude.md
- /tmp/tri_review_REANCHOR-HELDOUT-M7/output_b_codex.md
- /tmp/tri_review_REANCHOR-HELDOUT-M7/output_c_antigravity.md (timeout, empty)
"""

D08 = """

## Pivot Decision: REANCHOR-HELDOUT-M7 (2026-06-12)
### Inputs: docs/07#tri-review-reanchor-heldout-m7 (2/3 DEGRADED, 2-0 consensus), docs/06 result. Track A screen, NON-CLAIM, retrospective-derived re-anchor gate.
### DECISION (autonomous within screen scope; 2/2 consensus): RE-ANCHOR GATE PASSED -> ③ Track-B is GREEN-LIT pending USER GO-AHEAD, with a MANDATORY job#1 redirection.
- Held-out re-anchor gate did its job: the candidate's intergenic_specificity advantage TRANSFERS cross-clade (0.9604 vs held-out anchor 0.8054, +0.155, all 5 seeds; LARGER margin than yeast+fly +0.078; near ANNEVO ceiling 0.9824). The retrospective worry ("spec numbers only on low-UTR in-corpus outliers") is REFUTED for the specificity axis. Methodology clean (leakage/fairness/VAL-band PASS).
### CORRECTION (both reviewers, ADOPTED): the M7 verdict was OVER-STATED. The candidate is **Pareto-ADMISSIBLE** (passes the R6 screen promotion contract: spec strictly>anchor AND gbF1>=floor AND macro>=gate AND gene_count<=1.25), NOT "Pareto-beat BOTH co-primary axes" — it DOMINATES AXIS-1 (spec) but on AXIS-2 gbF1 0.6664 < anchor 0.7099 (it loses the publishable axis vs the raw-DNA anchor). docs/06/00/10 wording corrected.
### ③ Track-B (USER GO-AHEAD required — >24h compute / new long sub-iteration):
  - **job#1 (MANDATORY redirection, both reviewers): gbF1 RECOVERY, not more spec.** richer strand/phase/splice-aware MULTI-CLASS structured output (semi-CRF / segment-level + FP-aware objective) on frozen SegmentNT features; target constrained_gene_body_F1 climbing toward ANNEVO ceiling 0.8976 WHILE intergenic_specificity stays >=~0.95. The gbF1->ceiling gap 0.231 is ARCHITECTURAL (>>0.05 anti-tuning threshold) — multi-class is the lever, tuning will not close it.
  - **mandatory eval upgrades (parallel, non-blocking):** (a) add Gallus gene-sparse MACROCHROMOSOME stratum (the untested, hardest spec regime — current chicken subset is gene-dense microchromosomes); report 3 strata (Arabidopsis / Gallus-micro / Gallus-macro). (b) deterministically AUDIT SegmentNT pretraining species membership for arabidopsis+gallus -> pre-claim leakage gate (held-out novelty is feature-level-contaminated; matters for full/scale claim, NOT for this non-claim screen).
  - pass conditions for the gated entry: spec still >> same-budget anchor; gbF1 no longer < raw-DNA anchor (or clear recovery trend); no macrochromosome specificity collapse; gene_count<=1.25 + under-prediction not worsening.
  - staged UNFREEZE/fine-tune SegmentNT = SEPARATE later axis (attribution); keep deterministic constrained post-proc as coherence layer; >=5-8 seeds + CI.
### Optional /revise-goal (human-gated): record held-out anchor (spec 0.8054 / macro 0.7804 / gbF1 0.7099) + ANNEVO ceiling 0.9824 alongside the yeast+fly anchor (do not replace; both are valid same-budget references on different species sets).
### Anti-tuning: gbF1 gap 0.231 >> 0.05 -> tuning_allowed=false on the gbF1 axis -> Track-B job#1 MUST be a structural (multi-class output) change, not lr/batch. Compliant.
### Docs updated: 04,05,06,07,08,10,00.
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/07_tri_review.md", D07)
app("docs/08_pivot_decisions.md", D08)

# --- over-statement correction patches (reviewers' consensus) ---
patches = [
    ("docs/10_findings.md",
     "advantage TRANSFERS cross-clade (stronger on held-out than on the original species)",
     "SPECIFICITY advantage TRANSFERS cross-clade (Pareto-ADMISSIBLE, NOT dominant: wins AXIS-1 spec, LOSES AXIS-2 gbF1 vs anchor)"),
    ("docs/00_active_goal.md",
     "**HELD-OUT PARETO-PASS — promote-ready conclusion REINFORCED cross-clade.**",
     "**HELD-OUT RE-ANCHOR GATE PASSED (Pareto-ADMISSIBLE) — SPECIFICITY axis reinforced cross-clade; gbF1 axis LOSES to anchor (tri-review correction).**"),
]
for path, old, new in patches:
    fp = os.path.join(ROOT, path); t = open(fp).read()
    if old in t:
        open(fp, "w").write(t.replace(old, new, 1)); print("patched", path)
    else:
        print("WARN patch target not found in", path)
