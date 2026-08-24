"""docs/07 tri-review + docs/08 pivot for TA-FOUNDATION-DECODER-M4. python3 scripts/_pivot_m4.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D07 = r"""

# Tri-Review: TA-FOUNDATION-DECODER-M4  (2026-06-11)
## Mode: independent_parallel_cli | A·Claude success | B·Codex success | C·Antigravity(agy) success | Quorum: 3/3
## Subject: foundation features -> structured decoder (MAIN bet), 3 candidates x 5 seeds. FPLOSS PASSES dual co-primary (spec 0.9303>anchor 0.8710, gbF1 0.6157>0.5576, macro 0.8431>0.7978, Pareto). FUSION 0.8615<anchor. CRF 0.8298 high-var<anchor but best gene_count 0.90.

## A·Claude — iterate (synthesis; NOT promote-as-is) (Medium)
FPLOSS is the only dual-axis Pareto winner, BUT gene_count_ratio 2.25 = FRAGMENTATION -> the full/scale HARD guardrail (<=1.25) would BLOCK it; base-weighted spec + base gbF1 are both insensitive to fragmentation, only gene_count_ratio catches it (already alarming). Promoting a 2.25-fragmented winner sends up something that breaks at the promotion ruler. Discounts on the +0.059: (1) NOT same-n (FPLOSS 5 seeds vs anchor 3); anchor mean 0.871 is dragged by ONE collapse (0.773) — anchor's other 2 seeds = 0.923/0.917 ≈ FPLOSS median 0.921, so +0.059 is largely 'FPLOSS never collapses vs anchor collapsed once' -> MUST rerun anchor to 5 seeds + paired test. (2) anchor(from-scratch raw-DNA) vs FPLOSS(pretrained SegmentNT features) mixes in a pretraining dividend; the CLEAN comparison is FPLOSS vs M1 probe (same features+head): spec 0.842->0.930 (+0.088 from loss alone) = clean strong evidence. (3) FPLOSS has the LOWEST gbF1 (0.616) — trades recall for spec. CRF: KEEP — most informative; its only failure is ONE collapse (0.593); other 4 seeds mean ~0.889>anchor; gene_count_ratio 0.90 = ONLY structurally-correct candidate (decoder does its job). FPLOSS lacks coherence, CRF lacks anti-spillover -> textbook complementarity. FUSION: iterate (add FP-loss), not drop. Next: (1) explicit synthesis FP-aware-CE-as-CRF-aux + CRF decoder (specificity x coherence) + anchor-to-5-seeds + paired; (2) THEN richer Tiberius multi-class (phase/splice give the CRF transitions real meaning — strong synergy, likely the ceiling-approach). NOT chase ceiling 0.9917 (different regime). repro: verify FP-loss lambda NOT tuned on test (only real risk); train-only feature normalization; SegmentNT pretraining-species vs test-clade overlap for future claim.

## B·Codex — promote-to-track-b (FPLOSS primary) (Medium)
FPLOSS is a sound Track B candidate (passes the screen gate, all 5 seeds > anchor mean), but NOT claim-quality yet; +0.059 gap is small and anchor (3 seeds, one collapse) is high-variance — vs anchor's 2 non-collapsed seeds the advantage nearly vanishes -> 'promotable mechanism', not 'robustly proven superior'. CRF: keep secondary, do NOT promote naked (coherence 0.90 valuable but spec collapse + variance) -> iterate as FPLOSS+CRF with variance controls. FUSION: drop standalone (fails spec+macro, gene_count 3.40), maybe only with FPLOSS. Next: Track B scale FPLOSS as clean primary + small bounded Track A side-test FPLOSS+CRF (do NOT combine all 3 — muddies attribution). Mandatory Track-B diagnostics: CI/bootstrap; >=5-seed anchor rerun; per-species (yeast); FRAGMENTATION (gene_count 2.25 -> gene length / exon count / merged-span / transcript-span precision-recall); keep intergenic_FPR<=0.01 future guardrail (0.070 far from final). FP-loss legitimate UNLESS lambda tuned on test. Don't chase 0.9917 yet; next meaningful step = move promoted FPLOSS toward the real structured task (richer CDS/intron/intergenic/phase/splice labels). Confidence Medium (anchor seed count, variance, simplified labels, fragmentation unresolved).

## C·Antigravity(agy) — promote-to-track-b (FPLOSS) + concurrent Track A hybrid (High)
FPLOSS sound robust winner at this scale (all 5 seeds > anchor mean, worst 0.890 > anchor mean; only candidate passing all 3 gates). CRF: iterate not drop (gene_count 0.90 best coherence; stabilize via FP-aware aux + regularization to fix the 0.593 collapse). FUSION: DROP (fails spec+macro, gene_count 3.40, raw-DNA fusion adds params without a SOTA path). Next: promote FPLOSS to Track B (verify scalability) + concurrent Track A hybrid FPLOSS-loss + CRF-decoder (specificity + coherence 0.90). FP-loss legitimate (biological prior vs class imbalance, test unseen). gene_count 2.25 = fragmentation -> justifies CRF for structural constraints. 2-epoch smoke unreliable for structured decoders (reversed at 8 epochs). Confidence High (rigorous 5-seed protocol exposes base-metric-wins-FPLOSS vs structural-coherence-wins-CRF trade-off).

## Cross-reviewer agreement (3/3)
- FPLOSS is the lead: only candidate passing all 3 dual-co-primary gates; Pareto-beats anchor; robust 5/5 (never collapses). Direction (foundation features + FP-aware objective) VALIDATED at screen.
- CRITICAL shared concern: FPLOSS gene_count_ratio 2.25 = FRAGMENTATION -> would fail the full/scale HARD gene_count guardrail (<=1.25); base-weighted spec & base gbF1 are blind to it.
- CRF: KEEP + iterate (best gene_count coherence 0.90; failure = ONE collapsed seed / variance, not direction). FUSION: not standalone (B/C drop; A iterate-with-FP-loss).
- anchor only 3 seeds + one collapse (0.773) inflates the +0.059 -> rerun anchor to 5 seeds + paired test before strong claims.
- FP-aware loss is LEGITIMATE cost-sensitive learning, NOT cheating — IFF lambda was not tuned on test (it was hardcoded 1.0, not tuned — confirmed).
- Don't chase ceiling 0.9917 (different pretrained+full-data regime). Richer multi-class (phase/splice) = strong synergy with CRF transitions, next orthogonal axis.

## Split: B,C = promote-to-track-b (FPLOSS) + concurrent FPLOSS+CRF side test; A = iterate (synthesis + fragmentation fix + anchor-5-seed) FIRST. A's fragmentation->full/scale-guardrail point is a genuine promotion blocker B also flagged (mandatory fragmentation diagnostics).

## ⚠️ Agent correction the reviewers missed: the CRF candidate ALREADY = --loss fp_aware --decoder crf (= the 'FPLOSS+CRF synthesis' they recommend). So the synthesis WAS tried: FP-loss alone (FPLOSS) = spec 0.930 / fragmented 2.25; FP-loss + CRF (CRF cand) = coherent 0.90 / spec 0.830+collapse. The learned CRF decoder TRADES specificity for coherence + adds variance. So the next step is NOT naive 'combine' — it is (a) FPLOSS + cheap constrained-decode POST-PROCESSING (merge fragments, no learned-CRF instability; reuse the old CONSTR mechanism) and/or (b) stabilize the CRF (diagnose the collapse seed / regularize / warm-start emissions), plus anchor-to-5-seeds.

## Aggregated recommendation: FPLOSS = validated screen winner (lead). Pivot = ITERATE one cheap screen round to FIX FRAGMENTATION (FPLOSS + constrained-decode post-proc; and/or stabilized CRF) + rerun anchor to 5 seeds + transcript-level/fragmentation diagnostics, THEN promote. Track-B promote-as-is withheld due to the 2.25 fragmentation (full/scale gene_count guardrail blocker). Confidence: Medium (2 promote / 1 iterate; fragmentation is a real promotion blocker).
## Raw: /tmp/tri_review_TA-FOUNDATION-DECODER-M4/output_{a_claude,b_codex,c_antigravity}.md
"""

D08 = r"""

## Pivot Decision: TA-FOUNDATION-DECODER-M4 (2026-06-11)
### Inputs: docs/07#TA-FOUNDATION-DECODER-M4 (3/3 quorum), docs/06 result. Track A screen, NON-CLAIM.
### Result: FPLOSS (FP-aware specificity-targeted loss on frozen SegmentNT features) PARETO-beats the same-budget anchor on the dual co-primary: intergenic_specificity 0.9303 > 0.8710 (all 5 seeds > anchor mean), gene_body_F1 0.6157 > 0.5576, macro 0.8431 > 0.7978. FIRST candidate to strictly exceed the anchor on the new ruler -> MAIN architecture bet (foundation features + FP-aware objective) VALIDATED at screen. FUSION 0.8615 (just below anchor) no; CRF 0.8298 (high variance ±0.119, one seed 0.59) no but best gene_count coherence 0.90.
### DECISION (autonomous, reconciling 3/3 split 2 promote / 1 iterate): ITERATE one cheap screen round to FIX FRAGMENTATION before Track-B promotion. FPLOSS is the validated lead; promote-as-is WITHHELD.
### Rationale: ALL 3 reviewers flagged FPLOSS gene_count_ratio 2.25 = fragmentation; the full/scale HARD guardrail predicted_gene_count_ratio<=1.25 would BLOCK it -> promoting a 2.25-fragmented candidate sends up something that fails at the promotion ruler. base-weighted spec + base gbF1 are blind to fragmentation; only gene_count_ratio catches it. AGENT CORRECTION (reviewers missed): the CRF candidate ALREADY = FP-loss + CRF decoder (their recommended 'synthesis'), and it TRADED specificity for coherence (spec 0.830, coherent 0.90, +variance). So the learned CRF is not a free coherence fix.
### Next round (NEW goal; cheap Track A screen, reuse FEATCACHE, same protocol):
  - Direction A (MAIN): FPLOSS (FP-aware loss winner) + CHEAP constrained-decode POST-PROCESSING (reuse src/screen_anchor/decoders.py constrained_decode: merge small intergenic gaps / drop tiny CDS) to fix the 2.25 fragmentation WITHOUT the learned-CRF's specificity cost + variance. Target: keep spec ~0.93, pull gene_count_ratio 2.25 -> <=1.25.
  - Direction B (stabilize the structured decoder): diagnose the CRF collapsed seed (0.593); regularize / warm-start emissions from FPLOSS / stronger FP-aware weight inside CRF, to recover spec while keeping coherence 0.90.
  - Direction C (validity, REQUIRED): rerun the ANCHOR to 5 seeds (currently 3, one collapse 0.773) -> valid paired test vs FPLOSS on intergenic_specificity. Report transcript-span precision/recall + gene-length/exon-count distributions (fragmentation diagnostics, per B).
  - >=5 seeds + CI. Confirm FP-loss lambda was NOT tuned on test (it was hardcoded 1.0 — confirmed; document).
### THEN promote the coherence-fixed winner to Track B (scale-up = new long sub-iteration -> USER GO-AHEAD required per decision-autonomy >24h-compute exception). Track-B job#1: scale data/epochs/seeds + CI + the richer Tiberius-style multi-class output (CDS/intron/intergenic/phase/splice) — multi-class gives the CRF transitions real biological meaning (reviewer A: likely the ceiling-approach step).
### DEFER: chasing ceiling 0.9917 (different pretrained+full-data regime, not screen-comparable); unfreeze/fine-tune SegmentNT (Track B); GENERanno (parallel probe). Pre-claim guard: verify test clade not in SegmentNT pretraining.
### FUSION: dropped as a standalone mechanism (fails spec+macro, worst fragmentation 3.40); its data_view idea may return only fused with FP-loss later.
### Anti-tuning: gap FPLOSS->anchor is FPLOSS ABOVE anchor (it won); the iterate is a structural/coherence fix (decoder/post-proc axis), not lr/batch tuning -> compliant.
### Docs updated: 04,05,06,07,08,10,00.
"""

def append(p, t):
    with open(os.path.join(ROOT, p), "a") as fh:
        fh.write(t)
    print("appended", p)

append("docs/07_tri_review.md", D07)
append("docs/08_pivot_decisions.md", D08)
