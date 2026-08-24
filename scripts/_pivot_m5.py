"""docs/07 + docs/08 for TA-COHERENCE-FIX-M5. python3 scripts/_pivot_m5.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D07 = r"""

# Tri-Review: TA-COHERENCE-FIX-M5  (2026-06-11)
## Mode: independent_parallel_cli | A·Claude success | B·Codex success | C·Antigravity(agy) success | Quorum: 3/3
## Subject: FP-FRAGFIX-CONSTR (FPLOSS + deterministic constrained post-proc) vs 5-seed anchor. spec 0.9272 (paired +0.0836±0.037 vs 5-seed anchor 0.8436, all 5 positive, Claude computed t≈5.0 p<0.01), gbF1 0.6581, macro 0.8555, gene_count 2.25->1.28 (0.03 over the full/scale guardrail 1.25).

## A·Claude — promote-to-track-b (with a Track-B job#0 gating step) (Medium-high)
Arithmetic verified (means, paired diffs all positive, t≈5.0 p<0.01, n=5 significant). Promotable. ATTRIBUTION: the +0.0836 vs anchor mixes FP-loss + SegmentNT pretraining dividend + constrained; the CLEAN net contribution of constrained-decode = vs M4 FPLOSS (spec 0.930->0.927 KEPT, gbF1 0.616->0.658 UP, gene_count 2.25->1.28 FIXED) — deterministic, no CRF instability. CONSTR seed spread (0.091) < anchor (0.151) = more stable. gene_count 1.28: do a quick param re-run to clear <=1.25 BEFORE Track B main training (NOT 'promote then tune') — UNTESTED tradeoff: max_fill_gap↑/min_cds_len↑ merges more into gene-body -> may RAISE intergenic FP / lower spec (spec & gene_count are COUPLED, not 'trivially tunable' for free); cheap (deterministic, no retrain) so make it Track-B job#0: sweep on TRAIN/VAL requiring gene_count<=1.25 AND spec>=anchor. screen_anchor: update to 5-seed 0.8436 via /revise-goal (keep 3-seed 0.8710 recorded + 'high variance' note) but it does NOT change promotion (paired test bypasses the point estimate). LEAKAGE PRECONDITION (must confirm): constrained params (max_fill_gap/min_cds_len) chosen on train/val NOT test — else test leakage blocker. (Confirmed by agent: defaults 30/20, same as train_screen_ref, never touched test -> OK.) unfreeze SegmentNT = LATER Track-B step (not mixed with first scale run — attribution). Only 2 species -> macro 'stability' semantics limited at screen. Confidence Medium(-high): paired stats robust, downside limited (screen non-claim, Track B retests).

## B·Codex — run-sanity-check-first then promote (Medium-High)
CONSTR is a real Track A winner (all 5 paired deltas positive, +0.0836, F1 up, fragmentation mostly fixed — not a tuning artifact, directly fixes the M4 blocker). But do NOT promote the exact gene_count=1.28 config into expensive Track B: do one cheap deterministic post-proc sweep (max_fill_gap/min_cds_len, paired seeds) requiring gene_count<=1.25 while keeping spec>0.8710 AND gbF1>=0.5276, THEN promote. Don't spend Track B compute on a config known to violate the full/scale guardrail (the exact M4 failure mode). screen_anchor: do NOT immediately revise to 0.8436 unless the 5-seed protocol is the official frozen one; CONSTR passes even the stricter old 0.8710 so conclusion unchanged; keep reporting both, /revise-goal only after confirming the 5-seed rerun is same code/data/metric. No hard blocker: deterministic constrained decode = acceptable test-time inference (frozen before claim); FP-loss lambda 1.0 not test-tuned fine; clean internal comparison = M4 FPLOSS vs CONSTR (anchor mixes pretraining dividend). Unfreeze SegmentNT = separate Track-B axis / staged ablation. Confidence Medium-High (high the route deserves Track B; medium the exact config is ready due to the small guardrail violation on the route's main failure mode).

## C·Antigravity(agy) — promote-to-track-b (High)
Robust winner: paired +0.0836±0.037 (all 5 positive), de-fragmented 2.25->1.28, gbF1 up 0.616->0.658. Meets Track-B bar (scale/multi-class/unfreeze). gene_count 1.28 overshoot tiny (0.03) -> tune max_fill_gap/min_cds_len in Track B or a quick local run. screen_anchor SHOULD update to 5-seed 0.8436 (more representative, higher variance) via /revise-goal; does NOT change conclusion (CONSTR 0.9272 > both 0.8436 and 0.8710). No blockers: constrained_decode deterministic (no test leakage), lambda 1.0 not test-tuned, M4 FPLOSS = clean internal baseline. Confidence High.

## Cross-reviewer agreement (3/3)
- CONSTR is a robust, real Track-A winner: paired-significant Pareto over the 5-seed anchor (+0.0836, all 5 positive, p<0.01) + de-fragmentation (2.25->1.28) + F1 kept/improved. Worth Track B.
- BEFORE Track-B main compute: do ONE CHEAP deterministic constrained-param sweep (max_fill_gap/min_cds_len) on TRAIN/VAL to clear gene_count<=1.25 while keeping spec>anchor + gbF1>=floor (A=Track-B job#0, B=sanity-first, C=quick-tune — same action). Spec & gene_count are COUPLED -> verify, don't assume free.
- constrained params must be non-test-tuned (= defaults 30/20, confirmed) — else leakage blocker.
- clean attribution baseline = M4 FPLOSS (anchor mixes SegmentNT pretraining dividend); report it.
- screen_anchor 5-seed 0.8436 < 3-seed 0.8710 -> /revise-goal update candidate (A/C yes, B after-confirm); does NOT change promotion (paired test bypasses point estimate).
- unfreeze SegmentNT = staged separate Track-B axis, NOT in the first scale run (attribution).

## Aggregated recommendation: PROMOTE-READY pending one cheap constrained-param sweep (train/val) to clear gene_count<=1.25 + keep spec>anchor. Then Track-B promote (= user go-ahead, new long sub-iteration). Confidence: Medium-High (3/3 promote-route; only the 0.03 gene_count overshoot + spec-gene_count coupling to verify cheaply first).
## Raw: /tmp/tri_review_TA-COHERENCE-FIX-M5/output_{a_claude,b_codex,c_antigravity}.md
"""

D08 = r"""

## Pivot Decision: TA-COHERENCE-FIX-M5 (2026-06-11)
### Inputs: docs/07#TA-COHERENCE-FIX-M5 (3/3 quorum), docs/06 result. Track A screen, NON-CLAIM.
### Result: FP-FRAGFIX-CONSTR (FPLOSS + deterministic constrained post-proc) PARETO-beats the 5-seed anchor, PAIRED-SIGNIFICANT: intergenic_specificity 0.9272 vs anchor 0.8436 (paired +0.0836±0.037, all 5 seeds positive, t~5.0 p<0.01), gbF1 0.6581 > anchor 0.5768, macro 0.8555 > gate, gene_count 2.25 -> 1.28 (de-fragmented, 0.03 above full/scale guardrail 1.25). Clean net constrained contribution (vs M4 FPLOSS): spec kept (0.930->0.927), F1 up (0.616->0.658), fragmentation fixed.
### DECISION (autonomous, 3/3 consensus): PROMOTE-READY pending ONE cheap constrained-param sweep. CONSTR is the validated winner; the only blemish (gene_count 1.28 > 1.25) is cleared by a deterministic, no-retrain param sweep BEFORE spending Track-B compute.
### Next (NEW goal -> ③; the cheap sweep is its first step):
  - STEP 0 (cheap, deterministic, gate): save raw pre-constrained per-base predictions (small code add) OR re-run the 5 CONSTR seeds; sweep constrained_decode (max_fill_gap/min_cds_len) on TRAIN/VAL, require gene_count_ratio<=1.25 AND intergenic_specificity>=anchor AND gbF1>=floor. (spec & gene_count are COUPLED — verify, since merging more into gene-body can raise intergenic FP.) Confirm params chosen on train/val NOT test (current defaults 30/20 are non-test — OK).
  - THEN promote to Track B (= USER GO-AHEAD, new long sub-iteration, >24h-compute exception): scale data/epochs/seeds + CI; add richer Tiberius-style multi-class output (CDS/intron/intergenic/phase/splice) — gives the structure real meaning; UNFREEZE/fine-tune SegmentNT as a SEPARATE staged axis (NOT mixed into the first scale run — attribution). Keep CONSTR (deterministic post-proc) as the coherence layer.
### Optional /revise-goal (human-gated, does NOT change promotion): update screen_anchor 0.8710(3-seed) -> 0.8436(5-seed, more representative + higher variance); keep both recorded. CONSTR beats both via the paired test.
### Recorded: constrained params = defaults (min_cds_len 30 / max_fill_gap 20), NOT test-tuned (no leakage). FP-loss lambda 1.0 hardcoded, not test-tuned. Clean attribution baseline = M4 FPLOSS (anchor comparison mixes SegmentNT pretraining dividend). CRFSTAB deferred (CONSTR superseded the need this round). FUSION dropped (M4).
### Docs updated: 04,05,06,07,08,10,00.
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/07_tri_review.md", D07)
app("docs/08_pivot_decisions.md", D08)
