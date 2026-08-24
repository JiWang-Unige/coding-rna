"""result-log for TB-GBF1-MULTICLASS-M8. python3 scripts/_log_m8_result.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D06 = """

## Result: TB-GBF1-MULTICLASS-M8 (③ Track-B — multi-class structured output for gbF1 recovery, on CLEAN held-out plants)
### Meta
- Date (UTC) 2026-06-12. Track B scale-up, NON-CLAIM (M2 sota_benchmark pending). submit-and-handoff (shared-gpu AMPERE). promoted_from FP-FRAGFIX-CONSTR via REANCHOR-HELDOUT-M7.
- Goal: recover the gbF1 short-fall (M7: candidate gbF1 0.666 << ANNEVO ceiling 0.8976, an ARCHITECTURAL gap) via richer strand-aware MULTI-CLASS structured output (8-class: intergenic/CDS-phase0-2/intron/UTR/donor/acceptor) + CRF(8x8) transitions on frozen SegmentNT features. Evaluated on CLEAN held-out plants {arabidopsis, rice} (both-layer SegmentNT-clean — backbone excludes plants, segmentation fine-tune species {human,mouse,chicken,fly,zebrafish,worm} has no plant), because the M7 SegmentNT audit found chicken/fly are segmentation-FINE-TUNE CONTAMINATED.
- Multi-class label code: build_labels_multiclass + collapse_mc_to_3class (collapse→3class ZERO-mismatch vs the 3-class builder, IoU 1.0 — eval ruler unchanged). train_probe_head --label-scheme multiclass + --decoder crf. Data: rice GCF_034140825.1 subset 139Mb + arabidopsis full 119Mb, chromosome-level split, check_data PASS.

### Result (CLEAN held-out {arabidopsis, rice}, base-weighted; anchor n=3, 3c n=3, mc n=5)
| model | intergenic_specificity | gene_body_F1 | gene_count_ratio |
|---|---|---|---|
| raw-DNA anchor (tiberius_like, 3c) | 0.9045 +-0.018 | 0.6960 +-0.010 | 3.46 |
| 3c-candidate (FP-FRAGFIX-CONSTR, M7 cfg) | 0.9663 +-0.008 | 0.7392 +-0.006 | 0.936 |
| **mc-candidate (M8 multi-class+CRF)** | 0.9683 +-0.011 | **0.7189 +-0.022** | 0.6625 |
- per-species gbF1: anchor arab 0.783/rice 0.565; 3c arab 0.805/rice 0.645; mc arab 0.773/rice 0.640. per-species spec: 3c arab 0.959/rice 0.969; mc arab 0.945/rice 0.977.
- 2 of 5 3c seeds FAILED on a TRANSIENT beegfs/conda read error (numpy/_function_base_impl.py FileNotFoundError under concurrent env reads — not a code bug); n=3 signal is low-variance and the verdict is robust.

### Verdict — M8 PRIMARY BET FAILED (key negative result) + clean POSITIVE side-finding
- ❌ **Multi-class did NOT recover gbF1**: mc gbF1 0.7189 is NOT > 3c gbF1 0.7392 (−0.020, slightly worse), and mc gene_count 0.66 = SEVERE under-prediction (8-class CRF over-merges; constrained mcl=60 tuned for 3-class is wrong for mc). The gbF1->ceiling gap (~0.16) is NOT closed by richer decoder labels. The M8 hypothesis (structured multi-class output is the gbF1 lever) is REFUTED on clean species.
- ✅ **Clean POSITIVE (leakage-free)**: the 3c-candidate (frozen SegmentNT + FP-aware + constrained) PARETO-beats the raw-DNA anchor on CLEAN plants on BOTH co-primary axes: spec 0.9663 > 0.9045 (+0.062) AND gbF1 0.7392 > 0.6960 (+0.043), with far better coherence (gcount 0.94 vs anchor 3.46). And SegmentNT backbone NEVER saw plants — so this is a genuinely clean foundation-features-help signal (cleaner than M7's chicken-contaminated +0.155 headline). The honest cross-clade lead is the 3c-candidate on clean plants.
- IMPLICATION: gbF1 short-fall is structural and NOT addressed by frozen-feature + richer decoder. Next lever per protocol negative-result branch: staged UNFREEZE/fine-tune SegmentNT (frozen features likely cap gbF1) OR backbone-only self-trained head — a SEPARATE architecture axis. Multi-class output is NOT scaled.

### Component exp_ids (ledger)
mc(5): M8-MC-CAND-s0..4 (8559000-04). 3c(3 ok / 2 transient-fail): M8-3C-CAND-s0/s2/s4 ok, s1/s3 FAILED (8559005-09). anchor(3): SCREENREF-tiberius_like-m8clean-s0..2 (8558997-99). featcache: FP-SEGMENTNT-FEATCACHE-M8 rice (8558832). smoke: M8-MC-SMOKE.
"""

D10 = r"""

## Finding (Research) 2026-06-12 -- multi-class structured output does NOT recover gbF1 on clean species; but frozen-SegmentNT+FP-aware DOES cleanly Pareto-beat raw-DNA on plants the backbone never saw (TB-GBF1-MULTICLASS-M8)
Tested whether richer strand-aware MULTI-CLASS output (8-class CDS-phase/intron/UTR/splice + CRF 8x8 transitions on frozen SegmentNT features) recovers the gbF1 short-fall (M7: ~0.16 below the ANNEVO ceiling 0.8976). On CLEAN held-out plants {arabidopsis, rice} (both-layer SegmentNT-clean): mc-candidate gbF1 0.719 is NOT better than the 3-class candidate 0.739 (slightly worse + gene_count 0.66 severe under-prediction from over-merging). The gbF1 gap is NOT closed by richer decoder labels -> the M8 architecture bet is REFUTED. LESSON: structured multi-class output is not the gbF1 lever; the short-fall is more fundamental (likely the FROZEN features cap gene-level F1). Next axis = staged UNFREEZE/fine-tune SegmentNT or backbone-only self-trained head. SIDE-FINDING (clean, important): the 3-class candidate (frozen SegmentNT + FP-aware + constrained) Pareto-beats the same-budget raw-DNA anchor on BOTH co-primary axes on CLEAN plants (spec 0.966 vs 0.905 +0.062; gbF1 0.739 vs 0.696 +0.043) — and SegmentNT's backbone EXCLUDED plants, so this is a genuinely leakage-free "foundation features transfer to an unseen kingdom" result, replacing M7's chicken-contaminated +0.155 headline with an honest clean +0.06/+0.04 dual-axis win. The foundation-feature route is validated clean; multi-class is a dead end for gbF1.
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/06_results_log.md", D06)
app("docs/10_findings.md", D10)
app("docs/04_experiment_iterations.md",
    "\n## ITER-FP-006 -- TB-GBF1-MULTICLASS-M8 (2026-06-12)\n"
    "- Track B scale-up (NON-CLAIM), ③: multi-class structured output for gbF1 recovery on CLEAN held-out plants {arabidopsis,rice}.\n"
    "- NEGATIVE (primary): mc-candidate gbF1 0.7189 NOT > 3c-candidate 0.7392 (multi-class did NOT recover gbF1; gcount 0.66 under-pred). M8 bet REFUTED.\n"
    "- POSITIVE (clean): 3c-candidate Pareto-beats raw-DNA anchor on clean plants BOTH axes (spec 0.966>0.905, gbF1 0.739>0.696) — leakage-free (SegmentNT backbone excludes plants).\n"
    "- Next axis: staged UNFREEZE/fine-tune SegmentNT or backbone-only self-train (frozen features cap gbF1). Parent: REANCHOR-HELDOUT-M7. Pivot: docs/08.\n"
    "- Components: M8-MC-CAND-s0..4, M8-3C-CAND-s0/s2/s4 (s1/s3 transient-fail), SCREENREF-tiberius_like-m8clean-s0..2.\n")

p = os.path.join(ROOT, "docs/00_active_goal.md"); t = open(p).read()
mk = "## last_result_summary\n"; i = t.index(mk) + len(mk)
blk = (
"- exp_id: `TB-GBF1-MULTICLASS-M8` (③ Track-B scale-up, NON-CLAIM)\n"
"- date: 2026-06-12 UTC\n"
"- **M8 PRIMARY BET FAILED (key negative result)**: multi-class structured output did NOT recover gbF1 on CLEAN held-out plants {arabidopsis,rice}: mc-candidate gbF1 0.7189 NOT > 3c-candidate 0.7392 (−0.020, gcount 0.66 under-pred). gbF1->ceiling gap (~0.16) NOT closed by richer decoder labels -> structural, frozen features likely cap gbF1.\n"
"- **CLEAN POSITIVE side-finding**: 3c-candidate (frozen SegmentNT+FP-aware+constrained) PARETO-beats raw-DNA anchor on clean plants BOTH axes (spec 0.966 vs 0.905 +0.062; gbF1 0.739 vs 0.696 +0.043), leakage-free (SegmentNT backbone excludes plants) — honest replacement for M7's chicken-contaminated +0.155 headline.\n"
"- NEXT (pending tri-review/pivot): multi-class NOT scaled; next gbF1 axis = staged UNFREEZE/fine-tune SegmentNT OR backbone-only self-trained head (route-level, >24h, USER go-ahead).\n"
"- --- prior results kept below for trend ---\n")
open(p, "w").write(t[:i] + blk + t[i:]); print("docs/00 updated")
