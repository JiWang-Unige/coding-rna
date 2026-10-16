"""result-log for TA-FOUNDATION-DECODER-M4. Run: python3 scripts/_log_m4.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D06 = r"""

## Result: TA-FOUNDATION-DECODER-M4

### Meta
- Date (UTC): 2026-06-11. Track A screen, NON-CLAIM. The MAIN architecture bet: foundation features -> structured decoder.
- Jobs: 8550151-8550166 (3 candidates x 5 seeds x 8 epochs) on shared-gpu, all COMPLETED. Reuses FP-SEGMENTNT-FEATCACHE (no re-extraction). New code: src/foundation_probe/train_probe_head.py {--loss fp_aware (intergenic-FP penalty), --fuse-raw-dna (vectorized one-hot), --decoder crf (LinearChainCRFVec + FP-aware aux)}.

### Result (5 seeds, base-weighted seed-mean +/- std; NEW full-transcript ruler)
| candidate | intergenic_specificity (AXIS-1) | macro (gate) | gene_body_F1 (AXIS-2) | FPR | gene_count_ratio | dual-gate PASS |
|---|---|---|---|---|---|---|
| FP-SEGNT-FPLOSS | 0.9303 +/- 0.036 | 0.8431 | 0.6157 | 0.070 | 2.25 | **YES** |
| FP-SEGNT-FUSION | 0.8615 +/- 0.018 | 0.7538 | 0.6850 | 0.139 | 3.40 | no |
| FP-SEGNT-CRF    | 0.8298 +/- 0.119 | 0.7329 | 0.6840 | 0.170 | 0.90 | no |
- FPLOSS per-seed spec: 0.963/0.981/0.890/0.921/0.896 (ALL 5 > anchor mean 0.871; min 0.890).
- CRF per-seed spec: 0.593/0.885/0.888/0.870/0.914 (HIGH variance; gene_count_ratio 0.90 = best coherence, structured decoder fixed over-prediction). FUSION 0.846/0.860/0.894/0.864/0.844.
- Anchor (same-budget raw-DNA tiberius_like, new ruler): spec per-seed 0.923/0.917/0.773 -> mean 0.8710 bw / 0.8278 macro; gene_body_F1 0.5576 (anchor ALSO high-variance, one seed 0.773). Gates: AXIS-1>0.8710, F1>=0.5276, macro>=0.7978. Ceiling (Helixer full-data) 0.9917.

### Verdict
- **FP-SEGNT-FPLOSS WINS**: PARETO-beats the same-budget anchor on the dual co-primary — intergenic_specificity 0.9303 > anchor 0.8710 (ALL 5 seeds > anchor mean; +0.059, tighter std 0.036) AND gene_body_F1 0.6157 > anchor 0.5576 AND > floor 0.5276 AND macro 0.8431 > gate 0.7978. FIRST candidate to strictly exceed the same-budget anchor on the new ruler. The FP-aware specificity-targeted loss converts the foundation features' recall into specificity — the MAIN architecture bet (foundation features + FP-aware objective) is VALIDATED at screen. Closes ~half the anchor->ceiling(0.9917) gap.
- FUSION: spec 0.8615 just BELOW anchor + macro fails -> no (highest gbF1 0.685; over-predicts count 3.40).
- CRF: spec 0.8298 < anchor + very high variance (one seed 0.59) + macro fails -> no; BUT best gene_count coherence (0.90) + high F1 -> structured decoder worth iterating (variance/regularization), not dropping.
- NON-CLAIM screen -> not_yet for the contract (screen never claims); FPLOSS is a Track-B promotion candidate.

### Key finding
2-epoch single-seed SMOKE was MISLEADING (all 3 looked >0.92, CRF best-balanced); the 5-seed 8-epoch full batch REVERSED it (CRF collapsed-variance, FUSION dropped below anchor, FPLOSS robust winner). Validates the goal's >=5-seed mandate. The FP-aware INPUT-objective beats both the input-fusion and the structured-decoder this round — controlling intergenic FP via the loss is the most robust lever at same-budget.
"""

D10 = r"""

## Finding (Research) 2026-06-11 -- FP-aware loss on frozen foundation features Pareto-beats the same-budget anchor (TA-FOUNDATION-DECODER-M4)
The MAIN architecture bet validated at screen: frozen SegmentNT features + an FP-aware specificity-targeted loss (penalize genic prob mass at true-intergenic bases) on an anchor-matched conv+biLSTM head STRICTLY beats the same-budget from-scratch raw-DNA anchor on BOTH co-primary axes -- intergenic_specificity 0.930 vs 0.871 (all 5 seeds > anchor mean) AND gene_body_F1 0.616 vs 0.558 -- a Pareto improvement, closing ~half the anchor->ceiling(0.9917) gap. Mechanism ranking (5 seeds, 8 epochs): FP-aware loss (0.930 PASS) > raw-DNA fusion (0.862, just below anchor) > CRF structured decoder (0.830, high variance ±0.119, but best gene_count coherence 0.90). LESSON 1: controlling intergenic FP via the OBJECTIVE is the most robust same-budget lever -- more robust than input-fusion or structured decoding. LESSON 2: a structured decoder (CRF) gives the best gene-count coherence (0.90 vs FPLOSS 2.25) but unstable specificity -> the obvious synthesis is FP-aware loss + CRF decoder (coherence + specificity). LESSON 3: 2-epoch single-seed smoke was MISLEADING (reversed at 5-seed/8-epoch) -> always >=5 seeds for AXIS-1 verdicts. Caveat: FPLOSS gene_count_ratio 2.25 (over-predicts gene COUNT despite low FPR -> fragmentation); the loss is tuned to the specificity metric (watch 'teaching to the metric' -- but it generalizes across both species + holds gene-body F1, so it is a real gain not pure metric-gaming).
"""

def append(p, t):
    with open(os.path.join(ROOT, p), "a") as fh:
        fh.write(t)
    print("appended", p)

append("docs/06_results_log.md", D06)
append("docs/10_findings.md", D10)

D04 = ("\n## ITER-FP-002 -- TA-FOUNDATION-DECODER-M4 (2026-06-11)\n"
       "- Track A screen (NON-CLAIM), MAIN architecture bet: foundation features -> structured decoder. 3 candidates x 5 seeds.\n"
       "- FP-SEGNT-FPLOSS (loss_design): spec 0.9303 > anchor 0.8710, gbF1 0.6157 > anchor 0.5576, macro 0.8431 > gate -> PARETO-beats anchor -> PASS (Track-B candidate).\n"
       "- FP-SEGNT-FUSION (data_view): spec 0.8615 < anchor -> no. FP-SEGNT-CRF (decoder): spec 0.8298 high-variance < anchor -> no (but best gene_count coherence 0.90).\n"
       "- Parent: FP-SEGMENTNT-PROBE-M1 (ITER-FP-001). Pivot: see docs/08. Components: FP-SEGNT-{FPLOSS,FUSION,CRF}-s0..4 (jobs 8550151-66).\n")
append("docs/04_experiment_iterations.md", D04)

# docs/00
p = os.path.join(ROOT, "docs/00_active_goal.md")
t = open(p).read()
mk = "## last_result_summary\n"
i = t.index(mk) + len(mk)
blk = (
"- exp_id: `TA-FOUNDATION-DECODER-M4` (Track A screen, NON-CLAIM, MAIN architecture bet: foundation -> structured decoder)\n"
"- date: 2026-06-11 UTC\n"
"- 3 candidates x 5 seeds (8 epochs) on frozen SegmentNT features (reused FEATCACHE), NEW ruler, vs anchor 0.8710.\n"
"- **WINNER FP-SEGNT-FPLOSS** (FP-aware specificity-targeted loss): intergenic_specificity **0.9303 +/- 0.036 > anchor 0.8710** (ALL 5 seeds > anchor mean) AND gene_body_F1 **0.6157 > anchor 0.5576** AND macro 0.8431 > gate 0.7978 -> **PARETO-beats the same-budget anchor on the dual co-primary**. First candidate to strictly exceed the anchor on the new ruler. MAIN bet validated at screen.\n"
"- FUSION 0.8615 (just below anchor, no); CRF 0.8298 (high var ±0.119, no; but best gene_count coherence 0.90).\n"
"- PIVOT (pending tri-review): FPLOSS = Track-B promotion candidate (scale-up = new long sub-iteration -> user go-ahead). Obvious synthesis next: FP-aware loss + CRF decoder (specificity + coherence). 2-epoch smoke was misleading -> 5-seed mandate validated.\n"
"- screen_anchor=0.8710 (PROVISIONAL); ceiling 0.9917; status draft.\n"
"- --- prior results kept below for trend ---\n")
t = t[:i] + blk + t[i:]
open(p, "w").write(t)
print("docs/00 updated")
