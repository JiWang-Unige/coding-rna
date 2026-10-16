"""result-log for TA-FRAGFIX-SWEEP-M6. python3 scripts/_log_m6.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ids = " ".join(f"FP-FRAGFIX-CONSTR-rp-s{s}" for s in range(5))

D06 = f"""

## Result: TA-FRAGFIX-SWEEP-M6

### Meta
- Date (UTC): 2026-06-11. Track A screen, NON-CLAIM. STEP-0 promote-gate before Track B: clear FP-FRAGFIX-CONSTR gene_count 1.28->≤1.25 via a VAL-chosen (no test leakage) constrained-decode param sweep.
- Jobs: FP-FRAGFIX-CONSTR-rp-s0-4 (8552452-56, --save-raw-pred) COMPLETED. New code: train_probe_head --save-raw-pred (raw pre-constrained val+test per-seqid preds + val_eval_subsets); scripts/_sweep_constrained_m6.py (offline VAL grid sweep -> pick -> apply test, torch-free).

### Method (no test leakage)
Saved RAW (pre-constrained) per-seqid predictions for VAL+TEST. Swept constrained_decode (max_fill_gap in {{20,40,60,100,150}} x min_cds_len in {{30,60,90}}) on VAL only; chose max val_spec s.t. val gene_count<=1.25 -> max_fill_gap=20, min_cds_len=90 (val_spec 0.9349, val_gcount 0.966). Applied those params ONCE to TEST.

### Result (TEST, 5 seeds, VAL-chosen params mfg=20/mcl=90)
| metric | test mean ± std | per-seed |
|---|---|---|
| intergenic_specificity | 0.9262 ± 0.019 | 0.944/0.947/0.929/0.915/0.896 (ALL > anchor) |
| macro_intergenic_specificity | 0.8389 ± 0.042 | — |
| gene_body_F1_unconstrained | 0.6376 ± 0.015 | — |
| gene_count_ratio | 0.939 ± 0.281 | 1.348/1.038/0.703/0.553/1.053 |
- Anchor 5-seed: spec 0.8436 / macro 0.802 / gbF1 0.5768 / gene_count 2.89. 3-seed 0.8710. Ceiling 0.9917.

### Verdict — STEP-0 GATE CLEARED (all 4 on the mean)
- spec 0.9262 > anchor (0.8710 & 0.8436) PASS; gbF1 0.6376 >= 0.5276 PASS; macro 0.8389 >= 0.7978 PASS; gene_count 0.939 <= 1.25 PASS (from 1.28).
- FP-FRAGFIX-CONSTR is now PROMOTE-READY: paired-significant Pareto over the anchor (M5) + de-fragmented to within the full/scale gene_count guardrail, params chosen on VAL (no test leakage).
- CAVEAT: gene_count high seed variance (0.55-1.35); mcl=90 aggressive -> 2 seeds UNDER-predict (0.55/0.70<1.0, may merge/miss real genes); gbF1 slight drop vs M5 (0.658->0.638). A milder param (mfg=20/mcl=30, val_gcount 1.21) targets ratio≈1.0 — a Track-B tuning choice.

### Component exp_ids (ledger)
{ids} — jobs 8552452-56, all COMPLETED.
"""

D10 = r"""

## Finding (Engineering) 2026-06-11 -- VAL-chosen constrained params clear the coherence guardrail without test leakage (TA-FRAGFIX-SWEEP-M6)
To clear FP-FRAGFIX-CONSTR's gene_count 1.28 -> <=1.25 rigorously: saved RAW pre-constrained predictions (train_probe_head --save-raw-pred) for VAL+TEST, swept constrained_decode (max_fill_gap x min_cds_len) on VAL ONLY, chose params there (max val_spec s.t. val gene_count<=1.25 -> mfg=20/mcl=90), applied ONCE to TEST -> all 4 gates pass on the mean (spec 0.926, macro 0.839, gbF1 0.638, gene_count 0.939). Deterministic post-proc + VAL-only param selection = no test leakage (addresses the M5 reviewer concern). LESSON: the 'max-spec s.t. <=1.25' val rule lands AGGRESSIVE (mcl=90 drops CDS<90bp) -> high gene_count seed variance (0.55-1.35), some seeds under-predict (<1.0); a target-ratio≈1.0 rule (mfg=20/mcl=30) would be more biologically faithful. FP-FRAGFIX-CONSTR (frozen SegmentNT features + FP-aware loss + VAL-tuned constrained post-proc) is the promote-ready same-budget-screen winner: Pareto-beats the anchor on both co-primary axes AND within the full/scale coherence guardrail.
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/06_results_log.md", D06)
app("docs/10_findings.md", D10)
app("docs/04_experiment_iterations.md",
    "\n## ITER-FP-004 -- TA-FRAGFIX-SWEEP-M6 (2026-06-11)\n"
    "- Track A screen (NON-CLAIM), STEP-0 promote-gate: VAL-chosen constrained param sweep clears gene_count 1.28->0.939 (<=1.25).\n"
    f"- FP-FRAGFIX-CONSTR (rp, mfg=20/mcl=90): test spec 0.9262>anchor, gbF1 0.6376, macro 0.8389, gene_count 0.939 -> ALL 4 GATES PASS -> PROMOTE-READY.\n"
    f"- Parent: TA-COHERENCE-FIX-M5 (ITER-FP-003). Pivot: docs/08 (promote to Track B = user go-ahead). Component exp_ids: {ids}.\n")

p = os.path.join(ROOT, "docs/00_active_goal.md"); t = open(p).read()
mk = "## last_result_summary\n"; i = t.index(mk) + len(mk)
blk = (
"- exp_id: `TA-FRAGFIX-SWEEP-M6` (Track A screen, NON-CLAIM, STEP-0 promote-gate)\n"
"- date: 2026-06-11 UTC\n"
"- **STEP-0 GATE CLEARED -> FP-FRAGFIX-CONSTR is PROMOTE-READY.** VAL-chosen (no test leakage) constrained params mfg=20/mcl=90 -> TEST 5-seed: intergenic_specificity **0.9262** (>anchor 0.8710/0.8436, all 5 seeds), gene_body_F1 0.6376 (>floor/anchor), macro 0.8389 (>gate), **gene_count_ratio 1.28 -> 0.939 (<=1.25 guardrail)** -> ALL 4 gates PASS.\n"
"- CAVEAT: gene_count seed variance 0.55-1.35 (mcl=90 aggressive, some under-predict); milder mfg=20/mcl=30 (gcount~1.0) is a Track-B option.\n"
"- PIVOT (pending tri-review): PROMOTE FP-FRAGFIX-CONSTR to Track B = ③ (USER GO-AHEAD): scale data/epochs/seeds + Tiberius multi-class (CDS/intron/intergenic/phase/splice) + staged SegmentNT unfreeze.\n"
"- ladder (new ruler): FLOOR 0.8805 / anchor 0.8436(5s) / FP-FRAGFIX-CONSTR 0.926 / ceiling 0.9917. status draft.\n"
"- --- prior results kept below for trend ---\n")
open(p, "w").write(t[:i] + blk + t[i:]); print("docs/00 updated")
