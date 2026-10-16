"""result-log for REANCHOR-HELDOUT-M7. python3 scripts/_log_m7.py (run anywhere; only edits docs)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cand_ids = " ".join(f"FP-FRAGFIX-CONSTR-ho-s{s}" for s in range(5))
anch_ids = " ".join(f"SCREENREF-tiberius_like-ho-s{s}" for s in range(3))

D06 = f"""

## Result: REANCHOR-HELDOUT-M7 (held-out/UTR-rich cross-clade re-anchor + candidate re-test)

### Meta
- Date (UTC): 2026-06-12. Track A screen, NON-CLAIM. Retrospective-2026-06-11-derived re-anchor GATE before Track-B (③). Submit-and-handoff (shared-gpu AMPERE; private-teodoro full).
- Purpose: every prior spec number was on yeast+fly (low-UTR, gene-dense, in-corpus OUTLIERS). Re-derive screen_anchor + ceiling on held-out/UTR-rich cross-clade species and re-test the promote-ready candidate FP-FRAGFIX-CONSTR — does its intergenic_specificity advantage survive cross-clade?
- Species (held-out, UTR-rich; WebFetch-verified RefSeq): Arabidopsis thaliana GCF_000001735.4 (TAIR10.1 full, 7 seqids/119Mb, UTR 41.7% of exon) + Gallus gallus GCF_016699485.2 (bGalGal1 GRCg7b, SUBSET NC_<=20Mb = 30 seqids/182Mb for screen cost; UTR 62.3% of exon). vs yeast+fly UTR~0. check_data PASS (no seqid leakage, all split pairs). chicken subset rule: assembled NC_ chromosomes <=20Mb (drops 6 macrochromosomes>20Mb + 172 NW_ scaffolds; full preserved); gene-dense microchromosomes -> gene-sparse macrochromosome intergenic-spec is a Track-B/full concern.

### Method (same unified screen protocol, ONLY species changed; anchor & candidate use SAME data -> fair)
- Held-out anchor: random-init tiberius_like, 3 seeds, sample_fraction 0.3, eval new full-transcript intergenic ruler (--span-mode cds).
- Candidate: FP-FRAGFIX-CONSTR IDENTICAL config (frozen SegmentNT feats + FP-aware loss convlstm head + constrained post-proc), 5 seeds, --save-raw-pred. Constrained band RE-SELECTED on held-out VAL (two-sided [1.0,1.25], max val_spec) -> mfg=20/mcl=60 (SAME as M6 — coherence params transfer), applied to TEST once (no leakage).
- Ceiling: ANNEVO clade-matched (Magnoliopsida/Aves), eval on the SAME test chromosomes (anchor eval_subsets) + new ruler. Helixer/Tiberius ceiling deferred (weights need download; ANNEVO=published-SOTA candidate is the key reference).

### Result (TEST, base-weighted) — held-out ladder
| metric | anchor (3-seed) | candidate (5-seed) | ANNEVO ceiling |
|---|---|---|---|
| intergenic_specificity | 0.8054 +-0.027 | **0.9604 +-0.0076** (all 5 > anchor) | 0.9824 |
| macro_intergenic_specificity | 0.7804 +-0.034 | **0.9621 +-0.0077** | 0.9781 |
| gene_body_F1_unconstrained | 0.7099 +-0.022 | 0.6664 +-0.012 | 0.8976 |
| predicted_gene_count_ratio | 1.9685 | 0.9688 +-0.078 | 0.732 |
- per-species spec: anchor arab 0.892 / gallus 0.669; candidate arab 0.954 / gallus **0.970** (gallus +0.30 vs anchor — vertebrate cross-clade gain is huge); candidate per-species gbF1 arab 0.792 / gallus 0.533 (vertebrate big-intron CDS harder).

### Verdict — HELD-OUT PARETO-PASS (re-anchor gate CLEARED; promote-ready conclusion REINFORCED cross-clade)
- spec 0.9604 STRICTLY > held-out anchor 0.8054 (+0.155, all 5 seeds, std 0.0076 -> paired-significant); macro 0.9621 > anchor 0.7804 (+0.18); gbF1 0.6664 > floor 0.5276; gene_count 0.9688 <= 1.25 HARD guardrail (mild 3% under-prediction, VAL-selected no leakage — benign direction; the over-prediction guardrail is cleared with huge margin).
- The candidate's intergenic_specificity advantage is NOT a yeast+fly artifact: on held-out UTR-rich cross-clade species the margin over the same-budget anchor is LARGER (+0.155) than on yeast+fly (+0.078), and absolute spec is HIGHER (0.9604 vs 0.9218), nearly reaching the pretrained ANNEVO ceiling (0.982). foundation-features + FP-aware objective TRANSFERS cross-clade — retrospective concern positively refuted.
- CAVEAT (Track-B): (1) gene_count mild under-prediction on held-out (mcl=60 from yeast+fly slightly aggressive cross-clade) -> per-clade band calibration with more data. (2) chicken subset is gene-dense microchromosomes; gene-sparse macrochromosome intergenic-spec untested (Track-B/full). (3) candidate gbF1 (0.666) < anchor gbF1 (0.710) — the candidate trades CDS-F1 for specificity (its design); SOTA-comparable claim axis needs the multi-class output planned for Track-B.

### Component exp_ids (ledger)
candidate (5): {cand_ids} (jobs 8554530, 8555770-73). anchor (3): {anch_ids} (jobs 8554369, 8554520-21). ceiling: REANCHOR-CEILING-ANNEVO-M7 (8554546). featcache: FP-SEGMENTNT-FEATCACHE-M7 (8554368).
"""

D10 = r"""

## Finding (Research) 2026-06-12 -- foundation+FP-aware candidate's intergenic-specificity advantage TRANSFERS cross-clade (stronger on held-out than on the original species) (REANCHOR-HELDOUT-M7)
The promote-ready candidate (frozen SegmentNT features + FP-aware loss + constrained post-proc) was re-tested on held-out UTR-rich cross-clade species (Arabidopsis thaliana plant + Gallus gallus vertebrate; UTR 42%/62% of exon vs yeast+fly ~0). Result: intergenic_specificity 0.9604+-0.008 (5 seeds) vs same-budget held-out anchor 0.8054 (+0.155, ALL seeds positive) — a LARGER margin than on yeast+fly (+0.078) and a HIGHER absolute (0.9604 vs 0.9218), nearly reaching the pretrained ANNEVO ceiling (0.982). The per-species vertebrate (chicken) gain is the largest (+0.30 spec vs anchor). LESSON: the retrospective worry that "all spec numbers are on low-UTR gene-dense in-corpus outliers (yeast+fly)" is POSITIVELY REFUTED — the architecture bet generalizes cross-clade and the held-out same-budget anchor is LOWER (0.805 vs 0.844/0.871), so the candidate's true margin was UNDER-stated before, not over-stated. De-risks the >24h Track-B spend.

## Finding (Engineering) 2026-06-12 -- VAL-selected constrained band transfers but has a val->test gene_count gap on held-out (REANCHOR-HELDOUT-M7)
Re-selecting constrained_decode params on held-out VAL (two-sided band [1.0,1.25], max val_spec) picked mfg=20/mcl=60 — IDENTICAL to the M6 yeast+fly choice (coherence params transfer cross-clade). But VAL gcount 1.070 (in-band) -> TEST gcount 0.969 (3% under 1.0): a real val/test chromosome gene-density generalization gap, NOT fixable without test-set tuning (leakage). LESSON: the constrained post-proc mcl=60 mildly under-predicts on held-out test chromosomes; Track-B should do per-clade band calibration with more val data, and report gene_count with the val->test gap explicit. Eval pitfall fixed en route: pretrained-baseline ceiling eval MUST filter predictions to the test seqids before scoring against a test-chrom genome, else whole-genome predicted-genic bases are counted against the test-chrom intergenic denominator -> FPR>1 -> negative specificity (caught & fixed; ANNEVO ceiling 0.982 after filtering).
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/06_results_log.md", D06)
app("docs/10_findings.md", D10)
app("docs/04_experiment_iterations.md",
    "\n## ITER-FP-005 -- REANCHOR-HELDOUT-M7 (2026-06-12)\n"
    "- Track A screen (NON-CLAIM), retrospective-derived re-anchor GATE before ③ Track-B. Held-out/UTR-rich cross-clade {Arabidopsis thaliana, Gallus gallus(subset)}.\n"
    "- Held-out same-budget anchor (tiberius_like 3-seed): spec 0.8054 / macro 0.7804 / gbF1 0.7099. ANNEVO ceiling spec 0.9824.\n"
    "- Candidate FP-FRAGFIX-CONSTR (5-seed, IDENTICAL config, VAL-selected mfg=20/mcl=60): spec 0.9604+-0.008 (all 5 > anchor, +0.155), macro 0.9621, gbF1 0.6664, gcount 0.9688 (<=1.25). HELD-OUT PARETO-PASS.\n"
    "- Margin LARGER cross-clade (+0.155) than yeast+fly (+0.078); absolute HIGHER (0.9604 vs 0.9218). Retrospective concern REFUTED. Parent: TA-FRAGFIX-SWEEP-M6 (ITER-FP-004). Pivot: docs/08.\n"
    f"- Component exp_ids: {cand_ids} {anch_ids} REANCHOR-CEILING-ANNEVO-M7 FP-SEGMENTNT-FEATCACHE-M7.\n")

p = os.path.join(ROOT, "docs/00_active_goal.md"); t = open(p).read()
mk = "## last_result_summary\n"; i = t.index(mk) + len(mk)
blk = (
"- exp_id: `REANCHOR-HELDOUT-M7` (Track A screen, NON-CLAIM, retrospective-derived re-anchor GATE before ③ Track-B)\n"
"- date: 2026-06-12 UTC\n"
"- **HELD-OUT PARETO-PASS — promote-ready conclusion REINFORCED cross-clade.** On held-out/UTR-rich {Arabidopsis(plant), Gallus(vertebrate)} (UTR 42%/62% of exon vs yeast+fly ~0): candidate FP-FRAGFIX-CONSTR (5-seed, IDENTICAL config) intergenic_specificity **0.9604+-0.008** (all 5 > held-out anchor 0.8054, +0.155) > macro 0.9621 (>0.7804) > gbF1 0.6664 (>floor) ; gene_count 0.9688 (<=1.25, mild 3% under-pred, VAL-selected no leakage). ANNEVO ceiling 0.9824.\n"
"- KEY: margin over same-budget anchor LARGER cross-clade (+0.155) than yeast+fly (+0.078); absolute spec HIGHER (0.9604 vs 0.9218), near ceiling. Retrospective worry (numbers on low-UTR outliers) POSITIVELY REFUTED.\n"
"- CAVEAT (Track-B): mild held-out gene_count under-prediction (per-clade band calib); chicken subset gene-dense (macrochromosome intergenic untested); candidate gbF1<anchor (trades CDS-F1 for spec; multi-class output planned).\n"
"- NEXT: ③ Track-B promotion now de-risked on held-out — USER GO-AHEAD required (>24h compute). Optional /revise-goal to record held-out anchor 0.8054 / ceiling 0.9824 alongside yeast+fly.\n"
"- --- prior results kept below for trend ---\n")
open(p, "w").write(t[:i] + blk + t[i:]); print("docs/00 updated")
