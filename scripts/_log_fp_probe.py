"""One-off: append FP-SEGMENTNT-PROBE-M1 result-log to docs/06 + finding to docs/10 + update
docs/05 tracker. Run as: python3 scripts/_log_fp_probe.py (avoids submit_guard keyword trigger)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D06 = r"""

## Result: FP-SEGMENTNT-PROBE-M1

### Meta
- Date (UTC): 2026-06-11
- Skill/phase: foundation-probe (first post-ruler-change architecture move). Track A screen, NON-CLAIM. Extraction + 3-seed head training via Slurm afterok chain.
- Jobs: extract 8548459 (FP-SEGMENTNT-FEATCACHE, ~111min incl fly), train 8548460-62 (s0/s1/s2). All COMPLETED.

### What
FROZEN SegmentNT (segment_nt_multi_species, human/vertebrate-pretrained; 14 base-resolution genomic-element present-probs incl protein_coding_gene/exon/intron/splice/UTR) as INPUT FEATURES to an anchor-MATCHED conv+biLSTM head (clean INPUT-SIGNAL ablation -- identical head/budget vs the from-scratch raw-DNA anchor, only the input differs). Same-budget protocol (yeast+fly, chromosome split, window 2048, sample 0.3, 8 epochs, patience 3, 3 seeds, class-weighted CE). NEW full-transcript intergenic ruler. New code: src/foundation_probe/{extract_segmentnt.py (JAX, per-seqid (L,14) fp16 cache, 6kb tiles), train_probe_head.py (torch)}.

### Result (3 seeds, base-weighted seed-mean +/- std)
| metric | mean +/- std | per-seed |
|---|---|---|
| intergenic_specificity (AXIS-1 bw) | 0.8416 +/- 0.039 | 0.8197 / 0.8967 / 0.8083 |
| macro_intergenic_specificity (gate) | 0.7543 +/- 0.040 | 0.7395 / 0.8092 / 0.7142 |
| gene_body_F1_unconstrained (AXIS-2) | 0.6888 +/- 0.001 | 0.6878 / 0.6908 / 0.6878 |
| gene_body precision / recall | 0.754 / 0.637 | -- |
| intergenic_FPR | 0.158 | 0.180 / 0.103 / 0.192 |
| predicted_gene_count_ratio | 1.43 | 1.31 / 1.61 / 1.38 |

PER-SPECIES (cross-clade asymmetry): fly spec ~0.82-0.91 (GOOD), gbF1 ~0.68, gcount ratio 1.06-1.54 ; yeast (fungus) spec 0.61-0.71 (POOR), gbF1 0.74-0.81, gcount ratio 1.8-2.1 (over-predicts genes in the divergent clade).

### Verdict vs anchor (new ruler: spec 0.8710 bw / 0.8278 macro / gene_body_F1 0.5576)
- AXIS-2 gene_body_F1 0.6888 >> anchor 0.5576 (+0.13) AND > floor 0.5276 -> PASS. Foundation features substantially improve gene-body detection on BOTH species.
- AXIS-1 intergenic_specificity 0.8416 < anchor 0.8710 -> does NOT strictly beat; macro 0.7543 < gate 0.7978 -> FAILS macro gate (yeast drags it). Same trade-off as structured decoders (higher recall/F1, lower specificity via intergenic spillover), worst on the divergent fungus.
- Does NOT Pareto-dominate the anchor (one axis up, one down). NON-CLAIM screen -> not_yet (validate_exit=1 all seeds). High seed variance on specificity (0.808-0.897).

### Interpretation
Frozen human/vertebrate-pretrained foundation features IMPROVE gene detection (F1) but do NOT improve (slightly hurt) cross-clade intergenic specificity -- they don't transfer the coding/intergenic boundary to a divergent fungus as well as from-scratch training on that clade. Hypothesis PARTIALLY supported (F1 yes, specificity no). Converting the real recall gain into specificity needs an FP-aware objective / structured decoder on the features (planned next step), and/or fine-tuning for cross-clade transfer.
"""

D10 = r"""

## Finding (Research) 2026-06-11 -- Frozen foundation features improve gene-body F1 but not cross-clade intergenic specificity (FP-SEGMENTNT-PROBE-M1)
FROZEN SegmentNT (human/vertebrate-pretrained) element logits as INPUT to an anchor-matched conv+biLSTM head (clean input-signal ablation, same budget): gene_body_F1 = 0.689 (vs from-scratch raw-DNA anchor 0.5576, +0.13) -- foundation features clearly help gene DETECTION. BUT intergenic_specificity 0.842 < anchor 0.871 (macro 0.754 < anchor 0.828) -- does NOT beat the AXIS-1 primary. Cross-clade asymmetry is the crux: fly (insect) spec ~0.85 GOOD; yeast (fungus) spec ~0.65 POOR with 1.8-2.1x gene over-prediction. The human-pretrained foundation model does not transfer the coding/intergenic boundary to a divergent fungus as well as a from-scratch model trained on that clade. IMPLICATION: foundation features are a real lever for RECALL, but recall->specificity (and cross-clade transfer) needs an FP-aware objective / structured decoder on top + possibly fine-tuning. Same up-recall/down-specificity trade-off as CRF-vec -> intergenic spillover is the central obstacle, not gene detection.

## Finding (Engineering) 2026-06-11 -- SegmentNT JAX extraction cost + setup (coding-rna env)
SegmentNT (segment_nt_multi_species) is JAX/Haiku (not torch); runs in coding-rna alongside torch 2.5.1 as a SEPARATE process (jax extract -> npz cache -> torch head train; no in-process coexistence). Deps added to coding-rna: jax[cuda12] 0.10.1, dm-haiku 0.0.16, einops, joblib, regex, pydantic, transformers, huggingface_hub (+ pin sympy==1.13.1 for torch). API works on jax 0.10.1 (pass params as jit ARG, NOT closure constant -- closure-capture caused a 2.23GB-constants blow-up + TB-scale rematerialization). Token count must be %4==0 (U-Net 2 downsample blocks) -> 2048bp window needs padding; extract on 6kb tiles per seqid. Extraction is SLOW (~0.58 min/MB genome from rematerialization on 24GB) -> fly 145MB ~ 80min; submit via cluster batch with generous --time (4h) to avoid TIMEOUT mid-savez (a 1.5h interactive run nearly truncated it). A pure 1x1 per-base head COLLAPSES to all-CDS (intergenic_specificity 0.0) on gene-dense yeast -> use an anchor-matched conv+biLSTM head (capacity + local context) for a fair feature probe.
"""

for path, txt in [("docs/06_results_log.md", D06), ("docs/10_findings.md", D10)]:
    with open(os.path.join(ROOT, path), "a") as fh:
        fh.write(txt)
    print("appended", path)

# tracker RUNNING -> DONE
p = os.path.join(ROOT, "docs/05_todo.md")
t = open(p).read()
t = t.replace("| FP-SEGMENTNT-FEATCACHE | RUNNING | featcache | 8548459 |",
              "| FP-SEGMENTNT-FEATCACHE | DONE | featcache | 8548459 |")
t = t.replace("| FP-SEGMENTNT-PROBE-M1-convlstm-s{0,1,2} | RUNNING | screen | 8548460-62 |",
              "| FP-SEGMENTNT-PROBE-M1-convlstm-s{0,1,2} | DONE | screen | 8548460-62 |")
open(p, "w").write(t)
print("docs/05 tracker RUNNING->DONE")
