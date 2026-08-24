---
exp_id: M11-L12-SPEC-CALIBRATION
date: 2026-06-16
approach_family: pretrained-LM
parent_exp: M10-M9L12-CLEANPLANTS
motivated_by: "M10 pivot: M9-L12 strong gbF1/coherence but aggregate FPR 0.0174 > full/scale 0.01"
track: B
profile: screen
status: done
primary_metric: intergenic_specificity
value: 0.9913
vs_anchor: "+0.1203 vs screen_anchor 0.8710; FPR 0.0087 clears 0.01"
one_liner: "M9-L12 validation-only decode calibration clears FPR<=0.01 without stronger FP objective"
---

## Why / Motivation

M10 established NT-v2 top-12 unfreeze as the mainline on clean plants, but it remained blocked by the hard FPR `<=0.01` claim guardrail. M11 tested whether this was an operating-point problem rather than a need for a new backbone or stronger objective.

## Hypothesis

Saved raw emissions from the same M9-L12 model should allow no-leakage validation-only decode calibration to lower intergenic FPR while preserving gbF1 and gene-count coherence.

## Architecture

Same as M10: NT-v2-500m ESM backbone with top-12 layers trainable, 3-class intron-aware convLSTM head, FP-aware loss (`fp_lambda=1.0`). Added raw VAL/TEST score saving and an offline calibrator that sweeps intergenic logit bias plus constrained-decode parameters on validation seqids only.

## Data

Clean held-out plants `{arabidopsis_thaliana, oryza_sativa}` from `data/m1_screen/`, deterministic seqid/chromosome-aware split. Per seed: arabidopsis train=5 / val=1 / test=1; rice train=6 / val=1 / test=1.

## Config

Train config: `window=2046`, `sample_fraction=0.3`, `epochs=4`, `batch_size=4`, `lr=1e-3`, `backbone_lr=1e-5`, seeds `0/1/2`.

Calibration grid: intergenic bias `{0,0.5,...,4.0}`, `min_cds_len={60,90,120}`, `max_fill_gap={0,20}`. Selection requires validation `FPR<=0.01`, `gbF1>=0.70`, and `gene_count<=1.25`, then applies once to TEST.

## Result

All three seeds completed successfully. Seed mean: intergenic_specificity `0.9913`, FPR `0.0087`, macro_specificity `0.9909`, gbF1 `0.8178`, constrained_gbF1@0.01 `0.8178`, gene_count_ratio `1.003`. Selected points were `b2p5_mcl60_mfg20`, `b3p0_mcl60_mfg20`, and `b1p5_mcl60_mfg0`.

## Findings

Validation-only calibration fixes M10's FPR tail (`0.0174 -> 0.0087`) without collapsing gbF1 or gene count. Stronger FP objective is not the immediate next move. Per-species caveat: arabidopsis seed2 FPR is `0.0111`, so full/scale needs per-species sensitivity and more validation chromosomes.

## Decision

Pending tri-review/pivot. Likely direction: promote calibrated M9-L12 to full/scale/comparability work; defer stronger FP objective unless full/scale reopens the FPR blocker.

## Links

- result-log: docs/06_results_log.md#result-m11-l12-spec-calibration
- iteration: docs/04_experiment_iterations.md#iter-m11-l12-spec-calibration--2026-06-16
- sbatch: sbatch/M11-L12-SPEC-CALIBRATION.sbatch
- outputs: outputs/M11-L12-SPEC-CALIBRATION-s{0,1,2}
- logs: outputs/fp_segnt_logs/M11L12CAL_8934130_{0,1,2}.out
