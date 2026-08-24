---
exp_id: M23-NTV2-CLEAN-TRANSFER-s0
date: 2026-07-01
approach_family: pretrained-LM-clean-provenance
parent_exp: M10-M9L12-CLEANPLANTS-s0
motivated_by: "M22 negative; user requested clean-provenance NT-v2 transfer single-seed screen without M22 gb_tversky/CRF/calibration"
track: B
profile: screen
status: done
primary_metric: intergenic_specificity
value: 0.9833
vs_anchor: "screen-pass at FPR<=0.02, hard FPR<=0.01 fails; exactly reproduces M10 s0"
one_liner: "Clean-provenance NT-v2 direct transfer is valid but not frontier-moving"
---

## Why / Motivation

M22 refuted the GENERanno `gb_tversky` objective and the user explicitly stopped M22 `gb_tversky`/CRF/calibration continuation. M23 checks the parallel claim route: whether the clean public-provenance NT-v2 backbone, adapted with our known 3-class FP-aware recipe, can offer a better paper-facing route than the stronger but provenance-blocked GENERanno adaptation.

## Hypothesis

If direct NT-v2 transfer is the right clean-provenance route, this single seed should approach the useful M19 tradeoff while avoiding GENERanno's overlap blocker. If it simply reproduces M10, then direct NT-v2 remains a clean baseline but not a new claim route.

## Architecture

`InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`, top-12 backbone layers unfrozen, 3-class convLSTM head, `fp_aware` loss, constrained post-processing (`min_cds_len=60`, `max_fill_gap=20`). No raw-score saving, no validation-only calibration, no CRF decoder, and no M22 `gb_tversky`.

## Data

Clean plant panel `{arabidopsis_thaliana, oryza_sativa}` from `data/m1_screen/`, using deterministic chromosome/seqid train/val/test splits. NT-v2 public provenance remains cleaner for plants because public model/dataset cards exclude plants.

## Config

Config path: `configs/M23-NTV2-CLEAN-TRANSFER.yaml`. Sbatch: `sbatch/M23-NTV2-CLEAN-TRANSFER-s0.sbatch`. Seed `0`; window `2046`; sample_fraction `0.3`; epochs `4`; batch size `4`; head lr `1e-3`; backbone lr `1e-5`; class weighting `sqrt_inv`.

## Result

Semantic success passed. Slurm job `9854668` completed `0:0` on `gpu034` in `19:25:54`; metrics are finite; `STATUS=COMPLETED`; no OOM/traceback. Loss decreased `0.7453 -> 0.6161 -> 0.5487 -> 0.4823`; best validation macro-F1 was `0.8213`.

Aggregate TEST metrics:
- gbF1 `0.8427`
- FPR `0.01673`
- specificity `0.98327`
- macro specificity `0.98050`
- gene_count_ratio `0.867`
- constrained gbF1@0.01 `0.0`

Per species:
- Arabidopsis: gbF1 `0.8983`, FPR `0.02520`, gene_count_ratio `0.785`.
- Rice: gbF1 `0.7606`, FPR `0.01380`, gene_count_ratio `1.033`.

## Findings

M23 exactly reproduces historical `M10-M9L12-CLEANPLANTS-s0`, so it is a clean-provenance route checkpoint rather than a new performance improvement. Direct NT-v2 transfer remains coherent and screen-valid at FPR<=0.02, but it fails hard FPR<=0.01 and is weaker than M19 GENERanno calibrated seed1 (`gbF1=0.8815`, FPR `0.0065`).

## Decision

Do not rerun more direct M10/M23-style NT-v2 seeds. If clean-provenance NT-v2 remains the desired claim route, the next GPU must make a structural change or serve a claim-panel/comparability purpose. Recommended next framework step: combined M22+M23 `$tri-review`/`$pivot`.

## Links

- result-log: `docs/06_results_log.md#result-m23-ntv2-clean-transfer-s0`
- iteration: `docs/04_experiment_iterations.md#iter-m23-ntv2-clean-transfer-s0--2026-07-01`
- output: `outputs/M23-NTV2-CLEAN-TRANSFER-s0`
- log: `outputs/fp_segnt_logs/M23NTV2S0_9854668.out`
