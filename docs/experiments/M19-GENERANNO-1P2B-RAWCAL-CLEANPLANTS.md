---
exp_id: M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS
date: 2026-06-20
approach_family: pretrained-CDS-adaptation
parent_exp: M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0
motivated_by: "docs/08 pivot M17+M18: switch next primary challenger to GENERanno 1.2B CDS-preview raw-score/two-seed preflight"
track: B
profile: screen
status: done
primary_metric: intergenic_specificity
value: 0.9935
vs_anchor: "non-claim screen; best calibrated seed s1 FPR 0.0065 / gbF1 0.8815; released clean-plant callers still higher gbF1"
one_liner: "GENERanno 1.2B CDS-preview is stable across two seeds and calibratable, but remains adaptation/comparability evidence"
---

## Why / Motivation

M18 showed a single strong GENERanno 1.2B CDS-preview result under stronger FP objective. M19 tests whether that was seed luck and whether saved raw scores can support validation-only FPR calibration without test-label tuning. It also feeds the paper-facing same-evaluator comparison table against ANNEVO, Tiberius, and Helixer.

## Hypothesis

If the M18 1.2B result is real, two M19 seeds should both keep aggregate FPR<=0.01, sane gene counts, and non-collapsed gbF1. If raw-score calibration helps, it should improve gbF1/FPR/gene-count tradeoff on TEST after selecting only on VAL.

## Architecture

`GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` token-classification backbone, k=6, LoRA r=8 on q/k/v/o projections, our 3-class FP-aware convLSTM head, constrained decode, `fp_lambda=2.5`, `min_cds_len=90`, `max_fill_gap=20`. M19 adds `--save-raw-scores` and two seeds.

## Data

Clean plant panel `{arabidopsis_thaliana, oryza_sativa}` with chromosome/seqid-aware train/val/test splits inherited from M18/M12. GENERanno provenance for these species remains public-overlap-unknown, so this is non-claim mechanism/comparability evidence.

## Config

Config path: `configs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS.yaml`. Slurm array: `9141356_[0-1%2]`. Seeds: `0/1`. Training window caps: `train_windows=1536`, `val_windows=768`; batch size `1`; LoRA lr `2e-5`; head lr `8e-4`.

## Result

Semantic success passed. Both array tasks completed `0:0`; metrics are finite; raw score files exist for VAL/TEST across both species and seeds; no OOM/Traceback.

Raw decode:
- s0: gbF1 `0.8390`, FPR `0.0088`, specificity `0.9912`, gene_count_ratio `0.967`.
- s1: gbF1 `0.8593`, FPR `0.0059`, specificity `0.9941`, gene_count_ratio `0.805`.

Validation-only calibrated decode:
- s0 selected `b2p0_mcl60_mfg20`: gbF1 `0.8421`, FPR `0.0083`, specificity `0.9917`, gene_count_ratio `1.083`.
- s1 selected `b0p0_mcl60_mfg20`: gbF1 `0.8815`, FPR `0.0065`, specificity `0.9935`, gene_count_ratio `0.830`.

## Findings

M19 confirms the M18 1.2B positive was not a random seed. Calibration is useful, especially for seed1, but the route still trails released clean-plant callers on gbF1. The strongest paper-safe claim is adaptation/comparability: CDS-specialized GENERanno plus our 3-class/FP-aware adaptation reaches a low-FP coherent operating regime, while generic 0.5B base does not.

## Decision

Pending tri-review/pivot. Candidate next actions are: cleaner held-out species panel for claim hygiene, a segment/structured head to close the gbF1 gap, or freezing GENERanno as adaptation evidence while shifting the claim route.

## Links

- result-log: `docs/06_results_log.md#result-m19-generanno-1p2b-rawcal-cleanplants-s01`
- iteration: `docs/04_experiment_iterations.md#iter-m19-generanno-1p2b-rawcal-cleanplants--2026-06-19`
- report: `reports/M19-COMPARABILITY-EVIDENCE/comparison_tables.md`
- runs: `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0`, `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1`
