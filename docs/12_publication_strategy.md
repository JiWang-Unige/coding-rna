# Publication Strategy / 投稿推进计划

> 由 `/publication-plan` 维护。用于“已经有完整思路或已超越 SOTA，不是盲目迭代，而是把研究做成可投稿故事”的阶段。

## 0. Target positioning
- Target venue / journal tier: computational genomics / genome annotation method paper tier; exact venue pending after same-panel baseline evidence.
- Backup venues: benchmark+method or resource+pipeline venues if full SOTA claim remains blocked by published-anchor comparability.
- Article type: method paper with benchmark/evaluator component.
- Expected novelty bar: fixed-model cross-species ab initio gene annotation that improves practical false-positive control while preserving gene-level recovery.
- Audience: computational genomics method developers and biology users who need deployable genome annotation with fewer spurious genes.

## 1. Core story
- One-sentence paper claim (draft): A calibrated NT-v2-based intron-aware gene annotator achieves lower intergenic false-positive burden and robust gene-body recovery across held-out eukaryotic species than existing gene callers under a shared evaluation panel.
- Why now: M11 shows the lead M9-L12 route can hit strict intergenic FPR on clean plant screen data, but the paper value now depends on comparability and fixed-model generalization rather than more screen-scale metric polishing.
- Why existing work is insufficient: Tiberius, Helixer, ANNEVO, and pretrained annotation models are not yet compared under one raw-DNA panel, one span/evaluator rule, and one practical FP/gene-count utility view.
- Our key insight: pretrained DNA emissions plus an intron-aware 3-class head and validation-only FP calibration can make gene-body regions coherent without inflating intergenic predictions; M12 must prove this is not merely a generic pretrained-model effect.
- 2026-06-20 update: the current strongest route is no longer a broad fixed NT-v2 claim. M18/M19 make `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` a stable adapted pretrained-CDS backbone: calibrated M19 seeds reach gbF1/FPR/gene_count `0.8421/0.0083/1.083` and `0.8815/0.0065/0.830`. Public provenance still does not clear Arabidopsis/rice overlap, and released clean-plant callers still have higher gbF1, so write GENERanno results as adaptation/comparability evidence, not clean no-overlap held-out SOTA.

## 2. Contribution menu

| Contribution ID | Claim | Evidence needed | Current evidence | Risk | Keep? |
|---|---|---|---|---|---|
| C1 | Fixed calibrated M9-L12 works as a cross-species annotation model, not just a same-species screen model. | `M12A-FIXEDMODEL-CROSSSPECIES`: train/calibrate once, freeze checkpoint/calibration, evaluate unseen species/clades. | M11 clean plants screen: spec 0.9913 / FPR 0.0087 / gbF1 0.8178, but not fixed unseen-species proof. | Pretraining overlap or validation/test leakage could invalidate claim. | yes |
| C2 | Our method is more practical than Tiberius/Helixer/ANNEVO under the same panel because it reduces intergenic FP/gene-count inflation while retaining gene-body recovery. | Same evaluator/same panel tables with CDS/full-transcript span contract, FPR sensitivity, and gene-count reporting. | `reports/M19-COMPARABILITY-EVIDENCE`: clean-plant released callers have higher gbF1 (`~0.922-0.927`); M19 GENERanno 1.2B is FPR-valid after calibration (`0.0083`, `0.0065`) and gene-count sane (`1.083`, `0.830`), while Tiberius is the closest practical comparator but under-calls genes (`0.628`). | Current GENERanno gbF1 does not beat released clean-plant callers; story must be utility/adaptation/tradeoff, not blanket superiority. | yes, reframed |
| C3 | The gain is not just “any pretrained model works”; CDS-specialized pretraining plus our intron-aware/FP-aware adaptation matters. | Compare GENERanno 1.2B CDS-preview and 0.5B base under the same 3-class head/evaluator/objective; test seed stability and raw-score calibration. | M18 1.2B CDS-preview is strong (`FPR=0.0071`, gbF1 `0.8494`, gene_count `0.864`); M18 0.5B base is negative (`FPR=0.0967`, gbF1 `0.6561`, gene_count `1.617`). M19 1.2B is stable and calibratable across two seeds. | GENERanno overlap remains unknown; claim-clean species or explicit adaptation framing is required. | yes |
| C4 | Validation-only FP calibration is a practical deployment knob, not post-hoc test tuning. | M11/M12A sensitivity at FPR thresholds 0.005/0.01/0.02; checkpoint + calibration artifact contract. | M11 selects on validation and clears aggregate FPR<=0.01. | Small validation diversity may overfit. | yes |

## 3. Figure / table plan

| Fig/Table | Message | Required experiments / analyses | Owner docs | Status |
|---|---|---|---|---|
| Fig.1 | Task, split, evaluator, and practical utility metrics. | Frozen same-panel species list, span contract, metric definitions. | `docs/19_evaluator_contract.md`, `docs/14_validation_matrix.md` | TODO |
| Fig.2 | Pareto comparison: intergenic specificity vs gene-level F1 / gene-count ratio. | Clean-plant same-panel Tiberius, Helixer, ANNEVO, M18/M19 GENERanno 1.2B, and 0.5B base control. | `reports/M19-COMPARABILITY-EVIDENCE/comparison_tables.md`, `docs/14_validation_matrix.md`, `docs/06_results_log.md` | READY_NONCLAIM |
| Fig.3 | Fixed-model cross-species generalization. | M12A per-species and aggregate metrics for unseen species/clades. | `docs/14_validation_matrix.md` | TODO |
| Fig.4 | Mechanism/ablation: calibration and backbone challenge. | M10/M11 calibration comparison; M12C GENERanno fair challenger; optional stronger FP objective only if needed. | `docs/14_validation_matrix.md`, `docs/10_findings.md` | TODO |
| Table 1 | Main benchmark numbers and utility guardrails. | M12A/M12B final metrics with CI/seed handling. | `docs/06_results_log.md` | TODO |
| Table 2 | Runtime/resource/deployability. | Inference wall time, GPU/CPU requirements, tool versions, checkpoints. | `docs/16_artifact_registry.md` | TODO |

## 4. Validation burden by venue tier

| Evidence type | Minimal | Strong | Needed for target? | Planned run/analysis |
|---|---|---|---|---|
| Main benchmark | Clean same-panel comparison on at least two held-out species. | Multi-clade same-panel with ANNEVO-compatible published anchor freeze. | yes | M12B + benchmark contract |
| Downstream/generalization | Fixed checkpoint tested on unseen species from the same broad clade. | Held-out clade / broader eukaryote panel with no per-species tuning. | yes | M12A |
| Ablation | M10/M11 calibration and GENERanno fair challenger/base-vs-CDS control. | Backbone/head/objective/calibration ablations plus failure-mode analysis. | yes | M11, M18, M19 |
| Robustness/OOD | Per-species FPR/gene-count sensitivity at 0.005/0.01/0.02. | Macro specificity, chromosome/contig strata, annotation-density strata. | yes | M12A/M12B |
| Runtime/cost | Tool versions and wall time on same panel. | Reproducible containers/envs + checkpoint release contract. | yes | Artifact registry + M12B |
| Statistical test | Multi-seed only where training stochasticity matters. | Paired species/chromosome bootstrap + seed CIs. | likely | After fixed panel is frozen |

## 5. Rebuttal risk pre-mortem

| Likely reviewer criticism | Evidence to preempt | Where captured |
|---|---|---|
| “This only works because train/val/test come from the same species pool.” | M12A fixed-model unseen-species/clade evaluation with frozen calibration. | `docs/14_validation_matrix.md` |
| “Numbers are not comparable to Tiberius/Helixer/ANNEVO.” | M12B same-panel baseline contract and raw command/version artifacts. | `docs/20_baseline_reproduction.md`, `docs/14_validation_matrix.md` |
| “Pretraining overlap contaminates the held-out species.” | M12/M19 provenance audits; current rule: GENERanno Arabidopsis/rice rows are adaptation/challenger evidence unless training manifest clears overlap. | `refs/dossiers/m19_generanno_provenance_audit.md`, `docs/15_evidence_register.md` |
| “GENERanno or any pretrained model would do the same.” | M18 1.2B CDS-preview vs 0.5B base under the same stronger FP objective; M19 seed stability and raw-score calibration. | `docs/14_validation_matrix.md`, `reports/M19-COMPARABILITY-EVIDENCE/` |
| “FP calibration is test-set tuning.” | Persist validation-only selection rule, selected parameters, and held-out test application. | `docs/06_results_log.md`, run artifacts |

## 6. Manuscript readiness checklist
- [ ] Main claim has comparable full/scale result.
- [ ] Fixed-model cross-species protocol is frozen and run.
- [x] Tiberius/Helixer/ANNEVO same-panel baseline artifacts are reproduced or explicitly justified for clean plants and summarized in a paper-facing report.
- [x] NT-v2 and GENERanno pretraining/species overlap audit v0.1 is complete for M12 planning: NT-v2 plant panel is clean; GENERanno remains overlap-unknown/mechanism-only.
- [x] GENERanno fair challenger seed-stability/calibration is complete. M18 1.2B vs 0.5B is complete; M19 two-seed raw-score cohort and VAL-only calibration are complete.
- [ ] Data/split/metric provenance archived in `refs/dossiers/`.
- [ ] Figure/table evidence mapped to run IDs.
