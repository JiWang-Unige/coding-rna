# Validation Matrix / 下游任务与可靠性论证矩阵

> 由 `/publication-plan`、`/generalization`、`/sota-randomized` 共同维护。用于回答“模型/流程为什么可靠，哪些下游任务必须做”。

## 1. Main result
| Candidate | Dataset/split | Metric | Value | Comparator | Comparable? | Evidence path |
|---|---|---|---:|---|---|---|
| M11 calibrated M9-L12 | clean plants screen `{arabidopsis,rice}`, same panel used for train/val/test partitions | intergenic_specificity / FPR / gbF1 / gene_count | 0.9913 / 0.0087 / 0.8178 / 1.003 | M10-M9L12, same route before calibration | Partially; screen non-claim only | `outputs/M11-L12-SPEC-CALIBRATION-s{0,1,2}/metrics/metrics.json` |
| M12A fixed M9-L12 | frozen train/calibration species -> unseen species/clade test on `configs/m12_publication_panel.yaml` | same primary + per-species FPR sensitivity | TBD | Tiberius/Helixer/ANNEVO same-panel where possible | No, ready not submitted | `M12A-FIXEDMODEL-CROSSSPECIES` config/sbatch/gate ready |
| M12B same-panel baselines | same clean plant panel and evaluator for Tiberius/Helixer/ANNEVO | Pareto: specificity, FPR, gbF1, gene_count | Tiberius `0.9927/0.0073/0.9252/0.628`; ANNEVO `0.9883/0.0117/0.9269/0.726`; Helixer `0.9784/0.0216/0.9220/0.820` | direct external tools | Yes for clean-plant utility table | `outputs/M12B-SAMEPANEL-BASELINES-*`; `reports/M19-COMPARABILITY-EVIDENCE/` |
| M18/M19 GENERanno 1.2B challenger | same clean plant panel/head/evaluator; 1.2B CDS-preview + 0.5B base control | same primary + seed stability + raw-score calibration | M18 1.2B: `0.9929/0.0071/0.8494/0.864`; M19 calibrated s0 `0.9917/0.0083/0.8421/1.083`, s1 `0.9935/0.0065/0.8815/0.830`; M18 0.5B negative | Tiberius/ANNEVO/Helixer clean-plant rows | Mechanism/adaptation comparable; not clean held-out claim because provenance unknown and released callers remain higher-gbF1 | `outputs/M18-GENERANNO-*`; `outputs/M19-GENERANNO-*`; `reports/M19-COMPARABILITY-EVIDENCE/` |
| M21 GENERanno 1.2B CRF decoder | same clean plant panel/head/evaluator; M19 route plus trained CRF decoder | same primary + FPR sensitivity + gene count | seed0 `spec/FPR/gbF1/gcount=0.9727/0.0273/0.8544/0.956`; seed1 rescue `0.9808/0.0192/0.8744/0.690` | M19 non-CRF s1 `0.9935/0.0065/0.8815/0.830`; Tiberius `0.9927/0.0073/0.9252/0.628` | Negative ablation; CRF worsens FPR and does not beat M19 | `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s0`; `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt`; `docs/09_decisions_log.md#dec-001-generanno-12b-trained-crf-decoder-route` |

## 2. Downstream / external validation tasks
| Task ID | Task | Purpose for paper claim | Dataset | Metric | Baselines | Required seeds | Status | Output path |
|---|---|---|---|---|---|---:|---|---|
| D1 | Fixed-model cross-species evaluation | Prove deployable model behavior beyond same-species train/val/test pool. | Clean same-panel species, with held-out species/clade separated before training/calibration. | intergenic_specificity, macro_specificity, gbF1, constrained gbF1, gene_count ratio, FPR thresholds 0.005/0.01/0.02. | M9-L12 frozen; Tiberius/Helixer/ANNEVO as contextual comparators. | 3 seeds prepared for A2R. | READY_BLOCKED_BY_SLURM | `outputs/M12A-*` |
| D2 | Same-panel external gene-caller comparison | Make practical utility legible to biology users and reviewers. | Exact clean plant panel. | Same evaluator and span rule; FPR sensitivity; gene-count ratio; fixed-vs-adapted label. | Tiberius, Helixer, ANNEVO, GENERanno 1.2B adapted route. | released tools: 1; adapted GENERanno: 2 M19 seeds. | READY_REPORT_WRITTEN | `reports/M19-COMPARABILITY-EVIDENCE/`; `outputs/M12B-*`; `outputs/M18-*`; `outputs/M19-*` |
| D3 | GENERanno fair challenger | Test whether strong pretrained annotation backbones alone solve the task or whether CDS-specialized backbone + our intron-aware/FP-aware adaptation is doing specific work. | Exact clean plant panel; GENERanno overlap status remains `unknown`. | same primary metrics + FPR/gene-count diagnostics + raw-score calibration. | `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`, `GenerTeam/GENERanno-eukaryote-0.5b-base`, Tiberius/ANNEVO/Helixer context. | 2 M19 seeds for 1.2B; 0.5B single-seed ablation done. | DONE_NONCLAIM | `outputs/M18-GENERANNO-*`; `outputs/M19-GENERANNO-*` |

## 3. Robustness / OOD / sensitivity
| Check ID | Perturbation / OOD split | Expected behavior | Pass criterion | Status | Result |
|---|---|---|---|---|---|
| R1 | NT-v2 and GENERanno pretraining/species-overlap audit | No hidden overlap that invalidates clean held-out claim. | Species/accessions categorized as clean, contaminated, or unusable before GPU. | BLOCKER_OPEN for GENERanno: public sources do not expose full species/accession exclusion list | `refs/dossiers/m12_prereq_audit.md`; `refs/dossiers/m19_generanno_provenance_audit.md` |
| R2 | Per-species FPR sensitivity | Aggregate specificity should not hide one failing species. | Report FPR at thresholds 0.005/0.01/0.02 and macro specificity; full/scale hard guardrail remains FPR<=0.01 where applicable. | TODO | M12A/M12B |
| R3 | Held-out species/clade split | Model should generalize without species-specific test calibration. | Frozen checkpoint/calibration; no test-label parameter selection. | TODO | M12A |
| R4 | Annotation-density/chromosome strata | Performance should not be driven only by easy contigs or gene-dense regions. | Stratified metrics reported or justified as not feasible. | TODO | after M12A panel freeze |

## 4. Ablations
| Ablation ID | Removed/changed component | Hypothesis | Metric | Expected delta | Status |
|---|---|---|---|---:|---|
| A1 | M10 uncalibrated vs M11 validation-only calibration | Calibration lowers FPR below 0.01 without destroying gene-body recovery. | FPR, gbF1, gene_count | FPR down; gbF1 acceptable | DONE screen |
| A2 | GENERanno 1.2B CDS-preview vs 0.5B base | If generic eukaryotic pretraining were sufficient, 0.5B base should approach the 1.2B CDS-preview under the same objective/head. | specificity, FPR, gbF1, gene_count | 1.2B far stronger and stable across M19 seeds; 0.5B remains FP/coherence negative | DONE screen |
| A3 | Stronger FP objective | Only needed if fixed-model/full-panel evidence reintroduces FPR failure. | FPR at 0.01, gbF1 | FPR down with bounded gbF1 loss | DEFERRED |
| A4 | GENERanno trained CRF decoder | If M19 emissions are low-FPR but recall-limited, a CRF decoder should recover genes without breaking FPR. | gbF1, FPR, gene_count | NEGATIVE: best CRF seed `gbF1=0.8744/FPR=0.0192`, worse than M19 s1 `0.8815/0.0065`; route abandoned | DONE negative |
| A5 | Fixed calibration vs per-species calibration | Per-species tuning would inflate apparent generalization. | delta in specificity/gbF1 | fixed calibration remains competitive | TODO if reviewers require |

## 5. Randomized SOTA small-sample retraining
> `init` 列必填：random=随机初始化重训（同预算公平参考，建 screen_anchor）/ pretrained=载入官方权重（仅作上界对照，**不混入随机初始化的 mean**）。split_scheme/metric_impl 必须与我们的 Track A 一致（防泄漏可比）。
| Model | init | sample_fraction | seeds | split_scheme | Metric mean±std | Our comparable run | Verdict | Link |
|---|---|---:|---|---|---|---|---|---|
| tiberius_like | random | screen | 3 | M1 same-budget | gbF1 seed-mean 0.5576 | screen anchor family | coherent anchor | `outputs/SCREENREF-tiberius_like-s{0,1,2}` |
| helixer_like | random | screen | 3 | M1 same-budget | gbF1 seed-mean 0.5579 but fragmented | excluded as fragmented anchor | contextual | `outputs/SCREENREF-helixer_like-s{0,1,2}` |
| Tiberius/Helixer/ANNEVO official | pretrained | TBD | 1 | M12 same-panel pending | TBD | M12B | pending | `outputs/M12B-*` |

## 6. Statistical tests
| Comparison | Test | Paired? | n/seeds | p-value/CI | Status |
|---|---|---|---:|---|---|
| M11 vs M10 calibration | paired by seed/species | yes | 3 seeds | TBD if needed | screen evidence only |
| M12A M9-L12 vs external tools | paired bootstrap by species/chromosome or fixed panel units | yes | panel-dependent | TODO | pending panel freeze |
| M19 GENERanno 1.2B seed stability | paired by clean plant panel/species | yes | 2 seeds | calibrated gbF1 range `0.8421-0.8815`; FPR range `0.0065-0.0083`; gene_count ratio range `0.830-1.083` | DONE screen/non-claim |
