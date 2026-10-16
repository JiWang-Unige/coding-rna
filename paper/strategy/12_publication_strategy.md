# Publication Strategy / 投稿推进计划

目标是 Nature Communications 及以上档次的原始方法论文：提出从 raw FASTA 直接输出 gene 分组、完整 CDS 链、strand 和 phase 的 ab initio 注释模型，以更好的完整链性能或经独立验证的更广适用范围为核心。当前仍在研发，没有支持投稿级优越性或独立泛化的模型结论。

## 2026年10月9日 M28 全语法模型：当前证据与投稿主线

当前路线已从有限候选集合训练推进到全语法模型及严格 raw-FASTA 完整读出。工程接入成立不等于完整链性能领先；目标仍是 Nature Communications 及以上的原始方法论文，不能以负结果记录或流程复现代替该目标。下方10月7日及更早段落保留为阶段历史，其中的作业状态和授权边界不作为当前依据。

### 已完成的原 TRAIN 同尺对标

以下仅比较原两物种、六条完整 TRAIN 染色体，保留156779632 bp输入及22688条 primary-reference CDS链（拟南芥17739、水稻4949），不按候选可达性、编码长度或非标准结构删难例。主终点为同strand、相同有序CDS坐标的 exact-chain F1；phase／结构合法性另查，不等于完整UTR或全部isoform恢复。表中F1均为百分数，macro为两个物种F1的算术均值，worst为较低物种值；pooled按合并TP／预测／参考计数计算，不是物种均值。

| 固定模型／权重 | 拟南芥 | 水稻 | pooled | macro-species | worst-species |
| --- | ---: | ---: | ---: | ---: | ---: |
| M28 原始末4608步 | 50.487 | 28.278 | 45.590 | 39.383 | 28.278 |
| M28 LOCAL末9216步 | 44.160 | 33.382 | 41.729 | 38.771 | 33.382 |
| ANNEVO 2.3.3 / Magnoliopsida | 76.831 | 71.179 | 75.542 | 74.005 | 71.179 |
| Tiberius 2.0.8 / 固定植物资产 | 75.417 | 63.577 | 72.655 | 69.497 | 63.577 |
| Helixer 0.3.7 / land_plant_v0.3_a_0080 | 68.122 | 50.294 | 63.592 | 59.208 | 50.294 |

来源为已完成评价13443329／13451623／13461307及其[原4608步结果](../../experiments/M28/full_grammar/reports/full_grammar_raw_train_result.md)、[原生TRAIN对标](../../benchmarks/native_baselines/reports/native_train_baseline_result.md)和[LOCAL恢复诊断](../../experiments/M28/full_grammar/reports/full_grammar_train_recovery_local_result.md)。本次写入前由TP／预测／参考重算逐物种与pooled F1，并核对macro／worst及共同分母。完整读出只给模型原FASTA，不将参考链或训练缓存用作推理输入；上述评价没有DEV、封存访问或参数更新。TRAIN适应效果、上游资产暴露和独立泛化资格必须分开；这些数值也不能与不同面板的旧A31 DEV分数直接比较来判断性能回退。

两个已完成M28权重均未领先三个原生基线。LOCAL相对原4608步恢复2119条旧遗漏、丢失2542条旧命中，净TP减少423；pooled CDS-base F1从90.360%升至94.097%，但完整链F1下降，纯背景预测链从162增至491。水稻与worst-species改善不等于整体收益；纯背景按注释范围定义，不自动等同生物学假阳性。原4608步少4608次更新，尚不能将这组变化因果归于local目标，须等待匹配更新数的STRUCTURED_CONTINUE结果。

### 机制证据和创新性边界

[Tiberius正式论文](https://academic.oup.com/bioinformatics/article/40/12/btae685/7903281)已将可微HMM接入端到端训练；其24个标准HMM参数并未训练。这个实现差异不足以单独构成新颖性或收益证据；GLM／神经表示加结构解码、结构约束和局部／结构联合学习均有先例，见[现有先例核对](../manuscript/M28_complete_chain_competition_prior_work_20261003.md)。

2026年10月8日[Carbon-A官方发布](https://huggingface.co/blog/HuggingFaceBio/carbon-annotator-genbank-genome-annotation)及作者[技术报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)进一步限制宽泛创新主张：作者披露的完整流程已有transition TCN、Viterbi基因语法和confidence TCN，不能仅据概率输出示例将其排除为“只做分割”。本项目固定GENERanno b0483c23b6b63787b61a6d3a204a9b517d6ba345来自CDS注释监督checkpoint；拟南芥GCF_000001735.4和水稻GCF_034140825.1属于作者报告的训练快照参考assembly。实际优化器窗口、验证划分与同源隔离仍未核实，不称纯无监督GLM、clean zero-shot或已证实DEV标签泄漏。来源与限制见[预训练暴露](../supplement/pretrained_exposure.md)。完整Carbon-A调用资产尚未具备固定重跑资格，见[接入核对](../../benchmarks/carbon_a/reports/carbon_a_pipeline_readiness.md)；保留其方法先例地位，不冒充已执行对标，也不改变冻结基线家族。

当前待检验的是具体改动能否在完整输出中恢复强基线可命中的原缺口，同时保留旧命中并控制背景负担，而不是凭组合名称宣布机制创新。原目标对照、native CDS prior、新抽样和24576共同编码仍按各固定TRAIN方案推进；组合分支另作独立四格比较。编码器上下文与whole-contig语法范围不是同一因素，新抽样也不只是“覆盖率”单因素。完整[八权重恢复读出](../../experiments/M28/full_grammar/reports/full_grammar_train_recovery_plan.md)尚未完成，LOCAL局部诊断不代表全家族结果，不预选赢家。

### 接续与发表判断

先完成既定TRAIN终态核验和主处理，再按[冻结目标对照DEV方案](../../experiments/M28/full_grammar/reports/full_grammar_objective_contrast_dev_plan.md)、[先验DEV方案](../../experiments/M28/full_grammar/reports/full_grammar_native_cds_prior_dev_plan.md)及[抽样／上下文DEV方案](../../experiments/M28/full_grammar/reports/full_grammar_sampling_context_dev_plan.md)推进对应固定家族；不按TRAIN排名增删候选、调阈值或换checkpoint。组合与新增恢复CPU不成为额外DEV闸门。原重复使用DEV只能支持开发比较，不能升级为独立泛化；新独立面板、Setaria和封存数据仍未放行。

正文主张仍需更好的完整链性能，或固定模型在具备独立资格的更广范围上的实证收益，并有与收益相连的机制结果。当前未满足：raw-FASTA工程接入和数字可信已经成立，但性能优越、机制因果、独立泛化与投稿准备度均不能据此宣布成立。负结果、背景权衡、unsupported／长链分层及资产暴露保留为补充证据，不改变论文核心目标或授权范围。

## 2026年10月7日 M28 阶段记录（已由上节更新）

当前主终点是同 strand、相同有序 CDS 坐标链的 primary-reference exact-chain P/R/F1，phase 和结构合法性另行检查。gene-body／CDS-base 覆盖不能代替完整链恢复，primary CDS 链也不等于完整 UTR 或所有 isoform。原两条 DEV 已反复用于开发，不能升级为独立测试；共同 7728 条 primary 分母及全部原背景必须保留，不按候选池、窗口宽度或预测可达性删难例。

已闭合的 A31_A23 在这两条 DEV contig 上合并完整链 F1 为 67.8133%，高于固定 Helixer 的 62.8512%，仍低于 Tiberius 的 71.8047% 和 ANNEVO 的 75.7025%。这一进展不能写成对三者普遍优越、两个完整基因组测试或独立泛化；逐物种及背景差异见[A31 DEV 结果](../../archive/m28_pre_full_grammar/reports/a31_local_path_caller_dev_result.md)和[原生基线结果](../../benchmarks/native_baselines/reports/native_baseline_result.md)。

A43 已在原六条 TRAIN 染色体、1491 个唯一窗口和 1536 draw×3 的范围完成一次拟合及完整自然读出。实际七项 CPU 数值证明成立，但原五条件全部未通过，拟南芥／水稻完整正确链均为0，模型按原规则关闭，详见[A43 完整结果](../../archive/m28_pre_full_grammar/reports/a43_result.md)。19189 个有限监督配对、原 gold、unknown 保护、末 004608 和零操作点未改变；自然预测不读取监督池或标签，搜索闭合后才按原 R4 读出。该实验不是全语法似然训练、在线难例挖掘或匹配所有监督曝光与梯度尺度的纯单因素消融。

| 发表主张 | 当前证据边界 | 下一步需要回答的问题 |
| --- | --- | --- |
| 完整基因结构性能更好 | 已有两条重复使用 DEV 的同尺原生比较，但 A31_A23 尚低于两个较强对标；A43 原 TRAIN 完整自然读出 FAIL0/5，不是 DEV 或 raw 结果。 | 新模型收益能否到达 reference-free 的 raw-FASTA 完整输出，并在逐物种、worst-species、完整链及背景负担上成立？ |
| 结构机制有实质贡献 | [A36](../../archive/m28_pre_full_grammar/reports/a36_structured_emission_result.md)原五条件通过 2 项并关闭，[A42](../../archive/m28_pre_full_grammar/reports/a42_emission_lattice_train_protocol.md)及 A43 均通过 0 项并关闭；较高解码目标不等于更准确。 | 控制表示、容量、监督支持和曝光后，收益来自哪里，是否保留到最终注释？GLM 加语法层本身不是已证实的新贡献，见[先例核对](../manuscript/M28_complete_chain_competition_prior_work_20261003.md)。 |
| 更广适用范围 | [迁移候选](../../archive/m28_pre_full_grammar/reports/independent_transfer_candidates.md)尚非冻结且获准执行的面板；GEN 上游暴露仍 unknown，各对照资产暴露不同。 | 固定模型及规则后，能否在有独立资格的物种和组装上保持效果，并分别披露[资产暴露](../supplement/pretrained_exposure.md)与同源边界？ |
| 完整部署可用 | 现有原生结果仅覆盖两条 DEV contig；人工跨窗／短输入检查和运行收据不等于全核基因组性能或公平速度比较。 | 在事前确定的完整核 contig 范围、碎片化输入和长链上，是否保留全部输入与参考，并量化完整流程的时间、GPU／CPU 和峰值内存？ |

正文应以完整 raw-FASTA 比较、独立适用能力和与收益相连的机制证据为主；补充材料保留负结果、unsupported／长链／短 contig 分层、暴露差异及计算条件。未知链和 nonexact-primary 都不能直接称作生物学假阳性。原窗口—链配对不是独立基因样本，重叠窗口也不是独立重复；不据这些配对数制造显著性。重复训练和统计确认留到固定候选模型及评价协议之后。

[A44 匹配拟合](../../archive/m28_pre_full_grammar/reports/a44_paired_fit_protocol.md)已固定原局部位置均值＋有限结构NLL／valid_bases的新损失单位，local-only与local＋structure使用同366191容量／seed0／原1536draw×3，系数分别1／0与1／1、不扫参；新两项单位检查13425929通过后，两臂首次提交为13425946／13425945。原[A44 零更新准备](../../archive/m28_pre_full_grammar/reports/a44_local_structured_protocol.md)的候选均值probe和A43损失不改变。新尺度不保证梯度平衡或性能收益；两臂完整自然TRAIN读出尚未闭合，JOINT须同时相对较强A31_A3和匹配LOCAL通过原五条件才继续本结构机制。不能断言缺少局部监督导致A43失败，不挑checkpoint、调阈值或加seed／epoch救援；真实raw／DEV、新域、Setaria及公开发布仍未放行。已授权的常规实施继续自主推进，不新增审批步骤。下方2026年9月26日及更早状态、贡献表和readiness勾选均为历史，不作为本轮科学目标、发表证据或授权依据。

## 2026-09-26 superseding status

The September 1 direct-annotator route below is also historical. The current authority for scientific direction is [research_roadmap_20260926.md](../../docs/research_roadmap_20260926.md). Same-scope M26 comparison places B below all three cached baseline callers in exact-chain F1; semantic source pairing closes the tested DEV source-dispute line; fixed-budget global boundary selection loses 225 supported chains net and fails its predeclared follow-up gate.

Do not use the old C1–C4 contribution menu, backbone-size comparison, or manuscript-readiness ticks below as current evidence of novelty, causal pretraining effects, clean generalization, or submission readiness. The user still seeks a publishable research contribution and ultimately Nature Communications; neither the fixed-route NO-GO nor a negative-result report completes that objective. A new paper claim remains unestablished. No manuscript-level superiority claim is supported.

## 2026-09-01 superseding status

The M12-era strategy below is retained as history, but its original broad fixed-caller claim is no longer current. M12–M17 showed that the NT-v2 route does not transfer as a universal fixed caller. M18/M19 established useful low-FPR GENERanno 1.2B adaptation, M24 showed that coarse gene-body quality does not translate into exact gene structures, and M25R stopped on development because complete-gene count is only `0.095–0.331x` reference across the frozen grid.

The current paper route is a direct raw-FASTA structural annotator, with exact CDS intervals/chains/genes, strand, phase, validity, FPR and count evaluated together. The live state and staged Nature Communications route are in `docs/28_current_research_state_2026-09-01.md` and `docs/29_end_to_end_technical_roadmap_2026-09-01.md`. No manuscript-level performance claim is currently supported.

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
