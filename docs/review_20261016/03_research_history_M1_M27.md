# 研究历程与方案谱系（M1 – M27，衔接 M28）

编制日期：2026-10-16。性质：只读梳理，未修改任何既有文件，未读取 Setaria 注释，未运行实验。
所有数值均摘自仓库内文件，并在括号或表格中标注来源路径（见下方"引用缩写"）。凡属本文的归纳、推断或判断，都明确写作"归纳""推断"，不与文件原文混淆。

## 0. 阅读范围、引用缩写与证据缺口

### 0.1 实际阅读的材料

`docs/06_results_log.md`（全文 1719 行）、`docs/10_findings.md`（全文）、`docs/12_publication_strategy.md`（含未提交修改，即 M28 相关新增段落）、`docs/14_validation_matrix.md`、`docs/15_evidence_register.md`（E001–E078 的摘要列）、`docs/25`–`docs/30`、`docs/experiments/`（6 个文件）、`reports/PROJECT-STAGE-REPORT-2026-08-16/`（report.html 提取文本）、`reports/` 下 M12B/M13/M18/M19/M20/M21/M22/M24/M25R/M26 各子目录的 report/summary/result 文件、`configs/` 中 M10–M25R 的 yaml、`src/screen_anchor/data.py::assign_splits`。M28 仅读取 `docs/M28_*_working.md`、`docs/12` 中 M28 段落及 `reports/M28-METHOD-RESTART/` 的少量结果标题，用于衔接。

### 0.2 引用缩写

| 缩写 | 路径 |
|---|---|
| LOG | `docs/06_results_log.md`（后附行号 Lxxx） |
| FIND | `docs/10_findings.md` |
| PUB | `docs/12_publication_strategy.md` |
| VM | `docs/14_validation_matrix.md` |
| EV | `docs/15_evidence_register.md`（后附 Exxx 编号） |
| D25–D30 | `docs/25_…` 至 `docs/30_…` |
| STAGE | `reports/PROJECT-STAGE-REPORT-2026-08-16/report.html` |
| R-M19 / R-M20 / R-M24 / R-M25R / R-M26 / R-M28 | 分别指 `reports/M19-COMPARABILITY-EVIDENCE/`、`reports/M20-SOTA-ERROR-ANALYSIS/`、`reports/M24-DIRECT-STRUCTURE-DIAGNOSTIC/`、`reports/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/` 与 `reports/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/`、`reports/M26-SAME-SCOPE-MECHANISM/`、`reports/M28-METHOD-RESTART/` |
| CFG | `configs/` |

### 0.3 本地仓库缺失、因而无法核验的材料（读者须知）

- 多处文档引用但本地不存在：`docs/11_evaluator_contract.md`、`docs/19_evaluator_contract.md`（评价器合同原文）、`docs/20_baseline_reproduction.md`、`docs/07_tri_review.md`、`docs/08_pivot_decisions.md`、`docs/09_decisions_log.md`、`docs/21_code_review_log.md`、`ACTIVE_GOAL.json`、`refs/dossiers/*`、`research_roadmap_20260926.md`（PUB 称其为 2026-09-26 起科学方向的"当前权威"，全仓库检索不到该文件）。评价器口径因此只能从 LOG、FIND、STAGE 的转述中还原。
- 本地 `outputs/` 仅含 M28 的 67 个目录，M1–M26 的原始 metrics.json 不在本地；M1–M26 数值全部来自 LOG/FIND/EV/reports 的转录，未能回溯到原始 json。
- **未找到名为 M27 的任何实验、报告、sbatch 或配置**（对 docs/reports/sbatch/scripts/configs/src 全文检索，唯一命中是 STAGE html 内嵌 base64 数据中的偶然字符串，不是实验记录）。时间线上 M26 之后直接进入 M28。本文把"M27"视为空缺，不虚构内容。
- 2026-07-01（M23）与 2026-08-24（M24）之间，LOG/FIND 无研究记录；STAGE 自述"最新研究记录停在 2026-07-01"（STAGE，"当前不确定性"一节）。

### 0.4 编号与口径陷阱（先读）

1. **编号不唯一**：M20 同时指 "M20-STRUCTURED-DECODER-IMPL-SMOKE"（LOG L1516）、"M20-SOTA-ERROR-ANALYSIS"、"M20-CLAIM-CLEAN-PANEL-FREEZE"；M21 同时指 GENERanno-CRF 筛选与"clean backbone provenance scout"；M22 同时指"非 CRF 路线设计"与"gb_tversky 筛选"；M25R 之后另有 Stage 1 R1–R4、case correction、phase-only、block assembly、run-graph R1 等一长串子分析，均挂在 M25R 名下。M1 内部的"M3–M8"是 Track A 的编号，与后来的 M3x 无关。
2. **指标名称相同、分母不同**。同一个 "gbF1/FPR/F1" 在不同阶段对应不同口径，见 §6.5 的"分母阶梯"。尤其：
   - M1–M23：以 CDS-span 的碱基级 gene-body F1 + 完整转录本 span 补集的 intergenic FPR 为主；
   - M24 起：同 strand 的 exact CDS interval / exact CDS-chain F1；
   - M25R/M26：6,450 条 complete-primary 参考链（拟南芥 4,151 + 水稻 2,299）；
   - M28：7,728 条 primary 参考链（拟南芥 5,437 + 水稻 2,291）与 TRAIN 22,688 条。
   这些数值不能跨阶段直接比较。
3. 对 M1–M23 的"gbF1"，R-M20 报告自己指出：旧表中标为 gbF1 的列其实是 `constrained_gene_body_F1`（FPR 超阈值即置 0 的硬门值），不是普通准确率；Helixer 的 0 不能解读为"零准确率"（R-M20 report.md，"2026-09-26 reporting correction"）。本文表中凡写"gbF1"，均尽量标明是 unconstrained（普通）还是 constrained。

### 0.5 证据级别定义（本文自定，用于表格"证据级别"列）

| 级别 | 含义 |
|---|---|
| L0 | 工程/冒烟/复现性检验，不含科学结论 |
| L1 | 探索性 screen：单种子或少量种子，开发面板上 val/test 染色体被反复使用，`ACTIVE_GOAL.status=draft`，文件自身标注 non-claim |
| L2 | 开发期冻结门槛：事前写定成功/停止门槛与产物，一次性执行，仍在被反复使用的开发染色体上；可做 go/no-go，不可作确认性主张 |
| L3 | 确认性：封存盲测（Setaria）或经资格审查的独立面板、冻结模型 | 

**迄今没有任何里程碑达到 L3。** Setaria 注释一直未被放行（R-M25R `terminal_summary.md`；PUB 2026-10-09 段"新独立面板、Setaria 和封存数据仍未放行"）。

---

## 1. 总览：一页纸版研究弧线

| 阶段 | 时间 | 里程碑 | 核心问题与走向 |
|---|---|---|---|
| A. 基线与评价器 | 2026-06-10 | M1 | 复现 Tiberius/Helixer/ANNEVO；发现跨工具 UTR 口径不可比，统一 CDS-span；两层 anchor（预训练上限 vs 同预算随机初始化）|
| B. 同预算架构筛选 | 06-11 – 06-12 | M3–M8 | 随机初始化结构化解码（CONSTR/CRF）→ 换尺子（intergenic 为主轴）→ 冻结 SegmentNT 特征 + FP-aware 损失 + 确定性后处理；发现 SegmentNT 在鸡/果蝇上标签泄漏；多类 CRF 无收益 |
| C. NT-v2 解冻与校准 | 06-14 – 06-16 | M9–M11 | 解冻 NT-v2-500m 顶部 L 层 → gbF1 与特异性同升；验证集阈值校准把 FPR 压到 ≤0.01（同物种池）|
| D. 泛化失败 | 06-17 – 06-19 | M12–M18 | 固定调用器跨物种全面失败（拟南芥→水稻、近缘、动物）；释放权重基线远强于己方；GENERanno 1.2B CDS 版显著优于 0.5B base |
| E. GENERanno 路线与封路 | 06-19 – 07-01 | M18–M23 | GENERanno 1.2B 低 FPR 稳定；CRF（M21）、Tversky（M22）无增益，NT-v2 直接迁移（M23）只是复现 M10 |
| F. 精确结构转向 | 08-24 – 09-26 | M24–M26 | 发现 gbF1 高 ≠ 结构精确；结构头（M25/M25R）no-go；Stage 1 诊断（含勘误）；无参考图解码 B 仍未达门槛；同范围重评分显示三基线完整链 F1 更高 |
| G. 全语法模型 | 10-03 起 | M28（另有 agent 负责）| 完整基因语法 likelihood + raw-FASTA 全读出；至 10-13 仍低于三基线 |

---

## 2. 时间线表

说明：每行"关键指标"尽量标注口径；"结论"用 成功/失败/无效实现/部分 标注；"证据级别"见 §0.5。

### 2.1 阶段 A：基线复现与评价器（M1，2026-06-10）

| 里程碑 | 研究问题 | 方法/架构 | 数据与物种 | 关键指标（口径）| 结论 | 证据级别 |
|---|---|---|---|---|---|---|
| M1-a：Tiberius mini-smoke / PILOT（LOG L1–L117）| 评价器与基线 runner 是否可用 | Tiberius 官方 SIF；评价器 `eval_gene_body_mask.py`（初版 transcript-span）| Tiberius 自带 Panthera pardus mini；S. cerevisiae + D. melanogaster（RefSeq）| mini-smoke：unconstrained gene-body F1 0.9196，FPR 0.0187，exact CDS F1 0.8594，transcript-chain F1 0.3124。PILOT 聚合：unconstrained 0.7087，FPR 0.0287，D.mel recall 0.5172 | runner 成功；聚合因 FPR 守卫被硬置零 → 暴露"guardrail 硬置零被误判为 failed_run"（已修 `completed_poor`，L119–L158）| L0 |
| M1-b：Helixer / ANNEVO smoke + 跨工具 span 统一（L161–L284）| 三工具是否可比 | `--span-mode cds`：只用 CDS(+start/stop) 作 gene-body span | 同上 | 统一前 S.cer FPR：Helixer 0.6544 vs Tiberius 0.0164（UTR 口径差）；统一后 Helixer 0.0333/F1 0.9869。三工具 CDS 口径 base-weighted F1：Tiberius 0.8608 / Helixer 0.9213 / ANNEVO 0.9197（macro 0.9150/0.9494/0.9429）| 成功：确立"CDS-only 为跨工具公共层"（FIND 2026-06-10 [M1-SPAN-HARMONIZE-CDS]）| L0–L1 |
| M1-c：same-budget screen_anchor（L287–L332）| 公平的"同预算"基线应是多少 | 随机初始化 tiberius_like (CNN+biLSTM) / helixer_like（空洞卷积），窗口 2048、8 epoch、3 seed、class-weighted CE | yeast + fly，染色体级切分（`assign_splits`，≈60/20/20）| FLOOR(ORF) 0.3735 < screen_anchor 0.5576 < 预训练上限 0.9213（CDS-span F1）。helixer_like 虽 0.5579 但 gene_count_ratio 51/94/153（严重碎片化）→ 被排除 | 成功：纠正"拿预训练 SOTA 当同预算 bar"的错误（FIND 06-10 "screen_anchor SEMANTICS"）| L1 |

### 2.2 阶段 B：同预算架构筛选（M3–M8，2026-06-11 – 06-12）

| 里程碑 | 研究问题 | 方法/架构 | 数据/物种 | 关键指标 | 结论 | 证据级别 |
|---|---|---|---|---|---|---|
| M3：TA-DECODER（L334–L389）| 结构化解码能否修复碎片化 | 固定 tiberius_like 骨干；(a) CONSTR 约束 Viterbi 后处理；(b) 线性链 CRF（python）；(c) semi-CRF；(d) 向量化 CRF-vec（log 空间结合扫描）| yeast+fly，3 seed | 旧尺子（gene-body F1）：anchor 0.5576（gene_count_ratio 2.74）→ CONSTR 0.5791（1.12）→ CRF-vec 0.6186（0.88，seed 方差 0.58–0.66）。python CRF/semi-CRF 因 W=2048 不可行被放弃 | 旧尺子下"成功"，随即被 M3b 推翻 | L1 |
| M3b：REVISE-INTERGENIC-PRIMARY（L391–L420）| 尺子对不对 | intergenic 改为 genome−**完整转录本 span（含 UTR）**；intergenic_specificity 成为主轴，gbF1≥0.5276 为底线（双主轴 Pareto）| 同上，重评分（无重训）| 新尺子下 anchor spec 0.8710；CONSTR 0.8369；helixer_like 0.7954；**CRF-vec 0.7138（最差，FPR 0.2862）**；FLOOR(ORF) spec 0.8805（最高）但 gbF1 0.3735 | **排名反转**：结构化解码靠外溢进 intergenic 抬召回；纯 specificity 可被"少预测"刷分 → 需双主轴（FIND 2026-06-11 两条 Finding）| L1 |
| FP-SEGMENTNT-PROBE（L422–L454）| 冻结基模特征的价值 | 冻结 SegmentNT(multi_species) 14 元素概率 → anchor 匹配的 conv+biLSTM 头（JAX 抽取 + torch 训练）| yeast+fly，3 seed | gbF1 0.6888（anchor 0.5576）；spec 0.8416（anchor 0.8710），macro 0.7543 | 部分：召回显著升，特异性不升（酵母最差）| L1 |
| M4：TA-FOUNDATION-DECODER（L456–L483）| 用什么手段把召回转化为特异性 | 在冻结 SegmentNT 特征上比较：**FP-aware 损失** / raw-DNA 融合 / CRF 解码 | yeast+fly，5 seed×8 epoch | FPLOSS spec 0.9303±0.036、gbF1 0.6157、gene_count 2.25（Pareto 胜 anchor）；FUSION 0.8615/0.6850/3.40（负）；CRF 0.8298±0.119/0.6840/0.90（方差大）| 成功（FPLOSS）。2-epoch 单 seed 冒烟曾得相反结论 → "AXIS-1 判定须 ≥5 seed" | L1 |
| M5：COHERENCE-FIX（L485–L507）| 去碎片化 | FPLOSS + **确定性约束后处理**（合并小间隙/去小 CDS）| 同上；anchor 扩到 5 seed | CONSTR：spec 0.9272±0.036、gbF1 0.6581、gene_count 1.28；5-seed anchor 0.8436±0.066（3-seed 为 0.8710）；配对差 +0.0836±0.037（5/5 为正）| 成功；同时暴露"3-seed anchor 不可靠" | L1 |
| M6：FRAGFIX-SWEEP（L509–L536）| 后处理参数能否在**无测试泄漏**下达标 | VAL 上扫 max_fill_gap×min_cds_len，一次性应用于 TEST | 同上 | 采用 mfg20/mcl60：spec 0.9218±0.021、gbF1 0.6439、macro 0.8331、gene_count 1.037±0.312（seed 间 0.55–1.35）| 成功（门槛通过）；但文档对"mcl=90 还是 mcl=60"有前后修订（L516 vs L535）| L1 |
| M7：REANCHOR-HELDOUT（L538–L570）| 优势是否只是 yeast+fly 的低 UTR 假象 | 同候选在 Arabidopsis + Gallus（UTR 42%/62%）重测；ANNEVO 作上限 | A. thaliana 全；G. gallus 子集（NC_≤20Mb 共 30 seqid） | 候选 spec 0.9604（anchor 0.8054，ANNEVO 0.9824）；gbF1 0.6664（anchor 0.7099）| 原结论"跨 clade 成立"；**随后被 M8-CK4 审计推翻**：SegmentNT 分割头在鸡/果蝇等上被微调，鸡的 +0.30 是标签泄漏，仅拟南芥（+0.06）干净（L568；FIND 06-12 "CRITICAL"）| L1，后部分作废 |
| M8：TB-GBF1-MULTICLASS（L572–L595）| 8 类（CDS 相位/内含子/UTR/剪接位点）+ CRF(8×8) 能否补 gbF1 缺口 | 同冻结 SegmentNT 特征 | 干净植物 {拟南芥, 水稻} | 3 类候选 spec 0.9663 / gbF1 0.7392 / gcount 0.936；多类 0.9683 / 0.7189 / 0.6625；raw-DNA anchor 0.9045 / 0.6960 / 3.46 | 多类 **失败**（gbF1 更低、严重欠预测）；3 类对 raw-DNA 的 Pareto 增益是"干净"的（+0.062 spec，+0.043 gbF1）。结论：冻结特征封顶 gbF1 | L1（3c n=3，2/5 seed 因 beegfs 读错失败）|

### 2.3 阶段 C：NT-v2 解冻与 GENERanno 初探（M9–M11，2026-06-14 – 06-16）

| 里程碑 | 研究问题 | 方法/架构 | 数据/物种 | 关键指标（拟南芥/水稻，unconstrained）| 结论 | 证据级别 |
|---|---|---|---|---|---|---|
| M9 / M9-DEEP：TB-UNFREEZE（L597–L616，L655–L721）| 解冻 NT-v2-500m 顶部 N 层是否打破冻结特征的 gbF1 上限 | NT-v2-500m（ESM），顶部 L0/L2/L4/L6/L8/L12 层可训，3 类 convLSTM 头 + FP-aware 损失 + 约束后处理；`generanno` 环境（transformers 4.49）| 单物种拟南芥，单 seed | gbF1：L0 0.8284 → L2 0.8544 → L4 0.8759（FPR 0.0246，>0.02 被硬置零）→ L6 0.8843 → L8 0.8840 → **L12 0.9035，spec 0.9810，FPR 0.0190，gcount 0.792**；ANNEVO 参考 ≈0.898 | 成功（路线内）：解冻单调同升 gbF1 与特异性；FPR 悬崖来自硬门而非"bug"（FIND 06-14 指标审计）| L1（单 seed 单物种）|
| GENERANNO-NATIVE-PROBE（L618–L653）| 官方 GENERanno-1.2B-CDS 开箱行为 | 官方 binary CDS mask + 官方 stair-refine，或我方 constrained 后处理（无微调）| 拟南芥测试染色体 | spec 0.9947 / FPR 0.0053 / gbF1 0.7709 / gcount 4.219；换后处理 0.9960 / 0.7607 / 3.844 | 特异性极强但基因被碎片化：根因是 binary CDS mask 无内含子概念 → 后处理无法归并外显子为基因 | L1 |
| M10：M9L12-CLEANPLANTS + GENERanno LoRA smoke（L723–L884）| L12 能否跨种子、跨干净植物；GENERanno+LoRA+3 类头能否跑通 | NT-v2 L12（window 2046，lr 1e-3，backbone_lr 1e-5，fp_lambda 1.0）；GENERanno 1.2B + LoRA(r=8) 8 窗 smoke | {拟南芥, 水稻}，3 seed（同物种池内染色体切分）| L12 均值：spec 0.9826、FPR 0.0174、macro spec 0.9801、gbF1 0.8398、gcount 0.897（拟南芥 gbF1≈0.90，水稻≈0.73–0.77）。GENERanno smoke：spec 0.9491、gbF1 0.7525、gcount 4.43 | L12 路线验证；剩余阻塞为 FPR 尾部（0.0142–0.0212>0.01）。GENERanno smoke 仅证明工程路径通 | L1 |
| M11：L12-SPEC-CALIBRATION（L886–L970）| FPR≤0.01 是否只是工作点问题 | 保存 raw logits；VAL 上扫 intergenic bias {0…4.0}、min_cds_len {60,90,120}、max_fill_gap {0,20}，一次性施于 TEST | 同 M10，3 seed | 均值 spec 0.9913、**FPR 0.0087**、macro 0.9909、gbF1 0.8178、gcount 1.003；per-species 拟南芥 seed2 FPR 0.0111、水稻 gbF1 0.683–0.711 | 成功（同物种池内）。后被 M12A 证明"这不是固定跨物种模型"| L1 |

### 2.4 阶段 D：泛化与外部基线（M12–M18，2026-06-17 – 06-19）

| 里程碑 | 研究问题 | 方法 | 数据/物种 | 关键指标 | 结论 | 证据级别 |
|---|---|---|---|---|---|---|
| M12A：固定模型跨物种（L972–L1027）| M9-L12 + 校准能否作为固定调用器 | 只在拟南芥训练与校准，水稻不参与任何选择 | 拟南芥→水稻，3 seed | 水稻：gbF1 0.6556、spec 0.9689、**FPR 0.0311**、gcount 1.755、constrained@0.01=0 | **失败**；解释为低数据/同域适配模型，非通用固定调用器 | L1 |
| M12B：同面板释放权重基线 | 三基线在干净植物面板的位置 | ANNEVO Magnoliopsida、Tiberius angiosperms、Helixer land_plant_v0.3_m_0100 | 拟南芥+水稻（**整基因组**参考，见 §6.5）| ANNEVO 0.9269 / FPR 0.0117 / gcount 0.726；Tiberius 0.9252 / 0.0073 / 0.628；Helixer 0.9220 / 0.0216 / 0.820（R-M19 `comparison_tables.md`）| 释放权重显著强于己方固定模型，但各有 under-call 或 FPR 取舍 | L1（**范围不匹配**，见 §6.5）|
| M12C：GENERanno 公平挑战 smoke（L972 内）| CDS 特化基模 vs 通用 base | 1.2B-cds-annotator-preview vs 0.5B-base，同 LoRA+3 类头 | 干净植物 smoke | 1.2B：gbF1 0.7527、FPR 0.0432、gcount 4.405；0.5B：gbF1 0.5982、spec≈0.0001（塌缩）| 1.2B 有信号；0.5B 塌缩（后被 M15 修正为"过严"）| L1 |
| M13：距离扫描 + 失败机理（L1030–L1092；R-M13_FAILURE_SANITY）| 失败是否只因 rice 太远 | 拟南芥训练/校准，测近缘 A. lyrata 与水稻；先做失败机理（M12A vs M11）| A. lyrata（GCF_000004255.2，top-8 scaffold，仅诊断用）+ 水稻 | VAL 拟南芥 FPR 0.0085/gbF1 0.8965 完美；TEST 聚合 FPR 0.0340、gcount 1.616；A. lyrata FPR 0.0355/gcount 1.358；rice 0.0328/1.990。机理：水稻 pre-decode 把真 intergenic 判为 genic 的比例 0.2320（M11 为 0.0161）；test-oracle 网格 0/3 seed 可同时 FPR≤0.01 与 gcount≤1.25 | **失败**：近缘也失败，不是单纯距离问题，也不是可转移的校准 bug | L1（A. lyrata 为 scaffold 级）|
| M14：动物负控（L1095–L1143）| 固定拟南芥模型对动物 | 同 M13 | Gallus、Drosophila（诊断用）| 聚合 FPR 0.0406、gbF1 0.5448、gcount 2.254；gallus FPR 0.1156/gcount 4.987；fly FPR 0.0085 但 gbF1 0.5274 | 失败（负控）| L1 |
| M15：GENERanno LoRA 面板（L1147–L1195）| 在更大预算下 CDS 版 vs base | LoRA(q/k/v/o,r=8)+3 类 FP-aware 头 | 干净植物 | 1.2B：gbF1 0.8510、FPR 0.0258、gcount 1.177；0.5B：0.7623、0.0562、1.326 | CDS 特化显著更好；0.5B 不再塌缩（"base 完全不可用"表述过强）| L1 |
| M16：多物种训练（L1198–L1258）| "训练物种太少"假说 | NT-v2 L12，拟南芥+水稻训练/校准 | 测 A. lyrata、gallus、drosophila | VAL 点 FPR 0.0087；TEST 聚合 FPR 0.0197、gbF1 0.5615、gcount 1.326；gallus gcount 4.987→2.368 | 部分：多物种帮助一致性，不解决固定模型泛化 | L1 |
| M17：同面板释放基线（L1261–L1316）| 释放权重在宽面板表现 | 同 M12B 的三基线，加 Aves/Insecta/vertebrate 等权重 | A. lyrata、水稻、gallus、drosophila（整基因组）| ANNEVO 0.9115/FPR 0.0240/0.840；Tiberius 0.8791/0.0173/0.556；Helixer 0.8797/0.0526/0.931 | 己方 M16（0.5615）远不具竞争力；动物行有已知预训练重叠 | L1（与 M16 范围不一致，见 §6.5）|
| M18-a：多 clade 训练诊断（L1318–L1372）| 加 drosophila 训练是否足够 | NT-v2 L12，训练/校准拟南芥+水稻+果蝇 | 测 A. lyrata、gallus、yeast | 聚合 FPR 0.0444、gbF1 0.6170、gcount 1.572；gallus FPR 0.1359/gcount 3.994；yeast gcount 0.591 | 部分：近缘植物改善（A. lyrata FPR 0.0355→0.0184），广谱仍失败 | L1 |
| M18-b：gallus test-oracle 校准（R-M18 `report.md`；LOG L1357）| gallus 失败是否仅因工作点 | 用 test 标签扫 bias/mcl/mfg（诊断）| gallus | 满足 FPR≤0.01 & gcount≤1.25 的最优点退化：FPR 0.0、gbF1 0.0060、gcount 0.183；最高 gbF1 点 FPR 0.7347 | 失败不是阈值问题 | L1（诊断，含 test 标签）|

### 2.5 阶段 E：GENERanno 路线与封路（M18–M23，2026-06-19 – 07-01）

| 里程碑 | 研究问题 | 方法 | 数据/物种 | 关键指标 | 结论 | 证据级别 |
|---|---|---|---|---|---|---|
| M18-c：GENERanno 0.5B 加强 FP 目标（L1375–L1411）| base 能否被目标救活 | 0.5B-base + LoRA，`fp_lambda=2.5`，`min_cds_len=90` | 干净植物 | FPR 0.0967、gbF1 0.6561、gcount 1.617；rice gbF1 0.4981/FPR 0.1239 | **失败**；仅作消融 | L1（单 seed）|
| M18-d：GENERanno 1.2B 加强 FP 目标（L1413–L1457）| CDS 版在同目标下 | 1.2B-CDS + LoRA + 3 类头，`fp_lambda=2.5`，`min_cds_len=90` | 干净植物 | FPR 0.0071、gbF1/constrained@0.01 0.8494、gcount 0.864、macro spec 0.9943 | **成功**（首个 guardrail-valid 的 GENERanno LoRA）| L1（单 seed）|
| M19：RAWCAL（L1460–L1514）| M18-d 是否偶然；raw-score VAL 校准 | 同上，两 seed，保存 raw logits，VAL 校准 | 干净植物 | s0：raw 0.8390/FPR 0.0088/0.967 → 校准 0.8421/0.0083/1.083；s1：raw 0.8593/0.0059/0.805 → 校准 **0.8815/0.0065/0.830** | 成功（稳定；内部最强的低 FPR 点）；落后释放基线 gbF1 ≈0.92–0.93（但范围不匹配）；GENERanno 来源 overlap_unknown 阻塞 claim | L1（2 seed）|
| M20：CLAIM-CLEAN-PANEL-FREEZE / SOTA-ERROR-ANALYSIS / CRF smoke（FIND 06-21；L1516–L1544；R-M20）| 能否声称 GENERanno 干净 held-out；错误集中在哪 | 来源核查；同面板误差分析；CRF 实现烟测（8 窗，gbF1=0 属设计内）| Arabidopsis/rice | 结论：GENERanno 的 Arabidopsis/rice 来源 overlap_unknown，claim 被阻塞；**R-M20 聚合排名后被标注"不可作同范围排名"** | 来源 gate 关闭；CRF 工程通 | L0–L1 |
| M21：GENERANNO-CRF-SCREEN（L1546–L1590）| 训练 CRF 解码器能否恢复召回而保持低 FPR | GENERanno 1.2B LoRA + CRF 头（`crf_aux_ce=1.0`）| 干净植物，s0；s1 因共享队列 TIMEOUT 后由 `s1-opt` 救援 | s0：gbF1 0.8544/FPR 0.0273/0.956；s1-opt：0.8744/0.0192/0.690（M19 s1 为 0.8815/0.0065）| **失败**：更差且丢失低 FPR；3/3 tri-review 放弃（FIND 06-22）。全窗 CRF 训练 ≈8–10 s/batch（A100）| L1（2 seed）|
| M22：非 CRF 路线设计 + gb_tversky（L1592–L1645）| 校准/后处理空间是否已耗尽；标量辅助损失 | (a) 在 M19 校准网格上做 test-oracle 重放；(b) GENERanno 1.2B + gene-body Tversky 辅助损失（α=0.8, β=0.2, λ=2.0），`decoder=none` | 干净植物 | (a) s1 校准点即 oracle 最优（0.8815/0.0065）；s0 仅有小幅隐藏收益 0.8421→0.8603；(b) gbF1 0.8635、FPR 0.02865、gcount 0.648，promotion gate=false | (a) 校准空间耗尽；(b) **失败**，不加 seed、不调参（seed1 未提交）| L1（(b) 单 seed）|
| M23：NTV2-CLEAN-TRANSFER-s0（L1647–L1701）| 干净来源 NT-v2 路线是否仍是 claim 方向 | 与 M10 完全相同配方、同 seed、无校准/CRF/Tversky | 干净植物 | gbF1 0.8427、FPR 0.01673、gcount 0.867——**与 M10 s0 数值逐位相同** | 工程成功，**无新信息**（确定性复现）；弱于 M19 | L0 |

### 2.6 阶段 F：精确结构转向（M24–M26，2026-08-24 – 09-26）

| 里程碑 | 研究问题 | 方法 | 数据/范围 | 关键指标 | 结论 | 证据级别 |
|---|---|---|---|---|---|---|
| M24：DIRECT-STRUCTURE-DIAGNOSTIC（LOG L1703；R-M24）| 已保存的 M19/M8 产物与释放基线在**同 seqid** 上的结构精确度；SegmentNT 缓存有无剪接信号 | 纯 CPU 评价；exact CDS-run interval、coordinate pseudo-chain；基线按 primary transcript（最长 CDS）| M19 held-out seqid：拟南芥 NC_003075.7（18,585,056 bp；3,161 条完整 primary）、水稻（28,605,474 bp；1,981 条）| M19：exact interval F1 0.0531–0.1498，pseudo-chain 0.0082–0.0123；M8 3 类 ≤0.0207；ANNEVO interval 0.8614/0.8882、chain 0.7479/0.7324；Tiberius 0.8547/0.8631、0.7254/0.6660；Helixer 0.8117/0.8121、0.6339/0.5850。SegmentNT 6-kb 缓存：exon AUCPR 0.5866–0.6569，donor/acceptor 仅 0.0314–0.0443 | **改变目标**：gbF1 高并不意味结构精确；M19 的 strand/phase 是占位符。SegmentNT 结论仅限该 6-kb 缓存 | L1（范围明确）|
| M25（首次）（EV E072）| 结构头能否输出完整模型 | GENERanno 1.2B + 区域/边界/相位多头 | 同上 | 4/1,536 训练窗结构 mask 全空 → 边界损失空均值非有限 | **无效实现**（checkpoint 不可用）| L0 |
| M25R（L1711–L1719；D26/D27；R-M25R）| 修复后同一冻结实验的开发期准入 | GENERanno 1.2B（官方 6-mer，6,144 bp 窗，LoRA r=8）+ 区域头(fp_aware, λ=2.5) + 4 边界头(focal, coef 2.0) + 4 相位头(CE, coef 0.5) + 受约束 HMM 式状态解码（I/C/G，GT-AG，±6 bp 吸附）；3 epoch；seed 0 单次拟合；5,625 个（epoch×3、region 阈值×3、边界阈值 5⁴）验证网格 | 训练：拟南芥+水稻 train seqid；验证：NC_003074.8 + NC_089041.1（共 6,450 条参考链）；Setaria 封存 | 训练 loss 0.8205→0.6669→0.6177；**0/5,625 可准入**：全部 FPR≤0.020，但 gene-count ratio 仅 0.0953–0.3305（要求 0.80–1.20）；最佳 tuple（epoch1）interval 0.1204、chain 0.3250、FPR 0.01247、count 0.3253；epoch 越晚越差（count 0.3253→0.1829→0.1034）| **有效 no-go**（针对此冻结组合）；不能据此判定 GENERanno 表示本身不足 | **L2**（冻结门槛）|
| M25R-Stage1：R1–R3（EV E076–E077；R-M25R-DEV `terminal_summary.md`、`r3_terminal_summary.md`）| 在不重训的前提下定位丢基因的阶段 | 重放 3 个 epoch 的各一个冻结 tuple；逐参考链归因；独立 GFF3 合法性审计 | 同上 | R1 因"非法即致命"的诊断控制流 bug、R3 因变量遮蔽 `KeyError: 'epoch'` 失败；R3 部分证据：earliest-stage 计数 region 105、无覆盖 block 1,922、CDS-run 数不符 586、motif 601、donor/acceptor 473、边界选择 229、phase 1,145 等 | **无效实现（诊断代码缺陷）**，部分账本有效 | L0 |
| M25R-Stage1：R4 + 大小写勘误（R-M25R-DEV `r4_terminal_review.md`、`r4_case_correction_report.md`）| 完成全部 3 epoch 诊断并勘误 | R4 全量重放；随后 CPU-only 勘误 | 同上 | 原结论"≈36% emitted 结构非法"（1338/2098 合法）**撤回**：`read_fasta` 不转大写，独立审计未大写化；大写后 2098/2098、1180/1180、667/667 全部合法。canonical-compatible 参考 2666→**6065/6450（94.03%）**；motif 可达 1693→3616。truth-assisted 897/577/320 → `NOT_RECOVERABLE_WITH_SAVED_ARTIFACTS`。主归因仍为 mixed attrition（no single stage ≥50%）| 诊断完成；撤回错误解释；NO-GO 不变 | L2 |
| M25R-E1 phase-only 上界 / block assembly（R-M25R-DEV `phase_only_feasibility_bound.md`、`block_assembly_result.md`）| 仅取消 phase 硬拒绝能否达门槛；允许跨 block 的图是否有空间 | 纯 CPU 上界推导与可达性 oracle | epoch-1 / row601 | phase-only 最优上界 F1=2×2538/(6450+3247)=**0.5235<0.55**（关闭）；跨 block 的 canonical 完整路径存在 3,585/6,450，oracle 上界 0.7145 | phase-only 关闭；组装机制未被排除（值得一次无参考实测）| L2 |
| M25R-E1 run-graph R1（B 臂）（R-M25R-DEV `run_graph_decoder_contract.md`、`run_graph_terminal_20260926.md`；EV E078）| 固定 epoch-1 权重，用无参考的候选 CDS 有向图 DP 替代"整 block + phase 硬拒绝" | 候选 CDS 图（run 支持、≥4 bp 内含子、完整 ORF）+ 固定无拟合评分 S(P)=Σ(zC−zI)+Σ(zG−zI)+Σ边界 logit 差；全局最大化跨度不重叠基因集 | 同上，6,450 参考 | A→B：exact interval 0.1204→**0.7345**；chain 0.3250→**0.5303**；count ratio 0.3253→0.9594；**FPR 0.01247→0.05606**；独立合法 6,188/6,188；strand 0.9996；phase 0.9962。门槛：interval≥0.80 ✗、chain≥0.55 ✗、FPR≤0.020 ✗（且触发 >0.030 立即停止线）| **NO-GO 结案**（固定 M25R 直接结构注释路线）。关闭范围限定为当前 checkpoint+heads+已测解码，非"GENERanno 架构普遍不可能"| **L2** |
| M26：SAME-SCOPE-MECHANISM R1–R3（R-M26）| 在同一坐标/参考/评价器下，内部改进是否仍落后于既有工具 | 对 5 套已保存预测重评分（无新推理）；R2 背景归属；R3 实际图支持分解（仅有协议，**本地未见结果文件**）| 同 6,450 参考（NC_003074.8 + NC_089041.1）| 完整链 TP/预测数/F1：A 1,389/2,098/0.3250；B 3,351/6,188/0.5303；ANNEVO 5,153/7,182/0.7560；Helixer 4,449/7,899/0.6201；Tiberius 4,956/7,427/0.7143。B 的 2,837 条 nonexact：精确匹配其他 isoform 158、intron 链相同端点不同 59、同链 span 重叠 1,537、仅反链 6、无重叠 1,077。**背景政策敏感性**：拟南芥 B FPR 12.7439%→1.3113%（改全注释 gene 背景）；Arabidopsis B 旧 FP 碱基 83.19% 落在被标 partial 的编码基因区，水稻仅 0.196% | 内部改进**不是**基线优势；特异性解释强烈依赖参考背景政策；"拟南芥高 FPR 源于重复/短 ORF"不成立 | L2（开发范围，探索性诊断）|

### 2.7 阶段 G：M28 衔接（简要；详见另一 agent 的报告）

| 项 | 要点 | 来源 |
|---|---|---|
| 范围转换 | 从"有限候选集训练"转为"全语法模型 + 严格 raw-FASTA 全染色体读出"；冻结 GENERanno-1.2B-CDS-preview（revision b0483c23…）+ 366,191 参数的 GRU 头；完整基因语法 likelihood（部分观测、coarse channel ε=0.1、未知位置边缘化）| `docs/M28_full_grammar_methods_working.md` |
| 分母 | DEV：拟南芥 NC_003074.8 + 水稻 NC_089041.1（即 M25R 的 val 染色体），primary 5,437+2,291=7,728；TRAIN：6 条染色体（拟南芥 NC_003070.9/003071.7/003076.8，水稻 NC_089043.1/089044.1/089046.1），22,688 条 primary | `R-M28/native_baseline_result.md`、`full_grammar_raw_train_result.md` |
| DEV 比较（完整链 F1，pooled）| A31_A23 67.8133%（Helixer 62.8512% < A31 < Tiberius 71.8047% < ANNEVO 75.7025%）；joint 66.623%；fresh-joint 70.107%、uniform-coarse 69.795%、codon-energy 64.579%；均高于 Helixer，低于 Tiberius/ANNEVO | `R-M28/a31_local_path_caller_dev_result.md`；`docs/M28_full_grammar_abstract_working.md` |
| TRAIN 同尺（pooled F1）| M28 原 4608 步 45.590；LOCAL 末 9216 步 41.729；ANNEVO 75.542；Tiberius 72.655；Helixer 63.592 | PUB（2026-10-09 段表格）|
| 负结果闭环（A 系列）| A6/A7、A9、A11、A13、A14、A15、A21、A36（2/5 通过）、A37（2/5）、A42/A43（0/5）等均在事前门槛下 FAIL/关闭；FULL Q vs C0 冻结 DEV 门槛 FAIL | `R-M28/*_result.md` |
| 暴露与独立性 | GENERanno 固定权重经 Carbon-A 来源映射到 RefSeq CDS 注释监督；两项目 assembly 属于其训练快照参考子集，实际优化窗口/验证/同源隔离未核实 → 不可称 clean zero-shot | `R-M28/pretrained_exposure.md`；`docs/M28_full_grammar_methods_working.md` |
| 当前判断（PUB）| "工程接入与数字可信成立，性能优越、机制因果、独立泛化与投稿准备度均不成立" | PUB（2026-10-09 段，"接续与发表判断"）|

---

## 3. 方案谱系（文字树）

图例：`[×]` 放弃/关闭（附原因）；`[→]` 被继承的发现或组件；`[◇]` 仍开放或有条件。

```
M1 基线与评价器
├─ [→] 评价器：CDS-span 为跨工具公共层（Helixer/Tiberius UTR 口径差）                   → 全程沿用
├─ [→] 两层 anchor 思想：预训练上限 ≠ 同预算 bar（0.9213 vs 0.5576）                    → 后续所有"同预算"比较
└─ [×] 用 yeast+fly 做 anchor：低 UTR、基因密、FPR 分母小，仅 provisional（M7 重推）

M3–M8 同预算架构筛选（随机初始化 / 冻结 SegmentNT）
├─ M3 结构化解码（旧尺子 gene-body F1）
│   ├─ [×] python 线性链 CRF / semi-CRF：W=2048 下不可行（工程）
│   ├─ [×] CRF-vec "旧尺子赢家"：换尺子后特异性最差（外溢 intergenic）→ Track-B 晋升在花 GPU 前被叫停
│   └─ [→] CONSTR 确定性约束后处理：以最小代价修复碎片化   → M5/M6 及此后所有 constrained_decode
├─ [→] 发现：尺子选择会翻转结论；需双主轴（specificity + gbF1 底线 + macro 稳定）
├─ FP-SEGMENTNT-PROBE → M4：冻结特征 + FP-aware 损失 Pareto 胜 anchor
│   ├─ [→] FP-aware 损失（对 true-intergenic 的 genic 概率质量惩罚）     → M9 起所有 3 类头（fp_lambda 1.0 → 2.5）
│   ├─ [×] raw-DNA 融合（低于 anchor）
│   └─ [×] 学习型 CRF（方差大，specificity 下降；但 gene_count 最佳）→ 后在 GENERanno 上又验证（M21）仍负
├─ M5/M6：[→] "学习型 CRF 不如确定性后处理 + VAL-only 选参"（无测试泄漏）
├─ M7：[×] "跨 clade 优势成立"的原结论 → 被 M8-CK4 审计推翻（鸡/果蝇 SegmentNT 分割头标签泄漏）
│   └─ [→] 保留的干净信号：拟南芥 +0.06；并产生"干净物种清单"原则（植物排除于 NT-v2 与 SegmentNT 监督）
└─ M8：[×] 8 类多类 + CRF(8×8)：gbF1 不升（0.7189<0.7392），严重欠预测
    └─ [→] 诊断"冻结特征封顶 gbF1" → 引出 M9 解冻

M9 NT-v2 顶部 N 层解冻
├─ [→] 解冻 L12：gbF1 0.9035、spec 0.9810（拟南芥单 seed）；M10 三 seed 0.8398
├─ M11 VAL-only 校准 [→] 沿用到 M19（raw-score 校准）；[×] 单靠校准无法跨物种（M13 oracle 0/3）
└─ M12A–M18-a 固定跨物种调用器
    ├─ [×] 拟南芥→水稻（M12A）、近缘 A. lyrata（M13）、动物（M14）：全部失败
    ├─ [×] 多物种训练（M16、M18-a）：仅部分改善；gallus 失败非阈值问题（M18-b oracle）
    └─ [→] 结论：该路线只能算"低数据/同域适配"，不是通用调用器 → 论文主张降格

GENERanno 路线（来源 overlap_unknown）
├─ 开箱（NATIVE-PROBE）：特异性极高，但 binary CDS mask 无内含子 → gcount 4.2×；后处理无解
├─ [→] CDS 特化 > 通用 base（M12C、M15、M18-c/d、M19）                          → 选择 1.2B-CDS 作为后续唯一 GENERanno 基模
├─ [×] 0.5B-base：加强 FP 目标也不行（M18-c）
├─ [→] 1.2B + LoRA + 3 类 FP-aware 头 + VAL 校准：低 FPR 稳定（M19 FPR 0.0065–0.0083）
├─ [×] 训练 CRF（M21）：丢失低 FPR；全窗 CRF 8–10 s/batch
├─ [×] gb_tversky 标量辅助损失（M22）：gbF1、FPR 双差
└─ [×] 后处理/校准继续挖（M22 oracle 重放：校准空间已耗尽）

M23 NT-v2 干净来源直接迁移：[×] 与 M10 s0 数值逐位相同，无新前沿（确定性复现）

M24 结构精确度（转折点）
├─ [→] 发现：gbF1/FPR 高 ≠ exact interval/chain（0.05–0.15 vs 0.81–0.89）
├─ [→] 评价主终点改为 exact CDS interval / chain / gene + strand + phase + 合法性
├─ [×] SegmentNT 6-kb 缓存作为首个 GPU 分支（donor/acceptor AUCPR 0.03–0.04；仅限该缓存）
└─ M25/M25R GENERanno 结构头 + 约束 HMM 式解码
    ├─ [×] M25：空 mask 非有限损失（无效实现）
    ├─ [×] M25R：0/5,625 可准入，count 0.095–0.331（有效 no-go）
    ├─ Stage1 诊断
    │   ├─ [×] R1/R3 诊断控制流 bug；R4 大小写 bug → "36% 非法"撤回
    │   ├─ [×] phase-only 修正：上界 0.5235<0.55（关闭）
    │   └─ [→] 保留：mixed attrition（无单一阶段 ≥50%）；motif/transition 可达率高、phase 检查与 block 组装损失大；跨 block 图有 oracle 空间
    └─ run-graph R1（B）：[→] 联合图解码释放大量结构（chain 0.3250→0.5303），但 [×] interval/chain/FPR 三项门槛失败 → 固定 M25R 路线 NO-GO 结案
        └─ [→] M26：同范围下三基线 chain F1（0.62–0.76）均高于 B（0.5303）；FPR 结论对背景政策敏感
            └─ [◇] 背景/可评价区域定义（partial 编码基因）成为独立研究问题

M28 全语法模型（衔接）
├─ [→] 继承：raw-FASTA 全染色体读出；exact-chain F1 主终点；primary 分母不删难例；"事前冻结门槛 + 闭合即停"
├─ [→] 继承：GENERanno 1.2B-CDS 表示（暴露未核实）、同尺原生基线对标
└─ [◇] 开放：M28 各固定权重在 TRAIN/DEV 上仍低于 Tiberius/ANNEVO；独立面板与 Setaria 未放行
```

放弃原因归纳（本文归纳）：
1. **口径/尺子问题**（M3 旧尺子、M20 聚合、R4 大小写）——先修尺子再谈模型。
2. **泄漏/来源问题**（M7 鸡/果蝇、GENERanno overlap_unknown）——使"干净 claim"路线受阻。
3. **科学无增益**（M8 多类、M21 CRF、M22 Tversky、M25R、phase-only、R1 B 门槛）——均有事前门槛，失败后不加 seed/不调阈值。
4. **工程可行性**（python CRF、全窗 CRF 吞吐）。
5. **无新信息**（M23）。

---

## 4. 已确立的"可靠正面结论"与"可靠负面结论"

分级说明：A = 多 seed/多物种/多口径一致；B = 单口径但冻结门槛或同范围测量；C = 方向可信但证据薄或有明确限定。

### 4.1 可靠正面结论

| # | 结论 | 强度 | 关键数值 | 证据位置 | 限定 |
|---|---|---|---|---|---|
| P1 | **CDS 特化的 GENERanno-1.2B 明显优于通用 0.5B base**（同 LoRA + 3 类 FP-aware 头 + 同目标）| B（限定后）| 1.2B：FPR 0.0071 / gbF1 0.8494 / gcount 0.864；0.5B：0.0967 / 0.6561 / 1.617（M18-c/d）；M15：1.2B 0.8510/0.0258 vs 0.5B 0.7623/0.0562；M19 两 seed 稳定 | LOG L1375–L1514；FIND 2026-06-17/18/19；R-M19 `comparison_tables.md` | **与模型规模混淆**（1.2B vs 0.5B），0.5B 仅单 seed；D30 明确"不能把 1.2B vs 0.5B 全部差异归因于 CDS 后训练"；基模来源 overlap_unknown（M28 进一步得知两 assembly 在 Carbon-A 训练快照参考子集）|
| P2 | **冻结基模特征 + FP-aware 损失 + 确定性约束后处理**在同预算筛选中 Pareto 胜 raw-DNA anchor；在**干净植物**上亦成立 | B | M5 配对差 +0.0836±0.037（5/5 seed）；M8 干净植物 spec 0.9663 vs 0.9045、gbF1 0.7392 vs 0.6960 | LOG L485–L507、L572–L595；FIND 2026-06-11/12 | 仅 screen；anchor 为随机初始化小预算，非 SOTA；M8 的 3c 臂 n=3 |
| P3 | **FP-aware 目标（惩罚 true-intergenic 的 genic 质量）是控制外溢的稳健杠杆，优于融合与学习型解码** | B | M4：FPLOSS spec 0.9303 vs FUSION 0.8615 vs CRF 0.8298±0.119 | LOG L456–L483 | 同预算 screen；后续 fp_lambda 2.5 在 GENERanno 有效 |
| P4 | **NT-v2 顶部层解冻能同时抬高 gbF1 与特异性，解除冻结特征的 gbF1 上限** | B | L0 0.8284→L4 0.8759→L12 0.9035（拟南芥单 seed）；M10 三 seed 干净植物均值 gbF1 0.8398、spec 0.9826 | LOG L597–L721、L793–L884 | 深度扫描单 seed 单物种 |
| P5 | **VAL-only 校准能把同物种池内的 FPR 压到 ≤0.01 而保持基因数**（工作点问题）| B（限定）| M11：FPR 0.0087、gbF1 0.8178、gcount 1.003（3 seed）；M19 s1 校准 gbF1 0.8593→0.8815 | LOG L886–L970、L1460–L1514 | 仅同物种池；跨物种失败（N1）；M22 重放显示收益已耗尽 |
| P6 | **评价尺子决定结论**：intergenic 必须是完整转录本 span 的补集；纯 specificity 可被"少预测"刷分；因此需双主轴 | A | ORF floor spec 0.8805（>anchor 0.8710）而 gbF1 0.3735；CRF-vec 在新尺子下由最优变最差（0.7138）| LOG L391–L420；FIND 2026-06-11 | 属方法学结论，在 yeast+fly 上发现，后以拟南芥/水稻再验（M9 起）|
| P7 | **gbF1/FPR 优秀并不意味精确基因结构** | A | M19 exact interval 0.0531–0.1498、pseudo-chain 0.0082–0.0123，释放基线 interval 0.8117–0.8882、chain 0.5850–0.7479（同 seqid）| R-M24 `report.md`；LOG L1703；FIND 2026-08-24 | M19 的 strand/phase 为占位符，不是公平的结构预测 |
| P8 | **释放基线在精确完整链上一致强于己方，排序稳定：ANNEVO > Tiberius > Helixer** | A（三个不同范围一致）| M24（测试 seqid）、M26（dev 染色体：0.7560/0.7143/0.6201）、M28 DEV（75.703/71.805/62.851）、M28 TRAIN（75.542/72.655/63.592）| R-M24；R-M26 `result.md`；PUB | 暴露不同：Helixer 训练列表含拟南芥/水稻；Tiberius 含水稻；ANNEVO 作者声明排除二者（`R-M28/pretrained_exposure.md`）|
| P9 | **联合、无参考的图解码能显著恢复结构信号**（相对于"整 block + phase 硬拒绝"）| B（解码机制不可单因素归因）| chain 0.3250→0.5303，interval 0.1204→0.7345，count 0.3253→0.9594；strand 0.9996、phase 0.9962、合法 6,188/6,188 | R-M25R-DEV `run_graph_terminal_20260926.md` | 同时 FPR 升至 0.0561，仍未达门槛；改进不能归因于 phase/跨 block/ORF/评分中的某一单项 |
| P10 | **统计与流程卫生**：≥5 seed 才可判 AXIS-1；3-seed anchor 不可靠；2-epoch 冒烟会误导；VAL-only 选参无泄漏 | A（经验结论）| anchor 0.8710→0.8436；冒烟 vs 5-seed 反转 | LOG L456–L507；FIND 2026-06-11 | — |
| P11 | **M25R 的训练在修复后数值有效，且 NO-GO 在大小写勘误后不变** | B | 3 个 checkpoint 张量有限；原 count 0.3253/0.1829/0.1034 未变 | R-M25R `terminal_summary.md`；`r4_case_correction_report.md` | — |

### 4.2 可靠负面结论

| # | 结论 | 强度 | 关键数值 | 证据位置 | 限定 |
|---|---|---|---|---|---|
| N1 | **固定的（单物种/少物种）NT-v2 L12 调用器无法跨物种**：拟南芥→水稻、→近缘 A. lyrata、→动物全部失败 | A | M12A 水稻 FPR 0.0311 / gcount 1.755（3 seed）；M13 A. lyrata FPR 0.0355；M14 gallus FPR 0.1156 / gcount 4.987；M16 聚合 FPR 0.0197 / gbF1 0.5615；M18-a 聚合 FPR 0.0444 | LOG L972–L1372；FIND 2026-06-17/18/19 | A. lyrata 为 scaffold 级；动物行存在预训练重叠；与基线的比较存在范围不一致（§6.5），但方向被 M24/M26 同范围测量印证 |
| N2 | **阈值/解码校准不能救跨物种失败，且在 GENERanno 上已耗尽收益** | A | M13 失败机理：水稻 oracle 网格 0/3 seed 可同时满足 FPR≤0.01、gcount≤1.25；M18-b gallus 最优有效点 gbF1 0.0060；M22 oracle：M19 s1 选点即最优 | R-M13_FAILURE_SANITY；R-M18 `report.md`；R-M22 `report.md` | oracle 诊断使用 test 标签，仅作上界性质判断 |
| N3 | **训练型 CRF 解码器无增益**（旧尺子下"赢"，新尺子下输）| B | M3b：CRF-vec spec 0.7138（最差）；M4：CRF spec 0.8298±0.119；M21（GENERanno）：best seed gbF1 0.8744/FPR 0.0192 vs M19 s1 0.8815/0.0065 | LOG L391–L420、L456–L483、L1546–L1590 | M21 仅 2 seed；"CRF"指本项目具体实现；3/3 tri-review 放弃（FIND 06-22）|
| N4 | **gene-body Tversky 辅助损失（固定设置）无增益** | C | gbF1 0.8635、FPR 0.02865、gcount 0.648，promotion gate=false | LOG L1592–L1645 | 单 seed、单固定设置；按预声明 gate 关闭 |
| N5 | **8 类多类 + CRF 不能补 gbF1；冻结 SegmentNT 特征封顶 gbF1** | B | M8：0.7189 < 0.7392，gcount 0.6625 | LOG L572–L595；FIND 06-12 | 仅干净植物、多类 n=5/3c n=3 |
| N6 | **GENERanno 开箱的 binary CDS mask 无法形成基因一致性，后处理无解** | B | gcount 4.219（官方）/3.844（我方后处理）| LOG L618–L653 | 仅拟南芥一条测试染色体 |
| N7 | **GENERanno 0.5B-base 不能靠加强 FP 目标救活** | C | FPR 0.0967、gbF1 0.6561（rice 0.4981）| LOG L1375–L1411 | 单 seed |
| N8 | **M25R 固定"结构头 + 约束 HMM 式解码"组合是开发期 no-go** | A（冻结门槛、完整 5,625 网格）| 0/5,625 可准入；count 0.0953–0.3305；epoch 越晚 count 越低 | R-M25R `terminal_summary.md`；D28 | 不隔离基模表示、监督、解码各自的贡献 |
| N9 | **仅取消 phase 硬拒绝不能达到 chain F1≥0.55** | A（严格上界）| 上界 0.5234608642 | R-M25R-DEV `phase_only_feasibility_bound.md` | 仅限固定 epoch-1/row601 |
| N10 | **run-graph R1（B）未达冻结开发门槛并触发 FPR 停止线** | A（冻结）| interval 0.7345<0.80；chain 0.5303<0.55；FPR 0.05606>0.020（>0.030 停止线）| R-M25R-DEV `run_graph_terminal_20260926.md` | 关闭范围仅此 checkpoint/评分/已测解码；不否定其他方法 |
| N11 | **同范围下 B 不是相对三基线的优势**（所有基线 chain F1 更高）| B | 0.5303 vs 0.6201/0.7143/0.7560 | R-M26 `result.md`、`comparison.tsv` | 基线为 2026-06-17 保存产物，暴露/上下文/预算未匹配 |
| N12 | **增加训练物种数不足以获得广谱固定调用器** | B | M16、M18-a 见 N1 | LOG L1198–L1372 | 面板小、单 seed |
| N13 | **（M28，简要）固定 M28 权重在 TRAIN/DEV 上均未领先三个原生基线，且多条候选机制按事前门槛关闭**（A6/A7、A9、A11、A13–A15、A21、A36、A37、A42、A43、FULL Q vs C0）| B | 见 §2.7 | R-M28 各 result；PUB | 另一 agent 负责，仅作衔接 |

### 4.3 需要限定或已被撤回的结论（读者勿再引用旧说法）

| 旧说法 | 状态 | 依据 |
|---|---|---|
| screen_anchor = 0.9213（三工具最大）| 改为"预训练上限，非门控"；同预算 anchor 0.5576；后来 intergenic 主轴的 anchor 0.8710 在 5-seed 复测降为 0.8436 | LOG L287–L332、L485–L507；FIND 06-10/06-11 |
| M3 "CRF-vec 获胜并可晋升 Track-B" | 被尺子修正推翻（改尺子后最差）| LOG L391–L420 |
| M7 "候选 specificity 优势跨 clade 成立（+0.155，鸡 +0.30）" | 鸡/果蝇受 SegmentNT 分割头微调污染；干净边际仅拟南芥 +0.06 | LOG L568；FIND 06-12 CRITICAL |
| M9 L4 "constrained postproc 管线 bug" | 实为 `constrained_gbF1 = gbF1 if FPR≤thr else 0` 的悬崖硬门；与双主轴 Pareto 哲学冲突（审计建议改软截断，**未见实施记录**）| FIND "指标审计 2026-06-14" |
| M12B "同面板基线" | 实为同物种、整基因组参考；与己方 test-seqid 分母不一致（M20 报告自述）；M24 起才在同 seqid 上重评 | R-M20 `report.md`；R-M24 |
| M12C "0.5B base 完全不可用" | M15 修正：较大预算下不再塌缩，但仍显著差于 1.2B-CDS | LOG L1147–L1195 |
| R-M20 `interval_overlap` 与聚合排名 | 已标注隔离：未按 seqid 过滤、"cds_base_*" 实为含内含子的 span | R-M20 `report.md` |
| R4 "约 36% emitted 结构非法 / 大多数参考与 canonical 文法不兼容" | 撤回：大小写 bug；大写后全部合法，canonical-compatible 94.03% | R-M25R-DEV `r4_case_correction_report.md` |
| M25R 验证网格中 `structurally_valid_complete_fraction = 1.0` | 无信息（"只要有预测即 1.0"的同义反复）| R-M25R `terminal_summary.md`；D30 |
| M26 前"拟南芥高 FPR 来自重复区/短 ORF" | 不能成立：83.19% 旧 FP 碱基落在 partial 编码基因区 | R-M26 `background_result.md` |

---

## 5. 评价体系与已知问题

### 5.1 数据与物种（配置层面可核实）

| 阶段 | 物种/accession | 用途 | 来源 |
|---|---|---|---|
| M1 pilot / Track A | S. cerevisiae GCF_000146045.2 (R64)；D. melanogaster GCF_000001215.4 | 工具复现与同预算 screen（低 UTR、基因密）| CFG `m1_data_manifest.yaml`；LOG L49–L117 |
| M7 | A. thaliana GCF_000001735.4 (TAIR10.1)；G. gallus GCF_016699485.2 (GRCg7b，子集 NC_≤20Mb) | 重锚；鸡后被标污染 | LOG L538–L570 |
| M8–M23 | A. thaliana、O. sativa GCF_034140825.1 (AGIS1.0)（"干净植物面板"）| 主开发面板 | CFG `m12_publication_panel.yaml` |
| M13/M16/M17 | A. lyrata GCF_000004255.2（v.1.0，scaffold 级，top-8 scaffold）；G. gallus；D. melanogaster；S. cerevisiae | 距离/动物/真菌诊断 | FIND "M13 diagnostic species frozen"；CFG `M18-MULTICLADE…yaml` |
| M24–M26 | 同拟南芥/水稻；开发 seqid NC_003074.8、NC_089041.1 | 结构精确度与解码 | CFG `M25R…yaml`；R-M26 `protocol.md` |
| 封存盲测 | **Setaria viridis** GCF_005286985.2 (Setaria_viridis_v4.0)，RefSeq 注释 RS_2025_03；仅 FASTA 已获取，注释 `reference_annotation: null` | M25 起冻结盲测物种 | CFG `M25R…yaml`；D26 |

### 5.2 切分与同源隔离

- **切分规则**：`assign_splits(seqids, val_mod=3, test_mod=4, mod=5)`：对排序后的 seqid 取 i%5，3→val，4→test，其余 train，≈60/20/20，**按整条 seqid 切**（`src/screen_anchor/data.py` L218–L225）。拟南芥 M10 起为 train 5/val 1/test 1；水稻 train 6/val 1/test 1（LOG L802–L805）。
- **验证染色体**：拟南芥 NC_003074.8；水稻 NC_089041.1（R-M25R-DEV `r4_terminal_review.md`；`R-M26 protocol.md`）。**M10–M19 的拟南芥 test 染色体为 NC_003075.7**（LOG L733，M24 范围 18.6 Mb）。
- **物种级切分**：M12A/M13/M14 为"拟南芥训练+校准，其他物种仅测试"；M16 加水稻；M18-a 加果蝇训练，测 A. lyrata/gallus/yeast（CFG `M18-MULTICLADE…yaml`）。
- **同源隔离**：本人在读到的 configs、results log、findings 中**没有发现任何序列同源（如 homolog 聚类、BLAST 去冗余）隔离**；M1 `check_data` 门槛只验证 seqid 不重叠（LOG L296）。M28 方法稿明确写出："DEV 染色体与训练染色体是**坐标不相交**，不是物种不相交，也未证明同源隔离"（`docs/M28_context_chain_methods_working.md` 第 101 行；`docs/M28_full_grammar_methods_working.md`）。拟南芥/水稻内的串联/片段重复使染色体级切分的泄漏风险未被量化（**这是本文的归纳，而非文件的明确声明**）。
- **验证/测试重用**：M9→M24 在同一组 val/test 染色体上反复做路线选择；M25R 起以 val 染色体作"开发范围"，M28 继续重复使用（PUB："原两条 DEV 已反复用于开发，不能升级为独立测试"）。
- **推断**：M10–M19 的"test"染色体（NC_003075.7 等）在 M25R/M28 中既不是 val 也不在 M28 的 6 条 TRAIN 列表中，可能处于闲置状态；但**文件未明确说明其被保留为封存集**，故不作结论。

### 5.3 sealed 数据规则

- **Setaria 注释禁令**（D26 §2）：训练、验证、阈值校准、分支选择均不得使用 Setaria；注释只能在 config、checkpoint、验证期解码参数、成功门槛、FASTA 哈希及 full/ablation 预测哈希全部不可变并记录后，一次性获取；`SETARIA_EMBARGO_RELEASED.json` 标记 + `eval_m25_structure.py` 校验；sbatch 在训练/推理前若发现任何 Setaria GFF/GTF/GBFF 即失败。
- 现状：M25R 无可准入 tuple，故无标记、无 Setaria 预测，**注释仍封存**（R-M25R `terminal_summary.md`）；Stage 1 的 `resolved_inputs.json` 记录 `setaria_files_read=false`；M28 文件亦称未放行。
- **盲测物种自身的暴露风险**（`R-M28/pretrained_exposure.md`）：Helixer 补充表把 *Setaria viridis* v2.0 列为作者 test；Tiberius 训练列表含 Sviridisv4.1。即使盲测数值未被打开，Setaria 对部分基线并非"物种未见"。
- 版本提示：docs 称 Setaria "RS_2025_03 / Setaria_viridis_v4.0"，而 Tiberius 列表是 v4.1，为同物种不同组装/注释版本，不能由此推断独立。

### 5.4 冻结阈值与门槛（按时间）

| 阶段 | 门槛 | 来源 |
|---|---|---|
| M1 合同 | screen：intergenic_FPR 仅作顾问（smoke/screen 阈值 ≤0.02）；full/scale：FPR≤0.01 硬门 + gene_count_ratio≤1.25 | LOG L38、L346–L353、L704（引述合同）；STAGE |
| 双主轴（06-11）| intergenic_specificity 严格 > anchor；gene_body_F1 ≥ floor 0.5276；macro_spec ≥ 0.7978 | LOG L415、L473 |
| M11 校准选择 | VAL：FPR≤0.01、gbF1≥0.70、gcount≤1.25，选 val gbF1 最大 | LOG L902 |
| M22 promotion gate | FPR≤0.01 且 gbF1>M19 s1(0.8815) 且 gcount≤1.25 | LOG L1630–L1635 |
| M25/M25R 准入 | FPR≤0.020、gene-count ratio 0.80–1.20、结构合法比例≥0.99；排名：pooled exact CDS-chain F1，tie-break：interval F1、低 FPR、枚举序 | CFG `M25R…yaml`；D27 |
| M25 Setaria discovery gate | interval≥0.80、chain≥0.55、gene≥0.50、strand≥0.98、phase≥0.90、合法≥0.99、FPR≤0.020、count 0.80–1.20、full−ablation interval/chain ≥+0.10、FPR 差≤+0.005；立即停止线：合法<0.95、interval<0.60、chain<0.30、gene<0.25、FPR>0.030、count 超出 0.50–1.50 | D26 §7；D27 |
| Stage 1（D29）| 硬诊断完整性：重放误差≤1e-5、6,450 参考与全部候选/预测对账、100% 合法性审计、模型/解码器不变、零 Setaria 读取；"dominant"=单一阶段解释 ≥50% 缺失链为**操作性启发式**；truth-assisted 增益 ≥0.15 同为启发式；曾建议的 interval≥0.60/chain≥0.50 被撤回 | D29 "Threshold taxonomy"；D30 |
| M28 | 事前冻结 DEV/TRAIN 门槛，闭合即停（如 A36/A37/A43"原五项继续条件"，FULL 的冻结 DEV 门槛）| R-M28 各 result |

### 5.5 基线对齐（same-scope）方式与"分母阶梯"

| 对齐步骤 | 做法 | 来源 |
|---|---|---|
| M1 | 三工具统一 `--span-mode cds`；per-species + base-weighted + macro 三视角 | LOG L203–L236 |
| M7 | 预训练基线评测**必须先把预测过滤到测试 seqid**，否则全基因组预测被计入测试染色体 intergenic 分母 → FPR>1（已修）| FIND 06-12 "Engineering" |
| M12B/M17 | 在拟南芥/水稻/A. lyrata/gallus/fly **整基因组**上评测释放权重；R-M20 指出参考基因数：己方 5,007/2,489，基线 33,461/13,324（A. lyrata：M16 测 3,609，M17 为 31,750）| R-M20；LOG L1286–L1304 |
| M24 | 基线按 M19 held-out seqid 过滤；同一 primary-transcript 策略、同一 FASTA 长度断言；不适用指标记 `not_applicable` 而非 0 | D25；R-M24 |
| M26 | 5 套**已保存预测**在同一开发染色体、同 6,450 参考链、同评价器上重评分；B 的四项冻结指标与旧 JSON 绝对误差 0；逐物种提供 frozen FPR 与"全注释 gene 背景"敏感性 | R-M26 `result.md`、`protocol.md` |
| M28 | 原生 raw-FASTA 重跑固定版本基线：ANNEVO 2.3.3/Magnoliopsida、Tiberius 2.0.8（在 2.0.7 依赖 SIF 内隔离执行）、Helixer 0.3.7/`land_plant_v0.3_a_0080`；primary 分母 7,728（DEV）/22,688（TRAIN）；exact-chain F1=2TP/(预测唯一链+参考 primary 链)；phase 另诊断，不入 exact-key | R-M28 `native_baseline_result.md`、PUB |

**分母阶梯（不可混用）**：
1. M1–M23：碱基级 CDS-span gbF1 / 完整转录本补集 FPR；
2. M24：M19 test seqid 上的 3,161+1,981 条 complete primary；
3. M25R/M26：6,450 条 complete-primary（拟南芥 4,151 + 水稻 2,299，dev 染色体）；
4. M28 DEV：7,728 条 primary（5,437+2,291，同两条 dev 染色体，不同 primary 定义——文件称"旧 6450 条 parent-primary 指标不是本表口径"，`R-M28/a31_local_path_caller_dev_result.md`）；
5. M28 TRAIN：22,688 条（拟南芥 17,739 + 水稻 4,949）。

### 5.6 指标定义要点

| 指标 | 定义 | 注意 |
|---|---|---|
| gene-body F1（CDS-span，unconstrained）| 碱基级 CDS(+start/stop) 转录本 span 的 F1 | M1–M23 主指标 |
| constrained gene-body F1 | = gbF1 if FPR≤阈值 else 0 | 悬崖硬门，审计建议软截断，未见实施（FIND 06-14）|
| intergenic_FPR / specificity | genome − **完整转录本 span（含 UTR）**；specificity=1−FPR；base-weighted 为 headline，macro 防跨物种掩盖 | 06-11 修订；M26 起显示对参考背景高度敏感 |
| gene_count_ratio | 预测 gene 数 / 参考 gene 数 | 06-10 修正分母（原 transcript 比被降为 diagnostic）|
| exact CDS interval F1 / exact CDS-chain F1 | 同 strand、坐标完全相同；chain 要求有序 CDS 链完全相同 | M24 起；phase 另算；"coding-gene F1"在冻结实现中同 chain 定义，非独立证据（R-M25R-DEV `run_graph_terminal`）|
| frozen background FPR（M25R/M26）| 预测 CDS span 落在 complete-primary 参考 exon span 补集内的碱基比例 | 非 gene-level 分类 FPR，不认证生物学假阳性 |

### 5.7 已知评价问题清单

1. **GENERanno 预训练物种重叠（overlap_unknown）**：公开资料没有完整物种/accession 排除清单（D28、FIND 06-20/06-21、R-M21 provenance scout）。2026-10-13 的 M28 来源更新：Carbon-A 披露 GENERanno-1.2B-CDS-preview（revision b0483c23…）的 RefSeq CDS 监督，并把两个项目 assembly（拟南芥 GCF_000001735.4、水稻 GCF_034140825.1）列为训练快照参考；实际优化窗口、验证划分、同源隔离仍未知（`R-M28/pretrained_exposure.md`；PUB 2026-10-09 段）。"overlap_unknown"不得改写为"已确认重叠"或"已确认排除"（D25 §1）。
2. **NT-v2 "干净来源"仅基于模型卡声明**（植物被排除于 850 基因组预训练）；SegmentNT 分割头在人类 + {mouse, chicken, fly, zebrafish, worm} 微调 → 鸡/果蝇不是干净 held-out（FIND 06-12 CRITICAL）。
3. **M12B/M17 范围不一致**（§5.5）：M16/M18 与 M17 的"己方 vs 基线"比较、M19 与 M12B 的 gbF1 差距（0.88 vs 0.925）不是同范围比较；M24 同 seqid gbF1（如 Tiberius 拟南芥 0.9571、水稻 0.8842）仍显示基线更强，故方向不变，但数值不应直接引用旧表。
4. **验证/测试重用与 anchor 不稳**：同一 val/test 染色体被重复用于选路；3-seed anchor 0.8710 vs 5-seed 0.8436；ACTIVE_GOAL 中 screen_anchor 仍为 3-seed 估计（LOG L502）。
5. **半成品基准**：`ACTIVE_GOAL.status=draft`，`sota_benchmark.value=0.0` 为占位，M2（ANNEVO 兼容的完整评价）从未完成（STAGE "当前不确定性集中在 claim 合同"）。
6. **参考背景政策**：拟南芥 dev 染色体上 83.19% 的旧 FP 碱基落在 gene 层标 partial 的编码区（面积 2,300,148 bp vs 水稻 29,821 bp），故原 FPR 的物种差异含参考政策效应（R-M26 `background_result.md`）。complete-primary 参考选择决定了"负区域"。
7. **工具链产生的诊断误导**：R4 大小写 bug、M25R 同义反复的合法比例、R1/R3 诊断 bug、M10 sbatch `--run-status COMPLETED` 误判 failed_run、validate_goal profile 回归（FIND 06-14、06-15）。
8. **环境事故**：2026-06-13 Write/Edit 工具回显成功但未落盘，导致两批 job（8575441-49、8623290-92）以错误配置 TIMEOUT，且 scancel 回显为假（FIND 2026-06-13 Engineering Finding）。
9. **基线版本/暴露不均**：基线是 2026-06-17 保存产物，上下文与训练预算未匹配；Helixer 含两物种于训练；Tiberius 实际 2.0.8 源码在 2.0.7 SIF 中隔离运行（`R-M28/native_baseline_result.md`、`pretrained_exposure.md`）。
10. **M26 R3 与研究路线图缺失**：R3 仅见协议；PUB 所称权威文档 `research_roadmap_20260926.md` 不在仓库。

---

## 6. 发表策略与审稿风险（文档中曾记录的内容）

### 6.1 目标与叙事的演变（PUB、D28–D30）

| 时间 | 核心叙事 / 主张 | 状态 |
|---|---|---|
| 2026-06-17（M12 时代）| "经校准的 NT-v2 内含子感知注释器在共享评估面板上以更低 intergenic FP 与稳健 gene-body 召回优于现有调用器"；贡献菜单 C1（固定校准 M9-L12 跨物种）、C2（同面板更实用）、C3（CDS 特化 + FP-aware 适配不是任意预训练都行）、C4（VAL-only 校准是部署旋钮，不是测试集调参）| C1 被 M12A/M13/M14/M16/M18 否定；C2 重框为"实用/适配权衡"；C3 保留但受规模混淆与来源限制；C4 保留（PUB §2）|
| 2026-06-20 | 最强路线从 NT-v2 转向 GENERanno 1.2B 适配（M19：0.8421/0.0083/1.083、0.8815/0.0065/0.830），写作为 adaptation/comparability，不写 clean held-out SOTA | PUB §1 |
| 2026-09-01 | 路线改为 raw-FASTA→直接结构注释器；目标 Nature Communications 及以上；"当前没有支持投稿级优越性或独立泛化的模型结论" | PUB 开头；D28 |
| 2026-09-26 | 固定 M25R 路线 NO-GO；M26 同范围比较 B 低于三基线；"仍寻求可发表贡献，但 NO-GO 与负结果报告都不完成该目标"；旧 C1–C4、骨干规模比较、投稿准备度勾选不得再作为新颖性/因果/干净泛化/投稿准备的证据 | PUB "2026-09-26 superseding status" |
| 2026-10-09 | M28 全语法模型；原 TRAIN 同尺对标：M28 两权重均未领先三基线；创新性边界受 Tiberius 可微 HMM、Carbon-A 等先例约束；"性能优越、机制因果、独立泛化、投稿准备度均不能据此宣布成立" | PUB 2026-10-09 段 |

### 6.2 ChatGPT Pro 系统评审要点（D30，两轮；2026-09-01）

- **范围声明**：Pro 报告 ACCESS_VERIFIED（私有库 `JiWang-Unige/coding-rna`，commit `fe38a675…`），但检查点与完整 5,625 网格不在 GitHub 快照中，Pro 只能审计代码/配置/摘要，**不能独立重跑**；未读 Setaria 注释。
- **联合判决表**：积极信号仅限"粗粒度低 FP 适配"；M24 改变科学目标；M25 无效；M25R 有效 no-go；项目处于"structural feasibility diagnosis"，**不具备 Nature Communications 性能主张**。
- **关键解释**：5,625 个 tuple 全部 FPR 合格而 count 全部失败，说明丢基因而非过预测；`validation_candidates()` 先在边界阈值为 0 时跑完整解码，再由 `_filter_models()` 只能"删"不能"补"，故阈值网格不可能救回召回；训练 loss 下降而复原数下降，提示目标与完整基因脱节；有效合法比例统计是同义反复。
- **纠偏（第二轮）**：区分硬完整性要求与操作性启发式；撤回 interval≥0.60/chain≥0.50 的发展期门槛；"强基线 +2 个百分点"仅为规划目标，非期刊规则。
- **决策树（Stage 2）**：按 Stage 1 归因选一条——候选缺失→SegmentNT 原生长上下文发射探针；转移/motif 高而边界/相位头弱→冻结 backbone/LoRA/region 头只训边界/相位头；训练强验证弱→物种均衡拟合；候选与 truth-assisted 上界高而终端低→只改一条解码规则；无单一主因→重新设计任务表示，不再给 M25 打补丁。**之后的 M25R-DEV 勘误与 run-graph R1（B）实际沿着"解码/组装"分支推进并以 NO-GO 收尾。**
- **允许的 claim 措辞（来源未知时）**：采用"适配公开发布的 CDS 特化 GENERanno checkpoint，并评估对任务训练/校准中被排除物种的下游泛化；上游预训练与 CDS 后训练语料的物种/accession 构成未完全披露，因此本实验不确立 pretraining-held-out、no-overlap 或 zero-shot 泛化"。需删除的说法：clean pretraining-held-out、zero-shot、no-overlap unseen species、Setaria unseen to the backbone、"通用基因语法无先验暴露"、把 1.2B vs 0.5B 全部差异归因 CDS 后训练。
- **投稿所需证据（D29 Stage 9 / D30）**：开发期通过的精确结构模型（事前门槛）；冻结候选 + 一次性 Setaria 盲测；Helixer/Tiberius/ANNEVO 在相同范围与评价器；最终架构 3 seed + ≥6 个 held-out 物种；匹配机制消融（backbone/heads/decoder/context）；exact interval/chain/gene/strand/phase/FPR/count/独立合法性；长基因、tile 边缘、非 canonical 剪接、重复区分层；独立 RNA/蛋白同源/共线性证据；运行时间与复现；来源闭合或收窄为 adaptation claim。
- **No-go 条件（D29）**：只提升粗粒度 gbF1；精确结构仍远低于释放基线；收益依赖目标标签调参；来源阻碍预期迁移主张。允许的替代：跨物种配对优势并有不确定性支撑且无不可接受的 FPR/count 代价，或同等强度的预注册 Pareto 价值（更低 FP、长基因恢复、运行时间、独立支持的参考修正）；否则只可能是"较窄的 adaptation/methods paper"，可能对应不同期刊。
- **Pro 建议删除的低价值工作**：继续粗 gbF1 阈值调参、GENERanno 0.5B 性能开发、通用 CRF/Tversky 扫描、过早加 seed、过早 MoE/长上下文、Setaria 访问、细胞器与 unplaced scaffold、替代 isoform、以 BUSCO 代替精确结构。

### 6.3 审稿风险预演（PUB §5，原文要点）与当前状态

| 预期审稿意见 | 文档中的预防证据 | 当前状态（本文评估）|
|---|---|---|
| 仅因 train/val/test 同物种池才有效 | M12A 固定模型跨物种评价 | M12A 已失败 → 该主张不可用（N1）|
| 与 Tiberius/Helixer/ANNEVO 不可比 | M12B 同面板协议、原始命令/版本 | M12B 范围不一致（§5.5），已由 M24/M26/M28 同范围重评分修正，且基线更强 |
| 预训练重叠污染 held-out | M12/M19 来源审计 | GENERanno 仍 overlap_unknown，M28 来源显示 assembly 属 Carbon-A 训练快照参考；"clean"主张不可用 |
| 任何预训练模型都行 | M18 1.2B vs 0.5B、M19 稳定性 | 与规模混淆（D30）；仅为适配证据 |
| FP 校准是测试集调参 | VAL-only 选择规则、一次性应用 | 成立，但 M22 表明收益已耗尽、M26 表明 FPR 对背景政策敏感 |

### 6.4 本文补充的潜在审稿风险（归纳，非文件原文；均有文件依据）

1. **新颖性**：PUB（2026-10-09 段）与 `docs/M28_complete_chain_competition_prior_work_20261003.md` 记录——可微 HMM 端到端训练（Tiberius）、GLM/神经表示 + 结构解码、结构约束、判别式结构竞争（mGene 2009）、Carbon-A（transition TCN + Viterbi 语法 + confidence TCN）均有先例；R-M26 亦指出"深度学习加结构解码…已被包含"。
2. **负结果压倒正面**：M25R/R1 B/M28 各固定权重均低于三基线；PUB 要求论文核心仍须"更好的完整链性能，或在具备独立资格的更广范围上的实证收益"。
3. **开发数据重复使用**：DEV 不是独立测试，任何 DEV 进展不能写成泛化（PUB）。
4. **基线暴露不均与 Setaria 暴露**：见 §5.3、§5.7。
5. **评价口径敏感性**：M26 R2 显示同一预测在不同参考背景下 FPR 变化一个数量级（拟南芥 12.74%→1.31%），工具排序在背景改变后会反转（Tiberius 变得更优）→ 审稿人会质疑 FPR 作为主终点。
6. **单 seed/少 seed**：多数正向 screen 为单 seed；M4 本身示例说明 2-epoch 单 seed 会反转结论。
7. **计算与可复现**：Pro 无法重跑检查点；本地仅保留 M28 outputs，M1–M26 原始 metrics 不在本地。
8. **来源声明基于模型卡**（NT-v2 排除植物）而非独立审计。

---

## 7. 对 M28 的衔接要点（供接手者）

- M28 把 M24–M26 的结论变成设计前提：主终点为 exact-chain F1，分母不剪枝，phase 与合法性另查，"更高 CDS-base F1 ≠ 更高完整链 F1"（如 LOCAL 9216 步：pooled CDS-base F1 90.360%→94.097%，但完整链 F1 45.590→41.729；PUB）。
- 继承的方法学纪律：事前冻结门槛，失败即关闭，不加 seed/不调阈值；诊断脚本的数值核验（M25R-DEV 经验）；撤回错误解释时保留原文件并出勘误报告。
- 未继承/已放弃：M25R 的"结构头 + HMM 式阈值网格"；GENERanno 上的学习型 CRF/Tversky；NT-v2 与 SegmentNT 路线（不再作为主线）。
- 仍悬而未决：GENERanno 暴露与同源隔离；独立物种面板与 Setaria 的放行；背景/可评价区域定义对 FPR 的影响（M26 R2 留下的问题）；M28 全语法模型相对 Tiberius/ANNEVO 的差距。

---

## 附录 A：关键里程碑—日期—作业/ID 速查

| 里程碑 | 日期 | Slurm job / 标识（取自 LOG/EV）|
|---|---|---|
| M1 PILOT | 2026-06-10 | 8528176；Helixer 8530344；ANNEVO 8537422 |
| M1 same-budget anchor | 06-10 | 8538949–8538954 |
| FP-SEGMENTNT-PROBE | 06-11 | 抽取 8548459；训练 8548460–62 |
| M4 | 06-11 | 8550151–8550166 |
| M5 / M6 | 06-11 | 8551173–77、8551128–29 / 8552452–56 |
| M7 | 06-12 | 8554530、8555770–73；ANNEVO 上限 8554546 |
| M8 | 06-12 | 8559000–09、8558997–99 |
| M9 / M9-DEEP | 06-14 | 8667188–90 / 8751498–8751500 |
| M10 | 06-15 | 8833071[0–2]；GENERanno smoke 8833070 |
| M11 | 06-16 | 8934130[0–2] |
| M12A / M12B / M12C | 06-17 | 8974902 / 8975814 / 8974903 |
| M13 / M14 / M15 / M16 | 06-17–18 | 9019532 / 9022700 / 9023295 / 9065776 |
| M17 / M18 / M19 | 06-19–20 | 9119473 / 9123661、9131867、9122868 / 9141356 |
| M20 smoke / M21 / M22 | 06-21–06-30 | 9249721 / 9259965、9298703 / 9570144、9570173 |
| M23 | 07-01 | 9854668 |
| M25 / M25R | 08-27–08-31 | 12094731（无效）/ 12116383 |
| M25R Stage1 R1–R4 | 09-02–09-08 | 12128521（R1 失败）/ 12166448（R3 失败）/ 12270926（R4 完成）/ 勘误 12496928 |
| phase-only / block assembly | 09-08 | 12499507 预检；12520954 |
| run-graph R1 | 09-26 | GPU 12522183（评价工程失败）；CPU 恢复 13225219 |
| M26 R1/R2 | 09-26 | 13225696 / 13225728 |

## 附录 B：证据源文件清单（最常引用）

- `docs/06_results_log.md`（L49 PILOT、L203 span 统一、L287 anchor、L334 M3、L391 换尺子、L456 M4、L538/L568 M7 与泄漏标注、L572 M8、L597 M9、L618 GENERanno 开箱、L655 M9-DEEP、L793 M10、L886 M11、L972 M12、L1030 M13、L1095 M14、L1147 M15、L1198 M16、L1261 M17、L1318–L1457 M18、L1460 M19、L1546 M21、L1592 M22、L1647 M23、L1703 M24、L1711 M25）
- `docs/10_findings.md`（06-10 至 09-01 各条；"指标审计 2026-06-14"；"CRITICAL SegmentNT 泄漏 2026-06-12"）
- `docs/12_publication_strategy.md`；`docs/25`–`docs/30`
- `reports/M24-DIRECT-STRUCTURE-DIAGNOSTIC/report.md`
- `reports/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/terminal_summary.md`
- `reports/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/`：`r3_terminal_summary.md`、`r4_terminal_review.md`、`r4_case_correction_report.md`、`phase_only_feasibility_bound.md`、`block_assembly_result.md`、`run_graph_decoder_contract.md`、`run_graph_terminal_20260926.md`、`nonexact_forensics_result.md`
- `reports/M26-SAME-SCOPE-MECHANISM/`：`result.md`、`background_result.md`、`protocol.md`、`comparison.tsv`
- `reports/M13_FAILURE_SANITY/report.md`；`reports/M18-MULTICLADE-ORACLE-CALIBRATION/report.md`；`reports/M19-COMPARABILITY-EVIDENCE/comparison_tables.md`；`reports/M20-SOTA-ERROR-ANALYSIS/report.md`；`reports/M21-CLEAN-BACKBONE-PROVENANCE-SCOUT/report.md`；`reports/M22-NONCRF-CLAIM-ROUTE-DESIGN/report.md`
- `reports/PROJECT-STAGE-REPORT-2026-08-16/report.html`
- `configs/M10…`、`M18-GENERANNO-1P2B-SPEC-OBJECTIVE.yaml`、`M19…`、`M21…`、`M22…`、`M23…`、`M25R…`、`m12_publication_panel.yaml`、`m1_data_manifest.yaml`
- `src/screen_anchor/data.py`（`assign_splits`）
- M28 衔接：`docs/M28_full_grammar_{abstract,methods}_working.md`、`docs/M28_context_chain_methods_working.md`、`reports/M28-METHOD-RESTART/{pretrained_exposure,native_baseline_result,a31_local_path_caller_dev_result,full_grammar_raw_train_result,full_train_result}.md`
