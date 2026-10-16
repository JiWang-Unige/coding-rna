# M28 当前阶段（方法重启）现状梳理

梳理日期：2026-10-16。只读梳理，未改动任何既有文件；本文件是唯一新增产物。
证据快照：最新项目文件为 2026-10-14 07:11 的 `docs/M28_full_grammar_manuscript_working.md`，之后没有新实验结果。
所有数值均注明来源文件（相对项目根 `/Users/jiwang/Desktop/coding RNA`）。F1 均为"同 strand、全部有序 CDS 区间完全一致"的 exact-chain F1，单位 %，pooled 指合并 TP/预测/参考计数后计算。

---

## 0. 一页结论

1. **M28 的方法主线**：先后经历两个阶段。
   - 阶段一（A3–A44，2026-09-29 至 10-07）：冻结 GENERanno 编码器，上面接 query/pointer 候选生成器加链评分头，在"有限候选池"里选链。
   - 阶段二（full-grammar，10-07 起）：放弃候选池，改为"15 角色逐碱基打分 + 完整 CDS 语法下的全集合配分函数 CRF 训练 + 整条 contig 的精确 MAP 解码"。
2. **当前最好的结果，没有任何一个同时满足"独立数据、单模型、超过 Tiberius/ANNEVO"**：
   - DEV（拟南芥 + 水稻各 1 条染色体，7,728 条 primary 参考，两条染色体已被反复使用）：最佳是 fresh-joint，pooled 70.107，低于 Tiberius 71.805 和 ANNEVO 75.703，高于 Helixer 62.851。
   - TRAIN（6 条训练染色体，22,688 条参考，也是头部的拟合数据）：rank-4 q/v LoRA 为 72.958，略高于 Tiberius 的 72.655，低于 ANNEVO 的 75.542。LoRA 没有跑过 DEV。
   - Carbon-A（同一 GENERanno 主权重的外部完整流程）在 TRAIN 上为 78.172（过滤）/ 76.063（不过滤），高于我们所有变体，也高于三个原生基线。
3. **没有 Setaria 盲测，也没有任何独立物种或独立面板评测**。Setaria viridis 的注释至今未被读取（`docs/28_current_research_state_2026-09-01.md`，各 M28 报告的"sealed/Setaria 未访问"）。
4. **没有多 seed，没有置信区间**。所有结果都是单 seed、单次拟合。消融只有 TRAIN 上的单 seed 小因子设计；DEV 上只有 6 个固定候选的单次比较。
5. **论文草稿措辞极为克制**（几乎每段都有"不建立泛化/优越性"的限定），明确过度的地方不多，但有结构性问题（见 §4.2）。
6. **最后停在**：`rank-4 q/v + o_proj LoRA` 这一对照（plan 与 7 个 sbatch 已写好，2026-10-14 04:59 UTC 因集群队列拥堵——运行 102 个、等待 1038 个，超过预设并发上限 80——尚未提交任何作业）。此外，E/M 能量池与最小能量池只做过 TRAIN 缓存重解码，没有整染色体新编码评估。

---

## 1. 数据划分与口径（先把"TRAIN/DEV"讲清楚）

| 名称 | 内容 | 用途 | 来源 |
|---|---|---|---|
| TRAIN | 拟南芥 NC_003070.9 / NC_003071.7 / NC_003076.8 + 水稻 NC_089043.1 / NC_089044.1 / NC_089046.1，共 156,779,632 bp；22,688 条 primary 参考（拟南芥 17,739、水稻 4,949） | M28 头部的拟合数据（1,491 个唯一窗口、1,536 draw × 3 遍）；同时作为"原 TRAIN 同尺对标" | `docs/M28_full_grammar_manuscript_working.md` Results；`docs/12_publication_strategy.md` |
| DEV | 拟南芥 NC_003074.8（23,459,830 bp）+ 水稻 NC_089041.1（29,936,421 bp），共 53,396,251 bp；7,728 条 primary（5,437 + 2,291） | 开发评测；从 M20/M26 起反复使用，**不是独立测试** | `reports/M28-METHOD-RESTART/full_grammar_dev_baseline_result.md` |
| 盲测 | *Setaria viridis* GCF_005286985.2 | 冻结的盲测物种，**注释从未被读取** | `docs/28_current_research_state_2026-09-01.md` L17–20、L74 |
| 独立迁移候选 | Brachypodium、Brassica rapa、番茄、Amborella（仅登记元数据） | **未冻结、未执行**；Brassica 已知在 Tiberius 训练集中 | `reports/M28-METHOD-RESTART/independent_transfer_candidates.md` |

**参考口径**：M26 的 "CDS-assessable 最长 eligible CDS，平分取 transcript ID"，primary 链 7,728（DEV）和 22,688（TRAIN）。
这与 M20/M26 旧口径（6,450 条 complete-primary）不同，**不可混用**。

**基线**：ANNEVO 2.3.3/Magnoliopsida、Tiberius 2.0.8/angiosperms、Helixer 0.3.7/land_plant_v0.3_a_0080，均为 2026-10-03 的原生运行，同 FASTA、固定默认操作点。

**暴露边界**（`reports/M28-METHOD-RESTART/pretrained_exposure.md`）：

- GENERanno 1.2b-cds-annotator-preview 是带 CDS 注释监督的 checkpoint，不是纯无监督基座。
- Carbon-A 技术报告把拟南芥/水稻 assembly 列入其训练快照参考子集，但实际优化器窗口与同源隔离未核实。
- Helixer 当前权重明确含这两个物种。
- Tiberius 含水稻（Kitaake），不含拟南芥。
- ANNEVO 作者声明排除了拟南芥和水稻。
- 因此"输入只有 raw FASTA"不等于"所有方法都没见过该物种"。

---

## 2. 任务 1：M28 子实验谱系表

判定符号：**go** = 通过事前门槛并推进下一步；**no-go** = 事前门槛失败并关闭；**诊断** = 无门槛、只解释机制；**对标** = 与外部方法比较；**工程** = 接线或数值核验，不是科学判定。
TRAIN 门槛里"五项继续条件"通常是：两物种正确链各自严格增加、两物种已知非匹配各自不增加、合并 ≥5-exon 正确链严格增加。

### 2.1 阶段〇：重启与基线（2026-09-29 至 10-03）

| 项 | 假设 / 改动 | 结果 | 判定 | 来源 |
|---|---|---|---|---|
| M25R/M26（前史） | 固定结构头 + 约束解码 | M25R B 在 6,450 同尺上 exact-chain F1 0.5303，ANNEVO 0.7560、Tiberius 0.7143、Helixer 0.6201；基因数恢复到 0.959，但 interval/chain/FPR 三门槛均未过，FPR=5.6% 触发 >3% 停止条件 | **no-go**（固定路线） | `reports/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_terminal_20260926.md`；`reports/M26-SAME-SCOPE-MECHANISM/result.md` |
| M26 背景敏感性 | 把背景从"参考 exon 补集"改为"全部已注释 gene span 补集" | B 的拟南芥 FPR 从 12.74% 降到 1.31%；各工具排名受参考政策强烈影响 | 诊断 | `reports/M26-SAME-SCOPE-MECHANISM/result.md`、`background_result.md` |
| M20 报告（有未提交修改） | 把表头改成"混合评价范围的历史数据，不是同尺比较"，并说明 LoRA 行参考基因数 5,007/2,489 对基线 33,461/13,324 | 旧 gbF1 表（ANNEVO 0.9269 等）**不能作为同尺结论** | 口径更正 | `reports/M20-SOTA-ERROR-ANALYSIS/report.md`（git diff +5/-3） |
| A3 | 冻结 GENERanno（四块 6,144 bp 编码）+ query/pointer 生成器，649,457 参数，1,536 draw × 3 遍 = 4,608 步，贪心 + beam4 产生自然候选池 | DEV pooled **52.320**（拟南芥 51.757、水稻 53.524） | 基线 | `reports/M28-METHOD-RESTART/a31_local_path_caller_dev_result.md` 表 |
| A5/A8/A9/A11/A12 | 在 8 个压力窗、27 对参考/窗上定位失败：正确前缀排名、完成池保留、前缀 critic、组件均值、质量评分 TRAIN 诊断 | A8 少命中 4 个；A9 critic 使自由整链可达性 early 4 对 late 15 下降；A11 15/27 不变；A13 表显示 A6 评分下正确链仅 1237/2916 正 gain | **no-go / 诊断** | `a5_true_prefix_diagnostic.md`、`a8_quality_complete_train.md`、`a9_query_prefix_result.md`、`a11_component_mean_result.md`、`a12_quality_train_diagnostic.md` |
| A6 / A7 | 在 A3 生成器上加 165,635 参数整链质量评分头（A7 为跨窗 max 聚合） | DEV pooled **47.925 / 48.533**；precision 高（A6 约 78%）但 recall 仅 34.5% | **no-go** | `a6_a7_quality_result.md` |
| FULL（Q vs C0） | 完整 TRAIN 25,512 窗、三遍 76,536 步：Q（A3 目标）对 C0（局部头 + 原生 HMM 解码，366,191 参数） | DEV pooled Q 56.28 / C0 55.66；Q 在水稻更高、拟南芥更低且背景增加，原门槛失败 | **no-go** | `full_train_result.md` |
| A10 | 四块 6,144 bp 编码改为一次 24,576 bp 整窗编码 | DEV micro 57.160（A3 52.320），两物种 F1 都升，拟南芥背景 84,489 > 54,809 bp 上限 | 门槛 FAIL，但是"整窗上下文有益"的第一个信号 | `a10_context_result.md` |

### 2.2 阶段一：候选池 + 链评分（A13–A44，2026-10-03 至 10-07）

| 项 | 假设 / 改动 | 关键数值 | 判定 | 来源 |
|---|---|---|---|---|
| A13 | 完整链竞争目标（null + 已知严格冲突竞争） | 8 条件 6 过 2 败：正确链正 gain 1237→2039，但 known 非匹配正 gain 311→1018 | **no-go** | `a13_chain_competition_result.md` |
| A14 | 兼容集合精确似然替代 BCE | 五条件 3 过 2 败 | **no-go** | `a14_chain_set_result.md` |
| A15 | 同上，限定自然候选域 | 3 过 2 败（正确链和长链增加，非匹配也增加） | **no-go** | `a15_natural_chain_set_result.md` |
| A16 | 诊断：合法事件概率质量 | start 97.8%、next-start 99.6%、end 联合 82.9%（≥5 exon 仅 71.2%）；入口非法概率不是主因 | 诊断 | `a16_legal_support_result.md` |
| A17 | 局部注意力生成器 | 三项生成覆盖门槛全败（拟南芥正确候选 1751→1690） | **no-go** | `a17_local_context_result.md` |
| A18/A19/A20 | 诊断/准备：事件均值目标审计、局部监督覆盖、部分观测链集合似然（未知链边缘化）的数学可行性 | 无性能结果；A20 指出"不当 0 标签"≠"零梯度" | 诊断 | `a18_*`、`a19_*`、`a20_*` |
| A21 | 剪接连接整链头（248,068 参数） | 3 过 2 败 | **no-go** | `a21_junction_chain_result.md` |
| A22 | 诊断：链选择位移 | 选择器无错；536 个正确链在自然池中根本不存在，下一优先级是候选生成 | 诊断 | `a22_chain_displacement_result.md` |
| A23 | 相对边界生成器（+160 参数） | 三项生成覆盖全过（正确候选 2380→2415，≥5 exon +44），但 known 非匹配 +2461 | 通过覆盖门槛，**代价是候选池膨胀** | `a23_relative_boundary_result.md` |
| A24/A27/A28/A29 | 诊断：536 个遗漏链为何缺失 | 缺首起点 246/536；局部必要条件排除 315、局部通过却仍丢 221；搜索淘汰分散在入口池/延伸池/完成池/分数界；其中 320/1708 的 cutoff 不能当负例 | 诊断 | `a24_*`…`a29_*` |
| A25 | 首起点提示 | 三项全败（拟南芥正确候选 1751→1698，非匹配 +4470） | **no-go** | `a25_start_weighted_result.md` |
| A26 | 相位证据 | 3 过 2 败 | **no-go** | `a26_phase_evidence_result.md` |
| A30 | 冻结搜索保留损失 | 五条件 **0/5** | **no-go**（A32 因此不执行） | `a30_search_retention_train_result.md` |
| **A31** | **不再训练新评分头，直接对 A3 的 15 角色局部 logit 沿候选链求和（单位权重、null=0），候选池取 greedy+beam4 联合池，再做同链向非重叠选择** | TRAIN 五条件全过（拟南芥 proxy F1 59.8→75.6）；**DEV A31_A3 66.167，A31_A23（A23 生成器）67.813**，高于 Helixer，低于 Tiberius 3.99、ANNEVO 7.89 | **go**，阶段一最佳 | `a31_local_path_caller_train_result.md`、`a31_local_path_caller_dev_result.md` |
| A31 诊断 | 评分与候选池因子分解（DEV） | greedy+query 52.32 → greedy+local 56.40 → union+local 66.17 | 诊断（有序条件，不是因果份额） | `a31_score_pool_dev_diagnostic_result.md` |
| A32/A33 | A30+A31 组合；发射分数驱动搜索 | A32 因 A30 失败未执行；A33 在指针重放检查 FAILED（浮点差），五条件未评价 | 未完成（工程失败） | `a32_search_local_caller_protocol.md`、`a33_emission_beam_execution.md` |
| A34/A35 | 诊断：把正确几何直接给固定评分 | 536 遗漏中 190 个 gain≤0、179 个正分未选、167 个被选；非正链的最负贡献来自 coding（91）和 intron（99）项 | 诊断 | `a34_missing_path_emission_result.md`、`a35_saved_role_contributions_result.md` |
| A36 | 从头训练整链监督的发射模型 | 五条件 2/5（拟南芥正确链 1531→1138） | **no-go** | `a36_structured_emission_result.md` |
| A37 | 联合上下文整链模型（612,324 参数） | 2/5（1531→1396） | **no-go** | `a37_context_chain_result.md` |
| A38 | 同架构，逐候选空集 BCE | 2/5（1531→1494） | **no-go** | `a38_logistic_context_protocol.md` |
| A39 | 同架构，BCE + 兼容集合 NLL 各半 | 2/5（1531→1496；合并长链 653） | **no-go** | `a39_null_anchored_set_protocol.md` |
| A40 | 局部角色证据输入 | 3/5（水稻 +1 不抵拟南芥 −47） | **no-go** | `a40_role_evidence_context_result.md` |
| A41 | 冻结 A3 发射上加零初始化线性修正（1,806 参数） | 2/5 | **no-go** | `a41_emission_residual_protocol.md` |
| A42 | **冻结 A3 发射，用完整 CDS 语法的精确 DP 搜索替换 query/beam 候选池** | 五条件 0/5：正确配对 2009→981，≥5-exon 700→81 | **no-go**，也是转向的关键证据 | `a42_emission_lattice_train_protocol.md` 开头 |
| A43 | 自然竞争监督 | **0/5：两物种完整正确链均为 0** | **no-go** | `a43_result.md` |
| A44 | 局部监督 + 有限结构项配对拟合（local-only 对 local+struct，366,191 参数，N 为 valid_bases） | 本地没有结果文件；后续文档统一称"旧 dense/A43/A44 失败""已停止的 A44 链" | 关闭（**结果文档缺失**，见 §5） | `a44_paired_fit_protocol.md`；`full_grammar_first_fit_contract.md` L29；`full_grammar_raw_train_contract.md` L9 |

**阶段一的整体教训**（`docs/12_publication_strategy.md` 10-07 段落）：
- A31 之后的 A36–A43 没有一个超过 A31_A3。
- A42 说明：对 A3 现有发射做"精确语法解码"反而崩溃（2009→981）。训练目标（有限候选 BCE）和推理空间（整条 contig 的语法）不一致。
- 这是转向"训练与推理同一语法"的 full-grammar 的直接动机。该因果链是文档推断，不是实验证明。

### 2.3 阶段二：full-grammar 变体（2026-10-07 起）

先说共同底座：C0 头 + 全语法配分函数训练 + 整 contig 精确 MAP。下表每行是"与某个对照比，只改一个东西"。**所有 TRAIN 结果是头部拟合数据上的结果，所有 DEV 结果是单 seed。**

| 变体 | 假设 / 改动 | TRAIN（22,688 参考）pooled F1 | DEV（7,728 参考）pooled F1 | 判定 | 来源 |
|---|---|---|---|---|---|
| **原始首训** | 全语法观测似然：精确 fine 角色硬约束 + coarse 区域观测通道（ε=0.1，正确 0.9/其余各 0.05）+ 未知位置边缘化；4,608 步，seed0，lr 3e-4，batch 1 | 45.590（拟南芥 50.487/水稻 28.278） | 44.608（50.844/30.096） | 工程成立，**性能远低于三基线**（在线 loss 降 63.8%） | `full_grammar_first_fit_result.md`、`full_grammar_raw_train_result.md`、`full_grammar_dev_baseline_result.md` |
| 结构续训（STRUCTURED） | 再 4,608 步，仍是结构似然 | 55.745 | 54.751 | 对照 | manuscript Results L61–78；DEV L360–376 |
| LOCAL 续训 | 换成逐碱基 15 类局部似然（同解码器） | 41.729 | 42.181 | **更差**；说明 CDS-base F1 高（94.1）≠ 链 F1 | 同上 |
| MIXED 续训 | 0.5/0.5 | 48.628 | 48.943 | 水稻升、拟南芥降 | 同上 |
| native CDS prior | 加 28 个零初始化系数读冻结原生 CDS 头的 CDS−NONCODING 分数 | 64.339 | 63.776 | 比结构续训高；水稻不及 Helixer | `full_grammar_native_cds_prior_fit.md` |
| fresh uniform draws | 换成 4,608 个新抽的均匀 draw（seed1，4,219 个唯一窗） | 63.073 | 62.685 | 有增益，但 72.4% 的净增益来自"两个计划都没覆盖完整链"的参考 | `full_grammar_fresh_draws_plan.md`；manuscript L313–317 |
| complete-core draws | 重要性校正的混合抽样（4,107 个唯一窗） | 63.676 | 未进 DEV | 覆盖更多完整链，净命中反降 243（水稻 −414） | `full_grammar_complete_core_draws_plan.md` |
| **joint24（joint 24,576 编码）** | 四块 6,144 bp 改为一次整窗编码 | 66.963 | **66.623** | 相对结构续训 +11.2/+11.9，**最清晰的单因素收益**；仍低于 Tiberius 5.18、ANNEVO 9.08（DEV） | `full_grammar_joint24_fit.md`；manuscript L96、L399–401 |
| joint + native prior | 两者叠加 | 62.451 | 未进 DEV | 低于任一单项，非加性（交互 −13.1） | `full_grammar_joint24_native_prior_plan.md` |
| **fresh-joint** | fresh 抽样 + joint24 | 70.303 | **70.107** | 四格交互 −3.99；DEV 最佳 | `full_grammar_fresh_joint24_result.md`；`outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-FAMILY-R1/summary.json` |
| uniform-coarse | 额外学习阶段 ε 0.1→2/3（coarse 变成无信息） | 70.645（比 informative 高 0.34） | 69.795（比 informative 低 0.31） | **coarse 通道无一致收益**，两物种方向相反 | `full_grammar_uniform_coarse_result.md` |
| coarse 信息探针 | 固定头下 informative 对 uniform 的 logit 梯度 | 首个含 coarse 的 draw 差异为 0；最大 coarse 视图差异 1e-3~5e-3 | — | 诊断：仅在 coarse 行多的视图才有局部信号 | `full_grammar_coarse_information_probe.md` |
| **codon energy** | 64 个中心化密码子能量系数（含 ATG/终止子） | 65.469（比 fresh-joint −4.83，6 条染色体全降） | 64.579（−5.53） | **no-go，负结果** | `full_grammar_codon_result.md` |
| **joint-LoRA rank-4 q/v** | 最后 4 层 q/v 的 rank-4 FP32 残差，102,400 参数，原生 lr 3e-4 | **72.958**（+2.66，6/6 染色体都升） | **未跑 DEV** | TRAIN 最佳单模型；背景和非匹配上升（wholly-bg 407→454） | `full_grammar_joint_lora_result.md` |
| full-rank q/v | 同位置全秩 18.87M 参数，原生 lr 1e-5 | 72.200 | 未跑 DEV | 与 LoRA 方向一致 | `full_grammar_dense_qv_result.md` |
| **原生 lr 四格** | rank × native lr | rank4@1e-5 65.929；rank4@3e-4 72.958；full@1e-5 72.200；full@3e-4 **58.982** | — | **交互 −20.2**，排序随 lr 反转，不能说"哪种 rank 更好"，优化敏感 | `full_grammar_qv_lr_controls_result.md`；`outputs/M28-FULL-GRAMMAR-QV-LR-FOUR-CELL/summary.json` |
| E/M 强一致交集 | 保留 rank-4 与 full-rank 都预测的链 | 74.619（P 81.3 / R 68.9） | — | **事后 TRAIN 探索**；召回降、CDS-base F1 降到 90.4 | `full_grammar_e_m_consensus_train.md` |
| E/M 等权能量池 | 两模型 15 角色 logit 0.5/0.5 平均后再解码 | 73.857（R 72.8，CDS-base 95.5） | — | 事后 TRAIN；只在有限 view 做过 raw 资格检查，**未做整染色体新编码** | `full_grammar_e_m_energy_pool_train.md`；`full_grammar_e_m_raw_native.md` |
| 最小能量（局部分歧惩罚） | 取两分支较低中心化能量 | 72.773（比均值池 −1.08） | — | **no-go** | `full_grammar_minimum_centered_energy_train_plan.md` |
| 基因后验排序 | 精确 CDS 事件边缘概率 vs 原始 gain 对已输出预测排序 | 后验 AP 0.91 对 raw gain 0.75；AUROC 0.81 对 0.58 | — | 排序有效，但**不是校准概率**（≥0.9 档实际 primary 命中仅 85–86%）；未生成过滤版 GFF | `full_grammar_gene_posterior_train_plan.md`；`full_grammar_energy_pool_posterior_train.md` |
| 参考辅助上界 | 22,688 条中 22,479 条单链语法可支持，22,460 条可联合共存 | oracle F1 上限 99.495 | — | 结构容量不是缺口来源 | `full_grammar_train_energy_diagnostic.md` |
| 得分竞争回放 | 7,762 条漏检中 2,486 条能量为负、5,067 条为正但被竞争链压过、209 条语法不支持 | — | — | 诊断，描述性，非因果 | `full_grammar_matched_energy_diagnostic.md` |
| Carbon-A（外部对标） | 同一 GENERanno 主权重 + 转移 TCN + 30 状态 Viterbi + 置信 TCN | 过滤 **78.172**；不过滤 76.063；24,576 窗口版 77.562/75.585 | **未跑 DEV** | 外部完整流程，**全面高于我们** | `carbon_a_original_train_result.md`、`carbon_a_context24_train_result.md` |
| **q/v/o LoRA** | 在 rank-4 q/v 之外再加最后 4 层 o_proj 的 rank-4 残差，共 167,936 参数，head+native 总 534,127 | **未运行** | — | **停在接线阶段** | `full_grammar_joint_qvo_lora_plan.md`（2026-10-14） |

---

## 3. 任务 2：最终采纳的方法

先说明一点：**项目内没有"冻结的最终方法"**。草稿里的主线模型是 fresh-joint，TRAIN 最强单模型是 rank-4 q/v LoRA。以下是当前 manuscript 描述的"基础模型"，加上它的可选升级。

**模型结构**（`manuscript` Methods L510–L526）
- 冻结 GENERanno-eukaryote-1.2b-cds-annotator-preview（revision b0483c23…）。
- 24,576 bp 窗口：fresh-joint 一次整窗编码（4,096 个六碱基 token）。最初的版本是四个 6,144 bp 块。
- 输出 2,048 通道的对齐六碱基表示，BF16 编码后转 FP32。
- 头 C0（366,191 参数）：2048→128 线性投影 → 双向 GRU（每向 64）→ 与 32 通道 DNA 卷积分支（kernel 9 和 5，GELU）拼接 → 线性/GELU/LayerNorm → 每个碱基输出 15 个角色分数。
- 15 角色：背景、编码相位、内含子相位、donor/acceptor 边界、start、stop。
- 代码：`src/m28/core.py`（`SharedContext` 与 `M28Core('C0',2048,128)`）。**本地 `src/m28/` 里没有这个文件**，只在 `origin/main`（commit d77d225 一线）和 Baobab 上；本地的 58 个 m28 模块与 `scripts/experiments/M28-METHOD-RESTART` 都是未跟踪文件。
- 可选升级：最后 4 层 q/v 的 rank-4 LoRA（`src/m28/joint_encoder_lora.py`），原生 lr 3e-4。

**训练目标**（Methods L528–L552）
- 一个基因集合 S 的能量是单位权重的"角色分数 − 背景分数"逐碱基之和，空集能量为 0，没有学习到的 null 偏置、长度归一或阈值。
- 训练损失是 `log Zall − log Zobs`，其中 Zall 对完整语法下**所有**兼容基因集合（含空集）求和：
  - fine 标签作为硬约束（允许角色掩码）；
  - 只有 coarse 标签的位置按 0.9/0.05/0.05 的区域观测通道加权；
  - 完全未知位置不进入似然（边缘化，而不是当背景）。
- 按已知观测行数 `max(1, F+C)` 归一，乘窗口 importance 权重。
- 以有限视图（24,576 bp 核心加两侧各一个 12,288 bp stride 邻窗）为单位的复合似然，**不是整条染色体联合 MLE**。
- FP64 log 空间 sum-product，256 碱基检查点反向求导。
- AdamW lr 3e-4、weight decay 0.01、clip 1.0、batch 1，fresh seed0 起始 4,608 步，再续 4,608 步（总 9,216，不选中间 checkpoint）。
- 抽样：fresh-joint 用 seed1 的 4,608 个新均匀抽样（4,219 个唯一窗）。
- 代码：`src/m28/full_grammar_coarse_channel.py`（观测通道）；`full_grammar_boundary_messages.py`、`full_grammar_boundary_compiled.py`、`full_grammar_partition*.py`（配分函数）；训练入口 `scripts/experiments/M28-METHOD-RESTART/fit_full_grammar_coarse_channel.py`（首训）、`fit_full_grammar_fresh_joint24.py`（fresh-joint）、`fit_full_grammar_joint_lora.py`（LoRA）。

**语法与解码**（Methods L520–L526、L602–L610）
- 语法由 DNA 自动机定义：ATG 起始；TAA/TAG/TGA 终止且无内部同框终止；GT 或 GC 供体加 AG 受体；内含子至少 5 bp；每基因至多 128 个外显子；每条链上基因不重叠，反链可重叠；编码区不允许 N，内含子可含 N。
- 生产解码器是"最大分数 + exon 预算 skyline"（`full_grammar_skyline.py`，`full_grammar_viterbi.py` 为参考实现）：在同位置、同阶段、同实际 codon carry 下，分数不更差且 exon 用得更少的状态支配其他状态。FP64，8,192 bp 回溯检查点。

**推理流程**（Methods L602–L610）
1. 只输入 FASTA 和固定的最终 checkpoint（不读注释、标签、特征缓存或候选池）。
2. 24,576/12,288 的固定窗口网格，含末端窗口。
3. 每个碱基取"最近包含它的窗口中心"的 15 角色分数，不对重叠预测做平均。
4. 对每个有向 contig 做一次完整基因集合 MAP，基因可以跨窗口边界。
5. 导出 gene/mRNA/CDS 的 GFF3，并核对坐标和相位。
6. 代码：`contig_emission_lattice.py`（`stitch_contig_emissions`）、`full_grammar_skyline.decode_full_grammar_skyline`、远程的 `c0_decode.genomic_chain`、`export.write_gff3`；入口 `infer_full_grammar_fresh_joint24_train.py`、`infer_full_grammar_joint_controls_dev.py`。

**什么不是"最终方法"**：
- 原生 CDS prior、complete-core 抽样、codon energy、最小能量池：结果为负或无增益，不进主线。
- E/M 能量池和强一致交集：TRAIN 上的事后探索，不是已实现的、通过整染色体新编码评估的调用器。
- 论文自己的定位（Methods "Relation to prior work"、`docs/M28_complete_chain_competition_prior_work_20261003.md`）：配分函数/CRF、冻结 DNA 表征加结构解码、局部加结构联合训练、密码子能量**都有先例**（Conrad、GeneDecoder、Tiberius、Carbon-A、GENATATOR、AUGUSTUS）。草稿并不声称新颖，而是要求"竞争性能或独立证明的更广适用范围"，两者目前都没有达到。

---

## 4. 任务 3：当前最佳结果与基线

### 4.1 同尺对比（DEV：7,728 条 primary，单 seed，反复使用的开发染色体）

| 方法 | 拟南芥 F1 | 水稻 F1 | pooled P / R / F1 | 来源 |
|---|---:|---:|---|---|
| **M28 fresh-joint（DEV 最佳，非后选）** | 72.858 | 63.725 | 70.966 / 69.268 / **70.107**（TP 5,353/预测 7,543） | `outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-FAMILY-R1/summary.json` |
| M28 uniform-coarse | 73.318 | 61.641 | 70.30 / 69.29 / 69.795 | 同上 |
| M28 codon energy | 69.038 | 54.330 | 63.94 / 65.23 / 64.579 | 同上 |
| M28 joint24（单独） | 70.499 | 57.502 | 67.67 / 65.61 / 66.623 | manuscript L399 |
| 阶段一最佳 A31_A23 | 69.421 | 64.142 | 69.94 / 65.81 / 67.813 | `a31_local_path_caller_dev_result.md` |
| **ANNEVO 2.3.3** | 76.265 | 74.432 | 78.57 / 73.03 / **75.703**（TP 5,644/7,183） | 同 FAMILY-R1 summary |
| **Tiberius 2.0.8** | 74.177 | 66.582 | 73.26 / 70.41 / **71.805**（TP 5,441/7,427） | 同上 |
| **Helixer 0.3.7** | 67.456 | 53.524 | 59.86 / 66.15 / **62.851**（TP 5,112/8,539） | 同上 |

- **fresh-joint 对 Helixer**：两物种 F1 都高，P/R 也都高（pooled F1 +7.26，多 241 个命中、少 996 条预测）。
- **fresh-joint 对 Tiberius**：pooled 低 1.70，对 ANNEVO 低 5.60，两物种都低。
- **Helixer 的暴露口径**：当前权重明确训练过拟南芥和水稻，因此"超过 Helixer"比"超过未见过该物种的工具"更容易，也更不能外推。
- CDS-base F1：fresh-joint 95.065，高于 Tiberius 和 ANNEVO，但完整链 F1 反而低，这是草稿反复强调的"碱基覆盖≠链恢复"。
- 背景：fresh-joint 的 wholly-background 链 145（Tiberius 165、ANNEVO 146、Helixer 563），背景 span 覆盖 0.921%（Tiberius 0.721%）。

### 4.2 TRAIN 同尺（22,688 条参考；头部拟合数据，不是泛化）

| 方法 | 拟南芥 | 水稻 | pooled F1 | 来源 |
|---|---:|---:|---:|---|
| M28 fresh-joint（冻结编码器） | 73.368 | 59.636 | 70.303 | manuscript L106–L113 |
| M28 rank-4 q/v LoRA | 75.316 | 64.908 | **72.958** | `full_grammar_joint_lora_result.md` |
| M28 full-rank q/v | 74.415 | 64.550 | 72.200 | `full_grammar_dense_qv_result.md` |
| M28 E/M 等权能量池（事后，缓存重解码） | 75.439 | 68.413 | 73.857 | `full_grammar_e_m_energy_pool_train.md` |
| M28 E/M 强一致交集（事后） | 76.241 | 68.540 | 74.619 | `full_grammar_e_m_consensus_train.md` |
| ANNEVO | 76.831 | 71.179 | 75.542 | `native_train_baseline_result.md` |
| Tiberius | 75.417 | 63.577 | 72.655 | 同上 |
| Helixer | 68.122 | 50.294 | 63.592 | 同上 |
| **Carbon-A 过滤（置信 ≥0.1 或 NA）** | 78.958 | 75.259 | **78.172** | `carbon_a_original_train_result.md` |
| Carbon-A 不过滤 | 76.745 | 73.530 | 76.063 | 同上 |

**读法**：
- 即使在**它自己的拟合数据**上，M28 单模型也只是在 pooled 上勉强越过 Tiberius（+0.30，拟南芥反而低 0.10），低于 ANNEVO 2.58，低于 Carbon-A 过滤版 5.2。
- 水稻差距比拟南芥大：rank-4 LoRA 对 ANNEVO，水稻低 6.3、拟南芥低 1.5。
- 把"超过 Tiberius"写在 TRAIN 上不应作为竞争性证据，因为 M28 的头部在这些染色体上拟合过。草稿本身已承认这一点，但摘要和讨论仍用"narrowed the gap"。
- 对比：fresh-joint 的 TRAIN（70.303）和 DEV（70.107）几乎相同，说明头部（366k 参数、1,491 窗）没有明显过拟合，但也说明容量/训练量是瓶颈而不是泛化。

### 4.3 回答你的具体问题

| 问题 | 答案 | 依据 |
|---|---|---|
| 是否在 Setaria 盲测上跑过 | **没有**。注释从未被读取；也没有独立物种面板被冻结或执行 | `docs/28_current_research_state_2026-09-01.md` L74、L83；`independent_transfer_candidates.md`；各 M28 报告"sealed/Setaria 未访问" |
| 多 seed / 置信区间 | **没有**。所有结果均是单 seed、单次拟合；`summary.json` 的 uncertainty 字段写 "single_fit_historically_reused_two_species_DEV_not_clean_zero_shot_or_model_selection"；没有 bootstrap、没有基因级重采样 | `outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-FAMILY-R1/summary.json`；manuscript 全文反复声明 |
| 是否有消融 | **有，但只是 TRAIN 上的单 seed 因子对照**：结构 vs 局部 vs 混合目标（confounded，同解码器）；sampling × encoding 2×2；context × native prior 2×2；coarse informative vs uniform；codon on/off；LoRA vs full-rank × 两种原生 lr；complete-core 对 fresh。DEV 上只有 6 个固定候选（fresh-joint、uniform、codon）和 10 方法家族。**没有**：不同解码器的消融（A42 在阶段一做过，但对象不同）、"未知填背景 vs 边缘化"、有无 15 角色头对 3 类头、去掉 GRU/DNA 卷积分支、不同 ε、训练量扩大（只用了 1,491/25,512 个窗口） | manuscript Results；`full_grammar_objective_contrast_contract.md` |
| 最佳结果在哪个口径 | DEV：fresh-joint 70.107；TRAIN：LoRA 72.958（单模型）。阶段一最佳 A31_A23 67.813（DEV） | §4.1、§4.2 |

---

## 5. 任务 4：论文草稿核心主张的证据状态

状态标记：**【确认】**已有确认性证据；**【探索】**仅开发集或 TRAIN 的探索性证据；**【无证据】**没有证据；**【矛盾】**与结果矛盾。
注意：本项目的 DEV 被多个实验反复使用，所以"确认"只用于计算性事实或工程验证，不用于性能主张。

### 5.1 主张表

| # | 草稿主张（位置） | 状态 | 证据与问题 |
|---|---|---|---|
| 1 | 部分观测下，对所有兼容基因集合（含空集）求和的似然训练是明确定义的，且未知位置被边缘化而非当背景（abstract、Methods） | **【确认】**（实现层面） | 对人工序列用穷举验证了 logZ 和梯度（`full_grammar_partition_preparation.md`、`full_grammar_gene_posterior_preflight.md`）；这只是实现正确，不证明有益 |
| 2 | 推理只用 FASTA，不枚举参考信息候选，整 contig 精确 MAP（abstract） | **【确认】**（工程） | caller DNA 与原 DNA 逐碱基相等、GFF 回读核验、skyline 对 dense 一致；但"精确"仅指在该语法域内 |
| 3 | joint 24,576 编码使 DEV 链 F1 从 54.751（结构续训）升到 66.623（abstract） | **【探索】** | 单 seed、复用 DEV、仅一对；方向在 TRAIN（55.745→66.963）和 DEV 一致，两物种皆升。这是全项目证据最强的一项，但仍无重复 |
| 4 | 单独 fixed family 达到 70.107/69.795/64.579，三者均超过 Helixer 的 62.851 且低于 Tiberius 71.805/ANNEVO 75.703（abstract） | **【探索】** | 数字已与 `summary.json` 核对一致；超过 Helixer 的两物种 P/R/F1 均更高（fresh-joint）。Helixer 训练过这两个物种；单 seed；DEV 复用。**"低于 Tiberius/ANNEVO"是确认的负结果** |
| 5 | 结构化似然训练优于逐碱基局部似然（Results："structured > local"） | **【探索】** | TRAIN 55.745 对 41.729，DEV 54.751 对 42.181，同一解码器。但目标、系数、有效梯度均不匹配，且 LOCAL 的 CDS-base 更高；草稿自己写"not a decoder ablation"。不能写成"complete-set competition 是原因" |
| 6 | 部分观测通道（coarse 观测）带来收益 / 标题所述"partially resolved supervision" | **【无证据】且部分**【矛盾】 | uniform-coarse 对照：TRAIN 信息版比均匀版 pooled −0.34，拟南芥 −0.91、水稻 +1.54；DEV pooled +0.31，拟南芥 −0.46、水稻 +2.08，符号随面板翻转。**"未知边缘化 vs 填背景"从未被检验**；ε=0.1 是固定值，没有扫描或估计。草稿 Discussion 对此措辞诚实，但标题和 abstract 第 2–3 句仍把它当作方法的核心 |
| 7 | 密码子能量（coding content）改善完整链 | **【矛盾】** | TRAIN −4.83、DEV −5.53，各染色体全降（草稿承认是负结果） |
| 8 | 编码器适配（rank-4 q/v）提高 TRAIN 链 F1，使其达到 Tiberius 水平 | **【探索】**（TRAIN-only） | 72.958 对 70.303，6/6 染色体提升；但 (a) 未跑 DEV；(b) 原生 lr 四格显示对 lr 极敏感（full-rank@3e-4 降到 58.982，交互 −20.2）；(c) 参数化、lr、clip、优化器年龄全部与"秩"混杂 |
| 9 | 7,762 条漏检中 5,067 条"可单链表示、能量为正但被竞争链压过"，区分"得分竞争"与"结构排除"（abstract） | **【探索】**（描述性） | 回放计算本身可复现（TRAIN，单模型）。但草稿已写明不能因果解释；"competition 是主要缺口"没有改善实验支持——阶段一 A13–A15、A36–A43 一系列竞争/集合目标全部未过门槛（A34/A35 只是诊断） |
| 10 | 参考辅助 oracle 上限 99.495%（语法容量不是瓶颈） | **【确认】**（计算事实） | 确定性最大兼容集合计算；不是模型准确率 |
| 11 | 结构后验比原始 gain 更适合排序已输出的预测（AP 0.91 对 0.75） | **【探索】**（TRAIN） | 仅在 TRAIN 且为事后；pooled AUROC 低于两个父模型；不是校准概率 |
| 12 | 与 Tiberius/ANNEVO/Helixer 的竞争性能（intro 中"intended contribution"） | **【矛盾】（目标未达成）** | DEV 上一直低于 Tiberius 和 ANNEVO；TRAIN 上 LoRA 才在 pooled 上勉强越过 Tiberius；Carbon-A 全面领先。草稿已如实写出，但"contribution under test"仍悬空 |
| 13 | 更广适用范围/跨物种泛化 | **【无证据】** | 只有两个物种；Setaria 盲测未做；独立候选未冻结 |
| 14 | 计算效率优势 | 草稿**已明确不主张** | CPU/节点/计时口径不一致；仅有 allocation 记录 |
| 15 | E/M 强一致交集或等权能量池"追平 Tiberius/逼近 ANNEVO"（TRAIN） | **【探索】**（事后） | 在父模型已知后才设计规则；两个父模型本身是同一 TRAIN 的拟合；没有整染色体新编码评估 |

### 5.2 草稿中的过度表述与结构问题

1. **标题与摘要的"partially resolved supervision"为主线，但没有任何实验把它与"全监督"或"未知填背景"对照**（主张 6）。摘要前半段把似然设计写成方法贡献，后半段的数据却说明收益来自编码上下文、抽样顺序和编码器适配。建议降格为"设计选择"，或补做"未知填背景"对照。
2. **"narrows the gap / exceeded Tiberius"建立在 TRAIN 上**。TRAIN 是头部拟合数据，ANNEVO 的训练排除了这两个物种，Tiberius 不含拟南芥，Helixer 含两者。这种比较对基线不利。草稿虽有限定，但 Discussion 的"further narrowed the original-TRAIN gap"和 Abstract 中 LoRA/E-M 的叙述容易被读成竞争性证据。LoRA 的 DEV 从未跑过。
3. **"demonstrates method improvement within the evaluated setting"**（abstract 末句）：被改进的对象是本项目自己的早期模型（结构续训/初始模型），不是已有方法。对外部强基线，结论是"仍低于 2/3"。
4. **复用 DEV 的"predeclared / frozen family"措辞**：这两条 DEV 染色体在 M20/M26/A10/A31 等数十个实验中被反复读数，"frozen"只对某个家族内部成立；Methods 已承认 fresh-joint 设计是在 DEV 结果出来之后才定的（L564）。应在摘要里同时写"DEV 为开发面板"。
5. **时效不一致**：abstract 文件（Oct 12 22:30，标题行写 Oct 13 快照）没有包含 10-13/14 的 lr 四格、E/M 池、最小能量池和后验排序；manuscript 顶部段落已更新。`docs/12_publication_strategy.md`（有未提交修改）最后一次实质更新停在 10 月 9 日，缺少此后全部 fresh-joint/LoRA/DEV 六方法的结果。
6. **可读性/可投稿性**：manuscript 约 170 KB，混有 Slurm job ID、监控、恢复流程等执行记录语言，不是论文体例；Results 里"后实验"类、纯诊断类、阴性类与主结果混排，缺少一个主表加一个明确的"主结论"。
7. **代码可复现性缺口**：核心模型文件 `src/m28/core.py`（以及 `training.py`、`c0_decode.py`、`export.py`）不在本地工作树，manuscript 和图 caption 都指向它；本地 `main` 与 `origin/main` 已分叉，全部 M28 文件均为未跟踪。若不迁移，外部读者无法复现。
8. **阶段一的失败（A13–A44）几乎完全没有进入论文**。草稿只有"Context chain methods"的旁支文档；这些负结果（尤其 A42 的"精确语法解码使 A3 发射崩溃"）是转向 full-grammar 的关键证据，建议放入补充材料或讨论。

---

## 6. 任务 5：悬而未决的问题与停止位置

### 6.1 最后几个实验

| 时间 | 事件 | 状态 |
|---|---|---|
| 2026-10-14 07:11 | manuscript 最后一次编辑 | 已含四格 lr 对照 |
| 2026-10-14 07:02 | `full_grammar_minimum_centered_energy_train_plan.md`、`full_grammar_minimum_centered_raw.md`、**`full_grammar_joint_qvo_lora_plan.md`** | 前两者完成（最小能量池 TRAIN 负结果；raw 资格检查过）；**q/v/o LoRA 只有计划和 7 个 sbatch**，plan 中写明 04:59 UTC 未满足提交条件（运行 102、等待 1038，要求并发 <80）；四组 CPU 接线检查 → 一个 20 min GPU native 预检 → 20 h 正式 fit → 6 个 raw → 评价 → 配对，均**尚未开始** |
| 2026-10-14 06:29–06:36 | `outputs/M28-FULL-GRAMMAR-QV-LR-FOUR-CELL`、`full_grammar_qv_lr_controls_result.md` | 完成 |
| 2026-10-13 | E/M 一致、能量池、后验排序、缺失参考 R2（部分任务曾超时，已恢复） | 完成 |
| 2026-10-12 | 三个固定候选 DEV 六方法家族 R1（首次启动因 source-key 命名失败，schema-only R1 完成） | 完成 |

**关于你提到的 `sbatch/M28-A15-*RECOVERY*`**：这些**已经完成**，不是未完成项。A15 自然链集合的 R1 恢复拟合/读出/核验在 10-05 全部 COMPLETED，结果见 `a15_natural_chain_set_result.md`（3 PASS/2 FAIL，关闭）。
真正"没有终态"的是：

- `M28-FULL-GRAMMAR-JOINT-QVO-LORA-*`（7 个 sbatch，10-14 创建，未提交）。
- A33 的 emission-beam（指针重放检查 FAILED，五条件从未评价）。
- A44 的 paired raw / 机制诊断：有 sbatch 和协议，本地没有结果；后来的文档把 A44 归入"已停止"。
- A36/A37/A38/A39 的 DEV 推理/评价 sbatch：均已"准备"，但按事前规则因 TRAIN FAIL 不执行（不是遗漏）；A40/A41 没有 DEV sbatch。

### 6.2 悬而未决的科学问题

1. **是否有任何模型在"独立数据"上超过 Tiberius/ANNEVO**：目前没有。DEV 最佳比 Tiberius 低 1.7，比 ANNEVO 低 5.6；LoRA/E-M 的 DEV 缺失。
2. **LoRA 是否能迁移到 DEV 和 Setaria**：LoRA 在 TRAIN 上最强，但 lr 敏感、未进 DEV、未做多 seed。是否应该为 LoRA 追加 DEV，与草稿"不按 TRAIN 排名选候选进 DEV"的原则存在张力，需要用户决定。
3. **"partially resolved supervision"到底有没有贡献**：需要至少一个对照（未知填背景；有无 coarse；不同 ε）才能保住标题，否则要改定位。
4. **Carbon-A 的存在**：Carbon-A（同一 GENERanno 主权重、含转移网络与置信网络、2026-10-08 发布）在 TRAIN 上全面领先，并使"冻结 GENERanno 加语法解码"的新颖性大幅受限（草稿已承认其为先例）。需要决定论文的差异化点是什么，例如更广范围、机制、效率，而不是"更高 F1"。
5. **训练量**：全部 M28 头部只用 1,491 个窗口（约占 25,512 个自然窗的 6%），3 遍，约 4,608–9,216 步。A 阶段的 FULL（76,536 步，全窗）已存在，但 full-grammar 目标还没有放大到全窗。容量/数据是否为主要瓶颈未检验。
6. **解码和训练之间的不一致**：整染色体联合 MLE 未实现，目前是有限视图复合似然。
7. **未知链保护（-1）在 full-grammar 中如何处理**：通过边缘化，但 coarse-only 行大量重复曝光、"不是独立基因"的统计问题尚未处理。
8. **假阳性/背景负担**：LoRA 等提升 F1 的同时增加了 wholly-background 链和背景 span，是否可以用后验过滤换取精度，尚未生成过滤版 GFF 评测。
9. **多 seed、置信区间、重采样检验**：完全缺失；基因级 bootstrap 可以用现有预测事后补做，不需要重新训练。
10. **近期能做、成本低的补充**：
    - 用已有 TRAIN/DEV 预测做基因级 bootstrap，给表格加 CI。
    - 在 DEV 上给 LoRA（不选择，仅作为已预声明的第 4 个候选）补一次评测。
    - 做"未知填背景"对照。
    - 整染色体新编码评估 E/M 能量池。

---

## 附录 A：已读核心文件索引

- 论文草稿：`docs/M28_full_grammar_manuscript_working.md`（含 Abstract/Intro/Results/Discussion/Methods）、`docs/M28_full_grammar_{abstract,introduction,results,discussion,methods}_working.md`、`docs/M28_context_chain_methods_working.md`、`docs/M28_complete_chain_competition_prior_work_20261003.md`、`docs/figures/M28_full_grammar_model*.{svg,png,md}`
- 战略：`docs/12_publication_strategy.md`（含未提交修改）、`docs/28_current_research_state_2026-09-01.md`
- M28 报告：`reports/M28-METHOD-RESTART/`（约 180 个 .md，已读 a5–a44 的 result/protocol 要点、full_grammar_*、carbon_a_*、native_*、pretrained_exposure、independent_transfer_candidates 及上面各行所引文件；未逐字通读全部 protocol）
- 输出（本地仅有汇总）：`outputs/M28-FULL-GRAMMAR-*`（共 67 个目录），核对了 `JOINT-CONTROLS-DEV-FAMILY-R1`、`JOINT-LORA-PAIRED`、`QV-LR-FOUR-CELL` 的 `summary.json` 数值与 manuscript 一致
- 前史：`reports/M25R-*`、`reports/M26-SAME-SCOPE-MECHANISM/`、`reports/M20-SOTA-ERROR-ANALYSIS/report.md`
- 代码：`src/m28/*.py`（58 个模块 docstring）、`scripts/experiments/M28-METHOD-RESTART/`（360 个文件）、`origin/main:src/m28/core.py`

## 附录 B：本次梳理的局限

- A 系列只读了各自的 result/protocol 开头与关键表，没有逐条核对 Slurm 日志和 Baobab 上的原始输出（本地只有 67 个 full-grammar 汇总目录，没有 A10–A44 的 outputs）。
- A44 的最终数值未在本地找到，表中只记录"后续文档称失败/已停止"。
- 对 A42→full-grammar 动机的解释是我的推断，来自 A42 结果与 `full_grammar_partition_preparation.md` 的前后文，不是文档明示。
- 没有重新运行任何评价，也没有独立核对 TRAIN 数值；DEV 六方法数值已与 `summary.json` 核对。
