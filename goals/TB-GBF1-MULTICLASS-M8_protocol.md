# TB-GBF1-MULTICLASS-M8 · Protocol

> ③ Track-B scale-up, promoted_from FP-FRAGFIX-CONSTR via REANCHOR-HELDOUT-M7 held-out re-anchor gate (docs/08 pivot 2026-06-12, tri-review 2/3 DEGRADED scale-to-track-b). USER GO-AHEAD given. submit-and-handoff (>24h). NON-CLAIM (M2 sota_benchmark pending).
> CORE: M7 confirmed the architecture (frozen SegmentNT feats + FP-aware loss + constrained post-proc) DOMINATES AXIS-1 intergenic_specificity cross-clade (held-out 0.9604 vs anchor 0.8054) but LOSES AXIS-2 gbF1 (0.666 < raw-DNA anchor 0.710; gap to ANNEVO ceiling 0.8976 = 0.231 = ARCHITECTURAL, not tunable). job#1 = recover gbF1 via richer multi-class structured output, NOT more spec.

## Permissions
- 内置工具不受限（Bash/Edit/Write/Read/WebFetch/Glob/Grep/TodoWrite/BashOutput/KillShell）。
- 并行 read-only subagent（Explore/code-plan-reviewer）勘查 + 训练前 code/metric/split review；写入型仅限本 exp_id scope。
- sbatch 经 `ssh baobab` 提交（remote_ssh，submit_guard 拦本机）；已有可复用脚本只 read-only review + Phase1 guard。
- 改 ACTIVE_GOAL（screen_anchor/ceiling/success_criteria/sota_benchmark）只经 `/revise-goal` 人闸；docs/03/09 只 draft patch，不直接 Edit。
- 首个全量训练是 >24h 承诺：USER GO-AHEAD 已给；但 smoke 通过前不提交全量。

## Final goal (milestone)
- Milestone M8, Track B scale-up. profile=full（更多 data/epochs/seeds），**NON-CLAIM**（M2 sota_benchmark 未冻结 + human gate 前不 claim）。
- 北极星：在 intergenic_specificity 严格超 published SOTA 且 **不牺牲 gene-F1**、跨支系泛化、可发表。M8 攻 AXIS-2 gbF1 短板（M7 已 de-risk AXIS-1 spec）。
- 3-layer gate：
  - **primary_progress_gate**：constrained_gene_body_F1 **不再 < raw-DNA anchor 0.710**（或明确朝 ANNEVO ceiling 0.8976 恢复的趋势）**AND** intergenic_specificity 维持 >= ~0.95（held-out base-w）。
  - sota_claim_gate：N/A（M2 pending；NON-CLAIM）。
  - review_decision_gate：3-strata 结果 → /tri-review → /pivot。
- 参照（held-out, base-w, M7）：same-budget anchor spec 0.8054 / gbF1 0.7099；候选 spec 0.9604 / gbF1 0.6664；ANNEVO ceiling spec 0.9824 / gbF1 0.8976。yeast+fly anchor spec 0.8436。

## Track + resource
- Track B scale-up，promoted_from FP-FRAGFIX-CONSTR（M7 held-out gate 通过：spec Pareto-admissible）。**非多候选正交 batch** → orthogonality N/A（单架构轴 = multi-class structured output）。
- 架构轴 major_axis=decoder/training_signal；mechanism_delta = strand-aware multi-class labels (CDS/intron/intergenic + phase 0/1/2 + splice donor/acceptor border) + semi-CRF/segment-level decoder with biologically-meaningful transitions + FP-aware objective, on FROZEN SegmentNT features.
- semi-CRF 需 vectorize segment DP（M3 deferred）；若 intractable → 退回 linear-chain CRF on multi-class labels + structural transition constraints（记录为 tractability fallback，不算失败）。
- scale：更大 sample_fraction 和/或 + 物种；epochs 增；**>=5-8 seeds + mean±CI + 两侧 gene_count band**。
- 资源：private-teodoro-gpu（若空）或 shared-gpu + `--constraint=COMPUTE_TYPE_AMPERE`（VRAM>=24GB）+ exclude 3080；expected>12h → submit-and-handoff。多分类 + semi-CRF 可能更吃显存/时间，/smart-sbatch 据此定 --time/--mem。
- 保留确定性 constrained post-proc（mfg=20/mcl=60，per-clade band 在 held-out VAL 重选）作 coherence 层。staged UNFREEZE/fine-tune SegmentNT = **单独后续轴**，不与本轮混（归因）。

## Execution mode details
- submit-and-handoff：训练 handoff 远程；决策逻辑全程前台可见（每步 inline）。job_id/output_dir/resume 写 docs/05 Run tracker；后台 monitor poll squeue，不前台阻塞、不假设成功。
- 第 1 epoch 出后确认 loss 降 / eval 正常 / 多分类 head 不 collapse，再降监控频率（30/60/120 三档）。

## Pre-submit gate (HARD)
1. **多分类标签正确性**：新 build_labels 从 GFF 提 CDS/intron/intergenic + phase(0/1/2, 从 CDS phase 列) + splice donor/acceptor border（intron 边界）；inline 校验每类 base 比例合理（phase 三类近均衡、splice border 稀疏但非零）+ 与旧 3-class 一致性（CDS∪intron∪intergenic 覆盖）。
2. **check_data 无泄漏闸**（不可跳）：chromosome-level held-out split 无 seqid 跨 split 泄漏（沿用 M7 split，新增物种重跑 check_data）。
3. **sanity smoke**（小规模真跑，srun/sbatch smoke profile）：多分类 + semi-CRF/CRF 解码端到端跑通、产出 metrics、head 不 collapse（多分类易 collapse 到 dominant 类——必查 per-class recall）。smoke 失败 → repair_advisor 有界修复 ≤3 次。
任一未过 → 不提交全量；final pivot = fix_eval / fix_data。

## Pre-submit code review
多分类 build_labels + semi-CRF/CRF multi-class decoder + 3-strata eval 是非平凡新代码 → 提交前 1 个 read-only subagent（code-plan-reviewer）核对：label 提取正确性（phase/splice 最易错）、decoder transition 合法性、metric 实现与 SOTA 可比、split 契约、无 target 泄漏。CRITICAL 未修不提交。建议非平凡架构先进 plan 模式勘查 refs/repos（Tiberius/Helixer 的 multi-class label 实现）再落码。

## Subagent fan-out
- read-only：勘查 refs/repos 的 multi-class/phase/splice label 实现（Tiberius 5-state、Helixer phase）作参考；勘查 src/foundation_probe/train_probe_head 改多分类的最小改动点；3-strata eval 实现勘查。
- 写入型：仅本 exp_id scope（TB-GBF1-MULTICLASS-M8-* 命名）。主 agent merge + 最终 sbatch + 最终写 docs。subagent 不再 spawn。

## MANDATORY eval upgrades (reviewer 共识，并行非阻塞)
(a) **3-strata 报告**：Arabidopsis / Gallus-microchromosome（≤20Mb，M7 已有）/ Gallus-MACROCHROMOSOME（>20Mb，用保留的 genome.full.fa 取大染色体子集，重做标签+split+特征缓存）。每 stratum 报 intergenic_specificity / macro_specificity / gene_body_F1 / gene_count_ratio。macrochromosome 是 spec 最吃紧、M7 未测的 regime。
(b) **SegmentNT 预训练物种重叠核查**：确定性查 SegmentNT(multi_species) 预训练语料是否含 arabidopsis/gallus（查 model card / 论文 supplement / HF repo 物种表）→ 写入 pre-claim leakage gate（docs/10 + ACTIVE_GOAL pre-claim guard）。held-out 在特征层可能被污染——关系 full/scale claim，非本 NON-CLAIM 轮致命，但必须留痕。

## Pivot decision menu (track-B scale-up)
- **gbF1 恢复 + spec 维持**（primary gate met）→ 强证据：foundation+multi-class structured 路线在两轴都站得住 → 下一步 M2 freeze sota_benchmark + /revise-goal status active + 走 claim 路径（human gate + pre-claim SegmentNT 审）。
- **gbF1 部分恢复但仍 < ceiling**（趋势向上）→ continue Track-B：scale 更多 data / 调 decoder 结构（仍架构轴）。
- **gbF1 在 spec 约束下结构封顶**（多分类也救不回，spec↔gbF1 张力不可破）→ **关键负结果**：pivot 回架构重审（换 decoder family / 重新权衡 co-primary / 怀疑 frozen-feature 表达力上限 → 考虑 staged unfreeze 轴），写 docs/08；**不白烧后续算力**。
- **macrochromosome specificity collapse**（spec 在 gene-sparse 区崩）→ 这削弱 M7 的 spec 结论 → 回 FP-aware emissions 在 gene-sparse 区的失效分析。

## Skill invocation chain
| step | skill | 说明 |
|---|---|---|
| (done) | /goal-prompt | 本文件 |
| 1 | /implement | 多分类 build_labels + multi-class semi-CRF/CRF decoder + 3-strata eval + SegmentNT-overlap audit；self-review + check_data + smoke |
| 2 | /smart-sbatch | Phase1 guard + Phase2（多分类/semi-CRF 吃显存→定 --mem/--time）；>=5-8 seed batch matrix |
| 3 | submit-and-handoff | ssh baobab 提交 + job_watch + 后台 monitor |
| 4 | /result-log | 3-strata + 5-8 seed mean±CI + gbF1 恢复判定 + spec 维持 |
| 5 | validate_goal.py → /tri-review → /pivot | gbF1 恢复? spec 维持? macrochromosome? 负结果? |
| 6 | /revise-goal (人闸) | 若进 claim 路径：M2 sota_benchmark + status active；或记 held-out anchor/ceiling |
| 7 | /exp-log + iter_ledger | 链路闭合 |

## Constraints (full)
- **NON-CLAIM**：profile=full 但 M2 sota_benchmark 未冻结 + human gate 前永不 claim。
- **反调参硬闸**：gbF1 gap 0.231 ≫ 0.05 → tuning_allowed=false on gbF1 轴 → 本轮必须结构性（multi-class 输出 + decoder）改动，禁 lr/batch/dropout 调参。
- 保留确定性 constrained post-proc 作 coherence 层；band per-clade 在 held-out VAL 重选（VAL-only TEST-once 防泄漏）。
- staged UNFREEZE SegmentNT = 单独后续轴，不与本轮混。
- macrochromosome stratum + SegmentNT-overlap audit 是 MANDATORY eval 升级（不阻塞训练，但 result-log 必含）。
- 多分类 head 易 collapse → smoke 必查 per-class recall；不可只看聚合 F1。
- 复用 M7 {arabidopsis,gallus} held-out + 原 {yeast,fly}；新增 macrochromosome 数据走标准 build_labels+check_data。
- 不动 docs/03/09；ACTIVE_GOAL 只经 /revise-goal 人闸。非 abandoned cousin（docs/09 空）；与 M4-M7 差异：首次 strand-aware 多分类结构化输出主攻 AXIS-2 gbF1。
