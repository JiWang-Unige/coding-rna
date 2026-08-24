# REANCHOR-HELDOUT-M7 · Protocol

> 稳健派重锚 (retrospective 2026-06-11, docs/08 `## Retrospective 2026-06-11`). Track A screen, NON-CLAIM.
> 核心目的：当前所有 spec 数字（含 promote-ready FP-FRAGFIX-CONSTR 0.9218）全在 yeast+fly（低-UTR、基因密集、in-corpus 离群物种）上；北极星 cross-species intergenic stability 从未跨支系评估。本 goal 在 held-out / UTR-rich 物种集上重导 anchor+ceiling 并复测候选，**扛住了才进 ③ Track-B**。

## Permissions
- 内置工具不受限（Bash / Edit / Write / Read / WebFetch / Glob / Grep / TodoWrite / BashOutput / KillShell）。
- 可并行 read-only subagent（Explore / general-purpose）做数据可获得性核实 / 代码勘查；写入型仅限本 exp_id scope。
- 允许下载公开 RefSeq genome+annotation（Arabidopsis thaliana, Gallus gallus）；下后 inline 报告 path/version/accession/md5。
- sbatch 经 `ssh baobab` 提交（remote_ssh，submit_guard 拦本机 sbatch）；已有可复用 sbatch 脚本只 read-only review + Phase 1 guard，不重写。
- 关键 source-of-truth（CLAUDE.md / docs/03 / docs/09 / ACTIVE_GOAL.json）只可 draft patch text，**不直接 Edit**。改 ACTIVE_GOAL（anchor/ceiling）必须走 `/revise-goal` 人闸——本 goal 只产数据 + 提议 diff。

## Final goal (milestone)
- Milestone: M7 (retrospective-derived re-anchor gate before ③ Track-B). Track A screen. Profile screen. **永不 claim**。
- 终极北极星不变：在 intergenic_specificity 上严格超越 frozen published SOTA、达可发表（M2 sota_benchmark 仍 pending）。本 goal 是其前置去风险，不直接推进 claim。
- 3-layer gate（沿用 ACTIVE_GOAL R6 dual co-primary）：
  - primary_progress_gate：候选在 **held-out** 物种集上 intergenic_specificity（base-weighted headline）严格 > **新 held-out screen_anchor**（本 goal 重导）。
  - sota_claim_gate：N/A（screen，永不 claim）。
  - review_decision_gate：Pareto 判定（见下）→ /tri-review → /pivot。
- **Pareto 判定（held-out 复测候选 FP-FRAGFIX-CONSTR mfg=20/mcl=60）**：
  1. intergenic_specificity（base-weighted）严格 > 新 held-out anchor（seed-mean），AND
  2. macro_intergenic_specificity ≥ 新 held-out macro gate（= held-out anchor macro，重导得到），AND
  3. gene_body_F1_unconstrained ≥ floor（screen floor 0.5276，沿用；若 held-out anchor F1 更高则取 max），AND
  4. predicted_gene_count_ratio ∈ [1.0, 1.25]（两侧 band，per-species 报告）。
- 参照值（旧 yeast+fly ruler，仅供对比，**不作 held-out 判定锚**）：anchor spec 0.8710(3s)/0.8436(5s)，candidate spec 0.9218，ceiling 0.9917。

## Track + resource
- Track A screen，NON-CLAIM。这**不是**多候选正交架构 screen batch——候选架构已固定（FP-FRAGFIX-CONSTR，M4-M6 已验证）；本轮是 **re-anchoring + single-candidate held-out 复测**。Orthogonality 表 = N/A（无竞争候选；anchor refs / pretrained baselines 是参照点重导，非被筛架构）。
- 统一 frozen 协议（与既有 same-budget screen 一致，仅物种变）：sample_fraction 0.3，chromosome-level held-out split（`src/screen_anchor/data.py` assign_splits i%5：==3 val, ==4 test），window 2048，epochs/patience 沿用既有 screen 配置，seeds：anchor refs 3、候选 5。
- Partition：优先 `private-teodoro-gpu`（RTX3090 24GB，免费 7d）；特征缓存（鸡基因组大、rematerialization-heavy）给足 `--time`（≥4h，参 fly 经验 + 鸡 ~7× 体量），idempotent skip-if-exists 防 TIMEOUT mid-savez。VRAM ≥20GB。
- pretrained baselines（Tiberius/Helixer/ANNEVO）用既有 env / SIF（Helixer SIF refs/containers；annevo 专用 conda env；Tiberius refs/repos）。ANNEVO 用 `annevo` env（sanctioned exception）。

## Execution mode details
- Mode = **submit-and-handoff**：新数据 + 特征缓存（可能 >12h，尤其鸡）+ 多训练 job → 训练 handoff 远程；**决策逻辑全程前台可见**（每步 inline）。
- submit-and-handoff marker：每个 sbatch 提交后写 job_id/output_dir/resume 到 docs/05 Run tracker；用后台 monitor（`run_in_background` poll squeue）等 job，不前台 sleep 阻塞、不中途下线假设成功。
- 第 1 个训练 job 出第 1 epoch 后确认 loss 降 / eval 正常，再降监控频率（30/60/120 min 三档）。
- 特征缓存先于训练（依赖）；anchor refs / pretrained baselines / 候选复测彼此独立，可并行提交（≤3 concurrent directions，受 cluster_config）。

## Pre-submit gate (HARD)
**check_data 无泄漏闸（不可跳）**：两物种数据下载 + 标签生成后，跑 `scripts/check_data.py`，确认 chromosome-level split **无 seqid 跨 split 泄漏**（held-out 的核心正确性）。`status=leakage` → **停**，不提交任何训练。`pass` 才继续。
- 另：held-out 严格性自检——确认 val/test 染色体与 train 不重叠；记录每物种 CDS/intron/intergenic/UTR 标签比例（验证 UTR-rich：Arabidopsis/Gallus 的 UTR 占比应显著 > yeast）。

## Pre-submit code/data review
非平凡新物种接入（鸡基因组大、植物 contig 多）→ 提交前用 1 个 read-only subagent（code-plan-reviewer 或 Explore）核对：(a) extract_segmentnt.py 对大基因组的 tile/内存路径；(b) data.py 标签生成对植物多 contig / 鸡 scaffold 的处理；(c) eval 新 full-transcript ruler 对含 UTR 注释的 intergenic 补集计算。CRITICAL 未修不提交。

## Subagent fan-out
- read-only：数据可获得性 + accession 核实（RefSeq Arabidopsis thaliana `GCF_000001735.x`, Gallus gallus `GCF_*` 当前 RefSeq release）；extract/data/eval 代码对新物种的适配勘查。
- 写入型（如需并行起 config）：每 subagent 只写自己 exp_id scope（REANCHOR-HELDOUT-M7-* 命名空间）。
- 主 agent merge + 最终 sbatch + 最终写 docs。subagent 不再 spawn subagent。

## Pivot decision menu (track-A-screen / re-anchor gate)
- **候选 Pareto-beat 新 held-out anchor** → promote-ready 结论**加固** → 下一步 ③ Track-B（USER GO-AHEAD，>24h compute 例外）；同时 `/revise-goal`（人闸）把 held-out anchor/ceiling 落进 ACTIVE_GOAL。
- **候选不扛**（spec 掉到 held-out anchor 下 / macro 崩 / gene_count 出 band）→ 这正是 retrospective 要抓的失效：**不进 Track-B**，pivot = change-objective-or-loss / replace-component，把方向拉回「FP-aware emissions 在 UTR-rich 跨支系上为何失效」的分析（per-species 诊断：哪个物种崩、UTR 区是否被误判、SegmentNT 特征是否跨支系退化）。
- **anchor refs 自身在 held-out 上 fragmented / 退化**（如 helixer_like 在 yeast+fly 曾 fragment）→ 重导锚需排除 fragmented ref，记录原因（沿用既有 coherent-ref 选锚纪律）。
- 不论结果，held-out anchor/ceiling 重导本身是 binding 地基修正，但**只经 /revise-goal 落盘**。

## Pivot autonomy 边界（本 goal 特定）
- anchor 重导后的「选哪个 ref 当 anchor」「候选是否 Pareto」判定 → 决策自治（列选项→tri-review→共识→写 docs/08）。
- **例外仍 pause**：启动 ③ Track-B（>24h compute / 新长 sub-iteration）；改 ACTIVE_GOAL（/revise-goal 人闸）；route abandon（docs/09）；scancel 他人 job。

## Skill invocation chain
| step | skill | 说明 |
|---|---|---|
| (done) | /retrospective | 已跑，本 goal 是其 advisory 的稳健派落地 |
| (done) | /goal-prompt | 本文件 |
| 1 | /implement | 新物种数据接入 + 特征缓存 + 复测脚本（复用 harness，仅物种适配）；self-review + check_data + smoke |
| 2 | /smart-sbatch | Phase 1 guard + Phase 2（分配感知）；特征缓存 / anchor refs / 候选复测 batch matrix |
| 3 | submit-and-handoff | 远程提交 + job_watch 对账，后台 monitor |
| 4 | /result-log | 新 held-out anchor/ceiling + 候选复测 + Pareto 判定 + per-species/CI |
| 5 | validate_goal.py | 确定性闸（screen，非 claim） |
| 6 | /tri-review | held-out 重锚有效性 + 候选是否真扛跨支系 |
| 7 | /pivot | 单一决策（promote-ready 加固→③ / 或失效分析） |
| 8 | /revise-goal (人闸) | 落 held-out anchor/ceiling 进 ACTIVE_GOAL |
| 9 | /exp-log + iter_ledger | 链路闭合 |

## Constraints (full)
- screen profile，**永不 claim**（held-out 也不行；这是 direction-selection + 地基重导）。
- check_data 泄漏闸不可跳；held-out split 必须 chromosome-level 且 val/test 不与 train 共享 seqid。
- 复用现有 harness（src/screen_anchor, src/foundation_probe, scripts/eval_gene_body_mask.py 新 full-transcript ruler）；不重造、不改 shared 训练码除非必要且 inline 声明。
- 候选复测必须用**完全相同**的候选配置（--loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20）——只换物种/特征缓存，否则不可比。constrained band 若按 Track-B 协议在新物种 VAL 上重选，必须 VAL-only、TEST apply once（防泄漏），并 inline 声明。
- anchor/候选用**同一套** SegmentNT 特征（同尺子公平比较）；SegmentNT 见过 ~850 物种，screen 不 claim 故可接受，但 future-claim 前须查 held-out 物种是否在 SegmentNT 语料（pre-claim leakage guard）。
- 不动 docs/03（roadmap 主体）/ docs/09（无 abandon）；ACTIVE_GOAL 改动只经 /revise-goal。
- 与既有 FP-001..006 差异化：同固定架构候选，但**首次**在 held-out / UTR-rich 跨支系物种上验证（之前全 yeast+fly）——非 abandoned cousin（docs/09 空）。
