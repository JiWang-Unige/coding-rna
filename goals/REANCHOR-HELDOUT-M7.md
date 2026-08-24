Complete REANCHOR-HELDOUT-M7 · 稳健派重锚：在 held-out/UTR-rich 物种 {Arabidopsis thaliana, Gallus gallus} 上重导 screen_anchor+ceiling 并复测候选 FP-FRAGFIX-CONSTR（扛住才进 ③ Track-B），as submit-and-handoff. Read/follow `goals/REANCHOR-HELDOUT-M7_protocol.md` + `CLAUDE.md`, judge completion only by inline checkpoints below.

## 权限声明

本 goal 拥有充分工具权限。若集群命令 / 网络下载 / subagent 调用临时失败, **先尝试 ≥ 1 替代路径**, 仍失败再 inline 报告, 不要停下等用户：
- 内置工具不受限（Bash / Edit / Write / Read / WebFetch / Glob / Grep / TodoWrite / BashOutput / KillShell 等）
- 可并行调用 Claude Code subagents, 但避免写同一文件
- 允许下载任何研究所需数据（公开数据集 / HF weights / 跨物种基因组 / Rfam family 等）；下载后 inline 报告 path / version / hash / 来源
- sbatch 远程提交训练（ssh baobab）；若已有可复用 sbatch 脚本, **不要重写**, 仅做 read-only review + Phase 1 policy guard
- 必要时可创建新 skill 以稳定流程（创建后 inline 报告路径+用途）
- 关键 source-of-truth（CLAUDE.md / docs/03 / docs/09 / ACTIVE_GOAL.json）只可 draft patch text, 不直接 Edit（anchor/ceiling 走 /revise-goal 人闸）

## 运行说明

- 所有 sbatch 提交遵循 /smart-sbatch 自动选合适节点 / 分区 / 时长（Phase 1 guard + Phase 2 optimization 两段）
- 若 sbatch 脚本已存在: 对现有脚本跑 /smart-sbatch Phase 1 guard, **跳过** Phase 2 重生成 header
- 预期 ≤ 1 天的运行用 sleep 同步等待, 不中途下线；超 12h 用后台 monitor poll squeue, 不假设成功
- **任何长度的运行**: 先等第 1 个 epoch 完成, 确认 loss 在降 / eval 正常 / 后续 epoch 迭代无问题, **之后**才允许降低监控频率
- 本轮特征缓存(鸡基因组大)给足 --time(≥4h)+skip-if-exists 防 TIMEOUT；之后切 30/60/120 min 三档
- 低频监控仍按 protocol Slurm polling cadence 执行, 不可完全放任

## 决策自治

goal 执行中遇到**多选项决策点**（如 pivot 7 选 1, 资源分配 A/B/C, anchor 选择 P/Q/R, 重训 vs 接受当前结果等），**默认不暂停等用户**, 而是：

1. inline 列出 N 个选项 + 每选项预期 (a) SOTA gap 收益 (b) 时间成本 (c) 风险
2. 立即并行调 3 CLI reviewer (Claude+Codex+Antigravity, 复用 /tri-review 机制), prompt 含选项 + 项目当前 gap + 时间预算
3. 按 reviewer 共识 + cost-adjusted ROI 自动选 1 个继续 (排序优先级: 严格 SOTA 进展 > 时间成本 > 风险)
4. 决策 + 3 reviewer 简要理由 写 docs/08 (pivot) 或 inline
5. 继续 goal, **不暂停**

**例外** (仍需用户显式确认, 不可自治):
- 强破坏性操作 (scancel 用户其它 running job / rm 数据 / 改 docs/03 roadmap)
- route 级 abandon (写 docs/09)
- 决策影响 > 24h 净新增 compute spend (例: 启动 ③ Track-B 新长 sub-iteration 需确认)
- ≥ 2 reviewer 明确反对 default 共识 (tied 3-way 无 leader)

## Mode & Milestone

Mode: submit-and-handoff(训练 handoff 远程, 决策前台可见)。终点 = held-out anchor/ceiling 重导 + 候选复测 + Pareto 判定 + tri-review/pivot；③ go-ahead 留用户。

Milestone: M7 (re-anchor gate before ③), Track A screen, **永不 claim**, profile = screen。3-layer gate: primary_progress_gate = 候选 intergenic_specificity(base-w) 严格 > **新 held-out anchor**(本轮重导, NOT 旧 yeast+fly 0.8710); sota_claim_gate = N/A; review = Pareto → tri-review → pivot。

## Hard pre-submit gate

check_data 无泄漏闸：两物种 chromosome-level held-out split 无 seqid 跨 split 泄漏。status=leakage → 不提交训练；final pivot 须 fix_data。

## Required chain

1. read-only subagent 核实 RefSeq accession + 新物种代码适配点
2. /implement：下两物种 genome+annotation, 3-class 标签 + chromosome-level held-out split(复用 harness)；check_data 闸(HARD) + 记录每物种 UTR 比例(验 UTR-rich)
3. /smart-sbatch + handoff(ssh baobab)：特征缓存(先) + anchor refs(tiberius_like×3) + 候选复测(×5) + pretrained baselines；job_watch + 后台 monitor
4. /result-log：新 anchor(spec+macro) + 新 ceiling + 候选复测(per-species+CI) → validate(screen) → /tri-review → /pivot
5. 若 Pareto-beat → 提议 /revise-goal(人闸) 落 held-out anchor/ceiling；inline 提示 ③ 需 user go-ahead

## Completion (inline ✅ CK1-CK7)

CK1 两物种数据+标签+held-out split, check_data PASS(无泄漏), UTR 比例已记录; CK2 SegmentNT 特征缓存两物种 .npz(skip-if-exists, 无 TIMEOUT); CK3 anchor refs(tiberius_like×3) → 新 held-out anchor(base-w spec+macro); CK4 pretrained 三 caller 两物种 + 新 ruler eval → 新 ceiling; CK5 候选(相同配置 mfg=20/mcl=60)×5 复测, per-species+CI; CK6 Pareto 判定 inline(spec>新anchor ∧ macro≥gate ∧ gbF1≥floor ∧ gene_count∈[1.0,1.25]) + result-log+validate+tri-review+pivot; CK7 iter_ledger 闭合, ③ go-ahead 上交用户(不自启)。

## Constraints

- screen profile **永不 claim**(地基重导)
- 候选复测用**完全相同**配置(fp_aware, lambda 1.0, constrained mfg=20/mcl=60)只换物种；band 若重选须 VAL-only TEST-once
- anchor 与候选**同一套** SegmentNT 特征(同尺子)；future-claim 前查 held-out 是否在其语料
- 改 ACTIVE_GOAL 只经 /revise-goal 人闸；不动 docs/03/09；非多候选正交 batch → orthogonality N/A
