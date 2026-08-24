Complete TB-UNFREEZE-BACKBONE-M9 · ④ staged-unfreeze PREFLIGHT：解冻 NT-v2-500m backbone 顶 N 层 + 3c head，bounded screen 实测 frozen 是否封顶 gbF1, as submit-and-handoff. Read/follow `goals/TB-UNFREEZE-BACKBONE-M9_protocol.md` + `CLAUDE.md`, judge by inline checkpoints below. ⚠️ FOCUSED ARCH BATCH on training_signal（frozen control + unfreeze N2/N4）

## 权限声明

本 goal 拥有充分工具权限。若集群命令 / 网络下载 / subagent 调用临时失败, **先尝试 ≥ 1 替代路径**, 仍失败再 inline 报告, 不要停下等用户：
- 内置工具不受限（Bash / Edit / Write / Read / WebFetch / Glob / Grep / TodoWrite / BashOutput / KillShell 等）
- 可并行调用 Claude Code subagents, 但避免写同一文件
- 允许下载任何研究所需数据（公开数据集 / HF weights / 跨物种基因组 等）；下载后 inline 报告 path / version / hash / 来源
- sbatch 远程提交训练；若已有可复用 sbatch 脚本, **不要重写**, 仅做 read-only review + Phase 1 policy guard
- 必要时可创建新 skill 以稳定流程（创建后 inline 报告路径+用途）
- 关键 source-of-truth（CLAUDE.md / docs/03 / docs/09）只可 draft patch text, 不直接 Edit

## 运行说明

- 所有 sbatch 提交遵循 /smart-sbatch 自动选合适节点 / 分区 / 时长（Phase 1 guard + Phase 2 optimization 两段）
- 若 sbatch 脚本已存在: 对现有脚本跑 /smart-sbatch Phase 1 guard, **跳过** Phase 2 重生成 header
- 预期 ≤ 1 天用 sleep 等待不下线；超 12h 用后台 monitor poll squeue, 不假设成功
- **任何长度的运行**: 先等第 1 个 epoch 完成, 确认 loss 降 / backbone 顶层 grad 流动 / head 不 collapse, **之后**才降频
- 本轮 bounded screen 预计 <12h, 第 1 epoch 出后切 30/60/120 min 三档; 低频仍按 protocol cadence, 不放任

## 决策自治

goal 执行中遇到**多选项决策点**（如 pivot 7 选 1, 资源分配, freeze 深度 N, 重训 vs 接受等），**默认不暂停等用户**, 而是：

1. inline 列出 N 个选项 + 每选项预期 (a) SOTA gap 收益 (b) 时间成本 (c) 风险
2. 立即并行调 3 CLI reviewer (Claude+Codex+Antigravity, 复用 /tri-review 机制), prompt 含选项 + 项目当前 gap + 时间预算
3. 按 reviewer 共识 + cost-adjusted ROI 自动选 1 个继续 (排序优先级: 严格 SOTA 进展 > 时间成本 > 风险)
4. 决策 + 3 reviewer 简要理由 写 docs/08 (pivot) 或 inline
5. 继续 goal, **不暂停**

**例外** (仍需用户显式确认, 不可自治):
- 强破坏性操作 (scancel 用户其它 running job / rm 数据 / 改 docs/03 roadmap)
- route 级 abandon (写 docs/09)
- 决策影响 > 24h 净新增 compute spend (全量 scale-up 需确认; bounded screen 自治内)
- ≥ 2 reviewer 明确反对 default 共识 (tied 3-way 无 leader)

## Mode & Milestone

Mode: submit-and-handoff（训练 handoff 远程, 决策前台可见）。终点 = 3-arm 对比 + result-log + tri-review/pivot。

Milestone: M9 Track B PREFLIGHT (promoted_from M8), **NON-CLAIM**, profile=screen(bounded)。primary_progress_gate = unfreeze gbF1(干净植物 base-w) > frozen control AND > 14-elem 3c 0.7392, spec 不塌(~0.93+), gcount 不欠预测; sota_claim_gate=N/A; review=tri-review→pivot。

## Hard pre-submit gate

兼容 smoke(NT-v2-500m trust_remote_code load+.esm forward) + token-base 对齐(6-mer vs per-base 验长度) + check_data 无泄漏 + sanity smoke(顶层 grad 非零, 不 collapse, 不 OOM)。任一未过 → 不提交全量, final pivot = fix_eval/fix_data。

## Required chain

1. /implement：新 torch-backbone-unfreeze trainer(.esm 丢 LM head+不碰泄漏 U-Net; freeze→解冻顶 N; token-base 对齐; bf16+grad ckpt; 复用 3c FP-aware constrained head)；先 plan 勘查 modeling_esm; 兼容 smoke+check_data+sanity smoke
2. /smart-sbatch(SOFT_WARN focused training_signal) + handoff：3 arm(frozen/N2/N4)×3 seed; job_watch+monitor
3. /result-log：3-arm vs 14-elem 3c 0.739, 干净植物 base-w spec/gbF1/gcount+CI → validate → tri-review → pivot
4. unfreeze 抬 gbF1 → full scale-up(USER 硬闸)；不抬 → 关键负结果 pivot(换 foundation/domain-adapt/重审)

## Completion (inline ✅ CK1-CK6)

CK1 兼容 smoke PASS(NT-v2-500m load+.esm forward)+token-base 对齐验证; CK2 unfreeze trainer+check_data PASS+sanity smoke(grad 流动, 不 collapse, 不 OOM); CK3 3 arm×3 seed 全量 COMPLETED; CK4 3-arm(frozen/N2/N4)+vs 3c 0.739, 干净植物 base-w spec/gbF1/gcount+CI; CK5 gate inline(unfreeze gbF1 > control 且 >0.739 ∧ spec 不塌 ∧ gcount sane)+result-log+validate+tri-review+pivot; CK6 iter_ledger 闭合, 下一步(full scale-up USER 硬闸/负结果重审)上交。

## Constraints

- **NON-CLAIM**；**LEAKAGE**：只在 segmentation-clean 植物 {arabidopsis,rice} 评估(chicken/fly 是 SegmentNT fine-tune 物种、永远污染、不进 M9); test 不进 early-stop/decode; raw-DNA ab-initio; 同 3-class ruler
- **反调参**：gbF1 gap 0.16≫0.05 → training_signal 架构轴; token-base 对齐是关键正确性点
- 显存 bf16+grad ckpt+冻结层 no-grad; 3090 24GB; fallback generanno 若 5.11 炸
- bounded screen 自治内; **全量 full-unfreeze scale-up = USER 硬闸**; 改 ACTIVE_GOAL 只经 /revise-goal; 不动 docs/03/09
