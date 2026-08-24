Complete TB-GBF1-MULTICLASS-M8 · ③ Track-B（promoted_from FP-FRAGFIX-CONSTR）：strand-aware 多分类输出恢复 gbF1（不补 spec），as submit-and-handoff. Read/follow `goals/TB-GBF1-MULTICLASS-M8_protocol.md` + `CLAUDE.md`, judge completion only by inline checkpoints below.

## 权限声明

本 goal 拥有充分工具权限。若集群命令 / 网络下载 / subagent 调用临时失败, **先尝试 ≥ 1 替代路径**, 仍失败再 inline 报告, 不要停下等用户：
- 内置工具不受限（Bash / Edit / Write / Read / WebFetch / Glob / Grep / TodoWrite / BashOutput / KillShell 等）
- 可并行调用 Claude Code subagents, 但避免写同一文件
- 允许下载任何研究所需数据（公开数据集 / HF weights / 跨物种基因组等）；下载后 inline 报告 path / version / hash / 来源
- sbatch 远程提交训练；若已有可复用 sbatch 脚本, **不要重写**, 仅做 read-only review + Phase 1 policy guard
- 必要时可创建新 skill 以稳定流程（创建后 inline 报告路径+用途）
- 关键 source-of-truth（CLAUDE.md / docs/03 / docs/09 / ACTIVE_GOAL.json）只可 draft patch text, 不直接 Edit

## 运行说明

- 所有 sbatch 提交遵循 /smart-sbatch 自动选合适节点 / 分区 / 时长（Phase 1 guard + Phase 2 optimization 两段）
- 若 sbatch 脚本已存在: 对现有脚本跑 /smart-sbatch Phase 1 guard, **跳过** Phase 2 重生成 header
- 预期 ≤ 1 天的运行用 sleep 同步等待, 不中途下线；超 12h 用后台 monitor poll squeue, 不假设成功
- **任何长度的运行**: 先等第 1 个 epoch 完成, 确认 loss 在降 / eval 正常 / 多分类 head 不 collapse, **之后**才允许降低监控频率
- 本轮预计 >24h, 第 1 epoch 出后切 30/60/120 min 三档
- 低频监控仍按 protocol Slurm polling cadence 执行, 不可完全放任

## 决策自治

goal 执行中遇到**多选项决策点**（如 pivot 7 选 1, 资源分配 A/B/C, decoder family 选择, 重训 vs 接受当前结果等），**默认不暂停等用户**, 而是：

1. inline 列出 N 个选项 + 每选项预期 (a) SOTA gap 收益 (b) 时间成本 (c) 风险
2. 立即并行调 3 CLI reviewer (Claude+Codex+Antigravity, 复用 /tri-review 机制), prompt 含选项 + 项目当前 gap + 时间预算
3. 按 reviewer 共识 + cost-adjusted ROI 自动选 1 个继续 (排序优先级: 严格 SOTA 进展 > 时间成本 > 风险)
4. 决策 + 3 reviewer 简要理由 写 docs/08 (pivot) 或 inline
5. 继续 goal, **不暂停**

**例外** (仍需用户显式确认, 不可自治):
- 强破坏性操作 (scancel 用户其它 running job / rm 数据 / 改 docs/03 roadmap)
- route 级 abandon (写 docs/09)
- 决策影响 > 24h 净新增 compute spend (本轮首训已 go-ahead; 但再启一条新长 sub-iteration 需确认)
- ≥ 2 reviewer 明确反对 default 共识 (tied 3-way 无 leader)

## Mode & Milestone

Mode: submit-and-handoff（训练 handoff 远程, 决策前台可见）。终点 = gbF1 恢复判定 + 3-strata result-log + tri-review/pivot。

Milestone: M8 Track B scale-up（promoted_from FP-FRAGFIX-CONSTR via M7 gate）, **NON-CLAIM**(M2 pending), profile=full。3-layer gate: primary_progress_gate = constrained_gene_body_F1 **不再 < raw-DNA anchor 0.710**(或明确朝 ANNEVO ceiling 0.8976 恢复趋势) AND intergenic_specificity 维持 >=~0.95; sota_claim_gate = N/A; review = 3-strata → tri-review → pivot。

## Hard pre-submit gate

多分类标签正确性(phase/splice border) + check_data 无泄漏 + sanity smoke(多分类+semi-CRF 跑通、head 不 collapse、per-class recall 非零)。任一未过 → 不提交全量, final pivot = fix_eval/fix_data。

## Required chain

1. /implement：多分类 build_labels(CDS/intron/intergenic+phase 0/1/2+splice donor/acceptor) + multi-class semi-CRF/CRF decoder(intractable 退 linear-chain CRF+structural transitions) + 3-strata eval(arab/gallus-micro/MACRO) + SegmentNT 语料核查；self-review+check_data+smoke(per-class recall)
2. /smart-sbatch + handoff(ssh baobab)：多分类/semi-CRF 吃显存→定 --mem/--time；>=5-8 seed；job_watch + 后台 monitor
3. /result-log(3-strata + CI + gbF1 恢复 + spec 维持 + macro 无 collapse + SegmentNT overlap) → validate → /tri-review → /pivot
4. gbF1 恢复&spec 维持 → 提议 M2 freeze + /revise-goal(人闸)进 claim；gbF1 结构封顶 → 关键负结果 pivot 回架构重审(不白烧)

## Completion (inline ✅ CK1-CK6)

CK1 多分类标签+macrochromosome 数据+check_data PASS+比例校验; CK2 多分类 decoder smoke 跑通(head 不 collapse, per-class recall 非零); CK3 >=5-8 seed 全量 COMPLETED; CK4 3-strata(arab/gallus-micro/gallus-MACRO)+CI + SegmentNT 重叠结论; CK5 primary gate inline(gbF1≥0.710或恢复趋势 ∧ spec≥~0.95 ∧ macro 不崩 ∧ gcount≤1.25)+result-log+validate+tri-review+pivot; CK6 iter_ledger 闭合, 下一步(M2/claim 或 架构重审)上交用户。

## Constraints

- **NON-CLAIM**(M2 未冻结+human gate 前不 claim)；**反调参硬闸**：gbF1 gap 0.231≫0.05 → 必须结构性(多分类+decoder)改动, 禁 lr/batch
- 保留确定性 constrained post-proc 作 coherence 层；多分类 head 易 collapse 必查 per-class recall
- macrochromosome stratum + SegmentNT 语料核查 = MANDATORY(result-log 必含); staged unfreeze = 单独后续轴
- 改 ACTIVE_GOAL 只经 /revise-goal 人闸; 不动 docs/03/09; orthogonality N/A
