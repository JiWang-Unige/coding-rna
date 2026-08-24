Complete TA-FRAGFIX-SWEEP-M6 · clear the FP-FRAGFIX-CONSTR gene_count_ratio 1.28->≤1.25 via a deterministic constrained-decode param sweep chosen on VAL (no test leakage), confirming promote-readiness, as run-and-evaluate. Read and follow `goals/TA-FRAGFIX-SWEEP-M6_protocol.md` and `CLAUDE.md`, but judge completion only by the inline checkpoints below.

⚠️ 便宜 STEP-0 gate（promote 前置，M5 tri-review 3/3）。① DONE(R6)，不重跑 /revise-goal。Retrospective advisory(~7 iter，非阻塞)。非 arch batch → 无 orthogonality 表。

## 权限声明

本 goal 拥有充分工具权限。若集群命令 / 网络下载 / subagent 调用临时失败, **先尝试 ≥ 1 替代路径**, 仍失败再 inline 报告, 不要停下等用户：
- 内置工具不受限（Bash / Edit / Write / Read / WebFetch / Glob / Grep / TodoWrite / BashOutput / KillShell 等）
- 可并行调用 Claude Code subagents, 但避免写同一文件
- 允许下载任何研究所需数据（公开数据集 / HF weights / 跨物种基因组 / Rfam family 等）；下载后 inline 报告 path / version / hash / 来源
- sbatch 远程提交训练；若已有可复用 sbatch 脚本, **不要重写**, 仅做 read-only review + Phase 1 policy guard
- 必要时可创建新 skill 以稳定流程（创建后 inline 报告路径+用途）
- 关键 source-of-truth 文件（CLAUDE.md / docs/03 / docs/09）只可 draft patch text, 不直接 Edit

## 运行说明

- 所有 sbatch 提交遵循 /smart-sbatch 自动选合适节点 / 分区 / 时长（Phase 1 guard + Phase 2 optimization 两段）
- 若 sbatch 脚本已存在: 对现有脚本跑 /smart-sbatch Phase 1 guard, **跳过** Phase 2 重生成 header
- 预期 ≤ 1 天的运行用 sleep 同步等待, 不中途下线
- **任何长度的运行**: 先等第 1 个 epoch 完成, 确认 loss 在降 / eval 正常 / 后续 epoch 迭代无问题, **之后**才允许降低监控频率
- 本轮预计 < 半天（5 seed 重训快 + sweep 离线零成本），第 1 epoch ≤ 数分钟, 之后 30/60 min polling
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
- 决策影响 > 24h 净新增 compute spend (例: 重训本身的资源决定可自治; 但启动一条新长 sub-iteration 需确认)
- ≥ 2 reviewer 明确反对 default 共识 (tied 3-way 无 leader)

## Mode & Milestone

Mode: run-and-evaluate；goal 在 result-log+tri-review+pivot 后结束。Milestone: Track A screen, NON-CLAIM, STEP-0 promote-gate（不 scale-up；Track-B=下个 goal 需 go-ahead）。
gate（在 M5 双 co-primary 已过基础上，补 coherence）：选定的 constrained 参数（**在 VAL 上选**）使 **gene_count_ratio ≤ 1.25** AND 仍 intergenic_specificity(bw) 严格>锚（3-seed 0.8710 与 5-seed 0.8436 两个都过）AND gbF1>=0.5276 AND macro>=0.7978；sota_claim=N/A。

## Hard pre-submit gate

--save-raw-pred 正确保存**未经 constrained 的** per-seqid 预测（供离线 sweep）；参数 sweep **只在 VAL 上选**（绝不看 test 的 gene_count/spec），选定后才 apply 到 test。违反 → test 泄漏，停。

## Required chain

1. /implement：train_probe_head 加 --save-raw-pred（存 pre-constrained per-seqid 数组 .npz，含 val+test seqids）+ 自审 + 1-seed smoke。
2. 重跑 5 CONSTR seed（--postproc constrained --save-raw-pred）；shared-gpu；job_watch。
3. 离线 sweep：(max_fill_gap × min_cds_len) 网格在 raw preds 上跑 constrained_decode → **VAL** eval（新尺子）→ 选 gene_count≤1.25 且 spec 最高的参数。
4. 选参 apply 到 **TEST** → 5-seed mean+CI；确认 gene_count≤1.25 AND spec>锚 AND gbF1>=0.5276 AND macro>=0.7978。
5. /result-log→/tri-review→/pivot；inline 选参 + test 双轴+gene_count + promote-ready 判定。

## Completion (inline ✅ CK1-CK5)

CK1 --save-raw-pred 实现+smoke（存 raw + val/test seqids）；CK2 5-seed 重跑 COMPLETED（job_watch）；CK3 VAL param sweep 选定 (max_fill_gap,min_cds_len)（gene_count≤1.25 & spec 最高，非 test 选）；CK4 选参 apply test，5-seed mean+CI，确认 gene_count≤1.25 + spec>锚 + gbF1>=0.5276 + macro>=0.7978；CK5 result-log+tri-review+pivot + promote-ready 判定。

## Constraints

- NON-CLAIM，不 scale-up（Track-B=下个 goal 需 go-ahead）。复用 FEATCACHE + 同 split/协议。不重跑 /revise-goal。
- 参数 **VAL 选、TEST apply**（防泄漏，HARD gate）；constrained=确定性后处理。
- gate：gene_count≤1.25 AND spec 严格>锚 AND gbF1>=0.5276 AND macro>=0.7978（spec↔gene_count 耦合，实测）。
- post-proc 校准非 lr/batch 调参；架构 CONSTR 已固定。env 只 coding-rna；sbatch 经 ssh baobab。
