Complete TA-COHERENCE-FIX-M5 · fix FPLOSS fragmentation (gene_count 2.25->≤1.25) without losing the specificity win + 5-seed anchor for a valid paired test, as run-and-evaluate. Read and follow `goals/TA-COHERENCE-FIX-M5_protocol.md` and `CLAUDE.md`, but judge completion only by the inline checkpoints below.

⚠️ ① DONE (intergenic_specificity=AXIS-1 primary + full-transcript complement in R6; 0.9303/0.8710 是 NEW-ruler) → 不重跑 /revise-goal。FOCUSED ARCH BATCH on decoder。Retrospective advisory (~6 iter，非阻塞)。

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
- 本轮预计 < 半天（cached 特征训练快 + post-proc 几乎零成本 + 2 个锚 seed），第 1 epoch 应 ≤ 数分钟出, 之后切 30/60 min polling
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

Mode: run-and-evaluate；goal 在 result-log+tri-review+pivot 后结束（不 scale-up；promote=③ 需 go-ahead）。Milestone: Track A screen, NON-CLAIM, FPLOSS 碎片化修复 + 锚补 5 seed。候选 FP-FRAGFIX-{CONSTR(主),CRFSTAB(次)}，正交性 SOFT_WARN focused on decoder。
gate（双 co-primary + coherence, 新尺子）：promotable iff intergenic_specificity(bw) 严格>5-seed 锚 AND gbF1>=0.5276 AND macro>=0.7978 AND **gene_count_ratio 趋向≤1.25**（vs FPLOSS 2.25）；sota_claim=N/A。

## Hard pre-submit gate

每候选 1-seed smoke 不塌缩 + 出 metrics；--postproc constrained 确实降 gene_count（smoke 即见 2.x→更低）。FEATCACHE 对齐 guard 过。崩/无改善 → 修再进全 batch。

## Required chain

1. /implement：train_probe_head 加 --postproc constrained（pred→constrained_decode→GFF，镜像 train_screen_ref）+ CRF 稳定化；自审 + 每候选 1-seed smoke。
2. 锚补 5 seed：run_screen_ref tiberius_like s3 s4，新尺子重 eval → 5-seed 锚 mean+per-seed。
3. /smart-sbatch 提交 FP-FRAGFIX-CONSTR×5(+CRFSTAB×5)（shared-gpu）；job_watch 对账。
4. eval 新尺子+aggregate，每候选 seed-mean+CI；**paired test vs 5-seed 锚**；fragmentation 诊断（gene_count、gene-length 分布、transcript-span P/R）。
5. /result-log→/tri-review→/pivot；inline 每候选 vs-锚双轴 + gene_count 2.25→? + promote/iterate 决策。

## Completion (inline ✅ CK1-CK5)

CK1 候选实现 + 各 1-seed smoke 过（不塌缩、--postproc 见 gene_count 降）；CK2 锚 s3/s4 COMPLETED + 5-seed 锚重算；CK3 候选 5-seed COMPLETED（job_watch）+ 新尺子 eval seed-mean+CI；CK4 paired test vs 5-seed 锚 + fragmentation 诊断；CK5 result-log+tri-review+pivot + inline vs-锚双轴 + gene_count 改善 + 决策。

## Constraints

- NON-CLAIM，不 scale-up（promote=③ 需 go-ahead）。复用 FEATCACHE + 同 split/协议。**不重跑 /revise-goal（①已完成）**。
- ≥5 seed/候选 + CI + **paired test vs 5-seed 锚**（现锚仅 3 seed、一个崩 0.773）。
- 门：spec 严格>锚 AND gbF1>=0.5276 AND macro>=0.7978 AND **gene_count_ratio 趋向≤1.25**（HARD 诊断）。
- λ=1.0 未在 test 调（文档化）；constrained=确定性后处理。env 只 coding-rna；sbatch 经 ssh baobab。
