Complete TA-FOUNDATION-DECODER-M4 · foundation features -> structured decoder, convert recall into intergenic specificity, beat the new-ruler anchor (dual co-primary), as run-and-evaluate. Read and follow `goals/TA-FOUNDATION-DECODER-M4_protocol.md` and `CLAUDE.md`, but judge completion only by the inline checkpoints below.

⚠️ Retrospective advisory: ~5 iterations since (never-run) retrospective; recommended (non-blocking).

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
- 本轮预计 < 半天（cached 特征上训轻 head，无需抽取；15 run = 3 候选×5 seed），第 1 epoch 应 ≤ 数分钟出, 之后切 30/60 min polling
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

Mode: run-and-evaluate；goal 在 batch 的 result-log+tri-review+pivot 后结束（不 scale-up）。Milestone: Track A screen, NON-CLAIM, 主架构赌注（foundation→structured decoder）。3 候选 FP-SEGNT-{FPLOSS,FUSION,CRF}，正交性 **verdict=PASS**（3 不同 major_axis：loss_design/data_view/decoder；详表见 protocol）。
3-layer gate（DUAL co-primary 新尺子）：primary_progress = intergenic_specificity(bw) 严格>锚 0.8710 AND gene_body_F1>=0.5276 AND macro_spec>=0.7978；sota_claim=N/A；review_decision=tri-review 后定 promote-Track-B/iterate/abandon。Bracket：FLOOR 0.8805(blocked)/anchor 0.8710/ceiling 0.9917。

## Hard pre-submit gate

每候选先 1-seed sanity smoke：出 metrics 且**不塌缩**（per_class 非单类），FEATCACHE 对齐 guard(feat L==genome L) 过。任一候选 smoke 崩/塌缩 → 修该候选再进全 batch，不盲投 15 run。

## Required chain

1. /implement 三候选（复用 src/foundation_probe + FEATCACHE，新代码仅 FP-aware loss / raw-DNA⊕特征融合 / CRF 接线）；自审 + 每候选 1-seed sanity smoke。
2. /smart-sbatch 提交 3 候选×5 seed（15 run，private-teodoro-gpu）；job_watch 对账。
3. predict→CDS GFF→eval --span-mode cds→aggregate bw+macro，每候选 seed-mean+CI。
4. validate_goal --profile screen 每候选。
5. /result-log→/tri-review→/pivot；inline 每候选 vs-anchor 双轴 + paired test + per-species(yeast) + gene_count + 决策。

## Completion (inline ✅ CK1-CK5)

CK1 三候选实现 + 各 1-seed sanity smoke 过（不塌缩、出 metrics）；CK2 15 run COMPLETED（job_watch 对账，非假设）；CK3 新尺子 eval 每候选 seed-mean+CI（spec bw+macro、F1、gcount、per-species）；CK4 validate_goal 每候选判定；CK5 result-log+tri-review+pivot 完成 + inline 每候选 vs-anchor(0.8710) 双轴对照 + paired test + 决策。

## Constraints

- NON-CLAIM screen，不 scale-up。复用 FEATCACHE（不重抽取）+ 与锚同 split/协议。
- **≥5 seed + CI + paired test vs 锚**（AXIS-1 方差脆弱：M1 spec 0.808-0.897，单 seed 已超锚）。
- 双门：spec 严格>0.8710 AND F1>=0.5276 AND macro>=0.7978（不达→not_yet）。
- **3-class 输出**（复用 harness）；多类(phase/splice/strand)=紧接 follow-up，非本批。CRF 候选带 FP-aware aux（非裸 semi-CRF）。
- env 只 coding-rna 绝不 base；sbatch 经 ssh baobab。
