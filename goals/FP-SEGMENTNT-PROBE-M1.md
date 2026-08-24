Complete FP-SEGMENTNT-PROBE-M1 · foundation-probe: frozen SegmentNT element logits → light head, same-budget screen, vs the recalibrated anchor under the new dual co-primary ruler, as run-and-evaluate. Read and follow `goals/FP-SEGMENTNT-PROBE-M1_protocol.md` and `CLAUDE.md`, but judge completion only by the inline checkpoints below.

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
- 本轮预计 < 半天（特征抽取一次性缓存 + 轻量 head screen）, 第 1 epoch 应 ≤ 30min 内出, 之后切 30/60 min polling
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

Mode: run-and-evaluate；goal 在 screen 结果的 result-log+tri-review+pivot 后结束（本轮不 scale-up）。Milestone: Track A screen, NON-CLAIM foundation-probe（首个 post-ruler 架构动作），profile=screen，exp_id 命名空间 FP-SEGMENTNT-*。
3-layer gate（DUAL co-primary, 新尺子）：primary_progress = intergenic_specificity(bw) 严格 > 锚 0.8710 AND gene_body_F1 >= 0.5276 AND macro_spec >= 0.7978；sota_claim = N/A；review_decision = tri-review 后定 promote/iterate/abandon。

## Hard pre-submit gate

ENV 共存 + 泄漏 smoke 必过：JAX 见 GPU 且与 torch 2.5.1 同 env 共存、SegmentNT 在极小输入上产出 (B,seq*6,14,2) logits；特征缓存经 check_data 无 train/test seqid 泄漏（与 anchor 同 split）。共存失败 → JAX-CPU fallback；仍失败或检出泄漏 → 不提交，pivot=fix_env_or_data。

## Required chain

1. 装 jax[cuda12]+dm-haiku+einops+huggingface_hub 进 coding-rna（非 base），smoke JAX+torch 共存 + SegmentNT 加载（优先 multi_species）。
2. /implement 特征抽取：SegmentNT 冻结 14-logits（'N' mask）→缓存，复用 anchor 窗口+同 split；check_data 泄漏闸。
3. /implement 轻量 head（linear/MLP）14-logits→3-class，与 train_screen_ref 同预算（2048/0.3/8ep/pat3/3seed/CW-CE）；sanity smoke。
4. /smart-sbatch 提交 3-seed screen（private-teodoro-gpu）；job_watch 对账。
5. predict→CDS GFF→eval --span-mode cds→aggregate bw+macro→validate_goal --profile screen。
6. /result-log→/tri-review→/pivot；报 seed-wise spec(bw+macro)+F1+gene_count+per-species，对照 ladder 0.8710/0.8369/0.7138。

## Completion (inline ✅ CK1-CK6)

CK1 env+JAX/torch 共存 smoke 过 + SegmentNT 加载（报变体）；CK2 特征缓存+check_data=pass；CK3 light-head sanity 出 metrics；CK4 3-seed COMPLETED（job_watch 对账）；CK5 新尺子 eval 出 spec(bw+macro)+F1+gene_count + validate_goal 判定；CK6 result-log+tri-review+pivot 完成 + inline vs-anchor 双轴对照 + 决策。

## Constraints

- NON-CLAIM；不 scale-up。Same-budget：复用 anchor 同 split+同 test subsets+同冻结协议，**不改 split**。
- 环境纪律：只装 coding-rna，**绝不 base**；sbatch 经 ssh baobab，脚本内 conda activate coding-rna。
- 新尺子 eval --span-mode cds；报 spec bw+macro + F1 地板 0.5276 + macro 门 0.7978（不达→not_yet 非 failed_run）。
- 'N' 预 mask；SegmentNT 人类训练，跨 yeast/fly 迁移本身是 probe 问题，优先 multi_species。
