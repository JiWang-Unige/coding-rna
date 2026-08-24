# FP-SEGMENTNT-PROBE-M1 · Protocol

> 执行手册（goal 运行期间读）。短 prompt 是 evaluator 唯一可见契约；本文件是操作细节。

## Permissions
- 全工具可用（Bash/Edit/Write/Read/WebFetch/Glob/Grep/TodoWrite/BashOutput/KillShell + Claude Code subagents）。
- 写 scope：本 exp_id 命名空间 `FP-SEGMENTNT-*` 的 configs/sbatch/src/outputs；新增的特征抽取 + light-head 代码放 `src/foundation_probe/`（新模块，不动 `src/screen_anchor/` 公共训练码，除非必要且声明）。
- 关键 source-of-truth（CLAUDE.md / docs/03 / docs/09）只 draft patch，不直接 Edit。ACTIVE_GOAL.json 不动（本轮非 revise-goal）。
- 装依赖只进 conda env `coding-rna`，**绝不 base**（CLAUDE §12 环境纪律）。

## Final goal（北极星对齐）
终极目标：在 primary（DUAL co-primary：AXIS-1 intergenic_specificity + AXIS-2 gene-level F1）上严格超越 frozen published SOTA 且可发表。本轮是**非 claim screen probe**：回答"用 SegmentNT 预训练元件信号作输入特征，能否在同预算下提升 intergenic_specificity 而不丢 gene_body_F1"——即 R6 finding 的核心假设（更好的输入特征是降 intergenic FP 的真正杠杆，而非结构化 decoder）。
- screen_anchor（PROVISIONAL，新尺子）：intergenic_specificity base-weighted **0.8710** / macro **0.8278**；anchor gene_body_F1 = 0.5576。
- same-budget ladder（新尺子，对照）：tiberius_like(锚) 0.8710 / CONSTR 0.8369 / helixer_like 0.7954(frag) / CRF-vec 0.7138。FLOOR(ORF) 0.8805 但 F1 0.3735 被地板挡。
- 双轴门（screen direction-selection）：intergenic_specificity(bw) 严格 > 0.8710 **AND** gene_body_F1_unconstrained >= 0.5276 **AND** macro_intergenic_specificity >= 0.7978。
- screen **永不 claim**；本 goal 不 scale-up（promote 留给后续 goal）。

## Track + resource
- Track A screen，单候选（SegmentNT logits-as-features）。profile=screen。
- 冻结协议（与 `src/screen_anchor/train_screen_ref.py` 同预算，保证 same-budget 可比）：window 2048、sample_fraction 0.3、epochs 8、patience 3、3 seeds（s0/s1/s2）、class-weighting sqrt_inv。
- 物种：yeast (S.cerevisiae) + fly (D.melanogaster)，**复用 anchor 的 chromosome-level held-out split + 同 eval_subsets test 区**（`outputs/SCREENREF-tiberius_like-s*/eval_subsets/<sp>/{genome.fa,reference.gff3}`），不改 split。
- 资源：private-teodoro-gpu，RTX3090 24GB。SegmentNT 500M ~2.2GB（fp32）；light-head 训练极小。expected walltime < 半天 → run-and-evaluate。
- per-seed 输出目录：`outputs/FP-SEGMENTNT-PROBE-M1-s{0,1,2}/{logs,checkpoints,metrics,predictions}`；特征缓存 `outputs/FP-SEGMENTNT-FEATCACHE/<sp>/`（seed 无关，共享）。

## Execution mode details (run-and-evaluate)
- 特征抽取是一次性（seed 无关）→ 先抽取缓存，再 3 seed 各训 light-head。
- 先等 seed s0 第 1 epoch 出（应 ≤ 30min；light-head 上 cached features 极快），确认 loss 降 + eval 正常，再降频 30/60min polling 其余 seed。
- job_watch 对账（sacct/squeue → COMPLETED/FAILED/TIMEOUT/OOM/STALE，不假设成功）；幽灵 run 按 failed_run。

## ENV 安装 + 共存 smoke（CK1，Hard pre-submit 的一半）
1. 经 ssh baobab 进 coding-rna 装（inline 报装什么/装哪）：
   `pip install "jax[cuda12]" dm-haiku einops huggingface_hub`（jaxlib 随 jax[cuda12] 带 cuda wheel；numpy 须 <2.0，SegmentNT setup 要求——若冲突，pin `numpy<2`）。
   - 注意 jax cuda wheel 与 torch 2.5.1+cu121 共存：两者各自带 CUDA runtime，通常可共存；若 XLA 抢显存，设 `XLA_PYTHON_CLIENT_MEM_FRACTION=0.4` 或 `XLA_PYTHON_CLIENT_PREALLOCATE=false`。
2. 共存 smoke（在 GPU 节点经 srun 或 sbatch smoke，**不在登录节点**）：import torch + import jax 同进程；`jax.devices()` 见 gpu；`get_pretrained_segment_nt_model(model_name="segment_nt_multi_species", max_positions=...)` 加载成功；对一条 ~2046bp 随机 ACGT 序列 forward 出 `logits` shape `(1, seq*6, 14, 2)`。
3. **fallback（可自治）**：若 GPU 下 JAX/torch 共存崩（CUDA 库冲突/OOM），改 JAX-**CPU** 一次性抽取（`JAX_PLATFORMS=cpu`）——500M 模型在 2 物种 test 区窗口上 CPU 抽取慢但一次性可接受；light-head 训练仍用 torch GPU。CPU 也跑不动 → inline 报 blocker，pivot=fix_env。

## Pre-submit gate（CK1+CK2，全 Hard）
- ENV 共存 smoke 过（上节）。
- 特征缓存 check_data 泄漏闸：`python3 scripts/check_data.py`，确认特征/标签的 train/val/test 按 anchor 同 chromosome-level split，无同 seqid 跨 split；`status=leakage` → 停，修 split。
- 数据契约：cached logits 与逐碱基标签**对齐**（6-mer unfold → base resolution；窗口边界对齐；'N' 位置 mask 且标签忽略）；schema/dtype 匹配；target 分布合理（intergenic/CDS/intron 比例记入 notes）；预处理只在 train fit。任一 ❌ → 不提交。

## 实现要点（CK2-CK3，/implement）
- **特征抽取器** `src/foundation_probe/extract_segmentnt.py`：
  - 加载 `segment_nt_multi_species`（跨支系更稳；同时记录是否也试 `segment_nt`）。
  - 输入：anchor 的 frozen 窗口（2048bp）。'N' → 预 mask（替换为占位 + mask 数组，标签处 ignore_index）。
  - 6-mer tokenize → forward（jax.jit/pmap）→ output `logits (B, seq*6, 14, 2)`；取 softmax 后 feature-present 概率（axis=-1 的 idx1）→ 14 维/碱基特征向量。可选同时存 layer-29 1280-d embedding 作 ablation（默认只存 14-logits，省空间）。
  - 缓存为 `.npy`/`.npz` per 物种 per split，含 base-resolution 特征 + 对齐的标签 + mask。
  - 关注 feature index：CDS=protein_coding_gene[0]、exon[2]、intron[3]、splice_donor[4]、splice_acceptor[5]、5UTR[6]、3UTR[7]（其余 regulatory 也带上，让 head 自己学权重）。
- **Light head** `src/foundation_probe/train_probe_head.py`：
  - 输入 = 14 维（或 +embedding ablation）/碱基；模型 = linear 或 1-hidden MLP（小，per-base 分类）；输出 3-class（intergenic/CDS/intron），class-weighted CE（sqrt_inv），mask 掉 'N'/pad。
  - 复用冻结协议超参（2048/0.3/8ep/pat3）；3 seeds。产出 metrics JSON（含 primary_metric 键，eval 会覆盖）。
  - sanity smoke：sample 极小 + 1-2 epoch 跑通到出 metrics；失败用 `repair_advisor.py` 有界修复（≤3 次）。
- **predict→eval**：复用 `src/screen_anchor/gff_io.py` 出 CDS GFF；`scripts/eval_gene_body_mask.py --span-mode cds`（新 full-transcript intergenic 尺子）per 物种 → `scripts/aggregate_gene_body_metrics.py`（bw+macro）→ `scripts/validate_goal.py --profile screen`。

## sbatch（CK4，/smart-sbatch）
- 经 ssh baobab 提交；脚本内 `source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh && conda activate coding-rna`。
- 若特征抽取需 GPU+JAX，可单独一个 extract job（或在同脚本先抽取再训）；3 seed 可 array 或 3 条 sbatch，各独立 output_dir。
- /smart-sbatch Phase 1 policy guard 必跑（VRAM≥20GB、output_dir 唯一、partition time、维护窗口）；CPU-only 抽取走 fast path。

## Orthogonality declaration
单候选（非 ≥2 并行 batch）→ 正交性表 N/A。mechanism delta vs 既有 ladder 清晰：**输入信号轴**（major_axis = training_signal/data_view）——输入从 raw one-hot DNA（tiberius_like/CRF-vec backbone 学的）换成 **SegmentNT 预训练元件 logits**，head 保持轻量 per-base，**非 decoder 改动**。why structural：新的输入张量来源（外部预训练 foundation 特征），forward 路径与梯度源全变。

## Subagent fan-out
- 可选 read-only `code-plan-reviewer`：提交前审 extract/head 代码的标签对齐 + 'N' mask + split 复用 + metric 口径（防 ground-truth 用错）。
- 不让 subagent 写同一文件；主 agent 负责 merge + 最终 sbatch + 写 docs。

## Pivot decision menu（track-A-screen，CK6）
- **promote-direction**：intergenic_specificity(bw) 严格 > 0.8710 AND F1>=0.5276 AND macro>=0.7978 → foundation 特征确实降 intergenic FP 不丢 recall → 下一 goal 接 semi-CRF + FP-aware objective（仍非 claim，screen→full 规划）。
- **iterate-probe**：接近但未过（如 spec 略高但 F1 跌破地板，或 multi_species 不如 human 变体）→ 调输入特征用法（logits vs embedding、加 regulatory、人类变体 vs multi_species）再来一轮 screen。
- **abandon-probe（写 docs/09 需人闸）**：SegmentNT 特征对 yeast/fly 跨支系完全不迁移（spec ≈ raw-DNA 或更差）→ 记 finding，转 GENERanno probe 或纯 from-scratch 路线。
- **fix_env_or_data**：env 共存/泄漏 gate 失败 → 先解 blocker。
- 多选项决策点按短 prompt 决策自治：列选项→3 CLI reviewer→共识+ROI 选→写 docs/08→不暂停（破坏性/abandon/>24h/tied 例外才停）。

## Skill invocation chain
| step | skill | 状态 |
|---|---|---|
| env+impl | /implement（含 check_data + sanity smoke + 有界 debug） | 用 |
| submit | /smart-sbatch（Phase1 guard + Phase2） | 用 |
| 对账 | scripts/job_watch.sh | 用 |
| 评估 | eval_gene_body_mask --span-mode cds + aggregate + validate_goal | 用 |
| 记录 | /result-log（写 06+04+05+00） | 用 |
| 复核 | /tri-review（3 CLI：claude+codex+agy） | 用 |
| 决策 | /pivot | 用 |
| 实验档 | /exp-log（写 docs/experiments + build_atlas） | 用 |
| /reproduce-baselines | — | 跳过（segmentnt-2025 dossier 已验证） |
| /retrospective | — | 跳过（未触发） |

## Constraints (full)
- NON-CLAIM：screen profile 永不 claim SOTA；本 goal 不 scale-up（promote 留下个 goal）。
- Same-budget apples-to-apples：复用 anchor 同 chromosome-level split + 同 test subsets + 同冻结协议；**不改 split**，否则与 ladder 不可比。
- 环境纪律：只装 coding-rna，**绝不 base**；装前 inline 说明装什么/装哪；sbatch 脚本显式 conda activate coding-rna（计算节点不继承登录态）。
- 新尺子：eval --span-mode cds（intergenic = full-transcript 含 UTR 补集）；报 intergenic_specificity bw+macro + gene_body_F1 + gene_count_ratio + per-species + paired-seed spread vs anchor。F1 地板/macro 门不达 → not_yet（非 failed_run）。
- 'N' 预 mask（6-mer tokenizer N→unk 会坏）；SegmentNT 主体人类 GENCODE 训练，对 yeast(真菌)/fly(昆虫) 跨支系迁移本身是 probe 问题——优先 segment_nt_multi_species，必要时对比 human 变体。
- 报告必须含 seed-wise（非只均值）+ per-species（base-weighted 会掩盖单物种崩溃）+ paired spread（CRF-vec 上轮方差极大的教训）。
