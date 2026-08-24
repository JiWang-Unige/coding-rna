# TB-UNFREEZE-BACKBONE-M9 · Protocol

> ④ staged-unfreeze PREFLIGHT (promoted from TB-GBF1-MULTICLASS-M8 pivot 2026-06-12, tri-review 2/3 DEGRADED unfreeze-finetune-backbone). USER GO-AHEAD given for the BOUNDED screen. Track B, NON-CLAIM.
> CORE diagnostic: M8 confirmed multi-class does NOT recover gbF1; both reviewers diagnose the gbF1->ANNEVO-ceiling(0.8976) ~0.16 gap as FROZEN features capping (ANNEVO is end-to-end-trained; frozen-head caps ~0.74). This bounded screen MEASURES whether unfreezing the backbone lifts gbF1 — before any >24h commitment.

## Permissions
- 内置工具不受限（Bash/Edit/Write/Read/WebFetch/Glob/Grep/TodoWrite/BashOutput/KillShell）。
- 并行 read-only subagent（Explore/code-plan-reviewer）勘查 + 训练前 review；写入型仅本 exp_id scope。
- 下载 HF weights（NT v2-500m ~2GB 到 HF cache，写共享 FS）；下后 inline 报告。
- sbatch 经 `ssh baobab` 提交；已有可复用脚本只 read-only review + Phase1 guard。
- 改 ACTIVE_GOAL 只经 /revise-goal 人闸；docs/03/09 只 draft，不直接 Edit。
- bounded screen compute 在自治内；**全量 >24h full-unfreeze scale-up 是单独 USER 硬闸**。

## Final goal (milestone)
- Milestone M9, Track B PREFLIGHT (bounded screen), promoted_from M8. profile=screen→full-ish (bounded), **NON-CLAIM** (M2 sota_benchmark pending + human gate).
- 北极星：intergenic_specificity 严格超 published SOTA 且不牺牲 gene-F1、可发表。M9 测「unfreeze backbone 能否抬 gbF1 越过 frozen 天花板」。
- 3-layer gate：
  - **primary_progress_gate**：unfreeze-treatment gene_body_F1（AXIS-2，干净植物 base-w）**方向性 > frozen-backbone control AND > 14-elem frozen 3c 0.7392**，同时 intergenic_specificity 不塌（对照 frozen FPLOSS spec ~0.93），gene_count 不欠预测（不像 M8 mc 0.66）。
  - sota_claim_gate：N/A（NON-CLAIM）。
  - review_decision_gate：tri-review → pivot。
- 参照（干净植物 base-w）：14-elem frozen 3c gbF1 0.7392 / spec 0.9663；raw-DNA anchor gbF1 0.696 / spec 0.905；ANNEVO ceiling gbF1 0.8976 / spec 0.9824。

## Track + resource
- Track B preflight。**FOCUSED ARCH BATCH on training_signal**（同 major_axis，不同 mechanism_delta = freeze depth）：
  - arm1 backbone-FROZEN-embeddings + head（control = backbone-only 自训；隔离 frozen 容量 vs unfreeze）。
  - arm2 backbone-UNFREEZE top-N=2 + head。
  - arm3 backbone-UNFREEZE top-N=4 + head（treatment）。
  - 3 arms × 3 seeds = 9 runs。干净植物 {arabidopsis, rice}，复用 src/screen_anchor/data.py chromosome-level split（same-budget 可比），window 2048，screen 短 epoch。
- 输入 = NT v2-500m backbone hidden states (1024-dim)，**NOT** SegmentNT 14-elem segmentation logits（那是泄漏 head 的输出）。这隔离「可训练 backbone embeddings」vs M8「冻结 segmentation 14-elem」。
- 架构: `AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True).esm` 取裸 backbone（29 层/1024 hidden/16 heads/rotary；丢 LM head；**不碰泄漏的 SegmentNT segmentation U-Net**）→ 冻结全部 → `esm.encoder.layer[-N:].requires_grad_(True)` → 接现有 train_probe_head 的 3c FP-aware constrained head（per-base 3-class + intergenic-FP penalty + constrained post-proc mfg=20/mcl=60）→ loss.backward() 进顶 N 层 + head。
- 显存: bf16 AMP + `esm.gradient_checkpointing_enable()` + 冻结层 requires_grad=False；RTX3090 24GB 够（顶 N 层 + window 2048）。首跑下 NT v2-500m ~2GB 到 HF cache。env=coding-rna（transformers 5.11.0 + torch 2.5.1，依赖就绪）；若 trust_remote_code 在 5.11.0 炸 → fallback generanno(4.49.0, 有 peft/accelerate)。
- 资源: private-teodoro-gpu(空了)或 shared-gpu `--constraint=COMPUTE_TYPE_AMPERE`（≥24GB）+ exclude 3080。

## Execution mode details
- submit-and-handoff：训练 handoff 远程；决策前台可见。job_id/output 写 docs/05；后台 monitor poll squeue，不假设成功。第 1 epoch 出后确认 loss 降 + backbone grad 流动 + head 不 collapse，再降频。

## Pre-submit gate (HARD)
1. **兼容性 smoke（先做，5min）**：在 coding-rna `python -c` 验 `AutoModelForMaskedLM.from_pretrained(NT_v2_500m, trust_remote_code=True)` 能 import+load + 取 `.esm` + forward 出 (B,L,1024) hidden。失败 → fallback generanno 或装 peft；仍失败 → 停报 blocker。
2. **tokenizer 对齐**：NT 是 6-mer tokenizer（1 token≈6bp）。head 输出是 per-base（1bp）——必须处理 token→base 对齐（6× 上采样/重复 或 per-token label 聚合）。这是关键实现点，smoke 必须验证输出长度与 label 长度对齐（per-base 3-class）。
3. **check_data 无泄漏**：干净植物 split（沿用 M8 split）无 seqid 跨 split 泄漏。
4. **sanity smoke**（小规模真跑）：1 arm（unfreeze N=2）1 seed 2 epoch，确认端到端跑通 + backbone 顶层 grad 非零 + head 不 collapse + gbF1 sane + 显存不 OOM。失败 → repair_advisor 有界修复 ≤3 次。
任一未过 → 不提交全量；final pivot = fix_eval/fix_data。

## Pre-submit code review
新 trainer（torch backbone 加载 + freeze/unfreeze + token-base 对齐 + AMP + grad checkpoint）是非平凡新代码 → 建议先 plan 模式勘查 `modeling_esm.py`（HF cache 里）接口 + train_probe_head head 模块复用点，再落码。提交前 1 个 read-only subagent（code-plan-reviewer）核对：token-base 对齐正确性（最易错）、freeze mask、grad 流动、metric/split 可比、无 target 泄漏。

## Subagent fan-out
- read-only：勘查 NT modeling_esm 接口（hidden states 提取、gradient_checkpointing、6-mer tokenizer 对齐）；勘查 train_probe_head head 模块复用点。
- 写入型：仅本 exp_id scope（TB-UNFREEZE-BACKBONE-M9-*）。主 agent merge + 最终 sbatch + 写 docs。subagent 不再 spawn。

## Pivot decision menu (preflight)
- **unfreeze 方向性抬 gbF1**（treatment > control > 14-elem 3c，spec 不塌，gcount sane）→ **frozen 确实封顶**确证 → 下一步 full-unfreeze scale-up（>24h，单独 USER 硬闸；更多层/epoch/seed + 全量数据 + 多 clean 物种）→ 朝 ANNEVO ceiling。
- **unfreeze 抬 gbF1 但 spec 塌** → spec↔gbF1 张力在 unfreeze 下重现 → 加 FP-aware 权重 / 分阶段（先 head 后 unfreeze）调和；仍架构轴。
- **unfreeze 不抬 gbF1**（≈ frozen control）→ **frozen 不是瓶颈**，关键负结果 → pivot 回：换 foundation model（更高分辨率/植物友好）/ self-sup domain-adapt / 重审 co-primary 是否可兼得 / 怀疑 head 容量。写 docs/08。
- **OOM/兼容性炸** → fallback generanno env / LoRA(peft) / 减 window / 减 N。

## Skill invocation chain
| step | skill | 说明 |
|---|---|---|
| (done) | /goal-prompt | 本文件 |
| 1 | /implement | 新 torch-backbone-unfreeze trainer（加载 .esm + freeze/unfreeze + token-base 对齐 + AMP + grad ckpt + 复用 3c head）；兼容性 smoke + check_data + sanity smoke |
| 2 | /smart-sbatch | Phase1 guard(orthogonality SOFT_WARN focused training_signal) + Phase2；3 arm × 3 seed batch matrix |
| 3 | submit-and-handoff | ssh baobab + job_watch + 后台 monitor |
| 4 | /result-log | 3 arm 对比（frozen control vs unfreeze N2/N4）+ vs 14-elem 3c 0.739 + spec/gcount，干净植物 base-w + CI |
| 5 | validate → /tri-review → /pivot | unfreeze 抬 gbF1 否？spec 塌否？下一步 full scale-up(user gate) 或负结果 pivot |
| 6 | /exp-log + iter_ledger | 链路闭合 |

## Constraints (full)
- **NON-CLAIM**（M2 未冻结 + human gate 前不 claim）。
- **LEAKAGE 纪律（HARD）**：只在 segmentation-clean 物种（植物 {arabidopsis,rice}）评估——chicken/fly 是 SegmentNT segmentation fine-tune 物种、永远 contaminated，不进 M9。test labels 不进 early-stop/decode-tuning；保持 raw-DNA ab-initio（不引入 target-derived/evidence features）；同 3-class collapse ruler（eval 不变）；gene_count 不用 test truth 调。
- **反调参**：gbF1 gap 0.16≫0.05 → training_signal（可训练 backbone）是架构轴，非 lr/batch 调参 → 合规。
- **token-base 对齐**是关键实现正确性点（NT 6-mer vs per-base label）——smoke 必验。
- bounded screen 在自治内；**全量 >24h full-unfreeze scale-up = 单独 USER 硬闸**。
- 改 ACTIVE_GOAL 只经 /revise-goal；不动 docs/03/09。非 abandoned cousin（docs/09 空）；与 M1-M8 差异：首次让 foundation backbone 可训练（之前全 frozen）。
- FOCUSED ARCH BATCH on training_signal（frozen control + unfreeze N2/N4）——不是 diverse batch，如实标注。
