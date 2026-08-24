# TA-FRAGFIX-SWEEP-M6 · Protocol

> STEP-0 promote-gate（M5 tri-review 3/3 一致要求）：把 FP-FRAGFIX-CONSTR 赢家的 gene_count_ratio 1.28 清到 ≤1.25，参数**在 VAL 上选**（防 test 泄漏），确认 promote-ready。便宜（5 seed 重训 cached 特征 + 离线 CPU sweep）。然后 Track-B 主战是**下一个 goal**（需 go-ahead）。

## Permissions
- 全工具。写 scope：`FP-FRAGFIX-CONSTR-s*`（重跑覆盖，带 raw-pred）+ `src/foundation_probe/train_probe_head.py` 加 `--save-raw-pred` + 一个离线 sweep 脚本 `scripts/_sweep_constrained_m6.py`。复用 `src/screen_anchor/decoders.constrained_decode`、eval/aggregate。
- 不动 ACTIVE_GOAL（①已完成；可选的 5-seed anchor /revise-goal 是 human-gated 另议）。装依赖只 coding-rna（本轮无需新依赖）。

## Final goal（北极星对齐）
M5：FP-FRAGFIX-CONSTR（FPLOSS + constrained 后处理）paired-显著 Pareto 超 5-seed 锚（spec 0.9272 vs 0.8436，paired +0.0836±0.037，gbF1 0.6581，macro 0.8555），唯一尾巴 **gene_count_ratio 1.28 > full/scale guardrail 1.25**（差 0.03）。本轮用 VAL 选定 constrained 参数把 1.28→≤1.25 且不丢 spec，确认赢家 promote-ready。
- gate（promote-ready）：选定参数下 test 上 gene_count≤1.25 AND intergenic_specificity(bw) 严格>锚（两个锚都过：3-seed 0.8710 / 5-seed 0.8436）AND gbF1>=0.5276 AND macro>=0.7978。
- bracket：FLOOR 0.8805 / 锚 0.8436(5seed) / CONSTR 0.9272(待清碎片) / ceiling 0.9917。

## Track + resource
- Track A screen，单一机制（CONSTR + constrained 参数校准）。profile=screen。无 arch batch → 无 orthogonality 表。
- 冻结协议同 M5（window 2048, sample 0.3, 8 epoch, patience 3, CW-CE, conv+biLSTM head, --loss fp_aware --postproc constrained, yeast+fly, 同 split）。复用 FEATCACHE。
- 资源：shared-gpu（private 被 M2/zebrafish 占；head 极小）。重跑 5 seed（~10min/seed P100）。sweep 离线 CPU（olympus /usr/bin/python3 或 baobab，纯 stdlib eval）。run-and-evaluate。

## 实现要点（CK1，/implement）
- **--save-raw-pred**：在 train_probe_head 的 predict 里，对 **val + test** seqids 都做 raw（pre-constrained）per-base argmax 预测，存 `outputs/<exp>/raw_pred/{split}_{species}.npz`（{seqid: int8 array}）。注意：当前 predict 只跑 test；加一个 val 预测分支（用 val seqids，同 batched predict）。constrained 仍照常 apply 到 test 出 GFF（不变）。
- 自审：raw pred 是 argmax(emissions)（CRF 时是 viterbi，但本候选 decoder=none → argmax），未经 constrained；val/test seqids 都存；对齐正确。1-seed smoke 确认 npz 存出 + 含 val/test。

## 重跑 5 seed（CK2）
- 重跑 FP-FRAGFIX-CONSTR s0-4（覆盖现有；脚本加 --save-raw-pred）。或新 exp-id `FP-FRAGFIX-CONSTR-rp-s*` 存 raw（避免覆盖 M5 已 result-log 的产物——推荐新 exp-id，保 M5 可复现）。shared-gpu。job_watch。

## 离线 param sweep（CK3，核心，防泄漏）
- `scripts/_sweep_constrained_m6.py`：对每 seed 的 **raw val preds** 跑 constrained_decode 网格，例如 max_fill_gap ∈ {20,40,60,100,150} × min_cds_len ∈ {30,60,90}，每组 → val GFF → eval_gene_body_mask --span-mode cds（VAL 的 reference/genome subset，需生成 val eval_subsets）→ 取 5-seed 平均 gene_count_ratio + spec。
- **选参规则（在 VAL 上）**：在 gene_count_ratio≤1.25 的参数中，选 intergenic_specificity 最高的一组（5 seed 一致）。打印完整网格表（透明）。**绝不看 test。**
- val eval_subsets：用 D.write_subset_fasta/gff 对 val seqids 生成（assign_splits 的 val），或扩展现有 eval 流程。

## apply 到 test + 终评（CK4）
- 用 VAL 选定的 (max_fill_gap, min_cds_len) 对 5 seed 的 **test raw preds** 跑 constrained_decode → test GFF → eval → aggregate（bw+macro）→ 5-seed mean+CI。
- 确认 gate：gene_count≤1.25 AND spec 严格>锚（0.8710 & 0.8436）AND gbF1>=0.5276 AND macro>=0.7978。spec 与 gene_count 耦合——若清碎片把 spec 压到 ≤锚，记录 tradeoff，pivot=iterate（调网格/换策略）。

## Pre-submit gate（CK1+CK3，HARD）
- --save-raw-pred 存出 raw + val/test seqids（smoke 验证）。
- 参数 **只在 VAL 选**（test 仅最后 apply 一次）。任何"看 test 调参"= 泄漏 blocker。

## Pivot decision menu（CK5）
- **promote-ready**：选定参数 test 上 gene_count≤1.25 + spec>锚 + gbF1>=floor + macro>=gate → FP-FRAGFIX-CONSTR 正式 promote-ready → **下个 goal = Track-B scale**（go-ahead）：扩数据/epoch/seed + 多类(CDS/intron/intergenic/phase/splice) + 解冻 SegmentNT（分阶段）。
- **iterate**：清碎片把 spec 压到 ≤锚（耦合）→ 调网格 / 试软 merge / 接受略高 gene_count 但记录；或 CRFSTAB。
- **fix**：raw-pred bug / val subset 问题 → 修。
- 多选项点按短 prompt 决策自治（破坏性/abandon/>24h/tied 例外才停）。**Track-B 启动需 go-ahead**。

## Skill invocation chain
| step | skill | 状态 |
|---|---|---|
| impl | /implement（+ sanity smoke） | 用 |
| submit | /smart-sbatch（Phase1 guard）| 用 |
| 对账 | scripts/job_watch.sh | 用 |
| sweep | scripts/_sweep_constrained_m6.py（VAL 选参，离线）| 用 |
| 评估 | eval_gene_body_mask --span-mode cds + aggregate（test apply）+ CI | 用 |
| 记录/复核/决策 | /result-log → /tri-review → /pivot | 用 |
| 实验档 | /exp-log + build_atlas | 用 |
| /revise-goal | — | 跳过（①已完成；5-seed anchor 更新另议）|
| /retrospective | — | advisory（~7 iter，非阻塞）|

## Constraints (full)
- NON-CLAIM screen，不 scale-up（Track-B=下个 goal，go-ahead）。复用 FEATCACHE + 同 split/协议。不重跑 /revise-goal。
- 参数 **VAL 选、TEST apply**（防泄漏，本轮 HARD gate）；constrained_decode 确定性。
- gate：gene_count≤1.25 AND spec 严格>锚（两锚都过）AND gbF1>=0.5276 AND macro>=0.7978。spec↔gene_count 耦合，实测不假设。
- 这是 post-proc 校准（清 coherence guardrail）非架构/lr 调参；架构 CONSTR 固定。
- 推荐新 exp-id（如 FP-FRAGFIX-CONSTR-rp-s*）保 M5 产物可复现，不覆盖。
- env 只 coding-rna 绝不 base；sbatch 经 ssh baobab，脚本显式 conda activate。
- 预 claim 闸（未来）：核实测试 clade 不在 SegmentNT 预训练语料。
