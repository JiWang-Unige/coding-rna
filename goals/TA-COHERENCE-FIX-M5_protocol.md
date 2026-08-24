# TA-COHERENCE-FIX-M5 · Protocol

> 执行手册。短 prompt 是 evaluator 唯一可见契约。本轮 = TA-FOUNDATION-DECODER-M4 的 pivot follow-up：把 FPLOSS 赢家的碎片化(gene_count 2.25)修到 ≤1.25 而不丢 specificity，并把锚补到 5 seed 做合法 paired test，让赢家真正可 promote。

## Permissions
- 全工具。写 scope：`FP-FRAGFIX-*` + 锚 seed `SCREENREF-tiberius_like-s{3,4}` 的 outputs；`src/foundation_probe/train_probe_head.py` 加 `--postproc` + CRF 稳定化选项（本模块自有）。复用 `src/screen_anchor/{train_screen_ref.py,run_screen_ref.sbatch,decoders.py}`（只读复用，不改）。
- 不动 ACTIVE_GOAL（①已完成，0.9303/0.8710 是新尺子——已实证：FPLOSS metrics.json `intergenic_definition=complement_of_full_transcript_span_incl_UTR`）。CLAUDE/docs03/09 只 draft。装依赖只进 coding-rna（本轮无需新依赖）。

## Final goal（北极星对齐）
TA-FOUNDATION-DECODER-M4 证明 foundation 特征 + FP-aware loss(FPLOSS) 在 screen 上 Pareto 超锚（spec 0.9303>0.8710, gbF1 0.6157>0.5576, macro 0.8431>0.7978），**唯一硬伤 gene_count_ratio 2.25**（full/scale HARD guardrail ≤1.25 会挡）。本轮修这个硬伤 + 建立合法 paired 对照，使赢家可 promote 到 Track B（=③，需 go-ahead）。
- 新尺子 bracket：FLOOR 0.8805 / **锚 0.8710（待 5-seed 重算）** / FPLOSS 0.9303(碎片 2.25) / ceiling 0.9917。
- 关键洞察（M4 tri-review 后蕾姆补）：CRF 候选**本就 = FP-loss + CRF decoder**，它拿 spec 换 coherence（spec 0.830/连贯 0.90/方差大，一个 seed 崩 0.593）。所以"合成"已试过——本轮改用**便宜的确定性 constrained-decode 后处理**修碎片（避开学习型 CRF 的不稳定），CRF 稳定化作次选。

## Track + resource
- Track A screen，2 候选（FP-FRAGFIX-CONSTR 主 / FP-FRAGFIX-CRFSTAB 次）× 5 seed + 锚补 2 seed。profile=screen。
- 冻结协议同 M4/锚（window 2048, sample 0.3, 8 epoch, patience 3, CW-CE, conv+biLSTM head, yeast+fly, 同 chromosome split）。复用 FEATCACHE。
- 资源：shared-gpu（private 被用户 M2 job + zebrafish 占；head 极小 ~2GB VRAM，任何 GPU 足；--time ≤12h，每 run <1h）。锚 seed 用 run_screen_ref（raw-DNA，也轻）。run-and-evaluate。
- 输出：`outputs/FP-FRAGFIX-{CONSTR,CRFSTAB}-s{0..4}` + `outputs/SCREENREF-tiberius_like-s{3,4}`。

## Parallel batch matrix
| exp_id | major_axis | mechanism_delta | seeds | flags |
|---|---|---|---|---|
| FP-FRAGFIX-CONSTR-s{0..4} | decoder/post-proc | FPLOSS + 确定性 constrained_decode 后处理（merge 小 intergenic gap / drop 小 CDS run）修碎片 | 5 | --loss fp_aware --postproc constrained |
| FP-FRAGFIX-CRFSTAB-s{0..4} | decoder | 稳定化学习型 CRF（warm-start / 正则 / 更强 fp-lambda）恢复 spec 保 coherence 0.90 | 5 | --decoder crf --loss fp_aware + stabilization |
| SCREENREF-tiberius_like-s{3,4} | (anchor validity) | 补锚到 5 seed，新尺子重 eval，做 paired test | 2 | run_screen_ref tiberius_like {3,4} |

正交性 verdict = **SOFT_WARN（focused-arch-batch on decoder）**：2 候选同 major_axis=decoder，mechanism 不同（确定性后处理 vs 学习型 CRF 稳定化），均 structural。锚 seed 非候选（validity job）。

## 实现要点（CK1，/implement）
- **--postproc constrained**（主）：train_probe_head predict 后、labels_to_cds_gff 前，对 per-seqid pred 数组跑 `from src.screen_anchor.decoders import constrained_decode`（已 import 过 LinearChainCRFVec，加 constrained_decode）。镜像 train_screen_ref 的 constrained 路径（min_cds_len/max_fill_gap 参数）。仅后处理，训练不变（仍 --loss fp_aware）。
- **CRF 稳定化**（次）：现 CRF 候选 spec 0.830 + 一个 seed 崩 0.593。稳定化选项：(a) warm-start——先 N epoch 只训 head(emissions) 再开 CRF；(b) CRF transition 正则/init；(c) 更强 fp-lambda。先实现最简单的 warm-start + 看是否消崩。若复杂，本轮可只交 CONSTR，CRFSTAB 留下轮（inline 说明）。
- 自审：label 对齐、constrained_decode 用对 per-seqid 数组、metric 口径同锚。每候选 1-seed smoke（2 epoch）：出 metrics、**确认 --postproc 使 gene_count_ratio 从 ~2.x 降**、不塌缩。

## 锚补 5 seed（CK2）
- `ssh baobab` 提交 run_screen_ref.sbatch tiberius_like 3 / tiberius_like 4（softmax，同冻结协议）。
- 完成后用 scripts/recompute（_recompute 思路 / eval_gene_body_mask --span-mode cds）对 s0-s4 全 5 seed 重 eval 新尺子 → 锚 5-seed mean + per-seed intergenic_specificity。当前 3-seed: 0.923/0.917/0.773 (mean 0.871)；补 s3/s4 后取 5-seed mean 作 paired 对照基准。

## Pre-submit gate（CK1，HARD）
- 2 候选各 1-seed smoke 过（不塌缩、出 gene_count_ratio、--postproc 见下降）。
- FEATCACHE 对齐 guard（feat L==genome L）。数据契约：constrained_decode 输入是 int8 per-seqid 数组（同 gff_io 期望）。任一 ❌ 不提交。

## sbatch（CK3，/smart-sbatch）
- 经 ssh baobab；脚本内 conda activate coding-rna。复用 run_TA-FOUNDATION-DECODER-M4.sbatch 模式（加 --postproc 透传）或新建 run_TA-COHERENCE-FIX-M5.sbatch <cand> <seed>。
- /smart-sbatch Phase 1 guard 必跑。shared-gpu（private 拥挤）。10 候选 run + 2 锚 run，注意并发，分批可。

## 评估 + 统计（CK4）
- 每 run：predict(+postproc)→GFF→eval --span-mode cds→aggregate(bw+macro)。
- 每候选 5-seed mean + **CI**（±std / bootstrap）on spec(bw+macro)、gbF1、**gene_count_ratio**；per-species(yeast)。
- **paired test vs 5-seed 锚** on intergenic_specificity（per-seed 配对，Wilcoxon/paired-t；注明 n）。
- **fragmentation 诊断**（本轮 HARD）：gene_count_ratio、predicted gene-length 分布、exon/CDS-run count 分布、transcript-span precision/recall（不只 base-level）。确认 constrained 把 2.25 拉下来。
- validate_goal --profile screen 每候选。

## Pivot decision menu（track-A-screen，CK5）
- **promote-ready**：CONSTR(或 CRFSTAB) 保 spec 严格>5-seed 锚 AND gbF1>=0.5276 AND macro>=0.7978 AND gene_count_ratio→≤1.25(或显著趋近) AND paired test 显著 → 赢家 promote-ready（③ Track-B，需用户 go-ahead；不在本 goal 启动）。
- **iterate**：spec 保住但 gene_count 仍>1.25（post-proc 不够）→ 调 constrained 参数 / 试 CRFSTAB / 组合；或 spec 因 post-proc 掉了 → 权衡。
- **fix**：post-proc bug / 锚补 seed 失败 → 修。
- 多选项点按短 prompt 决策自治（列选项→3 CLI reviewer→共识+ROI→docs/08→不暂停；破坏性/abandon/>24h/tied 例外）。**Track-B promote 启动需 go-ahead（>24h 新 sub-iteration 例外）**。

## Skill invocation chain
| step | skill | 状态 |
|---|---|---|
| impl | /implement（check_data 同 split 继承 + sanity smoke + 有界 debug） | 用 |
| 锚 seed | run_screen_ref.sbatch tiberius_like 3/4 | 用 |
| submit | /smart-sbatch | 用 |
| 对账 | scripts/job_watch.sh | 用 |
| 评估 | eval_gene_body_mask --span-mode cds + aggregate + 5-seed-anchor paired + fragmentation 诊断 | 用 |
| 记录 | /result-log | 用 |
| 复核 | /tri-review（3 CLI） | 用 |
| 决策 | /pivot | 用 |
| 实验档 | /exp-log + build_atlas | 用 |
| /revise-goal | — | **跳过（①已完成，0.9303/0.8710 新尺子已实证）** |
| /retrospective | — | advisory triggered（~6 iter，非阻塞）|

## Constraints (full)
- NON-CLAIM screen，不 scale-up（promote=③ 需用户 go-ahead，>24h 新 sub-iteration 例外）。
- 复用 FEATCACHE（不重抽取）+ 与锚同 chromosome split + 同冻结协议。**不重跑 /revise-goal**（①已完成；FPLOSS metrics.json 已证 intergenic 用 full-transcript 补集）。
- ≥5 seed/候选 + CI + paired test vs 5-seed 锚（现锚 3 seed、一个崩 0.773 拉低均值，paired 才能定 +0.059 是否真）。
- 门：spec 严格>5-seed锚 AND gbF1>=0.5276 AND macro>=0.7978 AND gene_count_ratio 趋向≤1.25（本轮 HARD 诊断，base-weighted spec 对碎片不敏感，必须独立看 gene-level）。
- FP-loss λ=1.0 硬编码、未在 test 上调（文档化，回应 reviewer "teaching to the metric" 关切）。constrained_decode 是确定性后处理（非学习型 CRF 的方差/崩塌）。
- env 只 coding-rna 绝不 base；sbatch 经 ssh baobab，脚本显式 conda activate。
- 预 claim 闸（未来 full/scale）：核实测试 clade 不在 SegmentNT 预训练语料。
