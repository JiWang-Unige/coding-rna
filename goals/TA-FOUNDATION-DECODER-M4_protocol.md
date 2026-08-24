# TA-FOUNDATION-DECODER-M4 · Protocol

> 执行手册。短 prompt 是 evaluator 唯一可见契约；本文件是操作细节。主架构赌注：foundation 特征 → 结构化解码器，把 recall 转成 intergenic specificity，严格超越新尺子锚。

## Permissions
- 全工具可用。写 scope：`FP-SEGNT-{FPLOSS,FUSION,CRF}` 命名空间的 configs/sbatch/outputs + `src/foundation_probe/` 内的新代码（FP-aware loss / fusion input / CRF wiring）。可改 `src/foundation_probe/train_probe_head.py`（本模块自有，非公共 screen_anchor 码）。
- 不动 ACTIVE_GOAL.json（本轮非 revise）。CLAUDE.md/docs03/09 只 draft patch。装依赖只进 coding-rna，绝不 base（jax+torch 已装齐，本轮无需新依赖）。

## Final goal（北极星对齐）
终极：在 DUAL co-primary 上严格超越 frozen published SOTA 且可发表。本轮 = 把 FP-SEGMENTNT-PROBE-M1 的发现（foundation 特征 gbF1 0.689>>锚 0.5576，但 intergenic_specificity 0.842<锚 0.871，yeast 溢出）转化为**击败 AXIS-1** 的架构。
- 新尺子 bracket（intergenic_specificity，test 子集）：FLOOR(ORF) 0.8805（gbF1 0.3735 被 F1 地板挡）< **screen_anchor 0.8710 bw / 0.8278 macro**（PROVISIONAL）< **pretrained_ceiling 0.9917**（Helixer 全数据）。anchor gbF1=0.5576。
- 双门（screen 方向选择，NON-CLAIM）：候选 promotable iff intergenic_specificity(bw) 严格>0.8710 AND gene_body_F1_unconstrained>=0.5276 AND macro_intergenic_specificity>=0.7978。screen 永不 claim。

## Track + resource
- Track A screen，3 候选并行 batch，各 5 seed = 15 run。profile=screen。
- 冻结协议（与锚/probe 同预算 same-budget）：window 2048、sample_fraction 0.3、epochs 8、patience 3、class-weighted CE(sqrt_inv)、anchor-matched conv+biLSTM head。yeast+fly，与锚同 chromosome split（复用 FEATCACHE 内 split）。
- **复用 `outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species/{yeast,fly}.npz`（per-seqid (L,14) fp16），不重抽取。**
- 资源：private-teodoro-gpu RTX3090 24GB。head 训练在 cached 特征上极快（~数分钟/seed）。run-and-evaluate。
- per-run 输出：`outputs/<cand>-s<seed>/{predictions,eval_subsets,metrics}`，<cand>∈{FP-SEGNT-FPLOSS,FP-SEGNT-FUSION,FP-SEGNT-CRF}。

## Parallel batch matrix
| exp_id | major_axis | mechanism_delta | seeds | config |
|---|---|---|---|---|
| FP-SEGNT-FPLOSS-s{0..4} | loss_design | 不对称 intergenic-FP / specificity-targeted loss（focal/precision-biased/FP-weighted）加到 CW-CE，惩罚落在真 intergenic(full-transcript 补集) 的预测-genic 碱基 | 5 | --head convlstm --loss fp_aware |
| FP-SEGNT-FUSION-s{0..4} | data_view | raw-DNA one-hot(5ch) ⊕ frozen SegmentNT logits(14ch) gated-concat → 同 conv+biLSTM head | 5 | --fuse-raw-dna |
| FP-SEGNT-CRF-s{0..4} | decoder | LinearChainCRFVec(reuse src/screen_anchor/decoders.py) on foundation emissions + FP-aware aux CE（非裸） | 5 | --decoder crf --loss fp_aware_aux |

提交：Option B 多 sbatch（3 候选脚本或一个带 --cand 参数的脚本 × 5 seed），各独立 output_dir。

## Orthogonality declaration
**verdict = PASS**（diverse batch，非 focused）。3 候选 major_axis 互不相同（loss_design / data_view / decoder），mechanism_delta 各异，均 structural（新 loss term / 新输入张量+fusion gate / 新 transition 参数+Viterbi）。共同点仅"输入用 foundation 特征"（继承自 probe，非本轮变量）+ 共同目标"recall→specificity"。无 HARD_FAIL（非全调 lr/batch/dropout；mechanism 不重复；各有 why-structural）。

## 实现要点（CK1，/implement）
- **train_probe_head.py 加三个开关**（默认 off，保持 probe 行为不变）：
  - `--loss {ce, fp_aware}`：fp_aware = CW-CE + λ·intergenic-FP penalty。FP penalty = 对预测为 genic(class>0) 但真值 intergenic(class 0) 的碱基加权（或 focal on class-0 precision）。λ 暴露为 arg（默认中等，可 sanity 调）。
  - `--fuse-raw-dna`：在 (B,14,W) 特征前 concat one-hot(5ch)→(B,19,W) 进 head（或 gated：learnable gate 融合）。需在 dataset 同时取 cached 特征 + 即时 one-hot（复用 D.one_hot_window，对齐窗口；窗口 seq 从 genome 取）。
  - `--decoder crf`：在 head 输出 emissions 上接 LinearChainCRFVec（reuse decoders.py）；训练用 nll + FP-aware aux CE；predict 用 batched viterbi（reuse train_screen_ref 的 decode_batch 思路）。
- 自审（不可跳）：label 对齐、loss 与目标一致、metric 口径同锚、无泄漏。CRITICAL 未修不进 smoke。
- **每候选 1-seed sanity smoke**（srun GPU，epochs 2-3）：出 metrics、per_class 非单类（不塌缩）、FEATCACHE 对齐 guard 过。崩/塌缩 → repair_advisor 有界修复 ≤3 次或修该候选；不盲投全 batch。

## Pre-submit gate（CK1，HARD）
- 3 候选各 sanity smoke 通过（不塌缩、出 spec+F1）。
- FEATCACHE 完整（trainer guard feat L==genome L）。
- 数据契约：fusion 候选的 one-hot 与 cached 特征**同窗口对齐**（同 seqid/s/e）；CRF 候选 emissions→viterbi shape 对。任一 ❌ 不提交。

## sbatch（CK2，/smart-sbatch）
- 经 ssh baobab 提交；脚本内 `source .../conda.sh && conda activate coding-rna`。
- /smart-sbatch Phase 1 guard 必跑（VRAM≥20、output_dir 唯一、time、exclude 3080、并发≤配额）。15 run 注意并发上限，必要时分批。
- 提交后 job_watch.sh 对账（COMPLETED/FAILED/TIMEOUT/OOM/STALE，不假设成功）；幽灵 run 按 failed_run。

## 评估（CK3-CK4）
- 每 run：predict→CDS GFF(gff_io)→`eval_gene_body_mask.py --span-mode cds`（新尺子）per 物种→`aggregate_gene_body_metrics.py`（bw+macro）。
- 每候选：5 seed seed-mean + **CI（±std 或 bootstrap）** on intergenic_specificity(bw+macro)、gene_body_F1、gene_count_ratio；**per-species**（yeast 单列）。
- **paired test vs 锚**：用锚的 5 seed（若只有 3 seed，注明）做 per-seed paired 比较 intergenic_specificity。
- `validate_goal --profile screen` 每候选（screen→not_yet/progress，claim_gate 报是否>锚）。

## Pivot decision menu（track-A-screen，CK5）
- **promote-to-Track-B**：某候选 intergenic_specificity(bw) 严格>0.8710 AND gbF1>=0.5276 AND macro>=0.7978 AND paired test 显著 → 该机制把 recall 转成了 specificity → Track-B scale（≥8 seed + 数据/epoch 扩 + 仍 NON-claim 直到 full）。
- **iterate**：接近未过（如 spec 升但未超锚、或 yeast 仍拖 macro）→ 调 λ / fusion gate / 组合两机制（FPLOSS+FUSION）再来一轮。
- **abandon-candidate**：某候选明确无效（spec 不升或更差）→ 记 docs/06，留其余候选。
- **abandon-route（写 docs/09 需人闸）**：三候选全部无法把 foundation recall 转成 specificity → foundation→decoder 路线受质疑，转 GENERanno probe 或 from-scratch 加 FP-aware。
- 多选项点按短 prompt 决策自治：列选项→3 CLI reviewer→共识+ROI→写 docs/08→不暂停（破坏性/abandon/>24h/tied 例外才停）。

## Skill invocation chain
| step | skill | 状态 |
|---|---|---|
| impl | /implement（check_data 同 split 继承 + sanity smoke + 有界 debug） | 用 |
| submit | /smart-sbatch（Phase1 guard + Phase2；15 run 并发分批） | 用 |
| 对账 | scripts/job_watch.sh | 用 |
| 评估 | eval_gene_body_mask --span-mode cds + aggregate + validate_goal + 跨 seed CI/paired | 用 |
| 记录 | /result-log | 用 |
| 复核 | /tri-review（3 CLI） | 用 |
| 决策 | /pivot | 用 |
| 实验档 | /exp-log + build_atlas | 用 |
| /reproduce-baselines | — | 跳过（已完成） |
| /retrospective | — | advisory triggered（~5 iter）；非阻塞，可在本轮前后跑 |

## Constraints (full)
- NON-CLAIM screen；不 scale-up（promote 留 Track-B 新 goal）。
- 复用 FEATCACHE（不重抽取）+ 与锚同 chromosome split + 同冻结协议（same-budget，apples-to-apples）。
- ≥5 seed + CI + paired test vs 锚（AXIS-1 方差脆弱：M1 probe spec per-seed 0.808-0.897，单 seed 已超锚）。
- 双 co-primary 门：intergenic_specificity(bw) 严格>0.8710 AND gene_body_F1>=0.5276 AND macro>=0.7978。不达→not_yet（非 failed_run）。
- 3-class 输出（复用 harness 干净对锚）；richer 多类(CDS/intron/intergenic/phase/splice/strand)=紧接 follow-up，**非本批**（需 label-schema+GFF-output 改）。
- CRF 候选**带 FP-aware aux 训练**，非裸 semi-CRF（M1：裸结构化 decoder 在 FP-heavy emissions 上伤 specificity）。
- per-species + gene_count_ratio 本轮作 **HARD 诊断**（yeast 过预测是主问题）。
- 环境纪律：只 coding-rna 绝不 base；sbatch 经 ssh baobab，脚本显式 conda activate。
- 预 claim 闸（未来 full/scale）：核实测试 clade 不在 SegmentNT 预训练语料（它见过脊椎动物）；pilot yeast+fly 对 screen 干净。
