# 31 P6 面板基线：设定核对、评分协议、评分器验证、暴露审计、运行记录与 C1 判定

来由：
- 大纲 `docs/paper_outline_20261027.md` 的 P6（面板基线），用于 Fig. 4 / Table 1。
- C1 面板判据见报告 28 §6（11-04 写定）；补跑 seed 的规则见报告 28 §7（11-07 写定）。
- 用户 11-07 指示：P6 不等 Baobab（12-01 credit 重置），ont 与 win3090 只要没有 TE 项目的任务在跑就可以用。

代码在 `experiments/PANEL-P6/`（ont 副本 `~/coding-rna-reorg/experiments/PANEL-P6`，输出 `~/coding-rna/outputs/PANEL-P6`）。

**本文件 §1–§4 写于 2026-11-07，早于任何面板分数**：此时我方 T2 的面板推理（P5）已完成，但没有任何工具在面板物种上被打分。§2 的评分协议之后不再修改；§5 起是运行记录与结果。

## 0. 状态

| 项 | 状态（2026-11-17 更新） |
|---|---|
| 设定核对 | 完成（§1）；发现两处与旧文档不符，均按固定版本处理（§1.2） |
| 输入 | 10 个整核基因组 FASTA 已在 ont 建好，计分主序列与 T2 读到的逐字相同（§1.3） |
| 评分器 | 写好并验证：59/59 项与已发表的 DEV / set B 值逐位一致（§3） |
| 暴露审计 | 完成（§4，明细 `experiments/PANEL-P6/results/exposure_audit_p6.md`） |
| GPU 运行 | 完成：11-09 15:45 至 11-15 09:22，全部在 ont；按"TE 优先"让行多次（§5.2、§5.3） |
| 面板打分与 C1 判定 | 11-17 07:50 打分一次。**C1 不成立**：T2 面板中位数 55.05，比 Carbon-A 过滤版低 12.60 点，且低于 Tiberius、Vipsania ft。补跑 seed 规则触发（§6） |

## 1. 设定核对

### 1.1 工具、版本、参数

除输入外，全部沿用 DEV 原生基线（M28，`benchmarks/native_baselines/reports/native_baseline_preparation.md`）、M29/M30 Carbon-A、M38 Vipsania 的同一资产与命令。资产从 Baobab 原样 rsync 到 ont `outputs/PANEL-P6/assets/`，没有重新拉取或重建镜像。

| 工具 | 版本与资产 | 模型 | 参数 | 判据中的角色 |
|---|---|---|---|---|
| Helixer | 0.3.7（源码 d17bb496，HelixerPost 4d4799b），SIF = DEV 同一文件 | `land_plant_v0.3_a_0080.h5`（官方优先级 80） | `--subsequence-length 64152 --batch-size 16 --deterministic --no-multiprocess` | 判据：T2 中位数须 ≥ 它 |
| Tiberius | 2.0.8 源码（871d497）隔离安装在 2.0.7 底座 SIF 中（DEV 同一组合；不是发布方的 2.0.8 镜像） | `angiosperms`（Mesangiospermae） | `--batch_size 4 --seq_len 200070 --no_softmasking`，无 ClaMSA | 判据：T2 中位数须 ≥ 它 |
| ANNEVO | v2.3.3（42c920f） | `ANNEVO_Magnoliopsida.pt`（blob 21487d87） | `annotation.py -l Magnoliopsida --batch_size 32 -t 8 --num_workers 8`，其余为默认（min_intron 20、min_prot 100） | 描述 |
| Carbon-A | HF `HuggingFaceBio/Carbon-A-1.2B` f6789ee，项目管线副本 `carbon_a_pipeline_copy` | 单一检查点；未过滤，与置信 ≥ 0.1 或 NA 的过滤版 | `model_inference.py --batch_size 4 --dispatch_chunk_size 1 --attn_implementation sdpa`；CPU `post_processing.py --viterbi_backend numba` | 判据：过滤版中位数 − T2 中位数 ≤ θ = 6.0 |
| Vipsania | 1.1.0，环境 `~/envs/vipsania_m38`（M38 同一环境） | Streptophyta `j9m0cdmk`（按 ID 固定，sha256 校验） | ft：`--finetune`，作者默认（10 epoch、lr 1e-4、有效批 64）；pre：只用预训练 | ft 进判据（T2 中位数须 ≥ 它）；pre 作描述 |
| 我方 T2 | P5 面板推理（报告 20 §14），T2 `model_002304.pt` | — | 不做后验过滤 | 被判定方 |

运行环境上的必要差异（都不改变预测）：
- ANNEVO 在 Baobab 上用的是 conda 环境 `annevo`，官方没有容器。ont 上新建 uv 环境 `~/envs/annevo_p6`，逐个包与 Baobab 环境同版本：Python 3.10（3.10.19 对 3.10.20）、torch 2.1.0+cu121、numpy 1.26.4、h5py 3.14.0、biopython 1.87、pandas 2.3.1、numba 0.65.1、tqdm 4.67.1、torchmetrics 0.8.2。
- Tiberius / Helixer 的容器挂载点改名为 `/PANEL-P6-*`，让 GPU 闸门能认出我们自己的容器进程（§5.1）。
- 两个 TensorFlow 工具加 `TF_FORCE_GPU_ALLOW_GROWTH=true`，不一启动就占满整张卡。
- 遇到显存不足时减半 batch，实际值记入每次运行的 `run.json`。batch 只影响速度。
- Vipsania 的 `--finetune_B` 规则同 M38：空闲 ≥ 23,040 MiB 时用 2，否则用 1，靠梯度累积保持有效批 64。新增一条：先最多等 2 h 争取 B = 2。
- Dryas、Eschscholzia 的 Vipsania pre 直接复制 M38 的注释，不重跑：M38 用的是同一个整核基因组（sha256 相同）、同一模型、同一环境（`run.json` 记 `reused_from`）。

### 1.2 与旧文档不符之处

1. **ANNEVO 参数**。
   - 报告 10 §4.2 写的"推荐用法 `--num_classes 15 --boundary-aware --min_intron_length 20`"在固定版本 42c920f 上不存在。
   - `annotation.py` 只有 `--min_intron_length`（默认 20）与 `--min_prot_length`（默认 100）；子代理查到这两个开关自 v2.3.0 起已删除，15 类边界感知写死在模型里。
   - 所以 DEV 时不带参数的命令就是该版本的官方用法（报告 deck 核对表 D20 的疑问到此解决），面板沿用。
2. **Carbon-A 现成注释（CADB）**。
   - 大纲 §5 第 3 条与报告 25 §9 曾建议先查 CADB。
   - CADB 只覆盖核心 10 种中的 6 种（缺 Eschscholzia、Hakea、Zizania、Themeda），且是否带置信分未核实。
   - 混用两种来源会让 Carbon-A 的物种分数口径不一，因此 10 种统一自己跑，与 DEV / set B 的 Carbon-A 同一管线（SDPA）。
   - SDPA 相对作者默认的 FlashAttention-2 在一条 TRAIN 染色体上约低 0.2 点（deck D16），面板上同样存在，写进局限。

GeneCAD 不在本次 P6 内：按报告 25 §9 与大纲 P6，它在 12-01 后在 Baobab 上跑，提交前先问用户。

### 1.3 输入：整核基因组 FASTA

脚本 `p6_genomes.py`，规则同 M38 `vip_make_genomes.py`：
- 来源：PANEL-PREP 冻结的 NCBI 下载（报告 25 §2，md5 已核对）。
- 去掉装配报告中 non-nuclear 单元（叶绿体、线粒体），其余全部序列保留。
- 两个版本：
  - `<key>.nuclear.fa`：保留 NCBI 软屏蔽，给 Vipsania（M38 同样做法）；
  - `<key>.nuclear.upper.fa`：转大写，给 Helixer、Tiberius、ANNEVO、Carbon-A（DEV 原生基线与 M29/M30 的 Carbon-A 都用大写输入）。
- 核对：每个物种的计分主序列（`prepared/<key>/lengths.json`）在整核基因组中长度相同，且大写后的序列 md5 与 T2 读到的 `prepared/<key>/genome.fa` 相同。
- Dryas、Eschscholzia 的软屏蔽版与 M38 用的整核基因组 sha256 相同。

| 物种 | 序列数 | 总长 Mb | 计分主序列（条 / Mb） | 其余核序列 Mb | 软屏蔽比例 |
|---|---:|---:|---|---:|---:|
| Dryas octopetala | 25 | 229.6 | 9 / 229.3 | 0.3 | 0.28 |
| Juncus effusus | 383 | 258.3 | 21 / 239.9 | 18.3 | 0.44 |
| Acorus gramineus | 137 | 398.6 | 11 / 289.1 | 109.5 | 0.42 |
| Eschscholzia californica | 329 | 400.7 | 6 / 379.2 | 21.4 | 0.34 |
| Rhynchospora breviuscula | 1,505 | 420.7 | 5 / 370.1 | 50.6 | 0.38 |
| Zizania latifolia | 289 | 595.8 | 17 / 571.7 | 24.1 | 0.44 |
| Akebia trifoliata | 913 | 652.8 | 16 / 605.6 | 47.1 | 0.42 |
| Amborella trichopoda | 25 | 700.3 | 13 / 699.5 | 0.8 | 0.39 |
| Hakea prostrata | 334 | 712.4 | 11 / 662.8 | 49.6 | 0.40 |
| Themeda triandra | 22 | 754.7 | 10 / 746.0 | 8.7 | 0.39 |
| **合计** | 3,962 | **5,123.9** | **119 / 4,793.3** | 330.6 | |

清单 `experiments/PANEL-P6/results/genomes_manifest.json`。

各工具都只在主序列上计分。只有 Vipsania 的微调会读整个基因组：它的微调抽样遍及全部核序列。其余工具逐条序列预测，主序列上的结果与是否输入其余序列无关。

## 2. 评分协议（2026-11-07 写定，早于任何面板分数）

### 2.1 范围与参考

- **物种**：核心 10 种（报告 25 §9）。
- **序列**：119 条主序列（与 T2 推理范围相同）。
- **参考**：冻结 ruler 的 M26 CDS-assessable 主参考（`common_ruler.select_reference`），由 `data/panel/prepared/<key>/reference.gff3` 构建。代码与 M29 `build_reference` 相同，只换路径。

### 2.2 指标

**主指标**：每个物种的 exact-chain F1。
- 每个预测基因取其主转录本（`primary_transcripts`，最长 CDS）的完整 CDS 链与链方向，唯一化后与参考精确匹配。
- **对所有工具和 T2 一律应用 CDS ≥ 96 bp 规则**，即丢弃总 CDS < 96 bp 的预测链。这与 M34 / P2b 我方各臂、M38 Vipsania 的主指标同口径。
- 不加此规则的值作为次要列一并报告。
- DEV 上此规则的影响已实测（§3）：ANNEVO +0.04，Tiberius 0，Helixer +0.25，Carbon-A 0（其管线自带 ≥ 96 bp 下限）。

**参评的 run**：
- T2：P5，未过滤；
- Helixer、Tiberius、ANNEVO；
- Carbon-A 的未过滤版与过滤版；
- Vipsania 的 ft（判据模式）与 pre（描述）。

### 2.3 面板统计量

- 每个 run 取 10 个物种 F1 的中位数（C1 用）、均值、最小值。
- **CI**：P2b 的 1-Mb 块 bootstrap，推广到多物种。
  - 块 =（序列，首个 CDS 起点 // 1 Mb），取参考和全部参评 run 的块之并集。
  - 2,000 次联合重抽，`default_rng(0)`。
  - 每次重抽得到每个物种的 F1、面板中位数、面板均值，以及 T2 与各工具在这些量上的配对差。
  - 只有 2 个物种时，这一代码路径与 P2b 的结果逐位相同（§3）。
- **物种层面配对检验**（报告 10 §4.1 预注册）：对 T2 − 工具的 10 个物种差值做 Wilcoxon 符号秩（双侧）和符号检验，对 6 个工具做 Holm 校正。

### 2.4 C1 判定（报告 28 §6；档次文字见大纲 §0）

记 med(X) 为 X 在 10 个物种上的 F1 中位数。

- **成立**：同时满足：
  - med(T2) ≥ med(Helixer)、med(Tiberius)、med(Vipsania ft)；
  - med(Carbon-A 过滤) − med(T2) ≤ 6.0。
- **不成立**：med(Carbon-A 过滤) − med(T2) > 10。
- **部分成立**：其余情况，即差距在 6 与 10 之间，或低于 Helixer、Tiberius、Vipsania ft 中任一个。
- 大纲"不成立"还有一条"面板上 T2 相对 NT-v2 + LoRA 的优势消失"。NT-v2 + LoRA 没有在面板上跑，这一条不评估，在结果里注明。

判定只用上述中位数的点估计。bootstrap CI 与物种检验作为描述一并报告。

### 2.5 补跑 seed 规则（报告 28 §7）

P6 打分后以下任一成立即触发：
- T2 对 Helixer、Tiberius、Vipsania ft 三者中位数的最小领先 < 2.5（含落后）；
- med(Carbon-A 过滤) − med(T2) 落在 [3.5, 8.5]。

触发后按 §7 在 ont 上对 T2_s2、T2_s3 各跑一次面板推理，同路径、同评分。C1 仍按 T2 判定。

### 2.6 盲态与一次打分

- 在全部基线对 10 个物种都完成之前，不在任何面板物种上计算任何分数。
- `p6_score.py panel` 在任一工具缺任一物种时拒绝运行；结果文件已存在时拒绝覆盖。
- 运行期间只看运行记录（基因数、GFF 合法性、耗时），不读参考注释。
- 评分器验证（§3）只用 DEV 与 set B 的序列。这两处的数字都已发表或已被看过。

### 2.7 敏感性与描述（不改判定）

- **Tiberius 的适用范围**：angiosperms 模型适用于 Mesangiospermae，Amborella 在其范围外（报告 10 §1.3）。
  - 判据按"核心 10 种"原文，Tiberius 的中位数用全部 10 种。
  - 另报去掉 Amborella 后 9 种的中位数，以及用它代入时的判定。
- **分层中位数**：TH / SH / AH 层、禾本科 / 非禾本科、RefSeq / GenBank 注释来源。
- **外显子数分层**：1 / 2–5 / ≥ 6 段。
- **set B 上的 Vipsania ft**：报告 30 §7 用户决定 B 集的 ft 分数取自面板运行。面板打分时，同时用 M38 `vip_score.py setb` 对面板 ft 的 Dryas、Eschscholzia GFF 在 set B 区段上评一次，补全报告 30 §8。

## 3. 评分器验证（2026-11-07，ont，`p6_score.py validate`，59/59 通过）

方法：用同一转换器（M38 `vip_score.convert`）和同一计分路径（`m30_score.Evaluator` + `load_gff_chains`，即 P2b），复算已发表的数字，要求逐位一致。结果文件 `experiments/PANEL-P6/results/validation.json`。

| 集合 | 对象 | 复算（tp / pred / ref，F1） | 已发表 |
|---|---|---|---|
| DEV（规则关，原口径） | ANNEVO | 5,644 / 7,183 / 7,728，75.7025 | 同（M28；两物种分项也相同） |
| | Tiberius | 5,441 / 7,427 / 7,728，71.8047 | 同 |
| | Helixer | 5,112 / 8,539 / 7,728，62.8512 | 同 |
| | Carbon-A 未过滤 / 过滤 | 5,730 / 7,503；5,687 / 6,966 | 同（M29 G0） |
| DEV（规则开） | T2（P5 用的模型） | 73.878 | 同（报告 20） |
| | Vipsania ft | 69.473 | 同（M38） |
| set B（面板参考） | 7 条序列的参考键集合 | 与 M30 参考（`data/m1_screen`）逐条相同 | — |
| | T2 / T10（池化与两物种） | 14,947 / 21,785 / 21,875（68.470）；15,058 / 22,200（68.329） | 同（P2b） |
| | 块 bootstrap：块数、各 run 的 CI、T10 − T2 的 CI | 276 块，CI 全部逐位相同 | 同（P2b） |
| | Carbon-A 未过滤 / 过滤（两物种） | 6,885 / 9,532、8,683 / 11,838；6,807 / 8,895、8,584 / 11,048 | 同（M30 / M33） |
| | Vipsania pre（池化与两物种） | 与 M38 表相同 | 同（M38） |
| 判定逻辑 | 5 组合成中位数 | 成立 / 部分成立 / 不成立与触发均按 §2.4–2.5 | — |

由此可确认：
- 转换器对 Helixer、Tiberius、ANNEVO 的 GFF 方言同样无损（M38 已验证 Carbon-A、GENERanno、Vipsania）。
- P2b 计分路径与 DEV 计分路径对外部工具给出相同计数。
- 面板参考的构建方式在 set B 上与 M30 完全一致。
- 加 CDS ≥ 96 规则后的 DEV 值：ANNEVO 75.743（pred 7,183 → 7,175）、Tiberius 71.805（不变）、Helixer 63.102（tp 5,112 → 5,110，pred 8,539 → 8,468）。

## 4. 暴露审计（补核报告 25 §5–6 的开放项）

明细与出处在 `experiments/PANEL-P6/results/exposure_audit_p6.md` 与同名 `.json`（子代理只读调研，2026-11-07）。

与 `PANEL-PREP/scripts/tool_exposure.json` 相比的变化：
- **Tiberius**：按预印本 v3 补充表 S3 / S32 与 `angiosperms.yaml` 逐格复核。
  - Eschscholzia 应记 **T**（Tiberius 测试种，训练近缘为目级），其余一致。
  - Acorus americanus 确在训练集，所以 Acorus gramineus 记 **G**。
  - Figshare 训练数据是一个 16.4 GB 的 tgz，没有下载，内部清单无法核实。
  - YAML 的训练表比 S3 少 Aquilegia coerulea，原因无法核实；这只影响 Eschscholzia、Akebia 的目级近缘。
- **ANNEVO 现行 Magnoliopsida 权重**：训练物种未公开，10 格全部记"无法核实"。论文版 Table 24 的码不适用于现行权重。
- **Vipsania j9m0cdmk**：
  - 预训练序列含面板的 Akebia、Hakea、Juncus、Amborella（同一装配或其配对 GCA），记 S；Acorus 记 G。
  - Juncus 同时是 Vipsania 的测试种。
  - 它不用注释，所以这是序列暴露，不是标签泄漏（报告 28 §1）。
- **GeneCAD**（12-01 后另跑）：现行默认模型在 8 个物种上微调，对面板只到同科或同目。骨干语料按 PlantCAD2 的 65 属推断。
- **Helixer a0080 与 Carbon-A**：未发现与现有码冲突的事实，沿用。

论文面板表按工具标注暴露码。Vipsania 的 S 与 Carbon-A 在 Amborella 上的 S 单独注明。

## 5. 机器路由与运行记录

### 5.1 路由（2026-11-07）

**win3090**：
- 11-07 07:27（北京时间）检查时，3090 被 TE 项目的 Vipsania 微调占满（`ab-initio-TE/scripts/experiments/PHASE1-VIPSANIA-20261103`，鸡 noMask 臂）。
- 每个 epoch 约 3.45 h，每臂约 35 h；`ft.tsv` 共 6 臂，按此估计约 11-15 前后才空。
- 按用户规则不抢占，本次不用 win3090。

**ont**：5 个工具全部放在 ont 的 4090 上。
- 此刻 GPU 上有我方 M34 k25U（续跑到 1,536 步）、TE 的 `PHASE1-WINOVERLAP` 推理，以及 ncRNA 的训练（11-07 00:27 检查时）。
- 插入时段需与维护会话商量（已发消息，11-07）。
- 估计 GPU 量（4090，按 DEV / set B / M38 实测外推到 5.12 Gb）：

  | 工具 | 估计 GPU 时 |
  |---|---|
  | ANNEVO | 约 3 h |
  | Tiberius | 约 5 h |
  | Helixer | 约 10 h |
  | Carbon-A | 约 27 h（另有 CPU 后处理） |
  | Vipsania | 约 38 h（10 个基因组各一次微调，每次约 3.5 h；B = 1 时约翻倍） |
  | **合计** | **约 80–90 h** |

**GPU 闸门**（`p6_lib.sh`，每个工具的每个基因组开始前检查），四个条件都满足才开始：
1. 存在 `$P6/GO`（与维护会话定好时段后手动写入）；
2. 没有 TE 项目的 GPU 进程。判定条件：
   - cwd 或参数含 `/work/te`、`ab-initio-TE`、`TE-GLM`、`te_benchmark`；
   - 或者是不属于我们的 Apptainer / Singularity 容器进程。

   这吸取了 ncRNA 会话 11-03 漏判 TE 容器化 Tiberius 的教训。11-07 实测：闸门正确把 TE 的 `winov_infer.py` 判为 TE，把 k25U 判为 coding-rna。
3. 没有其他 coding-rna GPU 作业（M34 / M39），除非维护会话同意同卡并行；
4. 空闲显存与磁盘足够。

条件不满足就等，每小时记一条 `gate_wait` 事件。

**11-07 更新**：
- TE 会话通知：TE 从 11-07 起占 ont GPU 约 2.5 天，用 `/home/ji/ab-initio-TE/outputs/TE_GPU_CLAIM` 声明，结束时删除。闸门已加一条：此文件存在即等待。
- 闸门只检查 nvidia-smi 的计算进程，TE 只用 CPU 的 Apptainer 作业（BRAKER）不会挡住它。
- 闸门还会等 M34/M39 的队列脚本全部退出，所以 P6 不会插进 M39 两臂之间的空档。
- M39 的续接脚本只看空闲显存，k25U 一结束就开跑。
- 维护会话（`local_f7d86278`，ont 队列与看板的维护者）回复：
  - P6 与 M39 的先后由用户定，看板决定卡 `ont-order-after-k25u`（推荐 P6 → P4 补 seed（若做）→ M39）；
  - M39 的续接脚本已加"无 TE GPU 进程"检查（提交 773665c）；
  - 在用户定之前不要自己往 ont 上排 GPU 作业，只做不占 ont GPU 的准备；win3090 一空出来可以先用。
- 据此 ont 上的 `GO` 已撤回（改名 `GO.withdrawn-20261107`）。队列仍在 ont 上，停在闸门处：没有 GO 就不会使用 GPU。

**win3090 准备**（不占 GPU）：
- 范围：win3090 只跑 ANNEVO、Tiberius、Helixer。
  - Vipsania 在 TE 的 win3090 运行里实测约 4 s/步，ont 4090 上是 0.185 s/步，所以 Vipsania 留在 ont。
  - Carbon-A 依赖 ont 上的 generanno 环境与权重，也留在 ont。
- 部署：
  - conda-forge 装 apptainer 1.5.4（`~/envs/apptainer_p6`，与 ont 同版本）；
  - ANNEVO 环境与 ont 同一组包版本；
  - SIF、权重、ANNEVO 代码从 Baobab 经 Mac 流式传入，实测约 1.4 MB/s，约 3 h，逐字节相同；
  - 大写基因组从 ont 流式传入，按清单 sha256 核对。
- 闸门：在 win3090 上，任何以 `scripts/experiments/` 启动、工作目录在 ab-initio-TE 下的进程都算 TE 任务。TE 的 Vipsania 链在 GPU 臂之间有 CPU 准备阶段，只看 GPU 进程会插进空档。
- 结果：GFF 传回 ont 的 `native/<tool>/<key>/`。ont 队列遇到已有 DONE 的物种会跳过，两台机器不会重复跑同一工具、同一物种。
- 11-07 核对：
  - win3090 上的 2 个 SIF、Tiberius 与 Helixer 权重、ANNEVO 代码与权重、Tiberius 2.0.8 隔离安装目录（32 个文件），sha256 与 ont 副本全部相同；清单见 `experiments/PANEL-P6/results/assets.sha256`。
  - 10 个大写基因组的 sha256 与 ont 清单相同。
  - 闸门在 win3090 上正确识别 TE 的 GPU 进程与 `vipsania_finetune.sh` / `vipsania_pretrained.sh`。
- 11-07 10:47（北京时间）：win3090 写入 GO，启动 `p6_win3090.sh`（冒烟 → ANNEVO → Tiberius → Helixer），停在闸门处等 TE 结束。
- TE 在 ont 上的 GPU 声明延长到约 11-10 中午（TE 会话 11-07 通知）。
- 11-08 约 11:12–11:41（北京时间），TE 的家蚕微调把 WSL 内存撑爆，整个 WSL 被停，随后由 TE 会话拉起。
  - 我们的 P6 队列当时还在闸门处、没用过 GPU，没有损失。11:44 重新挂上。
  - 为防同类情况，win3090 上 P6 的资源守护内存上限设为 32 GB（ont 为 64 GB）。
  - TE 作业现在跑在 MemoryMax=40G 的 systemd scope 里，闸门照样能认出。
- 11-08 用户定优先级（由论文 draft 会话转告）：TE 在 win3090、ont、CAU 上都排第一，P6 排在 TE 之后。
  - win3090：TE 的 Vipsania 剩余臂暂停，先跑 G2（约 25 GPU-h）。
  - ont：GENERanno 扫描之后接 G1（约 10 GPU-h）。
  - P6 两个队列继续停在闸门处，TE 空出卡时会通知。

**CPU 工作**：建基因组、评分器验证、Carbon-A 后处理都绑核 25–31、nice 10。TE 的 GeneMark-ETP 在 ont 上用核 2–31，我们的 CPU 进程优先级较低。

### 5.2 运行记录

**排期**：
- 11-09 用户在决定卡 `ont-order-after-k25u` 上选了 P6 → P4 补 seed → M39。
- 维护会话 14:38 停掉 M39（I10_k50U 在第 320 步存点，之后续训），14:39 写入 ont 的 `GO`。
- 当时 ont GPU 上只有 ncRNA 的训练（约 6 GB），TE 在 ont 的 GPU 工作已于 10:52 结束。
- win3090 上的 TE 还要占很多天（G2 之后接 Vipsania 剩余臂），所以五个工具全部在 ont 跑。win3090 的 `GO` 已作废（`GO.withdrawn-20261109`）。

**开跑时的三处修复**：都不改变计算，只影响能否跑起来。
1. **14:41，ANNEVO 冒烟失败**：`--tmp_path` 目录不存在。DEV 的运行入口会预建这个目录，我的驱动漏了。改为每次运行前为三个原生工具预建 `temporary/`。
2. **14:43–15:37，闸门死锁**：GPU 空了约 54 分钟。
   - 维护会话的 `m34_after_p6.sh` 只是在等 P6 结束，却被闸门当成 coding-rna 队列阻塞者，而它又在等 `p6_queue.sh` 退出。
   - 改为闸门忽略名字含 `after_p6` 的脚本。
3. **15:43，ANNEVO-Dryas 被资源守护杀掉**：M38 守护限制的是进程组的 RSS 之和，算出 73.7 GB > 64 GB。
   - ANNEVO 的 8 个 fork 出的 DataLoader 进程共享父进程的写时复制页，RSS 相加会重复计算；当时主机仍有 119 GB 可用。
   - 新写 `p6_guard.py`：接口与 M38 相同，但限制 PSS 之和（真实占用），另记 RSS 之和。重跑时整机已用内存约 6 GB。

**冒烟**（15:38–15:42，2 Mb 拟南芥训练切片 NC_003070.9:5,000,001–7,000,000，不计分）：五个工具全部通过，`SMOKE_OK`。

| 工具 | 基因数 | 墙钟 | GPU 峰值 |
|---|---:|---:|---:|
| ANNEVO | 538 | 6.5 s | —（运行短于 15 s 采样间隔） |
| Tiberius（apptainer --nv） | 550 | 37.6 s | 4.6 GB |
| Helixer（apptainer --nv） | 601 | 20.6 s | 4.6 GB |
| Carbon-A | — | 推理 22.6 s + 后处理 20.5 s | 6.7 GB |
| Vipsania pre | 通过 | — | — |

**正式运行**：15:45:20 起，ANNEVO → Tiberius → Helixer → Carbon-A → Vipsania，每个工具按基因组从小到大。逐物种记录在每个运行目录的 `run.json`，汇总表见 §5.3（全部完成后填）。

| 工具 | 时段 | 墙钟 | 每物种 | 内存峰值（PSS） | GPU 峰值 |
|---|---|---:|---|---|---:|
| ANNEVO | 11-09 15:45–17:08 | 1.3 h | 3–11 min | 13–21 GB | 5.4 GB |
| Tiberius | 11-09 17:08–22:26 | 5.3 h | 14–42 min | 22–36 GB | 4.6 GB |
| Helixer | 11-09 22:26 – 11-10 06:12 | 7.8 h | 20–70 min | 2–3 GB | 4.6 GB |
| Carbon-A | 11-10 06:12 – 11-12 15:19（其中约 30 h 按 TE 优先在闸门等待） | 约 26 h（推理 + 后处理） | | 后处理 7–17 GB | 6.7 GB（推理） |
| Vipsania | 11-12 15:19 起 | | | | 12.9 GB（B = 1 微调） |

**11-10 Carbon-A 的间歇性失败**：

失败情况：

| 输入文件 | 时间 | 现象 |
|---|---|---|
| Dryas NC_141366.1 | 06:19 | 推理已完成，"合并结果"之后、写 parquet 之前收到 SIGSEGV（rc −11） |
| Juncus unplaced | 08:12 | tokenizer 在纯 Python 处报 `TypeError: 'int' object is not iterable`（`tokens_trie.split`） |
| Juncus CM050892.1 | 08:33 | 同第一行，SIGSEGV |

- 三处都是同一输入重试即成功，原因未查明（像间歇性的内存或进程间问题）。M30 跑 set B 时没有出现过。
- 每次崩溃都留下一个 multiprocessing spawn 子进程（父进程死后被 init 收养），各占 6.7 GB 显存。闸门要 13 GB 空闲，于是 08:33–10:12 一直等待，GPU 空着（由维护会话发现）。

**修复**：
- `p6_guard.py` 在命令结束后（无论成败）按进程组清掉残余进程，结果文件记 `orphans_killed`；
- Carbon-A 每个文件最多重试 5 次（原 3 次）。

**11-10 16:03 Carbon-A 队列中止**：
- 经过：Rhynchospora 的 unplaced 文件（1,500 条短 scaffold）后处理其实成功写出了 GFF。但 Carbon-A 管线的密码子表启发式（只有 34% 的完整基因达到置信 0.1）把输出目录改名为 `unplaced_WARNING`，并写入 `WARNINGS.txt`。驱动在原目录找不到 GFF，判为后处理失败，中止了队列。
- 影响：
  - 该后处理只覆盖 ≥ 100 kb 的序列（`min_record_length=100000`，7 Mb），这部分不计分，预测本身不受影响；
  - GPU 空了 16:03–16:32（我没有及时处理等待器的返回）。
- 修复：驱动遇到 `*_WARNING` 时把目录移回原名，警告记入 `post.json` 的 `codon_table_warning` 和事件 `post_codon_warning`。本次已手工照此处理。
- 报告时列出所有触发此警告的文件；若有主序列触发，单独说明。

**11-10 17:08–17:25 磁盘闸门暂停**：
- ncrna-annot 项目写入大量文件（维护会话查明，如 83.9 GB 的 `mouse_chr19.npy`），ont 根分区余量从 242 GB 降到 91 GB，P6 在 120 GB 的下限处暂停。
- P6 本身只占 27 GB，后续再写不到 10 GB，所以把下限降到 30 GB 后重启。

**11-11 TE 重新占用 ont GPU**：按用户 11-08 的优先级，TE 排第一。
- 00:30–12:12：TE 的 GPU 作业在跑，P6 在 Amborella NC_142343.1 之前的闸门处等待。
- 12:12–12:26：TE 两臂之间 CLAIM 消失了 14 分钟，P6 跑完了两个各约 7 分钟的文件。
- 12:26 起：TE 臂 C 运行（TE 预计 11-12 约 01:15 结束），之后还有 PLANT-FT 约 9–12 GPU h。
- 为避免长步骤（如 Vipsania 微调约 3.5 h、22 GB）插进 TE 两臂之间的空档、压住 TE 下一臂，闸门新增安静期：最近一次看到 TE（CLAIM、GPU 进程、队列脚本、外来容器）满 30 分钟后才开始新步骤。
- 此时 Carbon-A 已完成 9 个物种中的 8 个（Amborella 只剩 unplaced 文件），还剩 Hakea、Themeda。

**11-12 段错误的规律**（dmesg）：
- P6 在 11-10 的三次 Carbon-A SIGSEGV（06:19、08:33、10:35）都发生在逻辑 CPU 8（core 16）。
- TE 在 11-12 的两次分别在 CPU 8 和 CPU 0。
- 11-05 至 11-11 其他项目的 python 段错误（Python 3.9 与 3.10 都有）都在 CPU 0。
- ont 是 i9-13900K（微码 0x133），CPU 8 是 5.8 GHz 的最高睿频核。不同项目、不同程序只在这两个 CPU 上崩溃，疑似 CPU 不稳定（这一代有已知问题），不是定论。
- P6 的 GPU 步骤 CPU 集合因此去掉了 8（原为 2,4,6,8,10,12,14,16–24），不影响计算。
- 已告知 TE 会话与维护会话；是否检查硬件由用户定。
- TE 看到的 `extract_p1_shard.py` 短 GPU 任务属于 ncrna-annot，不是 P6。P6 自 11-11 12:26 起在闸门处等 TE。

**Carbon-A 汇总**：
- 10 个物种全部完成。
- 4 次单次推理失败都在重试后成功：3 次 SIGSEGV、1 次 tokenizer 报错，全部发生在 11-10，之后去掉 CPU 8 就没有再出现。
- 只有 Rhynchospora 的 unplaced 文件（不计分）触发了密码子表警告，没有主序列触发。

**Vipsania**（11-12 15:19 起）：
- Dryas 的微调用 `--finetune_B 2` 第一次尝试时显存不足：XLA 一次申请 25 GB，超过了卡的 24 GB。按报告 30 §3 的规则降到 `--finetune_B 1` 重试，有效批仍是 64。
- M38 在拟南芥、水稻上用 B = 2 能跑（峰值 21.9 GB），Dryas 为什么不行未查明。
- B = 1 实测 177 ms/步，每个 epoch 约 19 分钟，10 个 epoch 约 3.2 h，与 M38 的 B = 2 用时相近。
- 由于驱动在一次 OOM 后就保持 B = 1，之后各基因组都用 B = 1，各自记在 `run.json`。

**11-13 Rhynchospora 注释阶段 OOM**：
- 经过：Vipsania 微调已跑满 10 个 epoch（06:23–09:44）；随后的全基因组注释于 09:45 因显存不足失败，原因是 TE 的两个 GPU 作业（OPOSSUM 冒烟、te-glm release 检查）在 P6 这一步运行中启动。
- 微调检查点 `finetuning_1/latest_checkpoint.weights.h5` 完整保留。
- 补救：驱动新增"微调完成则只补注释"：若失败前日志已有 10 个 epoch 的结束行且检查点存在，就写 `FT_COMPLETE`，重试时只运行 `vipsania annotate finetuning_<k> <genome> --model_dir <dir>`。这与一次性 `--finetune` 调用用的是同一个微调后模型，`run.json` 记 `annotated_from_checkpoint`。
- Rhynchospora 已手工照此标记，省去约 3.4 GPU h。

**11-13/14 Vipsania 进度**：
- Rhynchospora：用 `finetuning_1` 只补做注释，18:53 完成。实际微调发生在第一次尝试，用的是 `--finetune_B 1`。
- Zizania：B = 2 时显存不足，降到 B = 1 后 22:49 完成。
- 11-14 00:59 TE 起了一个约 1 h 的推理，与 P6 的 Akebia 微调同卡并行（12.9 + 4.5 GB），没有冲突。
- TE 随后排了 X2-BENCH（约 8 h）和 X1（约 22 h），接在 Akebia 之后。按"TE 优先"，P6 剩下的 Amborella、Hakea、Themeda（约 11 GPU h）约在 11-15 恢复。是否改变先后，已请维护会话交给用户定。

**11-14 晚 ont 磁盘写满**：
- 原因：TE 的 RepeatModeler `elements.db` 涨到 100 GB 以上（维护会话查明，已通知用户）。20:59 根分区余量 1.4 GB，随后降到 639 MB。
- 当时 Amborella 的 Vipsania 微调在第 8/10 个 epoch，第 7 个 epoch 的检查点（20:40）完整，401,637,584 字节，与 Akebia 的相同。
- 处理：
  - 为免第 8 个 epoch 结束时在满盘下覆盖写坏检查点，20:59:59 对该进程发 SIGSTOP 暂停；
  - ont 上另起小脚本，余量连续两次 ≥ 40 GB 时自动 SIGCONT 继续；
  - 驱动的"只补注释"判定收紧为：检查点 ≥ 400 MB 且能被 h5py 打开，否则从头微调。

**11-14 21:05–23:48**：
- 磁盘恢复（余量 123 GB）后，自动恢复脚本于 21:05 让 Amborella 继续；22:58 完成微调，15,392 个基因。
- 之后闸门把 TE 留下的空声明文件（22:09 起 0 字节，TE 已撤掉自己那一行）当成占用，GPU 空闲了 22:58–23:48。
- 修正：声明文件有非空白内容才算占用；"最后看到 TE 的时间"记为该文件清空的时刻。

**可信度核对**：P6 在 Dryas 四条 set B 染色体上的 Carbon-A 预测，与 M30 用同一模型、同一序列的输出逐条完全相同，共 2,844 / 2,480 / 2,220 / 1,988 条链，无任何差异。这只比较两份预测，不读参考。说明流程是确定的，重试成功的输出可信。

### 5.3 运行汇总

| 工具 | GPU 墙钟（不含等待） | 说明 |
|---|---:|---|
| ANNEVO | 1.3 h | 一次通过 |
| Tiberius | 5.3 h | 一次通过 |
| Helixer | 7.8 h | 一次通过 |
| Carbon-A | 约 26 h（推理 + CPU 后处理） | 4 次单次失败重试成功；1 个不计分文件触发密码子表警告 |
| Vipsania | 约 38 h（10 个基因组各一次微调 + 注释，外加 8 次 pre） | 微调都以 `--finetune_B 1` 完成（B = 2 时显存不足）；Rhynchospora 由检查点补做注释；Dryas、Eschscholzia 的 pre 复用 M38 |

逐物种的耗时、资源与基因数在各运行目录的 `run.json` 里（ont `~/coding-rna/outputs/PANEL-P6/{native,carbonA,vipsania}/`）。

## 6. 面板结果与 C1 判定

打分只做了一次：2026-11-17 07:50，在 ont 上运行 `p6_score.py panel`，用时 8 分钟，只用 CPU。打分前核对过，8 个参评对象在 10 个物种上全部完成。

结果文件：`experiments/PANEL-P6/results/p6_panel_result.json`，以及 `p6_panel_table.tsv`。

### 6.1 C1 判定（报告 28 §6，判定规则见本文 §2.4）

| 对象 | 面板中位数 [95% CI] | T2 − 对象 [95% CI] | 判据中的角色 | 是否满足 |
|---|---|---|---|---|
| **T2（我方）** | **55.05** [54.23, 55.78] | — | 被判定方 | — |
| Helixer a0080 | 42.10 [41.62, 42.75] | +12.95 [+12.06, +13.56] | T2 须 ≥ | ✓ |
| Tiberius angiosperms | 59.25 [58.51, 59.77] | **−4.20** [−5.06, −3.16] | T2 须 ≥ | ✗ |
| Vipsania ft | 57.41 [56.95, 57.96] | **−2.36** [−3.35, −1.56] | T2 须 ≥ | ✗ |
| Carbon-A 过滤 | 67.65 [67.05, 68.25] | **−12.60** [−13.65, −11.65] | 差距须 ≤ 6.0 | ✗（差距 12.60 > 10） |
| ANNEVO（描述） | 63.08 [62.50, 63.73] | −8.03 [−9.08, −7.14] | — | — |
| Carbon-A 未过滤（描述） | 65.73 [65.10, 66.35] | −10.68 [−11.75, −9.67] | — | — |
| Vipsania pre（描述） | 54.69 [54.30, 55.27] | +0.36 [−0.63, +1.05] | — | — |

CI 来自联合 1-Mb 块 bootstrap（2,000 次，§2.3），"T2 − 对象"是面板中位数之差。

**判定：C1 不成立**。
- 依据大纲 §0：与 Carbon-A 过滤版的差距超过 10 点。
- 另外两条也未满足：T2 的中位数低于 Tiberius 和 Vipsania ft。
- 大纲"不成立"里还有一条"面板上 T2 相对 NT-v2 + LoRA 的优势消失"。NT-v2 + LoRA 没在面板上跑，这一条没有评估，不影响判定。

**敏感性**（§2.7，不改判定）：
- Tiberius 只按其适用范围内的 9 种（去掉 Amborella）计算：T2 中位数 56.61，Tiberius 59.51，仍为不成立。
- 去掉 CDS ≥ 96 bp 规则：各中位数基本不变（T2 54.99，其余与上表一致或相差 < 0.6），判定不变。

### 6.2 补跑 seed（报告 28 §7）

**触发**：
- 条件 1 成立：T2 对 Helixer、Tiberius、Vipsania ft 的最小领先为 −4.20，小于 2.5。
- 条件 2 不成立：与 Carbon-A 过滤版的差距 12.60，不在 3.5–8.5 之间。

**§7 要求的做法**：在 ont 上对 T2_s2、T2_s3 各跑一次面板推理，同一路径、同一评分；C1 仍按 T2 判定。P5 的 T2 面板推理实测约 34 h 一次，两次合计约 68 GPU h。

**补跑不会改变 C1 判定**：
- 差距 12.60 比"不成立"线（10 点）高 2.6 点。
- 三个 seed 在跨物种评估上的极差约 2.5–3.0 点（报告 28 §7）。
- 即使三 seed 均值上移这么多，与 Carbon-A 的差距仍大于 10。

所以补跑的作用只是给论文报告 seed 离散度。是否、何时排这约 68 GPU h，由维护会话与用户决定；本报告没有启动。

### 6.3 逐物种 exact-chain F1（CDS ≥ 96 bp 规则）

| 物种 | 层 / 注释来源 | T2 | Helixer | Tiberius | ANNEVO | Carbon-A 未过滤 | Carbon-A 过滤 | Vipsania ft | Vipsania pre |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Dryas octopetala | TH / RefSeq | 69.47 | 54.04 | 69.65 | 72.49 | 72.95 | 74.66 | 64.93 | 62.47 |
| Eschscholzia californica | TH / RefSeq | 67.94 | 48.43 | 66.86 | 71.29 | 71.32 | 72.89 | 64.13 | 61.05 |
| Akebia trifoliata | TH / RefSeq | 53.49 | 40.49 | 60.87 | 67.56 | 68.59 | 70.58 | 59.90 | 57.05 |
| Hakea prostrata | TH / RefSeq | 49.03 | 40.28 | 58.99 | 65.87 | 68.44 | 70.69 | 56.78 | 53.86 |
| Zizania latifolia | SH / GenBank | 56.61 | 48.11 | 53.32 | 60.81 | 61.59 | 62.45 | 56.39 | 53.90 |
| Themeda triandra | SH / GenBank | 34.96 | 27.57 | 31.87 | 35.88 | 36.86 | 38.09 | 32.93 | 33.41 |
| Juncus effusus | SH / GenBank | 64.06 | 50.77 | 59.51 | 64.77 | 66.10 | 68.37 | 61.47 | 59.29 |
| Rhynchospora breviuscula | SH / GenBank | 60.96 | 43.71 | 59.80 | 60.46 | 62.32 | 64.45 | 58.03 | 55.48 |
| Acorus gramineus | SH / GenBank | 40.87 | 33.67 | 41.80 | 43.76 | 44.14 | 45.08 | 41.87 | 39.85 |
| Amborella trichopoda | AH / RefSeq | 38.91 | 32.78 | 51.97 | 61.40 | 65.36 | 66.94 | 52.64 | 49.03 |
| **中位数** | | **55.05** | 42.10 | 59.25 | 63.08 | 65.73 | 67.65 | 57.41 | 54.69 |
| 均值 | | 53.63 | 41.99 | 55.46 | 60.43 | 61.77 | 63.42 | 54.91 | 52.54 |
| 最小值 | | 34.96 | 27.57 | 31.87 | 35.88 | 36.86 | 38.09 | 32.93 | 33.41 |

一致性核对：T2 在 Dryas、Eschscholzia 整个物种上为 69.47、67.94，与此前 set B 的 68.57、68.39 一致。每条序列都有预测，未见工程故障。

### 6.4 物种层面配对检验（报告 10 §4.1；10 个物种，Holm 校正 6 个工具）

| T2 对 | 胜 / 负 | 物种差中位数 | Wilcoxon p（Holm） | 符号检验 p |
|---|---|---:|---|---:|
| Helixer | 10 / 0 | +10.87 | 0.0020（0.014） | 0.0020 |
| Tiberius | 5 / 5 | +0.45 | 0.85（1.00） | 1.00 |
| ANNEVO | 1 / 9 | −3.18 | 0.0039（0.016） | 0.022 |
| Carbon-A 未过滤 | 0 / 10 | −3.43 | 0.0020（0.014） | 0.0020 |
| Carbon-A 过滤 | 0 / 10 | −5.07 | 0.0020（0.014） | 0.0020 |
| Vipsania ft | 6 / 4 | +1.13 | 0.92（1.00） | 0.75 |

T2 与 Tiberius、Vipsania ft 逐物种互有胜负，配对差的中位数甚至略为正。但 T2 的物种分数分散更大，所以按 C1 的口径（各自取物种中位数再比）T2 落后。

### 6.5 描述：落后集中在长基因物种（未做受控分解）

| 物种 | T2 精确率 / 召回 | Carbon-A 过滤 精确率 / 召回 | T2 预测数 / 参考数 | ≥ 6 段基因召回：T2 / Carbon-A 过滤 / Tiberius / Vipsania ft |
|---|---|---|---:|---|
| Amborella | 35.6 / 42.8 | 70.8 / 63.5 | 1.20 | 29 / 58 / 51 / 48 |
| Hakea | 44.1 / 55.3 | 72.1 / 69.3 | 1.25 | 42 / 67 / 64 / 57 |
| Akebia | 51.7 / 55.4 | 75.2 / 66.5 | 1.07 | 49 / 67 / 65 / 60 |
| Dryas（对照） | 68.8 / 70.1 | 77.2 / 72.3 | 1.02 | 68 / 72 / 74 / 65 |
| Juncus（对照） | 60.4 / 68.2 | 66.2 / 70.6 | 1.13 | 68 / 71 / 65 / 71 |

- **在哪里落后**：T2 落后最多的三个物种（Amborella −28.0、Hakea −21.7、Akebia −17.1，对 Carbon-A 过滤版）正是报告 25 §3 指出的长基因、低基因密度物种。CDS 跨度超过 24,576 bp 的参考基因，Amborella 占 10.1%，Hakea 占 7.2%，Akebia 占 3.0%；Dryas 只有 0.3%。
- **表现形态**：在这三个物种上，T2 的预测数多于参考数，精确率低，而且 ≥ 6 段基因的召回大幅下降。在基因紧凑的物种上，这一层与 Carbon-A 接近。这像是长基因被切成几段，与报告 25 §3 预见的长基因结构性难例一致。
- **注释来源分层**（描述）：
  - SH 层（GenBank 注释，5 种）：T2 56.6，Tiberius 53.3，Vipsania ft 56.4，Carbon-A 过滤 62.4，T2 与 Carbon-A 的差距 5.8。
  - RefSeq 注释的 5 种：T2 53.5，Carbon-A 过滤 70.7，差距 17.2。
  - Carbon-A 的训练标签就是 RefSeq（Gnomon）风格，这与报告 10 §7.1 第 2 条的风险方向一致。
  - 这些只是描述，不改变判定。
- **禾本科**：禾本科 2 种 T2 45.8，非禾本科 8 种 T2 57.2；Carbon-A 过滤分别为 50.3 和 69.5。
- **注释质量**：Themeda 与 Acorus 所有工具都 ≤ 45，与报告 25 §6 的注释质量警示一致，对各工具的影响相近。
- **暴露**（§4）：
  - Vipsania 的预训练序列含 Akebia、Hakea、Juncus、Amborella；
  - Carbon-A 在 Amborella 上为 S（同物种）；
  - T2 落后最多的三个物种里有两个是 Vipsania 的预训练物种，一个是 Carbon-A 的训练物种。

### 6.6 set B 上的 Vipsania ft（补全报告 30 §8）

| 模式 | set B pooled | Dryas / Eschscholzia | 对 T2 | 对 Vipsania pre |
|---|---:|---|---|---|
| Vipsania ft | 64.36 | 64.03 / 64.61 | −4.11 [−4.69, −3.49] | +2.67 [+2.31, +3.02] |
| Vipsania pre（M38） | 61.68 | 61.86 / 61.55 | −6.79 [−7.41, −6.13] | — |

评分方式同报告 30 §4（`vip_score.py setb`，结果 `experiments/PANEL-P6/results/vip_setb_ft/vip_setb_result.json`）。

### 6.7 对论文档次的含义（按预注册规则，不在此改写大纲）

按大纲 §0，C1"不成立"的后果是：
- C1 回到 DEV 排序，零样本差距只在局限里写；
- 主体改为 C2 + C3 + 骨干分解。

面板上的稳定结论有三条：
- T2 在 10 个物种上全部高于 Helixer；
- 与 Tiberius、Vipsania ft 逐物种互有胜负；
- 在 10 个物种上全部低于 Carbon-A，差距集中在长基因物种。

大纲与看板的修改由维护会话与用户决定。
