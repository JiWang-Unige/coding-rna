# M38-VIPSANIA：Vipsania（无监督基因预测）在 DEV + 评估集 B 上的对比

状态（2026-11-04）：**CPU 准备完成**（环境、模型、基因组、评分器及其校验、CPU 烟测、GPU 作业脚本）。**GPU 步骤未运行**；`READY` 未创建；没有任何 Vipsania 数字在 DEV 或集 B 上产生。

协议与用户决定见 `docs/review_20261016/28_animal_C2_and_Vipsania_proposal.md` §1、§6（报告 30 将写完整协议）。要点：Vipsania 进入 C1 面板判据，判据模式是 Streptophyta 模型 `j9m0cdmk` + `--finetune`（作者默认参数）；只用预训练的模式同时跑，描述性，不进判据。本实验先测 DEV（拟南芥、水稻）与集 B（Dryas、Eschscholzia），每个整基因组各评一次。

## 1. 版本与路径（ont）

| 项 | 值 |
|---|---|
| 环境 | `~/envs/vipsania_m38`（uv 0.11.14 建，Python 3.12.12）。`~/envs/vipsania`、`~/envs/vipsania_ncrna` 是另一项目（ab-initio-TE / ncrna-annot）的，未触碰 |
| vipsania | **1.1.0**（`pip install vipsania`，PyPI） |
| tensorflow / keras / numpy | **2.19.1**（pip wheel，带 CUDA 12.5 pip 库，cudnn 9.3.0.75）/ 3.15.1 / 2.1.3；bricks2marble 0.1.3，hidten 1.0.0 |
| 作者容器 | NGC TF 2.17 + CUDA 12.8；我们用 pip 路线（TF 2.19.1），与另一项目 `vipsania_ncrna` 环境同版本（其 `~/envs/vipsania` 是 CPU 版 tensorflow-cpu） |
| 驱动 / GPU | 595.58.03（≥ 570）/ RTX 4090 24 GB |
| 模型 | clade `Streptophyta` → **`j9m0cdmk`**（`versions.json` 2026-08-20；作者 locus F1 0.602 微调 / 0.579 预训练）。脚本只用 ID。25,047,166 参数，T=20,000，B=8 |
| 模型缓存 | `~/coding-rna/outputs/M38-VIPSANIA/vipsania_cache/models/j9m0cdmk/`（`VIPSANIA_CACHE`，与另一项目的 `~/.cache/vipsania` 隔离）；sha256 见 `model_j9m0cdmk.sha256` |
| 版本记录 | `~/coding-rna/outputs/M38-VIPSANIA/env_versions.json`、`env_freeze.txt`（77 个包） |
| 输出根 | `~/coding-rna/outputs/M38-VIPSANIA/`（`genomes/`、`smoke/`、`validation/`；GPU 步骤将写 `runs/`、`results/`、`events.jsonl`、`DEVB_DONE`） |

GPU 可见性：按要求**没有**用 TF 列 GPU（TF 列设备会初始化 StreamExecutor，可能建 CUDA 上下文，占显存）。所有 CPU 步骤设 `CUDA_VISIBLE_DEVICES=""`，日志确认 `no CUDA-capable device`，没有 GPU 进程。旁证：12:05 观察到另一项目的 `vipsania_ncrna` 环境（同为 TF 2.19.1 wheel）在这块卡上占着约 7 GB 运行。真正的可见性检查是 `vip_devb.sh` 的第一步（GPU 烟测，要求日志出现 `Using N GPU`，否则 exit 20）。

## 2. 基因组（整套核基因组；每个 FASTA 保持 NCBI 软屏蔽 = 小写，作者的默认输入）

`~/coding-rna/outputs/M38-VIPSANIA/genomes/<key>.nuclear.fa`，清单 `genomes_manifest.json`（sha256、序列长度、与评分参考序列的 md5 比对）。由 `scripts/vip_make_genomes.py` 生成（原样复制序列与屏蔽，去掉细胞器；M30/M34 的工作 FASTA 不是整基因组：Dryas / Eschscholzia 的含 20-Mb 截取、水稻 `genome.fa` 只有 5 条染色体）。

| key | 组装 | 序列数 | 总长 | 文件 | 评分序列 | 来源 |
|---|---|---:|---:|---:|---|---|
| arabidopsis_thaliana | GCF_000001735.4（TAIR10.1） | 5 | 119,146,348 | 121 MB | **DEV** NC_003074.8 | `data/m1_screen/arabidopsis_thaliana/genome.fa`（去 NC_037304.1、NC_000932.1） |
| dryas_octopetala | GCF_963921425.1（drDryOcto1.1） | 25（9 染色体 + 16 未定位） | 229,579,301 | 232 MB | **B** NC_141367.1、NC_141369.1、NC_141371.1、NC_141373.1 | `data/m30/raw/GCF_963921425.1_…_genomic.fna.gz` |
| oryza_sativa | GCF_034140825.1（ASM3414082v1） | 12 | 385,710,679 | 391 MB | **DEV** NC_089041.1 | `data/m1_screen/oryza_sativa/genome.full.fa`（去线粒体、质体、质粒 B1） |
| eschscholzia_californica | GCF_041429995.1（dmEscCali1.0.p） | 329 个 scaffold | 400,665,884 | 406 MB | **B** NW_028238309.1、NW_028238461.1、NW_028238483.1 | `data/m30/raw/GCF_041429995.1_…_genomic.fna.gz` |

评分序列的序列在两处逐字相同（转大写后 md5，清单里 `scored_sequence_md5_equal_to_scoring_reference: true`）。

**与任务书的出入（需要确认）**：任务书把 DEV 写成拟南芥 NC_003071.7 + 水稻 NC_089044.1。代码里的 DEV 是 **NC_003074.8 + NC_089041.1**（`m29_common.DEV_SCOPE`、`m31_common.DEV_SEQS`、`m34_boot_pair`；已发表的 70.107 / 77.41 / 7,728 条参考都是这两条）。NC_003071.7 / NC_089044.1 是 M31/M34 的**内部验证**染色体（选点用，不是 DEV）。评分器按代码里的 DEV 评；两对染色体都在同一基因组里，Vipsania 整基因组注释时都会被注释，若要改评别的染色体只需重跑 CPU 评分。

## 3. Vipsania CLI 与输出

`vipsania annotate <model_id|clade> genome.fa -o out.gff3 [--finetune]`（输出格式由后缀定：`.gff3/.gff`、`.gtf`、`.gp`）。关键选项（`vipsania annotate --help`、`docs/annotate.md`、`docs/training.md`，仓库提交 397660e）：

| 选项 | 默认 | 说明 |
|---|---|---|
| `-T/--context` | 200,000 | 注释块长（nt）；不改 |
| `-B/--batch_size` | 按空闲显存推断 | 注释批大小（块数）。作者公式：每块约 2.25 GiB（25M 模型、200 kb），`floor(0.9 × 空闲 / 2.25)`，24 GB 卡约 9 |
| `--finetune` / `--finetune_epochs` / `--finetune_lr` | 关 / 10 / 1e-4 | 目标基因组上的无标签微调（MLM + HMM 损失），完成后用微调模型注释；检查点存在输出文件旁（`finetuning_1/`，约 400 MB） |
| `--finetune_B` | 按空闲显存推断 | 微调的显存批大小，必须整除 64；有效批大小恒为 64（梯度累积 64/B），每 epoch `train_steps = 100×(64/B)`。作者公式：每个 20-kb 窗口约 10.1 GiB，24 GB 空卡得 2，空闲不足 11.3 GiB 时直接崩（返回 None） |
| `-fo KEY=VALUE` | | 覆盖微调配置（我们只在烟测里用） |
| `--drop_repeats`、`--group_limit`、`--delta`（0.1）、`-p`、`--nojit` | | 重复序列取样 / 显存友好度 / 多小序列 / HMM 并行度 |

微调的计算量**与基因组大小无关**：10 epoch × 6,400 窗口 × 20 kb（共 64,000 个窗口，约 1.3 Gnt）。作者在 A100 80 GB 上约 160 min。显存：作者文档写明 24 GB 卡上 25M 模型训练 B=1。`vipsania annotate` 自己会设 `TF_GPU_ALLOCATOR=cuda_malloc_async`。

**输出 GFF3 结构**（CPU 烟测实测）：`##gff-version 3`（无 `##sequence-region`）；每个基因 `gene`（`ID=g1;Name=g1`）→ `mRNA`（`ID=g1_t1;Name=g1_t1;Parent=g1`）→ 每个编码块一行 `CDS` 与一行同坐标 `exon`（`ID=g1_t1.CDS.1;Parent=g1_t1`、`g1_t1.exon.1`）。无 UTR（exon = CDS）、无得分（`.`）、source 为 `Vipsania`。每个基因恰好一个转录本。CDS phase 为 0/1/2（按链 5' 端累计，负链也正确）。**CDS 含终止密码子**（GTF 版本里 stop_codon 取最后 3 个碱基且位于最后一个 CDS 内），与 RefSeq 参考一致，所以精确链匹配无需调整。序列名在第一个空白处截断。输出文件旁有同名 `.log`（含模型、批大小、用时、`peak memory` / GPU 峰值显存）。后处理自动删除 CDS < 9 bp 与框内终止密码子的基因。

## 4. 评分器 `scripts/vip_score.py`

不重新实现指标，只做格式转换并调用现有代码：

- **转换** `read_gff` → `write_contig`：流式读 Vipsania GFF3（gene / mRNA / CDS，Parent 链），按评分序列写成 `<run>/<contig>/prediction/full_grammar.gff3`（`m34_infer` 布局；gene / mRNA / exon / CDS，基因与转录本 ID 原样保留）。内置独立的 stdlib 链推导（每基因取最长 CDS 的转录本，键 = `(seqid, strand, ((start0, end), …))`），与项目 ruler 的链集合逐条核对。
- **DEV**：池化 exact-chain F1 与分区数（`a1_eval.load_chains` + `summarize`，CDS ≥ 96 bp 规则，M26 CDS-assessable 主参考，DEV 7,728 条）；配对 1-Mb 块 bootstrap 直接运行 `m34_boot_pair.py` 本身（2000 次、seed 0，Vipsania 作参考 run），对 T2 / T10 / k50U / k50F（有则加 T2_s2）以及它内置的 GENERanno fresh-joint、Carbon-A 过滤版。另报不加 CDS ≥ 96 规则的池化值（已发表外部行 70.107 / 77.41 走的是这条路径，见 §5）。
- **集 B**：`m34_p2b_score.run_scoring`（`m30_score.Evaluator` + `load_gff_chains`，同一 bootstrap），Vipsania 作 `new`，对 T2、T10（有则加 T2_s2、另一个 Vipsania 模式）；P2b 判据（T10−T2 ≥ +2.0）不适用，已从输出去掉。每个 run 的 provenance 取自 GFF 旁的 `run.json`。
- 结果文件 `vip_dev_result.json` / `vip_setb_result.json` 已存在时拒绝覆盖（DEV 与集 B 预先登记为每个 run 评一次），`--overwrite` 才重评。

```
# 转换（stdlib）：  python3.12 vip_score.py convert --gff F --out DIR --contigs A,B [--slice slice.json]
# DEV / 集 B（ont，generanno 环境，CPU）：
nice -n 10 taskset -c 18-23 ~/envs/generanno/bin/python -I vip_score.py dev  --out OUT --ft arabidopsis_thaliana=GFF oryza_sativa=GFF --pre arabidopsis_thaliana=GFF oryza_sativa=GFF
nice -n 10 taskset -c 18-23 ~/envs/generanno/bin/python -I vip_score.py setb --out OUT --ft dryas_octopetala=GFF eschscholzia_californica=GFF --pre dryas_octopetala=GFF eschscholzia_californica=GFF
# 校验：             … vip_score.py validate --out ~/coding-rna/outputs/M38-VIPSANIA/validation
# 烟测解析检查：     … vip_score.py smoke --gff smoke.gff3 --slice slice.fa.json --out-dir DIR
```

## 5. 评分器校验（2026-11-04，`validation/validation.json`，31/31 通过，28 s）

用已评分的外部 GFF 经**同一转换器**和同一评分路径复算，与已发表值逐位比较（浮点 1e-9）：

| 检查 | 复算 | 已发表 |
|---|---:|---:|
| GENERanno fresh-joint DEV，经转换器，不加 CDS≥96 规则 | **70.10673826206535** | 70.107（`boot_pair_*_DEV.json` 同一数） |
| 同上，原路径 `AC.gff_chains` | 70.10673826206535 | 同 |
| 同上，加 CDS≥96 规则（我方 arm 用的规则）= 原链集按规则过滤 | 70.17635874909854（+0.070） | — |
| Carbon-A 未过滤 DEV，经转换器 = `AC.carbon_dev` | **75.24128422296631** | 75.24 |
| Carbon-A 过滤 DEV（conf ≥ 0.1），经转换器 | **77.40574384102355** | 77.41 |
| k50U DEV GFF 经转换器 | 71.63115889453202 | 71.631（`boot_pair_k50U_DEV.json`） |
| bootstrap（共 54 块）：T10 F1 与 CI95 | 72.91357691789133 [71.0635, 74.5037] | 逐位相同 |
| T2 F1；T10 − T2 点估计与 CI95 | 73.87819414948277；−0.96462 [−1.51639, −0.28897] | 逐位相同 |
| T10 − Carbon-A 过滤（转换后）CI95 | [−5.28665, −3.68412] | 逐位相同 |
| 集 B：T10 / T2 的 B 集 GFF 经转换器 + `run_scoring`：池化、Dryas、Eschscholzia | 68.329 / 68.470；68.107 / 68.570；68.502 / 68.392 | `p2b_result.json` 逐位相同 |
| 集 B：T10 − T2 池化点估计与 CI95 | −0.14101 [−0.53061, +0.27747] | 逐位相同 |

两条路径的差别要知道：已发表的外部行（fresh-joint 70.107、Carbon-A）走 `AC.gff_chains`，**不加** CDS≥96 规则；我方 arm（T2/T10/k50…、集 B）走 `E.load_chains` / `SC.load_gff_chains`，**加**规则。评分器默认用后者作主数字（与要对比的 T2/T10 同规则），并另报前者。对 Vipsania，这条规则丢掉 CDS 9–95 bp 的预测（它自己的下限是 9 bp）。`tests/dryrun_score.sh` 用已评分的非 Vipsania 预测做 `dev` / `setb` 全流程演练（输出结构、覆盖保护）。

## 6. CPU 烟测（GPU 隐藏；训练区段，不是评估）

切片：拟南芥训练染色体 NC_003070.9 的 5,000,001–5,500,000（500 kb；小写比例 16.6%）。核 `smoke/`。

| 步骤 | 用时 | 峰值内存 | 输出 |
|---|---:|---:|---|
| 预训练注释（`-B 1`，7 核 nice 19，与其他任务共用 CPU） | 337 s（Vipsania 记 329 s；3 块初始预测 204 s，重预测 2 块 125 s） | 5.1 GB | 143 个基因，785 个 CDS 块，GFF3 171 KB |
| `--finetune` 管线检查（1 epoch，2 个 1 窗口微批，累积 2；只验证代码路径与检查点，不是真微调）+ 注释 | 799 s（XLA 编译约 3.5 min；训练步 227 s/步为 CPU 编译主导） | 26.9 GB | 143 个基因（权重几乎未变）；检查点 `ft/finetuning_1/latest_checkpoint.weights.h5` 401 MB |
| 去软屏蔽的预训练注释（旁证） | 375 s | 5.1 GB | 141 个基因；解析检查 tp 106、F1 71.86（屏蔽版 71.38） |

解析检查（**仅检查转换器与精确链匹配，训练区段，不是评估**；参考为落在切片内的 154 条主参考，预测 143 条均在切片内）：预训练 tp 106、precision 74.1、recall 68.8、F1 **71.38**；管线检查的 `ft` 输出（2 个微批后权重几乎未变）tp 105、F1 70.71；stdlib 链推导与项目 ruler 加载结果在三个输出上都一致。106 条精确匹配也证明了坐标约定（CDS 含终止密码子，1-based 闭区间）无需调整。

## 7. GPU 作业 `scripts/vip_devb.sh`（未运行）

由队列 `m34_after_t2s2.sh` 在 `M38-VIPSANIA/READY` 出现后调用（T2_s2 与 T2_s3 之间），成功标志 `DEVB_DONE`。顺序：

0. 门：没有 `READY` 就拒绝（exit 11）；`SKIP` / `DEVB_DONE` 已存在则安静退出；`flock` 单实例；核对模型 sha256 与基因组大小。
1. **GPU 烟测**：2-Mb 训练切片（NC_003070.9 5,000,001–7,000,000）上预训练注释，`TF_FORCE_GPU_ALLOW_GROWTH=true`，日志必须出现 `Using N GPU`，否则 exit 20（队列不阻塞）；再做 512 窗口的微调探针，选定 `--finetune_B`（空闲 ≥ 23,040 MiB 用 2，否则 1；OOM 则降到 1），测每窗口秒数，写事件 `gpu_smoke_done`（含 `projected_finetune_hours_per_genome_upper`）。若上界投影超过 `VIP_MAX_PROJ_H`（默认 6 h/基因组：用户批准 DEV + B 约 15–20 GPU h，4 个基因组的微调上界合计约 24 h；超出则停下问用户）则 exit 25，不浪费槽位。
2. 按基因组从小到大：拟南芥（DEV）→ Dryas（B）→ 水稻（DEV）→ Eschscholzia（B）。每个基因组先做 `pre`（只预训练），再做 `ft`（`--finetune`，作者默认：10 epoch、lr 1e-4、T 200,000、有效批 64）。模型按 ID `j9m0cdmk` 固定。**唯一由 24 GB 强制的偏离**：`--finetune_B`（有效批不变），以及注释 OOM 后 `-B` 降为 4 → 2 → 1。步骤哨兵 `runs/<g>/{pre,ft}/*_DONE` 与基因组哨兵 `runs/<g>/DONE`，可续跑（半途的微调会重来）。守卫 `vip_guard.py`：宿主 RSS > 48 GB 或超时（ft 16 h、pre 5 h）则杀进程组，并记录峰值 RSS 与峰值 GPU MiB。
3. 全部完成后用 `vip_score.py dev` 与 `setb`（CPU，`taskset 18-23`）评 `ft` 与 `pre`，写 `DEVB_DONE`（含各池化 F1）。

日志：`~/coding-rna/logs/m38-vipsania-<步骤>.{out,err}`；事件：`~/coding-rna/outputs/M38-VIPSANIA/events.jsonl`（`step_start` / `step_end`，含秒数、峰值 GPU MB、峰值 RSS）。查看进度：`tail -n 5 ~/coding-rna/outputs/M38-VIPSANIA/events.jsonl`。

**GPU 小时估计（外推，GPU 烟测会给实测值）**：微调与基因组大小无关，每基因组 64,000 窗口。作者 A100 约 160 min；4090 在 B=1–2 下我估 1.5–2.5 倍（带宽约为 A100 的一半，FP32 算力更高，微批更小）→ 约 4–6.5 h，保守范围 2.5–8 h。4 个基因组共 **10–32 h，中值约 18 h**。注释：作者 144 Mb 不到 6 min（≥ 24 Mb/min）；我们每个模式 1,135 Mb、共 8 次注释（2,270 Mb）→ 约 1.6–3.9 h。烟测与启动开销约 0.3 h。**合计约 12–36 GPU 小时，中值约 22 h**（只跑预训练模式的话仅约 1–2 h）。

## 8. 文件

| 路径 | 内容 |
|---|---|
| `scripts/vip_make_genomes.py` | 建四个整核基因组 FASTA 与清单 |
| `scripts/vip_score.py` | 转换器 + DEV / 集 B 评分 + 校验 + 烟测解析检查 |
| `scripts/vip_devb.sh` | GPU 作业（见 §7） |
| `scripts/vip_guard.py`、`vip_event.py`、`vip_slice.py`、`vip_env.sh` | 资源守卫、事件行、切片、环境 |
| `scripts/vip_smoke_cpu.sh` | CPU 烟测（`pre|ft|both|unmasked`） |
| `tests/test_vip_convert.py` | 转换器单元测试（本地即可：`python3 -m unittest experiments/M38-VIPSANIA/tests/test_vip_convert.py`） |
| `tests/test_devb_stub.sh` | `vip_devb.sh` 控制流测试（桩 vipsania，无 GPU：门、哨兵、续跑、事件、投影闸、失败退出码） |
| `tests/dryrun_score.sh` | `dev` / `setb` 全流程演练（已评分的非 Vipsania 输入） |
| `push.sh` | 推送到 ont（不要给两层 ssh 加 `</dev/null`） |

## 9. 已知风险 / 待你决定

1. **DEV 染色体的写法**见 §2 的出入。
2. **软屏蔽**：输入保持 NCBI 小写屏蔽（作者默认）。作者建议屏蔽质量不确定时去掉；500 kb 训练切片上去屏蔽几乎没有差别（F1 71.86 对 71.38，tp 同为 106，基因 141 对 143），所以保持默认。但这只是单个切片的旁证，拟南芥屏蔽比例仅 16.6%（水稻 34%、Dryas 28%、Eschscholzia 34%），不能外推到重复多的基因组；若想把“去屏蔽”预登记为选项，需在创建 `READY` 前定。
3. **微调显存**：B=1 是最小，若 24 GB 仍 OOM，没有进一步退路（不改 T 与数据设置）；GPU 烟测会先暴露。另一项目的 `vipsania_ncrna` 作业此刻占着这块卡约 7 GB，槽位期间若仍在跑，空闲显存和吞吐都会受影响（脚本会等到至少 11.5 GiB 空闲）。
4. **重叠说明**：Streptophyta 模型的训练集（仅序列）含面板 10 种中的 Juncus、Akebia、Hakea、Amborella，测试集含拟南芥；无监督，无标签泄漏，但对比时要写明。DEV 的拟南芥是模型的作者测试物种（微调后 locus F1 0.799）。
5. 面板 10 种不在本实验范围（随 12-01 的 P6 基线）；面板评分器将复用本评分器，规则（CDS ≥ 96 bp）需与 Helixer / Tiberius 的面板评分统一。
