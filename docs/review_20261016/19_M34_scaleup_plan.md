# 19 M34-NTV2-SCALEUP：多植物物种放大训练的数据准备与预注册协议草案（2026-10-20）

范围：为"NT-v2-500M（纯自监督）+ M31 B 配方（q/v LoRA）+ M28 完整 CDS 语法"从拟南芥/水稻 2 个物种扩大到 10 个物种训练（目标 DEV exact-chain F1 ≥ 70），并检验"注释量—性能曲线"，做**数据准备与方案设计**。本任务**未使用 GPU、未训练、未读取任何 DEV/sealed 序列**；Baobab 未提交计算、未使用 Modal；ont 上只用 CPU（核 24–31，nice，≤8 个并发进程），没有碰 M32/M33 的进程与文件，也没有改动任何冻结代码。

代码：ont `~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/{scripts,tests}`，本地副本 `/Users/jiwang/Desktop/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/`（含 `results/` 小结果；只新增，未 commit）。数据：ont `~/coding-rna/data/m34_train/{raw,prepared,catalog,homology}`；调度：ont `~/coding-rna/outputs/M34-NTV2-SCALEUP/schedules/`。

## 0. 结论先行

| 项 | 结果 |
|---|---|
| 新增训练物种 | **8 个**（真双子叶 5：杨 *Populus trichocarpa*、番茄 *Solanum lycopersicum*、菜豆 *Phaseolus vulgaris*、甜橙 *Citrus sinensis*、胡萝卜 *Daucus carota*；单子叶 3：二穗短柄草 *Brachypodium distachyon*（禾本科 Pooideae）、椰枣 *Phoenix dactylifera*（棕榈目）、几内亚薯蓣 *Dioscorea rotundata*（薯蓣目））；加上已有的拟南芥/水稻共 10 个物种。全部 RefSeq（NCBI EGAP/Gnomon）染色体级装配，基因组 271–833 Mb，注释蛋白编码基因 2.4–3.3 万，注释所用 RNA-seq 比对 12–97 亿 reads |
| 泄漏审计 | 8 个物种的属、accession 与独立面板（主 12 + 备选 15 + Lactuca + 项目 sealed 的 Setaria）**无一重叠**（`panel_exclusion_audit.json`，测试已写入）。科级重叠两处并已说明取舍：禾本科（Brachypodium 与已有的水稻 vs 面板 Zizania/Themeda/Paspalum）、豆科（Phaseolus vs 备选 Medicago sativa） |
| 数据 | 下载 **1.5 GB**（8 物种 FASTA+GFF3+蛋白，另 2 个旧训练物种蛋白、10 个面板蛋白用于同源审计），8 物种×3 个基因组文件 + 2 个旧物种蛋白 + 10 个面板蛋白的 md5 全部与 NCBI 校验文件一致，下载日期 2026-10-20 13:01–13:13 CEST。整个 `m34_train` 目录含全部派生数据 7.4 GB（< 10 GB） |
| 预处理产物 | 每个新物种：仅主核染色体的 `genome.fa`/`reference.gff3`/`lengths.json`/`splits.json`（m1_screen 同格式）；10 个物种全部生成 M28 格式的 `transcripts.jsonl`、`windows_train/val.jsonl`、`summary.json`、嵌套 25%/50% 标注基因集；7 个调度（arm）的确定性抽样序列 |
| 监督量 | 训练可用（eligible）完整基因：拟南芥+水稻 **16,609** → 10 物种 **203,737**（12.3×）。新物种每个留出 8 Mb 内部验证区段（565–955 条 CDS-assessable 参考，合计 5,793 条），拟南芥/水稻整条留出染色体沿用 M31（4,293 + 1,601 = 5,894 条，与 M31 逐条一致，见 §2.4） |
| 协议草案 | 7 个 arm、每个 3,072 步（batch 8，M31 B 配方），物种数曲线 2/4/6/10 + 部分标签曲线（新物种 25%/50% 基因，边缘化 U 对填背景 F）；**部分标签 arm 与全标签 arm 的标注基因总数匹配**（k25U ≈ T4：63.4k vs 62.5k；k50U ≈ T6：110.2k vs 110.8k），可直接比较"多物种浅标注 vs 少物种深标注"。DEV 每个 arm 只评一次（共 7 次，T10 为唯一的 ≥70 判据臂） |
| 时间（ont 单卡） | 每个 fit ≈ **22.5–24 h**（3,072 步 × 23.4 s + 验证 + DEV；乐观 ≈17 h，悲观 ≈33 h）；阶段一（T2 对照 + T10）≈ 47 h；全部 7 个 ≈ 166 h ≈ 6.9 天；排在 M32/M33 之后 |
| 不能保证的 | **DEV ≥ 70 不是承诺**。M31 B 是 63.61（1,536 步、仍在上升）；比 GENERanno fresh-joint 的 70.107 低 6.5。放大物种与步数能否补上这 6.5 点未知；本文给出的可能区间只是猜测（§3.9）。另外 DEV 仍是拟南芥/水稻留出染色体，不是独立物种检验 |

## 1. 训练物种选择

### 1.1 准则（先于任何结果固定）

1. 不在独立评测面板（`10_independent_plant_panel.md` 主面板 12 个 + 备选 15 个 + 被排除的 Lactuca）也不是项目 sealed 的 Setaria；**属级不同、accession 不同**。
2. RefSeq 当前版本、染色体级装配（Chromosome 级或 Complete Genome 级），避免 GenBank 提交者注释（质量与方法不一致）。
3. 注释证据充分：RefSeq EGAP（Gnomon + RNA-seq/蛋白/同物种 mRNA 比对），注释的 BUSCO 完整度 ≥ 92%，蛋白编码基因 ≥ 2 万。
4. 覆盖真双子叶多个目（蔷薇类 fabids/malvids、菊类 lamiids/campanulids）与单子叶（禾本科 Pooideae 亚科、棕榈目、薯蓣目）；基因组 ≤ 1 Gb（实际 271–833 Mb）。
5. 每个物种的 ≥ 90% 以上训练碱基都能进入采样（见 §2.2 的主染色体规则）。

### 1.2 物种表

物种按**固定的嵌套加入顺序**排列（T4 = 前 2 个，T6 = 前 4 个，T10 = 全部），顺序在看任何结果前按"真双子叶/单子叶交替、目不重复"定下。

| 序 | 物种（俗名） | 目 / 科 | accession（装配名） | 装配日期 | 大小 / 染色体数 | 蛋白编码 / 基因总数 | 注释 BUSCO（谱系） | 注释版本；RNA-seq 比对证据（NCBI 注释报告） |
|---|---|---|---|---|---|---|---|---|
| 1 | *Populus trichocarpa*（毛果杨） | Malpighiales / Salicaceae | `GCF_000002775.5`（P.trichocarpa_v4.1） | 2022-10-20 | 392 Mb / 19 | 29,617 / 34,621 | 98.8（eudicots_odb10） | Release 102（2022-12-02）；76 亿 reads，23.5 万 introns，同物种 mRNA/EST 比对 |
| 2 | *Brachypodium distachyon*（二穗短柄草） | Poales / Poaceae（Pooideae） | `GCF_000005505.3`（Brachypodium_distachyon_v3.0） | 2018-01-24 | 271 Mb / 5 | 25,535 / 31,335 | 99.3（poales_odb10） | Release 103（2018-03-16）；97 亿 reads，21.2 万 introns，同物种 mRNA/EST |
| 3 | *Solanum lycopersicum*（番茄） | Solanales / Solanaceae | `GCF_036512215.1`（SLM_r2.1） | 2024-09-03 | 833 Mb / 12 | 28,611 / 41,476 | 98.7（solanales_odb10） | RS_2024_09（2024-09-20）；92 亿 reads，16.6 万 introns，同物种 mRNA/EST |
| 4 | *Phoenix dactylifera*（椰枣） | Arecales / Arecaceae | `GCF_054393675.1`（Male_Medina_1_primary_PhDact_Kaust_v1.1） | 2026-01-07 | 690 Mb / 18 | 25,830 / 44,467 | 95.1（liliopsida_odb10） | RS_2026_09（2026-09-17）；69 亿 reads，21.5 万 introns；非编码基因占 33% |
| 5 | *Phaseolus vulgaris*（菜豆） | Fabales / Fabaceae | `GCF_000499845.2`（P._vulgaris_v2.0） | 2024-04-23 | 537 Mb / 11 | 29,679 / 34,165 | 98.4（fabales_odb10） | RS_2024_09（2024-09-16）；74 亿 reads，15.6 万 introns |
| 6 | *Citrus sinensis*（甜橙） | Sapindales / Rutaceae | `GCF_022201045.2`（DVS_A1.0） | 2022-02-11 | 299 Mb / 9 | 23,556 / 29,331 | 99.2（eudicots_odb10） | Release 103（2022-12-05）；78 亿 reads，17.3 万 introns |
| 7 | *Daucus carota*（胡萝卜） | Apiales / Apiaceae | `GCF_001625215.2`（DH1_v3.0） | 2023-10-20 | 441 Mb / 9 | 33,079 / 39,086 | 99.1（eudicots_odb10） | RS_2024_03（2024-03-08）；83 亿 reads，20.1 万 introns |
| 8 | *Dioscorea rotundata*（几内亚薯蓣） | Dioscoreales / Dioscoreaceae | `GCF_009730915.1`（TDr96_F1_v2_PseudoChromosome…） | 2019-11-25 | 584 Mb / 28 个伪染色体（20 条 ≥ 8 Mb） | 24,562 / 35,061 | **92.4**（liliopsida_odb10） | Release 100（2021-01-28）；**仅 11.6 亿 reads**，无同物种 mRNA/EST；标签质量观察对象（§5） |
| — | *Arabidopsis thaliana*（已有） | Brassicales / Brassicaceae | `GCF_000001735.4`（TAIR10.1） | — | — | 27,562 / 38,312 | — | Araport 注释（M28/M31 沿用） |
| — | *Oryza sativa*（已有） | Poales / Poaceae（Oryzoideae） | `GCF_034140825.1`（AGIS1.0） | — | — | 29,427 / 39,387 | 99.3 | RS_2024_06（M28/M31 沿用） |

来源：NCBI datasets API `dataset_report`（2026-10-20 查询）与各注释报告页（RNA-seq 汇总行）；BUSCO 为 NCBI 对注释蛋白的自报值（未重测）。

### 1.3 与面板的系统发育关系与取舍

分歧时间是文献量级的**近似值（未用 TimeTree 逐对核实）**；面板信息引自 10 号文档。

| 训练物种 | 与最近面板物种/类群的关系 | 取舍 |
|---|---|---|
| Populus（杨柳科，Malpighiales，fabids） | 距 *Dryas*（蔷薇科 Rosales，fabids）≈100–110 Mya；距 Akebia/Eschscholzia/Coptis（基部真双子叶）≈125–135 Mya | 最近的面板物种在目级以上；Dryas 是 M30 已用的零样本物种，距离大于拟南芥—Dryas（≈100–110 Mya 同量级），不构成"近缘泄漏" |
| Phaseolus（豆科，Fabales，fabids） | 距 Dryas ≈100 Mya；**与备选面板 *Medicago sativa* 同科**（菜豆族 vs 车轴草族，≈50–55 Mya） | 备选面板里的 Medicago 本来已被标为"同属泄漏"（M. truncatula 在 Carbon-A/Helixer 训练内），只作敏感性分析。若将来 Medicago 升入主面板，需在报告里标"科级近缘"，或把 Phaseolus 换成备选 *Cucumis sativus*（葫芦科，无面板同科） |
| Citrus（芸香科，Sapindales，malvids） | 与面板无同科；距拟南芥（malvids）≈100–110 Mya | 无泄漏问题 |
| Solanum（茄科，Solanales，lamiids） | 面板无菊类（asterids）物种 | 补上面板没有覆盖的 asterids 方向；对面板零样本是"远类群训练" |
| Daucus（伞形科，Apiales，campanulids） | 同上 | 同上 |
| Brachypodium（禾本科 Pooideae） | **与面板 Zizania（稻亚科，≈50 Mya）、Themeda（黍亚科 Andropogoneae，≈50 Mya）、备选 Paspalum（黍亚科）同科**；已有的水稻与 Zizania 同族（≈20–30 Mya） | 水稻已经让面板中的禾本科物种成为"近训练层"（10 号文档本就把 Zizania 标为"近开发物种"）。Brachypodium 来自第三个亚科，与 Themeda 的距离不小于水稻—Themeda（均 ≈50 Mya），因此**没有使任何面板物种变得更近**；保留它是为了"禾本科内部多样性"。禾本科面板物种在结论里必须与非禾本科分开报告 |
| Phoenix（棕榈科，Arecales） | 距面板 Acorus（石菖蒲目）≈130–140 Mya；距备选 Canna（姜目）≈110–120 Mya；距 Juncus/Rhynchospora（莎草目）≈100–110 Mya | 无同科 |
| Dioscorea（薯蓣科，Dioscoreales） | 距 Acorus ≈130–140 Mya；距备选 Stemona（露兜树目）≈110 Mya | 无同科；注释证据最弱（见下） |

**未纳入的类群及原因**：蕨类（面板 Dipteris 为蕨，无 ≤1 Gb 的 RefSeq 蕨类训练候选，且 Azolla/Adiantum 等在备选面板）、裸子植物（> 10 Gb 或在面板内）、苔藓/绿藻（备选面板含苔藓；绿藻与植物基因结构差异过大）、基部被子植物（Amborella 在面板，Tasmannia 只有 contig 级装配）。面板里的蕨类、基部被子植物、苔藓因此是"完全外推"层，写入协议的解读。

**替换规则（事前固定，防止事后挑选）**：若 §2 的审计或预检发现某物种不可用（如标签质量不达预设门槛），**按固定顺序用下列备选替换，不得改动其他物种**：真双子叶 → *Cucumis sativus*（`GCF_000004075.3`，225 Mb，注释蛋白编码 20,040）→ *Vitis vinifera*（`GCF_030704535.1`，495 Mb，25,187）→ *Sesamum indicum*（`GCF_026168435.1`，309 Mb，23,943）；单子叶 → *Ananas comosus*（`GCF_001540865.1`，382 Mb，22,251，BUSCO 93.9）→ *Asparagus officinalis*（`GCF_001876935.1`，1.19 Gb，26,460）。这些备选的 HEAD 已核实可下载（FASTA 71–360 MB），但**没有下载**。

## 2. 下载与预处理记录

### 2.1 下载清单（实际；NCBI FTP，ont，下载日期 2026-10-20）

| 项目 | 内容 | 大小 |
|---|---|---|
| 8 个新训练物种 | `*_genomic.fna.gz`、`*_genomic.gff.gz`、`*_protein.faa.gz`、`*_assembly_report.txt`、`md5checksums.txt` | FASTA 1,236 MB + GFF3 96 MB + 蛋白 78 MB（FASTA：杨 122、短柄草 87、番茄 252、椰枣 206、菜豆 165、橙 93、胡萝卜 129、薯蓣 184 MB） |
| 2 个旧训练物种蛋白 | 拟南芥 TAIR10.1、水稻 ASM3414082v1 的 `protein.faa.gz`（基因组/GFF 已在 `data/m1_screen`） | 22 MB |
| 10 个面板物种蛋白 | 主面板（除 M30 已有的 Dryas/Eschscholzia 外）的 `protein.faa.gz`，仅供同源审计 | 84 MB |
| 合计 | `raw/` 目录 | **1.5 GB**（限额 10 GB） |

下载前先用 HEAD 核实了全部 URL 与大小（200，见上），然后下载；全部在约 1.5 分钟内完成。校验：每个 `genomic.fna.gz`/`genomic.gff.gz`/`protein.faa.gz` 的本地 md5 与 NCBI `md5checksums.txt` 逐一比对，**24/24（8 物种×3）+ 2 + 10 全部 OK**。完整 md5 与下载时间在 ont `~/coding-rna/data/m34_train/raw/<assembly>/{local_md5.txt,md5checksums.txt,download_date.txt}`，汇总见本地 `results/download_manifest.tsv`。**GFF3 文件 md5 即注释版本冻结**（前 10 位）：

| 物种 | genomic.fna.gz | genomic.gff.gz | protein.faa.gz |
|---|---|---|---|
| Populus | c76a9756ea | 8de0782069 | 2f577263da |
| Brachypodium | cb4419f7cb | a02bf0b556 | dd3ebcd2a4 |
| Solanum | fbefc524fd | 38ace823cc | 2898310c86 |
| Phoenix | ed01ec76a8 | f231db2d5d | 31bef0c057 |
| Phaseolus | 3229368b33 | 6535df68d4 | 71aca82947 |
| Citrus | d514a81ba2 | fad032a3f7 | 5237742162 |
| Daucus | 88e9df3209 | 5bdfc49cbf | b9348bc240 |
| Dioscorea | abbab0225e | 55f6e7e81a | 2bd1a9362c |

### 2.2 预处理规则（`m34_prep.py`；先于任何模型输出固定）

1. **主核染色体**：装配报告中 `assembled-molecule` + 位置类型 `Chromosome` + `Primary Assembly` + 长度 ≥ 5 Mb。叶绿体/线粒体/质体、未定位/未锚定 scaffold 全部丢弃；Dioscorea 另有 8 条 < 5 Mb 的伪染色体（0.1–4.1 Mb，同一装配的未锚定部分）被长度门槛排除。
2. **内部验证区段**（每个新物种）：长度**下中位**（偶数条取较短的中位）的那条染色体，只保留其**前 8 Mb**（8 Mb 内完整的基因，坐标不变，伪 contig 名 `<染色体>__val8M`）；该染色体其余部分**弃用**（既不训练也不评估，同 M30 §2）。规则只依赖长度，不依赖基因密度或模型表现。弃用碱基合计 265.7 Mb（占 8 个新物种主染色体总长约 4.0 Gb 的 6.6%）。
3. 其余主染色体全部作为**训练染色体**。

| 物种 | 训练染色体（条 / Mb） | 内部验证染色体（长度）→ 8 Mb 区段 | 弃用 Mb |
|---|---|---|---|
| Populus | 18 / 371.4 | NC_037298.2（17.8 Mb） | 9.8 |
| Brachypodium | 4 / 211.9 | NC_016132.3（59.1 Mb） | 51.1 |
| Solanum | 11 / 764.0 | NC_090802.1（68.7 Mb） | 60.7 |
| Phoenix | 17 / 655.8 | NC_142554.1（34.4 Mb） | 26.4 |
| Phaseolus | 10 / 465.9 | NC_023756.2（48.0 Mb） | 40.0 |
| Citrus | 8 / 268.4 | NC_068563.1（30.6 Mb） | 22.6 |
| Daucus | 8 / 393.2 | NC_030385.2（47.6 Mb） | 39.6 |
| Dioscorea | 19 / 468.5 | NC_052485.1（23.4 Mb；20 条 ≥ 8 Mb 的下中位） | 15.4 |
| 拟南芥（沿用 M31） | 2 / 57.4（NC_003070.9、NC_003076.8） | NC_003071.7（19.7 Mb，整条） | 0 |
| 水稻（沿用 M31） | 2 / 55.7（NC_089043.1、NC_089046.1） | NC_089044.1（23.9 Mb，整条） | 0 |

拟南芥/水稻另有 **DEV 染色体 NC_003074.8、NC_089041.1**（M31 DEV，每个 arm 只评一次）和 **NC_003075.7、NC_089042.1**（项目 `split_test` 染色体，M34 不读取）；调度器的断言保证它们和所有 `__val` 伪 contig 永远不会被抽到。

### 2.3 产物（ont `~/coding-rna/data/m34_train/`）

- `prepared/<物种>/{genome.fa, reference.gff3, lengths.json, splits.json}`：新物种，m1_screen 同格式；`prepared/panel_exclusion_audit.json`。
- `catalog/<物种>/`（10 个物种全有）：`transcripts.jsonl`（主转录本目录：CDS_assessable_primary、CDS_complete_supported、eligible、overlapping_primary、序列结构失败原因、剪接基序）、`windows_train.jsonl`/`windows_val.jsonl`（M28 窗口行：24,576 bp，步长 12,288，两条链，`positive_ids`、未知区间、`no_annotated_gene_either_strand`、`motifs`、`label_policy=CDS_assessable_primary_R4`）、`summary.json`、`labeled_k25.json`/`labeled_k50.json`（嵌套标注基因集，seed 20261017，同 M30 协议 §3）。新物种窗口数：训练 34,490（短柄草）–124,344（番茄），验证 1,302（8 Mb 区段）；拟南芥/水稻训练 9,342 / 9,070。
- `homology/`：MMseqs2 审计输入与结果（0.5 GB；`tmp/` 已删）。
- 调度：`outputs/M34-NTV2-SCALEUP/schedules/schedule_<arm>.jsonl` + `schedule_summary.json`。
- 标签构建器 `m34_labels.Labeler` 与 M30 `Labeler` 同逻辑（M30 测试断言其与冻结 R4 标签在拟南芥窗口上 > 99.9% 一致，唯一已知差异是备选转录本冲突掩码；本任务另验证了目录计数，见 §2.4），支持处理 `FULL`/`U`（未选中基因位点标 unknown，边缘化）/`F`（填背景）。**不需要 GPU 编码缓存**：NT-v2 + LoRA 的编码器在训练中变化，特征不可预取（M31 B 同样如此）；如果以后加入冻结特征的 arm，其缓存体量与 M28 的 `cache_features` 同数量级，**此处只写方案、不执行**。

### 2.4 监督量与一致性检查

训练染色体上的主转录本统计（`catalog/*/summary.json`，汇总 `results/species_supervision.tsv`）：

| 物种 | 训练 Mb | 主编码基因 | CDS-assessable 参考 | **eligible（可作完整正例）** | eligible/Mb | 外显子 1 / 2–5 / ≥6 | CDS 跨度 > 24,576 bp | 非标准剪接占比 | 验证区段参考（assessable / eligible） |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|
| 拟南芥 | 57.4 | 13,487 | 13,446 | 13,325 | 232.1 | 3,026 / 5,856 / 4,564 | 2 | 0.5% | 4,293 / 4,269 |
| 水稻 | 55.7 | 3,377 | 3,348 | 3,284 | 58.9 | 853 / 1,479 / 1,016 | 0 | 1.7% | 1,601 / 1,561 |
| Populus | 371.4 | 27,933 | 27,745 | 26,842 | 72.3 | 5,101 / 12,069 / 10,575 | 192 | 2.3% | 821 / 795 |
| Brachypodium | 211.9 | 19,679 | 19,548 | 19,086 | 90.1 | 4,319 / 8,692 / 6,537 | 20 | 2.0% | 955 / 928 |
| Solanum | 764.0 | 25,589 | 25,473 | 25,126 | 32.9 | 5,395 / 11,681 / 8,397 | 284 | 1.0% | 570 / 562 |
| Phoenix | 655.8 | 23,901 | 23,787 | 23,116 | 35.2 | 4,032 / 10,571 / 9,184 | **1,134** | 2.0% | 630 / 619 |
| Phaseolus | 465.9 | 26,808 | 26,751 | 26,268 | 56.4 | 6,605 / 11,248 / 8,898 | 176 | 1.4% | 565 / 556 |
| Citrus | 268.4 | 21,248 | 21,104 | 20,481 | 76.3 | 4,308 / 9,425 / 7,371 | 153 | 1.7% | 832 / 815 |
| Daucus | 393.2 | 29,170 | 29,038 | 28,403 | 72.2 | 6,157 / 13,118 / 9,763 | 121 | 1.7% | 727 / 716 |
| Dioscorea | 468.5 | 20,349 | 20,067 | 17,806 | 38.0 | 3,567 / 9,233 / 7,267 | 115 | **10.8%** | 693 / 608 |
| **合计** | 3,712 | | | **203,737**（T2 = 16,609） | | | | | 新物种验证区段 5,793 条参考 |

一致性检查（证明新构建器与 M28/M31 的数据口径一致）：
- 拟南芥训练染色体 CDS-assessable 参考 13,446 = M29 `EXPECTED_PRIMARY`（7,135 + 6,311）；水稻 3,348 = 1,638 + 1,710；**两个内部验证染色体 4,293 与 1,601 与 M31 完全一致**（合计 5,894）。
- 测试 `tests/test_m34.py`、`tests/test_m34_labels.py`：主染色体规则、下中位规则、嵌套标注集、与 M30 同规则的基因窗口、面板排除审计、`splits.json` 一致性、U/F/FULL 语义（U 与 F 只差"背景→unknown"，FULL 与 k=1 逐位相等）、调度不含 DEV/验证/未用 contig——全部通过（ont：11 项通过，0 跳过，78 s）。

关注点：
- **Dioscorea 的非标准剪接占 10.8%**（其余 0.5–2.3%），注释证据最弱（无同物种 mRNA/EST，仅 11.6 亿 reads）。非标准剪接的基因被标为非 eligible（掩码为 unknown，不会作为错误正例），但这个比例本身提示该注释含较多低置信模型；列为"标签质量观察对象"（§5，设有事前替换门槛）。
- **Phoenix 有 1,134 个（4.8%）CDS 跨度 > 24,576 bp 的基因**，番茄 284、杨 192 个；这些基因无法完整落入一个 24 kb 核心窗口（M28 的上下文最长 49 kb），训练中只能作为部分信息，评估时则是难例。Phoenix 的"非编码基因 14,804（占 33%）"也值得注意：这些位点在我们的标签里是背景，其中可能含未被识别的编码基因。

### 2.5 面板排除审计与蛋白同源审计

- 排除审计（`m34_prep.py` + `panel_exclusion.py`）：10 个训练物种与 29 条面板/排除/sealed 记录比对，**属级与 accession 无违规**；科级重叠已列于 §1.3。
- 同源审计（`m34_homology_audit.py`，MMseqs2 `-s 6.0 -e 1e-5`，query 覆盖度 ≥ 0.7）：面板各物种蛋白对各嵌套 arm 训练蛋白的最佳命中一致度 ≥ 50/70/90% 的比例。结果：

每个面板物种随机抽 1,500 个蛋白（含 isoform，seed 0），对各嵌套 arm 的训练蛋白库（T2 = 拟南芥+水稻 9.1 万条；T10 = 全部 10 个物种 44.6 万条）做 MMseqs2（`-s 4.0 -e 1e-5`，query 覆盖度 ≥ 0.7），取最佳命中；表中数字 = 最佳命中一致度 ≥ 阈值的蛋白百分比（%）。

| 面板物种 | ≥50%：T2 → T4 → T6 → T10 | ≥70%：T2 → T4 → T6 → T10 | ≥90%：T2 → T10 |
|---|---|---|---|
| Dryas（RefSeq TH） | 69 → 81 → 83 → 85 | 28 → 48 → 49 → 53 | 1 → 4 |
| Eschscholzia（TH） | 61 → 70 → 73 → 74 | 24 → 31 → 33 → 35 | 2 → 4 |
| Akebia（TH） | 63 → 73 → 76 → 77 | 26 → 36 → 39 → 41 | 2 → 4 |
| Hakea（TH） | 68 → 77 → 80 → 81 | 27 → 38 → 41 → 43 | 2 → 4 |
| Chenopodium（TH） | 54 → 61 → 63 → 66 | 21 → 27 → 28 → 30 | 1 → 2 |
| Dipteris（TH，蕨） | 36 → 38 → 40 → 41 | 10 → 11 → 12 → 12 | 1 → 1 |
| Zizania（SH） | 84 → 85 → 85 → 85 | 74 → 74 → 74 → 74 | 24 → 24 |
| Themeda（SH） | 42 → 44 → 44 → 44 | 28 → 29 → 29 → 29 | 4 → 5 |
| Juncus（SH） | 65 → 68 → 72 → 73 | 27 → 30 → 33 → 34 | 3 → 5 |
| Rhynchospora（SH） | 62 → 65 → 68 → 69 | 27 → 29 → 32 → 33 | 5 → 7 |
| Acorus（SH） | 56 → 61 → 66 → 68 | 22 → 26 → 32 → 33 | 2 → 4 |
| Amborella（AH） | 63 → 69 → 72 → 73 | 22 → 28 → 31 → 33 | 2 → 3 |
| **12 个面板物种均值** | 60.3 → 66.1 → 68.5 → 69.8 | 28.0 → 33.7 → 36.1 → 37.6 | 4.2 → 5.5 |

解读：
1. 加入 8 个新物种使面板"拉近"：均值上，一致度 ≥70% 的蛋白占比 28.0% → 37.6%，蔷薇科 Dryas 最明显（28 → 53%，因为 Populus 同属 fabids）。这是**训练集更靠近面板**，不是面板物种进入训练；它是"更多物种带来迁移"的合理机制，但意味着 T10 与 T2 在面板上的差异混有"训练集与面板的蛋白距离缩短"，不能全部解释为"更好的泛化"。**因此面板结果必须按每个面板基因的最佳命中一致度分层（< 50% / 50–70% / ≥ 70%）报告**，并把距离最远的物种（Dipteris 12%、Chenopodium 30%、Themeda 29%）单独列出。
2. Zizania 在只有拟南芥+水稻时已达 74%（水稻同族），之后不再变化，符合 10 号文档把它标为"近开发物种"的定位；Brachypodium 几乎没有拉近 Themeda（28 → 29%），支持 §1.3 关于"禾本科取舍不增加面板泄漏"的判断。
3. 局限：每物种仅 1,500 个蛋白，`-s 4.0` 灵敏度对低一致度偏保守；一致度为比对区内一致度；本审计纯描述，没有据此过滤任何数据。

## 3. 预注册协议草案（M34；执行期间不得根据结果修改，修改须在 §3.12 追加记录）

### 3.1 要回答的问题与判据

- **P1（性能目标）**：T10 最终模型的 DEV exact-chain F1（拟南芥 NC_003074.8 + 水稻 NC_089041.1，7,728 条主参考，项目 ruler，CDS ≥ 96 bp 规则）**≥ 70.0**，即追平 GENERanno fresh-joint（70.107）。≥ 66 记为"实质进展"，< 63.61（M31 B）记为"无增益"。1 Mb 块自助法（2,000 次）报告 CI。
- **P2（数据效应，排除"只是训练更久"）**：T10 − T2 的 DEV 配对差（同一自助法重采样），判据 ≥ +2.0 且 CI 下界 > 0。T2 与 T10 步数、批、sampler、初始化完全相同，只差训练物种。
- **P3（注释量—性能曲线）**：Ath/Osa 内部验证 pooled F1 对 log10(训练标注基因数) 的曲线，两条轴合并拟合；趋势性结论，不设通过/失败门槛。
- **P4（部分标签，C2）**：k50U − k50F 在**新物种内部验证区段 pooled** 上 ≥ +2.0 点，且在 8 个新物种中 ≥ 6 个同向；"保留的全标签增益"= (k50U − T2)/(T10 − T2) 作为描述量。

### 3.2 数据划分

| 集合 | 内容 | 用途 |
|---|---|---|
| 训练 | §2.2 的训练染色体（10 物种共 3.71 Gb）；Ath/Osa 沿用 M31 的 4 条 | 抽样 |
| 内部验证 | Ath NC_003071.7、Osa NC_089044.1（整条）+ 8 个新物种各 8 Mb 区段 | checkpoint 选择（仅 Ath/Osa 两条，所有 arm 同一集合，防止物种多的 arm 在选择上占便宜）；新物种区段只作报告 |
| DEV | NC_003074.8 + NC_089041.1 | 每个最终 checkpoint 一次 |
| 弃用/未触碰 | 8 个新物种验证染色体的其余部分、NC_003075.7、NC_089042.1、sealed、Setaria | 不读取 |

### 3.3 模型与配方（沿用 M31 B，逐项）

NT-v2-500M-multi-species，FP32 权重 + BF16 autocast，融合 SDPA；特征 hidden_states **21 + 27** 拼接（2,048 维）；窗口 12,288 bp，步长 6,144；头 `M28Core('C0', 2048, 128)`（366,191 参数）；**从 M31 runA `model_002048.pt`（头，冻结编码器训练 2,048 步）初始化头，LoRA 零初始化**；LoRA：q/value，rank 8、alpha 16，块 13–27（0 起 12–26，491,520 参数，FP32）；AdamW（wd 0.01，β=(0.9,0.999)），头 lr 峰值 1e-4、LoRA lr 峰值 3e-4，cosine 到 1e-5，5% 线性 warmup，梯度裁剪 1.0；batch 8；损失为配分函数 CRF（精细标签硬约束 + ε=0.1 粗通道 + 未知位置边缘化 + 128 外显子）；DP 并行 12 进程；每 64 步存可恢复检查点；上下文几何同 M31（24,576 核心 + 两侧各一个 12,288 邻窗）。所有 arm 同一初始化——这意味着所有 arm 共享"A 阶段在 Ath/Osa 上训练头"的先验，在曲线中是常数偏移，并记入局限。

### 3.4 步数与 sampler

- **3,072 步**（= 24,576 次抽样，M31 B 的 2×）。依据：M31 B 的内部验证 56.1 → 61.7 → 63.3（512/1,024/1,536 步）在退火完成时仍在上升；放大数据后需要更长的训练。不设早停（cosine 已定义到 3,072），checkpoint 每 768 步保存。
- **sampler（一条规则用于所有 arm，seed 1）**：每次抽样先按 arm 的物种权重选物种（**等权重**，M28 的"物种均等"惯例）；然后以 50% 概率在该物种全部（训练 contig、链、24,576-bp 网格窗口）中均匀选窗口，以 50% 概率"基因优先"（均匀选一个已标注 eligible 基因，再在覆盖它的 ≤5 个窗口中均匀选，取该基因的链）。这是 M30 协议的 sampler；与 M31 B 的"层内均匀"不同（基因优先让稀疏基因组的有效监督更高，并使部分标签 arm 的批内有标注基因数不被 k 稀释）。**偏离已声明**：T2 控制臂用同一 sampler，所以 T2 vs M31 B 的差异混有 sampler 因素，T10 vs T2 则不混。
- **物种权重取舍**：基因数比例会让 Ath/Osa 只占 8.2%、新物种各占 8–14%；平方根比例 12.4%；等权重 20%。选等权重（简单、与 M28 一致、避免大基因组/基因多的物种主导）。风险是每个物种只分到约 2,458 次抽样：新物种各自约 2–7% 的窗口会被抽到，而拟南芥/水稻约 26%。**条件分支**：若 T10 的 Ath/Osa 内部验证比 T2 低 ≥ 1.0 点（稀释），追加一个"锚定"arm（Ath/Osa 合计 40%，其余 8 物种均分 60%），事前写明、不事后挑。

### 3.5 arm 表

训练标注基因数 = 各训练染色体上 eligible 基因，部分标签 arm 的新物种取嵌套的 25%/50%（拟南芥/水稻始终 100%）。

| arm | 训练物种 | 新物种标注 | 未选基因处理 | 标注基因总数 | 目的 |
|---|---|---|---|---:|---|
| **T2** | Ath, Osa | — | — | 16,609 | **算力匹配对照**（同步数/同 sampler，只有 2 个物种） |
| T4 | + Populus, Brachypodium | 100% | — | 62,537 | 物种曲线 |
| T6 | + Solanum, Phoenix | 100% | — | 110,779 | 物种曲线 |
| **T10** | + Phaseolus, Citrus, Daucus, Dioscorea | 100% | — | **203,737** | **P1/P2 主臂** |
| T10_k50U | 同 T10 | 50% | U（边缘化） | 110,174 | C2；标注基因数 ≈ T6 |
| T10_k25U | 同 T10 | 25% | U | 63,394 | C2；标注基因数 ≈ T4 |
| T10_k50F | 同 T10 | 50% | F（填背景） | 110,174 | C2 对照（边缘化 vs 填背景） |

标注基因数匹配的配对：**k25U vs T4**（63.4k vs 62.5k，10 物种浅标注 vs 4 物种深标注）、**k50U vs T6**（110.2k vs 110.8k）。若"多物种浅标注"优于"少物种深标注"，则**物种多样性本身有价值**；反之则只是总基因数起作用。这是"需要多少注释监督"这个问题里最干净的一对比较。

调度实际暴露量（`results/arms_exposure.tsv`）：| arm | 标注基因总数 | 24,576 次抽样的唯一窗口 | 含 ≥1 个完整标注基因的抽样占比 | 唯一的"被完整看到"标注基因（占总数） | 完整基因暴露总次数（次/基因） |
|---|---:|---:|---:|---:|---:|
| T2 | 16,609 | 12,739 | 79.6% | 15,250（91.8%） | 52,375（3.4） |
| T4 | 62,537 | 19,740 | 77.0% | 25,317（40.5%） | 42,087（1.7） |
| T6 | 110,779 | 22,000 | 70.2% | 25,560（23.1%） | 34,892（1.4） |
| T10 | 203,737 | 23,340 | 68.4% | 27,514（13.5%） | 32,100（1.2） |
| T10_k50U / k50F | 110,174 | 23,223 | 59.7% | 20,654（18.7%） | 25,095（1.2） |
| T10_k25U | 63,394 | 22,988 | 54.5% | 17,203（27.1%） | 21,778（1.3） |

（"完整看到"= 基因连同两侧 256 bp 完整落在被抽到的 24,576-bp 核心窗口内；部分落入窗口的基因仍提供局部监督，未计入。）**这张表揭示的设计事实**：在**恒定预算**（3,072 步 × 8）下，T4、T6、T10 被完整看到的唯一标注基因数几乎相同（25.3k / 25.6k / 27.5k），T2 则把 1.5 万个基因平均重复看 3.4 次；也就是说 T4→T10 的差别主要是**物种多样性**，而不是"看到的基因更多"。要检验"更多标注基因"本身，需要更多步数或更大批，这不在恒预算的主设计内，作为阶段四的可选 arm（见 §3.9）。k25U 与 T4 的配对是**标注池大小**匹配（63.4k vs 62.5k），不是暴露量匹配（17.2k vs 25.3k）。

### 3.6 内部验证与 checkpoint 选择

- 每 768 步（768、1,536、2,304、3,072）在 **该 arm 的验证集合**上推理 + 精确 skyline 解码，算 pooled exact-chain F1（沿用 `c1_infer.py`/`c1_val.py`，参考 = 区段内完整的 CDS-assessable 主参考，预测 CDS ≥ 96 bp）。
- **最终 checkpoint = {1,536、2,304、3,072} 中 Ath/Osa 两条验证染色体 pooled F1 最高者**（并列取步数大者）。所有 arm 同规则；选择不使用 DEV、不使用新物种区段。
- 选定后，对每个 arm 再在**全部 8 个新物种区段**上评一次（对没有该物种的 arm 即跨物种迁移，对有该物种的 arm 即同物种留出区段）——这是物种数曲线的另一半证据。

### 3.7 DEV 的使用

每个 arm 的最终 checkpoint 在 DEV 上**评一次**，合计 7 次；加上 M31 的 2 次，DEV 累计 9 次，必须在论文里如实说明。**只有 T10 的 DEV 结果参与 P1 判据**；其他 arm 的 DEV 值是曲线的辅助点，不做多重比较校正。DEV 是同物种留出染色体，不能代表对独立物种的泛化。

### 3.8 次要评估（独立性证据）

- **面板零样本**（单次）：最终 T2 与 T10 在 M30 评估集 B（Dryas S2/S4/S6/S8 102.8 Mb；Eschscholzia SCAF_2/4/6 170.5 Mb，共 273 Mb）和 M33 的时间留出面板上评估，与 M30 的 Carbon-A 原生/闸门结果同尺比较；约 2.4 h/模型。
- 按外显子数分层（1 / 2–5 / ≥6）、CDS 碱基 F1、背景负担，沿用 M31 §5.2 的口径。

### 3.9 预期（猜测，非承诺）

M31 B 在 1,536 步内部验证 63.3，退火后仍升；多物种可能主要改善多外显子链的连接精度（M31 的差距集中在 ≥6 段：53.98 对 GENERanno 65.81）。单靠"更多步数"可能给 T2 带来 +1–3 点；更多物种再带来几点是合理猜想，但没有数据支持"一定能达到 70"。另据 §3.5 的暴露表，恒预算下 T10 看到的唯一标注基因数只比 T4 多 8%，物种多样性之外的增益（更多基因）需要更多步数：可选的**阶段四**为 T10 的 6,144 步版本（≈48 h；仅当 T10 在内部验证上高于 T4 时才做，且须另行预注册）。若 T10 < 66，下一步是**另行预注册**的容量杠杆（更大 rank、更多目标矩阵、更长训练、后验过滤/温度校准——M31 §8 已列），而不是在本协议里临时加。

### 3.10 部分标签（C2）的方案细节

- 训练物种只给 k 比例的 eligible 基因（嵌套，seed 20261017，按（contig, 起点, id）排序后固定排列）；Ath/Osa 保持全标签（C2 关注的是"新物种只有不完整注释"）。
- **U**：未选中 eligible 基因的 CDS 外包区间在细/粗标签里标 unknown（IGNORE），损失边缘化；**F**：当作不存在，位置标背景（角色 0）。同 k 的 U/F 共用同一抽样序列（抽样只依赖物种集合与标注集合，不依赖 U/F）。
- 注意：U 假定"缺失注释的位点已被某种证据标成未知"，这是边缘化的上界设定（M30 同样的限定）；真实的不完整注释更嘈杂。
- 有效监督会随 k 下降：k25 时一个批次内仍有标注基因的窗口比例见 §3.5 的调度统计（基因优先 sampler 使其不至于被 k 稀释）。

### 3.11 统计

单 seed（头初始化 seed 0、抽样 seed 1）；若 T10 通过 P1，则追加 seed 2、3 的 T10（各 ≈24 h）估计 seed 方差。配对自助法：54 个 1 Mb 块式重采样（M31 的 DEV 做法）用于 DEV；新物种区段用按基因重采样。所有 arm 的失败分解照报。

### 3.12 修订记录
（空。修订须写：时间 / 原因 / 是否已看过任何 arm 的验证或 DEV 结果。）

### 3.13 实施清单（本任务**未**实现，需要另一个实施任务）

1. `m34_train_lora.py`：复制 `c1_train_lora.py`，把数据来源换成 `schedule_<arm>.jsonl` + `m34_labels.Labeler`（对 10 个物种统一构建 FULL/U/F 目标），保持 DP 并行/检查点/日志；注册新物种 contig 到 M29 的 `read_sequence`（按物种分派数据目录）。
2. 20 步 smoke（GPU 占用 ≈5 GB，≤5 分钟，只在 M32/M33 空闲时）：验证新窗口（含 ≥128 外显子、跨度 > 24 kb 的基因）不触发 `InfeasibleObservation` 过多、DP 步时与 M31 B 同量级。
3. `c1_val.py`/`c1_infer.py` 对新物种验证伪 contig 的区域推理（`NC_xxx__val8M`）。
4. 每个 arm 的 `c1_stage` 风格监督脚本（含 OOM 自动恢复、单实例哨兵、DEV 一次性哨兵）。

## 4. 时间估算（ont 单卡 RTX 4090，参考 M31 B 的实测）

M31 B 实测：23.4 s/步（范围 18–35；batch 8，12 个 DP 进程，ont CPU 被他人任务压到负载 ≈38），验证推理 B 模型 32.2 kbp/s（DEV 53.4 Mb = 1,660 s），显存峰值 5.15 GB。当前 ont 负载 ≈ 55–62，故以 23.4 s/步为"可能"值。

| 项（S = 3,072 步） | T2 | T4 | T6 | T10 / k50U / k25U / k50F |
|---|---:|---:|---:|---:|
| 训练 3,072 × 23.4 s | 20.0 h | 20.0 h | 20.0 h | 20.0 h |
| 验证 ×4（Ath/Osa 43.6 Mb + 新物种区段） | 4 × 0.38 = 1.5 h | 4 × 0.51 = 2.1 h | 4 × 0.65 = 2.6 h | 4 × 0.93 = 3.7 h |
| 全部 8 个新物种区段一次（迁移） | +0.55 h | +0.41 h | +0.28 h | 0 |
| DEV 一次（53.4 Mb） | 0.46 h | 0.46 h | 0.46 h | 0.46 h |
| **单个 fit 合计** | **22.5 h** | **22.9 h** | **23.3 h** | **24.2 h** |

- 7 个 arm 合计 ≈ **166 h ≈ 6.9 天**。乐观（18 s/步 = M31 B 实测下限、推理 75 kbp/s）：每个 16.5–17.5 h，总计 ≈ 5.0 天；悲观（35 s/步、推理 32 kbp/s）：每个 32–34 h，总计 ≈ 9.7 天。步数缩放：S = 2,048 时训练 13.3 h，S = 4,608 时 30.0 h。
- 建议分阶段：**阶段一**（T2 对照 + T10）≈ 47 h，之后看 P2；**阶段二**（T4、T6）≈ 46 h；**阶段三**（k50U、k50F，最后 k25U）≈ 73 h。阶段间只做一个门：T10 的 Ath/Osa 内部验证若不高于 T2，则暂停阶段二、三并复盘（数据多样性无效时曲线实验的边际价值低；C2 不受此限，可单独推进）。
- 排队：ont 只有一块 GPU，同一时间只跑一个我方 GPU 任务；M32（PlantCAD2）完整训练与 M33（后验校准 + 时间留出零样本）排在前面，M34 预计在它们之后开始；还有他人任务占 2–19 GB 显存，M31 A 训练曾因此 3 次 OOM 重启（自动从 ≤64 步前恢复，最多损失 64 步）。
- CPU 数据准备耗时（已完成，8 个并发进程以内）：下载 1.5 分钟；预处理 8 物种约 6 分钟；标签/窗口 10 物种各 2–6 分钟；调度 7 个 arm 各约 0.5–2 分钟；同源审计约 25 分钟（8 线程、nice 10，在被占满的 CPU 上，query 抽样 1,500 蛋白/物种；首版全量 46 万 query 在预过滤阶段过慢而放弃）。

## 5. 风险

1. **ont 排队与算力**：单卡、他人任务、CPU 负载 50–60；7 个 fit ≈ 7–10 天连续占用，且排在 M32/M33 之后。缓解：分阶段、门控、可恢复检查点；3,072 步方案可按 §4 缩放。DP 的 CPU 时间（≈25 s/步）而非 GPU 是瓶颈，CPU 被占满会直接拉长所有 arm。
2. **数据泄漏**
   - 面板属级/accession 无重叠；**科级**重叠（禾本科、豆科）见 §1.3；深度同源（保守蛋白家族）不可避免，同源审计只描述、不过滤。
   - **注释管线同源**：训练标签（RefSeq Gnomon/EGAP）与面板 TH 层（也是 RefSeq EGAP）出自同一管线，共享注释风格与系统偏差；SH 层（GenBank 提交者注释）才是"注释风格独立"的检验，两层结果必须分开报告。
   - **NT-v2 预训练不含植物**仅基于模型卡/数据卡声明（M21 溯源 scout），未独立审计。
   - 内部验证区段与训练染色体之间的旁系同源/片段重复可使验证偏高（M31 同样存在）。
   - **DEV 不独立**：仍是拟南芥/水稻留出染色体，累计 9 次使用；≥70 的结论必须配面板零样本才能写成泛化。
3. **标签质量**：Gnomon 注释含 ab initio 与证据驱动模型；Dioscorea 证据最弱且非标准剪接 10.8%；Phoenix 非编码基因占 33% 且 4.8% 的基因跨度超窗口；U 处理的"缺失位点已知为基因位点"是上界假设。缓解：非标准剪接/部分/重叠基因掩码为 unknown；**事前门槛**：若 Dioscorea 的 T10 内部验证区段 F1 低于其余新物种中位数 −15 点（说明标签噪声/难度异常），按 §1.3 的替换规则换 *Ananas comosus*（只重做 T10 及其后 arm，T2/T4/T6 不受影响）。
4. **实施风险**：训练器未实现（§3.13）；新物种 contig 需要登记到 M29 辅助函数；多物种数据路径的 bug 会污染所有 arm——靠 20 步 smoke 与 U/F/FULL 单元测试降低。
5. **科学风险**：T10 可能达不到 70（§3.9）；等权重稀释 Ath/Osa（已设锚定条件分支）；3,072 步仍可能欠训练；单 seed；sampler 与 M31 B 不同（§3.4）；arm 共享 A 阶段头初始化。
6. **多重检验**：7 次 DEV、多个次要比较；只有 T10 的 P1 是确认性终点，其余是描述性曲线。

## 6. 复现

```
# ont：cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts   (python = ~/envs/generanno/bin/python；CPU 核 24-31，nice)
bash m34_download.sh                              # 8 物种，md5 校验
bash m34_download_legacy_proteins.sh; bash m34_download_panel_proteins.sh
python m34_prep.py                                # 主染色体 / 验证区段 / prepared/<key>/ / 面板排除审计
python m34_labels.py build <key> ...              # 10 个物种（含 arabidopsis_thaliana、oryza_sativa）
python m34_schedule.py [arm ...]                  # 7 个 arm 的确定性抽样序列
python m34_homology_audit.py 4                    # MMseqs2 面板-训练同源审计
python m34_tables.py ../results                   # 汇总表
python -m pytest ../tests -q
```

## 7. 产物清单

| 内容 | 位置 |
|---|---|
| 脚本 | `experiments/M34-NTV2-SCALEUP/scripts/`：`species.json`、`m34_common.py`、`panel_exclusion.py`、`m34_download*.sh`、`m34_prep.py`、`m34_labels.py`、`m34_schedule.py`、`m34_homology_audit.py`、`m34_tables.py`；`wait_event.sh`、`push.sh` |
| 测试 | `experiments/M34-NTV2-SCALEUP/tests/`（`test_m34.py`、`test_m34_labels.py`） |
| 小结果（本地） | `experiments/M34-NTV2-SCALEUP/results/`：`species_supervision.tsv`、`arms_exposure.tsv`、`schedule_summary.json`、`download_manifest.tsv`、`homology_audit.json`、`prepared/panel_exclusion_audit.json`、`prepared/<物种>/splits.json`、`catalog/<物种>/summary.json` |
| ont 数据 | `~/coding-rna/data/m34_train/{raw,prepared,catalog,homology}` |
| ont 调度 | `~/coding-rna/outputs/M34-NTV2-SCALEUP/schedules/` |
