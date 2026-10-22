# 25 独立时间留出植物面板：数据准备记录（2026-10-21/22，仅 CPU 与网络）

目的：为论文大纲（`docs/paper_outline_20261021.md` §2 R4、§3 P5/P6）的 Fig. 4 / Table 1 把 10 号文档主面板 12 个物种的数据下载、冻结、整理成 M34/m1_screen 同格式，并完成暴露与蛋白同源审计。**本任务没有运行任何模型、没有做任何推理或评分**；面板序列在看任何预测前已冻结（§2 的 md5 与 §3 的序列清单）。代码在本地 `experiments/PANEL-PREP/{scripts,tests,results}`（未 commit），ont 副本 `~/coding-rna-reorg/experiments/PANEL-PREP`。规则来源：10 号文档 §2.2/§3/§6，19 号文档 §0–§2（M34 的下载、md5、主核染色体、`panel_exclusion`、同源审计做法）；准备脚本直接复用 M34 的函数（`m34_common.parse_assembly_report/main_chromosomes`、`m34_prep.read_fasta/write_fasta`、`m34_labels.Labeler/stats`、`panel_exclusion.PANEL/audit`、`m34_homology_audit.faa`），未改动任何 M34/冻结代码。

## 0. 结论先行

| 项 | 结果 |
|---|---|
| 下载 | 12 物种 × (FASTA + GFF3 + 蛋白 + assembly_report + assembly_stats + NCBI `md5checksums.txt` + datasets v2 `dataset_report`)，ont `~/coding-rna/data/panel/raw/<accession>_<assembly_name>/`，2026-10-21 23:13–23:20 CEST，**60/60 个文件 md5 与 NCBI 校验文件一致**，原始数据 2.7 GB |
| 版本冻结 | 每个 RefSeq 物种记录 annotation release 名（`RS_YYYY_MM`）、注释日期与 GFF md5；Dryas/Eschscholzia 的 fna/gff/faa md5 与 M30（2026-10-17 下载）逐一相同，10 个面板蛋白文件 md5 与 M34 `panel_proteins`（2026-10-20）相同，即 NCBI 在此期间未重新发布这些文件 |
| 主核序列 | 11 个物种用 M34 规则（assembled-molecule + Chromosome + Primary Assembly + ≥ 5 Mb），Eschscholzia（Scaffold 级、无 assembled-molecule）用 M30 规则（Primary Assembly 中 ≥ 10 Mb 的序列，6 条）；合计 **179 条序列、8.25 Gb**（NCBI 报告总长 8.62 Gb 的 95.7%），丢弃 4,574 条短序列/细胞器/未锚定序列（0.37 Gb，19,640 个注释基因，主要是 Acorus 109.5 Mb 与 Rhynchospora 50.6 Mb 的未锚定 scaffold） |
| 参考统计 | 冻结 ruler（`common_ruler.select_reference`，M34 `Labeler` 不改）：主核序列上蛋白编码基因 334,699，**CDS-assessable 主参考 331,106**，完整且不重叠的 323,075；外显子分层 1 / 2–5 / ≥ 6 = 67,513 / 149,793 / 113,800 |
| 暴露审计 | 12 物种与 M34 10 个训练物种**无同属、无同 accession**；同科只有禾本科两例（Zizania、Themeda vs 水稻 + Brachypodium），且 sealed 的 Setaria 也只在科级。**Dryas 与 Eschscholzia 是 M30 的适配物种**（其标签集 A 染色体曾用于 M30 适配训练）：面板评测必须用从未读过这两份标签的 checkpoint（M34/M35 系列满足），否则只能评 M30 的 B 集 |
| 蛋白同源 | 全蛋白组（464,568 条，含 isoform）对 M34 T10 训练蛋白库（445,715 条）MMseqs2（`-s 4.0 -e 1e-5`，query 覆盖 ≥ 0.7）最佳命中一致度分布见 §5；一致度 ≥ 70% 的占比 12%（Dipteris）–73%（Zizania），12 物种中位数 34%；与 M34 1,500 条抽样审计的数字逐物种相差 ≤ 2 点 |
| 磁盘 | ont `data/panel/`：raw 2.7 GB + prepared 9.8 GB + homology 1.4 GB（全量命中表 6.0 GB 已 gzip 为 1.13 GB）= **14 GB**（< 15 GB，§7）；准备过程不超过 1 个 Python 进程 + mmseqs 6 线程，均绑定核 25–31、nice 10（核 24 在 M35 的 CPU 集合里，故未使用） |
| 未核实 | SH 层 5 个 GenBank 注释的方法/证据（10 号文档 §2.4 的 BUSCO/RNA-seq 资格测试**未做**）；Themeda 蛋白 45.6% 无命中与 Acorus 非标准剪接 6.6% 提示注释质量问题（§6）；Tiberius Figshare 清单与 ANNEVO 现行权重清单仍未核对 |

## 1. 物种表（12 个；来源 NCBI dataset_report 2026-10-21 + 本次统计）

层：`TH` 严格时间留出；`SH` 来源留出（RefSeq 之外）；`AH` 版本留出。暴露字母同 10 号文档 §1.7（`S` 同物种训练、`V` 仅验证、`T` 仅测试、`G` 同属、`Tr` 同族/亚科、`F` 同科、`O` 仅同目、`–` 无）。"我方"= M34 T10 训练物种的最近关系（§4）。

| # | 物种 | 层 | Accession（装配名；级别） | 注释版本（日期；方法） | NCBI 总长 / 蛋白编码基因 / 注释 BUSCO | 主核序列（条 / Mb；丢弃条数 / Mb / 基因） | CDS-assessable 主参考（1 / 2–5 / ≥6 段） | 我方 | Carbon-A | Helixer | Tiberius | ANNEVO | PGA | PC2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | *Dryas octopetala* | TH | GCF_963921425.1（drDryOcto1.1；Chromosome） | RS_2026_07（2026-07-24；Gnomon） | 229.6 Mb / 21,617 / 98.4 | 9 / 229.3；19 / 0.8 / 25 | 21,467（3,786 / 9,463 / 8,218） | F（无同科；**M30 适配物种**） | F | F | F | F | F | F |
| 2 | *Eschscholzia californica* | TH | GCF_041429995.1（dmEscCali1.0.p；**Scaffold**） | RS_2026_07（2026-07-30；Gnomon） | 400.7 Mb / 26,729 / 99.1 | 6 / 379.2；325 / 22.5 / 5,225 | 26,580（5,122 / 11,797 / 9,661） | O（**M30 适配物种**） | Tr | – | O（测试种） | – | – | O |
| 3 | *Akebia trifoliata* | TH | GCF_017979445.1（ASM1797944v1；Chromosome） | RS_2026_07（2026-07-28；Gnomon） | 652.8 Mb / 25,221 / 97.8 | 16 / 605.6；898 / 47.3 / 1,544 | 23,926（3,891 / 10,985 / 9,050） | O | O | – | O | – | – | O |
| 4 | *Hakea prostrata* | TH | GCF_041222555.1（UWA_Hpros_2；Chromosome） | **RS_2026_10**（2026-10-01；Gnomon） | 712.4 Mb / 24,188 / 98.6 | 11 / 662.8；323 / 49.6 / 702 | 23,486（3,923 / 10,344 / 9,219） | – | F | – | – | O | – | – |
| 5 | *Chenopodium album*（六倍体） | TH（敏感性层） | GCF_948465745.1（dcCheAlbu1.1；Chromosome） | RS_2026_03（2026-03-13；Gnomon） | 1,593.8 Mb / 61,810 / 97.2 | 27 / 1,588.2；131 / 6.1 / 1,564 | 61,683（10,543 / 28,890 / 22,250） | – | **G** | V/G | F | G | – | G |
| 6 | *Zizania latifolia* | SH | GCA_048537475.1（ASM4853747v1；Chromosome） | GenBank 提交注释（武汉植物园，2025-03-10；方法未披露） | 595.8 Mb / 36,280 / 未提供 | 17 / 571.7；272 / 24.1 / 179 | 35,671（8,542 / 15,515 / 11,614） | **Tr**（水稻同族 Oryzeae；Brachypodium 同科） | Tr | Tr | Tr | Tr | Tr | Tr |
| 7 | *Themeda triandra* | SH | GCA_057662145.1（ThtriWBW_v1.0；Chromosome） | GenBank 提交注释（塔斯马尼亚大学，2026-05-21；方法未披露） | 754.7 Mb / 43,880 / 未提供 | 10 / 746.0；12 / 8.7 / 364 | 41,205（14,277 / 19,614 / 7,314） | F（水稻、Brachypodium 同科） | Tr | Tr | Tr | Tr | Tr | Tr |
| 8 | *Juncus effusus* | SH | GCA_027726005.1（MPIPZ_junEff1m_1；Chromosome） | GenBank 提交注释（MPIPZ，2023-01-12；方法未披露） | 258.3 Mb / 18,942 / 未提供 | 21 / 239.9；362 / 18.3 / 1,162 | 17,779（2,463 / 8,025 / 7,291） | O（Poales） | O | O | O | O | O | O |
| 9 | *Rhynchospora breviuscula* | SH | GCA_027562975.1（MPIPZ_rhyBre1m_1.0；Chromosome） | GenBank 提交注释（MPIPZ，2023-01-04；方法未披露） | 420.7 Mb / 24,339 / 未提供 | 5 / 370.1；1,500 / 50.6 / 4,361 | 19,976（3,258 / 9,098 / 7,620） | O（Poales） | Tr | O | O | O | O | O |
| 10 | *Acorus gramineus* | SH | GCA_030737835.1（ASM3073783v1；Chromosome） | GenBank 提交注释（福建农林大学，2023-08-16；方法未披露） | 398.6 Mb / 24,974 / 未提供 | 11 / 289.1；126 / **109.5** / 3,609 | 21,365（5,217 / 9,353 / 6,795） | – | – | – | **G** | – | – | – |
| 11 | *Dipteris shenzhenensis*（蕨） | TH | GCF_047967505.1（Dip_SZ_1.0；Chromosome） | RS_2026_08（2026-08-19；Gnomon） | 1,899.3 Mb / 21,776 / 90.1 | 33 / 1,868.9；593 / 30.4 / 749 | 21,501（3,817 / 9,739 / 7,945） | – | – | – | –（超范围） | –（超范围） | – | – |
| 12 | *Amborella trichopoda* | **AH（单列）** | GCF_054954075.1（AtrichopodaSantaCruz75HAP2_v2.0；Chromosome） | RS_2026_09（2026-09-10；Best-placed RefSeq + Gnomon） | 700.3 Mb / 16,734 / 98.2 | 13 / 699.5；13 / 1.0 / 156 | 16,467（2,674 / 6,970 / 6,823） | – | **S**（AMTR1.0） | V | –（超范围） | T | – | – |
| | **合计** | | | | 8,617 Mb / 346,490 | **179 / 8,250；4,574 / 370 / 19,640** | **331,106（67,513 / 149,793 / 113,800）** | | | | | | | |

说明：(a) "蛋白编码基因"列的 NCBI 数字是整个装配；主核序列上的蛋白编码基因为 334,699（`results/species_table.tsv` 的 `kept_protein_coding_genes`）。(b) GenBank 注释的 `dataset_report` 没有 BUSCO 与方法字段，NCBI 只记"Annotation submitted by …"；5 个 SH 物种的 GFF 结构都是 gene/mRNA/CDS 三层且全部 `gene_biotype=protein_coding`（Themeda 另有 276 个 pseudogene 与 isoform：48,517 mRNA / 43,880 gene），冻结 ruler 能直接解析。(c) Amborella 注释方法含 "Best-placed RefSeq"，即把旧注释的 RefSeq 转录本投到新装配，进一步说明它只能是版本留出。(d) Hakea 的注释发布于 2026-10-01，是面板中最新的，本次已冻结 release 名与 md5。

## 2. 下载与 md5 清单

- 脚本 `scripts/panel_download.sh`（与 `m34_download.sh` 同逻辑，增加了 `_assembly_stats.txt` 与 `ncbi_dataset_report.json`）；日志 `results/ont_download.log`；完整清单 `results/download_manifest.tsv`（60 行：每物种 5 个文件的本地 md5、NCBI md5、状态、下载时间、注释名/日期/方法/BUSCO）。ont 上每个目录另有 `local_md5.txt`、`md5checksums.txt`、`download_date.txt`。
- 全部 60 个文件 `OK`，无 MISMATCH。GFF md5（前 10 位）即注释版本冻结：

| 物种 | genomic.fna.gz | genomic.gff.gz | protein.faa.gz | 下载时间（CEST） | 与此前副本比对 |
|---|---|---|---|---|---|
| Dryas | 6c573c854c | 7f7d087150 | 67ea644936 | 2026-10-21 23:14:04 | 三者 = M30（13 号文档 §1） |
| Eschscholzia | 440dff3a33 | dc5f8bb6c9 | ef9b45b74b | 23:14:33 | 三者 = M30 |
| Akebia | bba2bb78f0 | 962732789a | c5edfd17f3 | 23:15:17 | 蛋白 = M34 panel_proteins |
| Hakea | 3042e35ac2 | 982dc142cf | 73b51cf38e | 23:15:36 | 蛋白 = M34 |
| Chenopodium | 4a507c214c | 2ed159510e | 6b9b432644 | 23:16:23 | 蛋白 = M34 |
| Zizania | c2a6b5e130 | 20bbb52d18 | 0f6a9aadae | 23:17:04 | 蛋白 = M34 |
| Themeda | 4941fab7bf | 512bfef93d | 8e420806da | 23:17:29 | 蛋白 = M34 |
| Juncus | 042a0a730a | 7e887bc5bf | c89afe90fd | 23:17:38 | 蛋白 = M34 |
| Rhynchospora | 9eae82b7d3 | 5bab1ada92 | be1b842b3e | 23:17:55 | 蛋白 = M34 |
| Acorus | caecd1780d | d153ae2a64 | 5b80b568e8 | 23:18:27 | 蛋白 = M34 |
| Dipteris | 81d83708ef | de8fa9ec57 | cb63e84ac1 | 23:19:40 | 蛋白 = M34 |
| Amborella | 64111d86fe | 9268feb7cc | 7f47fa64dd | 23:20:01 | 蛋白 = M34 |

压缩大小合计 FASTA 2,636 MB + GFF 107 MB + 蛋白 100 MB，与 10 号文档 §6.1 的 HEAD 核验一致。

## 3. 主核序列文件（m1_screen / M34 `prepared` 格式）

- 产物 ont `~/coding-rna/data/panel/prepared/<key>/{genome.fa, reference.gff3, lengths.json, splits.json, reference_stats.json}`；`splits.json` 记录规则、每条保留序列的 RefSeq/GenBank accession 与长度、丢弃序列按角色计数、丢弃的注释基因数、NCBI 报告摘要、gff/fna md5。`train=[]`、`val=[]`、`test=全部保留序列`：整个文件都是评测数据，无任何裁剪。
- 规则（`panel_common.main_sequences`，先于任何模型输出固定）：有 assembled-molecule 时用 M34 规则（assembled-molecule + Chromosome + Primary Assembly + ≥ 5 Mb）；没有时（Eschscholzia）用 M30 规则（Primary Assembly 中 ≥ 10 Mb）。GCA 装配的序列名用 GenBank accession（FASTA 头即如此）。GFF 只按序列名过滤，不改坐标，去掉 `region` 行。
- 与 10 号文档 §3.5 的差别：§3.5 曾建议"不按长度删序列"。本次按任务说明采用与 M34/M30 一致的客观规则，并把丢弃量写进 `splits.json`，供论文披露：丢弃最多的是 Acorus（126 条未锚定 scaffold，109.5 Mb，3,609 个基因，占其基因的 14.5%）、Rhynchospora（1,500 条，50.6 Mb，4,361 个基因，17.9%）、Eschscholzia（325 条，22.5 Mb，5,225 个基因，含大量 <1 Mb 的短 scaffold，13.8%）；其余物种丢弃基因 ≤ 6.1%（Juncus）。若评审要求全序列评测，可用同一脚本把阈值设为 0 重跑（无模型依赖）。
- 一致性测试 `tests/test_panel_prep.py`（ont，generanno env，7 项通过，含 prepared 数据检查：FASTA 头 = 保留序列、GFF 坐标 ≤ 长度、md5 状态 OK、`splits.json` 自洽）。

参考统计（`results/reference_stats.tsv`，冻结 ruler 的 `CDS_assessable_primary_R4` 口径，M34 `Labeler`/`stats` 不改）：

| 物种 | 主核 Mb | 解析的编码基因 | CDS-assessable 主参考 | 每 Mb | 完整且不重叠 | 中位段数 / 中位 CDS bp | CDS 跨度 > 24,576 bp | 非标准剪接占比 | 不可评估原因（非标准剪接 / 末端密码子跨段 / 重叠或部分） |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|
| Dryas | 229.3 | 21,602 | 21,467 | 93.6 | 20,735 | 4 / 1,182 | 62 | 2.8% | 609 / 37 / 86 |
| Eschscholzia | 379.2 | 26,650 | 26,580 | 70.1 | 26,042 | 4 / 1,170 | 76 | 1.5% | 397 / 59 / 82 |
| Akebia | 605.6 | 24,038 | 23,926 | 39.5 | 22,855 | 4 / 1,128 | 710 | 3.8% | 898 / 58 / 117 |
| Hakea | 662.8 | 23,679 | 23,486 | 35.4 | 22,400 | 4 / 1,185 | 1,687 | 4.2% | 992 / 36 / 60 |
| Chenopodium | 1,588.2 | 61,773 | 61,683 | 38.8 | 60,629 | 4 / 1,164 | 628 | 1.3% | 794 / 102 / 160 |
| Zizania | 571.7 | 36,101 | 35,671 | 62.4 | 35,570 | 3 / 999 | 144 | 0.2% | 57 / 30 / 14 |
| Themeda | 746.0 | 43,516 | 41,205 | 55.2 | 40,582 | **2 / 792** | 64 | 0.9% | 372 / 14 / 239 |
| Juncus | 239.9 | 17,780 | 17,779 | 74.1 | 17,720 | 4 / 1,140 | 23 | 0.2% | 34 / 20 / 5 |
| Rhynchospora | 370.1 | 19,978 | 19,976 | 54.0 | 19,856 | 4 / 1,161 | 147 | 0.3% | 52 / 33 / 35 |
| Acorus | 289.1 | 21,365 | 21,365 | 73.9 | 19,921 | 3 / 912 | 228 | **6.6%** | 1,418 / 28 / 0 |
| Dipteris | 1,868.9 | 21,587 | 21,501 | **11.5** | 20,811 | 4 / 1,272 | **5,095** | 1.9% | 402 / 40 / 254 |
| Amborella | 699.5 | 16,630 | 16,467 | 23.5 | 15,954 | 4 / 1,176 | 1,667 | 2.3% | 385 / 35 / 93 |

观察：Dipteris 有 5,095 个（23.7%）参考基因 CDS 跨度超过 24,576 bp 的核心窗口（蕨类大内含子），Hakea 1,687（7.2%）、Amborella 1,667（10.1%）——这些基因在我方 24 kb 窗口方案下是结构性难例，论文按段数/长度分层时需单列。Themeda 单外显子基因占 35%、中位 CDS 仅 792 bp，与 §5 的 45.6% 无同源命中一起指向注释含大量低置信/TE 样模型（§6）。Acorus 非标准剪接 6.6%（其余 ≤ 4.2%），这些基因被冻结 ruler 标为不可评估，不进入 F1 分母。

## 4. 暴露审计（`scripts/panel_exposure.py` → `results/exposure_audit.json`）

训练侧集合：M34 10 个训练物种（拟南芥 GCF_000001735.4、水稻 GCF_034140825.1、杨、短柄草、番茄、椰枣、菜豆、甜橙、胡萝卜、薯蓣）、M30 适配物种（Dryas、Eschscholzia）、项目 sealed 的 Setaria。规则：面板物种不得与 M34 训练物种同属或同 accession；同科只报告。结果：

- **同属 / 同 accession：0**；M34 方向（`panel_exclusion.audit`）的违规也为 0；12 个面板 accession 都在 M34 的排除名单里且属名一致（`test_species_registry_matches_m34_exclusion_list`）。
- 同科：Zizania、Themeda 与水稻 + Brachypodium（禾本科）；Setaria 也只到科级。结论与 19 号文档 §1.3 相同：禾本科面板物种必须与非禾本科分开报告。
- **M30 暴露**：Dryas、Eschscholzia 的标签集 A（Dryas S1/S3/S5/S7 + S9 前 20 Mb；Eschscholzia SCAF_1/3/5）在 M30 中用于适配训练。对论文的最终模型（M34/M35 路线，从未读这两份 GFF）这不是泄漏；但**任何由 M30 适配 checkpoint 延续的权重**只能在 B 集（Dryas S2/S4/S6/S8；Eschscholzia SCAF_2/4/6）上评，或在面板表里标注。`exposure_audit.json` 对这两个物种记录了 `m30_adaptation_exposure`。
- 各工具暴露（`scripts/tool_exposure.json`，转录自 10 号文档 §1.7，未新增核实）见 §1 表。本次没有核对 Tiberius Figshare 清单与 ANNEVO 现行权重清单（10 号文档 §7.3 待办 2/3 仍开放）。

## 5. 蛋白同源审计（`scripts/panel_homology.py`，ont 2026-10-21 23:55 → 10-22 08:40，mmseqs 6 线程 ≈ 8.7 h）

设置与 M34 §2.5 相同（MMseqs2 18.8cc5c，`easy-search -s 4.0 -e 1e-5 --max-seqs 300`，最佳命中要求 query 覆盖 ≥ 0.7，一致度为比对区内 pident），差别是查询为**全部蛋白记录**（464,568 条，含 isoform）而非每物种 1,500 条抽样，并给出完整分布。目标库复用 M34 的 `train_T10.faa`（445,715 条）；T2 = 其中拟南芥 + 水稻。纯描述，不据此过滤任何数据。

| 面板物种 | 蛋白数 | 无命中 | <30 | 30–40 | 40–50 | 50–70 | 70–90 | ≥90 | ≥50 | **≥70** | ≥90 | 有命中者中位一致度 | T2：无命中 / ≥70 | 最近训练物种（最佳命中归属前 3） |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| Dryas（TH） | 30,373 | 3.5% | 0.9% | 4.0% | 7.5% | 32.2% | 47.1% | 4.8% | 84.1% | **51.8%** | 4.8% | 71.3 | 9.6% / 29.1% | 杨 10,876；橙 10,631；菜豆 3,006 |
| Eschscholzia（TH） | 42,020 | 6.3% | 3.3% | 7.5% | 10.4% | 36.9% | 32.2% | 3.4% | 72.5% | **35.6%** | 3.4% | 64.7 | 12.7% / 24.5% | 杨 9,984；橙 9,690；椰枣 6,476 |
| Akebia（TH） | 43,352 | 8.9% | 1.2% | 5.8% | 8.9% | 34.7% | 36.5% | 4.0% | 75.2% | **40.5%** | 4.0% | 67.7 | 15.6% / 25.0% | 杨 10,583；橙 10,250；椰枣 7,250 |
| Hakea（TH） | 34,200 | 2.8% | 1.6% | 6.7% | 8.9% | 38.8% | 37.4% | 3.7% | 80.0% | **41.2%** | 3.7% | 66.9 | 8.9% / 26.2% | 杨 8,859；橙 8,546；椰枣 5,740 |
| Chenopodium（TH） | 77,290 | 6.3% | 3.4% | 10.8% | 12.5% | 36.5% | 28.3% | 2.2% | 67.1% | **30.6%** | 2.2% | 62.0 | 14.1% / 21.0% | 橙 19,224；杨 17,689；胡萝卜 9,224 |
| Zizania（SH） | 36,280 | 12.5% | 0.1% | 1.0% | 2.0% | 11.3% | 48.8% | 24.2% | 84.3% | **73.0%** | 24.2% | 84.7 | 13.6% / 72.3% | **水稻 29,338**；短柄草 1,770 |
| Themeda（SH） | 48,213 | **45.6%** | 0.7% | 3.1% | 4.9% | 15.8% | 25.2% | 4.6% | 45.7% | **29.9%** | 4.6% | 72.3 | 46.9% / 28.4% | 水稻 16,471；短柄草 8,655 |
| Juncus（SH） | 18,942 | 3.6% | 0.9% | 6.6% | 12.5% | 41.9% | 29.8% | 4.7% | 76.4% | **34.6%** | 4.7% | 63.9 | 6.4% / 27.6% | 椰枣 8,137；水稻 3,395；薯蓣 2,378 |
| Rhynchospora（SH） | 24,339 | 8.4% | 2.2% | 7.6% | 11.3% | 36.5% | 27.2% | 6.8% | 70.5% | **34.0%** | 6.8% | 64.3 | 11.9% / 27.2% | 椰枣 9,353；水稻 4,059；薯蓣 2,920 |
| Acorus（SH） | 24,974 | 14.6% | 0.9% | 5.6% | 11.3% | 35.5% | 28.8% | 3.3% | 67.6% | **32.1%** | 3.3% | 64.3 | 20.4% / 21.3% | 椰枣 9,807；薯蓣 3,380；杨 2,085 |
| Dipteris（TH，蕨） | 32,820 | 21.5% | 3.1% | 18.2% | 15.4% | 29.8% | 11.1% | 1.1% | 42.0% | **12.2%** | 1.1% | 51.6 | 27.3% / 9.5% | 椰枣 6,238；杨 3,596；橙 3,099 |
| Amborella（AH） | 51,765 | 6.3% | 1.2% | 6.4% | 12.5% | 40.1% | 30.8% | 2.8% | 73.7% | **33.6%** | 2.8% | 64.2 | 11.4% / 22.9% | 椰枣 18,114；橙 7,109；杨 7,026 |
| **12 物种中位数** | | 7.4% | | | | | | | 73.1% | **34.3%** | 3.9% | 64.5 | | |

解读：
1. 与 M34 §2.5 的 1,500 条抽样（≥70% 列：Dryas 53、Esch 35、Akebia 41、Hakea 43、Cheno 30、Dipteris 12、Zizania 74、Themeda 29、Juncus 34、Rhyncho 33、Acorus 33、Amborella 33）逐物种相差 ≤ 2 点，抽样审计的结论成立。
2. 分层口径：按 10 号文档 §3.2 的 H（≥70% 且覆盖 ≥70%）/ M（40–70%）/ L（30–40%）/ N（无命中）档，各物种 H 档占 12–73%，N 档 3–46%；论文的"同源剥离"表（§3.3 表 C，去掉 ≥50% 的参考基因）会去掉 42–84% 的参考蛋白，剩余最多的物种是 Dipteris（58%）、Themeda（54%）、Chenopodium（33%）。此表给的是蛋白记录级，不是基因级（isoform 未合并，Amborella 51,765 条蛋白对 16,630 个基因）；基因级剥离集合需在评测时用 GFF 把蛋白映射回基因后再算。
3. Zizania 对 T2 已 72.3%，加入 8 个物种后 73.0%，几乎全部来自水稻（29,338/36,280 最佳命中），是"近缘层"的定量依据；Themeda 的最佳命中仍以水稻为主（16,471）而非短柄草（8,655），Brachypodium 没有把它拉近。
4. 最近训练物种的分布显示 M34 的杨/橙是真双子叶面板物种的主要"近邻"，椰枣是单子叶非禾本科、蕨类与 Amborella 的主要近邻；这是 M34 把训练集向面板拉近的具体路径，与 19 号文档 §2.5 第 1 条一致。
5. 本审计只对我方训练集（E_ours = M34 T10）做；10 号文档 §3.2 第 6 条要求对每个基线的 E_base(X) 分别审计（需 201 个 RefSeq 物种蛋白库，约 1.5–2 GB），本次未做。

## 6. 风险与未核实项

1. **SH 层注释资格未测**：5 个 GenBank 注释的方法（是否含 Helixer/BRAKER 等 ab initio 预测）、RNA-seq 证据与 BUSCO 完整度本次均未核实（10 号文档 §2.4 门槛 2）。本次统计给出两个警示：Themeda 蛋白 45.6% 在 10 个训练物种中无任何 e ≤ 1e-5 命中（其余物种 3–22%），单外显子基因占 35%，中位 CDS 792 bp；Acorus 非标准剪接 6.6% 且 14.5% 的注释基因在未锚定 scaffold 上。两者在做 compleasm/BUSCO 资格测试前不应进入主表；备选替换顺序按 10 号文档 §2.3（Canna/Stemona/Paspalum/Aristolochia）。
2. **序列范围**：本次按 M34/M30 的客观长度规则丢弃了 4,574 条序列（0.37 Gb、19,640 个注释基因，占全部注释基因的 5.0%）；与 10 号文档 §3.5"不按长度删序列"的建议不同，已在 §3 说明并在 `splits.json` 记录丢弃明细，评审若要求可无模型地重跑。
3. **Dryas/Eschscholzia 的 M30 暴露**（§4）：评测前须确认所用 checkpoint 的训练谱系。
4. **蛋白级 vs 基因级**：§5 为蛋白记录级（含 isoform）；同源分层与剥离在评测时须折算到基因级。
5. **基线训练清单**：Tiberius Figshare、ANNEVO 现行权重、GeneCAD 训练名单仍未核对；PlantCAD2 预训练 65 个被子植物与本面板的关系只到属级（Chenopodium quinoa）。
6. **Hakea 注释极新**（2026-10-01）：若 NCBI 再发布 RS_2026_11 等版本，以本文 md5 为准不更新。
7. **Amborella** 注释含 "Best-placed RefSeq"（旧注释投射），与 Carbon-A 训练内的 AMTR1.0 注释高度相关，只能作版本留出敏感性分析。
8. 本任务没有做 BUSCO、RNA-seq 可得性、TimeTree 分歧时间、Mash/ANI 基因组级复核（10 号文档 §3.2 第 7 条）。

## 7. 磁盘与运行记录

- ont `~/coding-rna/data/panel/`：`raw/` 2.7 GB；`prepared/` 9.8 GB（12 个 `genome.fa` 共 8.25 Gb 未压缩 + GFF3）；`homology/` 运行结束时 5.9 GB（`best_hits_full.tsv` 6.0 GB + `panel_full.faa` 0.23 GB，现 1.36 GB），总计 18.4 GB，**超过任务的 15 GB 上限**；随后把 `best_hits_full.tsv` 改为 gzip 保存（`best_hits_full.tsv.gz` 1.13 GB，10-22 09:08 完成），**现总量 14 GB（`du -sh`），在 15 GB 内**；`homology_audit_full.json` 与 `homology_table.tsv` 已包含全部汇总数字，不需要解压即可引用。ont 根分区余 1.9 TB。
- 进程：下载 1 个 curl 串行；prep/stats 1 个 Python 进程（Dipteris/Chenopodium 峰值内存约 6–8 GB）；mmseqs 6 线程；全部 `taskset -c 25-31 nice -n 10`。未触碰 M35（核 2–24，GPU）与 M37 的任何文件/进程。
- 时间：下载 7 min；prep 8.5 min；stats 27 min；homology 8.7 h（prefilter 索引 2 min，主要是 46 万查询的比对）；总计约 9.3 h（10-21 23:13 → 10-22 08:42 CEST）。
- 日志：`results/ont_download.log`、`results/ont_run.log`；ont `~/coding-rna/outputs/PANEL-PREP/{download.log,run.log,results/}`。

## 8. 产物清单

| 位置 | 内容 |
|---|---|
| `experiments/PANEL-PREP/scripts/` | `panel_species.json`（12 物种登记）、`tool_exposure.json`（§1.7 暴露码）、`panel_common.py`（路径与序列规则）、`panel_download.sh`、`panel_prep.py`（prep / stats）、`panel_exposure.py`、`panel_homology.py`、`panel_tables.py`、`panel_run.sh`（ont 串行驱动） |
| `experiments/PANEL-PREP/tests/test_panel_prep.py` | 7 项（序列规则两分支、GCA/GCF 序列名、登记表与 M34 排除名单一致、暴露审计、同源分档、prepared 数据自洽） |
| `experiments/PANEL-PREP/results/` | `download_manifest.tsv`、`species_table.tsv`、`reference_stats.tsv`、`exposure_audit.{json,txt}`、`homology_audit_full.json`、`homology_table.tsv`、`disk_usage.txt`、`ont_download.log`、`ont_run.log` |
| `experiments/PANEL-PREP/{push.sh,wait_event.sh}` | 推送与后台事件等待器（经 baobab 跳板到 ont） |
| ont `~/coding-rna/data/panel/` | `raw/`、`prepared/`、`homology/`（见 §7） |

下一步（不在本任务内）：SH 层注释资格测试（compleasm odb10/odb12 + SRA RNA-seq 元数据）；基因级同源剥离集合；各基线 E_base(X) 审计；面板推理排在 M35 → M34 → M36 之后（11 月）。
