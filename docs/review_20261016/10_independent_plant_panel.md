# 10 独立植物评测面板设计（只读调研，2026-10-16）

产物说明：本文件是唯一新增产物。调研期间未修改、移动、删除任何已有本地/远端文件，未提交 Slurm；未下载任何大文件（只下载了以下小文件到会话临时目录：Carbon-A 技术报告 PDF 3.8 MB、Helixer 补充 PDF 4.9 MB、ANNEVO 补充 PDF/XLSX 约 15 MB、Tiberius 补充 PDF 0.5 MB、PlantCAD2/PlantGeneAnn 补充 XLSX 各 <0.3 MB）；基因组 FASTA/GFF 只做了 HEAD 请求。

来源标记：**[PDF]** Carbon-A 技术报告（`https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf`，本机已下载解析，Baobab `refs/carbon-a/` 下只有模型仓库快照、没有 PDF 与补充表）；**[NCBI]** 2026-10-16 当天对 `https://api.ncbi.nlm.nih.gov/datasets/v2/` 的 dataset_report / taxonomy 查询；**[推断]** 本人推断；**[未核实]** 未能验证。所有 F1 数字为各作者自报，不可跨面板拼接。

---

## 0. 结论先行（含对 05/06/07 的更正）

1. **推荐主面板 12 个物种，共 8.6 Gb**（5 个 RefSeq-2026 时间留出被子植物 + 5 个 RefSeq 之外的单子叶 + 1 个蕨类 + 1 个 Amborella 新装配体）。另有 15 个备选。主面板 FASTA+GFF3+蛋白合计约 **2.84 GB**（FASTA 2.64 GB）。
2. **Amborella `GCF_054954075.1` 不能当作 Carbon-A 的时间留出物种**。核实结果：新 GCF 的装配（JGI HiFi，2026-02-09）和 RefSeq 注释（`RS_2026_09`，2026-09-10）确实晚于 Carbon-A 权重，但同一物种旧装配 `GCF_000471905.2`（AMTR1.0，RefSeq 注释 2017-03-29，BUSCO 97.8）早已在 RefSeq，因此在 Carbon-A 训练快照内；Helixer 补充表 S18 也把 Atrichopoda v1.0 列为 train_val（验证），ANNEVO 论文 Table 24 把 `GCF_000471905.2` 列为 Test。**只能作"装配/注释版本留出"，不是物种留出**，且两代装配很可能来自同一个体 Santa Cruz 75（[NCBI] 新装配名 `AtrichopodaSantaCruz75HAP2`），序列近乎同一，必须单列。
3. **Carbon-A 的"时间留出"边界可以精确到证据级**：（a）14 个官方留出基因组的 RefSeq 注释日期为 2026-02-23 至 2026-07-24（[NCBI] 逐个查询，见 §1.1）；（b）与 GENERanno-preview 同权重的公开日期为 2026-02-10（`https://github.com/GenerTeam/GENERanno` README News）。所以"RefSeq 首次注释发布 ≥ 2026-02-23 且物种此前无任何 GCF"可以被证明不在 Carbon-A 训练数据里；快照确切日期未披露，按植物 182 个物种/accession 反推更像 2025 年中（§1.1，仅为推断）。
4. **Carbon-A 权重不是从 GENERator 预训练权重继续训练的**。[PDF] Table 4 明写 Carbon-A-mammal 为 "Random initialization"，Table 3 只说架构"derived from GENERator-1.2B backbone"（把因果解码器换成双向编码器）。因此 06 号文件 §0 第 1 条里"骨干来自 GENERator-1.2B"只应理解为架构血缘；GENERator 的 RefSeq 基因区预训练暴露（植物 7.28 M 基因/30.5 Gbp，arXiv 2502.07272 Table S13）对这份权重**不构成**数据暴露路径（除非另有未披露的初始化，作者原文没有）。
5. **05 §3.4 的 Setaria viridis 暴露已可以补齐**：Carbon-A（RefSeq `GCF_005286985.2`，注释 2025）、Tiberius angiosperms 训练（Sviridisv4.1）、PlantGeneAnn 训练、PlantCAD2 预训练均**训练过**；Helixer（a0080）与 ANNEVO 论文版是**测试**。结论：Setaria 对 Carbon-A/Tiberius 不是盲测，彻底放弃作新颖盲测。
6. **ANNEVO 现行植物权重的训练物种没有公开清单**。论文版（Embryophyta）清单在补充表 Table 24（40 Train + 9 Val&Test + 119 Test），而现行 Magnoliopsida 权重是 2026-04 发布的新模型，作者只说"排除模式生物（如拟南芥、水稻）并更新了未发表的数据预处理"（`docs/boundary_aware_model.md`），无清单。**对 ANNEVO 的物种暴露只能写"论文版已知、现行版未知"**，面板选物种时按"两版都未见过"的保守口径（即同时不在 Table 24 的 Train/Val 里）。
7. **面板的最大三个风险**：（i）RefSeq 之外的 5 个单子叶，其 GenBank 提交注释的证据类型/质量没有逐一核实（BUSCO 需自测）；（ii）RefSeq 参考注释（Gnomon）与 Carbon-A 训练标签同风格，对 Carbon-A 有"标签风格偏置"，需用非 RefSeq 来源的物种做分层对照；（iii）所有物种在 Carbon-A 快照里都有同科/同目近缘，同源审计不是可选项（§3）。
8. 算力：主面板我方推理约 **141 GPU-h（4090，17 kbp/s）**，Carbon-A 约 **72–96 GPU-h**（06 号文件的 15–20 min/30 Mb）；去掉两个大基因组（Chenopodium album 1.6 Gb、Dipteris 1.9 Gb）的"核心 10 物种"为 5.1 Gb、我方约 84 h、Carbon-A 43–57 h（§5）。

---

## 1. 暴露矩阵

### 1.1 Carbon-A（= GENERanno-eukaryote-1.2b-cds-annotator-preview 同权重）

| 项 | 内容 | 出处 |
|---|---|---|
| 训练数据 | 仅 GCF 开头的 RefSeq 装配，记录长度 ≥ 98,304 bp；2,055 个装配/2,047 个物种；植物 182 个物种、182 个装配、源 335.42 Gbp；植物滑窗步长 49,152 bp | [PDF] §5.2、Table 1/2 |
| 初始化 | 随机初始化，课程式累加类群，植物在第 4 阶段加入（~0.5 epoch，之后 fungi+protozoa 2 epoch） | [PDF] Table 4 |
| 训练物种清单 | **不逐个披露**（只给计数）。基准里 28 个"训练快照内"参考基因组中的 5 个植物为：拟南芥 `GCF_000001735.4`、莱茵衣藻 `GCF_000002595.2`、水稻 `GCF_034140825.1`、小立碗藓 `GCF_000002425.5`、毛果杨 `GCF_000002775.5` | [PDF] Table S1/S3/S4/S5 |
| 时间留出 | 14 个"2026 年起新增到 RefSeq"的基因组，**均为动物**（哺乳 4、其他脊椎 5、无脊椎 5），**没有植物** | [PDF] Table S1 星号行；§2.1 |
| 留出基因组的 RefSeq 注释日期 | 最早 2026-02-23（Giraffa tippelskirchi `GCF_054371585.1`，RS_2026_02），最晚 2026-07-24（Bos grunniens `GCF_055505025.1`）；其中 Ips typographus `GCF_016097725.1` 的 GenBank 装配是 2021 年的，但 **GCF 在 2026-04 才建立**，说明留出判据是"GCF/RefSeq 注释首次出现时间"，而不是装配发表时间 | [NCBI] 逐个 dataset_report |
| 同权重的公开日期 | GENERanno-eukaryote-1.2b-cds-annotator-preview：2026-02-10 | GENERanno README News |
| 辅助网络（转移 TCN 11,077 参数、置信 TCN 14,433 参数） | 训练集/划分未披露 | 06 号文件 §1.3；[PDF] §5.5.2/5.5.4 |
| 由 RefSeq 反推快照 | 当前 RefSeq 植物共 214 个当前装配（356 个含历史版本），按物种（二名法）首次 GCF 版本日期 < 2026-02-10 的有 201 个，< 2025-07-01 有 191 个，< 2025-04-01 有 181 个（[NCBI]）。哺乳类同法：236 个物种在 ≈2025-07-01 前达到（报告值 236）。**[推断] 快照大概率在 2025 年中前后；在没有作者确认前，把"2025-04-01 至 2026-02-10 之间首次出现的 GCF"视为暴露未知，不用作留出** | 本节脚本复现见 §7.2 |

**对 Carbon-A 的"留出资格"判据（本文使用）**
- 严格时间留出 `TH`：该物种此前**没有任何 GCF 版本**，且首个 RefSeq 注释发布 ≥ 2026-02-23（不早于官方留出集最早日期，晚于权重公开日 13 天）。
- 来源留出 `SH`：物种在 RefSeq 里从未出现（只在 GenBank/社区注释），Carbon-A 只用 GCF，故物种级未见；但**不是时间留出**（序列数据在权重发布前已公开）。
- 版本留出 `AH`：物种在快照内但装配/注释是新版（如 Amborella）；对 Carbon-A 视同训练物种。

### 1.2 GENERanno 预训练与 CDS 后训练

| 项 | 内容 | 出处 |
|---|---|---|
| GENERator 预训练 | 386 B 核苷酸真核 DNA，来自 RefSeq 注释的**基因区间**（含内含子/UTR）；植物 7,279,591 基因、30.55 Gbp；全部真核 RefSeq，快照日期未写 | arXiv 2502.07272v5 Table S13、§4；`https://arxiv.org/abs/2502.07272` |
| CDS annotator preview | README 表中 Data 记 1T；无物种/装配清单、无评估指标；模型卡 "More technical details are coming soon" | `https://github.com/GenerTeam/GENERanno`；HF 模型卡 |
| 与 Carbon-A 的关系 | 权重 sha256 与 Carbon-A 的 `model.safetensors` 一致（06 号文件已核）；Carbon-A 报告写明随机初始化，所以 GENERator 预训练不进入这份权重的数据路径（§0.4） | [PDF] Table 4；`provenance.json` |
| 对我方模型的含义 | 我方骨干 = Carbon-A 权重，故**对我方模型而言面板物种暴露 = Carbon-A 暴露 + 我方头训练物种（拟南芥/水稻，若扩训则加上扩训物种）** | 本文 |

### 1.3 Tiberius（2026 多谱系版；被子植物 `angiosperms` 模型）

出处：`https://doi.org/10.64898/2026.04.24.720536`（v3，2026-10-01）及其补充 PDF Table S3、S17、S32（`https://www.biorxiv.org/content/biorxiv/early/2026/10/01/2026.04.24.720536/DC1/embed/media-1.pdf`）；固定配置 `model_cfg/angiosperms.yaml`（commit `871d497…`，date 2026-04-15，tiberius_version 1.1.6，target `Mesangiospermae`）。

- **训练 45 个物种**（多数为 Phytozome 注释，部分来自 NCBI/GDR/Figshare），**验证 2 个**（Populus trichocarpa、Typha latifolia），**测试 7 个**（Arabidopsis thaliana、Brachypodium distachyon、Eschscholzia californica、Freycinetia multiflora、Medicago truncatula、Mimulus guttatus、Urochloa brizantha）。完整名单见附录 A。
- 部分训练基因组"发表时尚未公开"，名单在 Figshare `https://doi.org/10.6084/m9.figshare.32086728`（Fort Lauderdale 约束），**[未核实]** 我没有打开该 Figshare；选面板前须核对其中是否含面板物种。
- 适用范围：Mesangiospermae（被子植物除 Amborellales/Nymphaeales/Austrobaileyales）。蕨类、苔藓、Amborella 超出训练/适用范围，只能作 OOD 探索，不进入"公平主表"。
- 另有 `chlorophyta` 模型（2026-01-26，9 个训练物种：Ostreococcus lucimarinus CCE9901、Chlamydomonas reinhardtii、Volvox carteri、Picochlorum sp. BPE2、Astrephomene gubernaculifera、Nannochloris sp. desiccata、Dunaliella primolecta、Chloropicon primus、Tetradesmus obliquus；测试 Bathycoccus prasinos `GCF_002220235.1` 与 Edaphochlamys debaryana `GCA_016858145.1`）。注意 Bathycoccus 在 Carbon-A 训练快照里。
- 作者自述：33 个测试物种中 Helixer 已将 17 个纳入训练或验证、ANNEVO 2 个、Tiberius 0 个（Table S17）。

### 1.4 Helixer（`land_plant_v0.3_a_0080.h5`，v0.3.7 默认植物模型）

出处：Nature Methods 补充 PDF（`https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-025-02939-1/MediaObjects/41592_2025_2939_MOESM1_ESM.pdf`）Table S18（物种/版本/划分）、S23（训练集组成，我用渲染图核对了 seed_009_017 列）、S29（模型与 data_dir 对应）；README 默认模型 `land_plant_v0.3_a_0080.h5`（`https://github.com/weberlab-hhu/Helixer`）。

- a0080 对应 data_dir `seed_009_017`：**训练 32 个物种**（Acomosus、Athaliana、Bdistachyon、Carietinum、Cgrandiflora、Creinhardtii、Crubella、Esalsugineum、Fvesca、Gmax、Ljaponicus、Lusitatissimum、Macuminata、Mesculenta、Mguttatus、Mpolymorpha、MpusillaCCMP1545、Mtruncatula、Ncolorata、Olucimarinus、Osativa、Othomaeum、Pacutifolius、Ppersica、Pumbilicalis、Sbicolor、Slycopersicum、Sparvula、Taestivum、Tpratense、Vcarteri、Zmays）；**S18 中其余 train_val 物种每物种随机 800 条子序列用于验证/早停**（即约 45 个物种有验证暴露，含 Amborella、Chenopodium quinoa、Brassica rapa 等）。
- **测试 13 个植物**：Arachis hypogaea、Brassica napus、Cannabis sativa、Coffea arabica、Hibiscus syriacus、Nicotiana attenuata、Oryza brachyantha、Papaver somniferum、Phoenix dactylifera、Setaria viridis、Solanum pennellii、Triticum dicoccoides、Vitis riparia（S18 test (used)）。**它们全部是 RefSeq GCF**，因此全部在 Carbon-A 快照内——"在 Helixer 测试集上比较"对 Carbon-A 不成立。
- 其他 land_plant 变体（a0090/0100/0200/0300/0400，m0100/m0200）用不同 data_dir（S29）；若换用它们，暴露需重算，本文固定 a0080。

### 1.5 ANNEVO

出处：Nature Methods 补充 XLSX（`https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-026-03036-7/MediaObjects/41592_2026_3036_MOESM6_ESM.xlsx`）Table 24（Embryophyta 物种用途）、Table 33（亚谱系）；仓库 `https://github.com/xjtu-omics/ANNEVO`（v2.3.3，`saved_model/ANNEVO_Magnoliopsida.pt`，许可：ANNEVO Non-Commercial License）。

- 论文版 Embryophyta：**40 Train、9 Val&Test、119 Test**（名单与 accession 见 Table 24；摘要见附录 A）。例：Oryza sativa `GCF_034140825.1`、Nelumbo nucifera、Medicago truncatula、Musa acuminata 为 Train；Arabidopsis thaliana、Amborella、Papaver somniferum、Setaria viridis 为 Test。
- 现行 Magnoliopsida 权重：**训练清单未公开**；作者说明"排除模式生物如拟南芥、水稻"（`docs/boundary_aware_model.md`）。
- 范围：现行权重只有 Magnoliopsida（NCBI 口径的全体被子植物），**无蕨类/苔藓/藻类模型**。

### 1.6 PlantCAD2（备选骨干）与其他植物模型

| 模型 | 预训练/训练数据 | 出处 |
|---|---|---|
| PlantCAD2（bioRxiv 2025.08.27.672609 v4，已发表于 Cell Genomics，DOI 10.1016/j.xgen.2026.101329；期刊名由 DOI 前缀推断） | **65 个被子植物基因组（每属一个）**，Solanaceae 整科排除；名单见补充表 1（附录 A）；无注释监督（DNA 掩码预训练，采样偏重编码/调控区） | `https://www.biorxiv.org/content/10.1101/2025.08.27.672609v4`；补充 `…/DC1/embed/media-1.xlsx` |
| PlantGeneAnn（Mamba-MoE，可选基线） | 微调 40–42 个 RefSeq 物种（RefSeq 2025-04-10 下载），验证 8 个；预训练 PlantBiMoE 的 42 个物种；基准把与 Helixer 测试集重叠的 5 个物种换成 Arachis duranensis、Oryza glaberrima、Nicotiana sylvestris、Gossypium raimondii、Daucus carota | `https://doi.org/10.64898/2026.06.25.733695`（v2，2026-10-14）补充表 1/2 |
| GeneCAD（PlantCAD2+CRF，可选基线） | 训练物种清单**未核实**（本次未读补充）；5 个留出：Juglans regia、Coffea arabica、Phaseolus vulgaris、Nicotiana sylvestris、N. tabacum | 05 号文件 §1.1；`https://doi.org/10.1101/2025.10.31.685877` |
| AgroNT | 48 个植物基因组预训练，**[未核实]** 清单 | `https://doi.org/10.1038/s42003-024-06465-2` |

### 1.7 面板物种 × 方法 暴露矩阵（主面板，物种级与最近已见近缘）

字母：`S`=同物种已用于训练/预训练；`V`=仅验证；`T`=仅测试（该物种出现在其论文测试集，非训练）；`G`=同属有训练物种；`Tr`=同族/亚科；`F`=同科；`O`=仅同目；`–`=同目也没有。最近近缘取自 NCBI 分类（[NCBI] taxonomy lineage 求最近共同祖先）。ANNEVO 列只对应论文版 Table 24，现行权重未知。

| 面板物种 | Carbon-A / 我方骨干 | Helixer a0080 | Tiberius angiosperms | ANNEVO（论文版） | PlantGeneAnn | PlantCAD2 | 备注 |
|---|---|---|---|---|---|---|---|
| Dryas octopetala | F（Rosaceae：Fragaria/Prunus/Malus/Rosa/Argentina 等） | F（Fragaria、Prunus 训练） | F（Fragaria、Malus） | F（Prunus、Rosa 训练） | F（Prunus、Rosa） | F | 蔷薇科同科近缘，非同属 |
| Eschscholzia californica | Tr（Papaver somniferum，罂粟科同亚科） | –（Papaver 仅为 Helixer 测试物种，非训练） | O（Aquilegia）；**本物种是 Tiberius 测试物种（非训练）** | –（Papaver 在其论文里是 Test） | – | O（Aquilegia） | |
| Akebia trifoliata | O（Ranunculales：Papaver） | – | O（Aquilegia） | – | – | O | 木通科 |
| Hakea prostrata | F（Proteaceae：Macadamia、Telopea） | – | – | O（Nelumbo nucifera 训练，Proteales） | – | – | |
| Zizania latifolia | Tr（Oryzeae：Oryza 三种） | Tr（Oryza sativa 训练） | Tr（O. sativa 训练） | Tr（Oryza 训练） | Tr | Tr | **与水稻同族，近缘；用于"近缘"层** |
| Themeda triandra | Tr（Andropogoneae：Sorghum/Zea/Miscanthus） | Tr | Tr | Tr | Tr | Tr | 高粱族近缘 |
| Juncus effusus | O（Poales） | O | O | O | O | O | 灯心草科，非禾本科 Poales |
| Rhynchospora breviuscula | Tr（Cyperoideae：Carex rostrata，RefSeq 2025） | O | O | O | O | O | 莎草科 |
| Acorus gramineus | –（无 Acorales） | – | **G（Acorus americanus 训练）** | – | – | – | 菖蒲属，单子叶最早分支 |
| Chenopodium album（六倍体） | **G（C. quinoa）** | V（C. quinoa 仅验证）G | F（Beta vulgaris） | G（C. quinoa Val&Test） | –（C. quinoa 仅验证） | G（C. quinoa） | **同属泄漏风险，放入敏感性层** |
| Dipteris shenzhenensis（蕨） | –（无蕨类；仅 Selaginella/苔藓） | – | –（超出适用范围） | –（超出适用范围） | – | – | 全谱系留出 |
| Amborella trichopoda | **S（AMTR1.0 `GCF_000471905.2`）** | **V**（Atrichopoda v1.0 验证） | –（超出适用范围） | T（论文版 Test） | – | – | **只能当"装配版本留出"** |

我方模型：拟南芥 TAIR10.1 与水稻 GCF_034140825.1 加 Carbon-A 暴露（同上第一列）。

---

## 2. 候选面板

### 2.1 筛选过程与数字（可复现，[NCBI]）

1. 全部 RefSeq 植物（taxon 33090，`filters.assembly_source=refseq`）：214 个当前装配（356 条含历史版本），均有 NCBI EGAP 注释（`Gnomon` 等），**全部报告了 BUSCO**（`annotation_info.busco`）。
2. 对每个物种（二名法）取其所有 GCF 版本的最早注释/装配日期：首次出现于 2026-01-01 之后的物种只有 **16 个**：Lactuca serriola/saligna/virosa（2026-02-06/09/09，距权重发布日不足 5 天，**排除**）、Chenopodium album（2026-03-13）、Medicago sativa（2026-07-09）、Dryas octopetala（07-24）、Adiantum capillus-veneris（07-27）、Akebia trifoliata（07-28）、Coptis chinensis（07-29）、Eschscholzia californica（07-30）、Cibotium barometz（07-31）、Azolla caroliniana（08-13）、Dipteris shenzhenensis（08-19）、Lygodium microphyllum（08-24，`.2` 版）、Nelumbo lutea（09-24）、Hakea prostrata（10-01）。
3. 另有 ≥ 25 个"装配/注释版本更新"（如 Amborella、Theobroma cacao、Phoenix dactylifera、Triticum urartu、Carica papaya 等），其物种早已在 RefSeq，**只算 AH，不进主面板**。
4. RefSeq 中**没有任何物种级时间留出的单子叶**（所有 2026 年新增 GCF 都是已有物种的新版本）。单子叶只能走 `SH`（RefSeq 之外的 GenBank 注释）：对 taxon 33090 的 GenBank 有注释的染色体级装配（2022 年起，398 个）逐个去除了所有在上述暴露名单里的物种，剩 216 个；其中单子叶 42 条记录 → 按谱系/大小/小基因组优先选出 Zizania、Themeda、Juncus、Rhynchospora、Acorus。
5. Ensembl Plants（release 116 / e!63，2026-06；267 个基因组、130 个物种，`https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/summary.txt`）基本是已有注释的二次分发，且其物种几乎都在上述训练名单内（new_genomes 清单为燕麦属 Avena 25 份、棉/莴苣/睡莲等已有物种的新装配）。**不作为独立来源**，只作镜像下载备选。
6. Phytozome：下载需要 JGI 登录，且 Tiberius 补充说明未发布数据受 Fort Lauderdale 约束（Tiberius 补充表 S3 脚注）；Phytozome 数据政策页是 JS 应用，**本次无法读取，许可条款[未核实]**。**不推荐**作主来源。
7. NCBI 数据使用条款：NCBI 自身对数据使用/再分发"无限制"，但提交者可能主张知识产权（`https://www.ncbi.nlm.nih.gov/home/about/policies/` 的 Molecular Data Usage 一节，已读取）。

### 2.2 主面板（12 个物种）

分层：`TH`=严格时间留出（物种级新、RefSeq 2026）；`SH`=来源留出（RefSeq 之外）；`AH`=版本留出。大小为 NCBI 报告的装配总长（含少量细胞器，核序列范围需预注册，见 §3.5）。系统发育距离列的 Mya 是**近似值，未用 TimeTree 逐对核实（TimeTree 页面 API 无法程序化访问，Tiberius 补充 S32 给出了 E. californica 距其最近训练物种 108 Mya，已核实）**。

| # | 物种 | Accession（装配名） | 层 | 装配（日期，水平，技术） | 注释（来源，日期，证据） | 大小/PCG/BUSCO | 与拟南芥/水稻的最近共同类群（NCBI）与近似分歧 | 难度标签 |
|---|---|---|---|---|---|---|---|---|
| 1 | **Dryas octopetala**（山地仙女木，蔷薇科） | `GCF_963921425.1`（drDryOcto1.1） | TH | 2024-01-10，染色体 9 条，PacBio+Arima（Darwin Tree of Life） | RefSeq EGAP（Gnomon+RNA-seq/蛋白比对）`RS_2026_07`，2026-07-24 | 230 Mb / 21,617 / 98.4（eudicots_odb10） | rosids / Mesangiospermae；与拟南芥≈100–110 Mya，与水稻≈150–200 Mya [近似] | **小基因组**，成本最低 |
| 2 | **Eschscholzia californica**（花菱草，罂粟科） | `GCF_041429995.1`（dmEscCali1.0.p） | TH | 2024-08-22，Scaffold（N50 66 Mb），PacBio+Illumina | RefSeq EGAP `RS_2026_07`，2026-07-30 | 401 Mb / 26,729 / 99.1（embryophyta_odb10） | 基部真双子叶（毛茛目）；≈125–135 Mya 与核心真双子叶分开 [近似]；距 Tiberius 最近训练物种 108 Mya（Tiberius S32） | 基部真双子叶；Tiberius 测试物种 |
| 3 | **Akebia trifoliata**（木通，木通科） | `GCF_017979445.1`（ASM1797944v1） | TH | 2021-04-16，染色体 16 条，PacBio | RefSeq EGAP `RS_2026_07`，2026-07-28 | 653 Mb / 25,221 / 97.8（embryophyta_odb10） | 毛茛目；同上 | 基部真双子叶，第二目内物种 |
| 4 | **Hakea prostrata**（山龙眼科） | `GCF_041222555.1`（UWA_Hpros_2） | TH | 2025-02-18，染色体 11 条，Nanopore+Illumina | RefSeq EGAP `RS_2026_10`，2026-10-01（**仅 15 天前，发布后再变动风险**，须固定 release 名） | 712 Mb / 24,188 / 98.6（embryophyta_odb10） | 山龙眼目（Proteales）；≈120–125 Mya [近似] | 基部真双子叶；与 Carbon-A 的 Macadamia/Telopea 同科 |
| 5 | **Chenopodium album**（藜，苋科；**六倍体**） | `GCF_948465745.1`（dcCheAlbu1.1） | TH（有同属泄漏） | 2023-02-03，染色体 27 条（含 3 个亚基因组），PacBio+Arima | RefSeq EGAP `RS_2026_03`，2026-03-13 | 1,594 Mb / **61,810** / 97.2（eudicots_odb10） | 石竹目（Caryophyllales），≈110–120 Mya [近似] | **困难：多倍体**；同属 C. quinoa 为 Carbon-A/ANNEVO/PlantCAD2 暴露 |
| 6 | **Zizania latifolia**（茭白，禾本科稻族） | `GCA_048537475.1`（ASM4853747v1） | SH | 2025-03-10，染色体 17 条，PacBio | GenBank 提交注释（武汉植物园），2025-03-10；**证据类型[未核实]** | 596 Mb / 36,280 / BUSCO 未提供（须自测） | 与水稻同族 Oryzeae，≈20–30 Mya [近似]；与拟南芥≈150–200 Mya | **近缘（对水稻）**，用于"近开发物种"层 |
| 7 | **Themeda triandra**（红燕麦草，禾本科） | `GCA_057662145.1`（ThtriWBW_v1.0） | SH | 2026-05-20，染色体 10 条，PacBio | GenBank 提交注释（塔斯马尼亚大学），2026-05-21；**证据类型[未核实]**；PCG 43,880 提示可能含多倍体/重复基因 | 755 Mb / 43,880 / BUSCO 未提供 | 与高粱/玉米同族 Andropogoneae，≈15–20 Mya [近似] | 禾本科第二个，距水稻 ≈50 Mya [近似] |
| 8 | **Juncus effusus**（灯心草） | `GCA_027726005.1`（MPIPZ_junEff1m_1） | SH | 2023-01-12，染色体 21 条，PacBio HiFi（MPIPZ，Hofstatter 等 Cell 2022；DOI 10.1016/j.cell.2022.06.045，仅读摘要） | GenBank 提交注释，2023-01-12；证据类型[未核实] | 258 Mb / 18,942 / BUSCO 未提供 | 禾本目内非禾本科，≈100–110 Mya 与水稻 [近似] | **非禾本科单子叶，小基因组** |
| 9 | **Rhynchospora breviuscula**（莎草科） | `GCA_027562975.1`（MPIPZ_rhyBre1m_1.0） | SH | 2023-01-04，染色体 5 条（全着丝粒），PacBio HiFi | GenBank 提交注释，2023-01-04；证据类型[未核实] | 421 Mb / 24,339 / BUSCO 未提供 | 禾本目；与 Carbon-A 的 Carex rostrata（RefSeq 2025）同亚科 | **非禾本科单子叶；全着丝粒（Hi-C 染色体，结构特殊）** |
| 10 | **Acorus gramineus**（石菖蒲） | `GCA_030737835.1`（ASM3073783v1） | SH | 2023-08-16，染色体 11 条，PacBio+Illumina | GenBank 提交注释（福建农林大学），2023-08-16；证据类型[未核实]；同属 A. tatarinowii 有 Nat Plants 2022（DOI 10.1038/s41477-022-01187-x）与 Sci Data 2025（DOI 10.1038/s41597-025-05933-6，基因模型 BUSCO 97.2%，仅读摘要） | 399 Mb / 24,974 / BUSCO 未提供 | 单子叶最早分支（Acorales），≈130–140 Mya 与水稻 [近似] | **早期单子叶**；Tiberius 训练含 Acorus americanus（同属，须披露） |
| 11 | **Dipteris shenzhenensis**（蕨） | `GCF_047967505.1`（Dip_SZ_1.0） | TH | 2025-02-24，染色体 33 条，PacBio | RefSeq EGAP `RS_2026_08`，2026-08-19 | 1,899 Mb / 21,776 / 90.1（embryophyta_odb10） | 蕨类 vs 种子植物 ≈380–400 Mya [近似] | **困难：蕨类（全谱系留出，大内含子）**；仅 Carbon-A/Helixer/我方在适用范围内 |
| 12 | **Amborella trichopoda**（无油樟） | `GCF_054954075.1`（AtrichopodaSantaCruz75HAP2 v2.0） | **AH（单列）** | 2026-02-09，染色体 13 条，PacBio HiFi（JGI） | RefSeq EGAP `RS_2026_09`，2026-09-10（另有 JGI 注释 GCA_054954075.1：21,711 PCG） | 700 Mb / 16,734 / 98.2（embryophyta_odb10） | ANA 级基部被子植物，≈180–200 Mya [近似] | 基部被子植物；**物种已在 Carbon-A 训练内，只作版本留出/敏感性** |

覆盖检查：真双子叶 5 个"目/大类"（蔷薇目 Rosales、毛茛目 ×2、山龙眼目 Proteales、石竹目 Caryophyllales）；单子叶 5 个（禾本科 2、灯心草科、莎草科、菖蒲目）；基部被子植物 1（Amborella，AH）；蕨类 1；规模梯度 230 Mb–1.9 Gb（<500 Mb：Dryas、Eschscholzia、Juncus、Rhynchospora、Acorus，共 5 个）。

**主面板内的三个子集（论文里分别报告）**
- **TH 集**（物种级时间留出，对 Carbon-A 证明可用）：Dryas、Eschscholzia、Akebia、Hakea（核心 4 个）+ Chenopodium album、Dipteris（困难 2 个）。
- **SH 集**（来源留出，物种级未见，但非时间留出）：Zizania、Themeda、Juncus、Rhynchospora、Acorus。
- **AH 集**：Amborella，只作敏感性分析。

### 2.3 备选面板（15 个）

| 物种 | Accession | 层 | 大小/PCG/BUSCO | 取舍说明 |
|---|---|---|---|---|
| Coptis chinensis（毛茛科） | `GCF_015680905.1` | TH | 936 Mb / 32,840 / 95.9 | 毛茛目第三个，Tiberius 有同科 Aquilegia；Akebia 的替代 |
| Nelumbo lutea | `GCF_056320005.1` | TH（**同属泄漏**） | 838 Mb / 22,264 / 98.5 | 同属 N. nucifera 在 Carbon-A 与 ANNEVO 论文版训练内，仅作"同属泄漏敏感性"对照 |
| Medicago sativa（四倍体苜蓿） | `GCF_051527215.1` | TH（**同属泄漏**） | 776 Mb / 39,275 / 97.5 | 同属 M. truncatula 为 Carbon-A/Helixer/ANNEVO/PGA 训练，Tiberius 为测试；多倍体敏感性对照 |
| Azolla caroliniana（水蕨科） | `GCF_043228905.1` | TH | 523 Mb / 18,042 / **78.9** | 小蕨类，但参考注释 BUSCO 低，Scaffold 级，作低成本蕨类备选 |
| Adiantum capillus-veneris（铁线蕨） | `GCF_014529385.1` | TH | 4,823 Mb / 31,630 / 86.7 | 超大蕨类，仅可用染色体子集 |
| Lygodium microphyllum（海金沙） | `GCF_053477565.2` | TH | 4,748 Mb / 22,911 / 87.9 | 同上；`.2` 版需核对 `.1` 首次出现日期（[NCBI] 版本列表未见更早版本） |
| Cibotium barometz（金毛狗） | `GCF_041138505.1` | TH | 3,490 Mb / 32,172 / 83.8 | 同上 |
| Canna indica（美人蕉） | `GCA_034359265.1` | SH | 752 Mb / 29,421 / 未提供 | 姜目，Carbon-A 有同目 Curcuma/Zingiber/Musa |
| Stemona tuberosa（百部，露兜树目） | `GCA_055697745.1`（单倍型 2：`GCA_055697765.1`） | SH | 803 Mb / 25,211 / 未提供 | 2026-03 T2T 装配；露兜树目无任何方法训练；注释证据未核实 |
| Paspalum saurae（雀稗，Panicoideae） | `GCA_036689595.1` | SH | 557 Mb / 45,090 / 未提供 | PlantCAD2 预训练含同属 P. vaginatum；PCG 数偏多，注释证据不明 |
| Aristolochia fimbriata（马兜铃，木兰类） | `GCA_019845555.1` | SH | 257 Mb / 21,838 / 未提供 | 木兰类基部代表，小基因组；Carbon-A 有同属 A. californica（RefSeq） |
| Nymphaea thermarum（睡莲） | `GCA_011799765.1` | SH（同属泄漏） | 368 Mb / 25,461 / 未提供 | Scaffold 级（N50 275 kb，仅 Illumina），质量差，不推荐 |
| Sphagnum compactum（泥炭藓） | `GCA_965153455.1` | SH | 381 Mb / 21,050 / 未提供 | 苔藓类备选（ELIXIR Norway 注释，流程[未核实]）；Carbon-A 有 Physcomitrium |
| Ceratodon purpureus（小金发藓） | `GCA_014871385.1` | SH | 362 Mb / 31,482 / 未提供 | JGI 注释，苔藓备选 |
| Edaphochlamys debaryana（绿藻） | `GCA_016858145.1` | SH | 142 Mb / 19,227 / 未提供 | **注释由 BRAKER2 生成（Tiberius 补充 S2 说明），与 ab initio 方法存在循环风险**，Tiberius 报告该物种基因 F1 仅约 39–41%；只作可选探索 |

**未能纳入的类群**
- **裸子植物**：RefSeq 里只有 Cryptomeria japonica（2023，暴露于 Carbon-A/ANNEVO）；GenBank 有注释的 Pinus longaeva（23.8 Gb）、Taxus chinensis（10.2 Gb）均 >10 Gb，我方推理约 160–390 h，不可行；若必须，仅取 1–2 条染色体（按染色体编号事前注册）。
- **绿藻**：RefSeq 的 Chlorophyta 物种全部在 Carbon-A 快照内；Tiberius `chlorophyta` 测试物种之一（Bathycoccus）也在内。只有上表 Edaphochlamys 一个候选，且注释有 BRAKER2 偏置，**推荐不作主表**。
- **asterids（菊类）**：主面板缺少物种级时间留出的菊类（RefSeq 2026 新增的 Lactuca 三种距权重发布不足 5 天且与 L. sativa 同属）。GenBank 216 个候选里有若干菊类/其他科，**本次没有逐个核实注释证据**，列为后续筛选项。

### 2.4 主面板的资格门槛（冻结前须通过，写入预注册）
1. 物种与近缘暴露表（§1.7）更新到冻结日；若冻结日前出现新训练名单（如 ANNEVO 现行权重清单、Tiberius Figshare 清单）发现重叠，替换为备选，**替换规则事前固定**（同层同大小级的下一个）。
2. 对 SH 层和所有 GenBank 注释：用 compleasm/BUSCO（odb12 与 odb10 均测）测注释蛋白的完整度，**阈值 ≥ 90% 且不低于同谱系 RefSeq 注释的中位数 −5 点**；未达阈值者换备选。同时用 SRA 公开 RNA-seq（若有）做内含子支持率（§4.2）。
3. 每个 RefSeq 物种固定 annotation release 名（如 `GCF_963921425.1-RS_2026_07`）与 GFF 文件 md5，防止后续重注释。
4. 面板 accession、GFF md5、核序列 contig 范围清单在**看任何预测结果之前**冻结（§3.5）。

---

## 3. 同源隔离方案

### 3.1 三类暴露集合
- **E_ours**：我方头/微调的训练物种（现为拟南芥、水稻；若扩训则加入，同时从面板中移出）。
- **E_base(X)**：基线 X 的训练+验证+预训练物种（§1 和附录 A 的名单），Carbon-A 用 RefSeq 反推集合（201 个物种）。
- **E_union**：以上全部之并集。

### 3.2 蛋白同源审计（对每个面板物种逐一做）
1. 参考侧：面板物种的参考 CDS 翻译蛋白（最长转录本一条/基因；用面板冻结 GFF + FASTA 提取，或直接用 NCBI `_protein.faa.gz`）。
2. 训练侧：E_union 内各物种的蛋白组（RefSeq 的 `_protein.faa.gz` 约 5–19 MB/物种，201 个物种约 1.5–2 GB；Phytozome/GenBank 来源用其蛋白文件）。
3. 搜索：**MMseqs2 easy-search**（`-s 7.5 -e 1e-5 --max-seqs 50 --cov-mode 2 -c 0.5`）作为主检索，**DIAMOND `--very-sensitive`（e ≤ 1e-5）** 在随机 10% 基因子集上做一致性复核；输出最佳命中的序列一致性 id_max、查询覆盖 cov。
4. 分档（以查询基因为单位）：

   | 档 | 条件 | 解释 |
   |---|---|---|
   | H（同源强） | id_max ≥ 70% 且 cov ≥ 70% | 视为训练集中有直系同源，等同"见过" |
   | M | 40% ≤ id_max < 70%，cov ≥ 50% | 近缘同源 |
   | L | 30% ≤ id_max < 40%，e ≤ 1e-5 | 远缘同源 |
   | N | 无命中 | 家族级新颖 |
5. 基因家族拆分：对 面板蛋白 ∪ E_union 蛋白做 `mmseqs cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0`（或 OrthoFinder），记每个面板基因所在簇是否含 E_ours/E_base(X) 成员；报告"家族内已见（含训练成员）"与"家族新颖（簇内无训练成员）"两个子集。
6. **按 X 分别审计**：对 Carbon-A、Helixer、Tiberius、ANNEVO、PlantGeneAnn 各用其 E_base(X)；对我方用 E_ours 与 E_ours∪Carbon-A（Carbon-A 暴露也要计，因为骨干同权重）。
7. 基因组级复核：面板基因组与 E_union 中同属/同科基因组做 Mash（`-s 10^6`）/ 或 minimap2 asm5 的覆盖度，标记 ANI > 95% 的（如 Amborella 新旧装配，Nelumbo lutea 与 N. nucifera）；这类物种自动转入敏感性层。

### 3.3 报告口径（论文里并列两张主表 + 一张敏感性表）
- **表 A：时间留出口径（对 Carbon-A）**：只含 TH 核心 4 个（Dryas、Eschscholzia、Akebia、Hakea）+（困难）Chenopodium、Dipteris；表头明确"物种级未见，GCF 首次出现于 2026-03 之后"；全部 5–6 个基线同表。
- **表 B：物种留出口径（对全部方法，保守）**：TH ∪ SH 中**物种本身不在任何方法训练/验证名单里**的物种（即主面板去掉 Amborella），按 NCBI 分类标注各方法的最近已见近缘（§1.7），再按同源档 H/M/L/N 给出 F1 的**分层曲线**。
- **表 C：同源剥离口径**：去掉 id_max ≥ 50% 且 cov ≥ 50% 的参考基因，以及预测蛋白对 E_union 同样条件命中的预测基因，两侧同时剥离，**对所有方法使用同一并集剥离集合**（避免对某个方法有利）。
- **敏感性：** Amborella（AH）、Nelumbo lutea、Medicago sativa（同属泄漏）、Chenopodium album 的 F1 与其"同属暴露方法"的相对排序是否改变。
- 论文必须披露：RefSeq 标签风格偏置（Carbon-A 训练标签同源风格）、ANNEVO 现行权重未知、Tiberius 未公开训练基因组名单（Figshare）核对结果、GENERanno 预训练暴露对这份权重不适用的依据（Table 4 随机初始化）。

### 3.4 近缘不可避免的说明
Tiberius 自己的基准里"除 D. melanogaster 外测试物种与最近训练物种在科级以上"，分歧 14–281 Mya（Tiberius v3 §3.1）。本面板的 Dryas（同科）、Eschscholzia/Akebia（同目）与此一致；Zizania/Themeda（同族）更近，因此单独标成"近缘层"，不与"远缘层"混报平均值。

### 3.5 核序列范围与短序列
- RefSeq 装配里有叶绿体/线粒体（如 Dryas、Hakea 等），Primary Assembly 之外的细胞器序列不进入核基因评测；contig 范围用 `sequence_reports` 在冻结时导出（`/genome/accession/<acc>/sequence_reports`），**不按长度、基因密度或预测表现删序列**。
- Eschscholzia 与 Azolla 为 Scaffold 级（Eschscholzia scaffold N50 66 Mb；Azolla 1,904 contigs），含大量短序列；M28 已有"短 contig 预检"记录，按同一规则保留全部输入。

---

## 4. 评测协议

### 4.1 指标
| 层 | 指标 | 说明 |
|---|---|---|
| **主指标** | exact CDS chain F1（基因/位点级）：完整 CDS 外显子链与参考的某条转录本精确一致才算对；参考分"最长转录本"与"任一转录本"两种口径（ANNEVO Suppl. Table 22 的做法，任一转录本会抬高分数） | 与项目既有 ruler 一致，并同时给 `gffcompare -r ref --no-merge -e 3 --strict-match`（Tiberius 口径，允许起止端 ±3 bp）与 `gffcompare --no-exon-merge --strict-match`（ANNEVO 口径）；gffcompare v0.12.10（两份文献都用此版本） |
| 次指标 | CDS 外显子级 P/R/F1；内含子级、intron-chain 级；核苷酸级 CDS F1（含相位 phase F1，Helixer 口径） | |
| 完整度 | BUSCO（compleasm 或 BUSCO≥5）对预测蛋白：odb12 与 odb10 并报；谱系：eudicots / liliopsida（或 poales）/ embryophyta；同时报参考注释的 BUSCO 作为天花板（ANNEVO 做法：318/793 个物种 BUSCO 高于参考） | Carbon-A 用 odb12（[PDF] Table S4） |
| 蛋白水平 | 预测蛋白对 SwissProt/UniRef90 的 DIAMOND 命中率（e ≤ 1e-10），用于"无参考重叠的新颖预测"（Tiberius v3 的 novel-loci hit-rate）；参考注释外预测占比 | |
| 误差分解 | 融合/断裂基因率、按 CDS 长度/外显子数分层 F1、假阳性来源（碎片 mRNA、与非编码转录本重叠，ANNEVO 报告 60.6% / 18.3%） | |
| 可选 RNA-seq | 内含子剪接读段支持率（STAR/HISAT2 比对 + 唯一映射 ≥ N 条跨接读段）；预测基因 CDS 被 StringTie 转录本覆盖比例；对"参考外预测"区分有无转录证据。**仅对 SRA 有公共 RNA-seq 的物种做**，BAM 数据量大（每物种数十 GB 原始数据），先用 SRA 元数据估计再决定，**本次未核实各物种 RNA-seq 可得性** | |
| 同源分层 | §3.3 的 H/M/L/N 档 F1；同源剥离 F1 | |
| 置信度（C3 用） | 基因级 ECE/Brier/可靠性曲线（按物种分层），与 Carbon-A 置信分（阈值 0.10）同召回比精度 | 与 07 号文件 C3 对应 |

统计：每物种按 1 Mb 区块 bootstrap 基因（≥1000 次）给 CI；面板级报告**物种中位数**、均值、最差物种；配对检验用物种为单位的 Wilcoxon 符号秩 + 符号检验；我方 3 个训练种子分别报告并以种子均值进入配对检验；多重比较（对 4–5 个基线）用 Holm 校正。预注册：主终点 = 表 A 和表 B 上"我方 vs 每个基线的 exact CDS chain F1 物种中位数差"；Nature Methods 级声明所需（建议）：对 Carbon-A（未过滤与 ≥0.10 过滤两个工作点都比）的中位数差 > 0，CI 不含 0，且在 ≥ 2/3 物种上胜出；对 Helixer/Tiberius/ANNEVO 同。

### 4.2 基线工具与版本（均须在 Baobab/ont 上冻结 revision，并在论文里写明）

| 基线 | 版本/资产（已在项目中固定的优先） | 正确模型 | 适用范围与注意 |
|---|---|---|---|
| **Carbon-A** | HF `HuggingFaceBio/Carbon-A-1.2B` revision `f6789ee79e6ef397de333577ffec1cb734496620`（Baobab `refs/carbon-a/…`，Apache-2.0） | 单一检查点；**未过滤**与**置信 ≥ 0.10 过滤**两个工作点都报 | 全谱系；本机/远端官方 `annotation_pipeline` 已齐；转移/置信网络训练数据未披露 |
| **Tiberius** | v2.0.8 隔离包（项目 `native_baseline_preparation.md` 已固定；配置 commit `871d4979568ae09e56de7eca5d7489046ee833ea`）；许可 MIT（Tiberius README，见 05 号文件） | `angiosperms`（Mesangiospermae，权重 `https://bioinf.uni-greifswald.de/bioinf/tiberius/models/angiosperms.tar.gz`）；`chlorophyta`（藻类可选） | 蕨类/苔藓/Amborella 超范围，只能标 OOD；论文用 `--seq_len` 默认 200,070（angiosperms） |
| **Helixer** | 0.3.7（d17bb496…），`land_plant_v0.3_a_0080.h5`（GPL-3.0） | `--lineage land_plant` 默认 | 含藻类/苔藓/蕨类训练（Selaginella、Marchantia 等），范围最广；训练暴露见 §1.4 |
| **ANNEVO** | 源码 42c920f7…（v2.3.x），`saved_model/ANNEVO_Magnoliopsida.pt`（blob 21487d87…，非商业许可） | `-l Magnoliopsida`，并按 `docs/boundary_aware_model.md` 加 `--num_classes 15 --boundary-aware --min_intron_length 20`（现行植物模型的推荐用法） | 只有被子植物；蕨类/苔藓/藻类无模型；当前权重训练清单未知 |
| PlantGeneAnn（可选） | HF `qzzhang/PlantGeneAnn-multi-species`（MIT），需 mamba-ssm | multi-species | 预印本；暴露见 §1.6；Zizania 等不在其训练名单 |
| GeneCAD（可选） | `plantcad/genecad`（CC BY-NC） | | 预印本；训练名单[未核实] |
| AUGUSTUS ab initio（可选） | 3.5.x，用最近缘已有物种参数（拟南芥/水稻/玉米） | 与 Helixer NM 论文同做法 | 对 AUGUSTUS 只能给近似基线；Helixer 论文指出 AUGUSTUS 仅 4/13 植物有预训练模型 |
| BRAKER3（可选，证据型上界） | 3.0.x（Genome Research 34:769，2024） | 需 RNA-seq + 蛋白证据；Tiberius 用 VARUS 生成 RNA-seq，并用 Tiberius 训练物种蛋白做证据 | 单基因组均耗时 2,170 min（72 线程，Tiberius v3）；**仅建议在 3 个小基因组（Dryas、Juncus、Eschscholzia）做** |
| **我方** | 冻结 checkpoint + 固定操作点 | 3 个种子 | 先冻结再看面板预测，caller 只读 DNA |

公平性条款：所有基线用**官方发布权重与官方推荐参数**，不对面板调参；所有方法在同一冻结 FASTA（含全部核序列）上跑；基线是否使用软屏蔽按其官方要求（Tiberius angiosperms/Helixer 不需要，Carbon-A 不需要）；`gffcompare` 在同一参考 GFF（仅蛋白编码基因、取 CDS）上评分。

---

## 5. 算力估计

假设（均来自 06 号文件与公开文献）：我方 ≈ 17 kbp/s（单张 RTX 4090，推理）；Carbon-A ≈ 15–20 GPU-min / 30 Mb；ANNEVO/Tiberius/Helixer 按 ANNEVO README（`https://github.com/xjtu-omics/ANNEVO`）的 12 个模式物种平均（单张 RTX 4090：ANNEVO 12.18 min、Tiberius 43.59 min、Helixer 286.61 min）按基因组大小线性外推，我假设这 12 个物种平均约 1.2 Gb（**依据：README 未给平均大小，该值为我据其中已知物种的估计，[未核实]**），即每 Gb：ANNEVO ≈ 10 min、Tiberius ≈ 36 min、Helixer ≈ 4 h。

| 物种 | 大小 Mb | 我方（h） | Carbon-A（h） |
|---|---:|---:|---:|
| Dryas octopetala | 230 | 3.8 | 1.9–2.6 |
| Eschscholzia californica | 401 | 6.6 | 3.3–4.5 |
| Akebia trifoliata | 653 | 10.7 | 5.4–7.3 |
| Hakea prostrata | 712 | 11.6 | 5.9–7.9 |
| Zizania latifolia | 596 | 9.7 | 5.0–6.6 |
| Themeda triandra | 755 | 12.3 | 6.3–8.4 |
| Juncus effusus | 258 | 4.2 | 2.1–2.9 |
| Rhynchospora breviuscula | 421 | 6.9 | 3.5–4.7 |
| Acorus gramineus | 399 | 6.5 | 3.3–4.4 |
| Chenopodium album | 1,594 | 26.0 | 13.3–17.7 |
| Dipteris shenzhenensis | 1,899 | 31.0 | 15.8–21.1 |
| Amborella trichopoda | 700 | 11.4 | 5.8–7.8 |
| **主面板合计** | **8,618** | **≈141** | **≈72–96** |
| 核心 10（去掉 Chenopodium、Dipteris） | 5,125 | ≈84 | ≈43–57 |
| 5 个 <500 Mb 小基因组 | 1,709 | ≈28 | ≈14–19 |

- 其他基线（主面板 8.6 Gb，4090 估计）：ANNEVO ≈ 1.4 h；Tiberius ≈ 5 h；Helixer ≈ 34 h；三者合计 ≈ 40 h。基线要在各自范围内跑（ANNEVO/Tiberius 不跑蕨类 → 约再减 0.3–1 h）。
- 我方 3 个训练种子：**训练**成本不变于面板；**推理**若 3 个种子都要在面板全量跑，则 ≈ 3 × 141 ≈ 420 h，建议：3 种子只在 5 个小基因组 + Dryas 之外的 TH 核心 3 个（合计 ≈ 6 个物种、约 2.6 Gb）上跑，其余物种用最终选定的 1 个种子（把 3 种子推理降到约 3×43 + 1×98 ≈ 227 h；仍偏大，若紧张再把 Chenopodium/Dipteris 限制为"前 N 条染色体"，N 按染色体编号事前注册，不得按预测表现选）。
- BRAKER3（可选）三个小基因组：按 Tiberius v3 的平均 2,170 min/基因组估计约 36 h CPU/基因组 × 3，**需 RNA-seq 获取与比对额外成本，未估**。
- 重计算优先 ont；Baobab 单作业 > 50 credits 须先问用户（07 号文件规则）。
- 数据：主面板下载 2.84 GB，预测/中间文件（Carbon-A 的概率轨/特征缓存）另需估计，**未核实**。

---

## 6. 可下载性核验与下载清单

方法：对 NCBI datasets 返回的 assembly_name 构造 FTP 路径 `https://ftp.ncbi.nlm.nih.gov/genomes/all/{GCF|GCA}/xxx/xxx/xxx/{acc}_{assembly_name}/`，对 `_genomic.fna.gz`、`_genomic.gff.gz`、`_protein.faa.gz` 做 `curl -sIL`（HEAD），**全部返回 200**，大小取自 `content-length`。**没有下载任何文件**。GenBank（GCA）的 GFF 为提交者注释（NCBI 转换），**未核对 GFF 内容与 `gene`/`CDS` 层级**。

### 6.1 主面板（大小单位 MB，gz）

| 物种 | Accession | FASTA | GFF3 | 蛋白 FASTA | 基础路径（`…/{acc}_{name}/{acc}_{name}` + 后缀 `_genomic.fna.gz`/`_genomic.gff.gz`/`_protein.faa.gz`） |
|---|---|---:|---:|---:|---|
| Dryas octopetala | GCF_963921425.1 | 71.2 | 7.7 | 6.9 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/921/425/GCF_963921425.1_drDryOcto1.1/GCF_963921425.1_drDryOcto1.1` |
| Eschscholzia californica | GCF_041429995.1 | 121.1 | 11.4 | 8.3 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/041/429/995/GCF_041429995.1_dmEscCali1.0.p/GCF_041429995.1_dmEscCali1.0.p` |
| Akebia trifoliata | GCF_017979445.1 | 204.2 | 11.1 | 7.7 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/979/445/GCF_017979445.1_ASM1797944v1/GCF_017979445.1_ASM1797944v1` |
| Hakea prostrata | GCF_041222555.1 | 220.9 | 9.0 | 7.6 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/041/222/555/GCF_041222555.1_UWA_Hpros_2/GCF_041222555.1_UWA_Hpros_2` |
| Zizania latifolia | GCA_048537475.1 | 175.6 | 6.4 | 9.0 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/537/475/GCA_048537475.1_ASM4853747v1/GCA_048537475.1_ASM4853747v1` |
| Themeda triandra | GCA_057662145.1 | 221.5 | 6.3 | 8.5 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/057/662/145/GCA_057662145.1_ThtriWBW_v1.0/GCA_057662145.1_ThtriWBW_v1.0` |
| Juncus effusus | GCA_027726005.1 | 73.3 | 3.4 | 4.8 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/726/005/GCA_027726005.1_MPIPZ_junEff1m_1/GCA_027726005.1_MPIPZ_junEff1m_1` |
| Rhynchospora breviuscula | GCA_027562975.1 | 128.2 | 4.2 | 6.0 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/562/975/GCA_027562975.1_MPIPZ_rhyBre1m_1.0/GCA_027562975.1_MPIPZ_rhyBre1m_1.0` |
| Acorus gramineus | GCA_030737835.1 | 121.5 | 4.3 | 5.7 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/737/835/GCA_030737835.1_ASM3073783v1/GCA_030737835.1_ASM3073783v1` |
| Chenopodium album | GCF_948465745.1 | 487.5 | 20.1 | 18.8 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/948/465/745/GCF_948465745.1_dcCheAlbu1.1/GCF_948465745.1_dcCheAlbu1.1` |
| Dipteris shenzhenensis | GCF_047967505.1 | 589.8 | 7.9 | 7.3 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/967/505/GCF_047967505.1_Dip_SZ_1.0/GCF_047967505.1_Dip_SZ_1.0` |
| Amborella trichopoda | GCF_054954075.1 | 220.7 | 14.7 | 8.9 | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/954/075/GCF_054954075.1_AtrichopodaSantaCruz75HAP2_v2.0/GCF_054954075.1_AtrichopodaSantaCruz75HAP2_v2.0` |
| **合计** | | **2,635.6** | **106.5** | **99.5** | **总计 ≈ 2.84 GB** |

### 6.2 备选面板（HEAD 均 200，大小 MB）

| 物种 | Accession | FASTA | GFF3 | 蛋白 | 基础路径（同上规则） |
|---|---|---:|---:|---:|---|
| Coptis chinensis | GCF_015680905.1 | 294.2 | 14.8 | 8.6 | `…/GCF/015/680/905/GCF_015680905.1_ASM1568090v1/GCF_015680905.1_ASM1568090v1` |
| Nelumbo lutea | GCF_056320005.1 | 261.7 | 12.6 | 7.5 | `…/GCF/056/320/005/GCF_056320005.1_American_lotus-T2T/GCF_056320005.1_American_lotus-T2T` |
| Medicago sativa | GCF_051527215.1 | 211.4 | 15.7 | 10.9 | `…/GCF/051/527/215/GCF_051527215.1_Bolivia_mono/GCF_051527215.1_Bolivia_mono` |
| Azolla caroliniana | GCF_043228905.1 | 160.9 | 6.4 | 6.4 | `…/GCF/043/228/905/GCF_043228905.1_voAzoCaro1.0.p/GCF_043228905.1_voAzoCaro1.0.p` |
| Adiantum capillus-veneris | GCF_014529385.1 | 1,494.3 | 11.7 | 9.1 | `…/GCF/014/529/385/GCF_014529385.1_ASM1452938v2/GCF_014529385.1_ASM1452938v2` |
| Lygodium microphyllum | GCF_053477565.2 | 1,488.2 | 10.1 | 7.5 | `…/GCF/053/477/565/GCF_053477565.2_Lygodium_microphyllum_v1.1.10/GCF_053477565.2_Lygodium_microphyllum_v1.1.10` |
| Cibotium barometz | GCF_041138505.1 | 1,019.4 | 12.5 | 9.9 | `…/GCF/041/138/505/GCF_041138505.1_ASM4113850v1/GCF_041138505.1_ASM4113850v1` |
| Canna indica | GCA_034359265.1 | 225.3 | 5.5 | 7.4 | `…/GCA/034/359/265/GCA_034359265.1_ASM3435926v1/GCA_034359265.1_ASM3435926v1` |
| Stemona tuberosa（hap1） | GCA_055697745.1 | 254.0 | 4.1 | 5.5 | `…/GCA/055/697/745/GCA_055697745.1_ASM5569774v1/GCA_055697745.1_ASM5569774v1` |
| Paspalum saurae | GCA_036689595.1 | 171.6 | 8.0 | 11.1 | `…/GCA/036/689/595/GCA_036689595.1_IICAR_Pnotatum_1.0/GCA_036689595.1_IICAR_Pnotatum_1.0` |
| Aristolochia fimbriata | GCA_019845555.1 | 75.0 | 3.6 | 5.6 | `…/GCA/019/845/555/GCA_019845555.1_ASM1984555v1/GCA_019845555.1_ASM1984555v1` |
| Nymphaea thermarum | GCA_011799765.1 | 114.2 | 7.1 | 6.5 | `…/GCA/011/799/765/GCA_011799765.1_ASM1179976v1/GCA_011799765.1_ASM1179976v1` |
| Sphagnum compactum | GCA_965153455.1 | 121.1 | 5.8 | 5.9 | `…/GCA/965/153/455/GCA_965153455.1_csSphComp1.1/GCA_965153455.1_csSphComp1.1` |
| Ceratodon purpureus | GCA_014871385.1 | 112.4 | 7.2 | 7.2 | `…/GCA/014/871/385/GCA_014871385.1_CpurpureusR40_1_0/GCA_014871385.1_CpurpureusR40_1_0` |
| Edaphochlamys debaryana | GCA_016858145.1 | 42.6 | 5.5 | 7.9 | `…/GCA/016/858/145/GCA_016858145.1_ASM1685814v1/GCA_016858145.1_ASM1685814v1` |
| **备选合计（15 个）** | | **≈6,046** | **≈131** | **≈117** | **≈6.29 GB** |

（`…` = `https://ftp.ncbi.nlm.nih.gov/genomes/all`。Lactuca serriola 的 FASTA 为 777.7 MB，未入面板。）

可选下载 API：`https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<acc>/download`（需 datasets CLI/认证流程，本次未用）。许可：NCBI 自身无使用限制（见 §2.1.7），GenBank 提交者可能主张权利；Ensembl Plants 为公开 FTP（`https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/`），条款[未逐页核实]；Phytozome 条款[未核实]。

---

## 7. 风险、未核实项与待办

### 7.1 风险
1. **SH 层注释证据**：5 个物种的 GenBank 注释方法/RNA-seq/BUSCO 未核实，是否含 ab initio（含 Helixer/BRAKER）预测，决定其对各方法是否公平。门槛见 §2.4。
2. **RefSeq 标签风格偏置**：Carbon-A 训练标签就是 RefSeq（Gnomon）标签；TH 层全是 RefSeq，会系统性有利于 Carbon-A。应用 SH 层（非 RefSeq 风格注释）做对照，并在论文里分开报告。我方头也是用 RefSeq 拟南芥/水稻训练的，同样有此偏置。
3. **快照日期不确定**：Carbon-A 快照"≥ 2025-04 至 ≤ 2026-02-10"之间首次出现的 GCF 暴露未知；本面板的 TH 物种均晚于 2026-03，不受影响；但如果以后想加入 2025 年新 GCF，必须先向作者确认快照日期。
4. **同属/同科近缘**：Chenopodium album（C. quinoa）、Nelumbo lutea、Medicago sativa 等同属泄漏；Dryas/Eschscholzia/Hakea 同科或同亚科近缘，无法避免，靠分层报告处理。
5. **Tiberius Figshare 名单与 ANNEVO 现行训练名单未核对**。
6. **RefSeq 注释可能被更新**：Hakea（2026-10-01）等新注释，固定 release 名与 md5。
7. **蕨类参考质量**：Dipteris BUSCO 90.1、Azolla 78.9、其他蕨类 83–88%，参考注释本身就不完整；蕨类 F1 的上限受限，不建议作为"超越"的主证据，仅作 OOD 稳健性。
8. **基线适用范围**：Tiberius/ANNEVO 对蕨类、Amborella 超出训练范围；公平主表只在被子植物 10 个物种上并列；蕨类只比 Carbon-A/Helixer/我方。

### 7.2 复现方法（本次所用的查询）
- 全部 RefSeq 植物含历史版本：`GET https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/33090/dataset_report?filters.assembly_source=refseq&filters.assembly_version=all_assemblies&page_size=100`（分页取 `next_page_token`），字段 `assembly_info.release_date`、`annotation_info.release_date/name/busco/stats.gene_counts`。
- GenBank 候选：`…/taxon/33090/dataset_report?filters.assembly_source=genbank&filters.has_annotation=true&filters.assembly_level=chromosome&filters.assembly_level=complete_genome&filters.first_release_date=2022-01-01T00:00:00.000Z`（398 条）。
- 分类学最近共同祖先：`…/taxonomy/taxon/<names>` 取 `lineage`。
- 中间文件（会话临时目录，非项目产物，会话结束后丢失）：`…/scratchpad/work/{refseq_plants_all.json,gb_plants.json,exposure_sets.json,lca_table.json,dl_check*.json}`；如需长期保存，应由用户决定放到项目 `reports/`。

### 7.3 待办（按优先级）
1. 对 5 个 SH 物种做注释资格测试（BUSCO、RNA-seq 可得性、注释方法文献），不合格则换 Canna/Stemona/Paspalum/Aristolochia。
2. 打开 Tiberius Figshare `10.6084/m9.figshare.32086728`，核对面板物种。
3. 向作者（或看 HF 仓库更新）确认 Carbon-A 快照日期与辅助网络训练集；确认 ANNEVO 现行权重训练清单。
4. 用 TimeTree 5 逐对计算面板物种到拟南芥/水稻/最近已见近缘的分歧时间，替换本文的近似值。
5. 在 Baobab 上建立 E_union 蛋白库与 MMseqs2 审计流水线（CPU），在看任何预测之前完成并冻结。
6. 补一个菊类（asterid）物种：从 GenBank 216 个候选里筛选，需核实注释证据。
7. 冻结面板 accession/GFF md5/核序列 contig 范围，在登记文件里写明替换规则。

---

## 附录 A：各方法训练/验证/预训练植物物种名单（按本次读取的来源）

**Helixer a0080 训练（S23，seed_009_017，32 个）**：Ananas comosus, Arabidopsis thaliana, Brachypodium distachyon, Cicer arietinum, Capsella grandiflora, Chlamydomonas reinhardtii, Capsella rubella, Eutrema salsugineum, Fragaria vesca, Glycine max, Lotus japonicus, Linum usitatissimum, Musa acuminata, Manihot esculenta, Mimulus guttatus, Marchantia polymorpha, Micromonas pusilla CCMP1545, Medicago truncatula, Nymphaea colorata, Ostreococcus lucimarinus, Oryza sativa, Oropetium thomaeum, Phaseolus acutifolius, Prunus persica, Porphyra umbilicalis, Sorghum bicolor, Solanum lycopersicum, Schrenkiella parvula, Triticum aestivum, Trifolium pratense, Volvox carteri, Zea mays。

**Helixer S18 train_val 全集（77，含上述 32 个与验证物种）**：Amaranthus hypochondriacus; Amborella trichopoda; Ananas comosus; Arabidopsis lyrata; Arabidopsis thaliana; Asparagus officinalis; Beta vulgaris; Brachypodium distachyon; Brachypodium hybridum; Brassica oleracea; Brassica rapa; Capsella grandiflora; Capsella rubella; Carica papaya; Chenopodium quinoa; Chlamydomonas reinhardtii; Chromochloris zofingiensis; Cicer arietinum; Cinnamomum kanehirae; Citrus clementina; Citrus sinensis; Coccomyxa subellipsoidea; Cucumis sativus; Daucus carota; Dioscorea alata; Dunaliella salina; Eucalyptus grandis; Eutrema salsugineum; Fragaria vesca; Glycine max; Glycine soja; Gossypium raimondii; Helianthus annuus; Hordeum vulgare; Kalanchoe fedtschenkoi; Lactuca sativa; Linum usitatissimum; Lotus japonicus; Lupinus albus; Malus domestica; Manihot esculenta; Marchantia polymorpha; Medicago truncatula; Micromonas pusilla; Micromonas sp. RCC299; Mimulus guttatus; Musa acuminata; Nymphaea colorata; Olea europaea; Oropetium thomaeum; Oryza sativa; Ostreococcus lucimarinus; Panicum hallii; Panicum virgatum; Phaseolus acutifolius; Physcomitrella patens; Poncirus trifoliata; Populus trichocarpa; Porphyra umbilicalis; Prunus persica; Ricinus communis; Salix purpurea; Schrenkiella parvula; Selaginella moellendorffii; Setaria italica; Solanum lycopersicum; Solanum tuberosum; Sorghum bicolor; Spirodela polyrhiza; Theobroma cacao; Trifolium pratense; Triticum aestivum; Vigna unguiculata; Vitis vinifera; Volvox carteri; Zea mays; Zostera marina。

**Helixer 测试（13）**：Arachis hypogaea; Brassica napus; Cannabis sativa; Coffea arabica; Hibiscus syriacus; Nicotiana attenuata; Oryza brachyantha; Papaver somniferum; Phoenix dactylifera; Setaria viridis; Solanum pennellii; Triticum dicoccoides; Vitis riparia。

**Tiberius angiosperms 训练（45）**：Acorus americanus; Agave tequilana; Aquilegia coerulea; Arundo donax; Beta vulgaris; Brassica rapa; Brocchinia micrantha; Capparis spinosa; Chasmanthium laxum; Citrus limon; Cleome violacea; Coffea arabica; Corylus avellana; Daucus carota; Dendrobium catenatum; Dioscorea alata; Eucalyptus grandis; Fragaria vesca; Glycine max; Gossypium hirsutum; Hydrocotyle leucocephala; Joinvillea ascendens; Lactuca sativa; Lindenbergia philippensis; Malus domestica; Manihot esculenta; Mentha longifolia; Miscanthus sinensis; Moringa oleifera; Musa acuminata; Myrothamnus flabellifolia; Nicotiana tabacum; Oryza sativa; Osyris compressa; Panicum virgatum; Pharus latifolius; Phoenix dactylifera; Setaria viridis; Sorghum bicolor; Thinopyrum intermedium; Ulmus americana; Vitis vinifera; Yucca filamentosa; Zea mays; Zingiber officinale。验证：Populus trichocarpa、Typha latifolia。测试：Arabidopsis thaliana、Brachypodium distachyon、Eschscholzia californica、Freycinetia multiflora、Medicago truncatula、Mimulus guttatus、Urochloa brizantha。（来源：Tiberius v3 补充 PDF 表格；注意固定 `angiosperms.yaml` 的 training_species 用版本化名称，如 Lsativa.LsPKU06、Mdomestica WA_38，对应上述物种，另含 Mlongifolia、Cspinosa、Dendrobium_catenatum_GCF_001605985 等；yaml 版本 2026-04-15，比补充表可能多或少少数物种，**以 yaml 为准做最终审计**。）

**ANNEVO 论文版 Embryophyta（Table 24）**
- Train（40）：Alnus glutinosa; Arabidopsis lyrata; Arachis hypogaea; Asparagus officinalis; Camelina sativa; Capsella rubella; Carya illinoinensis; Coffea eugenioides; Cucumis sativus; Cucurbita maxima; Cucurbita moschata; Euphorbia lathyris; Gossypium hirsutum; Hevea brasiliensis; Ipomoea nil; Ipomoea triloba; Lycium barbarum; Mangifera indica; Manihot esculenta; Medicago truncatula; Musa acuminata; Nelumbo nucifera; Oryza brachyantha; Oryza sativa; Panicum hallii; Panicum virgatum; Phragmites australis; Physcomitrium patens; Pisum sativum; Prunus persica; Ricinus communis; Rosa chinensis; Solanum dulcamara; Solanum lycopersicum; Syzygium oleosum; Theobroma cacao; Triticum aestivum; Vicia villosa; Vitis vinifera; Ziziphus jujuba。
- Val & Test（9）：Brachypodium distachyon; Chenopodium quinoa; Daucus carota; Humulus lupulus; Lotus japonicus; Mercurialis annua; Rhododendron vialii; Salvia splendens; Sorghum bicolor。
- Test（119）：见补充表 Table 24（含 Arabidopsis thaliana、Amborella trichopoda、Papaver somniferum、Setaria viridis、Setaria italica 等）。

**PlantGeneAnn 微调训练（40）**：Alnus glutinosa; Arabidopsis lyrata; Arachis hypogaea; Brachypodium distachyon; Camelina sativa; Capsella rubella; Coffea arabica; Cucumis sativus; Daucus carota; Euphorbia lathyris; Glycine max; Gossypium hirsutum; Ipomoea triloba; Juglans regia; Lactuca sativa; Lycium barbarum; Mangifera indica; Manihot esculenta; Medicago truncatula; Musa acuminata; Nymphaea colorata; Oryza brachyantha; Oryza glaberrima; Panicum virgatum; Phragmites australis; Populus trichocarpa; Prunus persica; Ricinus communis; Rosa chinensis; Salvia splendens; Setaria viridis; Solanum dulcamara; Solanum lycopersicum; Syzygium oleosum; Theobroma cacao; Trifolium pratense; Triticum aestivum; Vitis vinifera; Zea mays; Ziziphus jujuba。验证：Brassica oleracea、Chenopodium quinoa、Henckelia pumila、Hibiscus syriacus、Lolium perenne、Lotus japonicus、Solanum tuberosum、Sorghum bicolor。

**PlantCAD2 预训练（65）**：Amaranthus hypochondriacus; Ananas comosus; Aquilegia coerulea; Arabidopsis thaliana; Arachis hypogaea; Asparagus officinalis; Beta vulgaris; Betula platyphylla; Boechera stricta; Brachypodium distachyon; Capsella rubella; Carya illinoinensis; Chenopodium quinoa; Cicer arietinum; Cinnamomum kanehirae; Citrus clementina; Coffea arabica; Corylus avellana; Corymbia citriodora; Cucumis sativus; Daucus carota; Dioscorea alata; Eleusine coracana; Eucalyptus grandis; Eutrema salsugineum; Fragaria vesca; Glycine max; Helianthus annuus; Hordeum vulgare; Lactuca sativa; Lens ervoides; Linum usitatissimum; Lotus japonicus; Lupinus albus; Malus domestica; Manihot esculenta; Medicago truncatula; Miscanthus sinensis; Musa acuminata; Nymphaea colorata; Olea europaea; Oropetium thomaeum; Oryza sativa; Panicum virgatum; Paspalum vaginatum; Poncirus trifoliata; Populus trichocarpa; Portulaca amilis; Prunus persica; Quercus rubra; Ricinus communis; Salix purpurea; Schrenkiella parvula; Setaria viridis; Sorghum bicolor; Spinacia oleracea; Spirodela polyrhiza; Thinopyrum intermedium; Triticum aestivum; Vaccinium darrowii; Vicia faba; Vigna unguiculata; Vitis vinifera; Zea mays; Zostera marina。

## 附录 B：Carbon-A 快照内的 RefSeq 植物物种（按 NCBI 当前数据反推，首个 GCF 版本日期 < 2026-02-10，共 201 个；**[推断]，真实快照为 182 个物种、日期未披露**）

Abrus precatorius; Actinidia eriantha; Aegilops tauschii; Alnus glutinosa; Amaranthus tricolor; Amborella trichopoda; Ananas comosus; Andrographis paniculata; Apium graveolens; Arabidopsis lyrata; Arabidopsis thaliana; Arachis duranensis; Arachis hypogaea; Arachis ipaensis; Arachis stenosperma; Argentina anserina; Aristolochia californica; Asparagus officinalis; Auxenochlorella protothecoides; Bathycoccus prasinos; Benincasa hispida; Beta vulgaris; Bidens hawaiensis; Brachypodium distachyon; Brassica napus; Brassica oleracea; Brassica rapa; Cajanus cajan; Camelina sativa; Camellia sinensis; Cannabis sativa; Capsella rubella; Capsicum annuum; Carex rostrata; Carica papaya; Carya illinoinensis; Castanea sativa; Chenopodium quinoa; Chlamydomonas reinhardtii; Chlorella variabilis; Cicer arietinum; Citrus sinensis; Citrus × clementina; Coccomyxa subellipsoidea; Coffea arabica; Coffea eugenioides; Cornus florida; Corylus avellana; Cryptomeria japonica; Cucumis melo; Cucumis sativus; Cucurbita maxima; C. moschata; C. pepo; Curcuma longa; Cynara cardunculus; Daucus carota; Dendrobium catenatum; Dioscorea cayenensis; Diospyros lotus; Durio zibethinus; Elaeis guineensis; Erigeron canadensis; Erythranthe guttata; Eucalyptus grandis; Euphorbia lathyris; Eutrema salsugineum; Fragaria vesca; Gastrolobium bilobum; Glycine max; G. soja; Gossypium arboreum/hirsutum/raimondii; Helianthus annuus; Henckelia pumila; Hevea brasiliensis; Hibiscus syriacus; Hordeum vulgare; Humulus lupulus; Impatiens glandulifera; Ipomoea nil/triloba; Jatropha curcas; Juglans microcarpa×regia/regia; Lactuca saligna/sativa/serriola/virosa（注：后三者按 NCBI 日期 2026-02-06/09 仍在 02-10 前，**实际是否在快照内未知**）; Lathyrus oleraceus; Lolium perenne/rigidum; Lotus japonicus; Lupinus angustifolius; Lycium barbarum/ferocissimum; Macadamia integrifolia; Magnolia sinica; Malania oleifera; Malus domestica/sylvestris; Mangifera indica; Manihot esculenta; Medicago truncatula; Mercurialis annua; Micromonas commoda/pusilla; Miscanthus floridulus; Momordica charantia; Monoraphidium neglectum; Morus notabilis; Musa acuminata; Neltuma alba; Nelumbo nucifera; Nicotiana attenuata/sylvestris/tabacum/tomentosiformis; Nymphaea colorata; Olea europaea; Oryza brachyantha/glaberrima/sativa; Ostreococcus lucimarinus/tauri; Panicum hallii/virgatum; Papaver somniferum; Phalaenopsis equestris; Phaseolus vulgaris; Phoenix dactylifera; Phragmites australis; Physcomitrium patens; Pistacia vera; Populus alba/euphratica/nigra/trichocarpa; Primulina eburnea/huaijiensis/tabacum; Prosopis cineraria; Prunus avium/dulcis/mume/persica; Punica granatum; Pyrus communis/×bretschneideri; Quercus lobata/robur/suber; Raphanus sativus; Rhodamnia argentea; Rhododendron vialii; Ricinus communis; Rosa chinensis/rugosa; Rutidosis leptorrhynchoides; Salvia hispanica/miltiorrhiza/splendens; Selaginella moellendorffii; Sesamum indicum; Setaria italica/viridis; Silene latifolia; Solanum dulcamara/lycopersicum/pennellii/stenotomum/tuberosum/verrucosum; Sorghum bicolor; Spinacia oleracea; Syzygium oleosum; Tarenaya hassleriana; Tasmannia lanceolata; Telopea speciosissima; Theobroma cacao/umbraticum; Trifolium pratense; Tripterygium wilfordii; Triticum aestivum/dicoccoides/urartu; Typha angustifolia/latifolia; Vicia villosa; Vigna angularis/radiata/umbellata/unguiculata; Vitis riparia/vinifera; Volvox carteri; Wolffia australiana; Zea mays; Zingiber officinale; Ziziphus jujuba。

> 注：其中 Lactuca saligna/serriola/virosa 与 Cucurbita pepo、Triticum urartu、Vigna umbellata 等的"首次 GCF"日期为 2026-02 或更早，若快照实际更早（本文推断 2025 年中），它们并不在 Carbon-A 训练内，故附录 B 对暴露的估计是**保守上界**。
