# P6 面板基线训练暴露核查（Tiberius / ANNEVO / GeneCAD / Vipsania / CADB）

读取日期：2026-11-07（本机日期）。性质：只读调研。没有下载基因组、权重或 GFF；没有运行任何计算；没有登录 ont/Baobab；没有修改已有文件；没有提交 git。
下载到会话临时目录的只有小文件：Tiberius 补充 PDF（0.5 MB）、ANNEVO 补充 XLSX（0.2 MB）、GeneCAD 补充 XLSX（2 个，0.03/0.23 MB）、若干仓库文本文件和 CADB 的 `metadata.json`（8 个，每个 5–135 KB）。
机器可读版：`exposure_audit_p6.json`（同一张表，`{物种: {工具: {code, status, evidence}}}`，另有 `_cadb_coverage`）。

码：S=同物种（Vipsania 里指同一装配的序列被预训练见过）；V=仅验证；T=仅测试；G=同属；Tr=同族或同亚科；F=同科；O=仅同目；–=无。
状态：(✓)已核实 = 读到原始来源并核对；(~)推断 = 由相关事实推出，无直接记载；(?)无法核实 = 来源未公开；(=)沿用 = 取自 `tool_exposure.json`，本次未重做（见 §2.6，未发现冲突）。

## 0. 结论摘要

1. Tiberius angiosperms：12 个格子全部用预印本 v3 的补充表 S3/S32 和 `angiosperms.yaml` 复核。与现有 `tool_exposure.json` 只有一处记法不同（Eschscholzia 应记 T，见 §3）。Acorus americanus 确在训练集（同属 G）。Figshare 里是一个 16.4 GB 的 tgz，没有下载，所以 Figshare 内部清单本身无法核实；物种清单以 S3 为准。
2. ANNEVO 现行 `ANNEVO_Magnoliopsida.pt`：训练物种清单**未公开**，12 个格子都只能记「无法核实」。它不是论文的 Embryophyta 模型（40 个 Train），论文版码不能代表它。
3. ANNEVO 参数：在 42c920f（v2.3.3）上**不需要也不能**加 `--num_classes 15 --boundary-aware`。这两个参数从 v2.3.0（2026-05-20）起就已删除，15 类边界感知是写死的默认；`--min_intron_length` 默认就是 20。所以 DEV 不带这些参数的命令，在 42c920f 上是对的（§2.2b）。
4. GeneCAD：现行默认模型（2026-11-05 刚替换）在 8 个物种上微调，论文版是 5 个；两者对面板物种的微调暴露只在 Eschscholzia、Akebia 上有差别（8 种版本有 Aquilegia）。骨干预训练语料模型卡里没写，按 PlantCAD2 的 65 属名单推断。
5. Vipsania：训练表含面板的 Akebia、Hakea、Juncus、Amborella，**还含敏感性物种 Dipteris**；accession 与面板相同或是面板 GCF 的配对 GCA。Juncus 同时是 Vipsania 的测试种。
6. CADB：核心 10 物种中 6 个在（Dryas、Akebia、Juncus、Rhynchospora、Acorus、Amborella），4 个不在（Eschscholzia、Hakea、Zizania、Themeda）；敏感性物种 Chenopodium 在，Dipteris 不在。GFF 无需登录即可下载。

## 1. 暴露码总表

| 物种（面板 accession） | Helixer a0080 | Tiberius angiosperms | ANNEVO 现行权重 | Carbon-A | Vipsania j9m0cdmk | GeneCAD 默认：整体 ｜ 微调(8 种) ｜ 骨干 |
|---|---|---|---|---|---|---|
| Dryas octopetala（GCF_963921425.1） | F (=) | F (✓) | 无法核实 (?) | F (=) | O (✓) | F (~) ｜ – (✓) ｜ F (~) |
| Eschscholzia californica（GCF_041429995.1） | – (=) | **T** (✓)；训练近缘 O | 无法核实 (?) | Tr (=) | Tr (✓) | O (~) ｜ O (✓) ｜ O (~) |
| Akebia trifoliata（GCF_017979445.1） | – (=) | O (✓)† | 无法核实 (?) | O (=) | **S** (✓) | O (~) ｜ O (✓) ｜ O (~) |
| Hakea prostrata（GCF_041222555.1） | – (=) | – (✓) | 无法核实 (?) | F (=) | **S** (✓) | – (~) ｜ – (✓) ｜ – (~) |
| Zizania latifolia（GCA_048537475.1） | Tr (=) | Tr (✓) | 无法核实 (?) | Tr (=) | F (✓) | Tr (~) ｜ Tr (✓) ｜ Tr (~) |
| Themeda triandra（GCA_057662145.1） | Tr (=) | Tr (✓) | 无法核实 (?) | Tr (=) | Tr (✓，同亚科) | Tr (~) ｜ F (✓) ｜ Tr (~) |
| Juncus effusus（GCA_027726005.1） | O (=) | O (✓) | 无法核实 (?) | O (=) | **S** (✓，也是其测试种) | O (~) ｜ O (✓) ｜ O (~) |
| Rhynchospora breviuscula（GCA_027562975.1） | O (=) | O (✓) | 无法核实 (?) | Tr (=) | Tr (✓，同亚科) | O (~) ｜ O (✓) ｜ O (~) |
| Acorus gramineus（GCA_030737835.1） | – (=) | **G** (✓) | 无法核实 (?) | – (=) | **G** (✓) | – (~) ｜ – (✓) ｜ – (~) |
| Amborella trichopoda（GCF_054954075.1） | V (=) | – (✓，范围外) | 无法核实 (?) | S (=) | **S** (✓) | – (~) ｜ – (✓) ｜ – (~) |
| Chenopodium album（敏感性，GCF_948465745.1） | V/G (=) | F (✓) | 无法核实 (?) | G (=) | Tr (✓，同亚科) | G (~) ｜ – (✓) ｜ G (~) |
| Dipteris shenzhenensis（敏感性，GCF_047967505.1） | – (=) | – (✓，范围外) | 无法核实 (?) | – (=) | **S** (✓) | – (~) ｜ – (✓) ｜ – (~) |

† Akebia 的 O（以及 Eschscholzia 的训练近缘 O）依赖 Aquilegia coerulea 确在 Tiberius 训练集：论文 S3、S32 写在，但 `angiosperms.yaml` 的 `training_species` 里没有（§2.1）。保守按含处理。
GeneCAD 整体码 = 微调码与骨干码取较高者。骨干语料未核实，所以整体一律记推断。

## 2. 逐项依据

### 2.1 Tiberius angiosperms（Mesangiospermae 模型）

来源（读取 2026-11-07）：
- Figshare 条目 `https://doi.org/10.6084/m9.figshare.32086728` → `https://figshare.com/articles/dataset/Tiberius_Mesangiospermae_Training_Data/32086728`；API `https://api.figshare.com/v2/articles/32086728`。发布 2026-04-23，许可字段 MIT，**只有 1 个文件** `tiberius_mesangiospermae_training_data.tgz`，16,394,213,247 字节，md5 `9aa8d554412190f492d02dd9184e13ed`。条目描述：「Genomes and annotations used to train, validate, and test the Mesangiospermae Tiberius model.」另一句：「use of these data is subject to the Fort Lauderdale Accord」。tgz 是单个 gzip 流，没有可读的文件目录，未下载，所以**内部清单无法核实**。
- 预印本 v3（2026-10-01）`https://doi.org/10.64898/2026.04.24.720536`，补充 PDF `https://www.biorxiv.org/content/biorxiv/early/2026/10/01/2026.04.24.720536/DC1/embed/media-1.pdf`（541,665 字节）。
- 配置 `https://github.com/Gaius-Augustus/Tiberius/blob/871d4979568ae09e56de7eca5d7489046ee833ea/model_cfg/angiosperms.yaml`（与 main `b3ea48f2`，2026-10-21 内容相同）。`target_species: Mesangiospermae`，`ncbi_tax_id: 1437183`，`date: 2026-04-15`。

**S3 的划分**：训练 45、验证 2、测试 7。
- 训练 45：Acorus americanus、Agave tequilana、Aquilegia coerulea、Arundo donax、Beta vulgaris、Brassica rapa、Brocchinia micrantha、Capparis spinosa、Chasmanthium laxum、Citrus limon、Cleome violacea、Coffea arabica、Corylus avellana、Daucus carota、Dendrobium catenatum、Dioscorea alata、Eucalyptus grandis、Fragaria vesca、Glycine max、Gossypium hirsutum、Hydrocotyle leucocephala、Joinvillea ascendens、Lactuca sativa、Lindenbergia philippensis、Malus domestica、Manihot esculenta、Mentha longifolia、Miscanthus sinensis、Moringa oleifera、Musa acuminata、Myrothamnus flabellifolia、Nicotiana tabacum、Oryza sativa、Osyris compressa、Panicum virgatum、Pharus latifolius、Phoenix dactylifera、Setaria viridis、Sorghum bicolor、Thinopyrum intermedium、Ulmus americana、Vitis vinifera、Yucca filamentosa、Zea mays、Zingiber officinale。
- 验证 2：Populus trichocarpa、Typha latifolia。
- 测试 7：Arabidopsis thaliana、Brachypodium distachyon、**Eschscholzia californica**、Freycinetia multiflora、Medicago truncatula、Mimulus guttatus、Urochloa brizantha。
- 正文：「No test species shares a species with any Tiberius training species.」S32 行：「E. californica order Ranunculales (order) 107.9 A. coerulea」，即 Eschscholzia 的最近训练种是 Aquilegia coerulea。

**YAML 与 S3 的差异**：`training_species` 44 项，逐项能对上 S3 的训练表，唯独**缺 Aquilegia coerulea**。YAML 的 `comment` 是「Angiosperm training weights trained by Tomas Bruna (JGI).」。仓库里同类文件出现过漏项（提交 871d497 的标题就是 "vertebrates.yaml: add 3 missing training species"），所以更可能是 YAML 漏写，但没有证据，实际权重是否含 Aquilegia 无法核实。影响：Eschscholzia 的训练近缘、Akebia 的 O 取决于它；若不含则为 –。

**逐物种**（面板物种自身都不在训练/验证里；Eschscholzia 是测试种）：
| 物种 | 码 | 依据 |
|---|---|---|
| Dryas | F | 训练含 Fragaria vesca、Malus domestica（Rosaceae） |
| Eschscholzia | T（训练近缘 O） | 测试种；最近训练种 Aquilegia coerulea（目级，107.9 Mya） |
| Akebia | O | 同目只有 Aquilegia coerulea |
| Hakea | – | 训练/验证/测试都没有 Proteales |
| Zizania | Tr | Oryza sativa（Kitaake）同族 Oryzeae |
| Themeda | Tr | Sorghum bicolor、Zea mays、Miscanthus sinensis 同族 Andropogoneae |
| Juncus、Rhynchospora | O | Poales 训练 12 种，无 Juncaceae、Cyperaceae |
| Acorus | G | Acorus americanus 训练（S3 与 YAML 都有） |
| Chenopodium | F | Beta vulgaris（Amaranthaceae） |
| Amborella、Dipteris | –（范围外） | NCBI 谱系核对：Amborella 不在 Mesangiospermae（1437183）内，蕨类更不在；S3 无此二者 |

另：Tiberius 预印本 S28 写明基准里 ANNEVO 用的是 `ANNEVO_Embryophyta.pt`，即 ANNEVO 的论文版，不是现行 Magnoliopsida 权重。

### 2.2 ANNEVO 现行 Magnoliopsida 权重

来源（读取 2026-11-07）：仓库 `https://github.com/xjtu-omics/ANNEVO`，固定 commit `42c920f70f52f96c7c06c665a639ec47d9304fff`（main 当前就是它，2026-08-11，v2.3.3）。

#### a. 训练了哪些物种：无法核实；证据表明不是论文的 40 个 Train 物种

1. 权重文件：`saved_model/ANNEVO_Magnoliopsida.pt`，git blob `21487d877430abf8e552c509560fd1218f19305e`，60,520,287 字节。同一 blob 在 v2.2.3（2026-04-23，路径当时是 `ANNEVO_model/`）、v2.3.0、v2.3.1、v2.3.2、v2.3.3 都在，没有变过。
2. 论文版的植物模型是 `ANNEVO_Embryophyta.pt`（v2.0/v2.1，blob `881caadd31…`），v2.2.3 之后不再随仓库发布。v2.3.0 发布说明：「For plants, we continue to use the latest model from v2.2.3」。Issue #36 中作者（2026-10-12）说植物模型是实验分支：「we did not retrain it in subsequent updates」。
3. 它和论文版不是同一类模型：分类头是 15 类（Issue #26 的报错：「copying a param with shape torch.Size([15, 64]) from checkpoint」），论文版是 5 类；v2.2.3 文档说明它用了「using a new strategy, which has not yet been published」的数据处理。
4. 作者的话（`docs/boundary_aware_model.md`，commit 745b3f9；该文件在 v2.3.0 的 3fdecd5 被删，所以 42c920f 上是 404）：「deliberately trained with all model organisms excluded, such as Arabidopsis thaliana and rice」。URL：`https://raw.githubusercontent.com/xjtu-omics/ANNEVO/745b3f9/docs/boundary_aware_model.md`。42c920f 上的 `docs/notes.md` 另有：「ANNEVO was intentionally trained without including model species」。论文 Table 24 把 Oryza sativa（GCF_034140825.1）列为 Train，与「排除水稻」矛盾，进一步说明两者训练集不同。
5. Issue #30（MarioStanke 问 2.3 的训练物种，作者 2026-07-13 答复）：「ANNEVO does not train on complete genomes from the selected species」，名单只承诺私下邮件发送（2026-07-15 对方邮箱已满，未见公开）。这条针对 2.3 系列整体，植物模型沿用 v2.2.3，没有单独的说明。
6. 其他公开渠道：HF 搜索 ANNEVO 无相关仓库；Zenodo 只有 `16749615`（论文时期对 566 个 RefSeq 装配的预测，不含权重说明）。README、docs、release notes 均无名单。
7. 论文补充 XLSX（`https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-026-03036-7/MediaObjects/41592_2026_3036_MOESM6_ESM.xlsx`）的 Table 24 是 Embryophyta 论文版：40 Train、9 Val&Test、119 Test；Table 33 给亚谱系（Embryophyta 含 Sapindales、Malvaceae、Brassicales、…、commelinids、Asparagales、Physcomitrium）。仍是论文版，不适用于现行权重。

结论：现行权重的 12 格全部记「无法核实」。**不能写成「ANNEVO 没见过 X」**。作者的邮箱（README 联系人）可以索取名单，本次没有联系。

**论文版 Table 24 的参考码**（12 格逐一复核，与 `tool_exposure.json` 的 `annevo_paper` 列一致，无误）：Dryas F（Prunus persica、Rosa chinensis）、Eschscholzia –（Papaver somniferum 仅 Test）、Akebia –、Hakea O（Nelumbo nucifera）、Zizania Tr、Themeda Tr（Panicum hallii、P. virgatum 同亚科）、Juncus O、Rhynchospora O、Acorus –、Chenopodium G（C. quinoa 为 Val&Test）、Dipteris –、Amborella T。

#### b. 参数问题：42c920f 上不需要、也不能加 `--num_classes 15 --boundary-aware`

- v2.2.3 时这个 15 类模型是「附加模型」：`annotation.py` 有 `--num_classes`（默认 5）和 `--boundary-aware`，`--min_intron_length` 默认 1。文档要求：「Simply add `--num_classes 15 --boundary-aware`.」另加 `--min_intron_length 20`。Issue #26：不加就会因分类头形状不符而报错。
- 在 42c920f：
  - `annotation.py` 的参数只有 `-g -m -l -o -t -s --tmp_path --batch_size --num_workers --overlap_pred --comp --min_intron_length(默认 20) --min_prot_length(默认 100) --show_log`，**没有** `--num_classes`、`--boundary-aware`；加上会被 argparse 当作未知参数报错。
  - `src/predict_nucleotide.py` 的 `base_pred(..., num_classes=15, ...)`，15 是写死的默认；`src/HMM.py` 的解码状态为 start/end/DSS/ASS/CDS/intron 的边界感知版（1+3+3+3+3+1+1=15 类），只有这一条解码路径。
  - README 的植物示例就是：`python annotation.py -g genome.fa -m saved_model/ANNEVO_Magnoliopsida.pt -l Magnoliopsida -o out.gff -t 48 --show_log`。`-l Magnoliopsida` 只决定窗口 `[30720, 5120, 32]`。README 说 `--overlap_pred` 对哺乳类、Actinopteri 强烈推荐，植物不需要（性能表里只有动物用了）。
- v2.3.0、v2.3.1、v2.3.2、v2.3.3 的 `annotation.py` 我逐个读过，都没有这两个参数，`--min_intron_length` 默认都是 20。
- 所以：照搬 DEV 不带这些参数的命令，**在 42c920f 上正确**。需要核对的只有一点：DEV 用的确实是 ≥v2.3.0 的代码。可以看 DEV 输出 GFF 的第一行，解码器从 v2.3.0 起会写「# This output was generated with ANNEVO (v2.3.x).」，x 即所用版本。若 DEV 其实是 v2.2.3 以前的代码，不带参数不会得到结果，而是报形状错误。
- 「Magnoliopsida.pt 本身是 15 类边界感知模型」：已核实（间接）。证据是 Issue #26 的报错（检查点的 `task_layer.fc.weight` 为 [15, 64]）、v2.2.3 文档和 42c920f 的 15 类默认；没有直接打开权重。

### 2.3 GeneCAD

来源（读取 2026-11-07）：
- 预印本 2026-05-12 版（bioRxiv API 把 v3、v4 都标在这一天，取到的 JATS 内部标签是 v3）JATS `https://www.biorxiv.org/content/early/2026/05/12/2025.10.31.685877.source.xml`；补充表 1 `https://www.biorxiv.org/content/biorxiv/early/2026/05/12/2025.10.31.685877/DC1/embed/media-1.xlsx`。
- 仓库 `https://github.com/plantcad/genecad`（main `f641abde`，2026-10-29）的 `train.sh`、`src/config.py`、README。
- HF：`https://huggingface.co/plantcad/genecad_plant`（main `21c98c7f895a7f3c5b65f9e7c4c76b7ac7229963`，2026-11-05）、`https://huggingface.co/plantcad/genecad_revision_models`（`training_inputs/README.md`、`eight_species.yaml`）、骨干 `https://huggingface.co/emarro/pcad2-200M-cnet-baseline`。

**现行默认模型**（`genecad predict -m plant`）：模型卡写明它是修订版 M seed 11，「Training: eight angiosperm genomes, 1,753 steps, effective batch 384」。8 个微调物种见 `eight_species.yaml`：Arabidopsis thaliana、Oryza sativa、Glycine max、Hordeum vulgare、Populus trichocarpa、Aquilegia coerulea、Nymphaea colorata、Triticum aestivum。验证窗口：「A. thaliana, O. sativa and T. aestivum provide the validation windows.」
**旧版**：`genecad_plant` 在 2026-07-03 提交（标签 `released-2026-07-03`）里是旧模型，模型卡为空，训练集未记载。论文（2026-05-12 版）和 `train.sh` 的默认是 5 个物种（`SPECIES_IDS="Athaliana Osativa Gmax Hvulgare Ptrichocarpa"`）。以后引用时要写清 revision，不要只写 `main`。
**留出/测试**（论文）：Juglans regia、Coffea arabica、Phaseolus vulgaris、Nicotiana sylvestris、Nicotiana tabacum；另有 Zea mays（「Maize was excluded from the GeneCAD and ANNEVO training sets」）；动物测试 Gallus gallus。**面板物种均不在**。
**骨干**：`emarro/pcad2-200M-cnet-baseline`（HNetForCausalLM，258M 参数，revision `5574ff5d`）的模型卡是空模板，预训练语料无法核实。论文称 PlantCAD2 在 65 个植物基因组上预训练；补充表 1 列了 64 个基因组的调查表。这里按报告 10 附录 A 的 PlantCAD2 65 属名单推断骨干暴露，结果与 `tool_exposure.json` 的 `plantcad2` 列一致。注意骨干名单里有 Zea mays，而微调里排除了它。

微调码（8 种）与骨干码（推断）：Dryas –/F（Fragaria、Malus、Prunus）；Eschscholzia O/O（Aquilegia）；Akebia O/O；Hakea –/–；Zizania Tr/Tr（Oryza）；Themeda F/Tr（微调只有 Poaceae 的 Oryza、Hordeum、Triticum；骨干有 Sorghum、Zea、Miscanthus）；Juncus O/O；Rhynchospora O/O；Acorus –/–；Chenopodium –/G（C. quinoa）；Dipteris –/–；Amborella –/–（Nymphaea colorata 属 Nymphaeales，不同目）。
对旧版（5 种）：Eschscholzia、Akebia 的微调码变为 –，其余不变；整体码取骨干推断值，不变。

### 2.4 Vipsania Streptophyta（j9m0cdmk）

来源（读取 2026-11-07）：仓库 `https://github.com/Gaius-Augustus/Vipsania` 的 `docs/training_species.tsv`、`docs/test_species.tsv`（HEAD `397660e`，2026-10-21）。会话笔记目录里的两个副本与 HEAD 逐字节相同。README：模型表 Streptophyta `j9m0cdmk`，训练 200 个、测试 13 个物种；「Since Vipsania never learns from annotations」，所以训练集里有某种不是标签泄漏。`--finetune` 会再在目标基因组序列上无监督继续训练。

**accession 核对**（面板 accession 用 NCBI Datasets 取得配对 GCA）：
| 物种 | 面板 accession | Vipsania 训练表 | 关系 |
|---|---|---|---|
| Akebia trifoliata | GCF_017979445.1 | GCA_017979445.1 | 面板 GCF 的配对 GCA，同一装配 |
| Hakea prostrata | GCF_041222555.1 | GCA_041222555.**2** | 面板 GCF 的配对 GCA 就是 .2 版，同一装配 |
| Juncus effusus | GCA_027726005.1 | GCA_027726005.1 | 相同；同时在其 Streptophyta 测试表 |
| Amborella trichopoda | GCF_054954075.1 | GCA_054954075.1 | 配对 GCA，JGI v2.0 装配 |
| Dipteris shenzhenensis（敏感性） | GCF_047967505.1 | GCA_047967505.1 | 配对 GCA，同一装配（报告 28 未列） |
| Acorus gramineus | GCA_030737835.1 | 无；有 Acorus calamus GCA_030737845.1 | 同属，相邻 accession 但不是同一装配 |

其余物种的近缘取自 NCBI Taxonomy 谱系（逐个核对到属/族/亚科/科/目）：
- Dryas：Rosales 内只有 Ficus carica（Moraceae），无 Rosaceae → O。
- Eschscholzia：Papaver somniferum GCA_010119995.1（同亚科 Papaveroideae）、Corydalis saxicola → Tr。
- Zizania：Poaceae 有 Holcus mollis、Panicum virgatum、Pharus latifolius，没有 Oryzoideae → F。
- Themeda：Panicum virgatum GCF_016808335.1 同亚科 Panicoideae → Tr（亚科）；Andropogoneae 无训练种，Sorghum、Zea 只是 Vipsania 的测试种。
- Rhynchospora：Carex hirta GCA_965112215.1（Cyperaceae、Cyperoideae）→ Tr（亚科）。
- Chenopodium：Dysphania ambrosioides GCA_026401315.1（同亚科 Chenopodioideae）→ Tr；Amaranthus hybridus 同科。
- Rhynchospora、Juncus 同在 Poales，Typha、Tillandsia、Ananas、Carex、Holcus 等都是同目。

这次也确认了报告 28 §1 的说法：核心 10 物种里 4 个在训练集（Juncus、Akebia、Hakea、Amborella），且补充：Dipteris 也在，Juncus 还是测试种。Vipsania 其他 16 个类群模型（含 `cg6grhms` 的「以上都不是」）不在本次范围内。

### 2.5 Carbon-A 的 CADB

来源（读取 2026-11-07）：
- 博客 `https://huggingface.co/blog/HuggingFaceBio/carbon-annotator-genbank-genome-annotation`：「The current release contains 48,167 assemblies from 22,617 taxa」；「We have now annotated roughly half of the genomes in our target GenBank set」，计划约三周后再发一批。所以 CADB 是不完整、仍在增长的库。
- 集合 `https://huggingface.co/collections/HuggingFaceBio/carbon-annotation-database`（2026-10-28 更新）；数据在公开存储桶 `https://huggingface.co/buckets/HuggingFaceBio/genbank-annotations`（71 TB、766,839 个文件，匿名可读）；浏览 Space `https://huggingface.co/spaces/HuggingFaceBio/genbank-annotation-explorer`。
- 桶的布局：`annotations/<division>/…parquet` 是逐碱基 CDS 概率（不是 GFF）；**基因模型在 `cds/<division>/<GCA 号>/`**，每个装配有 `annotations.gff3.gz`、`cds.fna.gz`、`proteins.faa.gz`、`metadata.json`。`cds/plant` 目前 5,838 个装配。
- 下载方式（匿名，HEAD 实测 200，没有下载正文）：`https://huggingface.co/buckets/HuggingFaceBio/genbank-annotations/resolve/cds/plant/<GCA 号>/annotations.gff3.gz`；也可用 `hf` CLI。
- **只按 GenBank 的 GCA 号索引**。面板里的 GCF 要换成配对的 GCA 来查。

覆盖（用 `cds/<6 个 division>/<GCA>` 逐个探测；探测结果与 Space 的 `data/coverage.json`（2026-10-05 的清单）一致）：
| 物种 | 面板 accession | CADB 目录 | 在否 | 基因数 | 文件上传日 | GCF/GCA 配对比较 |
|---|---|---|---|---|---|---|
| Dryas octopetala | GCF_963921425.1 | GCA_963921425.1 | 在 | 21,654 | 2026-10-06 | different |
| Eschscholzia californica | GCF_041429995.1 | GCA_041429995.1 | **不在** | – | – | – |
| Akebia trifoliata | GCF_017979445.1 | GCA_017979445.1 | 在 | 23,445 | 2026-10-06 | different |
| Hakea prostrata | GCF_041222555.1 | GCA_041222555.2 / .1 | **不在** | – | – | – |
| Zizania latifolia | GCA_048537475.1 | 同 | **不在** | – | – | – |
| Themeda triandra | GCA_057662145.1 | 同 | **不在** | – | – | – |
| Juncus effusus | GCA_027726005.1 | 同 | 在 | 20,805 | 2026-10-07 | na（无 RefSeq） |
| Rhynchospora breviuscula | GCA_027562975.1 | 同 | 在 | 23,265 | 2026-10-12 | na |
| Acorus gramineus | GCA_030737835.1 | 同 | 在 | 22,229 | 2026-10-07 | na |
| Amborella trichopoda | GCF_054954075.1 | GCA_054954075.1（v2.0）；另有 GCA_000471905.1（AMTR1.0） | 在 | 16,086（AMTR1.0 为 17,122） | 2026-10-13 | v2.0：na（CADB 抓取时 GCF 尚未发布） |
| Chenopodium album（敏感性） | GCF_948465745.1 | GCA_948465745.1 | 在 | 78,367 | 2026-10-06 | different |
| Dipteris shenzhenensis（敏感性） | GCF_047967505.1 | GCA_047967505.1 | **不在** | – | – | – |

核心 10 物种：6 个在，4 个不在。这 4 个仍要自己跑 Carbon-A（或等下一批 CADB，到 12-01 前再探测一次）。

**版本与日期**（`metadata.json`）：概率模型写的是 `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`，revision「Not recorded」；Carbon-A-1.2B 仓库的 `provenance.json` 说明其 `model.safetensors` 与该预览版（b0483c23）字节相同，所以可以视为同一权重，但 CADB 没有记录具体 revision（推断）。后处理是 `HuggingFaceBio/GENERanno-cds-annotation-pipeline` rev `6ffa159`，`minimum_cds_length_bp` 96、`maximum_boundary_intron_length_bp` 20000，带 `transition_network.pt`、`confidence_network.pt`。文件上传日 2026-10-06 至 10-13。Carbon-A-1.2B 仓库的流程在此后又改过（2026-10-09 到 11-04），CADB 不一定对应最新流程。
**GFF 是否带置信分**：推断为带。依据是 metadata 里有置信网络权重、博客写「a confidence score」、本地 Carbon-A 流程的说明是「Every gene row of prediction.gff … carries confidence=<score>」。没有读 GFF，属性名未直接核实。
**与面板 FASTA 的坐标对应**（只用 NCBI 序列报告和 metadata，不读 GFF）：CADB 注释的是 GenBank 序列，面板的 TH 物种用 RefSeq 序列。Dryas、Akebia、Chenopodium 的配对比较是 different。逐 contig 核对：公共 contig 长度完全一致（Dryas 24/24、Akebia 890/890、Chenopodium 156/156，0 个长度不符）；差别只是少量小 contig 在一边有、另一边没有（Akebia 的 RefSeq 多 24 条 <7 kb 的未定位序列；Dryas 的 RefSeq 少 3 条 GenBank contig、多 1 条 2 kb 的）。可用，但需要用序列报告把 GenBank 号映射到 RefSeq 号，序列碱基级一致性只比了长度，没有比序列。
**许可**：HF 桶与集合页没有找到 CADB 的许可条款，未核实。

### 2.6 Helixer a0080 与 Carbon-A：沿用，未见冲突

本次读到的事实没有与 `tool_exposure.json` 的这两列冲突。（Tiberius S17 与 Helixer 的训练集信息一致；Carbon-A 的训练快照仍是 2,055 个 RefSeq GCF，面板的 TH 物种都是之后才有 RefSeq 的。）CADB 与 Vipsania 的发现不影响这两列。

## 3. 与现有 `tool_exposure.json` 不同的格子

1. Tiberius × Eschscholzia：现为「O (test species)」，建议记 **T**（本种是测试种），并注明训练近缘 O（Aquilegia）。实质不变。其余 11 格全部复核一致。
2. ANNEVO：`annevo_paper` 列对论文版复核无误，但它不是现行权重；现行权重 12 格全部应记**无法核实**。凡是引用「ANNEVO 对 X 为盲测」的地方都要改。
3. 新增两列（原文件没有）：Vipsania（§2.4）、GeneCAD（§2.3）。
4. 报告 28 §1「训练集含 4 个核心种」属实；补充：Dipteris（敏感性）也在训练集，Juncus 同时是其测试种。
5. 报告 10 §1.5「现行权重训练清单未公开」属实；§1.6「GeneCAD 训练物种清单未核实」现已核实（微调 8 种，旧版 5 种）。

## 4. 仍无法核实

- ANNEVO 现行权重的训练物种（只有作者能提供）。
- Tiberius：Figshare tgz 内部清单（16.4 GB，未下载）；`angiosperms.yaml` 缺 Aquilegia coerulea 是漏写还是实际未训练。
- GeneCAD：骨干 `pcad2-200M-cnet-baseline` 的预训练语料；旧版（2026-07-03）的微调集；修订稿（Bioinformatics BIOINF-2026-2172）的正文未见公开。
- CADB：GFF 的置信分属性名（未读 GFF）、概率模型的具体 revision、许可条款、缺失 4+1 个物种何时补齐。
- Vipsania 预印本全文（bioRxiv 10.64898/2026.08.26.747235）本次未读，只用了仓库 docs 和 README。

## 5. 对后续工作的影响

- ANNEVO 用现行权重跑面板时，暴露码写「未披露」，不要写「盲测」；若要做分层曲线，只能把它放在不带暴露假设的总体表里。
- 面板里 Vipsania 的 S 物种（Akebia、Hakea、Juncus、Amborella、Dipteris）是序列级见过，要在分层结果里单独标记；它无监督，不构成标签泄漏，但作者推荐的 `--finetune` 本来就会在目标基因组上训练。
- GeneCAD 在 12-01 后跑时固定 revision（建议记 `genecad_plant@21c98c7f`，或明确选旧版标签）。
- Carbon-A：核心 10 物种里 6 个可从 CADB 直接取 GFF（需映射 contig 名，并按 §2.5 说明版本）；Eschscholzia、Hakea、Zizania、Themeda 要自己跑。
