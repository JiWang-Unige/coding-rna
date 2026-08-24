# Deep Research Report — claude · 2026-06-09

> 占位文件，由 `$research-interview` 在 2026-06-09 自动生成。
>
> **使用方式**: 把外部平台（claude）跑出来的完整 deep research 报告整段粘贴到下面
> "## Report body" 之下。可以保留这段说明也可以删除，由 `$research-synthesize` 自动跳过。
>
> 对应的 prompt: `docs/inputs/deep_research_prompt_2026-06-09.md`

---

## Source platform
- Platform: claude
- Mode: Deep Research
- Run date:
- Number of retries / refinements:
- Subjective quality (filled later): low / medium / high

## Report body

# 跨物种 ab initio 深度学习用于真核生物蛋白编码基因注释：研究现状报告

## TL;DR

- **当前最强的 ab initio 真核基因注释模型是端到端深度学习方法 Tiberius（CNN+LSTM+可微分 HMM）、2025 年版 Helixer（CNN/LSTM+HMM 后处理），以及 2026 年发表于 Nature Methods 的 ANNEVO（MoE 基因组语言模型+Viterbi 解码）**；在哺乳动物上 Tiberius 平均达到 89.7% 外显子级 F1、55.1% 基因级 F1，人类基因级 F1 达 62%（次优 ab initio 方法仅 21%），已可与使用 RNA-seq 与蛋白证据的 BRAKER3 媲美——这是 Stanke 等所述"25 年来由 HMM 主导"的基因预测领域最大的一次范式转变。
- **以 SegmentNT、DNABERT-2、GENERanno 为代表的"基础模型+微调"路线在单核苷酸级元素分割上表现强劲，但大多止步于碱基级标注、缺乏结构化解码，难以直接产生符合生物学约束（剪接位点、阅读框、ORF）的完整合法基因模型**；GENERanno 真核注释器目前仍是无技术文档、无公开真核指标的 preview 版本。
- **面向 Nature Methods 级别贡献的最大空白在于：(a) 长上下文骨干（Mamba/S4/扩张卷积达 100kb+）、(b) 可微分结构化解码（semi-CRF/可微分 HMM）以保证基因级一致性、(c) 真正的跨物种表示学习与生物学约束感知解码**；GeneDecoder、GeneCAD 等近期工作已开始用 CRF/latent-CRF 把约束嵌入解码层（GeneCAD 在 5 个留出被子植物上 transcript 级 F1 较 Helixer/BRAKER3 平均提升 8–10%），但尚未在大规模跨物种基准上系统验证。

---

## Key Findings

1. **范式已从"浅层 HMM"转向"深度学习+结构化解码"。** 经典 GHMM 基因发现器（AUGUSTUS、SNAP、GeneMark 系列）统治了约 25 年（Tiberius 摘要："For more than 25 years, learning-based eukaryotic gene predictors were driven by hidden Markov models"）；2020 年 Helixer 首次证明深度学习层+独立 HMM 后处理可超越 AUGUSTUS，2024 年 Tiberius 进一步把 HMM 做成端到端可微分层，人类基因级 F1 从次优 ab initio 方法的 21% 跃升至 62%。

2. **跨物种泛化仍是核心痛点。** Tiberius 最初仅训练于哺乳动物，作者明确不建议未经重训练用于非脊椎动物；2026 年其多 clade 扩展版为 Mesangiospermae、Fungi、Vertebrata、Insecta、Chlorophyta、Bacillariophyta 分别训练 lineage-specific 模型。Helixer 同样分 clade 训练脊椎动物/无脊椎动物/陆生植物/真菌模型。没有任何模型实现真正"一个模型通吃所有真核生物"。

3. **基础模型路线"标注准、成模型难"。** SegmentNT、GENERanno、DNABERT-2 等做的是碱基/token 级多类分割，输出的是概率而非合法转录本；正如 GeneCAD 作者指出"nucleotide-level labeling alone does not constitute gene prediction; without structure-aware decoding, outputs can violate annotation rules and cannot be assembled into valid transcript models"。

4. **生物学约束的整合程度差异巨大。** Tiberius 与 Helixer 通过 HMM 状态结构显式编码相位、剪接基序（GT-AG/GC-AG/AT-AC）；ANNEVO 用 Viterbi 解码强制内含子状态群对应三种剪接模式；而纯分割模型（SegmentNT）几乎不施加这些约束，依赖后处理。

5. **可复现性整体良好但有梯度。** Tiberius、Helixer、ANNEVO、SegmentNT、DNABERT-2 均开源且提供预训练权重；GENERanno 真核注释器权重可下载但无论文指标、无真核评测子集，复现风险最高。

---

## Details

### 第 1 节 方法族（Method Families）

#### (a) 经典 HMM/GHMM 基因发现器
- **核心思想**：用（广义）隐马尔可夫模型把基因结构（外显子/内含子/UTR/基因间区）建模为状态序列，状态发射概率刻画密码子、剪接位点、长度分布等信号，Viterbi 解码求最可能的状态路径。
- **代表工作**：
  - AUGUSTUS：Stanke & Waack, "AUGUSTUS: a web server for gene finding in eukaryotes"（2003，Bioinformatics）；Stanke et al., "AUGUSTUS: ab initio prediction of alternative transcripts"（2006，NAR，含可变剪接预测，DOI 10.1093/nar/gkl200），基于 GHMM。
  - SNAP：Korf, I., "Gene finding in novel genomes"（2004，BMC Bioinformatics）。
  - GeneMark 系列：GeneMark-ES/ET/EP+，自训练 GHMM。
  - GenScan：Burge & Karlin（1997，J Mol Biol 268:78–94）。
- **优点**：无需 GPU、可解释、对紧凑基因组（如酵母）精度尚可、产生符合语法的完整基因模型。
- **缺点**：参数需逐物种训练/调参；精度随基因组增大、基因密度降低而显著下降；对长内含子、复杂结构建模能力有限。
- **典型性能**：在 G3PO 多物种基准（Scalzitti et al. 2020，BMC Genomics）上，AUGUSTUS 核苷酸级 F1 约 0.52（最高），SNAP 约 0.39；外显子级灵敏度仅约 0.27。在 Tiberius 论文的哺乳动物测试中 AUGUSTUS 平均外显子 F1 67.3%、基因 F1 仅 12.4%。

#### (b) 混合 HMM+证据/神经方法
- **核心思想**：在 GHMM 之上整合外部证据（RNA-seq 比对、蛋白同源比对）或把 GHMM 嵌入更大流水线。
- **代表工作**：
  - GeneMark-ETP：Brůna, Lomsadze, Borodovsky, "GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes"（2024，Genome Research 34(5):757–768，DOI 10.1101/gr.278373.123）。整合基因组、转录组、蛋白证据，先确定"高置信"基因再迭代训练 GHMM。
  - BRAKER3：Gabriel et al.（2024，Genome Research），结合 GeneMark-ETP、AUGUSTUS、TSEBRA，使用 RNA-seq+蛋白库。
  - AUGUSTUS-CGP（comparative gene prediction）：基于多基因组比对的比较基因预测。
- **优点**：在有证据物种上精度最高，是社区标准。
- **缺点**：依赖 RNA-seq/蛋白证据的可得性与质量；流水线复杂、运行慢。
- **典型性能**：BRAKER3 较 BRAKER1/2 转录本级 F1 平均提升约 20 个百分点；GeneMark-ETP 在大型 GC 不均一基因组（鸡 G. gallus、小鼠 M. musculus）上较旧工具有两位数 Sn/Pr 提升。

#### (c) CNN/LSTM/U-Net 序列到序列模型
- **核心思想**：直接以 one-hot DNA 为输入，用卷积/循环网络做逐碱基多类分类，再以 HMM 后处理或可微分 HMM 解码出基因模型。
- **代表工作**：
  - Helixer（初版）：Stiehler, Steinborn, Scholz, Dey, Weber, Denton, "Helixer: cross-species gene annotation of large eukaryotic genomes using deep learning"（2020，Bioinformatics 36(22-23):5291–5298，DOI 10.1093/bioinformatics/btaa1044）。
  - Helixer（2025 Nature Methods 版）：Holst, Bolger, Kindel, Günther, Maß, Triesch, Kiel, Saadat, Ebenhöh, Usadel, Schwacke, Weber, Bolger, Denton, "Helixer: ab initio prediction of primary eukaryotic gene models combining deep learning and a hidden Markov model"（Nature Methods，2025，DOI 10.1038/s41592-025-02939-1）。
  - Tiberius：Gabriel, Becker, Hoff, Stanke, "Tiberius: end-to-end deep learning with an HMM for gene prediction"（2024，Bioinformatics 40(12):btae685，DOI 10.1093/bioinformatics/btae685；bioRxiv 2024.07.21.604459）；端到端整合 CNN+双向 LSTM+可微分 HMM 层。
- **优点**：GPU 上极快（据 Tiberius 作者测试，注释人类基因组约 1 小时 39 分，Helixer 约 8 小时 54 分）；无需逐物种证据；端到端可微分 HMM 保证基因级一致性。
- **缺点**：固定窗口限制长上下文；跨 clade 需重训练；训练依赖高质量参考注释。
- **典型性能**：见第 2 节。

#### (d) 基础模型+微调
- **核心思想**：先在大规模无标注 DNA 上自监督预训练得到通用表示，再加分割头/token 分类头微调到注释任务。
- **代表工作**：
  - Nucleotide Transformer / SegmentNT：de Almeida, Dalla-Torre, Richard et al.（2024，bioRxiv 2024.03.14.584712，DOI 10.1101/2024.03.14.584712）。
  - DNABERT-2：Zhou et al., "DNABERT-2: Efficient Foundation Model and Benchmark For Multi-Species Genome"（2023，arXiv 2306.15006；ICLR 2024），BPE tokenization、ALiBi、Flash Attention；权重 zhihan1996/DNABERT-2-117M。
  - GENERanno：Li, Wu, Zhu, Feng, Ye, Wang, "GENERanno: A Genomic Foundation Model for Metagenomic Annotation"（2025，bioRxiv，DOI 10.1101/2025.06.04.656517）。
  - GeneCAD（PlantCAD2 基础模型+ModernBERT+CRF）：bioRxiv 2025，DOI 10.1101/2025.10.31.685877。
- **优点**：可迁移、少标注下泛化、捕获进化保守信号。
- **缺点**：输出多为碱基级概率、缺结构化解码；上下文窗口受限（NT/SegmentNT 约 30kb、GENERanno 8kb、DNABERT-2 更短）；推理成本高。
- **典型性能**：见第 2 节。

#### (e) 新兴 transformer/attention/MoE ab initio 方法
- **代表工作**：ANNEVO（Ye et al., "Highly accurate ab initio gene annotation with ANNEVO"，Nature Methods 2026，DOI 10.1038/s41592-026-03036-7；Research Square 预印 DOI 10.21203/rs.3.rs-6402260/v1）——MoE 基因组语言模型建模长程依赖与跨物种联合进化关系，配合 Viterbi 结构化解码。配套 News & Views "Learning genes deeply"（DOI 10.1038/s41592-026-03035-8）。
- **优点**：建模远程依赖与进化关系；566 物种大规模基准；无需外部证据即可超越 ab initio 基线、媲美证据流水线。
- **缺点**：模型大、推理成本高；非 OSI 开源许可（ANNEVO Non-Commercial License，仅限学术非商业）。

#### (f) 事后基因级精修 / 结构化解码
- **核心思想**：在神经网络输出的碱基级 logits 之上，用 CRF/semi-CRF/HMM/latent-CRF 强制状态转移合法性（相位一致、特征顺序、剪接规则）。
- **代表工作**：
  - GeneDecoder：Marin, Pultz, Boomsma, "Gene finding revisited: improved robustness through structured decoding from learned embeddings"（2025，arXiv 2505.03377）——扩张卷积+残差连接+BiLSTM 特征模型 + latent CRF，可手动设定转移势以编码合法状态转移。
  - GeneCAD：用染色体级 CRF 强制 splice-phase 与 feature order，BILOU 标签捕获 UTR/内含子/编码段边界，再用蛋白语言模型筛除 repeat 驱动的假 ORF。据论文，在 5 个留出被子植物（含异源四倍体烟草 Nicotiana tabacum）上 transcript 级 F1 较 Helixer 和 BRAKER3 平均提升 8–10%，并锐化起止密码子与剪接连接边界。
  - HelixerPost / Tiberius 可微分 HMM：把结构化解码内嵌网络。
- **优点**：保证输出合法基因模型、提升边界精度。
- **缺点**：仍依赖底层表示质量；大规模跨物种验证不足。

### 第 2 节 当前 SOTA 候选

#### Tiberius
- 论文：Bioinformatics 2024，DOI 10.1093/bioinformatics/btae685（bioRxiv 2024.07.21.604459）；多 clade 扩展预印本 "Accurate ab initio gene prediction in eukaryotes with Tiberius in multiple clades"（bioRxiv 2026，benchmark 覆盖 33 物种，对比 Helixer 与 ANNEVO）。
- GitHub：https://github.com/Gaius-Augustus/Tiberius 。预训练权重随仓库提供，并有 Web 界面与 Galaxy ToolShed。
- 架构：CNN + 双向 LSTM + 可微分 HMM 层，自定义 gene prediction loss；训练模式输出后验状态概率并最小化误分类损失。
- 指标（哺乳动物测试集 H. sapiens、B. taurus、D. leucas 平均，Supplementary Table S6）：外显子级 F1 **89.7%**、基因级 F1 **55.1%**（Helixer 对应 72.9%/19.3%，AUGUSTUS 67.3%/12.4%）；人类基因级 F1 **62%**（次优 ab initio 方法 21%）；BUSCO 完整度平均 96.0%（Helixer 92.1%、AUGUSTUS 74.2%）。de novo 模式下约 2/3 人类基因外显子-内含子结构完全正确；正确预测的最长人类基因 **328,931 bp**。
- 数据/划分：37 个哺乳动物基因组，RefSeq 注释，RepeatModeler2/RepeatMasker/TRF 软屏蔽；物种级留出测试（held-out species）。
- 跨物种 vs 同物种：跨物种（在未见物种上评测），但限于训练 clade（哺乳动物）。

#### Helixer（2025 Nature Methods 版）
- 论文：Nature Methods 2025，DOI 10.1038/s41592-025-02939-1（早期 Bioinformatics 2020 版 DOI 10.1093/bioinformatics/btaa1044）。
- GitHub：https://github.com/weberlab-hhu/Helixer 。权重见 Zenodo（DOI 10.5281/zenodo.17414354；v2 记录 https://zenodo.org/records/17850139 ，2025-12-07）。
- 架构：CNN+双向 LSTM 逐碱基预测 + HMM（HelixerPost）后处理输出 gff3；HMM 把 CDS 拆为按相位与密码子类型（start/stop/regular）的 10 个子状态、内含子拆为 60 个子状态（按起始/延续、剪接基序、外侧状态）。
- 指标：跨真菌、植物、脊椎动物、无脊椎动物；论文 Table 2 报告 exon/intron/intron chain/transcript 级 F1（均值），与专家注释多指标接近。在与 Tiberius 直接对比的哺乳动物设置中（H. sapiens、B. taurus、D. leucas），Helixer 落后 Tiberius（基因级 F1 约 19.3% vs Tiberius 55.1%）。在植物测试集上 Helixer 召回高于参考但精度低于参考；装配最碎的物种（Papaver somniferum、Triticum dicoccoides）表现最差。
- 跨物种：明确以跨物种为设计目标，分 clade 预训练模型。

#### ANNEVO
- 论文：Nature Methods 2026，DOI 10.1038/s41592-026-03036-7（Research Square 预印 DOI 10.21203/rs.3.rs-6402260/v1）。
- GitHub：https://github.com/xjtu-omics/ANNEVO （ANNEVO Non-Commercial License，学术免费、商业需授权；联系人 Pengyu Zhang、Kai Ye，西安交大）。
- 架构：MoE 基因组语言模型建模远程依赖+联合进化关系；Viterbi 解码，内含子状态群对应 GC-AG/GT-AG/AT-AC 三种剪接模式（如从 CDS0 进入内含子则从 CDS1 退出，强制相位一致）。
- 指标：在 566 个系统发育多样物种上基准；12 个模式物种（跨 Fungi、Embryophyta、Invertebrates、Mammalia、Vertebrate_other）上，**较 BRAKER3 基因级 F1 绝对提升 4%、较深度学习基线 Helixer 提升 11%**；某编码区识别任务 mean F1 0.92、recall 0.922（最高，体现更完整的编码区识别）。较 Augustus 在核苷酸级 F1 提升 7.5–22.3%、基因级 F1 提升 9.9–38.5%、BUSCO 提升 9–34.5%。在 793 物种（566 RefSeq + 227 Ensembl）中，有 **252 个（31.8%）BUSCO 超过参考注释**，其中真菌、被子植物中分别有约 70%、71% 物种超过 Ensembl 注释。在 ANNEVO 自身基准下，Tiberius 在 Vertebrate_other clade 上性能骤降。
- 跨物种：以进化关系建模实现跨物种，覆盖五大 clade。

#### SegmentNT
- 论文：de Almeida, Dalla-Torre, Richard et al.，bioRxiv 2024，DOI 10.1101/2024.03.14.584712（"Annotating the genome at single-nucleotide resolution with DNA foundation models"）。
- 权重/代码：https://huggingface.co/InstaDeepAI/segment_nt ；https://github.com/instadeepai/nucleotide-transformer 。许可 CC-BY-NC-SA-4.0。
- 架构：NT v2 500M multi-species encoder + 1D U-Net 分割头（2 下采样+2 上采样卷积块，分割头 53M 参数，总计 562M），预测 14 类基因组元素（蛋白编码基因、lncRNA、5'UTR、3'UTR、外显子、内含子、剪接供体/受体位点、polyA 信号、组织不变/特异启动子与增强子、CTCF 结合位点）。各元素独立预测概率（允许重叠）。输入最长 30kb（训练，5001 tokens），可零样本泛化到 50kb（用 YaRN/Rotary 重标定）。在 23B tokens 上训练 3 天（8×H100）。
- 数据/划分：人类全染色体，chr20/chr21 留作测试集、chr22 验证集——染色体级划分；启动子/增强子来自 ENCODE cCRE。
- 指标：碱基级与剪接位点检测强。据独立基准（bioRxiv 2026.02.22.707219），在精简训练数据类上剪接位点 F1 约 0.90，对外显子/基因综合注释 F1bw 约 0.624–0.660、F1seg 约 0.664–0.724，随注释复杂度（如 TE 衍生外显子、可变剪接外显子）下降。注意 SegmentNT 存在与序列内位置相关的系统性偏差（bioRxiv 2025.04.09.647946）。
- 跨物种：主体为人类同物种（染色体级留出）；论文另把最佳 30kb 模型微调到多物种，泛化到未见动植物。

#### GENERanno-eukaryote-1.2b-cds-annotator-preview
- HuggingFace：https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview （2026-02-10 发布）。论文 DOI 10.1101/2025.06.04.656517（v3，2025-08-14）——但论文内容几乎全是**原核 0.5B 模型**，真核注释在论文中被列为 future work/局限。
- GitHub：https://github.com/GenerTeam/GENERanno 。
- 架构：GENERanno 家族为 Transformer encoder（类 Llama-for-MLM/ModernBERT，全层全局注意力），**单核苷酸 tokenization**（非 k-mer/BPE），上下文 **8192 bp**；该模型 GitHub 标为 1.2B 参数（HF 元数据约 1B，F32）；token classification，**双链双头**（TokenCLS+/TokenCLS− 对应正/负链）。
- 标签系统：已公开文档仅描述**原核**任务为每链二分类（CDS=1 / non-coding=0，连续 1 为一个 CDS 区）；真核 preview 模型的具体标签集（是否含 intron/UTR/intergenic 多类）**未公开**（model card 写 "More technical details are coming soon..."）。
- 指标：**无任何公开的真核定量指标**；model card 仅有定性声明"superior accuracy compared to Augustus, Helixer, SegmentNT"，无数字、无评测物种。cds-annotation 数据集当前仅含 bacteria（57 行）与 archaea（42 行）子集，无真核评测子集（仅以 fly GCF_000001215.4 作示例输入）。
- 可得性：权重 MIT 许可可下载（发布后一月内约 3.3 万次下载）；代码公开；但真核评测复现风险高。

### 第 3 节 数据集与基准

- **GENCODE（人/鼠）**：来源 GENCODE/EMBL-EBI。最新人类 Release 49、小鼠 M36（参考构建 GRCh38、GRCm39）。Release 47/M36（2024-10）人类蛋白编码基因约 19,433 个（较早期版本因剔除早年 ab initio 误判而减少）。标签丰富（biotype、attribute；CDS/UTR/外显子/内含子/可变剪接/假基因/lncRNA/sRNA/免疫球蛋白与 TCR 片段）。下载：https://www.gencodegenes.org 。已知问题：可变剪接导致转录本 collapse 不一致；以参考为 ground truth 在性能接近参考时失真；近两年人鼠合计有超过 3.7 万模型被修改。GENCODE 论文：Frankish et al.（NAR 2023，53(D1):D966 为 2025 版）。
- **Ensembl 多物种**：来源 Ensembl/Ensembl Genomes，多 clade、多 genome build。Helixer、ANNEVO 大量使用。
- **RefSeq（NCBI）**：Tiberius 训练/评测所用（37 哺乳动物）。
- **Tiberius 训练数据**：37 个哺乳动物基因组（RefSeq 注释，软屏蔽），物种级留出 H. sapiens/B. taurus/D. leucas。Processed RefSeq 注释随论文提供。
- **Helixer 训练数据**：脊椎/无脊椎/陆生植物/真菌分 clade；初版即用 186 动物基因组+51 植物基因组做碱基级预测。2025 版测试物种含约 13 植物/11 脊椎/15 无脊椎，物种级留出评测。
- **SegmentNT 基准**：人类 GRCh38；chr20/21 测试、chr22 验证（染色体级留出）；元素注释来自 GENCODE 与 ENCODE cCRE（790k 增强子、34k 启动子）。
- **AUGUSTUS/BRAKER 基准集 & G3PO**：G3PO（Scalzitti et al. 2020，BMC Genomics，DOI 10.1186/s12864-020-6707-9）含 147 物种、人工校验真实基因，设多个子集评估装配质量、基因结构复杂度、蛋白长度等。
- **DNABERT-2 GUE 基准**：36 个数据集、9 类任务、4 物种，输入长度 70–10000 bp（偏分类而非完整基因结构）。
- **已知数据问题**：注释质量参差（参考本身是流水线产物）；同源泄漏（训练/测试物种系统发育过近）；染色体级泄漏（同物种染色体划分时 paralog 泄漏）；转录本 collapse 不一致；装配碎片化显著降低精度。

### 第 4 节 评估指标

- **碱基级 P/R/F1（逐类）**：对 CDS/intron/intergenic 等每类算 precision/recall/F1。Helixer 报告 base-wise F1（F1bw）。
- **外显子级 F1**：可为精确边界匹配（两端坐标全对）或重叠式。Tiberius/Helixer 用精确边界（TP=两端都对）。
- **基因级/转录本级 F1**：需明确——重叠阈值、链一致性、转录本/基因 collapse、CDS-only vs 全基因体、partial gene 处理。Tiberius 用 recall=TP/(TP+FN)、precision=TP/(TP+FP)，并指出领域内"specificity/sensitivity"有时与"precision/recall"混用。BRAKER 圈常用 transcript-level F1。
- **intron chain F1**：整条内含子链是否完全一致（Helixer 2025 Table 2 含此项）。
- **segment/interval F1**：SegmentNT 用 F1seg。
- **intergenic 假阳性率**：ab initio 方法在基因间区误报基因的比率（AUGUSTUS 典型问题，如 EGASP 中反向链误报基因）。
- **边界精度**：剪接供受体、起始/终止密码子定位准确率；GENERanno 原核任务报告 Start Accuracy/End Accuracy/Boundary Accuracy/Exact Match Rate。
- **macro-F1**：跨类平均，受少数类（内含子/剪接位点）影响大（如酵母 <1% 基因含内含子，F1 普遍偏低）。
- **BUSCO 完整度**：常用于新注释质量评估（Tiberius 96.0%、Helixer 92.1%、AUGUSTUS 74.2%）。
- **SOTA 特定评分**：Tiberius/BRAKER 用 TSEBRA/BRAKER 评测脚本；评测常忽略可变剪接（只取代表转录本）。

### 第 5 节 当前 SOTA 的已知局限

- **跨物种泛化差距**：Tiberius 不建议未重训练用于非脊椎动物；ANNEVO 指出 Tiberius 在 Vertebrate_other clade 上性能骤降。各模型需分 clade 训练。
- **长上下文建模**：多数模型用固定窗口（SegmentNT 30–50kb、GENERanno 8kb、DNABERT-2 更短）；人类基因可跨数十万 bp（Tiberius 正确预测最长 328,931 bp 已属罕见），超长内含子/基因仍困难。
- **基因间区假阳性**：纯分割模型缺结构化解码时易在基因间区误报。
- **边界精度**：剪接位点/起止密码子精确定位仍是难点；SegmentNT 存在与序列内位置相关的系统性偏差。
- **基因级一致性**：碎片化预测（一个基因被拆成多个、或多个基因被合并）。
- **罕见结构**：超长内含子、极短外显子、重叠基因、转剪接基因（transpliced genes）。
- **低质量装配**：碎片化装配显著降低精度（Helixer 在装配最碎的物种 Papaver somniferum、Triticum dicoccoides 表现最差）。
- **训练数据噪声**：参考注释本身含错误，限制以参考为 ground truth 的评测上限（ANNEVO 甚至在 31.8% 物种上 BUSCO 超过参考）。
- **推理成本**：基础模型（NT/GENERanno）推理昂贵；Helixer 比 Tiberius 慢约 5 倍。
- **对外部预训练/基础模型的依赖**：SegmentNT/GeneCAD/GENERanno 依赖大型预训练骨干，受其上下文与领域覆盖限制。

### 第 6 节 面向 Nature Methods 级贡献的开放机会

- **架构创新**：结构化解码器、可微分 HMM 头、semi-CRF（Tiberius 已做可微分 HMM，semi-CRF 在该任务仍少见）。
- **长上下文骨干**：Mamba/S4/扩张卷积建模 100kb+ 上下文，覆盖超长内含子与远程调控；PlantCAD2（8,192bp）、Evo2 等已探索单核苷酸长上下文。
- **基因级一致性的结构化解码**：把 CRF/latent-CRF/semi-CRF 与基础模型表示结合（GeneDecoder、GeneCAD 已起步，GeneCAD transcript F1 较 Helixer/BRAKER3 提升 8–10%，但缺大规模跨物种基准）。
- **跨物种表示学习**：进化感知预训练（ANNEVO 的联合进化建模、PlantCAD 的多基因组保守信号）。
- **生物学约束感知解码**：把 GT-AG 剪接、ORF/阅读框约束、起止密码子一致性写入解码层或损失。
- **基础模型适配策略**：参数高效微调、分割头设计、长序列重标定（SegmentNT 的 YaRN rescaling）。
- **多任务/多物种课程学习**：分 clade curriculum、多任务（剪接位点+外显子+基因体联合）。

### 第 7 节 生物学约束与后处理对比

| 约束 | Tiberius | Helixer(2025) | ANNEVO | SegmentNT | GENERanno(真核) | GeneDecoder/GeneCAD |
|---|---|---|---|---|---|---|
| GT-AG 剪接基序 | HMM 状态结构隐式编码 | HMM 内含子子状态显式（按剪接基序，60 子状态） | Viterbi 内含子状态群（GC-AG/GT-AG/AT-AC） | 仅预测剪接位点类，无强制 | 文档未述 | CRF 转移约束强制 |
| 起止密码子一致性 | 可微分 HMM 编码相位/密码子类型 | HMM CDS 子状态含 start/stop/regular | Viterbi 状态转移 | 无 | 原核仅 CDS 二分类 | CRF 特征顺序约束 |
| ORF/阅读框相位 | HMM 相位状态 | 10 个相位/密码子类型子状态 | 状态机相位 | 无 | 无（真核未述） | CRF latent graph |
| 链感知解码 | 是 | 是 | 是 | 各元素分别预测 | 双链双头 | 是 |
| 最小外显子/内含子长度 | HMM 长度约束 | HMM 子状态隐含 | 状态机 | 无 | 无 | CRF 可编码 |

- **整合层次**：Tiberius/Helixer/ANNEVO 把约束放在**模型架构/解码层**（可微分 HMM/Viterbi）；GeneDecoder/GeneCAD 放在**CRF 解码器**；SegmentNT/GENERanno 主要靠**后处理或根本不施加**。
- **改进边界过滤的空间**：用剪接基序与相位约束在解码阶段过滤非法外显子边界（GeneCAD 的蛋白语言模型 ORF 筛选即此思路，可锐化起止密码子与剪接连接），可显著降低基因间区假阳性与边界误差。

### 第 8 节 可复现性状态

- **Tiberius**：代码公开可用、权重随仓库提供、训练数据（processed RefSeq）有文档、评测脚本可得。难度：**中等**（TensorFlow>2.13、需 GPU）。**强烈建议先复现作为 baseline**。
- **Helixer**：代码公开、Zenodo 权重、训练数据多 clade 有文档、Web 界面。难度：**中等**。值得作为跨 clade baseline。
- **ANNEVO**：代码公开（非商业许可）、权重可下载。难度：**中等**。值得作为进化建模 baseline，但许可限制商业用途。
- **SegmentNT**：HuggingFace 权重+Colab notebook、代码公开（CC-BY-NC-SA）。难度：**易—中等**。适合作为基础模型分割 baseline，但需注意 30kb 上下文与位置偏差。
- **DNABERT-2**：HuggingFace 权重（zhihan1996/DNABERT-2-117M）、GUE 基准公开。难度：**易**。但非专门基因结构注释模型。
- **GENERanno（真核）**：权重可下载（MIT）、代码公开，但**无真核论文指标、无真核评测子集、真核标签系统未公开**。难度：**硬**（文档缺失）。**不建议作为首要复现对象**，可作为前沿对照关注后续正式发布。

---

## Recommendations

**阶段一（1–2 周，建立 baseline）**：先复现 Tiberius（哺乳动物）与 Helixer（跨 clade）作为两条主线 baseline，因二者代码/权重/数据/脚本最完整。用 RefSeq/GENCODE 物种级留出（held-out species）而非仅染色体级划分，以暴露真实跨物种差距。基准阈值：能在人类上复现 Tiberius 基因级 F1 ≥55%、外显子 F1 ≥88% 即视为环境正确。

**阶段二（3–6 周，定位空白）**：在统一基准（建议覆盖 ANNEVO 的五大 clade 子集 + 物种级留出）上同时跑 Tiberius、Helixer、ANNEVO、SegmentNT+分割头，统一用 transcript-level F1 + intron chain F1 + BUSCO，量化跨物种差距、基因间区假阳性、边界误差三大痛点。若发现长基因（>100kb）召回显著低于平均，则确认长上下文为最高价值方向。

**阶段三（核心贡献）**：押注"长上下文骨干 + 可微分结构化解码"。具体：用 Mamba/S4 或扩张卷积把上下文扩到 100kb+，接 semi-CRF/可微分 HMM 头显式编码 GT-AG/相位/ORF 约束，做跨物种 curriculum 预训练。决策阈值：若新方法在物种级留出上 transcript F1 较 Tiberius/ANNEVO 绝对提升 ≥5 个百分点且对长基因召回提升 ≥10 个百分点，即具备 Nature Methods 级新意。

**何时改变路线**：若结构化解码带来的提升 <2 个百分点，则转向"数据质量"路线（用 masked-motif score 等清洗参考注释，参照 GeneCAD）或"基础模型适配"路线（参数高效微调大型 DNA LM，如 GENERanno/NT/PlantCAD2）。

## Caveats

- 多处 F1 数字来自各模型自报，评测设置（是否忽略可变剪接、collapse 规则、留出方式）不完全一致，跨论文直接比较需谨慎；ANNEVO 对 Tiberius 的对比是在 ANNEVO 自身基准下，可能不利于 Tiberius（其 clade 范围限制）。
- GENERanno 真核注释器为 preview，定性"超越 Augustus/Helixer/SegmentNT"的声明无公开数字支撑，未经同行评审验证；其论文（DOI 10.1101/2025.06.04.656517）内容几乎全为原核模型，真核标签系统与指标均未公开——本节相关结论以"未公开"明确标注，而非推断。
- SegmentNT 的"gene annotation F1"与 Tiberius/Helixer 的"基因级 F1"定义不同（前者偏元素分割、后者偏完整基因模型），数值不可直接比较；SegmentNT 的 F1bw/F1seg 数字来自一项第三方基准（bioRxiv 2026.02.22.707219），非原论文。
- Helixer 2025 Nature Methods 版的具体 per-clade、per-species F1 表（Table 2）本报告仅得其结构与部分数值，精确逐物种数字需查原文 Table 2。
- 部分预印本（Tiberius 多 clade 版、GeneCAD、SegmentNT 偏差分析）尚未经同行评审，结论可能调整。Tiberius/Helixer/ANNEVO 的速度对比（约 1h39m vs 8h54m）来自 Tiberius 作者设置，会因硬件而异。
