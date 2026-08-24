# Deep Research Report — chatgpt · 2026-06-09

> 占位文件，由 `$research-interview` 在 2026-06-09 自动生成。
>
> **使用方式**: 把外部平台（chatgpt）跑出来的完整 deep research 报告整段粘贴到下面
> "## Report body" 之下。可以保留这段说明也可以删除，由 `$research-synthesize` 自动跳过。
>
> 对应的 prompt: `docs/inputs/deep_research_prompt_2026-06-09.md`

---

## Source platform
- Platform: chatgpt
- Mode: Deep Research
- Run date:
- Number of retries / refinements:
- Subjective quality (filled later): low / medium / high

## Report body

# 跨物种真核蛋白编码基因从头注释深度学习研究现状报告

## 总体判断

如果把任务严格限定为“只输入原始基因组 DNA、输出逐碱基蛋白编码基因结构标签、并且要求跨物种泛化”，那么截至 2026 年中，公开证据最强、最接近你所说的 *Nature Methods* 级比较基准的候选，主要是 **ANNEVO** 与 **Tiberius**；**Helixer** 仍然是必须复现的强公开基线；**SegmentNT** 很有研究价值，但它的原始任务定义更接近“多标签基因组语义分割”，而不是传统 gene caller 的全基因结构重建；**GENERanno-eukaryote-1.2b-cds-annotator-preview** 很值得关注，但目前公开的真核专项评测与复现文档还不够完整，尚不足以像 ANNEVO/Tiberius/Helixer 那样做严格可审稿比较。这个判断来自论文、官方仓库、模型卡和当前公开 benchmark 的综合比对，而不是单看某一篇论文的单一数字。citeturn9view0turn41view0turn11view3turn18view0turn30view0turn35view1

对你的具体研究目标而言，**最关键的设计选择不是“要不要用 foundation model”本身，而是“是否能在 whole-genome、held-out species、严格控制 intergenic false positives 的协议下，把局部逐碱基分类提升为全基因一致的结构解码”**。当前文献已经很清楚地显示：只做局部 per-base segmentation 不足以代表 gene annotation SOTA；真正拉开差距的是长程上下文、跨物种表示、以及把 phase / splice / start-stop / ORF 约束并入 decoder 或训练目标的能力。ANNEVO 把这件事做到了“MoE evolutionary model + Viterbi”，Tiberius 做到了“CNN/LSTM + differentiable HMM”，Helixer 做到了“deep basewise prediction + HelixerPost HMM”，而 SegmentNT 则展示了 foundation model 做高分辨率分割的上限，但其原生输出并不是传统意义上的完整 gene model。citeturn26view0turn44view1turn41view0turn18view0turn20view0turn30view0

换句话说，如果你的最终目标是“在可比 benchmark 上严格超过已发表 deep-learning SOTA”，那么**最可信的主赛道**应该瞄准 ANNEVO/Tiberius/Helixer 这条线：跨物种、whole genome、带结构解码、输出 gene model；而 **SegmentNT 更适合作为 backbone / pretraining / auxiliary head / long-context 表示学习的技术来源**，不太适合作为唯一最终对手；**GENERanno preview 则更像一个潜在强对手，但当前公开证据链仍不够闭环**。citeturn11view3turn41view0turn17view0turn30view0turn35view0turn35view1

## 方法谱系

这个问题大致可以分成六个方法家族。为了避免把“传统 gene finder”和“现代深度学习 segmentation/decoding”混在一起，下面我按**核心建模假设**来分，而不是按年代来分。每个家族我都给出代表性工作、长处、短板，以及你真正需要关心的“典型可达到水平”。

**自训练概率图模型家族**的核心思想，是把基因结构写成状态机或 HMM/GHMM，由“信号传感器”去识别起始密码子、终止密码子、剪接位点、密码子统计特征，再用动态规划求最优路径。这个方向里，仍然有历史和方法学价值的代表作包括 Alex Lomsadze 的 *Gene identification in novel eukaryotic genomes by self-training algorithm*（2005，*Nucleic Acids Research*）和 V. Ter-Hovhannisyan 的 *Gene prediction in novel fungal genomes using an ab initio algorithm with unsupervised training*（2008，*Genome Research*）。它们的长处，是**强可解释性、很强的生物学先验、对缺数据物种仍可工作、输出天然结构一致**；短板，是**特征工程重、跨远缘物种时泛化差、对长程依赖和复杂替代剪接受限**。在今天的公开比较里，它们仍然是重要 baseline，但在跨物种 held-out benchmark 上，尤其是长内含子真核生物中，已经更多扮演“必须打败的传统基线”而不是“最终 SOTA”。citeturn18view0turn26view0turn30view0

**早期跨物种深度逐碱基分类家族**的代表是 Felix Stiehler 的 *Helixer: Cross-species gene annotation of large eukaryotic genomes using deep learning*（2020，*Bioinformatics*）。这一家族的核心是：把 gene annotation 先视为 per-base labeling，再尝试利用跨物种训练获得更好的表示。它相对传统 HMM 的进步在于**自动学习局部序列模式**，对训练物种之外的 genome 更稳健；但当年的模型还没有很好解决全基因一致性与长程依赖问题，因此更像是“把信号识别从 hand-crafted sensor 换成神经网络”。它的重要性在于开创了**跨物种 ab initio deep learning gene annotation**这个子方向。citeturn17view0

**深度逐碱基预测加独立 HMM/后处理家族**的成熟形态是 Felix Holst 的 *Helixer: ab initio prediction of primary eukaryotic gene models combining deep learning and a hidden Markov model*（2025，*Nature Methods*）。这一家族把 deep net 和结构化后处理明确分层：神经网络输出 UTR/CDS/intron/intergenic 等逐碱基概率，再交给 HMM/HelixerPost 生成最终 gene model。它的优势是**工程上可用、标签体系丰富、开放模型多、跨多个谱系有现成权重**；而缺点也很明确：**模型和输入长度按谱系分开调、速度较慢、在长基因和困难植物 benchmark 上明显掉点、HMM 部分仍依赖手工设计 penalty**。从当前公开比较看，Helixer 仍然是最重要的公开 deep baseline 之一，但在 2026 年的跨谱系头部比较中，已经落后于 ANNEVO 和当前版本 Tiberius。ANNEVO 官方仓库在 12 个模式物种上的复评给出 Helixer 平均 exon recall/precision 为 86.1/75.3，locus recall/precision 为 50.2/47.0，BUSCO 为 92.5。这个数字不是 Helixer 论文原文的主表数字，而是 2026 年最新版本之间的横向比较，因此我建议你把它视作“当前工程版本参考水平”，而不是论文原始结果。citeturn18view0turn17view0turn20view0turn11view3

**端到端深度模型加可微分结构解码家族**的代表是 Lars Gabriel 的 Tiberius。官方仓库把它概括为“CNN/LSTM 与 differentiable HMM layer 的端到端整合”，并在 2026 年扩展到了 Mammalia、Vertebrates、Insecta、Angiosperms、Fungi、Diatoms、Chlorophyta 等多个 clade。这个家族相较 Helixer 的关键创新，是**把结构一致性更深地并入模型本体，而不是单纯后接一个外部 HMM**。优点是**gene-level consistency 更强、对复杂基因结构更友好、已有多谱系现成模型**；弱点是**依然依赖 clade-specific 权重，且当前最公开、最完整的原始论文首先是在哺乳动物上建立声誉**。按 ANNEVO 仓库 2026 年的统一复评，Tiberius 在 12 物种上的平均 exon recall/precision 为 89.8/88.7，locus recall/precision 为 74.0/68.8，BUSCO 为 96.3，明显强于 Helixer。citeturn41view0turn42view0turn11view3

**foundation model 高分辨率语义分割家族**的代表是 B. P. de Almeida 的 *Annotating the genome at single-nucleotide resolution with DNA foundation models*（2025，*Nature Methods*），对应模型 **SegmentNT**。它把问题显式表述为 multilabel semantic segmentation，使用预训练 DNA foundation model（NT）加 1D U-Net，对 14 个 genomic elements 做逐碱基分割。优势是**局部标签分辨率高、可自然兼容多任务、很适合迁移 backbone、对 splice site / exon / intron 这类局部结构很强**；缺点是**原生输出不是完整 gene model、全基因一致性依赖额外后处理、最初 benchmark 是人类染色体切分而非传统 whole-genome held-out species gene finder 协议**。它在 human 14-element benchmark 上的 best model SegmentNT-30kb 平均 MCC 为 0.45，在 50 kb 推理时可到 0.47；其 multispecies 版本在 human-distant 动物上由 0.49 提升到 0.57，在植物上由 0.34 提升到 0.45。这非常强，但**不应直接当作 ANNEVO/Tiberius/Helixer 那种 gene caller 指标的同协议 SOTA**。citeturn30view0turn32view0turn32view1turn29view1

**跨谱系进化建模与专家混合家族**的当前代表是 Pengyu Zhang 的 *Highly accurate ab initio gene annotation with ANNEVO*（2026，*Nature Methods*）。它的核心是把 gene annotation 视为**“长程序列依赖 + 跨谱系进化关系 + 结构解码”**的联合问题：卷积塔形成 bins，Transformer 建模远距关系，MoE 建模不同亚谱系，multi-task 头预测 category / phase / state，再用 Viterbi 得到符合生物学规则的 gene structure。优点是**跨 566 物种 benchmark 的证据最强、在不使用外部实验数据时已接近甚至匹配 evidence-assisted pipelines、当前开源版本更新很快**；弱点是**仍然采用 clade-level 建模、解码与工程复杂度高、许可证并非 OSI 开源**。2026 年当前仓库版本在 12 个模式物种上给出平均 exon recall/precision 91.4/90.2、locus recall/precision 76.3/74.3、BUSCO 97.8；而论文版在 566 RefSeq 物种上相对优化版 Augustus 的平均提升为 nucleotide-level F1 提升 7.5–22.3%、gene-level F1 提升 9.9–38.5%、BUSCO 提升 9–34.5%。这两组数字来自不同 benchmark 和不同软件版本，不应该直接混算，但足以说明它是当前最强公开 ab initio 候选之一。citeturn9view0turn26view0turn44view1turn11view3

**大规模 genomic foundation model 下游 CDS 注释家族**目前公开可见的代表是 Qiuyi Li 的 *GENERanno: A Genomic Foundation Model for Metagenomic Annotation*（2025，bioRxiv）及其后续发布的 **GENERanno-eukaryote-1.2b-cds-annotator-preview**。这一路线的吸引力在于：把大规模 eukaryotic DNA pretraining 与 token-classification 式 CDS annotation 结合，潜在上能以统一 FM 做跨物种迁移。优点是**模型规模大、预训练数据大、对你想做“foundation-model adaptation + structured decoder”的方向很有启发性**；短板是**截至目前，真核 preview 模型的公开评测、分割脚本和可审稿比较表还明显不如 ANNEVO/Tiberius/Helixer 完整**。因此我把它列为“高潜力方向”，而不是“当前最稳妥的主比较基线”。citeturn35view0turn35view1turn35view2

## 主要 SOTA 候选

**Tiberius**  
论文（原始哺乳动物版本）：`https://doi.org/10.1093/bioinformatics/btae685`。多谱系训练预印本入口由官方仓库指向 2026 clade-training preprint：`https://doi.org/10.64898/2026.04.24.720536`。官方代码：`https://github.com/Gaius-Augustus/Tiberius`。权重入口：`https://github.com/Gaius-Augustus/Tiberius/tree/main/model_cfg`；仓库明确说明每个 YAML 都包含公开 `weights_url`，并列出当前可用 clade：Mammalia、Vertebrates、Insecta、Eudicotyledons、Monocotyledonae、Fungi、Diatoms、Chlorophyta。官方当前接口可直接通过 `--model_cfg` 自动解析下载权重。就可核验的公开数字而言，我本次检索里最稳妥的是 2026 年 ANNEVO 仓库的统一复评：12 个模式物种平均 exon recall/precision 89.8/88.7，locus recall/precision 74.0/68.8，BUSCO 96.3。性质上它显然是**cross-species、clade-aware**，而不是 same-species 训练。需要特别提醒的是：Tiberius 2024 论文版、2026 多谱系预印本版、以及当前仓库版不是同一个冻结系统，因此你做对比时必须明确版本、clade 和是否 softmasking。citeturn41view0turn42view0turn11view3

**Helixer 2025 版**  
论文：`https://doi.org/10.1038/s41592-025-02939-1`。官方代码：`https://github.com/usadellab/Helixer`。当前官方权重总入口：`https://zenodo.org/records/17850139`；仓库推荐用 `scripts/fetch_helixer_models.py` 自动下载。仓库 README 列出的经典 lineage 模型文件包括 fungi、land_plant、vertebrate、invertebrate；Nature Methods 论文与后续 Zenodo 记录又补充了 mammal 模型。Helixer 的原生标签体系是 **intergenic / UTR / CDS / intron**，随后经 HMM 产出主转录本 gene models。对于“当前 SOTA”这一问题，最可比较的公开数字有两套：其一是它自己论文补充材料中各测试物种的 Genic F1、Subgenic F1 和 exon-level gffcompare 指标；其二是 2026 年 ANNEVO 仓库基于 12 个模式物种、最新版本软件的横向复评，得到平均 exon recall/precision 86.1/75.3、locus recall/precision 50.2/47.0、BUSCO 92.5。两套数字协议不同，所以不能混成一个总分，但都支持同一个结论：**Helixer 仍是强公开 baseline，但已不是最强公开 ab initio 结果**。训练/测试方面，它是**同 lineage 内的 cross-species model**，不是单物种拟合。citeturn17view0turn18view0turn16search2turn20view0turn24view0turn24view1turn11view3

**ANNEVO**  
论文：`https://www.nature.com/articles/s41592-026-03036-7`；开放预印本 PDF：`https://doi.org/10.21203/rs.3.rs-6402260/v1`。官方代码：`https://github.com/xjtu-omics/ANNEVO`。官方权重目录：`https://github.com/xjtu-omics/ANNEVO/tree/main/saved_model`，其中已公开 `ANNEVO_Mammalia.pt`、`ANNEVO_Insecta.pt`、`ANNEVO_Aves.pt`、`ANNEVO_Actinopteri.pt`、`ANNEVO_Magnoliopsida.pt`、`ANNEVO_Fungi.pt`。论文版在 566 个 RefSeq 物种、五大 clade 上，相对优化版 species-specific Augustus 的平均提升为 nucleotide-level F1 提升 7.5–22.3%，gene-level F1 提升 9.9–38.5%，BUSCO 提升 9–34.5%；在 *Sus scrofa* 上，论文图注给出 CDS nucleotide F1 为 0.89、gene completeness 为 78%、BUSCO 为 91.6。与此同时，仓库当前 v2.3.x 版又进一步宣称人类 BUSCO 从论文版的 95.7 提高到 98.3，并在 12 个模式物种上给出平均 exon recall/precision 91.4/90.2、locus recall/precision 76.3/74.3、BUSCO 97.8。这里存在一个必须明确标注的**版本冲突**：论文版和当前 GitHub 发布版不是同一冻结模型，因此如果你要和 ANNEVO 比，最好比“论文冻结版”和“当前工程版”两种设置。训练/测试上，ANNEVO 是**cross-species、clade-level**，不是 same-species。citeturn9view0turn26view0turn44view1turn44view2turn11view3turn13view0turn14view0

**SegmentNT**  
论文：`https://www.nature.com/articles/s41592-025-02881-2`。代码与官方文档入口：`https://github.com/instadeepai/nucleotide-transformer` 与 `https://github.com/instadeepai/nucleotide-transformer/blob/main/docs/segment_nt.md`。权重：`https://huggingface.co/InstaDeepAI/segment_nt`；跨物种版权重：`https://huggingface.co/InstaDeepAI/segment_nt_multi_species`。它的原始任务是**14 类 genomic elements 的逐碱基多标签分割**，并不等同于传统 gene caller 的 GFF3 重建；人类版训练在除 chr20/21/22 之外的所有染色体上，test 为 chr20/21，val 为 chr22，并移除了与训练/验证基因同源的 test chunks；multispecies 版在 human + mouse/chicken/fly/zebrafish/worm 的 genic labels 上再微调，并对另外 10 个动物和 5 个植物做 held-out species 测试。核心数字方面，SegmentNT-30kb 在 human 14-element benchmark 上平均 MCC 为 0.45，在 50 kb 推理时可达 0.47；在跨物种上，human model 对 human-close species 平均 MCC 为 0.62、human-distant 为 0.49，而 multispecies 版分别为 0.64 和 0.57；植物上 human model 为 0.34、multispecies 为 0.45。作者还报告在 whole-chromosome、all-isoform 的 human gene benchmark 上，SegmentNT-30kb 在各 gene elements 上显著优于 AUGUSTUS。我的判断是：**如果你要做“从 FM 到 gene caller”的论文，SegmentNT 是极好的 backbone 参照；但如果你要做“严格可比的 ab initio protein-coding gene annotation SOTA”，它不是最标准的主比较对象。**citeturn30view0turn32view0turn32view1turn33view0turn29view1turn31search0

**GENERanno-eukaryote-1.2b-cds-annotator-preview**  
当前可核验的主要公开入口不是一篇真核专项完整论文，而是官方仓库和 Hugging Face 模型卡。仓库：`https://github.com/GenerTeam/GENERanno`。模型卡：`https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`。相关通用论文入口：`https://www.biorxiv.org/content/10.1101/2025.06.04.656517`。官方仓库说明该 preview 是一个**1.2B 参数、面向 eukaryotic CDS annotation 的 expert model**，并把它列入 `GENERanno` 系列；模型卡声称它在综合评估中优于 Augustus、Helixer 和 SegmentNT，但同时也明确写着 “More technical details are coming soon”。仓库还给出 eukaryotic CDS annotation 的 CLI 与测试数据集入口 `GenerTeam/cds-annotation`，以及一个受访问条款限制的 `annotation_data_eukaryote` 数据集；然而**我在本次检索中没有找到像 ANNEVO/Tiberius/Helixer 那样完整的、可直接重算的真核公开数值表、冻结 split 说明和官方 metric script**。因此，当前最稳妥的结论只能是：**它是高潜力候选，但公开可审计证据尚不足以把它排在 ANNEVO/Tiberius 之前。** 关于“cross-species 还是 same-species”，官方定位显然是跨物种真核注释，但公开 benchmark 文档尚不足以做细粒度审计。citeturn35view0turn35view1turn36view2turn37search2turn37search6

综合来看，如果你今天就要选三个“最值得先复现并超过”的公开对象，我会把优先级排成：**ANNEVO（当前最强公开证据） > Tiberius（结构化 deep caller、极强基线） > Helixer（标准公开 deep baseline）**。SegmentNT 适合作为架构来源与 supplementary baseline；GENERanno preview 适合持续跟踪，但不建议把它作为你稿件里唯一的头部对手，除非其作者后续公开了真核专项完整评测。这个排序是我基于公开证据链做的推断，不是任何单一论文的原文结论。citeturn11view3turn41view0turn17view0turn30view0turn35view1

## 数据集与基准协议

这一领域最大的现实问题之一，不是模型本身，而是**benchmark 并不统一**。表面上大家都在做“gene annotation”，但实际上存在三类明显不同的协议：**多物种 whole-genome/whole-chromosome gene caller 协议**、**human-centric chromosome split 语义分割协议**、以及**foundation-model 下游 CDS 标注协议**。如果你将来要做可被 *Nature Methods* 接受的主 benchmark，最好自己构建一个明确冻结、可复现、跨物种 held-out、whole-genome 的比较体系，并且保留 Helixer/Tiberius/ANNEVO 所采用的那种“intergenic 区域也要计分”的设置。citeturn20view0turn26view0turn30view0turn35view0

**Helixer/Tiberius 这条主线的多物种 benchmark**，公开信息最完整的是 Helixer 2025 补充材料。它给出了 assessed species 的总量：fungi 298、plants 77、vertebrates 314、invertebrates 201；并给出 seeds per split 与 training species per split，说明它不是简单的随机 chunk 切分，而是**以物种为单位的 train/val/test 组织**。具体数据源为：fungi training/validation/test 来自 RefSeq（2022-03-04）；plants 的 train/val 来自 Phytozome13（2021-06-07），plant test 来自 RefSeq（2022-07-14）；vertebrates 和 invertebrates 来自 RefSeq（2022-05-06）；mammals 则来自 RefSeq（2025-03-13），并明确写明“采用与 Tiberius mammal model 相同的 species selection and partition”。补充表 S17–S21 还列出了物种名、assembly accession、版本与 split。这个协议的最大优点，是**真 held-out species、whole-genome 标注、真正对应 gene finder 任务**；但已知问题也很多：作者承认为了与其它工具可比，test species 部分是从“有接近 Augustus 模板的物种”和“额外每 N 个选一个物种”中构造的，这会带来**benchmark convenience bias**；plants 的 train 和 test 来自不同数据库，也可能引入**源域偏移**；而 RefSeq 标注本身并非无噪声金标准。citeturn20view0turn25view0turn25view1turn25view2turn25view3turn25view4

如果你要复现这一主线，最实用的“入口 URL”不是某一个固定压缩包，而是数据库入口页与补充表清单：RefSeq/NCBI 官方框架页可从 `https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process/` 进入；Helixer 论文的冻结物种/版本清单则在补充 PDF 中；Tiberius 当前公开的模型配置和 clade 权重入口位于 `https://github.com/Gaius-Augustus/Tiberius/tree/main/model_cfg`。这里我要明确说明不确定性：**论文并没有提供一个一键冻结所有训练/测试 genome 的单一官方下载脚本**，因此复现必须依赖补表中的 accession/version 清单自行重建。citeturn34search13turn20view0turn41view0turn42view0

**ANNEVO 的主 benchmark**使用 566 个 RefSeq 物种，覆盖 Fungi、Embryophyta、Invertebrates、Vertebrate_Mammalia、Vertebrate_other 五个 clade；作者在方法部分明确写到，他们先比较了 RefSeq 和 Ensembl 的 BUSCO 表现，最终选择 RefSeq 作为主训练/主评测来源，而 Ensembl 仅用于补充性 completeness 分析。对 metazoan lineages，他们又从 candidate test sets 中随机抽取 50 个物种，以兼顾系统发育多样性和计算可行性。更重要的是，ANNEVO 把任务定义为**最长 coding transcript 的 sequence-based prediction**，这意味着它在标签定义上更接近“主 CDS 基因结构”的统一目标，而不是完整 alternative isoform 生态。这个 benchmark 的强项是**规模大、跨谱系、whole-genome 意义明确**；但它同样有已知问题：物种选择中显式借助了 BUSCO 质量筛选，这会提升参考注释质量，却也可能**低估现实世界低质量注释/低质量组装的难度**。可复现的数据库入口可从 RefSeq 与 Ensembl 官方站点进入，而作者的冻结协议仍以论文与补充说明为准。citeturn26view0turn27view0turn27view3

**SegmentNT 的协议**与上述 gene caller 主线显著不同。它的 human segmentation dataset 来自 GENCODE 与 ENCODE，其中 gene elements 包括 protein-coding genes、lncRNAs、5′UTR、3′UTR、exon、intron、splice acceptor/donor，regulatory elements 包括 poly(A) signal、tissue-invariant / tissue-specific promoters / enhancers 以及 CTCF-bound sites。作者按染色体切分：human test 用 chr20/21，validation 用 chr22，其余用于训练，并通过 Ensembl BioMart 去掉 test 中与 train/val 基因同源的 chunks，以减少泄漏；但作者也明确承认，这并**不能移除同源 distal regulatory elements 的潜在泄漏**。它的跨物种协议则是 human + 五个训练物种，再用十个动物和五个植物做 held-out species。下载入口方面，GENCODE 人类注释总入口是 `https://www.gencodegenes.org/human/`，ENCODE cCRE 注册表入口是 `https://screen.encodeproject.org/`。这个体系非常适合作为“高分辨率 DNA foundation model segmentation”基准，但和 gene caller 文献中的 whole-genome gene-level F1 并不是同一类 benchmark。citeturn30view0turn33view0turn33view1turn28search8

**GENERanno 的真核数据公开程度目前最低**。仓库公开给出了下游 `cds-annotation` 数据集入口 `https://huggingface.co/datasets/GenerTeam/cds-annotation`，以及一个带访问条件的 `annotation_data_eukaryote` 数据集；另外还有开放的 eukaryotic pretraining corpus `https://huggingface.co/datasets/GenerTeam/pretrain_data_eukaryote`。但对你要做的严格 benchmark 来说，关键问题是：**真核 preview 模型的官方 split、species manifest、label schema、metric script 和 frozen baseline numbers 目前还没有像 Helixer/ANNEVO 那样公开完整说明**。因此我建议把 GENERanno 相关数据用于“探索性补充实验”而非“主 benchmark 锚点”。citeturn35view0turn36view2turn37search0turn37search2turn37search6

如果你要自己设计一个更好的 benchmark，我的建议是：**以 species-held-out 为主切分单元，以 whole chromosome / whole genome 为评测单元，以 protein-coding gene body + exon/intron/CDS 作为主标签层级；对所有方法统一 transcript collapsing 规则，并同时给出 main-isoform 与 collapsed-gene 两个版本。** 这不是现成论文文字，而是我基于现有 benchmark 差异做的归纳性建议。支持这一建议的直接原因，是现有文献已经反复显示：仅在局部 windows、只看 BUSCO、或只在 human chromosome split 上比，是不够支撑 gene annotation SOTA 结论的。citeturn26view0turn30view0turn20view0turn27view3

## 评价指标与评分规则

在这个任务里，**最容易“看起来很高、其实没用”的指标是 BUSCO 和只在 gene-rich windows 上计算的局部 F1**。ANNEVO 明确指出 BUSCO 不能反映 false positives，而 gene-level F1 能提供全基因组范围内、同时考虑 false positives 的补充视角。SegmentNT 也因此不仅报告局部 segment 指标，还专门在 whole test chromosomes 上比较 gene elements。对你这种想控制 intergenic 假阳性的工作，**真正该当主指标的，是 whole-genome / whole-chromosome 条件下的 genic F1、gene/locus precision-recall-F1，以及 exon/intron/CDS 的精确结构指标**。citeturn27view3turn33view0

当前文献中最常见的指标可以分成五层。第一层是**逐碱基二分类或多分类指标**，例如 base-level precision / recall / F1、genic F1、subgenic F1；Helixer 补充材料给出了 Test species 的 Subgenic F1 与 Genic F1 表。第二层是**结构单元级指标**，典型是 exon、intron、CDS 的 precision / recall / F1，以及 splice donor/acceptor 的识别准确度；SegmentNT 特别强调 splice site 元素，Helixer 和 ANNEVO/Tiberius 类工作则更常以 exon/intron/gene/locus 为单位。第三层是**gene/locus 级指标**，ANNEVO 当前官方仓库用 gffcompare 的 exon recall/precision 与 locus recall/precision 做 12 物种的统一对比。第四层是**区域级重叠指标**，例如 SegmentNT 使用的 SOV。第五层是**辅助完整度指标**，如 BUSCO。citeturn24view0turn24view1turn11view3turn30view0

这里最需要精确定义的是 **gffcompare 风格的“gene/locus matching”**。根据 GffCompare 论文：在 **transcript level**，多外显子转录本的 TP 被定义为“full exon chain match”，要求所有内部外显子完全匹配，终端外显子边界默认最多允许 100 bp 的轻微差异；单外显子转录本则要求和参考单外显子转录本有“显著重叠”，默认是长者长度的 80% 以上。到了 **locus level**，GffCompare 把 locus 定义为**基于 exon overlap 聚成的一簇转录本**，并认为一个预测 locus 和一个参考 locus 匹配，当且仅当该 locus 中至少一个预测 transcript 在 transcript level 上匹配到相应参考 locus 中的某个参考 transcript。源码文档进一步显示 locus 聚类是在**同一链上的重叠 mRNAs**上完成的。也就是说，**strand consistency、内含子链一致性、终端边界容差、单外显子显著重叠阈值**，都是 gene-level / locus-level F1 的实际组成部分；你在写论文时必须把这些规则写清楚，否则“gene F1”四个字没有统一含义。citeturn45search1turn45search4

需要特别强调的是，**不同 SOTA 论文对“gene-level”并不完全同义**。ANNEVO 论文正文把 gene-level F1 与 gene completeness 并列使用，并在部分图中比较 `Gene precision / Gene recall`；当前仓库又使用 `locus recall / locus precision`。Helixer 的长期传统是 `Genic F1 / Subgenic F1` 加 gffcompare 的 exon-level 或其他结构结果。SegmentNT 则在 human gene benchmark 中对“main isoform only”“all confident isoforms”“whole test chromosomes”分别给 F1、MCC 和 SOV，但这本质上还是**per-base multi-label segmentation**的评价，而不是输出最终 GFF3 后用 gene matching 规则计分。因此，你未来如果要做“严格超过 SOTA”，**必须主动统一 scoring rule，而不是沿用各文献自家最有利的定义。** 这是我基于现有论文差异做的归纳判断。citeturn26view0turn11view3turn24view0turn33view0turn32view3

对你的任务，我建议把 scorecard 排成这样：主指标用 **whole-genome genic F1** 与 **gffcompare locus/transcript F1**；并列报告 **exon/intron/CDS precision-recall-F1**、**splice donor/acceptor precision-recall**、**boundary error distribution**、**fragmentation / fusion counts**、**BUSCO completeness**；所有结果至少给出 **collapsed-gene** 与 **main-isoform** 两个版本。这个建议不是直接来自某一篇论文，而是综合了 Helixer、ANNEVO、SegmentNT 与 gffcompare 的共同经验：只有这样，你才能既控制 intergenic 假阳性，又衡量结构边界精度，还能避免被 alternative isoform 定义牵着走。citeturn27view3turn24view0turn33view0turn45search1

## 当前瓶颈与高价值机会

先说**文献已经明确暴露出来的失败点**。第一，**跨物种泛化仍然强依赖 clade 划分**。Helixer 官方现成模型按 fungi / land_plant / vertebrate / invertebrate 分开，ANNEVO 按五大 clade 建模，Tiberius 2026 版也按 Mammalia / Vertebrates / Insecta / Plants / Fungi / Algae-Protists 分提供权重。ANNEVO 论文甚至直接批评 Helixer 的“clade-specific settings”“custom input lengths and architectures”限制了 unified cross-clade model 的发展。这意味着“一个统一真核 ab initio 模型”在公开 SOTA 中仍然没有被真正解决。citeturn17view0turn41view0turn26view0

第二，**长程上下文仍然是主要瓶颈**。ANNEVO 把可建模上下文写到 40 kb，并专门把 Mammalia 的改进归因于对长基因的更好建模；SegmentNT 通过 RoPE rescaling 把有效推理长度推到 50 kb，并观察到 longer context 对 gene elements 明显有益；但从真实真核基因结构看，超长内含子、超长基因、远端剪接线索、重复区上下文，很多时候仍然超过这些尺度。ANNEVO 用 human TRIO 365 kb 作为典型例子，说明长基因重建仍是拉开差距的关键场景。你的题目如果要做头部工作，**长上下文 backbone 仍然是最自然的高价值切入点**。citeturn44view1turn27view2turn27view4turn32view0

第三，**intergenic false positives 与 gene-level consistency 仍未被完全统一处理**。ANNEVO 明确强调 BUSCO 不反映 false positives；SegmentNT 也因此在 whole-chromosome 测试上单独比较；Helixer 通过 HMM 过滤 basewise 输出；Tiberius 则把结构化层嵌入模型本体。但从结果形式上看，当前方法仍普遍存在两个未彻底解决的问题：一是**把局部高概率片段错误地拼成碎片基因**，二是**对长 intergenic 区域的假阳性控制仍依赖解码器和阈值工程**。这也是为什么你若只做 segmentation head 改进，很难在真正的 gene-level benchmark 上稳赢。citeturn27view3turn33view0turn20view0turn41view0

第四，**参考注释噪声和组装质量问题会反向限制模型上限**。ANNEVO 专门在训练时对“pre-identified erroneous regions”做 loss masking，并展示它可以修正 Ensembl/RefSeq 中的结构错误；SegmentNT 用同源过滤减少 leakage，却也承认远端调控同源仍可能残留；ANNEVO 仓库还连续发布了针对**高度碎片化组装**的内存优化版本，并给出一个 59,693 条 contigs 的鸟类 RefSeq 组装作为例子。这说明今天的高分模型，至少部分性能，仍然受制于训练标签噪声与 assembly fragmentation，而不只是网络表达能力。citeturn26view0turn27view3turn33view0turn9view1turn9view2

基于这些公开局限，我认为最像 *Nature Methods* 级贡献的机会主要有五类。**第一类是长上下文 backbone 创新**：不是简单把窗口从 40 kb 提到 80 kb，而是把“局部 motif/phase + 中程 exon-intron 依赖 + 超长基因全局一致性”分层建模。SegmentNT 证明了 foundation model 表示对局部元素有效，ANNEVO 证明 40 kb 和 MoE 有用，Tiberius 证明结构层必须进模型；真正的新贡献，很可能是把这些优点统一到一个**长上下文、跨谱系、可解码**的架构里。这个判断是推断，但直接由三条技术线的公开结果支撑。citeturn30view0turn26view0turn41view0

**第二类是 decoder/head 的结构化创新**。当前最强公开方法并没有把所有经典约束都变成可学习的、统一优化的目标。HelixerPost 仍有手工 transition penalties；ANNEVO 虽然把 category/phase/state 学进去了，但解码逻辑仍是分步式；SegmentNT 则基本没有原生 gene-structure decoder。非常自然的机会，是构造一个**生物约束感知的 structured decoder**：显式控制 start/stop、reading frame、splice motif、minimum intron length、同链一致性、ORF 完整性、gene overlap 约束、以及 partial gene 规则，并把这些约束写进训练目标或可微分动态规划中。citeturn20view0turn26view0turn30view0

**第三类是专门针对 intergenic false positives 的目标函数设计**。SegmentNT 和 ANNEVO 都间接说明 whole-chromosome/whole-genome 评测的重要性，但现有模型多数仍以局部 label loss 为主。你完全可以把研究重点放在**whole-genome negative mining、intergenic calibration、FP-aware loss、fragmentation penalty、fusion penalty、boundary-aware curriculum** 上。只要协议设计得好，这种“不是更大模型，而是更对的 objective”反而更容易做出可审稿的强结论。这个机会更多是我的推断，但与现有方法的薄弱点高度对齐。citeturn27view3turn33view0turn20view0

**第四类是跨物种表示学习与 taxonomy-aware adaptation**。ANNEVO 的 MoE 已经表明“把进化关系显式纳入建模”有效；SegmentNT multispecies 也表明从 human-only 到 multispecies 的微调能显著改善远缘泛化；GENERanno/GENERator 则说明超大规模 eukaryotic pretraining 在表示层面很有潜力。因此，一个前沿但仍然可做的方向是：**foundation model 负责跨物种表示，轻量 taxonomy-aware adapter / MoE / routing 负责谱系差异，structured decoder 负责 gene consistency**。如果你能在 held-out species 上证明这个组合同时优于 ANNEVO 和 Tiberius，那会非常有说服力。citeturn26view0turn32view1turn35view0turn35view2

**第五类是 annotation-noise robustness 与 cross-reference consensus learning**。ANNEVO 已经给出一个重要信号：模型不仅可以拟合参考注释，还可以纠错参考注释。对一个 aiming high 的工作来说，最有价值的不是“更像 RefSeq/Ensembl”，而是“在 noisy annotations 下学到更稳定的 coding-gene grammar，并能在 orthogonal evidence 下更常为真”。如果你能系统地把 RefSeq/GENCODE/Ensembl 的冲突区域、低置信区域、isoform collapsing 不一致区域纳入训练和评测，这会比单纯追逐 1–2 个百分点更像高影响力方法学贡献。citeturn27view3turn26view0turn30view0

## 生物学约束、后处理与复现现状

在“生物学约束是否被真正用起来”这件事上，当前方法之间差异非常大。**Helixer/HelixerPost 是最明确把经典 gene finding 规则写进后处理的现代模型之一**。补充材料写得很具体：它先在每条链上做两遍扫描，先找 probable genic regions，再用一个 73-state HMM 生成 primary gene models；状态包含 intergenic、5′UTR、CDS、3′UTR、start/stop codon、coding phase；还显式建模 GT-AG、GC-AG、AC-AT 三类 intron 形式，并设最小 intron 长度；HMM transition penalties 是基于生物学知识手工设定的。换言之，Helixer 并不是“单纯黑盒深度学习”，而是**经典 gene grammar + 现代概率输出**的混合体。对边界过滤、假阳性抑制和 phase consistency，它依然很有研究启发。citeturn20view0

**Tiberius** 的核心价值则在于把结构化约束更深地并入模型本体。官方仓库把它定义为“deep learning + differentiable HMM layer”的 end-to-end gene predictor，并提供多 clade 权重、softmasked / unmasked / ClaMSA 等不同模式。这说明它不是把 gene grammar 完全放到后处理去做，而是更接近“**让模型在训练时就感知结构一致性**”。对你来说，这条路线的重要含义是：如果你的新方法只是在 backbone 上替换 CNN/Transformer，却没有回答“结构一致性在训练中如何进入”，那么要超过 Tiberius 会比较难。citeturn41view0turn42view0

**ANNEVO** 介于两者之间。它在神经网络端显式学习三类输出：category prediction、phase prediction、state prediction；然后在 gene structure decoding 阶段用 Viterbi 把这三类概率整合成符合生物学规则的 gene structures。作者还强调因为模型已经学到了 evolutionary/contextual learning，解码所需人工调参很少。也就是说，ANNEVO 的策略不是完全抛弃规则，而是把更多规则前移成可学习的 phase/state 表示。这个思路对你的课题很重要，因为它提示你：**可以把“phase、边界、状态迁移”从 decoder 里的硬编码，逐步改造成神经网络和 decoder 共享的中间变量。**citeturn26view0turn44view1

**SegmentNT** 在这方面则明显更“神经分割化”。它把任务写成 14 类或 7 类元素的多标签 per-base segmentation，使用阈值 0.5 做逐碱基预测，并在 human gene benchmark 上与 AUGUSTUS 比 F1/MCC/SOV；但它的原生工作流程没有像 Helixer/Tiberius/ANNEVO 那样公开强调 start/stop 一致性、ORF 完整性、相位约束和最终 gene-model 解码。因此，SegmentNT 非常适合做**高质量局部标签器或 backbone**，但如果你想把它变成真正的 ab initio protein-coding gene caller，几乎必然还要加一层结构化 decoding / post-processing。citeturn30view0turn33view0

**GENERanno preview** 在“生物学约束怎么合入”这件事上，当前公开信息最少。模型卡和仓库都强调它是 eukaryotic CDS annotator，并提供 end-to-end CLI，但没有像 Helixer 补充材料或 ANNEVO 方法部分那样，把 start/stop、phase、splice motif、strand-aware decoding 的机制完整公开出来。因此现阶段你可以把它当作“强潜在 FM annotator”，但还不能把它当成“约束集成策略的参考实现”。citeturn35view0turn35view1

从**复现角度**看，我建议这样排优先级。**最值得先复现的是 Tiberius**：代码、Docker/Singularity、模型配置、当前 clade 权重入口都公开，工程上有一定复杂度，但路径很清晰；而且它和你的目标任务最对齐。**第二是 Helixer**：代码、Zenodo 权重、自动下载脚本、web tool 都成熟，但因为模型/后处理/依赖链较长，实际 setup 难度略高，且速度明显慢于当前 ANNEVO 和 Tiberius。**第三是 ANNEVO**：它现在很可能是最强公开候选之一，代码与权重都在仓库中，但环境管理、CPU 解码资源、非商业许可证和“论文版 vs 当前版”的版本漂移都提高了复现成本；如果你的目标是发高水平方法学稿，仍然非常值得复现。**SegmentNT** 的复现门槛在“能跑”层面不高，因为 Hugging Face 权重和 notebook 都公开；但若要把它纳入可比 gene-caller benchmark，反而需要你自己补大量 protocol glue code。**GENERanno preview** 目前是最不建议作为“首个复现对象”的：虽然模型卡、代码和 CLI 已公开，但缺少足够透明的真核 benchmark 说明，容易把大量时间花在“猜作者 protocol”上。citeturn41view0turn17view0turn16search2turn9view1turn11view3turn29view1turn35view0turn35view1

最后，明确列出本次检索中的**公开信息缺口**。第一，Tiberius 当前各 clade 的**逐模型直链权重 URL**在仓库里是通过 YAML `weights_url` 解析的；我本次核验了这一机制与入口目录，但没有逐个拉取全部 YAML，因此如果你要做精确复现，仍需在本地把冻结 config 保存下来。第二，ANNEVO 的论文版和 GitHub v2.3.x 版已有明显性能和速度漂移，例如 human BUSCO 从 95.7 到 98.3、human runtime 从 82 分钟到 31 分钟再到 19 分钟，比较时必须显式冻结版本。第三，GENERanno 的真核 preview 仍缺少足够公开的数值表、split 描述与 metric scripts，这会直接影响你是否能把它放进主 benchmark。citeturn42view0turn9view2turn11view3turn35view1
