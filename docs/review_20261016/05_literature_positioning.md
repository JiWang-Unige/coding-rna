# 竞品与发表标准调研：GLM 微调 + 结构化解码做植物 ab initio 基因注释

调研日期 2026-10-16。范围：2023–2026 文献，重点是深度学习 ab initio 基因结构注释。先读过项目内 `docs/M28_complete_chain_competition_prior_work_20261003.md` 与 `docs/12_publication_strategy.md`，本文只补充其未覆盖的内容，不重复其中已核对的 Conrad、Tsuboi、BSO、GeneDecoder 等先例。

核查等级标注：
- **[原文]**：本次直接读到论文全文、扩展方法或官方模型卡。
- **[摘要]**：只读到摘要、预印本元数据或网页摘要，数字未在正文核对。
- **[二手]**：来自第三方转述，仅作线索。
- **[未核实]**：检索未能确认。

所有数字都是各作者在自己的评估面板和指标下报告的，不能互相拼接，也不能与本项目的 exact-chain F1 直接比较。

---

## 0. 一页结论（先看这里）

1. **本项目的骨干已经有"上游原作者的完整竞品"**。本项目使用的 `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` 的 config（hidden 2048、26 层、k=6、最大 16384 token）与 HuggingFaceBio 发布的 `Carbon-A-1.2B` 的 config 逐项一致。Carbon-A 技术报告写明它沿用 GENERanno 框架、骨干来自 GENERator-1.2B，并带转移网络、序列约束 Viterbi、基因置信度网络，发布了 GFF3 流程，许可 Apache-2.0。因此"GENERanno 微调 + 结构化解码"在审稿人眼里首先要和 Carbon-A 全流程比较。这是本次调研最重要的发现。依据：[GENERanno preview config](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview/raw/main/config.json)、[Carbon-A-1.2B config](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/main/config.json)、[Carbon-A 技术报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)。"二者为前驱/同源关系"是据 config 一致和技术报告措辞的推断，不是作者的明确声明。
2. **GLM + 序列语法解码的路线在 2025–2026 已经拥挤**：GeneCAD（PlantCAD2 + 染色体级 CRF）、PlantGeneAnn（Mamba-MoE + semi-Markov）、Carbon-A（GLM + 转移网络 + Viterbi + 置信度）、GeneDecoder（GPN/NT + 潜变量 CRF）、ANNEVO（MoE + Viterbi）、Tiberius（可微 HMM）。单凭"GLM + 解码"不构成新颖性。
3. **准确率门槛已被抬高**。同类工作被接受（Nature Methods）或预印本发布时，普遍是在 12–566 个物种上报告，并且对比 Helixer、Tiberius、ANNEVO 取胜或持平。本项目当前在拟南芥/水稻 DEV 上的完整链 F1 为 55.17%/58.90%，原生 ANNEVO 为 76.27%/74.43%（项目内记录）。在这个差距下，"新 SOTA 工具"这条叙事不成立。
4. **Setaria viridis 不能当作新颖的盲测卖点**。它本身就是 Helixer（HelixerPost 测试集）的测试物种，并被 PlantGeneAnn 保留在其泄漏控制后的测试集中；它同时是否进入 ANNEVO、Carbon-A 的训练集尚未核实（见 §3.4）。
5. **Nature Communications 对"准确率低于现有工具的工具论文"几乎不可能接受**。对"校准不确定性 / 低假阳性 / GLM 表征分析 / 负结果机制解释"这类角度，只有在做成跨多种 GLM、多种解码器、多物种且有明确可迁移结论的系统性分析，并且最好带来一个可验证的实用收益时才有可能；仅用单一骨干、两个开发物种、一个盲测物种做不到（§5）。该判断来自对同类论文证据体量的归纳，不是 Nature Communications 的官方标准（我未能读取其编辑标准原页，见 §6 的未核实项）。

---

## 1. 竞品逐项核对

### 1.1 总表

| 工具 | 发表 / 年份 | 骨干与解码 | 评估范围（作者报告） | 开源与易用性 |
|---|---|---|---|---|
| Helixer（NM 版） | Nature Methods 23(4):732–739，2025-11-24 在线，[DOI 10.1038/s41592-025-02939-1](https://doi.org/10.1038/s41592-025-02939-1) | CNN+BiLSTM 逐碱基分类 + HelixerPost HMM（神经网络与 HMM 分开训练） | 训练/评估基因组由 237 扩至 936；45 个测试物种（13 植物、11 脊椎、15 无脊椎、其余真菌）[原文] | GitHub 手动安装仅 Linux，推荐 Docker/Singularity，另有 Galaxy 与网页；需 Nvidia GPU（8–11 GB 可运行）；GPL-3.0 [原文：[Helixer README](https://github.com/weberlab-hhu/Helixer)] |
| Tiberius（2024） | Bioinformatics 40(12):btae685，[DOI 10.1093/bioinformatics/btae685](https://doi.org/10.1093/bioinformatics/btae685) | CNN+LSTM + 可微 HMM 层，端到端训练 | 37 个哺乳动物基因组训练与评估；测试 human 等 3 个基因组 [原文，[PMC11645249](https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/)] | 见下行 |
| Tiberius（多谱系） | bioRxiv 2026-04-28，[DOI 10.64898/2026.04.24.720536](https://doi.org/10.64898/2026.04.24.720536) | 同上，扩到 Mesangiospermae（被子植物）、Fungi、Vertebrata、Insecta、Chlorophyta、Bacillariophyta | 33 个测试物种（5 个谱系三工具同台）；gffcompare 精确匹配的 CDS-exon/转录本/基因 F1；不含 UTR；无 BUSCO 结果 [摘要+正文提取] | pip 安装 launcher + Singularity/Docker/源码/网页服务；TensorFlow；建议 GPU ≥8 GB；MIT [原文：[Tiberius README](https://github.com/Gaius-Augustus/Tiberius)] |
| ANNEVO | Nature Methods 23:740–748，2026-03-12，[DOI 10.1038/s41592-026-03036-7](https://doi.org/10.1038/s41592-026-03036-7) | 卷积分 bin + Transformer + MoE（8 专家）+ 分辨率恢复；Viterbi 解码；5 个谱系专用模型 | 566 个 RefSeq 测试物种（Embryophyta 与真菌全部候选物种入测，动物谱系各随机 50 个）；与 12 个模式物种上的 BRAKER3、GeneMark-ETP、Helixer、Tiberius 对比；793 个物种的 BUSCO 对照参考注释 [原文] | GitHub，conda 环境；Linux+CUDA；单张 RTX 4090 基准；人类基因组 GPU 峰值显存 3.8 GB；学术非商用许可 [原文 + [ANNEVO README](https://github.com/xjtu-omics/ANNEVO)] |
| Carbon-A | 技术报告 + 博客，2026-10-08 前后发布（未同行评议），[报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)、[博客](https://huggingface.co/blog/HuggingFaceBio/carbon-annotator-genbank-genome-annotation) | 1.2B 双向编码器，6-mer，GENERator 骨干改成 encoder；转移 TCN + Viterbi + 置信度 TCN | 42 个基因组（28 个训练快照内的参考基因组 + 14 个 2026 年后新增的时间留出基因组）；对比 AUGUSTUS、Helixer、Tiberius、ANNEVO、OrionGeno、SegmentNT、NTv3；Iso-Seq 湿实验验证；注释全部 GenBank 真核组装生成 CADB [原文] | HF 上单检查点 Apache-2.0；含 `run_annotation.sh` 流程，输出 GFF3 与逐基因置信度；需 CUDA、flash-attn 或 SDPA [原文：[模型卡](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B)] |
| GeneCAD | bioRxiv 2025-11-03，[DOI 10.1101/2025.10.31.685877](https://doi.org/10.1101/2025.10.31.685877) | PlantCAD2 表征 + 8 层 ModernBERT + 染色体级 CRF（Viterbi）+ ReelProtein 蛋白语言模型过滤；每位点输出一个规范转录本 | 5 个留出被子植物（Juglans regia、Coffea arabica、Phaseolus vulgaris、Nicotiana sylvestris、N. tabacum）；对比 BRAKER3、Helixer、Tiberius；转录本级 F1 平均比 Helixer 与 BRAKER3 高 8–10% [原文] | GitHub plantcad/genecad；CC BY-NC 预印本 |
| PlantGeneAnn | bioRxiv 2026-06-26，[DOI 10.64898/2026.06.25.733695](https://doi.org/10.64898/2026.06.25.733695) | PlantBiMoE（Mamba-MoE，116M）+ U-Net 解码头 + semi-Markov 解码 | 摘要称 20 个物种基准、对 Helixer、Tiberius、ANNEVO 在 exon/intron/locus 级平均 F1 为 0.856/0.880/0.707；我读到的 JATS 全文版本写 13 个留出物种、对 HelixerPost、ANNEVO、AUGUSTUS、NTv3，转录本级 F1 0.575 对 ANNEVO 0.569（两处不一致，可能是版本差异）[原文+摘要] | GitHub qzzhang0131/PlantGeneAnn；HF 权重 MIT；需 mamba-ssm [原文：[模型卡](https://huggingface.co/qzzhang/PlantGeneAnn-multi-species)] |
| OrionGeno | bioRxiv 2026（BGI），仅从聚合页和 PyPI 确认存在；作者称跨脊椎、无脊椎、绿色植物、真菌蛋白级 F1 相对最佳方法提升 37.2% | 多物种谱系感知深度学习框架，同时预测基因、UTR 与重复序列 | 一份 CNGB 数据报告了 1,249 个真核科的注释，BUSCO 中位 94.7% [二手：[聚合页](https://vitagri.org/evidence/IRmoq83umo-99eccd/)、[CNGB](https://db.cngb.org/data_resources/project/CNP0009228)] | PyPI `oriongeno` 1.0.2；AGPL-3.0 [二手] |
| GENERanno / GENERator | GENERator：[arXiv 2502.07272](https://arxiv.org/abs/2502.07272)（v5 2026-01-31）；GENERanno：[bioRxiv 2025.06.04.656517](https://www.biorxiv.org/content/10.1101/2025.06.04.656517v1)（只覆盖原核） | 真核 CDS annotator preview 是 token classification（NON_CODING/CDS）模型 | **GENERanno 仓库与模型卡不报告任何 CDS 注释指标**：README 明确没有给注释器的评估结果，模型卡只声称优于 AUGUSTUS、Helixer、SegmentNT 但没有数字，并写"More technical details are coming soon" [原文：[GitHub](https://github.com/GenerTeam/GENERanno)、[模型卡](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview)] | MIT；`pip install -r requirements.txt`，`cds_annotation.py`；预览版 |
| SegmentNT | Nature Methods 2025，[DOI 10.1038/s41592-025-02881-2](https://doi.org/10.1038/s41592-025-02881-2) | Nucleotide Transformer + 1D U-Net，14 类元件 | 剪接位点 MCC 0.75/0.76（SpliceAI 0.67/0.59）[摘要]；不提供完整 CDS 基因模型 | JAX 与 PyTorch 权重公开 |
| NTv3 | bioRxiv 2025-12，[DOI 10.64898/2025.12.22.695963](https://doi.org/10.64898/2025.12.22.695963) | U-Net 型，1 Mb 上下文，监督使用 24 个动植物物种的功能轨道与注释 | 摘要称注释达到 SOTA [摘要]；Carbon-A 报告称其不输出完整 CDS 基因模型且碱基/外显子 F1 明显更低 [原文] | 权重公开 |
| GeneDecoder | [arXiv 2505.03377](https://arxiv.org/abs/2505.03377)（项目文档另查到 2023 年 OpenReview 版本） | GPN/NT 表征 + 潜变量 CRF + Viterbi | 位置级 MCC/F1，片段输入，不是整染色体 [见项目既有记录] | 未核实 |
| PlantCaduceus / PlantCAD2 / AgroNT | PlantCaduceus：PNAS 2025，[DOI 10.1073/pnas.2421738122](https://doi.org/10.1073/pnas.2421738122)；AgroNT：Commun Biol 2024，[DOI 10.1038/s42003-024-06465-2](https://doi.org/10.1038/s42003-024-06465-2)；PlantCAD2：[bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.27.672609) | 植物 DNA 语言模型（PlantCaduceus 单核苷酸，16 个被子植物预训练） | 评估的是剪接位点、翻译起始、变异效应等位点级任务，不是完整基因结构 [摘要/二手] | 公开 |
| Evo 2 | Nature 2026，[DOI 10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5) | 40B 参数 StripedHyena 2 | 机制可解释性发现外显子-内含子边界等特征；仓库有基于嵌入的外显子分类 notebook；不是注释工具 [摘要/二手] | 完全开放 |
| BEND | ICLR 2024，[arXiv 2311.12570](https://arxiv.org/abs/2311.12570) | 基准，不是方法 | 人类基因发现任务（GENCODE 44，80% 序列同一性切分）：AUGUSTUS MCC 0.80，最好的 DNA LM（NT-MS）0.68，NT-v2 0.64，DNABERT-2 0.43 [原文] | GitHub |
| GENATATOR、Modern Gene Finders | ICLR 2025 研讨会 / ICLR 2026 | DNA LM 注释，后者聚焦哺乳动物转录本边界（ModernGENA） | 摘要级 [二手：[GENATATOR](https://iclr.cc/virtual/2025/36979)、[Modern Gene Finders](https://iclr.cc/virtual/2026/10013741)] | 未核实 |
| 经典基线 | BRAKER3：Genome Research 34:769，2024，[DOI 10.1101/gr.278090.123](https://doi.org/10.1101/gr.278090.123)；GeneMark-ES/ETP；AUGUSTUS：[Stanke 2003](https://gobics.de/mario/papers/GenePred2003.pdf) | 需 RNA-seq 与蛋白证据（BRAKER3） | BRAKER3：11 个物种，转录本 F1 比 BRAKER1/2 高约 20 个百分点；提供 Docker/Singularity [摘要] | 开源 |

### 1.2 与你们路线最相关的几项细节

**Helixer（NM）**[原文，[PMC13076211](https://pmc.ncbi.nlm.nih.gov/articles/PMC13076211/) 经 Europe PMC 全文 XML 读取]：
- 指标体系：逐碱基 phase F1、subgenic F1、genic F1；特征级 exon、intron、intron chain、transcript、gene 的 P/R/F1；BUSCO；植物的 orthogroup 占据度和 MapMan4 蛋白家族 P/R；拟南芥 Araport11 逐位点比对（约 70% 精确匹配）；消融（转移权重、LSTM 与混合架构）；运行时间。
- 植物数字（作者平均值，每组按物种平均）：exon F1 Helixer 0.714、GeneMark-ES 0.181、AUGUSTUS 0.504（AUGUSTUS 仅 4/13 个植物物种有预训练模型）；intron chain 0.434；transcript 0.462。植物 phase F1 为 0.810。
- 与 Tiberius 的比较只在哺乳动物：Tiberius 仍更好，基因召回与精度约高 20 个百分点。作者承认并以"提供更广谱系"为卖点。
- 作者自己承认的局限：消融没有随机种子重复；比较工具若按其他场景（含重复屏蔽与同源蛋白）配置本会更高。这是审稿人对"基线是否调优"的常见质疑点。

**ANNEVO（NM）**[原文，作者放在 xjtu 个人主页的 PDF 全文]：
- 训练与测试切分：每个主要谱系按亚谱系选 BUSCO 分数最高的物种，随机在较高分类阶元上划分训练/验证（基因组大小比约 5:1），未入训练的物种为测试候选；Embryophyta 与 Fungi 的全部候选物种进入测试。
- 指标：NT(CDS)-F1、基因级 F1、BUSCO；明确声明只注释最长转录本，评估目标取参考的最长转录本；承认放宽到任意 isoform 会提高成绩。
- 假阳性分析：60.6% 来自碎片化 mRNA 预测，18.3% 与已注释非编码转录本的外显子重叠。
- 资源：GPU 5× 于 Tiberius、24× 于 Helixer、83× 于 BRAKER3；CPU-only 仍快 10× 以上；人类基因组 GPU 内存 3 GB（Helixer 8.6 GB，Tiberius 超过 32 GB）。
- 额外证据：在 793 个物种中 318 个的 BUSCO 高于参考注释；在甘蓝（Brassica oleracea）上纠正 55 个有 RNA-seq 独立支持的 BUSCO 基因。
- 方法上与本项目最相关的事实：上下文每个核心区 30,720 bp、两侧各 5,120 bp；错误注释掩码；focal 损失；与项目记录一致。
- README 报告的 12 个模式物种平均：ANNEVO exon 召回/精度 91.4/90.2，locus 召回/精度 76.3/74.3，BUSCO 97.8；Tiberius 89.8/88.7、74.0/68.8、96.3；Helixer 86.1/75.3、50.2/47.0、92.5 [原文：[ANNEVO README](https://github.com/xjtu-omics/ANNEVO)]，是作者自己的汇总，且基于 2026-06-01 的最新版本。

**Tiberius（多谱系预印本）**[摘要+正文提取]：
- 在 33 个物种上，Tiberius 的基因级 F1 比 Helixer 高 12–37 个百分点，比 ANNEVO 高 10–22 个百分点；与 BRAKER3 相比，Mesangiospermae、Fungi、Bacillariophyta、Chlorophyta 上基因级 F1 差距小于 3 个百分点。
- 平均运行时间（72 线程 + GPU）：Tiberius 26 min、ANNEVO 30 min、Helixer 178 min、BRAKER3 2,170 min。
- **作者各自的比较互相矛盾**：ANNEVO 论文称其全面优于 Tiberius，Tiberius 预印本称 Tiberius 比 ANNEVO 高 10–22 个百分点；第三方的 Camelina sativa 单基因组幻灯片称两者相当、Tiberius 略优（[JOBIM 2026](https://pepi-ibis.inrae.fr/sites/pepi-ibis/files/PoleAnnotationGenome/JOBIM2026/05_JOBIM2026_Legeai-Chathuant.pdf)，我未能解出正文，以 [二手] 对待）。这意味着一个公开、冻结、可复现的独立比较面板本身有价值，也意味着你们用冻结版本并在同一面板上逐物种比较的做法是对的。

**Carbon-A（技术报告）**[原文，已下载并解析 PDF]：
- 汇总：42 个基因组上整体核苷酸、外显子、基因级 F1 最高；植物平均基因级 F1 为 0.776。
- 基因级 F1（Table S3，百分比，各工具的 P|R|F1 中取 F1）：
  - 拟南芥：AUGUSTUS 51.3、Helixer 74.5、Tiberius 83.0、ANNEVO 84.3、OrionGeno 83.0、Carbon-A 82.8（过滤后 85.4）。
  - 水稻：AUGUSTUS 30.6、Helixer 58.6、Tiberius 70.6、ANNEVO 78.1、OrionGeno 60.7、Carbon-A 79.9（过滤后 81.3）。
  - 拟南芥与水稻都属于 Carbon-A 的"训练快照内"参考基因组（Table S5 列出 GCF_000001735.4 与 GCF_034140825.1）。
- **置信度**：逐基因置信度预测"是否与参考 CDS 精确匹配"，AUROC 0.876；阈值 0.10 时过滤掉 1.3% 的精确匹配与 33.4% 的其他预测，过滤集的基因级精度 0.792、召回 0.691、F1 0.736。
- OOD：只用动物训练的检查点（Carbon-A-animal）在植物上核苷酸 F1 0.873；非标准遗传密码纤毛虫 Tetrahymena 上核苷酸 F1 0.960，而其他方法接近 0。
- 效率：人类 1 号染色体 2,027 s（NTv3 1,894 s；Helixer 6,030 s；AUGUSTUS 单线程 7,510 s）；Tiberius、ANNEVO、OrionGeno 比 Carbon-A 快约 3–4 倍。
- 湿实验：猫、仓鼠、鸡、拟南芥的 Iso-Seq；完整 CDS 的 Iso-Seq 支持比例约 0.62，与参考相当。
- 项目既有文档（10-09）写"完整 Carbon-A 调用资产尚未定位齐"。当前模型卡声称含 `annotation_pipeline/run_annotation.sh`，并说明输出 GFF3、CDS/蛋白 FASTA、逐基因置信度。我只读到模型卡摘要，**未验证仓库文件实际可运行**，建议重新核查再决定是否把 Carbon-A 作为固定基线。

---

## 2. "被接受档次"的证据标准归纳

### 2.1 同类论文的证据体量

| 论文（期刊/平台） | 评估物种数与跨度 | 盲测方式 | 指标 | 基线 | 消融 | 资源报告 | 软件可用性 |
|---|---|---|---|---|---|---|---|
| Tiberius 2024（Bioinformatics） | 37 个哺乳动物，3 个测试基因组 | 物种级留出 | 基因/转录本/外显子 P/R/F1 | Helixer、AUGUSTUS、GeneMark、BRAKER3（带证据） | 有限 | 人类基因组 <2 h | GitHub，后有 Docker/Singularity/网页 |
| Helixer（Nature Methods 2025） | 936 个基因组用于训练评估，45 个测试物种，4 大类 | 测试物种预留，仅用于最终评估 | phase/subgenic/genic F1、特征级 F1、BUSCO、orthogroup、MapMan、Araport11 逐位点比对 | GeneMark-ES、AUGUSTUS、Tiberius（仅哺乳动物） | 有（转移权重、架构），无种子重复 | 单线程速度，27 min（263 Mbp）与 8.5 h（人类） | GitHub、Zenodo、Docker/Singularity、Galaxy、网页 |
| ANNEVO（Nature Methods 2026） | 566 个 RefSeq 测试物种 + 793 个物种的 BUSCO 对照；5 个谱系 | 训练/验证/测试按物种划分，在较高分类阶元随机分割；动物谱系各抽 50 个 | NT(CDS)-F1、基因 F1、BUSCO | Augustus（用近缘物种逐物种训练）、Augustus-Evidence、GeneMark-ETP、BRAKER3、Helixer、Tiberius | 有（见补充方法） | GPU/CPU 时间、峰值显存、规模 | GitHub（非商用）、conda |
| Tiberius 多谱系（bioRxiv） | 33 个测试物种，6 个谱系 | 谱系专用模型，测试物种留出 | gffcompare 精确匹配 exon/转录本/基因 | Helixer、ANNEVO、BRAKER3 | 实现层变更对精度影响 | 全面运行时间 | 网页、Docker、pip |
| GeneCAD（bioRxiv） | 5 个留出被子植物 | 物种留出，训练物种不在测试集 | 核苷酸/外显子/转录本 | BRAKER3、Helixer、Tiberius | 有（CRF、过滤等） | 有限 | GitHub |
| PlantGeneAnn（bioRxiv） | 13（全文）或 20（摘要）个留出植物，覆盖蔷薇类、菊类、单子叶 | 在 Helixer 测试集基础上，**重新审计预训练/训练泄漏并替换 5 个重叠物种** | 碱基、外显子、内含子、intron-chain、转录本；配对显著性检验 | HelixerPost、ANNEVO、AUGUSTUS、NTv3 | 有（预训练先验与 U-Net 头） | 峰值内存、微调成本 | GitHub + HF 权重 |
| Carbon-A（技术报告） | 42 个基因组，6 个真核类群 | **14 个 2026 年后新增的时间留出基因组**；声明留出基因组不在所有被评模型的训练数据中 | 核苷酸/外显子/基因 P/R/F1、蛋白级 BUSCO、置信度 AUROC | AUGUSTUS、Helixer、Tiberius、ANNEVO、OrionGeno、SegmentNT、NTv3 | 上下文长度等 | 单检查点吞吐 | HF 开放 |

### 2.2 可归纳的"最低门槛"与"强证据"（我的综合判断）

下列是基于以上论文的经验性归纳，不是任何期刊的官方清单。

| 维度 | 近年顶刊/高水平工作的常见做法 | 本项目现状 |
|---|---|---|
| 物种数量与系统发育跨度 | 至少 12 个物种且跨多个谱系（ANNEVO 的头对头基准），或谱系内 13–45 个测试物种（PlantGeneAnn、Helixer、Tiberius 多谱系） | 两个开发物种 + 1 个盲测物种，仅限禾本科/十字花科 |
| 盲测物种独立性 | 物种级留出；更好的做法是做上游预训练/训练泄漏审计并替换重叠物种（PlantGeneAnn）、或用时间留出基因组（Carbon-A） | GENERanno 预训练暴露 unknown（项目既有记录）；Setaria 暴露状态待核 |
| 指标体系 | 核苷酸 + 外显子 + 内含子链 + 转录本/基因（精确匹配）+ BUSCO；加分项为蛋白家族/orthogroup、Iso-Seq 或 RNA-seq 独立验证 | 项目已有 exact-chain、CDS-base、背景负担，缺 BUSCO 和独立转录证据 |
| 基线 | Helixer、Tiberius、ANNEVO 必须有；BRAKER3/GeneMark-ETP 作为带证据的上界；AUGUSTUS 作为经典基线；GLM 类的 Carbon-A、GeneCAD、PlantGeneAnn 越来越成为必需对照 | 已有三者；缺 Carbon-A、PlantGeneAnn、GeneCAD、BRAKER3 |
| 消融 | 对解码层、预训练先验、损失、转移权重做消融（Helixer、PlantGeneAnn、GeneCAD） | 项目的消融记录大多是已关闭的负结果，尚无成功模型的消融 |
| 资源与速度 | 报告 GPU/CPU 时间、峰值显存、与基线同硬件 | 部署指标空缺（项目文档已列为待办） |
| 跨物种泛化 | 谱系专用模型 + 留出物种，或单模型跨谱系（Carbon-A） | 未满足 |
| 软件可用性 | GitHub + Docker/Singularity 或 pip/conda + 预训练权重 + Zenodo 存档；GPU 需求写明 | 未发布 |

要点：
- 评估广度近两年快速抬高。2024 年 Tiberius 只用 3 个哺乳动物测试基因组就能发 Bioinformatics，而 2026 年的多谱系版本已经是 33 个物种。
- 竞争对手对"泄漏"越来越敏感：PlantGeneAnn 明确替换了与预训练/训练集重叠的 5 个物种。你们的盲测设计需要达到同等严格度，且需要公开说明 GENERanno 预训练暴露未知这一事实。
- 作者各自报告的基线数字互相矛盾（见 §1.2）。审稿人会特别关注基线是否使用冻结版本、官方配置、同一面板。你们已有的"固定版本与固定资产"做法符合这一点。

---

## 3. 空白与风险：针对"GLM 微调 + 结构化语法解码"

### 3.1 已被占据的点（做主贡献有高风险）

| 点 | 占据者与依据 |
|---|---|
| GLM 表征 + 序列约束 Viterbi/CRF 完整注释 | Carbon-A（转移网络 + 序列约束 Viterbi，含剪接类别先验）；GeneCAD（PlantCAD2 + 染色体级 CRF + 蛋白语言模型过滤）；PlantGeneAnn（semi-Markov）；GeneDecoder（GPN/NT + 潜变量 CRF）[原文] |
| 可微/端到端结构层训练 | Tiberius（可微 HMM，标准 HMM 的 24 个参数不训练）[项目既有记录 + 原文] |
| 边界加权、相位辅助、MoE、错误注释掩码、focal/Dice 损失 | Helixer（转移权重）、ANNEVO（掩码、focal、Dice、MoE）[原文] |
| 基因级置信度 | **Carbon-A**：置信度 TCN，AUROC 0.876，阈值 0.1 过滤后精度/召回数字齐全 [原文] |
| 训练标签噪声治理 | GeneCAD（MMLR 序列分数 + 蛋白语言模型过滤）、ANNEVO（错误区域掩码）[原文] |
| 湿实验独立验证 | Carbon-A（四个物种的 Iso-Seq，含拟南芥）[原文]；ANNEVO（对甘蓝 55 个 BUSCO 基因的 RNA-seq 验证）[原文] |
| 纠正参考注释、补全注释 | ANNEVO（793 个物种 BUSCO 对照）、Carbon-A（CADB）、Helixer（拟南芥缺失基因）[原文] |
| 长基因、非标准密码、动物到植物迁移 | ANNEVO（621 kb GMDS 位点）、Carbon-A（Tetrahymena、动物到植物迁移）[原文] |
| 剪接位点级变异效应零样本分析 | PlantGeneAnn（Bt2、OsMADS1、8,000 个位点的饱和突变）[原文] |
| 全基因组大规模注释资源库 | Carbon-A（CADB，48,167 个组装）、Tiberius（2,948 个脊椎动物组装）、OrionGeno（1,249 个真核科）[原文/二手] |
| 完整 CDS 链竞争训练、有限集合目标、搜索保留损失 | 项目既有记录已核对：mGene、BSO/LaSO、Conrad、Tsuboi 等先例 |

### 3.2 可能的空白（需要自己再做一轮穷尽检索才能宣称，我这次不是穷尽性检索）

1. **校准的基因/外显子后验，而不是学习出来的置信分**。
   - 现状：Carbon-A 的置信度来自单独的神经网络、以"精确匹配参考"为监督、报告 AUROC；Tiberius、Helixer、ANNEVO 没有逐基因后验。我这次检索没有找到报告 ECE、Brier、可靠性曲线或跨物种校准漂移的基因注释工作，但这是"检索中未见"，不是"不存在"。
   - 你们的全语法部分观测似然/配分函数理论上能给出结构的边缘后验。若做成"全 DNA 语法路径配分函数下的边缘后验，经跨物种验证校准，并且给出可选择性注释的 FDR 控制"，与 Carbon-A 的区别是有理论基础的后验而非判别式置信分。
   - 前提：项目既有文档已指出，当前 FP64 角色 gain 和零操作点"不是校准后的 gene 概率"，需要先证明这一点能做到。
2. **在匹配监督下，GLM 预训练、模型规模与解码器各自贡献多少的受控分解**。Carbon-A 报告写明 Tiberius、ANNEVO、OrionGeno 参数量分别为 8M、15M、212M，比 1.2B 的 Carbon-A 小 1–2 个数量级；PlantGeneAnn 的消融声称预训练先验与 U-Net 头都重要。我没有看到在相同标签、相同上下文、相同解码器下对"从头训练 vs GLM 初始化 vs 不同 GLM"做系统比较的植物工作。这是可行的分析型贡献，但不是工具型。
3. **低背景假阳性和注释噪声感知的评估框架**。ANNEVO 的假阳性分析显示 60.6% 来自碎片化预测；Helixer 在植物上召回高于参考但精度更低。若能量化"参考背景"对 FPR 的贡献（项目已有 reference-background 敏感性工作），并提出对未知注释链的保护式评估，可能有方法学价值；但目前只在两个物种上。
4. **完整 CDS 链作为单位的结构级评估与互斥/共存建模**。多数工作评估的是最长转录本的基因级或转录本级，且 Carbon-A 与 ANNEVO 都明确是"单一最长转录本"或"合并同链异构体"。对"一个窗口内多条不重叠基因共存"的完整链竞争建模，在文献中的先例主要来自 mGene、Conrad 一类早期工作。

### 3.3 若完整基因结构精度仍低于 Tiberius/ANNEVO

可发表角度的现实评估（我的判断，依据 §2 的证据体量和 §1 的竞争状况）：

| 角度 | 可发表性 | 需要什么才可能达到 Nature Communications | 主要风险 |
|---|---|---|---|
| 校准的基因后验 + 选择性注释（给定 FDR 的高置信子集） | 中。与 Carbon-A 置信度直接重叠，需要明显更强的校准证据 | 多物种（至少多个谱系）跨物种校准；与 Carbon-A 置信度、Tiberius/ANNEVO 概率的同台比较；说明高置信子集的精度/召回超过基线工具在同召回下的精度 | 如果高置信子集仍不及 ANNEVO 在同召回下的精度，就没有实用价值 |
| 低 FPR 工作点 | 低到中。M19 之后项目已以 FPR 为卖点，但完整链成绩落后会使"低 FPR"被解读为"少预测" | 在相同召回下比较 FPR，并有 BUSCO 与 Iso-Seq/RNA-seq 的支持 | Tiberius 被报告为"保守预测、精度更高"（ANNEVO 论文：基因精度高 2.3%、召回低 4.5%），低 FPR 本身不新 |
| GLM 表征能力的系统性分析（探针、层级、上下文长度、tokenization、预训练 vs 随机初始化） | 中。属于分析型论文，近期有 [Nature Communications 2025 的 DNA 基础模型基准论文](https://doi.org/10.1038/s41467-025-65823-8)（Feng 等，我只读到摘要级转述）先例 | 需要多个 GLM（至少 NT-v2、GENERanno/Carbon-A、PlantCAD2/PlantCaduceus、Evo 2）、多个解码器、多个谱系，并给出可迁移的规律与实用建议 | 只用单骨干、两个植物物种会被认为证据面过窄；Carbon-A 与 PlantGeneAnn 已经做了部分消融 |
| 负结果的机制解释（为什么完整链低于逐碱基/局部指标） | 低。单独作为主论文不太可能；更适合作为补充或 Bioinformatics/Genome Biology 级别的分析 | 需要从机制上解释并在另一个独立面板上验证"解释可预测改进"，且改进后的模型至少追平基线 | 没有改进验证的解释被看作事后叙述 |
| 可解释性（哪些层/头编码剪接、相位） | 低到中。Evo 2 已在 Nature 上有 SAE 可解释性发现 [摘要/二手] | 需要与功能性验证联动（例如定向扰动改变预测，PlantGeneAnn 已做零样本变异实验） | 对"注释工具"论文是附属证据 |

**关于 Nature Communications 会不会接受**：
- 我的判断：若主张是"新的注释工具"，且完整链精度低于 Tiberius/ANNEVO，也低于 Carbon-A，则接受概率很低，因为同领域的方法论文头部落点是 Nature Methods（Helixer、ANNEVO、SegmentNT），且这些都是在广谱物种上取得领先或接近证据驱动流程；若主张是"系统性分析或实用的校准机制"，则有可能，但需要足够广的证据面并能给出可被他人复用的结论。
- 局限：我这次检索没有找到 2025–2026 年发表在 Nature Communications 的 ab initio 基因结构注释方法论文（见 §6），无法用该刊近期同类论文校准其门槛；Nature Communications 的官方编辑标准页面被重定向拦截，没有读到原文。第三方咨询博客的说法是"需要超出增量的进展，但不要求跨学科影响"，且这些博客可信度有限（[Manusights 概述](https://manusights.com/blog/is-my-paper-ready-for-nature-communications)）。
- 建议：用一次轻量的预投稿咨询（Nature Communications 接受 presubmission enquiry）验证定位，再决定投入方向。该刊是否接受 presubmission enquiry 我没有核实。

### 3.4 盲测物种 Setaria viridis 的具体风险

- **Setaria viridis 是 Helixer 的测试物种**。PlantGeneAnn 全文写明保留了来自 HelixerPost 测试集的 8 个物种，其中含 Setaria viridis（见 [PlantGeneAnn 预印本](https://doi.org/10.64898/2026.06.25.733695)，我读的是 JATS 全文）。因此 Helixer 公开的比较、以及后来使用同一测试集的工作，都已经覆盖这个物种，"盲测 Setaria"在文献里不是新场景。
- 它的参考注释来自 JGI 未发表流程：Ensembl Plants 页面显示 v2.0（A10.1）有 38,334 个基因模型和 14,125 个可变转录本，组装 395.1 Mb [二手：[Ensembl Plants](https://plants.ensembl.org/Setaria_viridis/Info/Annotation)]。参考本身带噪声会限制"精确链 F1"的上限，需要用 BUSCO 或转录组独立验证补强。
- **暴露状态待核实**：ANNEVO（补充表 23–27 列出各物种的训练/验证/测试用途）、Carbon-A（Table S5）、Tiberius（被子植物模型的训练物种）、PlantGeneAnn 是否把 Setaria viridis 放进训练，我这次没有核对。这决定它对各基线是否"盲"，直接影响公平比较。需要逐一核对后再冻结盲测面板。
- 建议：把盲测扩展到多个不与任何被比工具训练集重叠的植物（参考 PlantGeneAnn 的做法：为每个重叠物种替换为系统发育上最近的无暴露近缘种）。

---

## 4. 项目当前位置与外部基线的粗略对照（不可视为同尺比较）

| 来源 | 拟南芥 | 水稻 | 指标与面板 |
|---|---|---|---|
| 项目：M28 FULL Q（DEV，同尺） | 55.17 | 58.90 | 完整链 exact-chain F1，两条 DEV contig，已用于开发 |
| 项目：原生 ANNEVO（DEV，同尺） | 76.27 | 74.43 | 同上 |
| 项目：原生 ANNEVO（TRAIN 同尺） | 76.83 | 71.18 | 完整链 exact-chain F1，六条 TRAIN 染色体 |
| 项目：原生 Tiberius / Helixer（TRAIN） | 75.42 / 68.12 | 63.58 / 50.29 | 同上 |
| Carbon-A 报告：基因级 F1 | Carbon-A 82.8；ANNEVO 84.3；Tiberius 83.0；Helixer 74.5 | Carbon-A 79.9；ANNEVO 78.1；Tiberius 70.6；Helixer 58.6 | 全基因组，精确 CDS 匹配；拟南芥与水稻均在 Carbon-A 训练快照内 |

数字上看，项目内部的 ANNEVO/Tiberius/Helixer 完整链 F1 与 Carbon-A 报告中对同样三个工具的基因级 F1 在同一量级，排序也一致（ANNEVO ≳ Tiberius > Helixer，水稻上差距更大）。这说明你们的同尺口径没有明显偏离外部报告，但两边的面板、参考集、精确匹配定义都不同，不能合并成一个榜单。

---

## 5. 对路线的含义（给决策用）

1. **先把 Carbon-A 补成固定基线**，并核对其流程在项目环境能否复跑；否则审稿人一定会问"为什么不和上游作者的完整流程比较"。若暂时不能复跑，只能把它写作已发表的方法先例，不能填作者数字入榜。
2. **不要把"GLM + 语法解码"当成主贡献**。它已被 Carbon-A、GeneCAD、PlantGeneAnn、GeneDecoder 占据。
3. **若走校准后验路线，最小可发表证据集**：
   - 校准度量（ECE、Brier、可靠性曲线）按物种分层并跨物种漂移；
   - 与 Carbon-A 置信度和其他工具的同台比较；
   - 选择性注释的实用收益：在相同召回下比基线工具精度更高，或在相同精度下召回更高；
   - 至少 10 个以上植物物种、覆盖多个谱系，并做泄漏审计；
   - 独立转录证据（Iso-Seq、RNA-seq、蛋白质组）与 BUSCO。
4. **若走分析型路线**，需要把单骨干扩为多 GLM 对照，并加受控的"预训练 × 规模 × 解码器"设计，否则证据面不足。
5. **若两条都不成立**，更现实的目标是 Bioinformatics、Genome Biology、The Plant Journal 或 NAR 一类的方法/分析论文。该判断是我的推断，不是对这些期刊门槛的核实。

---

## 6. 未核实项与方法局限

- Nature Communications 官方编辑标准页面被重定向拦截，没有读到原文；§3.3 里的刊物接受度判断是基于同类论文的体量与第三方博客的推断。
- 没有找到 2025–2026 年在 Nature Communications 发表的 ab initio 基因结构注释方法论文，可能是检索工具局限，不能据此断言该刊不发表此类论文。
- Helixer NM 与 ANNEVO 的补充材料（逐物种表、物种使用表）没有读到，所以各工具的训练物种清单、Setaria 暴露状态未核对。
- ANNEVO 全文 PDF 来自作者主页托管文件，页面标注为 Nature Methods 排版版本；我解析了文本，但图中的逐物种数字没有提取。
- OrionGeno 只通过聚合页、CNGB 和 PyPI 确认存在，bioRxiv 的 DOI 与全文没有核对。
- PlantGeneAnn 的摘要与我读到的全文对基准规模（20 个与 13 个物种）、基线列表（含/不含 Tiberius）描述不一致，需要以最终版本为准。
- Carbon-A 为技术报告和博客，未经同行评议；其评估由作者自己执行，基线配置在其 Appendix Table S5。
- Carbon-A 与 GENERanno 预览版"同源"是根据 config 完全一致和技术报告措辞做出的推断，没有作者的明确声明。
- GeneDecoder 的数字、BEND 之外的 gene-finding 基准、GENATATOR 与 Modern Gene Finders 只读到摘要。
- 本次没有向任何检索服务发送项目数据、权重或未公开指标。

## 7. 主要来源清单

- Helixer，Nature Methods 2025：[DOI](https://doi.org/10.1038/s41592-025-02939-1)，全文经 [Europe PMC/PMC13076211](https://pmc.ncbi.nlm.nih.gov/articles/PMC13076211/) 读取；预印本 [bioRxiv 2023.02.06.527280](https://www.biorxiv.org/content/10.1101/2023.02.06.527280v2)。
- ANNEVO，Nature Methods 2026：[DOI](https://doi.org/10.1038/s41592-026-03036-7)；[README](https://github.com/xjtu-omics/ANNEVO)；[Behind the Paper](https://communities.springernature.com/posts/closing-the-annotation-gap-with-annevo)。
- Tiberius：[Bioinformatics 2024](https://doi.org/10.1093/bioinformatics/btae685)；[多谱系 bioRxiv](https://doi.org/10.64898/2026.04.24.720536)；[GitHub](https://github.com/Gaius-Augustus/Tiberius)。
- Carbon-A：[技术报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)；[博客](https://huggingface.co/blog/HuggingFaceBio/carbon-annotator-genbank-genome-annotation)；[模型卡](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B)。
- GENERanno / GENERator：[GitHub](https://github.com/GenerTeam/GENERanno)；[模型卡](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview)；[GENERanno bioRxiv（原核）](https://www.biorxiv.org/content/10.1101/2025.06.04.656517v1)；[GENERator arXiv](https://arxiv.org/abs/2502.07272)。
- GeneCAD：[bioRxiv](https://doi.org/10.1101/2025.10.31.685877)。PlantGeneAnn：[bioRxiv](https://doi.org/10.64898/2026.06.25.733695)、[模型卡](https://huggingface.co/qzzhang/PlantGeneAnn-multi-species)。OrionGeno：[CNGB](https://db.cngb.org/data_resources/project/CNP0009228)。
- SegmentNT：[DOI](https://doi.org/10.1038/s41592-025-02881-2)。NTv3：[bioRxiv](https://doi.org/10.64898/2025.12.22.695963)。BEND：[arXiv](https://arxiv.org/abs/2311.12570)。Evo 2：[DOI](https://doi.org/10.1038/s41586-026-10176-5)。PlantCaduceus：[DOI](https://doi.org/10.1073/pnas.2421738122)。AgroNT：[DOI](https://doi.org/10.1038/s42003-024-06465-2)。
- BRAKER3：[Genome Research 2024](https://doi.org/10.1101/gr.278090.123)。AUGUSTUS：[Stanke & Waack 2003](https://gobics.de/mario/papers/GenePred2003.pdf)。
- DNA 基础模型基准，Nature Communications 2025：[DOI 10.1038/s41467-025-65823-8](https://doi.org/10.1038/s41467-025-65823-8)（摘要级转述）。
- Setaria viridis 参考注释：[Ensembl Plants](https://plants.ensembl.org/Setaria_viridis/Info/Annotation)。
