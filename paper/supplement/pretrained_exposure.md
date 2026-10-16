# M28 当前预训练模型的训练暴露

截至2026-10-13。本记录把当前运行资产与作者训练说明对齐，用于论文方法和比较表述。已执行的DEV比较仍是开发比较，不是新的独立测试：ANNEVO作者声明排除了拟南芥和水稻，Helixer当前权重明确包含这两个训练物种，Tiberius声明使用水稻。Carbon-A来源记录把当前GENERanno固定权重接到RefSeq注释监督训练说明；具体训练窗口、验证划分及同源隔离尚未核实。完整Carbon-A流程已补齐transition/confidence网络及调用代码，但两份辅助权重的训练暴露须与主模型分别记录。输入只有raw FASTA与“所有方法都未见过该物种”是两个不同条件。以下日期明确的历史追溯保留当时判断；分阶段暴露见2026-10-11条目，最新来源表述见2026-10-13补充。

## 当前资产与可支持的判断

| 模型 | 本轮使用的固定身份 | 训练暴露证据与边界 |
| --- | --- | --- |
| GENERanno encoder | GENERanno-eukaryote-1.2b-cds-annotator-preview，revision b0483c23b6b63787b61a6d3a204a9b517d6ba345 | 发布CDS注释checkpoint的backbone，不是纯无监督base。Carbon-A作者明确声明保存了这一revision的权重，技术报告披露RefSeq六类真核生物的CDS监督训练。拟南芥与水稻的项目assembly出现在报告的训练快照参考子集；实际优化器训练窗口、划分及同源隔离未核实，不支持清洁zero-shot。[固定来源记录](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/f6789ee79e6ef397de333577ffec1cb734496620/provenance.json)、[技术报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf) |
| Carbon-A完整原生流程 | HuggingFaceBio/Carbon-A-1.2B，revision f6789ee79e6ef397de333577ffec1cb734496620；主模型沿用上述b0483c23源权重，另含两份TCN权重 | 共享已注释监督训练的主模型，并非独立编码器。固定README说明transition/confidence网络均为已训练网络；此次读取的来源映射与权重头未给出两者各自的训练／验证assembly、窗口或同源隔离。不能把主模型的2055 assemblies清单直接当作辅助网络的训练清单。[固定流程说明](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/blob/f6789ee79e6ef397de333577ffec1cb734496620/annotation_pipeline/README.md)；分阶段边界见文末。 |
| ANNEVO | 源码42c920f70f52f96c7c06c665a639ec47d9304fff，Magnoliopsida权重 | 当前权重与v2.2.3的官方Git blob相同。该版本作者声明训练排除拟南芥、水稻等模式物种；这是作者披露，不是本项目已审计的完整assembly、验证集或同源隔离。[固定版本说明](https://github.com/xjtu-omics/ANNEVO/blob/f1ccd9c36bd6dce3e32b25dafe50c3f5809acd0e/docs/boundary_aware_model.md) |
| Helixer | 0.3.7/d17bb496b8842542590644b2437798711c0125d5，land_plant_v0.3_a_0080.h5 | 补充表S29把该权重对应到seed_009_017；S23该列含Athaliana与Osativa，S18对应TAIR10与v7.0。两物种均是作者列出的训练物种，不是物种未见。[补充表S18、S23、S29](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-025-02939-1/MediaObjects/41592_2025_2939_MOESM1_ESM.pdf) |
| Tiberius | 2.0.8隔离包871d4979568ae09e56de7eca5d7489046ee833ea，原angiosperms模型 | 固定配置训练列表含OsativaKitaakev3.1与Sviridisv4.1，未列拟南芥。水稻已知作者声明物种暴露；未列拟南芥不扩写成独立审计过的完整未见证明。[固定配置](https://github.com/Gaius-Augustus/Tiberius/blob/871d4979568ae09e56de7eca5d7489046ee833ea/model_cfg/angiosperms.yaml) |

模型资产、运行参数与原生执行收据见[native_baseline_preparation.md](../../benchmarks/native_baselines/reports/native_baseline_preparation.md)。表中作者训练物种与本项目具体DEV序列是否有坐标、样本或同源重叠尚未核对；物种训练暴露不自动等于泄漏，也不能由同物种不同assembly反推独立泛化。

## ANNEVO 新植物模型不能套用旧论文清单

作者[2026-04-23的v2.2.3发布说明](https://github.com/xjtu-omics/ANNEVO/releases/tag/v2.2.3)新增植物模型；其固定boundary-aware说明披露更新了数据预处理，且新策略尚未发表。后续[v2.3.0发布说明](https://github.com/xjtu-omics/ANNEVO/releases/tag/v2.3.0)说明植物继续沿用v2.2.3权重。当前固定源码的[更新记录](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/docs/update_history.md)另记v2.3.3于2026-08公开训练代码。

本轮直接比较官方Git树：v2.2.3的ANNEVO_model/ANNEVO_Magnoliopsida.pt与当前saved_model/ANNEVO_Magnoliopsida.pt均为blob 21487d877430abf8e552c509560fd1218f19305e。这个已有对象身份将当前实际权重接到作者的物种排除说明，不新增或重算大权重摘要。不能把原论文的旧训练清单无条件移植到更新权重，也不能把源码版本、权重版本和论文版本混成一个字段。

## GENERanno 和后续迁移表述

本节是2026-10-03的历史判断；上游来源披露已由文末2026-10-09的Carbon-A补充更新，尚未核实的坐标与同源隔离边界不变。

本项目feature_smoke.py与cache_full_train.py加载原AutoModelForTokenClassification，再取full.model并冻结。移除分类头、冻结encoder、当前forward不读取标签，都不抹去上游checkpoint可能继承的注释知识。确切revision模型卡没有提供训练或验证物种/assembly清单；既有[方法暴露说明](../../archive/m28_pre_full_grammar/docs/M28_implemented_methods_B2.md)继续有效，不使用较小eukaryote-0.5b-base或prokaryote论文的数据规模来替代本1.2B checkpoint的真实暴露。

因此，正文可描述“基于发布CDS注释checkpoint的冻结表示，进行raw-FASTA直接注释的开发比较”，不能描述为“所有模型均从未接触过目标物种”或“纯无监督预训练下的跨物种泛化”。后续独立面板必须分别标注每种模型的暴露资格；GENERanno未知仍需保留，不靠下游head未见标签替代。

Helixer补充表S18将Setaria viridis v2.0列为作者test；这与Tiberius训练列表中的Sviridisv4.1不同。这里只读公开元数据，不访问本项目封存的Setaria DNA、参考或结果，也不据此自动放行、选新物种或修改冻结实验。

## 本轮来源核对与执行影响

Exa进行了三个定向检索，每次请求5条结果，sources_reviewed按Skill计数为15，不代表15份独立来源均已核验。最终依据是已读取的固定模型卡、官方版本文档和Git对象、固定Tiberius配置及Helixer原始补充PDF。补充PDF的S23（印刷页81）与S29（印刷页85）已渲染检查列对应关系，避免把train_val池或邻列误读成实际训练集。

当前ANNEVO固定commit中旧boundary_aware_model.md已不存在；首次读取返回404，随后从原v2.2.3对应固定commit读取，不把失败当作内容缺失证明。v2.3.3 release页面未能读取，未用作证据，改用当前固定update_history.md。

本核对不改变运行中的FULL Q/C0、A6或六个原生基线，不重新训练、扫描、换阈值或打开新数据。原生新表与历史缓存表分开、共同7728 primary分母及各原gate不变。继续按原协议进行raw-FASTA部署比较；清洁zero-shot与因果消融仍需各自证据，不能用一张性能表同时替代。

## GENERanno 1.2B 来源的定向追溯补充

2026-10-03，本轮再作两次公开定向检索、每次5条返回额度，并核对实际 `feature_smoke.py`／native caller；固定身份仍为表中1.2B preview、revision b0483c23b6b63787b61a6d3a204a9b517d6ba345，没有下载新权重或改变运行配置。

作者仓库当前 main 由官方 commits API 接到29a93db37ec8b9a44e41502de5dbb4e5cf630830；其[固定README](https://github.com/GenerTeam/GENERanno/blob/29a93db37ec8b9a44e41502de5dbb4e5cf630830/README.md)模型表将本1.2B preview的Data记为1T，并记录2026-02-10发布。但该字段没有分开预训练与注释监督、唯一序列量与训练曝光，亦没有给出物种／assembly及划分清单；不能据此将1T解释成已审计的独立预训练数据集或补齐未知暴露。

当前[固定config](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview/raw/b0483c23b6b63787b61a6d3a204a9b517d6ba345/config.json)为GenerannoForTokenClassification、k=6、hidden_size=2048、CDS/NON_CODING两个类别、两个prediction heads。作者[2025预印本v3](https://www.biorxiv.org/content/10.1101/2025.06.04.656517v3)摘要和第2.1至2.3节描述的是500M原核模型、单碱基tokenizer及715B原核DNA。这些架构与对象不同，不能将其RefSeq来源、采样或零样本结论移植到本1.2B真核权重；也不能用仓库0.5B真核base的386B数据替代。

作者README链接的[cds-annotation数据卡](https://huggingface.co/datasets/GenerTeam/cds-annotation)描述原核评价数据，README另有真核fly运行示例；二者都不是本1.2B的训练manifest。这里只读公开文字元数据，没有下载例子、基因组或标签。固定模型卡仍未补齐训练／验证及同源暴露；本轮结论保持unknown，不推断“未见”或“已经泄漏”。

这不阻止当前raw-FASTA开发性能比较，但限制清洁zero-shot主张。以后可以分别检验下游监督物种外迁移和实际适用能力，并明确上游未知；“所有模型未见目标物种”另需完整资产来源证据。本次追溯不自动选择新面板、替换encoder、释放test/Setaria、调整A9或其他作业，也不发消息联系作者或公开研究结果。

## 固定资产的公开文件与训练入口

2026-10-06。GENERanno固定revision的[Hugging Face文件目录](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview/tree/b0483c23b6b63787b61a6d3a204a9b517d6ba345)完整返回11个文件，无后续分页：模型权重、README、模型／FSDP配置、模型实现、分词器及vocab，没有训练／验证物种或assembly manifest。只读取目录元数据，未下载权重、vocab或数据。已有作者固定仓库的目录及CDS调用器也不能补齐本1.2B模型的真实上游输入；暴露仍为unknown，不改写成未见或泄漏。

ANNEVO固定42c920f70f52f96c7c06c665a639ec47d9304fff的model_training目录有16个Git记录文件，未包含物种／assembly清单。[训练说明](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/model_training/README.md)和实际generate_h5_data.py要求调用者提供genome／annotation路径，model_train_seq.py要求提供训练／验证H5路径；示例species_a／b／c与通用train／val文件名不是作者实际训练资产身份。对该目录定向查询Vitis／riparia、Hibiscus／syriacus、Nicotiana／attenuata及GCF／GCA标识未见对应记录，不能据此判定未暴露。

更重要的是，公开训练入口model_train_seq→model_train→[model_construction_seq](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/model_training/src/utils.py#L150)构建ANNEVO_seq.ANNEVO；当前[原生推理工厂](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/src/utils.py#L7)在Magnoliopsida分支构建ANNEVO_plant.ANNEVO。因此，“已公开训练代码”不等于已经给出当前植物权重的完整复现合同；不能把这个sequence-only入口直接作为植物matched-retraining对照，也不能从它推定当前植物权重的物种、seed、更新量或实际训练配方。原作者植物权重与已关闭的原生比较保持不变。

训练说明还要求最长转录本并过滤premature stop及低于最小内含子长度的转录本；这是作者披露的训练预处理，不是本项目共同评价分母的删例规则，也不是实际数据已审计的证明。现有7728条DEV primary及unsupported政策不变；将来独立迁移的primary政策须在读预测性能前确定，不能追随某个训练器的可预测性缩小分母。

本轮一次Exa定向检索请求5条，随后读取固定公开源码／目录，并以已安装ANNEVO固定checkout的Git文件表及源码核实训练工厂；未执行训练、加载新权重、读取H5／候选DNA／注释或Setaria。Hugging Face目录抓取工具不支持该JSON内容类型后，使用同一官方API直接GET成功；GitHub当前main API限流403，未以此声称当前main没有更新。ANNEVO公开目录网页读取失败后，转读既有固定checkout和固定raw源码，失败不作缺数据证据。后续不重复审阅这些未变的入口；只有出现明确对应实际权重的新manifest或发布证据时再更新资格。当前任务继续等待A31原生DEV整链并执行已准备诊断，上游unknown不成为冻结开发比较的额外阻塞，也不授权更换encoder或清洁zero-shot主张。

## GENATATOR 候选对照的阶段曝光与版本

2026-10-07。GENATATOR是公开候选对照，尚未成为本项目运行资产。[作者流水线](https://huggingface.co/AIRI-Institute/genatator-pipeline)支持raw FASTA到GFF3的完整注释，不能因其分割模型只用基因序列训练就排除完整流水线。公开Hugging Face API本次返回流水线revision 7976d5cf9bdcdc3e03327f5bf8a3e5b8487c7fce，lastModified为2026-09-19；其[固定配置](https://huggingface.co/AIRI-Institute/genatator-pipeline/raw/7976d5cf9bdcdc3e03327f5bf8a3e5b8487c7fce/config.json)默认采用下列四个阶段。表中revision是公开仓库快照身份，不表示已下载、校验或执行对应权重。

| 阶段 | 固定模型卡 revision | 作者披露的监督训练范围 |
| --- | --- | --- |
| ModernGENA边界发现 | [df900a722e3945b8986bf8c5c88923a2cc3a3fd1](https://huggingface.co/AIRI-Institute/genatator-moderngena-base-multispecies-edge-model/raw/df900a722e3945b8986bf8c5c88923a2cc3a3fd1/README.md) | 39种哺乳动物的全基因组序列，包含基因间区域及mRNA、lncRNA注释；预测双链TSS与PolyA。 |
| ModernGENA区域发现 | [180c8ffcff74a0779f63e21fbd4ad0705ecd1039](https://huggingface.co/AIRI-Institute/genatator-moderngena-base-multispecies-region-model/raw/180c8ffcff74a0779f63e21fbd4ad0705ecd1039/README.md) | 同样披露39种哺乳动物、含基因间背景的全基因组监督；预测双链基因内区域与边界。 |
| Caduceus PS转录本类型 | [956d8c6c8c8d893e68082d344530efbc9bb9ab12](https://huggingface.co/AIRI-Institute/genatator-caduceus-ps-multispecies-transcript-type/raw/956d8c6c8c8d893e68082d344530efbc9bb9ab12/README.md) | 多物种基因序列，每基因一个转录本，无基因间区域；区分mRNA与lncRNA。该卡未列完整物种／assembly清单。 |
| Caduceus PS内部结构分割 | [2716f2b2cdaf72c9cd2b4b8ed069f37314ada0a7](https://huggingface.co/AIRI-Institute/genatator-caduceus-ps-multispecies-segmentation/raw/2716f2b2cdaf72c9cd2b4b8ed069f37314ada0a7/README.md) | 多物种基因序列，每基因一个转录本，无基因间区域；输出5UTR、exon、intron、3UTR及CDS。该卡未列完整物种／assembly清单。 |

两个发现阶段均声明人类chr8、chr20和chr21留出，非人类基因组包含长度大于100 kbp的染色体。这是这两个微调模型的作者披露，不能移植到类型、分割或上游预训练，也不证明四阶段对某一物种、assembly或同源序列均未暴露。后两张卡只声明multispecies，不据此假定四阶段用了完全相同的39种及划分；上述阶段卡不足以独立核实本项目未来面板的完整训练、验证及同源隔离。

固定流水线配置的model_revisions为null，阶段路径也没有各自revision；wrapper revision本身不足以证明四阶段的不可变模型身份。将来若执行对照，应分别固定流水线、四阶段权重及运行参数，不能把2026-09的流水线快照默认视作2026-06预印本评测所用版本。[UPSTREAM_MANIFEST](https://huggingface.co/AIRI-Institute/genatator-pipeline/raw/7976d5cf9bdcdc3e03327f5bf8a3e5b8487c7fce/UPSTREAM_MANIFEST.json)列的是源码归档文件身份，不是物种、assembly或训练划分清单，不能补齐训练曝光资格。

候选对照的有效比较单位是完整发现与结构注释流程，而不是由参考提供边界的孤立分割器；应在读性能前确定共同full-grid primary CDS链与背景口径，UTR和lncRNA另列。以上公开元数据不自动新增基线、读取DEV／独立面板或改变A40的TRAIN门槛；当前固定三条外部基线及各冻结结果保持不变。GENATATOR的raw-FASTA能力是作者公开功能说明，不是本项目已完成的性能或清洁zero-shot证据。

## Carbon-A 发布与当前固定编码器来源

2026-10-09的历史核对。[官方发布文章](https://huggingface.co/blog/HuggingFaceBio/carbon-annotator-genbank-genome-annotation)日期为2026-10-08。当时Carbon-A-1.2B固定revision为725713317756b1726ca026c281b31b9e59c4b25a；其[provenance.json](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/725713317756b1726ca026c281b31b9e59c4b25a/provenance.json)明确列出source_repository=GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview、source_revision=b0483c23b6b63787b61a6d3a204a9b517d6ba345，并将model.safetensors列为保留的源文件。这是作者对当前原始权重的来源映射，不是另一个无关编码器；该次核对没有下载新权重或独立做二进制比较，也没有用Carbon-A的兼容性实现替换当前运行代码。

### 注释监督及目标assembly资格

[技术报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)第18–23页披露“annotation as pre-training”：模型在RefSeq六类真核生物的双链CDS二分类上训练，CDS合并各isoform、无CDS背景窗口保留；初始哺乳动物模型随机初始化，再逐步累计其他类群。报告列出2047 taxids、2055 assemblies，使用98304 bp窗口。这支持“注释监督预训练表示”，不支持把当前固定encoder称为纯无监督GLM。GENERATOR架构来源也不等于先加载了其因果语言模型预训练权重。

报告第3页区分28个来自训练快照的参考assemblies与14个2026年新增的时间留出assemblies。第36页Table S3中拟南芥GCF_000001735.4与水稻GCF_034140825.1均属于前者、没有New标记；这两个accession与本项目configs/m1_data_manifest.yaml一致。因此本面板不能宣称目标assembly未进入上游来源体系或为新时间留出。表中的训练快照归属仍不足以证明哪些具体窗口参与了梯度更新，也不能直接推断本项目DEV标签已泄漏。

[公开训练数据说明](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-training-data)将数据接到RefSeq GBFF/FASTA，说明H5不包含逐窗口accession元数据，而由独立catalog提供。[catalog入口](https://huggingface.co/datasets/HuggingFaceBio/Carbon-A-RefSeq-training)本次返回401，未尝试登录或绕过权限。精确窗口／参考版本、训练与验证划分、同源隔离仍为未核实；下游head的TRAIN／DEV分离不能补齐这些上游资格。

### 原始模型与完整结构流程不能混为一谈

[固定模型卡](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/725713317756b1726ca026c281b31b9e59c4b25a/README.md)及[infer_genbank.py](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/725713317756b1726ca026c281b31b9e59c4b25a/infer_genbank.py)的发布示例输出逐碱基CDS概率NPZ，没有实现完整GFF调用器；但报告第24–25页的完整Carbon-A流程明确包含transition TCN、Viterbi基因语法、confidence TCN及GFF3/CDS/protein FASTA。不能据模型卡就将Carbon-A排除为“只做分割”。该流程和我们的语法解码存在实质先例关联，单凭冻结表示、边界网络或Viterbi组合不足以主张新颖性。

报告流程窗口重叠6144 bp，而发布示例stride=49152 bp；二者不是已核实的同一运行配置。本轮未定位并固定完整流程的可执行资产或复现其结果；它是后续匹配比较的候选，不加入或替换当前冻结的ANNEVO／Tiberius／Helixer比较。作者各assembly性能、置信度过滤与本项目最长CDS primary、unsupported保留及exact-chain分母也未完成口径对齐，不并入本项目结果表。

### 下一项TRAIN诊断的决策依据

报告第40页Table S8的人类chr1消融中，98304 bp与6144 bp上下文的未过滤gene F1分别为0.782与0.501。这是作者特定流程的结果，不是本项目水稻差距的因果证明或同口径性能。当前encode_window按6144 bp独立编码；已有preflight_full_grammar_glm_context.py把多个窗口输出拼接到通常49152 bp的语法视图，仍使用6144 bp编码块，不能让编码器获得共同长程注意力。因此下一项隔离诊断优先检验同一固定权重的原生长上下文，在原TRAIN见证上比较表示和冻结head输出；不将“拼接后的长语法视图”称为已检验原生长上下文，也不重复原拼接预检。

当前full_grammar_coarse_channel.py对fine角色仍作硬约束，仅在fine缺失而region存在时使用固定粗观测通道。它不是对错误fine标签的软化处理；“普遍抵抗注释噪声”尚无实验支撑。部分注释的任务收益可以继续作为机制假设，但不能仅以局部／结构化／混合目标三臂对比证明该主张。

本次读取公开来源文字、模型代码和43页技术报告的相关方法／表格，并渲染核对Fig 12、Table S3及Table S8。三次定向检索均请求5条结果，其中一次传输失败后转用公开搜索；不把请求额度当作完整已审阅文献数。没有读取训练H5、候选DNA／参考、DEV或sealed/Setaria数据，没有执行新训练、改阈值或改作业。本来源更新只改变方法归属、暴露表述及下一诊断优先级；三条已提交目标函数训练、其TRAIN闭环与已准备的冻结DEV方案保持不变。

同日后续已完成[原生上下文TRAIN诊断](../../experiments/M28/full_grammar/reports/full_grammar_native_context_diagnostic.md)：原八见证、四条件、固定4608步权重，短块基线精确重放；共同编码改善物种汇总NLL和部分精确链命中，但每条件仍有两个窗口NLL变差，且收益不随长度单调增加。此处为本项目新执行的TRAIN结果，不是Carbon-A作者评价或新独立测试。下一步核验原1491个TRAIN视图的完整效果，不改变上述上游暴露资格或冻结DEV规则。

2026-10-09 14:46 UTC的资产追溯当时只定位了概率推理入口及数据库概率segment检索；“完整流程仍缺可固定重跑的transition/confidence资产和调用配置”是该时点的历史结论。随后同日17:10 UTC发布的新仓库快照已补齐完整原生流水线，当前接入状态见[完整Carbon-A资产与比较资格追溯](../../benchmarks/carbon_a/reports/carbon_a_pipeline_readiness.md)及下一节。数据库结构缓存的输入、覆盖与口径仍未获得匹配资格；补齐可执行资产不等于本项目已经复现作者性能，也不改变冻结DEV家族。

## 完整Carbon-A流程的分阶段训练暴露

2026-10-11的历史接入状态：固定release为[f6789ee79e6ef397de333577ffec1cb734496620](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/tree/f6789ee79e6ef397de333577ffec1cb734496620)，已取得完整调用代码、两份辅助权重及配置；主权重复用原b0483c23资产。当时原两物种／六条TRAIN的独立原生对照已准备，隔离环境与人工序列原生预检已提交，尚无这六条TRAIN的Carbon-A性能结果；现已完成的固定比较见[原TRAIN结果](../../benchmarks/carbon_a/reports/carbon_a_original_train_result.md)。运行收据继续在readiness及结果记录维护，不在此重复作业台账。

此次实际读取了固定release的README、provenance.json、两份网络实现及safetensors文件头。它们支持下列分阶段判断，不补写尚未核实的训练数据：

| 组件 | 已知监督或来源 | 尚未核实的暴露 |
| --- | --- | --- |
| 主模型及原双链CDS概率头 | 来源映射仍指向原b0483c23；技术报告的RefSeq双链CDS监督说明适用于该主模型。新实现只是作者声明的Transformers兼容性修改。 | 实际参与优化的窗口、训练／验证划分与同源隔离仍未知；项目两目标assembly在训练快照参考子集中的事实不变。 |
| transition TCN | README明确为已训练网络；其11077个参数按相邻碱基边界预测stay、进／出gene、进／出intron。权重头给出channels、dilation_count与outcome_names。 | 本次读取的发布来源映射及权重头没有列出该网络的训练／验证assembly、窗口、标签构造或同源隔离，不能沿用主模型清单替代。 |
| confidence TCN | README明确为已训练网络；其14433个参数读取CDS概率及整链解码角色，为完整gene输出结构完全正确的概率。权重头给出网络结构、输入和角色名称。 | 同样未由此次固定发布元数据核实实际训练／验证数据及标签来源；输出叫confidence也不意味着校准已在本项目面板验证。 |

固定[provenance.json](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/f6789ee79e6ef397de333577ffec1cb734496620/provenance.json)的preserved_source_files映射主模型、配置及分词器；annotation_pipeline仅列为additional_files，没有为两份新辅助权重提供对应源checkpoint或训练catalog。这是当前已读证据的限制，不是断言作者在所有来源都没有披露；本轮没有重新完整审阅技术报告，也没有读取受限catalog或训练H5。

因此，下一张TRAIN表可比较“共享主模型的不同完整注释流程”，不能称为纯解码器因果消融：原生上下文、原CDS头、transition监督、语法及confidence处理均与本项目流程有差异。未过滤分支仍使用已训练transition网络，不是无监督对照。仅raw FASTA的推理接入不需要当前样本标签，但不抹去这些上游监督来源。

Carbon-A固定输出全部genes，后续共同口径分别报告confidence≥0.1或NA与未过滤结果，不依性能选择阈值。辅助网络暴露尚未核实不阻塞原TRAIN诊断，但禁止据其结果宣称清洁zero-shot、独立泛化或“只靠新解码器”的机制。此更新不放行DEV／封存数据，不新增物种，不改变已有冻结评价家族。

## 当前来源表述的精确边界

2026-10-13。固定[来源映射](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/f6789ee79e6ef397de333577ffec1cb734496620/provenance.json)仍将主权重接到原b048 revision。[技术报告](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)5.1／5.2及第36页Table S3与上述记录一致；表中两个项目assembly均无New星号，布局已核对。应写“已披露注释监督及目标assembly快照归属，具体优化窗口／验证／同源隔离尚未核实”，而非笼统的“GENERanno暴露未知”。该归属既不证明具体DEV标签泄漏，也不支持未见assembly、未见物种或同源隔离主张。
