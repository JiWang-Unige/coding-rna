# M27：等位序列干预的准备阶段

日期：2026-09-26。状态：PREPARATION_FROZEN；仅 CPU 参考材料准备与接口核实，**不是模型实验结果，也不冻结尚未解决的完整推理协议**。

## 为什么转向这里

Pro 第八轮已读公开固定提交 0b51d7a 的 R8 合同与完整结果，建议结束当前实现的修补，测试“内部提前终止变异是否诱发剪接结构绕行；神经分数与序列解码约束分别造成什么”。我们接受有限可证伪试验的方向，不接受把它当成已经具有 Nature Communications 新意的结论。

ACE 2017 已指出传统 gene finder 可为维持完整编码结构而错误处理失活变异；后续 SGRF 2018 直接用 exon-definition 信号避免编码保守假设。因此“ORF 约束可能掩盖损伤”不是我们的新发现。潜在新增贡献只能来自现代神经工具中可重复的分量干预、对实际等位序列的可靠验证及实用纠正原则，目前三者都未成立。

来源：[ACE](https://pubmed.ncbi.nlm.nih.gov/28011790/)、[SGRF](https://pubmed.ncbi.nlm.nih.gov/29701825/)。
Pro 是协作建议；代码与实验证据由本项目独立核验，没有把浏览器建议视为独立复现。

## 本阶段的 starting_point_assessment

- 已观察：ANNEVO 源码固定 37bdd9aa62ddf24fa55941fb827061f7ed49ce53（v2.3.2），Mammalia 权重已有 56 MiB；prediction.py 写入双链 HDF5，decoding.py 独立接收 DNA 与 HDF5。
- 已观察：Tiberius 源码固定 8c49fd029909eaa6098160247363e68e5ecde698；存在 tiberius_2.0.5.sif，源码有 lstm_prediction 和 hmm_prediction 两参数路径。不能仅凭源码目录认定容器相同，准备作业将读取实际版本与接口。
- 已观察：Tiberius 主流程还会按预测边界重新推理，并按 DNA 删除内部 stop 的输出。简单替换一次 HMM 参数并不足以说明端到端因果归因；必须把重推理、后处理及缓存输入来源一起控制。
- 已观察：原项目物种清单中没有本次 human chr22 输入；本阶段从公开参考重新获取，不读取任何 Setaria、旧 sealed 数据或模型预测。
- 未确认：两模型的完整训练物种暴露、当前容器内 mammalia 权重位置、端到端交叉解码重放一致性及实测 GPU 时间。本阶段不宣称 zero-shot。
- 文献已明确：固定 GENCODE 49 是 GRCh38.p14 的历史版本，不是“最新注释”。[固定发布页](https://www.gencodegenes.org/human/release_49.html)。

## 本阶段的 analysis_plan

决策：参考上是否存在足够清楚的配对干预单位，以及现有软件是否能进入小规模 WT 重放验证。实验单位是参考基因座，不是独立受试者或物种。

数据与操作：

1. 下载 UCSC hg38 chr22 FASTA（应为 50,818,468 bp / NC_000022.11）和 GENCODE 49 comprehensive CHR GTF。保留原始压缩文件，提取 chr22；大写化只用于未来无掩蔽、sequence-only 权重，不生成新 repeat track。
2. GTF 的 CDS 不包含 stop_codon；显式合并并验证起止密码子、完整 ORF、CDS phase、正负链坐标。排除部分 CDS、含内部 stop、selenocysteine、歧义碱基。记录所有排除原因。
3. 基因必须是 protein_coding，所有 protein_coding 转录本均可评估且只有一个去重后的 CDS 链。这个定义**不等于所有 RNA isoform 唯一**，NMD/其他 transcript_type 不进入这组链。
4. 在内部 CDS exon 选择 TAT/TAC：同一个基因组碱基可分别产生同义 Tyr 与 TAA；位点在 CDS 的 20–80% 位置，整个密码子距所在 CDS exon 两端至少 30 nt。仅要求目标基因链上跨编辑碱基的 GT/AG 二核苷酸不增不减；不宣称未改变 splicing regulatory motif 或反链预测。
5. 每基因选择距 CDS 中点最近的合格密码子，再按基因组顺序的等分位规则最多选择 24 个；不读预测，不补选高分位点。少于 24 则保留实际数目并报告，不擅自更换染色体。
6. 只输出位点清单、排除台账和 WT FASTA。**不生成合并了多个编辑的突变染色体，不运行模型，不训练。**

后续设计边界（尚未进入执行）：

- 各位点需要 WT / 单独 synonymous SNV / 单独 PTC SNV。不可把多基因同时编辑后当成单变异效应，也不改变 native window/重推理边界以迎合结果。
- 必须先证明 WT 的分离式执行与固定 native 执行一致，再做神经输入源 × 显式序列源交叉解码。混合输入只是机制干预，不是自然样本或可比准确率。
- 所有选择位点都保留，分别报告两模型各自 WT 正确与共同 WT 正确分母；WT 失败不是 mutant 失败。正常终止/截短不得归作“错误绕行”。
- Pro 的 >=20/24 共同 WT 正确和每工具 >=5 绕行是建议性的投资门槛，并无样本量依据。本准备阶段不采用它们，也不事后看结果降低门槛；模型实验的可判断标准需在第一次 WT 运行前写定。
- Pro 建议的总计 <=8 GPU 小时 / <=32 GiB RAM / <=50 GiB 存储暂作为未来先导设计上限，需先核算是否能容纳完整上下文与所需缓存；本阶段不消耗 GPU 预算。
- 如果 Tiberius 原生动态重推理使有限缓存无法可靠分离，停止该工具机制归因，不能从最终 GFF 变化补做因果故事。是否采用仅端到端表型试验需另行明确科学目标。

## 执行范围与停止条件

输出：outputs/M27-ALLELIC-PREP-R1；数据：data/m27_grch38_chr22；源码：scripts/experiments/M27-ALLELIC-INTEGRITY。

Slurm：private-teodoro-gpu，2 CPU，8 GiB，45 分钟，0 GPU；预计数分钟至十几分钟（主要下载），上限 1.5 CPU 小时。新输入及台账预计 <1 GiB，上限 2 GiB。不下载 Tiberius 权重，不更新环境。

坐标/phase/配对变异、正负链和选择函数先通过小型单元测试。下载/实现错误保留原始失败记录，只允许一次不改变上述科学范围的工程修复；不能因为位点数少、WT 难以预测或无绕行信号而换数据/阈值/种子。参考准备完成不等于干预假说成立。下一步取决于合格位点数与实际容器接口，随后带具体事实回 Pro 商议。
