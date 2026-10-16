# M27：先行研究与当前问题的新颖性边界

2026-09-26。M27 配对原生推理运行期间完成的定向核查。此文件不是结果分析，不修改 c284cd0 预注册的输入、分类、顺序、预算或停止规则，也不根据未完成病例选择研究方向。

## 判断

**“基因预测器可能通过改变结构绕开致损变异”已经是明确的先行研究结论。** 不能把它重新命名为“完整性偏置”就作为新概念；“用人工早停测试注释器”和“神经注释器对序列扰动敏感”也不是空白。

M27 仍有一个有限用途：检验固定现代原生流程在配对同义/PTC干预下是否出现可严格归属的响应，为是否值得进行后续机制确认提供证据。它本身既不是新方法，也不能仅凭换了现代工具就获得 Nature Communications 层级的新意。当前没有完整配对结果，不预判阳性或阴性。

## 最接近的工作

| 先行工作 | 已有证据与适用范围 | 对本项目的直接限制 |
|---|---|---|
| ACE，2017 | 指出传统 de novo gene finder 在 ABO O 等位上会修改阅读框以避开内部 stop；用参考注释和个人变异重建等位基因，不假设所有等位均功能完整。 | 不能声称首次发现完整 ORF 假设掩盖损伤，也不能将“保留损伤”本身当新设计原则。 |
| SGRF / ACE+，2018 | Supplementary Methods S1 描述在约19,000个人类基因的原预测 CDS 中随机置入早停，再运行 HMM gene finder，检查是否改剪接以绕开早停；主体方法用 exon-definition 特征，不依赖保守编码框，输入含参考注释和 phased VCF。 | 不能声称首次人工 PTC 干预或首次将剪接与编码完整性分开。其原 WT 不要求匹配参考，且不是本项目的同坐标 syn/PTC、整染色体现代原生流程。 |
| Helixer，2020/2021 | 对一个木薯基因做 in-silico mutagenesis：遮蔽剪接/start/stop motif 或扰乱编码潜力，观察神经网络逐碱基响应；去掉 stop 后可转向邻近下游 stop。该早期论文尚将完整 transcript 后处理列为后续工作。 | 不能声称首次对神经基因注释器做扰动解释；逐碱基神经响应不等于完整工具链响应。 |
| Helixer，2025 | 神经分数与 HMM 后处理输出完整基因；原文明示逐碱基性能不自动保证后处理性能，并分别评价。 | “局部识别不等于完整基因准确”不能单独作为新贡献。 |
| ANNEVO，2026 | Extended Data Fig.6b 描述修正错误 splice site 所造成的参考 premature stop，并有 RNA-seq 支持。 | “修复早停”不天然等于隐藏真实损伤；参考错误修正与真实等位失活必须用独立证据区分。 |
| TOGA，2023 | 同源/整基因组比对驱动的注释与 orthology 联合分析，显式区分 intact、missing 和 inactivating mutation 等状态。 | 损伤感知的基因注释已存在；它利用参考与比对，不能与本轮纯序列工具混成相同输入基线。 |
| LiftOn，2025 | 结合 DNA 与蛋白比对，以蛋白一致性选择结构和 ORF；论文图1明确涉及 frameshift、早停、下游重新起始等情况。 | 蛋白完整性导向的后处理及 restart 已有方法语境；不能把所有结构变化都标成不合理绕开，也不能把同源辅助工具当无证据输入的直接公平对照。 |
| Brůna 等，2026 | 棉花和大豆泛基因组研究显示注释方法可主导推断的 PAV；相同 CDS 序列也可能被注释为不同结构，统一流程明显缓解但不完全消除，且已纳入 Helixer。 | “注释不一致影响泛基因组结论”同样已有系统研究，不能仅用这一宏观后果补足新意。 |

原始来源与定位：

- [ACE 作者机构原文](https://ucgd.genetics.utah.edu/wp-content/uploads/2017/01/btw799.pdf)，Introduction、ABO 示例及 Discussion。
- [SGRF 作者机构原文](https://ucgd.genetics.utah.edu/wp-content/uploads/2018/05/bty324.pdf)，Introduction、Methods 2.1–2.3；[原文及补充材料入口](https://pmc.ncbi.nlm.nih.gov/articles/PMC6198862/)，Supplementary Methods S1。本次 S1 内容来自该条目可检索的补充材料文本；主体原文另经作者机构 PDF 核对。不据此引用未经核实的绕开比例。
- [早期 Helixer 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8016489/)，Methods 2.9、Results 末段及 Discussion。
- [完整 Helixer 原文](https://www.nature.com/articles/s41592-025-02939-1)，Results 中 base-wise / postprocessing 比较。
- [ANNEVO 原文](https://www.nature.com/articles/s41592-026-03036-7)，Extended Data Fig.6。
- [TOGA 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10193443/)，基因状态判定及 assembly incompleteness / base-error 区分。
- [LiftOn 正式论文](https://doi.org/10.1101/gr.279620.124)，Genome Research 35:311–325，图1；[前版本条目](https://pubmed.ncbi.nlm.nih.gov/38798552/)明确链接正式发表版本，不能继续把整项工作标成仅有预印本。
- [泛基因组注释一致性正式论文](https://doi.org/10.1093/nargab/lqag011)，NAR Genomics and Bioinformatics 8:lqag011，Methods、Results、Conclusions；[2025预印本](https://www.biorxiv.org/content/10.1101/2025.08.14.670405v1)已明确标注正式发表，不能继续称“尚未同行评议”。

另核查了作者维护的 [TOGA2 数据轨道说明](https://hgdownload.soe.ucsc.edu/hubs/GCA/004/027/085/GCA_004027085.1/contrib/TOGAv2/TOGAv2.html)。它描述 exon-wise annotation、SpliceAI 等更新，但该页引用的 TOGA2 论文标为 In preparation；本项目不把工具说明升级为已发表实验证据。

## 等待完整结果后，真正值得问的问题

1. 若出现严格的 PTC=1 / syn=0 配对，先按既定规则只确认每工具登记坐标最前的一个。需要知道现代流程中的哪一环实际导致变化，而不是仅给最终 GFF 起一个机制名称。神经分数、DNA约束、动态请求和过滤是不同环节；当前原生输出只能给出整体响应，不能单独归因。
2. 即使原生响应重现，仍需判断是否存在超出旧 HMM 现象的、可操作的新结论：例如某一可辨识环节是否决定了何时保留或掩盖损伤，以及该结论能否指导工具选择或具体修正。这里只列待检验问题，不宣称已经存在这种机制。
3. 人工 PTC 不提供真实剪接真值。未来若讨论“掩盖真实失活”或生物学后果，需要独立的真实等位与转录/蛋白证据，并排除参考错误、组装错误、同源投射偏差。保持 canonical splice dinucleotide 不变，也不能证明剪接生物学不受影响。
4. 没有可归属的严格配对，就按合同关闭本固定先导；未知保留未知，不放宽锚、不补选位点。若只有现代工具重现已知行为，即使阳性也不能自动扩大为 Nature Communications 主线。

## 检索范围与执行边界

检索截至2026-09-26。Exa 当前连接可用：三个预先限定的方向各请求10条结果（现代神经注释器、等位感知结构预测、泛基因组/基因失活），随后针对 SGRF 原文定位补充请求5条，共35个返回位置，包含重复论文、工具页及不相关记录。该数量不是35篇独立研究，也不是系统综述的纳入样本数。

按是否直接涉及 gene structure、完整性假设、损伤或注释后果筛选，本文讨论8项研究及1份官方工具说明；采用原文、作者机构稿或官方资源，排除仅有聚合摘要和不相关临床变异分类工作。Exa 部分日期/文献类型元数据与原文矛盾，本文以原文发表身份为准，不沿用错误日期。部分 PMC/PubMed 直接访问触发站点验证，OUP直读不可用；使用公开索引和正常可访问的作者机构版本，未绕过访问控制。未搜索到更接近工作不等于证明不存在。

此阶段只新增文献说明。没有训练、没有新增模型运行、没有新生物学数据下载、没有 Setaria 访问；运行中的 M27 合同、代码和所有结果均不改。此说明与完整 M27 结果一起交给 Pro 第11轮商议，AI商议不代替独立复现或同行评审。
