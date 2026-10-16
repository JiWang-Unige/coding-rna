# 完整链竞争训练的先例与本项目边界

2026-10-03。结论：结构竞争训练、自然候选中的局部对比学习均已有先例；不能仅把BCE换成排序／对比损失就称为论文创新。当前值得检验的问题，是在DNA-only自由候选中，保持多基因共存和未知链保护的前提下，让正确完整CDS链胜过null及实际冲突结构，且收益保留到最终注释。

本轮Exa两个定向检索共10个返回结果额度，选取两个原始研究来源核对方法。10不是10篇全文审阅，也不是穷尽性新颖性检索。只查询公开方法名称与概念，没有发送项目数据、权重或未公开结果。事实与项目差异分开记录，不拼接原论文成绩做本项目SOTA排名。

## 原始方法事实

[Schweikert等，mGene，Genome Research 2009](https://genome.cshlp.org/content/genome/19/11/2133.full.pdf)，Methods中Layer 2及Training：用HSM–SVM整合边界、内容及长度特征，以结构化间隔使真实标签比其他合法标签得分更高；训练按当前参数寻找最大违例标签并扩充约束，预测通过动态规划求最高分合法结构。故判别式基因结构学习、正确／错误结构竞争及难例更新不是新概念。该结果来自其对应任务和评估，不能推论当前植物模型的性能。

[Gao、Dai与Callan，Rethink Training of BERT Rerankers，ECIR 2021](https://www.cs.cmu.edu/~callan/Papers/ecir21-Luyu-Gao-b.pdf)，第3节及第6节：Localized Contrastive Estimation从目标初级检索器的高排名候选采负例，以一条正例与该组负例的对比目标训练reranker，并检验候选来源变化的影响。它说明候选来源与竞争目标已有明确研究先例，不证明本项目class imbalance、null、长度或表示是单一根因；文本检索成绩也不能移植为gene准确率。

## 对本项目的具体约束

这是基于上述方法与当前代码的比较推断，不是原作者对本项目的判断。本项目一个窗口可以有多条不重叠正确基因，因此不能把全窗所有链softmax成一个互斥正确类别；主参考相对0也不是生物学真负，-1未知／保护链不能当作难负例。若研究竞争目标，应限定为实际span冲突的已知非匹配结构，并保留绝对null竞争与背景负担控制；若仅排序而没有绝对零操作点约束，正确链在非正gain处仍可能全部被拒绝。

现有ChainHead使用可加exon分、ordered-GRU非加性分和共享window null，A6采用已知候选均值BCE。A12正在原TRAIN上固定末权重诊断，不选择中间权重或替代标签。是否启动新目标，应由完整TRAIN诊断决定，另定一次固定seed／曝光／候选池的协议，不改已闭合A6/A7等实验。

即使新目标有益，也只是任务特定机制证据。论文贡献仍需完整raw-FASTA输出、与ANNEVO／Tiberius／Helixer同尺逐物种比较、足够效应与代价、独立迁移及暴露审计、机制消融共同支持。TRAIN读出改善不是泛化，DEV筛选不是独立测试。此记录不授权新数据访问、解除Setaria封存、公共发布或论文claim。

## ANNEVO 方法核对与下一机制选择

2026-10-06。本节更新研究判断，不改上面的历史记录或任何冻结实验。定向公开检索核对原论文和本项目实际运行的官方源码；Exa本次返回连接错误，随后使用公开网页检索。没有将项目代码、数据或未公开结果发送给检索服务，也不是穷尽性新颖性审查。

[ANNEVO 原论文](https://www.nature.com/articles/s41592-026-03036-7)的 Extended Data Fig. 1–2 明确包含卷积压缩、长程 Transformer、专家组合、碱基分辨率恢复及 Viterbi 结构解码。期刊正文在本次网页访问中仅显示订阅预览；这里的方法事实来自可公开读取的扩展图说明，不冒称已读取完整期刊 Methods。

实际原生比较固定源码为 42c920f70f52f96c7c06c665a639ec47d9304fff。[src/utils.py](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/src/utils.py)明确把Magnoliopsida路由到[ANNEVO_plant.py](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/model/ANNEVO_plant.py)：卷积压缩、sin/cos位置编码、六层Transformer、八专家中top-2加权、上采样、15类逐碱基输出。该植物分支不是本项目的递归query边界指针，也不是[ANNEVO_seq.py](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/model/ANNEVO_seq.py)里的RoPE架构。后来公开的sequence-only训练说明不能直接当成实际植物权重的完整训练设置；权重来源和暴露仍按[既有记录](../reports/M28-METHOD-RESTART/pretrained_exposure.md)分别报告。

[固定 gene_decoding.py](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/src/gene_decoding.py)先从逐碱基概率定位、扩展与合并可能的基因区间，再调用[HMM.py](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/src/HMM.py)的状态动态规划；局部概率区分 coding phase、intron、donor、acceptor、start、end。因此相位／边界证据、长程表示和结构解码已有直接对标先例，A26不能仅以加入这些输入宣称创新。HMM的精确性限于给定区间和状态模型，不等于任意候选图或历史依赖评分都获得全局精确解。

[固定训练说明](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/model_training/README.md)区分中心监督区和两侧上下文；[固定注释入口](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/annotation.py)的植物配置为中心30720 bp、两侧各5120 bp。不能把另一个更长配置视为当前植物运行的设置，亦不能只按本项目24576 bp窗口和对方中心区宽度声称上下文公平。此差异暂不授权变更当前窗口或重新运行对标。

本项目的当前完整输出仍不具竞争性：[FULL Q](../reports/M28-METHOD-RESTART/full_train_result.md)在同尺 DEV 的拟南芥／水稻完整链 F1 为55.17%／58.90%，合并56.28%；[原生 ANNEVO](../reports/M28-METHOD-RESTART/native_baseline_result.md)为76.2651%／74.4318%，合并75.7025%。两条DEV已用于研发，这些分数不是独立测试。旧A3 greedy自由池的理想上限只约束该旧池，不冒称当前FULL或beam联合池的上限。

据此作出的研发判断是：不把增添MoE、相位输入或放宽beam当成核心贡献；先用已启动的A23相对边界和A25起点监督判断生成瓶颈，用A26判断显式路径证据是否改善固定选择，再由A27全量遗漏链诊断区分局部必要条件失败与仍需查明的路径竞争。A27有局部通过query不等于自然搜索返回，亦不直接证明EOS是原因。下一独立机制须针对实际读出的失败位置，保留已闭合FAIL，不事后救参；只有完整raw-FASTA收益及独立适用范围才能支撑论文主张。这是本项目的研究选择，不是ANNEVO作者的结论。

## A30 搜索保留损失与既有搜索训练

2026-10-06。结论：训练时观察beam错误、在参考前缀丢失处更新、对搜索竞争者施加margin，均已有方法先例。A30不是BSO或LaSO的直接复现；其任务特定差异需要实验支持，不能仅凭实现不同就宣布新算法。本节比较两个原始研究的方法段落，不构成穷尽性新颖性审查。

[Collins与Roark，ACL 2004，第4.2节](https://aclanthology.org/P04-1015.pdf)的early update在正确部分解析离开过滤后的候选集时结束该例的搜索，使用部分解析做perceptron更新，再处理下一个训练例。“首次搜索失败处学习”不是新概念，其线性模型性质也不自动适用于本项目神经网络。

[Wiseman与Rush，EMNLP 2016，第4至4.3节](https://aclanthology.org/D16-1137.pdf)的BSO用去掉softmax的非概率前缀分数，对gold与beam第K名的margin违例施加损失；末步要求超过最高分错误完整序列。它按LaSO方式在违例后从gold历史重新扩展，可收集多个违例，而非early update的首次失败后结束该例。

下面的A30事实来自[当前训练实现](../scripts/experiments/M28-METHOD-RESTART/query_search_retention_training.py)及[路径分数回放](../scripts/experiments/M28-METHOD-RESTART/query_search_score_replay.py)，比较关系是本项目的推断。

| 比较点 | BSO方法 | 当前A30实现 |
| --- | --- | --- |
| 搜索分数 | 非概率前缀分数，无最终softmax | 保留原指针log-softmax与EOS log-sigmoid的累计路径分数 |
| 额外损失触发 | 搜索中的margin违例，末步另有完整序列条件 | 每个原Hungarian正query最多一个实际首次丢失点；已返回链的额外项为零，即使未满足额外margin |
| 违例后的搜索 | 从gold历史重启后续beam扩展 | 原自然搜索不注入、不保活参考，不重置beam；参考仅用于观察，随后回放目标与cutoff分数 |
| 竞争者与汇总 | 原论文任务成本及违例集合 | 仅相对原标签判为0的cutoff接受hinge；1及-1保护。按全部已分配正query平均，包括返回、保护和不支持的零项 |

A30的额外项为`relu(margin + cutoff_score - target_score)`；raw grammar不支持的参考仍保留原teacher监督，不作为新负例，也不从原监督中删除。局部top-k、跨历史候选保留和原完成分数界属于实际观察位置，不等于修改搜索预算或保证参考链最终保留。归一化分数未变，故不能借BSO的去softmax机制宣称A30消除了label bias，亦不能引用原论文NLP成绩支持基因注释性能。

目前可检验的差异是：在真实DNA约束搜索的丢失位置，把注释资格与未知链保护用于额外竞争监督，并检查收益是否到达最终完整链输出。它尚不是已证实的论文贡献。原A30五项TRAIN门槛、末权重、曝光与停止条件保持不变；A32仍以A30全部通过且实际结果核验闭合为前提。任何TRAIN支持门槛通过都不能替代同尺raw-FASTA比较、匹配机制消融和独立泛化。此次文献比较不引入BSO论文的预训练或beam课程，不改变冻结协议，也不放行sealed数据。

## Tiberius Helixer 与 GeneDecoder 的结构训练边界

2026-10-06。可微基因语法、神经表示与CRF联合学习、预训练DNA表示接结构解码，均有直接先例；不能仅以加入HMM／CRF或GLM表示宣称新方法。这一判断约束后续机制选择，不改变A36已经冻结的训练、五项继续条件或末权重。

[Tiberius原论文，Methods 2.3至2.4.2，DOI 10.1093/bioinformatics/btae685](https://www.ovid.com/journals/bioinf/fulltext/10.1093/bioinformatics/btae685~tiberius-end-to-end-deep-learning-with-an-hmm-for-gene)明确区分先训练CNN／LSTM、再通过可微HMM微调；损失从HMM posterior回传至神经输出，结合CCE与exon类F1损失。其标准HMM的24个参数没有训练，transition对应按数据均值设置的几何长度分布。因此“端到端训练”不等于HMM转移参数也学习，“把结构层纳入梯度”本身也不是本项目创新。这里核对的是同一原论文的Ovid刊载正文；期刊主站本次只返回导航，不冒称在那里读到了Methods。

[Helixer正式论文，Architecture、HelixerPost及loss小节](https://www.nature.com/articles/s41592-025-02939-1)的神经训练使用逐碱基genic／phase交叉熵，并对类别及重要边界加权；HelixerPost随后从预测中识别候选基因区间，在区间内按生物语法解码。故边界强调、相位辅助和HMM合法结构均是直接对标已有机制；需区分神经层训练与HMM后处理，不能由“end-to-end”概述推断HMM已参与训练梯度。

[GeneDecoder预印本v1，2025-05-06提交](https://arxiv.org/abs/2505.03377)，[正文3.1、4.1至4.4及Conclusion](https://ar5iv.labs.arxiv.org/html/2505.03377)，用潜在CRF编码方向／codon相位和合法转移，联合学习DNA表示与转移势，以NLL训练、Viterbi预测；也报告GPN＋CRF与Nucleotide Transformer embeddings结果。实验按物种分别训练，输入为带1至2000bp侧翼的基因片段并限制极端长基因，按80%蛋白相似性做同源划分；主要MCC／F1针对位置标签，不是本项目完整染色体exact-chain F1。作者也指出短侧翼训练可能影响高intergenic背景表现。其成绩不能移入当前同尺榜单，也不能作为未见物种零训练迁移证据。

据此形成的项目判断是：如果后续研究绕开A3已完成候选池，直接从DNA边界证据解码合法结构，确实能检验当前重评分无法恢复缺失几何的范围；但“GLM＋CRF／语法层”不能单独构成论文贡献。需要用完整原生输入与背景、真实完整链收益、匹配机制消融和独立物种暴露边界证明差异，而不是更换输出名称或引用片段／碱基指标。具体下一机制仍应依据A36完整TRAIN核验和相应原生结果定位失败位置，在既有授权范围内另定实验；此比较不选择新loss、阈值、首训或数据范围，不部署A33例外，不放行Setaria、独立面板或公开发布。本轮仅进行公开方法定向检索和上述原论文段落核对，没有传输项目研究材料，也不是穷尽性新颖性审查。

## A44 有限集合目标的先例与含义

2026-10-07。GeneDecoder的arXiv v1提交日期仍为2025-05-06，但不能把这个日期当作该思路首次公开时间。[更早的OpenReview页面](https://openreview.net/forum?id=Rn50hCOX9XX)列相同三位作者，显示2023-02-01公开、2023-02-13修改及Submitted to ICLR 2023；[对应投稿正文](https://openreview.net/pdf?id=Rn50hCOX9XX)已描述潜在CRF、NLL训练、Viterbi解码和冻结GPN表示的探索。这是公开投稿证据，不是ICLR录用证明，也不把两个版本默认为逐字相同。因此，DNA表示接精确结构解码的先例至少已公开到2023年。

[ANNEVO原论文补充方法1至3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-026-03036-7/MediaObjects/41592_2026_3036_MOESM1_ESM.pdf)还给出注释错误掩码，以及focal、CDS／intron Dice和专家负载平衡损失；训练先用含genic信息的窗口，再纳入全数据，以验证集genic F1早停，随后按预测势作Viterbi解码。由此不能把“防止错误注释成为负例”或局部监督加生物语法本身称为原创。这些是论文版方法，不补齐当前更新植物权重的完整训练物种、assembly及实际配方；[固定资产暴露记录](../reports/M28-METHOD-RESTART/pretrained_exposure.md)中的unknown保持不变。期刊网页仍未提供完整主Methods，本补充依据的是已读取的原始Supplementary Information。

当前A44的含义由实际[chain_set.py](../src/m28/chain_set.py)、[structured_emission.py](../src/m28/structured_emission.py)和[local_structured_emission.py](../src/m28/local_structured_emission.py)确定：额外项仅在固定有限支持中的0／1已知候选上，对包括empty的全部span兼容集合计算logZknown；−1链不进partition或gold。它不是全DNA语法路径的训练配分函数，也不是将未知链作为潜变量积分的部分观测似然。预测才使用原完整DNA语法，且不把有限支持或标签作为输入。N为实际有效碱基、M为原局部可监督位置，两者不能因都涉及位置就互换；FP64角色gain和零操作点也不是校准后的gene概率。这些是当前实现的范围解释，不是新颖性证明。

据此，A44能够检验的是在相同容量、局部标签和draw曝光下，额外的已知完整结构竞争监督能否改善自然完整链；两臂增加的总标签信息和梯度并未匹配。仍须同时通过相对A31_A3及LOCAL的两组冻结五条件，不能仅凭损失下降、局部恢复或工程接通宣称结构收益。任一失败按[原配对协议](../reports/M28-METHOD-RESTART/a44_paired_fit_protocol.md)关闭；成功也只获得TRAIN机制继续条件，论文仍需完整raw-FASTA同尺收益、独立适用范围及机制证据。此次先例补充不改变运行中的损失、候选支持、checkpoint、阈值或停止条件，不选择新面板、增加训练或公开研究材料。

## 全路径似然与部分观测的先例

2026-10-09。全路径条件似然、未知标签边缘化和潜在相位均已有先例。当前完整DNA语法实现的差异应作为待检验的任务设计，而不是单凭配分函数或CRF名称宣称算法原创。本节不更改既有实验；固定4608步[原TRAIN诊断](../reports/M28-METHOD-RESTART/full_grammar_train_energy_diagnostic.md)和[DEV共同尺](../reports/M28-METHOD-RESTART/full_grammar_dev_baseline_result.md)仍是当前性能证据。

[DeCaprio等，Conrad，Genome Research 2007，Methods的CML及MEA小节](https://pmc.ncbi.nlm.nih.gov/articles/PMC1950907/)将合法基因结构作为semi-Markov CRF的输出，用全结构条件最大似然和forward-backward梯度训练，再以最高概率结构预测。它还比较最大期望准确率训练，包含强调剪接正确性的目标。因此“把完整结构竞争放进训练分母”及“让训练目标更接近结构准确率”不是首次提出。Conrad的线性feature权重、真菌及比较序列实验不能证明当前神经head的收敛、植物迁移或效果，也不作为本项目今日基线成绩。

[Tsuboi等，COLING 2008，第3至4节、式2至5](https://aclanthology.org/C08-1113.pdf)把每个位置的观测写成允许标签集合，将未知标签积分掉，优化相容序列集合的条件概率，即相容路径配分函数与全部路径配分函数之比。缺失位置不强填背景；歧义标签也不强选一个gold。作者明确指出这种边缘似然一般不再是凹目标。该研究实验是NLP，不是基因注释；它确立一般部分观测CRF原理，但不提供本项目完整CDS链、噪声注释或生物语法上的收益证据。

[GeneDecoder原文第3节、式1至4及Training](https://arxiv.org/html/2505.03377v1)用潜在CRF边缘化方向／codon-index，学习表示与允许的转移势，以NLL训练及Viterbi解码。潜在相位不是每个位置都没有genic标签；应与本项目缺失fine、已知coarse和unknown三类观测区分。其片段／短侧翼实验及位置指标仍不能移作whole-contig exact-chain比较。“GLM表示＋潜在相位＋精确结构解码”已有直接方法先例。

本项目当前定义来自实际[coarse_channel_nll](../src/m28/full_grammar_coarse_channel.py)、[可微配分函数](../src/m28/full_grammar_boundary_compiled.py)与[DNA状态转移](../src/m28/full_grammar_boundary_messages.py)：在固定有限DNA视图的完整非重叠gene-set宇宙内求和；fine是精确角色指示，coarse使用correct0.9／其余各0.05的三类观测channel，unknown贡献1。观测分子为所有合法路径的模型权重乘观测似然后求和，分母是同一合法宇宙的模型配分函数。由这些公式作出的比较推断是：hard部分观测是已有原理，soft coarse channel是观测模型的具体选择；两者组合本身尚无已证实的新颖性或性能优势。epsilon0.1不是由文献重新选择的值。

| 当前设计 | 应承认的先例或边界 | 真正仍需验证的差异 |
| --- | --- | --- |
| 全DNA合法路径NLL | Conrad已有全结构CML；GeneDecoder已有表示与CRF联合NLL | 在实际原生基因组背景中，是否改善完整链而非仅位置损失 |
| unknown路径边缘化 | Tsuboi已有部分／歧义标签的相容路径似然 | 对真实annotation不确定性是否更稳健；不能把缺失标签当真负例 |
| soft coarse观测channel | 它是额外似然因子的具体参数化，不是新CRF原理 | 标签信息、上下文与训练量匹配时，是否有独立收益 |
| 实际codon carry及字面DNA约束 | 精确性只针对当前状态图；生物约束本身并非新概念 | 合法性收益与拒绝真实非标准结构的代价，必须同时报告 |
| 跨块可微消息与whole-contig MAP | 数值精确／内存可执行不等于模型已充分拟合；有限视图fit不是whole-contig联合MLE | 长链准确率、实测代价及与局部视图的机制差异 |

据此，下一阶段不能只把“再训同一全路径NLL”包装成新模型，也不默认换为已有MEA、加权CE或更大encoder。应先通过受控TRAIN对照区分同目标增加优化量、表示变化与观测目标变化的作用；保持固定输入、初始条件和明确预算，报告全部15角色与完整链读数，不以背景占主导的一致率或归一loss作为拟合成功。不同因素同时变化的两个模型只能作候选比较，不能作单一机制归因。随后才按预先约定的开发评价选择路线，独立数据与上游暴露仍另受冻结边界约束。

本次定向Exa检索使用5个查询、每个5个返回结果额度，并读取上述三份原始方法文本；25是返回额度，不是25篇全文审阅，也不是穷尽性新颖性结论。仅检索公开方法名称与一般概念，没有向检索服务传输项目数据、代码、权重或未发表指标。该次核对时“不改训练参数、不增训”边界仍有效，以上对照尚未实施；封存访问、DEV调参和公开发布均未授权。这是解除阶段限制前的历史范围，不是对后续已授权原TRAIN实验的禁止。既有完整链缺口及原论文goal未完成的结论不变。

## Carbon-A 与多谱系 Tiberius 对论文定位的约束

2026-10-09。完整编码结构的模型先例不只包括ANNEVO、Tiberius和Helixer；Tiberius的当前适用范围也不只限于最初的哺乳动物模型。这两项事实收紧论文定位，不能把“预训练表示接生物语法”或“可用于植物”单独作为新贡献。

[Carbon-A技术报告，第24–25页、5.5.2至5.5.4节](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf)描述逐链coding probability经过transition TCN、序列约束Viterbi和confidence TCN，输出GFF3及CDS／protein FASTA。Viterbi按实际序列的起止密码子、相容剪接位点和跨内含子阅读框构建编码结构，并加入GT-AG、GC-AG及AT-AC剪接类别先验。置信度网络还直接以预测结构是否精确匹配参考作为监督；报告的过滤分析使用0.1操作点。因此，DNA表示、字面序列语法、精确解码及整结构质量监督都已有直接相关流程，不能把该方法缩写成一个逐碱基分类器。

该报告允许record边缘的不完整结构，其剪接范围、概率与transition评分、整基因置信度过滤也不同于本项目的complete-set角色gain与部分观测似然。这个区别来自报告与[当前Methods定义](M28_full_grammar_methods_working.md)的比较，不是原创性证明，也不推出任一方法性能更好。[固定资产追溯](../reports/M28-METHOD-RESTART/carbon_a_pipeline_readiness.md)在2026-10-09 14:46 UTC尚未定位齐完整重跑所需的transition／confidence资产与调用入口；这不是全网不存在的结论。已发布coding-probability示例和作者报告的完整流程须分开，不能拿自建简化调用器冒充Carbon-A原生基线或把作者成绩填入本项目共同尺。

[Tiberius作者仓库的模型列表与How to Cite](https://github.com/Gaius-Augustus/Tiberius)现列出哺乳动物、脊椎动物、昆虫、被子植物、真菌、硅藻和绿藻模型，另引2026年bioRxiv的《Accurate ab initio gene prediction in eukaryotes with Tiberius in multiple clades》（DOI 10.64898/2026.04.24.720536），与2024年原方法论文分开。故后续“更广适用”主张必须针对实际模型／谱系和上游暴露成立，不能仅由原论文的哺乳动物范围推出。仓库所述覆盖也不自动等于这些谱系在本项目上的独立泛化证据。

当前研发仍先闭合已授权原TRAIN中的目标、上下文、原生角色证据与采样对照，并用完整链和背景负担判断是否有真实收益；文献差异不替代机制消融。固定原生基线和[既定DEV候选方案](../reports/M28-METHOD-RESTART/full_grammar_sampling_context_dev_plan.md)保持不变，不因本节添加、替换或删除候选，不更新运行权重，不读取DEV／封存数据，也不发布研究材料。若最终只有碱基拟合或工程可执行性的改善，仍不足以支撑竞争性完整注释或更广适用这一论文核心主张。

## 密码子能量与编码内容评分的先例

2026-10-12。编码内容、阅读框相关序列统计及可学习结构feature权重均已有直接方法先例。新增64参数的拼接CDS密码子能量是具体的机制候选，不能仅凭加入这一项宣称新算法。

[Stanke与Waack的AUGUSTUS原文，第5–6页的模型小节](https://gobics.de/mario/papers/GenePred2003.pdf)使用四阶插值Markov模型描述编码内容，初始及内部3′内容组件具有三周期／阅读框依赖，并分别给TAG、TGA、TAA以24%、48%、28%的停止密码子概率。因此，在合法结构中加入序列组成或停止类别偏好不是新概念。这是2003年原模型的定义，不代表当前全部AUGUSTUS版本或本项目植物上的性能；其生成式内容模型也不等同于这里的自由条件能量系数。

[Vinson等的Conrad早期原文，NIPS 2006，第2至3节](https://proceedings.neurips.cc/paper_files/paper/2006/file/65a31da7ede4dc9b03fb5bbf8f442ce9-Paper.pdf)先把训练数据确定的生成式feature作为CRF输入，再通过判别训练学习feature权重，并加入非概率的EST及比对gap指示。第2.1节区分feature内部参数与条件似然权重的两阶段训练；第3.1节的比较模型使用多序列比对。它确立混合序列证据及可学习结构权重的先例，但不是本项目64个自由codon系数与神经head联合更新的直接复现，也不把该文的真菌／比较序列成绩移入今日raw-FASTA共同尺。此处核对的是可访问的2006原文；2007扩展版的SMCRF／CML比较仍见上节，不把两版默认为同一文本。

当前[实现](../src/m28/full_grammar_codon.py)及[固定协议](../reports/M28-METHOD-RESTART/full_grammar_codon_factor_preparation.md)在不变的gene-set宇宙内增加`Σk n_k(S;x) β̄_k`：完整拼接CDS的每个密码子计一次，包括起始ATG和末端stop；实际carry允许普通密码子跨内含子闭合，内含子及background不计此项。64个原参数分别按61个sense／3个stop去均值，零初始化嵌套原模型。原fine／coarse／unknown通道、全部路径分母、有限视图、fresh采样与joint编码不变；与fresh-joint的配对包含新增参数导致的梯度流和共同clip变化，不能归因为只改变一种生物学统计。

由当前公式作出的解释是：去均值移除了两个组内统一偏移方向，不等于经验codon usage、概率归一化或长度／基因数中性。每个基因必有一份ATG起始项及一份终止项，普通ATG还会在内部重复；非均匀组成与CDS长度也影响总和。系数与神经角色势联合学习，不能单独把正系数解释成已证实的生物偏好或独立机制。相应数学及梯度定义已补入[Methods](M28_full_grammar_methods_working.md)。

本分支的人工枚举、原生可丢弃更新及严格raw接入检查已经完成，正式fit仍在运行；尚无新的最终权重准确率。研究判断仍取决于固定最终9216步读出、原六条TRAIN的完整链与背景负担、相同参考总体和逐物种结果，而非参数符号或训练loss。三个新分支的[独立DEV准备方案](../reports/M28-METHOD-RESTART/full_grammar_joint_controls_dev_plan.md)与旧八／十方法结果分开；新研究DEV尚未执行，仍禁止DEV反馈调参、封存访问和对外发布。此次仅对两个原始方法的相关段落做定向核对，不构成穷尽性新颖性审查；检索未传输项目数据、代码、权重或未发表结果。
