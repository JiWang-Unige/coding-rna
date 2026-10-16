# A37 联合上下文整链模型的 TRAIN 负结果

2026-10-07北京时间。工程与实际数值核验完成，原五项继续条件仅通过2项，科学判定FAIL并关闭A37。不运行已准备的原生DEV，不挑其他checkpoint、不追加曝光、不改loss、gain阈值或权重。此结果只关闭这个固定候选；raw-FASTA模型、独立泛化和论文总目标仍未达到。

## 原条件与同池计数

粒度是原1491唯一TRAIN窗的自然候选／窗口配对，不是去重基因或完整染色体F1。原16321自然链的几何、顺序、greedy／beam4来源、原标签、gain0及span非重叠DP均保持不变。基线仍是已核验A31_A3（readout13398521／CPU13398524），不退回较弱对照或重判旧FAIL。

| 原继续条件 | A31_A3 | A37 | 结果 |
| --- | ---: | ---: | --- |
| 拟南芥已选正确配对严格增加 | 1531 | 1396 | FAIL，减少135 |
| 拟南芥已选known_nonmatch不增加 | 347 | 125 | PASS，减少222 |
| 水稻已选正确配对严格增加 | 478 | 431 | FAIL，减少47 |
| 水稻已选known_nonmatch不增加 | 251 | 114 | PASS，减少137 |
| 合并已选至少5-exon正确配对严格增加 | 700 | 639 | FAIL，减少61 |

eligible reference/window分母仍为拟南芥2172、水稻744，合计2916；至少5-exon分母854／232，合计1086。自然池正确1751／629、known_nonmatch6305／6261、未知853／522保持完整，不能按新预测缩小分母。

拟南芥全部已选2193→1675、未知315→154、至少5-exon正确548→499；水稻全部已选861→600、未知132→55、至少5-exon正确152→140。合并已选3054→2275，正确2009→1827（−182），known_nonmatch598→239（−359），未知447→209（−238）；两臂各满足已选＝正确＋known_nonmatch＋未知，下降779与三个分项一致。known_nonmatch只是primary参考相对的已知非匹配，不直接等同生物学假阳性；未知不计入错误或新基因。

观测结论是固定零操作点下选择更少，已知非匹配下降同时丢失正确与长链配对，未满足预声明目的。不能以较少错误替代三项失败，或从本结果断言联合表示、整链目标或GEN本身无效。相较A36，结构、容量和监督候选域都不同；相较A31_A3还有额外上下文／头训练曝光，不是匹配容量的因果消融。逐候选转移见后文；阶段定位仍不能把较少选择直接归因于null或某个分量。

## 实际产物与数值证据

fit13410061、natural13410739、CPU13410742的root／batch／extern均实际COMPLETED0:0，分别1:15:53／0:16:59／0:00:27，全部gpu034。fit和natural实际private1RTX3090／2CPU／16G，CPU实际0GPU／2CPU／8G；batch MaxRSS分别2213712K／1590504K／835128K。三个STATUS／JOBID与源作业相符，科研源未被重提或覆盖。

源fit固定fresh612324参数、三遍原1536 draw及末004608；training.jsonl实际4608行，001536／003072／004608三个checkpoint各7389798字节。自然natural_scores.jsonl和windows.jsonl各1491行，全部16321自然链，无参考注入、新搜索、新GEN、新R4、DEV、Setaria或A33访问，零新optimizer；全自然分数先关闭后才连接原保存标签。

实际CPU四项证明均true：原draw／importance／counts及4608更新；三模型checkpoint的有限FP32 tensor与分支optimizer曝光；原自然身份及精确FP32保存分量；原labels／DP／较强基线／2916和1086分母及五条件重算。SharedContext及null实际参数step4608，exon／链GRU／additive／nonadditive为4602，junction两分支为4425；这是候选调用决定的分支曝光，包括未知多exon候选，不一概假定所有头参数4608。

保存组件按原FP32逐值严格核对exon＋junction＝additive、additive＋nonadditive−null＝gain，没有放宽容差、改FP64或重跑新模型前向。GPU与CPU计数／分母／门槛完全一致，主再由两个实际summary独立重算五个整数条件、物种与合并小计，结果仍为2／5。这里的独立整数算术不是独立神经前向复现或第三方算法评阅，工程完成、数值一致与科学FAIL分别成立。

完整证据根为outputs/M28-A37-CONTEXT-CHAIN-FIT、outputs/M28-A37-CONTEXT-CHAIN-TRAIN-READOUT、outputs/M28-A37-CONTEXT-CHAIN-RESULT-CHECK；原checkpoint、4608训练行、两个1491行JSONL及全部summary／日志保留。原科学合同见[a37_context_chain_protocol.md](a37_context_chain_protocol.md)和[a37_context_chain_readout_protocol.md](a37_context_chain_readout_protocol.md)。已准备的[DEV合同](a37_context_chain_dev_protocol.md)及接口只是工程产物，因TRAIN FAIL不执行；两个DEV输出根实际不存在。

## 监控与下一判断

固定watcher于16:15的定时检查看到natural完成、CPU仍运行，没有发整链终态交接。主随后读取实际三根终态及科学结果，已直接消费这一cohort；同heartbeat于16:21 UTC原生暂停，读回PAUSED，ID／固定对话／完整prompt／15分钟／通知字段及created_at均保留。没有虚称终态自动消息已发送，不让下一次定时运行重复交接，监控对话保留复用。

后续阶段诊断仅配对已闭合的原自然分数、选择与标签；不重算模型／GEN，不打开参考或新数据，不改变A37停止条件。下一独立模型须依据实际失败位置另定机制，不能把诊断当阈值补救或已证实因果。A36仍关闭，A33例外未部署，独立面板、Setaria封存和公开发布边界不变。

## 已保存分数的配对诊断

diagnose_context_chain_attrition.py仅读取闭合的A37与A31_A3自然分数／选择及A37保存的原标签，复用已通过三项人工检查的transitions纯计数逻辑，不重写分类或重复未变的人工测试。统计保留、丢失、新增及双方未选，分别区分新gain≤0（含恰为0）与正分未选，另列≥5-exon正确、known_nonmatch与未知。实际全人口须对齐1491窗／16321配对、原两臂已闭合计数及2916／1086分母，否则不输出完成诊断。

这是决定下一机制优先处理链对空集还是正分竞争位置的描述性证据，不是null／nonadditive的因果消融。入口不加载权重或feature tensor，不打开TRAIN catalog、参考、原始DNA、GEN、DEV或新标签，不创建optimizer，不改变任何原分数或选定模型。

SBATCH沿用generanno／private0GPU／2CPU／8G／5分钟／no-requeue及原排除节点，唯一输出M28-A37-CONTEXT-CHAIN-ATTRITION。同类A36只读人口诊断实际36秒；本轮只产小型summary，不建缓存。短诊断由主直接消费，不重启已暂停的长期心跳；原A37仍科学FAIL关闭。

实际新源码／SBATCH／结果记录已明确同步canonical，远端AST及bash语法通过；提交前五项活调度／预约检查无相关维护、并发或输出冲突。仅提交sbatch --parsable sbatch/M28-A37-CONTEXT-CHAIN-ATTRITION.sbatch一次，真实13412326；即刻RUNNING／gpu034／private0GPU／2CPU／8G／5分钟／Requeue0／Restarts0。最终root／batch／extern均COMPLETED0:0／0:16，batch MaxRSS763404K，STATUS／JOBID／summary相符，stderr为空。

全部1491窗／16321配对准确复现原两臂已闭合计数、五条件及2916／1086分母；主再核对每个物种和pooled的六个互斥阶段之和＝pool、旧已选＝保留＋两类丢失、新已选＝保留＋新增。没有新模型、feature tensor、GEN、原始参考或新标签访问，只复用已核验的保存标签，也没有新optimizer；不把未变的人工计数测试再跑一遍当新科学证据。

### 实际正确配对转移

“丢失”是同一候选身份在A31_A3已选而A37未选；“新增”是相反方向，不是生成了新几何或发现了新基因。新gain≤0包括恰为0。

| 正确配对范围 | 旧已选 | 保留 | 丢失：新gain≤0 | 丢失：正分未选 | 新增 | 新已选 | 净变化 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | 1531 | 1338 | 168 | 25 | 58 | 1396 | −135 |
| 水稻 | 478 | 392 | 64 | 22 | 39 | 431 | −47 |
| 合并 | 2009 | 1730 | 232 | 47 | 97 | 1827 | −182 |
| 其中至少5-exon | 700 | 607 | 69 | 24 | 32 | 639 | −61 |

合并丢失279条原已选正确配对，232／279＝83.2%在A37 gain下非正，47／279＝16.8%为正分未选；新增97条不足以抵消损失。≥5-exon丢失93条，其中69／93＝74.2%非正，24条正分未选；新增32条，净减61。双方均未选的正确配对中还有221条非正和53条正分，原正确自然池2380＝1730＋232＋47＋97＋221＋53；≥5-exon自然池788与原1086参考配对分母不能互换。

known_nonmatch保留142、非正丢失384、正分未选丢失72、新增97，598→239；未知保留183、非正丢失251、正分未选丢失13、新增26，447→209。标签0／−1仍各自保留，不变成生物学错误或新基因。

这支持下一独立机制优先处理链对空集的评分／学习位置，而不只扩大A3候选搜索或改最终span DP：上述232条正确链几何已经存在，却在固定零操作点被A37拒绝。它不证明null、nonadditive、长度、训练域或类别不平衡是单一原因，也不授权调低阈值、改权重、加训练或救A37。具体新机制与数据范围仍须另行固定，完整raw-FASTA收益、匹配机制控制及独立泛化证据仍缺。
