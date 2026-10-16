# 联合上下文整链模型的逐候选空集目标

2026-10-07北京时间，在首次新梯度与拟合前固定。A38检验同一联合上下文整链模型中，逐候选的绝对空集约束能否比兼容集合似然更好地保留正确链，同时控制已知非匹配。A37科学FAIL关闭，不续训、不救参；本轮使用fresh参数和独立输出。科研目标仍为raw-FASTA基因注释模型与最终论文，不放行独立物种、Setaria或公开发布。

实际结果已闭合：工程和数值核验通过，原五条件仅2／5，A38科学FAIL关闭。较A37改善不替代较强A31_A3门槛；不执行已准备的原生DEV，不救阈值、checkpoint或训练曝光。后文保留原方案及实际结果。

## 依据和待检验机制

[A37已闭合结果](a37_context_chain_result.md)只通过五条件中的两项。相对同池A31_A3，丢失279条正确窗口—链配对，其中232条新gain≤0、47条正分未选，新增97条不能抵消；这是阶段定位，不证明集合目标或null是单一原因。

现有集合目标对非匹配的gain导数是其兼容集合边际概率除以known数；与高分正确链冲突的非匹配，即使绝对gain>0，也可能只有很小惩罚。逐候选BCE则独立锚定每条known链对空集，正链的梯度为(sigmoid(g)−1)/|K|，非匹配为sigmoid(g)/|K|。这种变化同时影响两类梯度，可能进一步牺牲召回，不能预先声称能恢复上述232条或改善校准。

A6已用过BCE，但其上下文被冻结、头无显式连接且包含参考注入；A13另有冲突竞争和等正负组质量；A14、A15及A21用集合似然；A36是固定局部角色累积。A38检验的是A37联合上下文与有序exon／junction头上的逐候选空集目标，不把BCE或现有GRU重新命名为创新。它与A37只改变目标，模型容量、输入候选、初始化配方和优化曝光匹配；一次seed的结果仍不足以证明一般性因果或论文贡献。

## 唯一变化和冻结曝光

直接复用既有core.chain_loss，不传候选权重：K仅含原自然0／1，loss为所有K候选的BCE均值；-1完全排除，独立gain的直接梯度严格零，推理继续保留未知。全未知／空候选给可微零，空gold窗口仍监督原known非匹配，保留原draw及AdamW过程。不是正负等组质量、类别加权、辅助loss、margin、温度或新标签。推理仍为原gain0及span非重叠DP，不改全窗多基因语义。

保持A37 ContextChain612324参数＝SharedContext364256＋JunctionChainHead248068，feature2048／H128；直接使用相同new_model与seed0单随机流构造，不读取任何旧训练权重。原gain为exon additive＋junction additive＋nonadditive−window null，FP32。GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345冻结，只读原四块6144bp特征和DNA，不做新GEN前向。

仍为原六条TRAIN染色体、1491唯一窗、1536draw的顺序／importance权重、16321自然A3 greedy＋beam4候选，标签1／0／-1为2380／12566／1375；不使用536注入。失去注入gold的72窗、空／全未知和zero-known draw不删除。固定fresh seed0、三遍4608 AdamW更新，lr3e-4／wd.01／betas(.9,.999)／eps1e-8／clip1；保存001536／003072／004608，仅末004608主读出，不挑epoch、增加seed或追加曝光。

## 核验和执行

四项人工检查只覆盖新目标的科学语义：known均值梯度及unknown保护；正确链与可共存非匹配集合冲突时，BCE绝对null惩罚仍存在；BCE到联合上下文和七个头分支的梯度及未知不改变known梯度；空／全未知／无gold语义。不重复未变的模型索引、集合DP或raw调用测试。

零更新preflight要求A37 fit13410061、自然13410739、CPU13410742及配对诊断13412326实际闭合，再复用原1491人口核验和八窗metadata选择规则，不按新分数挑窗。四个上下文及七个头分支须有限非零梯度；没有known连接时两连接分支须无梯度或精确零。实际gain梯度须符合均值BCE公式；原输入无梯度、参数前后相同、optimizer0、无checkpoint。通过后只拟合一次；分支实际AdamW step按原候选调用结构核验，不假定全部为4608。

generanno环境，private1RTX3090／2CPU／16GiB，preflight5分钟、fit1小时45分钟、自然读出25分钟；CPU核验private0GPU／2CPU／8GiB／5分钟。依据A37实际fit1:15:53、自然0:16:59留运行余量，时窗不是方法吞吐或完成保证。输出独立M28-A38-LOGISTIC-CONTEXT-{PREFLIGHT,FIT,TRAIN-READOUT,RESULT-CHECK}，三checkpoint及JSON预计<64MiB，不创建新特征缓存。提交前据活调度／维护／存储路由，no-requeue、原排除节点、唯一日志。

拟合运行期间并行进行一次CPU三臂接线核验，0GPU／2CPU／8GiB／5分钟，独立READOUT-CHECKS输出。两个人工检查覆盖三臂、未知和完整分母及不以较弱A37替代门槛、缺失或乱序的对照拒绝；随后重放全部原人口metadata，A38槽用原基线作人工占位，不打开新权重、特征tensor或前向，不产生A38性能。核验完成是自然读出的afterok依赖。

自然读出先完成全部1491窗／16321候选分数，再连接原保存标签；逐物种记录A31_A3、A37、A38同几何选择、unknown及≥5-exon正确配对。CPU重算原draw／importance／4608行、三checkpoint和实际分支曝光、保存FP32组件、三臂原DP／标签／reference-window分母和五条件；不做独立神经前向。原eligible2172／744＝2916、≥5-exon854／232＝1086，不按新预测缩小。candidate/window数量不是独立基因或染色体F1。

## 原继续条件和停止规则

门槛仍为较强A31_A3：两物种已选正确各严格高于1531／478、known_nonmatch各不高于347／251，pooled≥5-exon正确严格高于700。A37配对对比不替代此门槛；unknown另列，0只是primary参考相对非匹配，不当生物学真负或新基因。

任一FAIL关闭A38，不改阈值／loss／权重／checkpoint／曝光。PASS只支持另外固定真实raw-FASTA执行方案，不等于SOTA、独立泛化或论文成功，也不能恢复池外缺失几何。首个工程失败保留原状态；仅明确证实optimizer0时允许至多一次同科学、同资源最小修复到独立R1目录。已有更新、TIMEOUT或不明更新不自动续训或重提。

长期等待只复用固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同heartbeat coding-rna-r5-light-monitor，原gpt-5.6-luna／low／15分钟，静默正常状态；整链闭合或可行动故障一次交接。A33例外、独立面板、DEV及Setaria封存和公开发布边界保持不变。

## 实际执行收据

科学目标、结构、数据域、优化曝光和五条件均在13412601首次新梯度前固定；三臂接线核验的实现于fit已运行、自然读出尚未执行时补齐，没有改变上述科学条件。原A37/A36结果、checkpoint和源码未修改。

| 阶段 | 实际作业与状态 | 资源和依赖 |
| --- | --- | --- |
| 零更新preflight | 13412601，root／batch／extern COMPLETED0:0，37秒 | private1RTX3090／2CPU／16G／5分钟，gpu034 |
| 固定fit | 13412630，root／batch／extern COMPLETED0:0，1:06:53，gpu034 | private1RTX3090／2CPU／16G／1:45，export A38_PREFLIGHT_JOB=13412601 |
| 并行三臂接线核验 | 13412682，root／batch／extern COMPLETED0:0，31秒 | private0GPU／2CPU／8G／5分钟，gpu034 |
| 全自然候选评分 | 13412701，root／batch／extern COMPLETED0:0，21:40 | private1RTX3090／2CPU／16G／25分钟，afterok:13412630:13412682，export A38_FIT_JOB=13412630,A38_CHECKS_JOB=13412682 |
| 实际结果CPU核验 | 13412702，root／batch／extern COMPLETED0:0，47秒 | private0GPU／2CPU／8G／5分钟，afterok:13412701，export A38_FIT_JOB=13412630,A38_READOUT_JOB=13412701 |

preflight四项新目标检查PASS／0.921秒；八个原metadata窗口的全部适用分支有限非零、无known连接窗的连接两分支严格零，实际均值BCE gain梯度符合公式，unknown直接梯度严格零、参数不变、optimizer0。原1491／1536／16321、原三类计数及2916／1086分母一致，GPU峰值165758464字节，八窗内部load／forward／backward共约11.1秒。batch MaxRSS1838972K、STATUS／JOBID／summary相符。这是接线准备，不是训练成功或性能改善。

固定fit的最终工程收据已收取：13412630实际root／batch／extern COMPLETED0:0，耗时1:06:53，batch MaxRSS2009612K；STATUS／JOBID／summary相符。training.jsonl实际4608行，summary记录三遍各1536步；001536／003072／004608三份checkpoint各7389798字节，仅末004608读出。配置metadata与上述冻结结构、原16321自然候选、mean-known BCE、fresh seed0和GEN revision一致，没有参考注入、DEV／Setaria选择或A33访问记录。首次收取只读完成metadata、文件大小和行数，未加载tensor／optimizer或自然分数；不能由summary或4608行直接断言实际参数曝光。后续13412702已完成该实际核验，见下节；拟合工程完成本身仍不代表性能改善或五条件PASS。

并行CPU核验两项新三臂检查PASS／0.002秒；全部1491窗／16321候选、1536draw、1491feature-index、原标签／三臂DP配对／较强基线／2916和1086分母相符；没有新权重、特征tensor、模型前向或optimizer访问。batch MaxRSS772728K，STATUS／JOBID／summary相符，已由主消费，不纳入长期监控。

提交前已完成活调度、预约、唯一输出与存储检查；private空闲RTX3090满足实际24GiB，原并发低于80，没有相关维护，实际156TiB余量覆盖<64MiB新产物。只提交每个当前作业一次，Requeue0／Restarts0；exact SubmitLine／export／canonical命令及fit实际allocation相符。下游请求资源不冒称已实际分配；Slurm移除已完成13412682依赖，仅留afterok:13412630合法。后续整链首次snapshot为fitRUNNING4:58／gpu034、两下游PENDING0:00／valid Dependency。

已更新同一原生heartbeat到上述三root；2026-10-06T16:55:05Z实际完整prompt／ACTIVE／heartbeat／原固定target／15分钟／failed_runs_only准确读回，created_at保留。手动启动核对实际16:55:26Z为fit RUNNING／gpu034、自然与CPU正常依赖等待；对应turn01a11223-e840-7792-862e-bed116a0b11d completed／error-null／28.951秒，实际turn_context仍gpt-5.6-luna／low／原目录。没有新监控对话、自动化、模型变化或正常状态主交接。当前cohort已观察真实定时追加：heartbeat触发17:10:14.331Z，同对话turn01a11231-c63d-7b23-80ea-75460504630d completed／error-null／18.492秒，实际turn_context仍gpt-5.6-luna／low；17:10:26Z的最小调度核对仍为fit RUNNING及两正常依赖，DONT_NOTIFY。主于18:20:58Z直接核实整链已完成并消费科学结果，没有自动终态消息；同心跳随后原生PAUSED，完整prompt／ID／kind／固定target／15分钟／通知字段／created_at准确保留并读回，旧cohort不再检查或重复交接。最初12个源码／测试／SBATCH与canonical逐字一致；完整TRAIN自然读出及五条件已闭合，原生DEV未执行。

## 实际科学结果及关闭

完整1491窗／16321自然候选已评分，natural_scores.jsonl与windows.jsonl各实际1491行；所有分数关闭后才连接原保存标签，没有参考注入、新搜索、新GEN、参考catalog、DEV、Setaria或A33访问，没有新optimizer。natural及CPU的STATUS／JOBID／summary相符，batch MaxRSS分别1606940K／856376K。

CPU四项实际证明均true：原draw顺序／importance／counts及4608更新；三checkpoint有限FP32张量和实际分支AdamW曝光；原自然候选身份与保存FP32分量算术；原标签／三臂DP／较强基线／完整reference-window分母及五条件。SharedContext和null实际参数step4608，exon／链GRU／additive／nonadditive为4602，junction两分支为4425，与原候选调用曝光一致，并非所有头参数均4608。CPU没有新模型或backbone前向；这是保存组件和决策核验，不冒充独立神经前向复现。

| 原继续条件 | A31_A3 | A37 | A38 | A38结果 |
| --- | ---: | ---: | ---: | --- |
| 拟南芥正确配对严格增加 | 1531 | 1396 | 1494 | FAIL，较基线−37 |
| 拟南芥known_nonmatch不增加 | 347 | 125 | 127 | PASS |
| 水稻正确配对严格增加 | 478 | 431 | 477 | FAIL，较基线−1 |
| 水稻known_nonmatch不增加 | 251 | 114 | 101 | PASS |
| 合并至少5-exon正确配对严格增加 | 700 | 639 | 665 | FAIL，较基线−35 |

原eligible分母仍为2172／744＝2916，至少5-exon为854／232＝1086。拟南芥总已选1792＝1494正确＋127已知非匹配＋171未知；水稻624＝477＋101＋46。合并正确1971、已知非匹配228、未知217、总已选2416，较基线分别−38／−370／−230／−638；较A37正确＋144、已知非匹配−11、未知＋8、总已选＋141、长链＋26。

主由完整自然读出和CPU summary重新计算五项严格整数条件、原分母及合并计数，与CPU结论逐项一致：2／5，A38科学FAIL关闭。较A37接近基线不能替代门槛，水稻只差1也不放宽严格增加条件。该一次seed的TRAIN候选配对结果不等于完整gene F1、机制因果或独立泛化；不挑checkpoint、不追加曝光、不改变gain0／loss／权重。已准备的DEV接口和执行方案保留但不运行，两个DEV输出根实际不存在。

## 完成评分后的三臂转移诊断

复用未改动的diagnose_structured_emission_attrition.transitions，分别统计A31_A3→A38与A37→A38。六种转移区分保留、新增、gain≤0丢失、正分未选丢失及两种未选状态；逐物种和pooled保留1／0／-1与≥5-exon正确子集。较弱A37上的新增不能替代较强A31_A3继续条件；这些保存评分的描述性配对不证明目标的单一因果、独立基因F1或泛化。

新诊断脚本只在fit13412630、自然评分13412701与CPU13412702实际完成，四项CPU证明成立且原五条件算术一致后读取保存JSON；科学PASS与FAIL均可诊断。完整1491窗／16321候选的三臂计数须重现CPU结果、原基线和已闭合A37，reference-window分母仍为2916／1086。它不打开模型或特征tensor，不做前向、optimizer、标签重建、参考catalog、DEV、Setaria或A33访问，也不修改gain0、checkpoint、候选或停止规则。

接线核验13413133实际root／batch／extern COMPLETED0:0，19秒；private0GPU／2CPU／8G／5分钟，gpu034，batch MaxRSS720608K，STATUS／JOBID相符，Requeue0／Restarts0。五项人工检查PASS／0.006秒，覆盖新双对照统计、unknown与完整分母、较强基线条件、PASS／FAIL诊断与陈旧gate／source拒绝、不完整CPU证明在读分前阻断，以及新四路配对缺失／乱序。四个新增诊断源码／测试／SBATCH与canonical逐字一致；它们不修改原训练、评分或CPU核验源码。

整链闭合后仅提交一次sbatch --parsable sbatch/M28-A38-LOGISTIC-CONTEXT-ATTRITION.sbatch，实际作业13414584。root／batch／extern均COMPLETED0:0，21秒，gpu034，实际private0GPU／2CPU／8G／5分钟／no-requeue，batch MaxRSS785468K、Requeue0／Restarts0；STATUS／JOBID／summary及canonical命令／SubmitLine相符，stderr为空。提交前五项活调度／预约检查无相关维护、并发限制或输出冲突，新诊断源码和SBATCH与canonical逐字一致，13413133完成标记及源作业相符。该短诊断由主消费，不纳入已暂停的旧整链心跳，也不改变A38关闭状态。

诊断的全部1491窗／16321自然配对重现实际三臂计数及2916／1086完整参考分母。主核对两个物种与pooled、两个对照、1／0／-1及长链子集共24组，均满足六种互斥转移之和＝pool、旧已选＝保留＋两类丢失、新已选＝保留＋新增。粒度仍是候选／窗口配对，不是新基因或去重基因。

| 正确配对范围与对照 | 旧已选 | 保留 | 丢失gain≤0 | 丢失正分未选 | 新增 | A38已选 | 净变化 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A31_A3→A38 拟南芥 | 1531 | 1431 | 69 | 31 | 63 | 1494 | −37 |
| A31_A3→A38 水稻 | 478 | 423 | 29 | 26 | 54 | 477 | −1 |
| A31_A3→A38 合并 | 2009 | 1854 | 98 | 57 | 117 | 1971 | −38 |
| A31_A3→A38 合并至少5-exon | 700 | 628 | 42 | 30 | 37 | 665 | −35 |
| A37→A38 合并 | 1827 | 1719 | 56 | 52 | 252 | 1971 | ＋144 |
| A37→A38 合并至少5-exon | 639 | 595 | 21 | 23 | 70 | 665 | ＋26 |

较强基线原已选正确共丢失155条，98／155＝63.2%为A38 gain非正，57条仍正但未选；新增117条未抵消。至少5-exon丢失72条，42条非正、30条正分未选，新增37条未抵消。双方均未选的正确配对还包括213条非正、41条正分未选，完整自然池2380＝1854＋98＋57＋117＋213＋41；长链自然池788不能代替1086参考分母。

较基线的known_nonmatch保留151、非正丢失394、正分丢失53、新增77，598→228；未知保留186、非正丢失254、正分丢失7、新增31，447→217。较A37新增252条正确配对，同时丢失108条，净增144；不是新搜索产生的几何，也不证明已恢复A37的某个特定丢失子集。

下一机制不能继续只把损失定位为非正拒绝：剩余155条基线正确损失中有57条正分未选，长链中有30／72；两种位置都需要解释。与A37匹配初始化配方、容量和曝光的单seed对比提供探索性目标差异证据，但不能证明null、BCE或竞争头的单一因果，也不能用调低阈值、挑checkpoint或混合旧分数救A38。下一独立模型需另行固定机制和停止条件；完整raw-FASTA、独立泛化和论文主张仍未成立。

## 正分未选链的实际集合竞争诊断

闭合后的下一问题是57条基线曾选正确链为何在A38仍为正分却未选，其中30条至少5-exon。本次同时统计全部98条自然正确正分未选链及其中48条长链，避免把旧基线曾选的子集当成完整漏选人口。只能读取已闭合13412630／13412701／13412702／13414584和原基线13398521的完成metadata及保存JSON，不读取权重、特征、原始DNA、catalog、新标签或DEV。

复用未改动的audit_chain_displacement.exchange，按原严格半开整链跨度找到实际选中冲突成员，记录其1／0／-1标签组合、数量和总gain；逐链及每窗所有正分未选正确链的兼容联合集合均核对插入／移除分数。原集合替换舍入余量仍为1e-6×max(1,两总量绝对值)，不是放宽A38保存FP32组件核验或改变选择阈值。统计多成员阻挡及“每个单独竞争者都严格弱于正确链，但集合共同胜出”的数量。它只解释实际选中集合的可行替换，不计算强制插入后的最佳集合optimum，不证明某个头或目标的单一因果；unknown不改为负类。

新三路adapter须复现每窗A31_A3／A38的原DP选择和完整计数；全1491窗／16321候选、2916／1086参考分母，以及两物种和pooled的全部正分／基线丢失／长链四个分层须重现已闭合转移诊断。输出独立M28-A38-LOGISTIC-CONTEXT-COMPETITION，保存1491行windows.jsonl、summary、STATUS和JOBID；不改变原分数、候选、gain0或A38关闭结论。

四项新adapter检查覆盖混合known0／unknown两条较弱链的实际集合阻挡、基线丢失与原未选子集分开及零gain不入正分组、缺失／乱序或陈旧选择拒绝、缺少CPU证明时在读分前停止；不重跑旧A22纯exchange检查。入口M28-A38-LOGISTIC-CONTEXT-COMPETITION.sbatch使用generanno／private0GPU／2CPU／8G／5分钟／no-requeue及唯一日志，先新检查再实际诊断。依据13414584完整四路人口诊断21秒留操作余量，只写小型JSON；不建新监控或重启已暂停旧心跳。首个失败保留，仅允许一次同科学／同资源、明确实现问题的最小修复到独立R1输出，不据诊断结果救参或启动DEV。

提交前活调度／预约五项检查无相关维护、并发或输出冲突；只同步新脚本、新adapter检查、SBATCH和本记录四个明确文件，逐字与canonical一致。仅提交一次sbatch --parsable sbatch/M28-A38-LOGISTIC-CONTEXT-COMPETITION.sbatch，实际13414782；root／batch／extern均COMPLETED0:0，47秒，gpu034，actual AllocTRES为0GPU／2CPU／8G，5分钟时限、Requeue0／Restarts0、canonical命令／SubmitLine相符，batch MaxRSS769180K。四项新adapter检查PASS／0.004秒；STATUS／JOBID／summary相符，windows.jsonl实际1491行，全16321候选、原计数及2916／1086分母一致。全部逐链和每窗联合集合替换核对通过，不存在本轮选择器算错该可行替换的证据。

| 正分未选正确配对范围 | 合计 | 竞争者仅known0 | 竞争者仅unknown | 被至少2条链阻挡 |
| --- | ---: | ---: | ---: | ---: |
| 拟南芥基线曾选子集 | 31 | 23 | 8 | 0 |
| 水稻基线曾选子集 | 26 | 23 | 3 | 0 |
| 合并基线曾选子集 | 57 | 46 | 11 | 0 |
| 合并基线曾选且至少5-exon | 30 | 24 | 6 | 0 |
| 合并全部正分未选 | 98 | 65 | 33 | 1 |
| 合并全部正分未选且至少5-exon | 48 | 30 | 18 | 1 |

98条中只有1条水稻长链被两个known0成员共同阻挡，未落在基线曾选子集；实际blocker成员99次＝66次known0＋33次unknown，不能把配对次数当去重竞争基因。没有label1或混合0／-1竞争组合，也没有“所有单个竞争者均严格弱于正确链而集合共同胜出”的实际案例。这不否定集合选择语义，只说明本批损失不是主要由该特殊多成员机制造成。

下一独立实验优先考虑同时约束链对空集的绝对gain与已知冲突链的相对排序，而不是修选择器或仅扩大束宽。A37的集合似然与A38的BCE各自未满足五条件；旧A13已有局部冲突／null及等正负组质量目标，不能重新命名为新方法。任何新结合须与其区分，另行冻结相同容量、自然候选、初始化配方和曝光；仍不把unknown当负例，也不以11条unknown竞争为由改标签。以上只定位下一检验，不证明新目标必定改善、单组件因果、完整raw-FASTA收益或论文创新。A38继续关闭，旧心跳保持PAUSED。
