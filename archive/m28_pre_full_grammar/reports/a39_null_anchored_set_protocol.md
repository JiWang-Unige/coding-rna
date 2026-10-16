# 绝对空集锚与兼容集合联合目标

2026-10-07北京时间，在A39首次新梯度和拟合前固定。A38已完成数值核验但五条件2／5，科学FAIL关闭；本轮是同一联合上下文架构上的独立fresh目标实验，不是续训或救A38。目的仍是完整raw-FASTA基因注释模型和论文，TRAIN准备结果不代替该总目标。

实际整链已闭合，工程与数值核验通过，原五条件仅2／5，A39科学FAIL关闭。不启动已准备的原生DEV，不改变目标系数、阈值、checkpoint或曝光；下文保留事前方案和按时间顺序的执行证据。

## 依据与唯一机制变化

[A38完整诊断](a38_logistic_context_protocol.md)显示，较强基线原已选正确丢失155条，其中98条新gain非正、57条正分未选；后57条有46条被known0阻挡、11条被unknown阻挡。全部98条正分未选正确配对中仅1条被多个成员阻挡，没有所有竞争者各自较弱但共同胜出的案例。逐候选BCE恢复了部分A37正确配对，却没有同时达到绝对保留和冲突排序目的。这是机制位置，不证明任何单一原因。

固定目标为(mean-known BCE＋mean-known compatible-set NLL)／2，直接组合未修改的core.chain_loss和chain_set.chain_set_loss。两项均除以同一窗原known数，不做正负组重配、标签重加权、margin、温度或系数扫描。BCE独立约束每条known链对空集的绝对gain；集合NLL通过原半开跨度约束的全窗sum-product分区约束可共存正确集合与known竞争集合，包含空集，不把多个合法共存基因当互斥单类。无known跨度冲突时两项相同，均值仍严格等价于原BCE尺度，而不是把梯度翻倍。

-1排除于两项监督和分区，独立gain直接梯度严格零，推理仍保留全部unknown。不因11条unknown竞争而改变标签或把它们当真负。空／全unknown给含null的可微零，原zero-known draw和AdamW过程保留；无gold窗仍监督全部原known0。known0仍是primary参考相对非匹配，不是生物学真负。

旧A13已有逐正链局部conflict／null logsumexp＋等正负组质量目标，但使用冻结上下文ChainHead165635参数和536参考注入，且没有这个同known归一的全窗兼容集合／逐候选复合目标。本轮不是重新命名A13，也不宣称相加损失、BCE、GRU或结构似然首创。A37、A38与A39匹配架构、输入、初始化配方和曝光；单seed探索性比较仍不足以证明一般性因果或论文创新。

## 保持不变的结构与曝光

fresh ContextChain612324参数＝SharedContext364256＋JunctionChainHead248068，feature2048／H128／seed0单随机流，FP32；与A37／A38同构造，不读取旧监督权重。GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345冻结，仅原四块6144bp特征和DNA，没有新GEN前向或生成器更新。

仍为原六条TRAIN染色体、1491唯一窗、原1536draw顺序与importance、原16321自然A3 greedy＋beam4候选及2380／12566／1375标签；无536注入、新搜索、新R4、DEV、独立物种、Setaria或A33访问。缺少注入gold的72窗、空／全unknown和zero-known draw不删。固定三遍4608 AdamW更新，lr3e-4／wd.01／betas(.9,.999)／eps1e-8／clip1，保存001536／003072／004608，仅末004608主读出，不挑seed或epoch。

score仍为exon additive＋junction additive＋nonadditive−window null，原gain0及span非重叠DP不变。原eligible2172／744＝2916、至少5-exon854／232＝1086不按预测缩小；候选／窗口配对不是独立基因或完整染色体F1。

## 检查和执行边界

当前先运行四项新目标合成检查：无冲突loss／gradient与BCE同尺度；known冲突时精确复合gain导数和unknown零直接梯度；强正确链冲突下绝对known0锚仍存在及输给known竞争者时排序梯度仍存在；空／全unknown／无gold语义。使用generanno／private0GPU／2CPU／8G／5分钟，独立M28-A39-NULL-ANCHORED-SET-CHECKS，只合成scalar，无研究模型、数据、权重、optimizer或GEN访问；不是拟合或科学PASS。

后续实际零更新preflight须先核实A38三根及13414584／13414782实际闭合，再复用原人口和八个metadata梯度样本规则，不按gain或性能挑窗。检查同一构造／参数未变、原输入无梯度、四个上下文及七个头分支适用梯度有限非零；没有known连接窗的连接两分支无梯度或精确零；实际gain梯度为两原目标gain梯度的等权均值，unknown严格零、optimizer0、无checkpoint。通过后仅拟合一次，CPU核对真实分支step，不假定所有参数4608。

拟合与自然评分沿用generanno／private1RTX3090／2CPU／16G，拟合1:45、自然30分钟，实际CPU核验0GPU／2CPU／8G／5分钟，no-requeue及独立A39输出／唯一日志。预算依据A37实际fit1:15:53、A38fit1:06:53与自然21:40留操作余量，不是新方法速度结论。实际preflight／fit／readout／CPU入口已接线；preflight和完整配对准备实际通过后直接提交一次研究拟合及其afterok下游，不提交原生DEV。

自然分数须完成全部1491窗／16321候选后才连接原保存标签；完整原DP、较强基线、三checkpoint／分支曝光、保存FP32组件和完整分母须由CPU核验，不冒充独立神经前向复现。

配对准备只读取既有metadata／JSON，不读取新checkpoint或特征tensor，核对A31_A3／A37／A38／A39占位四臂完整identity、候选顺序、protected标签、原DP和分母；A39占位为原较强基线，不能当新模型结果。两项新adapter合成检查覆盖四臂unknown及完整分母／五条件，缺失、重排或参考用于评分的比较臂必须拒绝。最终CPU除原四项证明外，核对真实训练保存的两分量、等权组合与原importance算术。

## 已完成的目标数学检查

实际13415040于2026-10-06 18:51:41 UTC结束，root／batch／extern均COMPLETED0:0，private gpu034／0GPU／2CPU／8G，5秒，batch MaxRSS468032K，Requeue0／Restarts0；独立输出STATUS COMPLETED／JOBID13415040。四项新合成检查全部PASS，0.163秒；无研究模型、输入、checkpoint、optimizer或GEN访问。该结果仅支持复合loss定义和零直接unknown梯度语义，不是实际分支准备、训练或科学PASS。

## 首个工程失败及唯一R1

实际原preflight13415110 root／batch FAILED1:0、16秒、batch MaxRSS492208K，原配对准备13415111 root／batch FAILED1:0、24秒、MaxRSS704036K；extern均COMPLETED0:0，独立原根STATUS FAILED exit=1／JOBID匹配，Requeue0／Restarts0。两个新四臂adapter检查已PASS0.003秒；随后两入口均在closed_a38同一条件失败，preflight尚未到population或new_model，配对准备尚未读取metadata新人口，未打开研究模型、特征tensor、新权重或optimizer，更新明确0、无checkpoint，原输出和日志不修改。

原因是A38竞争诊断只记录A31_A3／A38两臂，而自然读出和attrition记录A31_A3／A37／A38三臂；错误要求完整字典相等，虽实际共同两臂精确相同仍拒绝。唯一修复是逐物种精确比较共同A31_A3／A38两个完整counts子字典，未改变任何计数、配置、实际proof、分母、loss、阈值或科学门槛。不伪造缺失A37臂。独立PREFLIGHT-R1及READOUT-CHECKS-R1保留原相同资源／时间、no-requeue，使用同一已固定新模型与原8metadata窗口规则。这是本轮唯一允许的工程重跑；任何后续失败不再自动二次修复或拟合重试。

实际R1 preflight13415165 root／batch／extern均COMPLETED0:0，46秒，private gpu034／1RTX3090／2CPU／16G／5分钟，batch MaxRSS1542464K，Requeue0／Restarts0，STATUS COMPLETED／JOBID13415165／完整summary匹配。8个原metadata样本及完整1491窗／1536draw／16321自然候选闭合；四上下文和七头适用分支梯度有限非零，无known连接窗两个连接分支梯度零；实际gain导数为原两目标等权均值、unknown直接导数零、两物种known连接覆盖，参数未变、optimizer0、无checkpoint，原特征和DNA detached、没有新GEN或封存／DEV访问。

实际R1配对准备13415166 root／batch／extern均COMPLETED0:0，25秒，private gpu034／0GPU／2CPU／8G／5分钟，batch MaxRSS808016K，Requeue0／Restarts0，STATUS COMPLETED／JOBID13415166／summary匹配。两个新source回归检查PASS0.006秒，覆盖真实两臂vs三臂差异，并确认共同counts改变或实际CPU证明false仍拒绝；未重复原已PASS的两个adapter或四个数学检查。完整1491窗／16321候选／1536draw／1491特征index的四臂metadata、原protected标签／DP／较强基线及完整分母核对成立，没有打开新权重、特征tensor或研究模型，没有optimizer更新。准备通过不代替拟合或科学PASS。

## 实际研究执行链

已一次提交fit13415276→natural13415309→CPU13415319。fit实际2026-10-06 19:04:23 UTC开始，private gpu034／1RTX3090／2CPU／16G／1:45，no-requeue／Restarts0；SubmitLine准确包含afterok13415165:13415166和export A39_PREFLIGHT_JOB13415165。初次实际RUNNING0:19、后续2:17，STATUS RUNNING／JOBID13415276且training.jsonl已有64行、stderr为空，是部分训练证据而非4608完成或科学结果。

natural准确提交afterok13415276:13415166，export A39_FIT_JOB13415276／A39_CHECKS_JOB13415166，请求private1RTX3090／2CPU／16G／30分钟；Slurm移除已满足13415166后仅保留afterok13415276合法。CPU准确提交afterok13415309，export A39_FIT_JOB13415276／A39_READOUT_JOB13415309，请求private0GPU／2CPU／8G／5分钟。两者实际PENDING0:00、正常Dependency，尚未分配节点，不把请求资源写成实际下游运行。三个canonical命令／唯一日志／no-requeue／Restarts0均核对。

提交前实时五项调度检查和输出碰撞检查通过，157TiB可用空间覆盖本轮小额新产物；未修改、取消或解释任何旧无效依赖及其他项目作业。原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3与同heartbeat coding-rna-r5-light-monitor只接管本组三根，15分钟、既有gpt-5.6-luna／low；原A38闭合及所有短准备排除，不创建新监控或改变主模型。配置／单次启动／本轮定时执行／整链交接分别按实际证据记录，尚未观察到本轮定时运行或终态。

实际同heartbeat于2026-10-06 19:07:37 UTC切换ACTIVE，保存完整新提示、ID／name／kind heartbeat／固定目标对话／15分钟／既有通知字段／created_at1790981930511精确核对。A39首次人工启动2026-10-06 19:08:32 UTC检查fit RUNNING gpu034和两个正常依赖；同一对话turn01a1129d-d417-7a83-b46d-9e4f627047e9完成、error-null、27.087秒，实际turn_context19:08:15.817Z为gpt-5.6-luna／low／原cwd。没有向主对话发送终态交接；不是定时运行或完成证明。监控遵循固定安静等待，整链完成或可行动故障后向原主一次接续并暂停，主负责原五条件判定及后续已授权工作。

## 原继续条件与停止规则

较强A31_A3五条件全部保持：拟南芥正确严格>1531，水稻正确严格>478；known_nonmatch各≤347／251；合并至少5-exon正确严格>700。A37／A38配对改善不替代此门槛；任一FAIL关闭A39，不改目标系数、阈值、权重、checkpoint或曝光。PASS只支持另定完整raw-FASTA执行，不等于SOTA、独立泛化或论文成功。

首个工程失败保留，只有明确optimizer0时才允许一次同科学／同资源最小实现修复到独立R1；该唯一修复已用于准备阶段。TIMEOUT、已有更新或不明更新不自动续训或重提。原A38／A37／A36／A13关闭、A33例外未部署、Setaria封存和公开发布边界不变。当前已提交的研究三根仅由原固定对话及同heartbeat监控；短准备由主消费，不增加旧对象、新监控或重复等待循环。

原生raw-FASTA接线和完整评价按[a39_null_anchored_set_dev_protocol.md](a39_null_anchored_set_dev_protocol.md)条件准备，尚不执行真实DEV。只有本轮原三根闭合、五项真实CPU证明和原五条件全部PASS才继续；准备检查不改变TRAIN协议或授权封存数据。

## 完整结果的配对归因接线

后续结果处理复用未修改的diagnose_structured_emission_attrition.transitions；分别记录A31_A3、A37、A38到A39的全部原自然TRAIN配对。只在fit13415276／natural13415309／CPU13415319实际COMPLETED、markers和模型／目标身份相符、五项实际CPU证明及原五条件算术一致后读取分数。该描述性处理允许分析实际科学PASS或FAIL，不以FAIL拒绝诊断，也不以诊断代替原生DEV放行条件。

全1491窗／16321候选五路记录和四臂原DP数量必须重合，原2916／1086分母保留。逐物种与pooled按known1、known0、protected−1及known正确至少5-exon，报告保留、新增、丢失且gain非正、丢失且正分未选、两者未选且非正／正分六类。净变化由新增减两类丢失解释；window／candidate配对不是独立基因，不把known0或unknown解释成生物学真负，不做阈值、系数、标签、checkpoint或曝光救援。

两个新增人工接线测试覆盖四臂人口／原DP／unknown及归因净量，以及关闭的科学FAIL仍可分析但缺少实际复合算术证明必须拒绝。测试只在后续同一归因作业中先执行，不重复原统计函数已有测试；准备时尚未运行，实际执行见下节。归因使用generanno／private0GPU／2CPU／8G／5分钟／no-requeue，独立M28-A39-NULL-ANCHORED-SET-ATTRITION及唯一日志，复用A38对应实际21秒并留四臂统计余量，输出预计低于1MiB。只读取已闭合JSON metadata／scores／labels，不打开模型、feature tensor、TRAIN catalog或参考，不作GEN／神经前向、optimizer、新搜索、R4、DEV、独立物种、Setaria或A33访问。失败保留日志并停止使用该归因结果，不自动修复或重提，也不把归因工程失败改判原A39科学结果。该短后处理由主消费，不改当前三研究作业或固定heartbeat对象。

## 实际完整结果及关闭

fit13415276、natural13415309、CPU13415319的root／batch／extern均COMPLETED0:0，实际耗时1:04:31／13:17／0:29，全部gpu034。fit与natural实际private1RTX3090／2CPU／16G，CPU实际0GPU／2CPU／8G；batch MaxRSS分别2213492K／1660900K／879452K，三根STATUS／JOBID与summary对应，natural和CPU stderr为空。fit固定4608更新、三份末遍checkpoint，training.jsonl实际4608行；natural两份JSONL各1491行，全部16321原自然候选完成后才连接保存标签，没有新候选、参考注入、GEN／R4、DEV、Setaria或A33访问。

CPU五项实际证明均true：原draw／importance／counts和4608更新，原两损失分量及等权组合算术，三checkpoint实际FP32 tensor和分支AdamW曝光，自然身份与保存FP32 score分量，以及原标签／DP／基线／完整分母和五条件重算。SharedContext及null的实际step为4608，其余非连接头4602、两连接分支4425。主由完整natural与CPU summary重新核对配置、来源、计数、分母和五个严格整数条件，结果一致；这不是独立神经前向复现。

| 原继续条件 | A31_A3 | A38 | A39 | A39结果 |
| --- | ---: | ---: | ---: | --- |
| 拟南芥已选正确严格增加 | 1531 | 1494 | 1496 | FAIL，较基线−35 |
| 拟南芥已选known_nonmatch不增加 | 347 | 127 | 131 | PASS |
| 水稻已选正确严格增加 | 478 | 477 | 478 | FAIL，持平不等于严格增加 |
| 水稻已选known_nonmatch不增加 | 251 | 101 | 101 | PASS |
| 合并已选至少5-exon正确严格增加 | 700 | 665 | 653 | FAIL，较基线−47 |

完整eligible reference/window分母仍为2172／744＝2916，至少5-exon854／232＝1086。A39拟南芥总已选1805＝1496正确＋131非匹配＋178未知，水稻626＝478＋101＋47；合并2431＝1974＋232＋225。相对A31_A3正确−35、非匹配−366、未知−222；相对A38正确仅＋3、非匹配＋4、未知＋8、至少5-exon正确−12。候选／窗口配对不是去重基因，known0不是生物学真负，未知不计作错误或新基因。

固定watcher于2026-10-06T20:23:25Z取得整链终态和结构白名单，向原主一次交接，记录事件并原生暂停同heartbeat；主读取实际指标后关闭A39。当前保存的kind heartbeat、固定目标对话、15分钟及PAUSED均核对，没有新监控或模型变化。两个A39 DEV输出根实际不存在；这是未执行DEV，不是失败的DEV或新gene F1结果。

## 实际同池配对归因

短归因13417784实际root／batch／extern均COMPLETED0:0，2026-10-06 20:32:55–20:33:25 UTC，30秒，private gpu034／0GPU／2CPU／8G／5分钟，batch MaxRSS807160K、Requeue0／Restarts0。两个新人工adapter检查PASS0.002秒；随后完整1491窗／16321候选、四臂原DP计数和2916／1086参考分母与已闭合结果一致。主核对两物种及pooled的三个对照、四类共36组：六种互斥状态之和＝pool，旧已选＝保留＋两类丢失，新已选＝保留＋新增，净变化一致。实际STATUS／JOBID／summary对应；stderr仅为两项通过的测试日志，没有异常。源码、测试、SBATCH和事前追加协议与canonical逐字一致；该短作业由主消费，同心跳维持PAUSED。

| 正确配对范围及对照 | 旧已选 | 保留 | 丢失gain≤0 | 丢失正分未选 | 新增 | A39已选 | 净变化 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A31_A3→A39 拟南芥 | 1531 | 1425 | 66 | 40 | 71 | 1496 | −35 |
| A31_A3→A39 水稻 | 478 | 425 | 31 | 22 | 53 | 478 | 0 |
| A31_A3→A39 合并 | 2009 | 1850 | 97 | 62 | 124 | 1974 | −35 |
| A31_A3→A39 合并至少5-exon | 700 | 620 | 44 | 36 | 33 | 653 | −47 |
| A38→A39 合并 | 1971 | 1883 | 45 | 43 | 91 | 1974 | ＋3 |
| A38→A39 合并至少5-exon | 665 | 621 | 19 | 25 | 32 | 653 | −12 |

相对较强基线，丢失159条原已选正确，其中97条非正、62条正分未选；新增124条不能抵消。长链丢失80条，44条非正、36条正分未选，新增33条。双方均未选正确还包括185条非正及62条正分；完整正确自然池2380＝1850＋97＋62＋124＋185＋62，长链自然池788不能替代1086参考分母。known0保留143、非正丢失396、正分丢失59、新增89；未知保留184、非正丢失254、正分丢失9、新增41，均与原净计数一致。

固定等权联合目标在本次单seed下没有消除两类正确损失，相较BCE还减少长链选择；不能以较少非匹配、局部新增或TRAIN loss下降替代原门槛。三种目标的同构造比较不证明某一组件的普遍因果或联合学习整体无效。后续独立模型应检验不同的信息表示或结构机制，而不继续扫描这三个已关闭目标的混合系数、零点或曝光；仍须真实raw-FASTA、独立泛化及机制控制证据。

权威产物为outputs/M28-A39-NULL-ANCHORED-SET-{FIT,TRAIN-READOUT,RESULT-CHECK,ATTRITION}，原checkpoint、4608训练行、两份1491行JSONL及日志保留。归因没有模型／feature tensor、catalog或参考访问，没有前向、optimizer、新搜索、新R4或封存数据访问；A37／A38／A39均保持科学FAIL关闭，论文总目标仍未达到。
