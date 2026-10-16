# 候选角色证据进入有序整链表示

2026-10-07北京时间，在首次新梯度前固定。A37／A38／A39均科学FAIL关闭；A39相对较强基线丢失159条原已选正确配对，其中97条非正、62条正分未选，新增124条未抵消，长链净减47。本轮检验不同的输入表示，不扫描已关闭损失的系数、零点、checkpoint或曝光。目标仍为完整raw-FASTA模型、独立适用范围和论文。

## 新机制及旧实验差异

冻结A3末004608的SharedContext和15类segmentation分支，原GEN b0483c23b6b63787b61a6d3a204a9b517d6ba345不更新。根据每条候选的原ascending phase、半开exon／intron、donor／acceptor及start／stop优先级，计算逐base的logit(role)−logit(intergenic0)。角色只由自然候选几何决定，不由参考或1／0／−1标签决定；沿用candidate_roles，不新增生物学语法。

新头在每个exon token原mean／first／last隐藏表示和三项几何后，追加相应三项角色margin：exon内平均、首base、末base。junction token在原五个隐藏位置和两项几何后，追加intron平均、前exon末base、intron首末base、后exon首base的五项角色margin；零长度gap三项intron值精确0，两个exon端点保留。输入先经非线性exon／junction投影，再交错进入同一链GRU及原additive／nonadditive分支。窗口null仍仅读取新上下文，没有标签过滤或后处理。

A26把14状态mean／count压成28维线性标量后加到末端nonadditive；本轮不复用该权重、缓存或输出项，而让候选角色条件化每个局部token和有序聚合。A31是冻结单位局部势的span求和；本轮没有把A31分数与旧A39分数混合。A36重新学习15类发射并按单位路径求和；本轮冻结既有局部预测，只拟合fresh链表示及新上下文。A35的旧角色贡献和A39两类损失提供检验位置，不证明该表示必定改善或存在一般性因果。

fresh seed0单流按SharedContext再扩展头的原七分支顺序构造，613348参数＝364256上下文＋249092头；相对A39增加1024个输入权重，没有新独立评分头。原共有输入维度变化会改变后续随机流，不冒称初始tensor逐值匹配。学生不读取任何旧监督权重；冻结A3源649457参数及其既有监督曝光是显式输入，不称全模型从零或匹配容量／曝光的因果消融。局部源无梯度、无optimizer更新、不解码query、不更新生成器或GEN；尚无组件新颖性或论文性能主张。

## 数据目标与停止条件

仍为原六条TRAIN染色体、1491唯一窗、1536draw原顺序和importance、原16321自然A3 greedy＋beam4链，标签2380／12566／1375。无536参考注入、新搜索、新R4、DEV、独立物种、Setaria或A33访问；缺失gold、空／全unknown及zero-known draw保留。原四块6144bp／feature2048和DNA不变，FP32；每窗读取冻结A3本地分支和新学生上下文，复用同一GEN缓存，不建立全量新特征缓存。

损失原样使用(mean-known BCE＋compatible-known-set NLL)／2及同一known归一，−1不进两项监督或集合partition，独立gain直接梯度0，推理全部候选保留。原score＝exon additive＋junction additive＋chain nonadditive−window null，gain0与原span DP不变。固定fresh seed0、三遍4608 AdamW更新，lr3e-4／wd.01／betas(.9,.999)／eps1e-8／clip1，保存001536／003072／004608，只读末004608，不选择epoch或增加seed。

较强A31_A3五条件不变：拟南芥正确严格>1531、水稻严格>478，两物种known_nonmatch各≤347／251，pooled至少5-exon正确严格>700。eligible reference/window仍为2172／744＝2916，至少5-exon854／232＝1086；候选／窗口配对不是独立gene F1，known0不是生物学真负，unknown另列。任一FAIL关闭A40，不救阈值、目标、权重或曝光；PASS也只支持另行固定完整raw-FASTA方案，不放行封存数据或公开发布。

## 零更新准备与后续执行

四个新增人工检查仅覆盖角色margin的归一／逐base平移／padding与span外语义，角色量进入正确exon和junction输入列及相触gap0，新列和联合上下文的梯度但冻结局部源与unknown无直接梯度，以及参数量／空和全unknown可微零；不重复旧角色映射、集合DP或raw接口完整测试。随后读取原完整人口，使用未修改的八窗metadata规则，核对冻结A3重建的局部路径分数与已闭合A31同源结果、原身份／标签／完整分母、四上下文和七头实际适用梯度、新3／5列梯度、无known连接时零连接梯度及两物种连接覆盖。学生和冻结源参数前后相同，optimizer0，不保存checkpoint；不按新分数换案例或判效果门槛。

本次唯一PREFLIGHT使用generanno／private1RTX3090至少24GiB／2CPU／16G／5分钟／no-requeue及唯一日志。依据A39八窗准备46秒，增加冻结局部前向及四项人工检查留余量；不是随机计时或吞吐结论。失败保留原输出，不自动改变机制或重复预检；已关闭A39结果不受准备失败影响。

预检通过后，仅准备并提交一次固定fit及完整natural／CPU链，fit1:45、natural30分钟，与A39相同GPU／CPU／内存，并据真实准备耗时和活调度核实。自然评分必须全部1491／16321完成后才连接原保存标签；CPU须核验全部4608 draw、两损失组合、三checkpoint形状／有限FP32值和实际学生分支AdamW曝光、冻结源无更新、保存role-token来源及score组件、原DP／较强基线／完整分母／五条件。尚未实现或运行的阶段不能写成已完成。

实际局部输入收据计划在完整natural时保存1491份原FP32 15类logits，预计约2.2GB，加新学生／optimizer／冻结源三checkpoint及JSON，总新产物低于3GiB。该局部张量是评分来源核验收据，不参与拟合缓存；CPU只用其重建原候选角色输入、全1491窗A31_A3源路径分数与当前组件，不作新神经前向。三学生checkpoint另保存冻结A3源state，CPU直接逐tensor与原A3末权重比较，并核对学生optimizer没有源参数；不使用哈希。因新CPU需读取全量局部收据，采用0GPU／2CPU／8G／15分钟，而非A39无此收据的5分钟；依据A36全量局部张量核验实际8:33留余量，在首次A40拟合前固定。

独立输出M28-A40-ROLE-EVIDENCE-CONTEXT-{PREFLIGHT,FIT,TRAIN-READOUT,RESULT-CHECK}，不覆盖任何旧数据。首个拟合工程失败保留；已有更新、TIMEOUT或不明更新不续训、重提或延期。长等待仅复用固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同heartbeat coding-rna-r5-light-monitor，gpt-5.6-luna／low／15分钟；准备阶段短作业由主消费，当时保留心跳暂停及旧对象，整链提交后才切换到本轮。

## 已完成的零更新准备

实际13418097 root／batch／extern均COMPLETED0:0，2026-10-06 20:46:45–20:47:28 UTC，43秒，private gpu034／1RTX3090／2CPU／16G／5分钟，batch MaxRSS1949404K、Requeue0／Restarts0。STATUS／JOBID／summary相符，四项新增人工检查PASS0.787秒；stderr仅通过日志，无异常。新模型、测试、preflight、SBATCH和事前协议与canonical逐字一致，提交前五项活调度／预约与输出碰撞检查通过。

全1491／1536／16321及原三类标签／2916和1086分母相符。八个原metadata窗口均未替换，冻结A3重建的原A31_A3全部候选路径分数最大误差逐窗0；学生613348、上下文364256、头249092及冻结源649457参数与协议一致。四上下文与七头适用分支有限非零，新3／5列实际有梯度；原无known连接的水稻样本两连接分支和新5列均精确零，两物种known连接覆盖成立。unknown直接gain梯度零，局部源无梯度，学生与局部源权重前后相同，optimizer0、无checkpoint，没有新搜索、GEN／R4、DEV、Setaria或A33访问。

八窗内部load／两局部上下文前向／学生backward合计约12.33秒；首窗3.065秒、其余七窗0.558–2.410秒、平均约1.323秒，峰值allocated174315008 bytes。这不是随机总体计时或吞吐结果。结合A39实际64:31拟合和本次追加局部源前向，保留1:45拟合上限，不因准备通过追加训练曝光或宣称科学PASS。完整自然读出及CPU核验仍未运行。

## 一次研究链提交与等待交接

拟合／完整natural／CPU入口、两项新配对和gate检查、三SBATCH及协议已完成语法检查和local／canonical逐字比较；没有改原数据或冻结科学条件。首次sbatch指定已完成预检13418097的afterok，被Slurm以Job dependency problem拒绝，没有创建job或拟合输出根。实际sacct及STATUS／JOBID／summary仍证明预检COMPLETED0:0，而controller的scontrol已Invalid job id；同名squeue／sacct无新拟合。仅移除已满足的调度依赖后，第一次实际创建以下研究链，程序内仍核验预检终态、原配置和八窗准备。不是既有fit失败后的重提、续训或修改参数。

| 阶段 | 实际作业 | 固定资源及依赖 | 首次完整调度核实 |
| --- | --- | --- | --- |
| fit | 13418488 | private 1RTX3090／2CPU／16G／1:45，程序固定preflight13418097 | RUNNING 3:24，gpu034，实际分配与请求相同 |
| natural | 13418506 | private 1RTX3090／2CPU／16G／30分钟，afterok13418488，export A40_FIT_JOB=13418488 | PENDING 0:00，正常Dependency，尚无实际分配 |
| CPU | 13418507 | private 0GPU／2CPU／8G／15分钟，afterok13418506，export A40_FIT_JOB=13418488,A40_READOUT_JOB=13418506 | PENDING 0:00，正常Dependency，尚无实际分配 |

三个canonical Command／SubmitLine／WorkDir与预定脚本、参数、export相符，Requeue0／Restarts0。fit实际于2026-10-06 20:58:56 UTC开始；2026-10-06 21:02:41 UTC其STATUS RUNNING／JOBID13418488相符，training.jsonl当时192行。stderr只有两项新增配对／严格五条件检查PASS0.003秒，无异常；192行不是4608完成，不读取中途loss或据此判断效果。提交前新产物低于3GiB计划与实际可用157TiB相容，五项活调度／预约及输出碰撞检查通过，没有取消或修改其他作业。

同一固定监控及heartbeat已通过原生更新切换本cohort M28-A40-ROLE-EVIDENCE-CONTEXT-NATURAL-20261007，实际保存kind heartbeat／ACTIVE／固定target／15分钟／既有failed_runs_only及完整白名单prompt逐字读回一致，created_at1790981930511保留；未新建对话／cron或覆盖模型。只监控13418488→13418506→13418507，排除已消费预检、A39诊断和所有旧cohort。三根实际COMPLETED0:0、全部规定产物及七个literal CPU证明成立才交接，fit完成不单独交接；科学gate不在监控读取范围，FAIL不是工程故障。正常状态静默，整链完成或可行动故障回传一次后暂停同heartbeat。手动启动实际2026-10-06 21:04:16 UTC核实fit RUNNING／gpu034及两下游PENDING／0:0正常Dependency，没有主对话交接；固定对话turn01a11307-b84f-70b1-a259-53cdf2a6e30b completed／error null／38.514秒，主已消费cursor b7267a5e-89dc-4556-b8a7-fb7096bf58b7:86。定时运行现已首次验证，见下节；终态自动回传及最终科学五条件仍待完整链闭合，不继承A39证据。

本轮实际手动turn_context 2026-10-06 21:03:55.511 UTC核实既有gpt-5.6-luna／low及同cwd，没有依赖旧轮模型记录推定。配置、手动及首轮定时检查成功已分别验证；终态回传仍待实际发生。

## 与拟合并行的人工 raw FASTA 接口检查

在等待现有研究链时补充infer_role_evidence_context入口，只处理原自然A3 greedy＋beam4候选、同gain0及原assembly／GFF3语义。每方向窗口GEN仅编码一次；冻结A3的SharedContext仅前向一次，同hidden用于原自然生成与15类segmentation，再把实际原FP32 logits作为新学生第五个输入。学生重新计算其独立上下文，原角色映射及新3／5列不改；既有监督源及学生曝光须分别记录，不称全模型从零或匹配容量消融。加载时直接比较fit保存的冻结源state与原A3源，不增加哈希。

入口在加载模型／GEN／FASTA或创建输出前，要求现有三root终态与匹配配置、1491窗／16321链／1491局部收据、七项实际CPU证明、2916／1086完整分母及较强A31_A3原五条件实际全通过。旧A39改进、缓存true标志或缺失源证明都不替代该判断；还须另行固定raw-FASTA协议，文件存在不构成数据访问或native执行授权。本轮不建立该DEV协议，不提交真实DEV／独立物种或评价，也不读取训练中权重或评分。

六项新增人工检查只验证新A40身份和七个证明、强基线条件／缺证明在native前阻断、缺协议阻断、人工不一致冻结源在学生加载前拒绝、实际第五输入为同源FP32 logits且自然几何／来源／非正候选与score组件不变，以及共享一次的正负方向、ambiguity转N和负链坐标phase。人工query和小tensor不含研究DNA、features、labels或checkpoint，不重复旧生物学边界全集。一次RAW-CHECKS使用generanno／private0GPU／2CPU／8G／2分钟／no-requeue；根据A39五人工检查实际18秒留余量，独立输出与日志。失败保留并报告，不自动重复，且不改变正在运行的fit／natural／CPU或其监控对象。

实际RAW-CHECKS 13418751 root／batch／extern均COMPLETED0:0，2026-10-06 21:11:26–21:11:43 UTC，17秒，private gpu034／0GPU／2CPU／8G／2分钟，batch MaxRSS738476K、Requeue0／Restarts0。STATUS COMPLETED／JOBID13418751匹配，六项新增人工检查PASS0.255秒；stderr仅通过日志、stdout0 bytes。新入口／测试／SBATCH／事前补充协议与canonical逐字一致，提交前五项活调度／预约与独立输出碰撞检查通过；Command／SubmitLine／WorkDir／实际分配与固定脚本相符。短检查由主消费，不加入监控或重启研究链；DEV协议及DEV-INFER／DEV-EVAL根实际不存在，故真实native与完整评价未运行，不把它们记作失败。本轮检查只支持人工接口路由与科学阻断，不是实际预测、吞吐或泛化证据。

## 与拟合并行的完整评价入口准备

evaluate_role_evidence_context沿用原完整DEV尺，不读取研究结果准备新阈值。入口先要求原TRAIN三根实际通过、另行固定的DEV协议和完整真实推理作业；当前补充实现不建立该DEV协议或提交真实推理／评价。新identity固定A40配置、末004608、冻结A3局部源与实际fit一致且未更新、649457源／613348学生参数及分别4608监督更新；要求GEN和A3共享上下文各每方向窗编码一次、原8690窗／四scope／完整contig长度。逐窗保留自然几何／顺序／greedy和beam4来源及空候选窗，与较强A31_A3配对，核对原FP32 score组件和局部源标识，然后重建原assembly并与A40.gff3负链坐标／phase和预测数一致。保存预测与GFF闭合前不读取参考。

参考侧保持7728 primary链、原两物种和背景区域分母，保留native不支持的参考strata。对照仍为较强A31_A3、不同生成器A31_A23、A3／A6／A／R5_C0及ANNEVO／Tiberius／Helixer；复用旧纯函数comparison_table，报告per-species、pooled、macro和worst-species及原stage coverage，不把重复DEV称独立泛化或曝光／容量匹配因果、原生运行时间公平比较或论文PASS。该实现仍待真实推理与评价，未产生性能数值。

四项新增人工检查覆盖A40源／参数／曝光／旧arm／完整网格错误在参考访问前拒绝，逐窗局部源与FP32组件及空窗完整性，A40负链GFF重建和计数先于参考，以及完整对照名／原参考和背景分母／macro／worst-species计算。仅人工metadata、JSON和GFF，读取既有长度常量而非任何研究FASTA／参考／checkpoint／feature／labels。一次EVAL-CHECKS使用generanno／private0GPU／2CPU／8G／2分钟／no-requeue及独立输出日志；依据本轮RAW-CHECKS实际17秒留余量。失败保留并报告，不自动重复，不改研究链或固定监控对象。

实际EVAL-CHECKS 13418923 root／batch／extern均COMPLETED0:0，2026-10-06 21:18:44–21:18:59 UTC，15秒，private gpu034／0GPU／2CPU／8G／2分钟，batch MaxRSS737352K、Requeue0／Restarts0。STATUS COMPLETED／JOBID13418923匹配，四项人工检查PASS0.289秒，stdout0 bytes；stderr保留既有scripts/eval_structure_diagnostic.py第51行open(path)的unclosed file ResourceWarning，GFF重建及计数通过，没有把警告隐藏或当作生物学结果。评价器／测试／SBATCH／事前补充协议与canonical逐字一致，实际Command／SubmitLine／WorkDir／资源及提交前五项调度检查相符。短检查由主消费；真实DEV推理、参考评价和独立物种未运行，科学结论不变。

## 首轮定时监控证据

同固定对话实际heartbeat触发2026-10-06 21:18:45.972 UTC，automation_id coding-rna-r5-light-monitor、current_time_iso 21:18:45.720 UTC；turn_context 21:18:45.767 UTC核实gpt-5.6-luna／low／同cwd。turn01a11315-4dda-7733-9105-c0bbcaf3c91e completed／error null／18.456秒，主消费cursor b7267a5e-89dc-4556-b8a7-fb7096bf58b7:88；实际21:18:54 UTC快照为fit13418488 RUNNING／gpu034、natural13418506及CPU13418507 PENDING／0:0正常Dependency，没有产物／性能读取或主对话交接。此证据只证明本轮一次scheduled检查已追加到原监控对话；同heartbeat仍ACTIVE，终态自动回传未发生，不把它记为整链成功。

## 闭合结果的链级变化诊断准备

diagnose_role_evidence_attrition只在原fit13418488／natural13418506／CPU13418507实际终态、markers与配置相符且七项真实CPU证明成立后，读取已闭合的原JSON分数及保存标签。科学PASS或FAIL均可做描述性归因；不以FAIL拒绝分析，也不将归因用于改变原五条件或放行DEV。只比较较强A31_A3→A40和已关闭A39→A40，复用未改的diagnose_structured_emission_attrition.transitions，逐物种及pooled保留known1、known0、protected−1和正确至少5-exon四类，报告保留、新增、非正丢失、正分未选及其净变化。

完整1491窗／16321自然配对、三臂原DP计数和2916／1086参考分母必须与闭合结果相符。没有模型／feature tensor／局部logit收据／TRAIN catalog或参考访问，没有GEN／神经前向、optimizer、新候选搜索、新R4、独立物种、Setaria或A33。A40相对A39改变输入维度、初始化流和既有监督局部源；描述性变化不是角色输入的因果消融、匹配容量／曝光或去重gene F1。

两个新增人工检查仅覆盖新三臂配对、未知与长链状态和净量，以及科学FAIL仍可分析但任一七项CPU证明缺失必须在分数读取前阻断。它们计划在一次实际后处理作业中先运行，不新增重复人工检查作业，不重复旧transitions全集。后处理使用generanno／private0GPU／2CPU／8G／5分钟／no-requeue，独立M28-A40-ROLE-EVIDENCE-CONTEXT-ATTRITION输出及唯一日志，依据A39四臂归因实际30秒留余量，新输出预计低于1MiB。失败保留日志并停止使用该归因，不修参数或自动重提。当前仅准备入口，不提交或读取未闭合的研究分数；主消费后处理，不修改原研究链和固定heartbeat对象。
