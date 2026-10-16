# A31 固定局部路径 caller：执行记录

2026-10-06。按首次固定的[a31_local_path_caller_protocol.md](a31_local_path_caller_protocol.md)实施。全1491窗读出及CPU核验已完成，两variant均通过原五项TRAIN准备条件，详见[a31_local_path_caller_train_result.md](a31_local_path_caller_train_result.md)。这不是完整基因F1或独立泛化。A30独立运行，不改变其损失、曝光或原五gate。

## 已完成预检

唯一作业13398501，private／1RTX3090 24GiB／2CPU／16GiB／10分钟，实际generanno。提交前核对分区与GRES占用、账户作业、GPU待排、预计开始和维护预约；gpu034的8卡用5，gpu035的8卡用0，无相关GPU维护截断或输出目录碰撞，其他TE作业排除。

root／batch／extern均COMPLETED0:0，实际2:42，batch MaxRSS1655336K；STATUS与JOBID匹配，无失败stderr。

四项unittest均PASS，0.174秒：原标签正反方向ascending相位对应；logit差与直接log-probability差恒等；逐base平移、padding／span外不变；零分null、非重叠最优选择、空池、无梯度和非有限输入。八项原metadata案例未按新分数选择；双模型自然池均保持来源，原A3 greedy几何重放一致。16项模型／case直接评分差均0，全部冻结参数未变、无梯度、optimizer0，A3／A23参数分别649457／649617，无新GEN编码或DEV／sealed／Setaria。

summary内部23.773069秒，八项双模型计算分别6.6933、3.7613、1.2327、0.7986、1.3102、2.4815、2.4866、3.5270秒，包含原greedy来源重放和额外直接log-probability比较；GPU峰值161618432 bytes。这不是自然候选质量或最终caller性能通过。

## 全量实现和资源依据

readout先对原1491唯一窗关闭natural_scores.jsonl，再加载TRAIN catalog；不重做beam或注入参考。冻结A3与A23各自原联合池，原A3 greedy／query-logit为主基准，既有A6选出计数仅作辅助。未知候选照常参加选择；候选／窗口标签、eligible参考／长链分母、所有selected计数和known-label proxy均单独报告，五gate用整数比较。

CPU核验不做模型前向，重新读取原缓存DNA与R4，核对1491条自然／标签记录的顺序、原池和greedy子集、标签、原DP选择、参考分母、proxy和五gate；同时核对两末checkpoint结构与有限值。它复核保存分数而不独立生成局部emission，不冒称独立模型复现。

全量GPU计划private／1RTX3090／2CPU／16GiB／90分钟：八case平均约2.79秒／窗，简单外推1491窗约69分钟，加metadata和运行余量；正式读出不需要直接log-probability重复比较，但个别更深greedy及I/O仍可能更慢。CPU计划private／0GPU／2CPU／8GiB／30分钟，依据相同1491窗已有CPU约14分钟。输出目录独立，不写模型checkpoint，新增JSON为小量产物；未提交前不称已开始。

## 实际全量提交

两次提交各自完成五项实时调度检查；private gpu034／035有空闲3090、没有相关维护或并发限制，不需共享分区迁移。所有本机新runner AST／三SBATCH语法检查通过，窄范围上传完成；正式新目录不存在，未覆盖已有数据。

sbatch --parsable sbatch/M28-A31-LOCAL-PATH-CALLER-TRAIN-READOUT.sbatch实际返回13398521，提交后RUNNING gpu034。sbatch --parsable --dependency=afterok:13398521 --export=ALL,A31_READOUT_JOB=13398521 sbatch/M28-A31-LOCAL-PATH-CALLER-RESULT-CHECK.sbatch返回13398524；实际scontrol确认PENDING、有效afterok:13398521、30分钟、0GPU／2CPU／8GiB。GPU读出为90分钟／1GPU／2CPU／16GiB。日志在outputs/M28-PILOT-R5/logs/M28A31LOCALREAD_13398521.*和M28A31LOCALCHECK_13398524.*。

04:09 UTC取证：A31 GPU RUNNING、CPU有效Dependency；A30原fit13398266同时RUNNING、两个原下游有效Dependency。A31尚无全量效果结论。

同一原生heartbeat coding-rna-r5-light-monitor于04:08:52.926Z更新ACTIVE，目标仍原固定watcher01a0f71f-8845-7432-bd13-bf54c1808fb3。整个prompt／kind／target／15分钟／通知策略／原created_at精确读回PASS，A30字面结构白名单未改。当前cohort为M28-A30-A31-PARALLEL-TRAIN-20261006，只监控这五根、两条独立整链。任一整链终态或实际故障才一次交接并暂停同一心跳；gate及性能不是结构白名单条件。配置成功、一次手工检查、当前cohort定时验证与终态回传分别记录，不能互相代替。

手工启动实际turn01a10f67-2b59-7f32-86e2-3ce4d268154d完成、error null、32.072秒；04:09:56.110Z快照为A30 fit及A31 readout RUNNING，其余三根均有效Dependency。实际turn_context04:09:42.012Z仍gpt-5.6-luna／low、相同项目，没有覆盖模型或新建对话。当前combined cohort的定时执行与终态回传尚未观察，不借用旧A30-only定时记录。

## 全量 TRAIN 结果闭合

05:28 UTC主对话直接取证：readout13398521和CPU13398524的root／batch／extern均COMPLETED0:0，分别1:00:41／19:28；STATUS／JOBID匹配，natural_scores.jsonl和windows.jsonl各1491行。CPU重算实际R4标签、自然池和greedy子集、DP选择、2916 eligible配对分母、proxy及整数五gate，两个variant各五项PASS。未知选择不计已知假阳性但保留完整计数；结果不能当完整gene F1或与原生基线直接比较。

同一固定watcher在04:20、04:35、04:50、05:05、05:20实际scheduled检查过这五根；05:20记录A31读出完成／CPU运行、A30继续运行。终态由主对话直接消费，没有发生本批次自动终态交接。原生heartbeat于05:33:17.959 UTC由主对话暂停，保持原ID、kind heartbeat、target、15分钟、通知策略和created_at，避免已消费A31重复回传。A30三根仍未消费；正式DEV链提交后只更新同一个监控对象再启用，不能把暂停监控写成暂停研究goal或取消作业。

两variant的完整raw-FASTA条件已由实际TRAIN核验满足，首次原生执行前固定[a31_local_path_caller_dev_protocol.md](a31_local_path_caller_dev_protocol.md)。原caller／评价Python未改，先前八项合成工程检查不重跑；本轮只新增相称的两个SBATCH并检查bash语法。两个入口都要求输出目录未存在，不能沿用旧A6脚本预建目录，否则将产生实际FileExistsError。新SBATCH仅在入口成功返回后写JOBID／COMPLETED；入口创建目录后的失败保存FAILED，没有创建目录的失败保留Slurm及日志。

## 完整 DEV 实际启动

本机两个新SBATCH语法通过，连同固定DEV协议及TRAIN结果窄范围上传成功；远程相同bash检查通过，新阶段目录都未存在。BeeGFS实际可用155TB。四次正式提交各自完成五项实时Slurm检查，私有gpu034／035均有合适3090和CPU资源，无相关维护、并发限制或输出碰撞；其他TE作业不操作。

sbatch --parsable --job-name=M28A31A3DEV sbatch/M28-A31-LOCAL-PATH-DEV-INFER.sbatch A3返回13399434，private／1RTX3090／2CPU／16GiB／4:30。sbatch --parsable --job-name=M28A31A23DEV --time=05:30:00 sbatch/M28-A31-LOCAL-PATH-DEV-INFER.sbatch A23返回13399441，同硬件／5:30。两个CPU命令使用对应afterok:13399434／afterok:13399441及M28A31A3EVAL／M28A31A23EVAL名运行sbatch/M28-A31-LOCAL-PATH-DEV-EVAL.sbatch A3／A23，分别返回13399454／13399455。scontrol实际确认两个afterok仍有效、0GPU／2CPU／8GiB／20分钟，不附已结束的TRAIN依赖。

提交后七根快照中A30原fit仍RUNNING，两DEV推理均RUNNING gpu034，四个下游均有效Dependency；原生stdout已显示A3 160／1909、A23 96／1909的拟南芥正链窗进展，没有失败stderr。这是实际raw-FASTA编码及注释进展，不是8690全窗闭合、gene F1或竞争性结果。各阶段独立INFER／DEV-EVAL目录和M28A31A3DEV／A3EVAL／A23DEV／A23EVAL_ID日志，不写或覆盖模型checkpoint。

同一个原生heartbeat于05:46:31.117 UTC重启ACTIVE，当前M28-A30-TRAIN-A31-DEV-PARALLEL-20261006只包含七根、三条独立完整链；整个prompt及固定target／kind／15分钟／通知／created_at精确读回相符，原A30结构白名单字面保留。新A31白名单只取结构，A23评价matched_A6_natural_candidate_geometry_order_sources_unchanged须为存在的JSON null、A3须true，不猜字段别名或据科学指标拒绝交接。RUNNING时原生输出目录或STATUS暂缺按入口设计正常。

原watcher手工启动turn01a10fc0-889f-78f3-82fb-08ca8329d63e实际completed／error null／34.888秒，取证快照05:47:35.926 UTC确认三根运行和四根Dependency；实际turn_context05:47:18.613 UTC保持gpt-5.6-luna／low和相同项目，未覆盖模型。当前cohort的定时运行和终态交接尚未观察，不能由配置或本手工检查推定；本机休眠、应用关闭或连接失败时不能保证即时回传。短三格TRAIN诊断由主对话消费，不加入七根监控。

## 评分与候选池辅助诊断

短CPU作业13399561 root／batch／extern COMPLETED0:0／47秒，三项合成检查PASS，全部真实1491窗及原两格／2916参考分母精确复现。固定greedy池改用局部分数减少known_nonmatch239／95、正确改变−10／0；固定局部分数扩池增加正确406／53、known_nonmatch96／81，水稻precision下降5.86百分点。见[a31_score_pool_train_diagnostic_result.md](a31_score_pool_train_diagnostic_result.md)。这不是完整四格交互或生物学因果证明，没有模型／backbone前向、新标注、搜索、训练、DEV或Setaria访问，不改变原门槛或正在运行的两个DEV variant。短作业已消费，不加入心跳。

## 完整 DEV 及五状态诊断闭合

2026-10-06 09:19 UTC，主对话已消费原固定watcher的09:03:20.685 UTC终态交接，并直接核对Slurm。13399434／13399454及13399441／13399455四根及batch／extern均COMPLETED0:0；两variant均8690窗、完整7728 primary、对应STATUS／JOBID与非空GFF，原CPU已完成真实自然池、装配、GFF及参考评价核对。主对话由原整数复算全部对照指标一致，未改变模型、操作点或分母。

原A3合并完整链F1为52.3200%，A31 A3／A23为66.1671%／67.8133%，固定Helixer／Tiberius／ANNEVO为62.8512%／71.8047%／75.7025%。两个variant均高于Helixer，但仍低于Tiberius及ANNEVO；原pilot辅助门槛PASS不是论文通过。原始结果、逐物种／背景及限制见[a31_local_path_caller_dev_result.md](a31_local_path_caller_dev_result.md)，不是独立泛化或纯评分因果结果。

全部五source root闭合后、五项最新Slurm检查及新输出核对通过，实际提交sbatch --parsable sbatch/M28-A31-DEV-SCORE-ATTRITION.sbatch返回13402231，private0GPU／2CPU／8GiB／5分钟。root／batch／extern COMPLETED0:0／1:14，09:07:06至09:08:20 UTC，四项检查PASS／0.003秒，references.tsv共7729行含表头；saved自然分数／owner／选择与原评价计数全部复现，未作新前向、预测、标注、训练或阈值调整。正确链非正gain的合并损失A6／A31 A3／A23为2967／261／255，当前A23最大剩余损失是1654条未进入自由池，另有535条正分未选，不进行DEV阈值补救。

原同一heartbeat已PAUSED、固定target及15分钟保留，已消费四根不再重复检查或交接；研究goal仍ACTIVE。短诊断由主对话消费，不加入新监控。此前提出的三原生对照短人工FASTA预检没有提交，优先处理正式结果；没有新训练、A32实际运行、独立面板或Setaria开放。

## 固定评分与候选池 DEV 三格诊断闭合

2026-10-06 09:36 UTC。固定[a31_score_pool_dev_diagnostic_protocol.md](a31_score_pool_dev_diagnostic_protocol.md)后，新增脚本、三项合成检查及单个SBATCH，窄范围上传；五项实时调度、五个source root闭合及唯一输出检查通过。sbatch --parsable sbatch/M28-A31-SCORE-POOL-DEV-DIAGNOSTIC.sbatch返回13402574，private／0GPU／2CPU／8GiB／5分钟。root／batch／extern COMPLETED／0:0，实际09:31:19至09:32:09 UTC、50秒、batch MaxRSS898120K；三项合成检查PASS／0.047秒，STATUS／JOBID相符。

全部8690窗的旧greedy去重几何与新联合池greedy来源集合完全匹配，原重复query记录保留，两个原格完整stage、装配计数及整数指标复现。先关闭候选及三格参考无关装配，再读保存7728参考；实际中间分数8690行及四scope选择，主对话复算数值和边际一致。完整链合并F1为52.3200%→56.4037%→66.1671%；评分变化非匹配primary预测减少1060、TP净减少10，固定局部分数扩池TP增加1261、至少五CDS片段正确链增加705，同时非匹配预测增加600。见[a31_score_pool_dev_diagnostic_result.md](a31_score_pool_dev_diagnostic_result.md)。

只作有序条件解释，不称完整二因素消融或因果份额；缓存A3与原生A31的hidden数值未证明相同。未作模型／GEN前向、候选搜索、R4标注、feature cache、训练、阈值扫描、Setaria访问或生产GFF。本短作业已消费，不启用旧A31心跳；原TRAIN及科学门槛未改，没有提交新生成器训练。研究goal继续ACTIVE。
