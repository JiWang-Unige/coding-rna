# A15 超时取证与工程恢复

2026-10-04 12:58 UTC。原A15拟合因时限结束，未产生可读出的head；五项TRAIN科学条件尚未执行。原协议、A14已闭合FAIL及论文目标不变。本文记录一次针对具体执行失败的工程处理，不是新seed、模型选择或增加科学曝光的实验。

## 已确认的失败

原fit13367306运行11:53:03–12:38:30 UTC，root TIMEOUT／0:0，batch CANCELLED／0:15，extern COMPLETED／0:0，实际45分28秒。`outputs/M28-PILOT-R5/logs/M28A15FIT_13367306.err`明确time limit。原目录`outputs/M28-A15-NATURAL-CHAIN-SET-FIT`仅有训练日志和标记，未存任何checkpoint或summary；STATUS遗留RUNNING，Slurm终态为准。

日志最后持久化step896、epoch1、draw ordinal895。没有模型／AdamW状态，不能精确续训；未刷盘尾部未知，896不是全部物理更新的保证。该次失败训练的权重丢失，后续若恢复须从同A6 seed0 fresh初始化重启，保持原1536 draw／权重／顺序×3＝4608；明确区分最终模型保留的4608更新和失败尝试消耗的计算。

原下游13367308 DependencyNeverSatisfied、13367310 Dependency，未执行，不视为读出或核验FAIL。旧目录、日志及下游均保留，不覆盖或取消。固定监控已完成一次故障交接并实际暂停；消费故障后同一监控当前仅包含有效A17三段链。

## 零更新计时：已完成

计时13367598经Slurm提交，private／1RTX3090／2CPU／16GiB／10分钟，2026-10-04 12:56:59–13:00:25 UTC，root／batch／extern均COMPLETED／0:0，206秒。独立输出`outputs/M28-A15-NATURAL-CHAIN-SET-TIMING`，日志`M28A15TIMING_13367598.out/.err`，summary completed／optimizer0／head不变／生成器冻结。2026-10-05主对话读取实际终态；日期跨越不沿用旧RUNNING快照。本机AST／bash语法检查通过，两新增文件明确逐个scp并与本地直接字节比较一致。

固定原前32 draw，沿用原cache、冻结A3生成器、自然候选／标签、集合损失及importance权重，分别计时：cache读取、输入转换／传输、冻结生成器、head、set loss、backward／clip。无optimizer step，验证head权重不变、冻结生成器无梯度；记录实际hidden形状／dtype、预估全部唯一窗hidden大小及实际GPU容量。无新GLM编码器前向、候选搜索、DEV、sealed或Setaria。

六段均值依次为3.625629／0.015175／0.126496／0.017663／0.007802／0.052422秒；包括首个draw的CUDA初始化，后续冻结生成器约0.067秒。读取占六段计时之和约94.3%，不把它说成整个job的94.3%。population另30.325秒，脚本总165.465秒，调度总206秒。这是当前读取主导的直接证据，不是历史45分钟超时的因果证明。

每窗hidden实际24576×128／FP32／12,582,912字节，1491窗合计18,761,121,792字节（约17.47GiB）；本作业GPU容量25,351,356,416字节，单窗峰值115,272,192字节。可在同类3090上一次生成并重用冻结hidden，无须写出18.8GB中间文件。新作业仍会在分配的实际GPU上检查容量；不足时在optimizer前停止，不改变精度或科学配置兜底。

## 唯一工程恢复R1：事前固定

2026-10-05，保持原配置，从同A6 seed0 fresh初始化，不继承任何丢失或旧权重。独立`recover_natural_chain_set.py`，旧runner／checkpoint／候选池不变。1491原唯一TRAIN窗各读原cache并经冻结A3生成器一次，保留原FP32 hidden于显存；训练仍严格依原1536 draw及权重／顺序重复3遍，共4608步。编码器无新前向、候选不重搜、自然标签与loss不变。

optimizer前按原metadata选择八个真实窗口，重算live hidden，与复用hidden比较（FP32 rtol／atol均1e-6），并比较head分数／set loss／截断后的逐参数梯度；head权重不变，optimizer0、生成器／link无梯度。实际检查记录写独立fit的`cache_receipt.json`，不把检查通过预先当事实。训练每256步保存同格式head／AdamW状态，包含001536／003072／004608；主读出仍仅允许004608，不挑中间权重。

固定输出`outputs/M28-A15-NATURAL-CHAIN-SET-{FIT,TRAIN-READOUT,RESULT-CHECK}-RECOVERY-R1`，三段afterok链。资源保持generanno／private，fit与readout各1RTX3090／2CPU／16GiB，CPU核验0GPU／2CPU／8GiB。操作时限修正为3小时／2.5小时／10分钟：计时读取均值×1491约90分钟，加冻结生成器／head、准备与余量；读出原路径也需要重读1491窗，不能继续用旧20分钟预算。此为操作上限，不新增epoch或科学预算，不能保证实际耗时。

R1读出与只读核验直接复用原实现，仅路由到R1目录；核验仍检查实际4608更新／权重／标签、三epoch保存文件、全部16321自然候选和原五项gate，额外核对缓存等价的零更新receipt。原失败尝试896条持久记录及未知未刷盘尾部明确为丢弃计算，不声称总物理更新仅4608。

## R1实际提交与等待

2026-10-05 13:54:39 UTC，fit13387908实际启动gpu034；13:56–13:57 UTC即时squeue／scontrol核对两fit并行RUNNING，下游四roots全部有效Dependency，命令、依赖及实际资源／时限符合事前记录。

|阶段|实际job|实际依赖／状态|操作时限|
|---|---:|---|---|
|A15 R1 fit|13387908|RUNNING gpu034|3小时|
|末head自然TRAIN读出|13387947|afterok:13387908，PENDING Dependency|2.5小时|
|原CPU结果核验＋cache receipt|13387977|afterok:13387947，传入13387908／13387947，PENDING Dependency|10分钟|

对应SBATCH为`M28-A15-NATURAL-CHAIN-SET-RECOVERY-R1-{FIT,TRAIN-READOUT,RESULT-CHECK}.sbatch`，日志前缀M28A15R1FIT／M28A15R1TRAIN／M28A15R1CHECK加实际job ID。两R1脚本AST和六SBATCH bash语法通过，新增／更新12个精确目标逐文件scp，全部与本地直接字节比较一致；原A15／A17科学runner也实际与本地比较一致且未改动。A17逐draw循环另经AST归一化比较，仅新增保存分支与命名空间引用，原计算逻辑相同。本检查不冒充独立科学审阅。

复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同一heartbeat coding-rna-r5-light-monitor，当前cohort M28-A15-A17-RECOVERY-R1-20261005，六roots同时包括A17 R1 13387910／13387946／13387978。完整原生prompt及ID／kind／ACTIVE／固定目标／15分钟／原通知设置／created_at已实际回读一致，gpt-5.6-luna／low不变，不创建新对话、fork或cron。仅任一完整三段链结构齐全或新可行动故障才一次交接后实际暂停，主对话消费后为剩余链重新启用同一心跳。此次R1范围的定时检查／最终回传尚未验证；实际取证与触发时间分开记录，本机需保持开机和应用运行，不能保证休眠时即时回传，见[官方说明](https://learn.chatgpt.com/docs/automations?surface=app)。

同一监控对话的手工启动检查于2026-10-05 13:58:48.161 UTC实际取得六root squeue／sacct：两fit RUNNING gpu034，四下游有效Dependency，未读未完成结果、未向主对话发送健康消息。对应turn01a10c5b-ce25-7233-a54f-31b5a49b96a4已completed／error-null，实际turn_context仍gpt-5.6-luna／low，原生心跳ACTIVE。此为本R1 scope的单次手工启动证据，不是定时执行或最终自动回传证据。

未完成：1491窗hidden物化、八个cache/live零更新等价case、4608步末head、全部自然TRAIN读出、CPU实际记录核验和原五项科学判断。以上RUNNING／提交均不视为结果。无DEV、封存访问、其他作业修改／取消、Git提交／推送或公开发布。

## 不变的科学终点

仍仅允许末step004608、自然16321候选／1491 TRAIN窗、原标签／未知规则、原A6／A14配对和原A6五项继续条件，操作点0。不得从失败日志或旧A14继承权重，不改变loss、draw、学习率、epoch、seed、阈值或gate。工程恢复完成仍须真实4608记录及同原CPU只读核验；TRAIN gate闭合前不提交DEV。资源／保存频率等针对TIMEOUT的操作变化须在此先记录，保持科学配置不变。

## R1 fit 已完成 自然读出进行中

2026-10-05主对话实时核查：fit13387908的root／batch／extern均COMPLETED、ExitCode=0:0、Elapsed=00:26:53；实际训练日志到4608，summary的optimizer_steps=4608、末head为R1 step_004608，原CONFIG未变。cache_receipt实际记录1491窗／18,761,121,792字节，物化1386.390秒；八个metadata witness的hidden差均0，分数／loss／clipped逐参数梯度等价检查通过，检查期间optimizer0、head未变、generator冻结。三遍实际训练各1536步，耗时32.379／32.917／32.120秒。缓存物化是主要时间；这不是模型科学收益。

下游13387947已实际RUNNING，自然TRAIN评分日志从64继续到128／1491；CPU13387977仍有效Dependency。A17 fit13387910仍RUNNING、日志到1344／4608。原上节未完成列表中的hidden／cache检查／末head现已取得工程证据；A15全自然读出、原CPU核验及五项继续条件仍未闭合，不提前判断PASS。各模型继续按原协议，不进行新DEV／sealed／Setaria读取或参数修改。

固定监控对话中的turn01a10c69-4298-7762-a99e-b22fa5740c26已completed／error-null，内容报告完成同一R1 cohort的健康检查，未向主对话发送健康交接；这确认后续检查追加到了原监控对话，本次未单独核验其定时触发来源。它取证时两fit仍运行，不能用其旧状态覆盖上述新Slurm终态。完整定时终态自动回传仍未验证，心跳对象不变。
