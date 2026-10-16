# 三种原生基线的人工短 contig 预检

2026-10-06，首次真实执行前固定。A33完整TRAIN链继续独立运行，本预检只补充ANNEVO、Tiberius、Helixer接收短多记录FASTA的工程证据，不接入独立物种DNA／注释、DEV输出、参考或特征缓存，不做训练、参数选择、模型排名或性能评价。

输入复用已闭合13401642的人工两记录，核对其确定性DNA再统一大写：artificial_332长332bp，artificial_7536长7536bp且人工ORF跨6144bp位置。保留两条完整输入，不补长、拼接、复制、删短contig或把人工ORF视为准确率标准。新副本仅写入各方法独立输出；原输入、所有权重和冻结结果不变。

直接复用run_native_baseline.commands的完整原生操作点、GPU探测和容器绑定，仅替换显式FASTA参数。ANNEVO固定v2.3.3源码42c920f70f52f96c7c06c665a639ec47d9304fff／Magnoliopsida权重／batch32／8worker；Tiberius固定2.0.8源码871d4979568ae09e56de7eca5d7489046ee833ea，在原实际2.0.7依赖SIF内运行angiosperms／batch4／请求seq_len200070／no_softmasking；Helixer固定0.3.7／a_0080／64152／batch16／deterministic／no-multiprocess。Helixer species参数只作为人工样本标签，不更换其显式模型。Tiberius已有的clamp_seq_len_to_genome原生逻辑保留；实际缩短后的长度须从日志记录，不能假称请求窗口就是实际窗口。默认min coding／min protein／min sequence规则均不调，过滤与零预测如实保留。

每方法先执行两项新入口检查，覆盖只替换FASTA且保持全部其他参数、原生命令完成标记／真实空GFF／多contig坐标和CDS strand／phase／Parent。随后真实执行原GPU探测及原生命令，记录每方法一张真实GPU、实际backend、精确命令、日志与原始GFF；核对全部预测seqid属于两条人工记录、1-based坐标不越界。必须存在原生完成标记，特别是Helixer仅exit0不能掩盖其后处理失败。零预测在原生完成且GFF实际存在时是有效空输出，不救回或用于排名。

统一array三任务0=ANNEVO、1=Tiberius、2=Helixer，最多并发3；private-teodoro-gpu／各1 RTX3090 24GiB／8CPU／32GiB／15分钟，annevo环境、原容器和固定权重，零新下载或安装。15分钟是原完整DEV入口89至360秒之外留出的首次短输入操作余量，不是公平速度基准或预测。输出outputs/M28-NATIVE-SHORT-CONTIG-PREFLIGHT/{ANNEVO,Tiberius,Helixer}，各有人工FASTA／commands.json／GPU收据／command.log／prediction.gff3／summary／STATUS／JOBID；日志M28NATIVESHORT_ARRAYID_TASKID.*。提交命令sbatch --parsable sbatch/M28-NATIVE-SHORT-CONTIG-PREFLIGHT.sbatch。预计每方法临时输出不超过100MiB；不覆盖既有原生数组、DEV预测或A33。

现有失败要分别记录：工具缺失／非零退出／缺完成标记／坐标错误为工程失败，真实默认短序列过滤或零预测不是准确率FAIL。首次执行失败不自动重提、延长、降batch、补padding、修上游版本或调整过滤；保留实际错误并据其性质决定下一步。若只剩长等待，复用原固定监控／同一心跳，不创建新chat、fork、cron或第二监控。

成功只支持这两个人工长度的原生命令完成与已产生GFF结构接通；空输出不证明非空gene的phase／Parent。两个contig不能代表37194条碎片化基因组的整体规模、所有序列实际被模型考虑、独立泛化或科学性能。参考门槛、面板资格、未知上游暴露、A30负结果、A32未执行及A33冻结五条件都不因本预检改变。

## 实际执行结果

提交array13403567，各任务独立目录和JOBID alias对应原方法。三项root／batch／extern均COMPLETED／0:0，实际STATUS／JOBID／summary与原生命令完成标记齐全；各自两项入口检查PASS（0.003至0.012秒），远程Python AST及SBATCH语法通过。输入仍为完整两记录332／7536bp，原人工源及原生runtime身份核对通过；没有生物基因组、注释、特征cache、DEV输出、Setaria、训练、下载或安装。

| 方法／实际任务 | Slurm用时 | batch MaxRSS | 实际backend／可见GPU | 原始GFF字节／实际feature数 |
| --- | ---: | ---: | --- | ---: |
| ANNEVO／13403567_0 | 1:42 | 1157476K | Torch2.1.0／1 RTX3090 | 755／0 |
| Tiberius／13403567_1 | 13:17 | 3014460K | TensorFlow2.17.0／1 RTX3090 | 16／0 |
| Helixer／13403567_2 | 5:02 | 1415764K | TensorFlow2.15.1／1 RTX3090 | 215／0 |

独立读取三份实际GFF，非注释feature计数均0，与各summary的总计及逐contig计数一致；原生完成标记均实际存在。因此是原默认条件下的真实空预测，不是缺文件、后处理失败、删除短输入或筛除预测来救回。零feature使坐标与CDS phase／Parent检查为空；不能宣称观察到非空gene结构或所有输入都被模型考虑。

Tiberius实际日志明确把请求seq_len200070按最长contig7536缩至7524，保留原生clamp规则。其原生打印“Tiberius took 2.1714 minutes to execute.”只覆盖内部计时；Slurm实际13:17，驱动内部总计791.772秒，二者不混用。ANNEVO／Helixer驱动内部总计分别96.720／296.427秒。容器／启动／加载开销和输入处理范围不等，以上不是公平速度比较，不外推碎片化规模。

当前交付已闭合这两条人工输入的三种原生入口工程预检；没有人工ORF正确性标注、准确率、模型排名、独立泛化或适用范围优胜结论。原15分钟上限内自然结束，没有重提、延期、参数调整、padding或取消。所有实际输出保留在原独立目录，不改变A33失败状态或其五项科学条件。
