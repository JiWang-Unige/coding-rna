# A33 局部角色分数驱动的全局搜索保留

2026-10-06，首次真实执行前固定。A31完整DEV与三格诊断支持进一步研究完整链生成／保留；A30原五条件失败、A32未执行的结论保持。A33是冻结A3上的独立零更新搜索组件，不救A30系数、checkpoint或margin，也不预判论文新颖性或性能。

## 搜索差异与固定条件

原beam按累计指针log-softmax／EOS log-sigmoid排序，最终A31 caller按15类局部角色相对intergenic0之和选择。A33仅将全局entry池、活跃continuation池和完成池的排序主键改为同一局部角色gain；原累计指针／EOS分数只用于平分次序，不新增混合系数。entry的未闭合预览为首ATG三base的START gain，或已有闭合前缀加待接内含子和acceptor首base gain；这是有界保留启发式，不宣称精确最优。

每个history的START／NEXT_START仍按原几何指针保留局部top4；所有原DNA合法DONOR／STOP事件、真实codon carry、GT／GC—AG、至少5bp内含子、连续START／STOP、Q16、beam4、completed4及max128不变。A27中全部query局部排除的315个原遗漏配对不由此自动解决；A33不宣称覆盖全部生成错误。

闭合DONOR前缀分数按原fine角色优先级累计，跨外显子保留coding phase；完成STOP后必须与原local_path_scores逐链一致。START覆盖短首exon的donor，STOP覆盖三base末exon的acceptor，一base内部exon的donor覆盖acceptor。FP32模型不变，FP64区间求和，不加长度归一化、温度、先验、阈值或角色权重。未来role增量可能为正，不能复用原非正log增量停止条件；只在全部保留live前缀的“已累计gain＋未访问位置最大非负role gain之和”严格低于完成池末名时提前停止。该上界不恢复已剪枝history，不构成全局精确解。

主方法保留原A31 A3 greedy＋beam自然池全部几何和来源，再加入新搜索的完成链，独立标记emission_beam4，按几何去重后用原局部分数和null0非重叠选择。旧池成员不能因新搜索没有重返而删除；本轮新增搜索预算／运行代价及池规模全部报告，不能称与旧方法计算量严格相等。

## TRAIN 预检与后续支持条件

唯一模型为A3末outputs/M28-A3-FIRST-FIT/step_004608.pt，原649457参数与4608更新不变，optimizer0、新参数0。只读原六条TRAIN染色体的cached特征，GEN revision和24576／真实尾窗几何不变；不做GEN编码、训练、DEV、sealed或Setaria访问。

先五项合成检查，再用原metadata固定八个TRAIN时序预检案例，不能按新效果换案例。核对增量role和直接logit／log-probability差求和、原自然池greedy+beam的实际几何／顺序／来源、原A31缓存局部分数、全部新完成链的语法和最终gain；记录新搜索耗时、深度及剪枝、池规模和GPU峰值，要求模型参数逐tensor不变、无梯度。TRAIN catalog仅用于原metadata案例选择，不进入候选API或新标签统计。八窗只判断工程接线及运行规模，不判科研效果门槛。

预检闭合后才准备全1491唯一TRAIN窗的自然读出和独立CPU重算；全部自然池／分数关闭后再开原R4参考标注，unknown仍参与选择并单列。新主对照是较强的A31_A3而不是原greedy。后续五项支持条件在这里预先固定：两个物种known-label window/chain proxy F1各自严格增加，两个物种proxy precision各自不下降，合并至少五CDS片段selected correct严格增加。沿用原eligible参考／unknown口径和整数交叉乘法，不称gene F1。任一FAIL即关闭本配方，不救阈值或预算；全部通过也只支持另定完整raw-FASTA协议，不自动开放独立测试或产生论文claim。

## 执行与失败处理

本轮仅提交八窗预检：private-teodoro-gpu／1 RTX3090 24GiB／2CPU／16GiB／30分钟，generanno。相较旧A31八窗，角色上界可能更松、访问深度更多，30分钟是首次观测的操作上限而非速度预测；不由此延长旧作业。唯一输出outputs/M28-A33-EMISSION-BEAM-PREFLIGHT，日志沿用outputs/M28-PILOT-R5/logs/M28A33EMITPRE_JOBID.*；命令sbatch --parsable sbatch/M28-A33-EMISSION-BEAM-PREFLIGHT.sbatch。

产物summary.json、八行natural_scores.jsonl、STATUS、JOBID及Slurm日志。实现／合成失败只允许一次保持上述科学条件的最小工程修复，原失败保留。数值、原池身份或语法不复现时停止解释；TIMEOUT不是工程通过，不自动改搜索范围、case或预算补救。没有首次新拟合、模型选择、GFF发布、Git提交或goal状态更改；若仅剩长等待，复用原固定轻量监控及同一心跳，不创建新对话。

## 一次表示类型修复

首次作业13402858 root／batch FAILED／1:0、extern COMPLETED／0:0，实际09:54:41至09:55:44 UTC、1:03、batch MaxRSS1616692K。五项合成检查PASS／0.888秒；真实首窗在新搜索之前的旧池相等检查停止，natural_scores.jsonl为0行。失败不是A33科学FAIL，也不是预检通过。

原completed_union的坐标为tuple，JSON反序列化为list；canonical轻量复现证明同坐标在该直接Python比较中仍不相等，规范为list后相等。R1只在比较前规范坐标容器表示，仍逐值／顺序／来源核对完整旧池，不忽略差异、不删窗或放宽数值检查。新搜索组件、参数、八例、语法、评分、预算和所有科学条件不变。保留原FAILED目录及日志，R1使用独立outputs/M28-A33-EMISSION-BEAM-PREFLIGHT-R1；原五项已通过的合成检查不重跑。若R1仍有真实身份／数值差异或其他失败，停止本轮，不做第二次修复。

## 预检闭合与完整 TRAIN 执行冻结

R1作业13402927 root／batch／extern均COMPLETED／0:0，实际10:00:06至10:01:29 UTC、1:23、batch MaxRSS1773576K；summary、八行自然分数和完成标记齐全。八窗旧几何／顺序／来源与A31分数全部复现；新完成链的增量gain和直接log-probability差误差均0，原指针／EOS重放最大误差2.0796884200535715e-05，满足原FP32相对／绝对容差而不是零误差。128个query深度2至24，cap存活0；模型权重不变、无梯度。八窗新增唯一链分别5、7、6、12、8、16、7、23，不知其正确性。内部八例总43.422秒，新搜索单例0.526至3.465秒；选定metadata案例不是代表性速度基准。

在任何A33新效果标注前固定完整读出：沿用1536抽样首次出现次序的1491唯一TRAIN窗，按ordinal modulo 4分到四个自然读出array任务（373／373／373／372），合并时恢复唯一原次序。各任务private-teodoro-gpu／1 RTX3090／2CPU／16GiB／2小时；这是根据八例运行规模设置的完整执行操作上限，不是预测时间，不调整Q16／beam4／completed4／cap128。实现检查先以CPU／2CPU／8GiB／10分钟执行，afterok放行array；全部array成功后CPU／2CPU／8GiB／40分钟合并标注；再以同资源独立CPU核验。任何缺窗、身份／数值／语法差异或作业失败均不解释效果、不自动重提／延长预算；原记录保留。

每窗原A31_A3候选直接复用已闭合原池，不重跑旧beam；旧局部分数在新冻结forward上逐值复现。新搜索的每条完成链仍进行原指针／EOS和raw grammar重放、增量／直接局部分数核验。保存FP32局部logits用于后续CPU独立直接求和，预计全体约2.2GiB，不存新GEN特征或参考注入。四份自然文件及模型不变检查全部关闭后，合并程序才打开原TRAIN R4 catalog；重算旧／新标签、unknown选择、2916原eligible配对分母、每窗计数和五项固定支持条件。独立CPU再从原TRAIN序列、catalog和保存logits重算几何／来源、原始语法、直接log-probability分数、标签、选择、分母及整数gate，无模型／GEN forward。不能由仅GPU完成或summary的PASS字段替代独立核验。

根路径outputs/M28-A33-EMISSION-BEAM-TRAIN-NATURAL/shard_0至shard_3、outputs/M28-A33-EMISSION-BEAM-TRAIN-READOUT、outputs/M28-A33-EMISSION-BEAM-RESULT-CHECK；各有STATUS／JOBID／summary。自然文件含ordinal以恢复原次序；读出natural_scores.jsonl及windows.jsonl各1491行。此后若只需等待，更新原固定监控对象与同一心跳，旧A31已消费作业不重新监控。提交和预检通过不等于完整TRAIN gate通过，更不等于raw-FASTA或独立泛化结果。
