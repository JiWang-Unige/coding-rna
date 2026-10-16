# A31 完整 raw FASTA DEV 注释协议

2026-10-06，首次原生DEV执行前固定。A31 A3和A31 A23的1491窗TRAIN读出13398521及CPU核验13398524均已COMPLETED0:0，两个variant独立满足原五项准备条件。因此两者都进入完整raw-FASTA测量，不按DEV结果选择模型、checkpoint或评分配方。此阶段仍是反复使用的DEV，不能证明独立泛化或论文竞争性通过。

## 模型和完整输入

分别固定A3末outputs/M28-A3-FIRST-FIT/step_004608.pt和A23末outputs/M28-A23-RELATIVE-BOUNDARY-FIT/step_004608.pt。参数、原4608更新及原训练曝光不变；不新增评分头、优化器更新或其他训练。沿用[a31_local_path_caller_protocol.md](a31_local_path_caller_protocol.md)的局部15类角色logit相对intergenic0之和，FP32模型、FP64求和、全部角色权重1、null0，无阈值、长度先验、归一化或温度调整。

输入是既有固定DEV原始FASTA：拟南芥NC_003074.8，长度23459830；水稻NC_089041.1，长度29936421。两个contig均处理正负方向，沿用原自然网格共8690窗、真实尾窗长度、非ACGT转N和同一四块原生GEN上下文，revision b0483c23b6b63787b61a6d3a204a9b517d6ba345。不读取标签或feature cache，不开放sealed或Setaria。原生编码及自然greedy＋beam4／completed4／cap128联合池从DNA重新计算；不把TRAIN的缓存池读出称为本次原生执行。

## 注释与完整评价

所有自然候选及来源，包括非正局部分数，先写入各窗JSONL。保持原ORF、5bp内含子、相位、nearest-owner、重复几何最大分数及span非重叠装配规则，不注入参考链，不按R4标签筛选或补回候选。每个variant保存四个方向scope和独立GFF3。原生入口仍先重算本variant的TRAIN门槛与核验身份，再加载模型、FASTA或创建输出；本协议文件是科学合同，不是人类授权凭据。

CPU评价先核对全部8690窗、候选到装配计数、实际GFF坐标／相位及末checkpoint，再打开参考。沿用全部7728条完整primary参考和全部背景分母，逐物种及合并报告exact-chain precision／recall／F1、CDS碱基指标、结构诊断、长链及free→owner→selected覆盖损失。未知或不支持参考层不能通过缩小分母消失。比较同输入的ANNEVO、Tiberius、Helixer、A3、A6、A和R5 C0；原pilot gate只作辅助，不替换本轮完整测量或原TRAIN五门槛。

A31 A3与A6必须逐窗保持同一冻结A3自然池的几何、顺序和greedy／beam来源；A31 A23是另一个生成器，只核对相同原始网格和上下文，不要求或宣称同池。A3 greedy主对照同时改变池与分数；A6有额外head训练曝光，A23改变生成器和拟合后的共享权重。本轮不是匹配容量、多seed因果实验，也不能与cached A3运行直接作公平耗时比较。另行三格TRAIN诊断只作机制解释，不改本协议、门槛或模型。

## 资源和停止规则

沿用实际generanno环境。相同8690窗原生A6执行内部11614.275秒，约3小时14分钟，GPU峰值2455940608 bytes；A31不加载质量头，但局部评分及A23较深搜索可能改变耗时，尚无本次原生速度实测。A3请求private／1RTX3090 24GiB／2CPU／16GiB／4:30，A23首次同输入原生推理保守请求5:30，仍为1GPU／2CPU／16GiB。两者独立并行；各自CPU依赖afterok，仅0GPU／2CPU／8GiB／20分钟。提交前分别检查实时Slurm、维护、并发、存储和唯一目录，不覆盖任何旧产物。

命令分别为sbatch --parsable --job-name=M28A31A3DEV sbatch/M28-A31-LOCAL-PATH-DEV-INFER.sbatch A3和sbatch --parsable --job-name=M28A31A23DEV --time=05:30:00 sbatch/M28-A31-LOCAL-PATH-DEV-INFER.sbatch A23；CPU为对应afterok及job-name=M28A31A3EVAL／M28A31A23EVAL的sbatch/M28-A31-LOCAL-PATH-DEV-EVAL.sbatch A3／A23。输出各自为outputs/M28-A31-A3-LOCAL-PATH-DEV-INFER／DEV-EVAL和A23对应目录；日志沿用outputs/M28-PILOT-R5/logs/%x_%j.*。入口要求输出未存在，SBATCH不能预建该目录；JOBID／COMPLETED仅在原生入口成功返回后写入，入口未创建目录的失败只保留Slurm和日志证据。

工程失败保留原产物，按现有授权仅允许一次最小明确修复，不改变科学配方；实际科学劣势不是工程错误，不救阈值、另换checkpoint或补训练。完整DEV注释及评价闭合后，主对话再判断性能差距与下一机制；不得由TRAIN proxy或单个推理完成宣称独立测试、竞争性或论文PASS。等待只复用原固定watcher和同一个15分钟心跳，沿用其模型设置，不新建对话或自动化。
