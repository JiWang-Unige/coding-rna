# A37 完整 raw FASTA DEV 注释协议

2026-10-07北京时间，原生DEV访问前固定的条件合同，当前未执行并关闭。实际训练13410061及自然13410739／CPU13410742已闭合，但原五项TRAIN条件仅通过2项，详见[a37_context_chain_result.md](a37_context_chain_result.md)。按原停止条件不提交本轮DEV、不救阈值或模型；以下原参数与计划保留为未执行记录，合同存在及工程PASS不能放行科学FAIL。论文总目标仍在推进，本轮没有独立泛化或原生性能结论。

## 源模型与执行前提

只有原fit、自然读出和CPU均实际COMPLETED0:0，且STATUS／JOBID／summary与各源作业匹配后，主对话才重算[a37_context_chain_protocol.md](a37_context_chain_protocol.md)和[a37_context_chain_readout_protocol.md](a37_context_chain_readout_protocol.md)的五条件。两物种已选正确须严格超过A31_A3的1531／478，known_nonmatch不得超过347／251，pooled至少5-exon正确须严格超过700；原2916和1086参考配对分母不变。实际四项CPU证明、4608更新、三个checkpoint及原自然16321候选身份均须相符，保存分量核验不是独立神经前向复现。

原生入口infer_context_chain.py在任何科研模型、GEN、FASTA或输出目录访问前再次执行这一核验。工程准备来自人工raw检查R1 13411198和完整评价检查13411914，均已实际COMPLETED0:0／五项PASS，但不替代上述科学条件。此合同不扩大数据或发布授权，不更改源TRAIN任何选择。

生成器固定A3末outputs/M28-A3-FIRST-FIT/step_004608.pt，649457参数／原4608更新；新评分器固定outputs/M28-A37-CONTEXT-CHAIN-FIT/step_004608.pt，612324参数／原4608更新。两者eval、全部参数冻结及无gradient，没有新optimizer。只使用原FP32 exon additive＋junction additive＋有序链nonadditive−窗口null的gain，空集gain0、操作点0，不调整尺度、温度、长度先验、阈值或权重，不选其他checkpoint或追加曝光。

## 固定原始输入与注释

输入仅为data/m1_screen/arabidopsis_thaliana/genome.fa中的NC_003074.8（23459830bp）和data/m1_screen/oryza_sativa/genome.fa中的NC_089041.1（29936421bp），不物化同一FASTA中的其他染色体。沿用原24576bp宽／12288bp stride及固定尾窗规则：拟南芥各方向1909窗，水稻各方向2436窗，共4个scope／8690窗。保持非ACGT转N、方向化、valid-length、原四块6144bp GEN编码及revision b0483c23b6b63787b61a6d3a204a9b517d6ba345，BF16 materialization后转FP32，与原TRAIN缓存口径相同。

每个方向窗口只运行一次GEN，两组件复用同一编码和DNA，各执行自己的上下文。原A3 greedy＋beam4／completed4／cap128及原ORF、5bp内含子和phase规则生成自然候选；A37仅给这些链评分，旧A3 query存在分或segmentation分不进入新gain。不读取reference、标签或feature cache，不注入参考、不补回缺失几何、不按标签筛选。保存全部自然链及来源，包括非正gain、四个FP32分量数组和scalar null；原nearest-owner、几何去重最大分数及span非重叠装配不变。

四个方向JSONL、完整summary及A37.gff3写入唯一outputs/M28-A37-CONTEXT-CHAIN-DEV-INFER。已有A31_A3实际同网格保存93195自然候选，新推理须逐窗复现其几何、顺序和来源；这是可检验的同池要求，不以该数量预填新结果。模型前后逐tensor一致、无梯度／新增更新，窗口与最终GFF完整才允许COMPLETED。

## 完整评价与解释边界

CPU入口evaluate_context_chain.py先核对新推理实际终态、完整网格、源模型及源TRAIN身份，逐窗配对实际A31_A3记录，再以原严格FP32口径核对exon＋junction＝additive与additive＋nonadditive−null＝gain。全部保存链、装配计数、负链半开CDS坐标／phase及最终GFF重建关闭后，才打开原参考；空池窗口也必须保留，缺失部分记录不能变成较小评价。

沿用全部7728条primary CDS-chain参考、每物种原背景bp和生成器不支持的参考层，不通过筛选缩小分母。输出逐物种及pooled exact-chain precision／recall／F1、物种macro F1、最差物种F1、CDS碱基／背景错误、结构及长链诊断，并单列free→owner→selected覆盖损失；selected阶段TP须等于最终exact-chain TP。候选／窗口配对不当作去重基因或完整染色体F1。

十项比较为A37、A31_A3、A31_A23、A3、A6、A、R5_C0、ANNEVO、Tiberius、Helixer。A31_A3共享提议，A31_A23改变生成器，不能作为同池因果消融；A37还有额外容量、监督曝光和计算。既有native基线的来源限制继续保留，不与cached方法作公平原生耗时排名。完整DEV只判断固定模型是否缩小性能差距及后续机制方向，不以单seed、重复DEV、TRAIN proxy或工程PASS宣布SOTA、创新、独立泛化或投稿成功。本轮不打开独立物种、Setaria／sealed或A33待答复例外，也不公开上传未发表材料。

## 资源与停止规则

沿用实际generanno及原conda初始化，不使用旧cluster_config默认环境。拟请求private／1RTX3090 24GiB／2CPU／16GiB／7:30／no-requeue。实际A31同网格Slurm 3:13:47，加A37零更新七个热窗load＋forward＋backward均值1.630秒×8690约3:56:05，得到约7:10的保守操作估计，再留约20分钟余量。热窗含缓存I/O和反向、仅为八个既定metadata窗口的部分测量，原生候选负荷也不同；该加法不是原生基准或完成保证，不支持速度主张。原A31实际GPU峰值2455143424 bytes；新增评分器preflight峰值165758464 bytes，二者不冒称本轮原生峰值，24GiB为原实际执行资源。

CPU拟请求private／0GPU／2CPU／8GiB／5分钟／no-requeue；旧A31 A3／A23完整评价实际43／69秒，本轮增加逐窗配对和FP32分量核对，5分钟只作操作余量。新保存记录估计低于64MiB，不保存GEN隐藏或新feature cache。实际提交前仍检查活Slurm／预约、必要存储及唯一输出／日志；当前只准备，没有原生作业ID或执行结果。

五条件实际全PASS后，GPU命令为sbatch --parsable --dependency=afterok:13410742:13411198:13411914 sbatch/M28-A37-CONTEXT-CHAIN-DEV-INFER.sbatch；不能仅凭afterok（它不检查科学门槛）预先提交。CPU命令绑定实际新GPU ID的afterok，并以--export=ALL,A37_DEV_INFER_JOB=<实际ID>传同一源ID，调用sbatch/M28-A37-CONTEXT-CHAIN-DEV-EVAL.sbatch。SBATCH不预建入口要求不存在的输出目录，JOBID／COMPLETED只在入口成功后写入；未创建目录的失败保留Slurm及日志证据。

工程失败保留产物，不自动重提、取消、延长、续接或覆盖；若有具体实现错误，先按实际错误及既有授权判断，不改变冻结科学配方或把性能劣势当工程修复。完整DEV及评价闭合后主负责解释结果。长期等待仅复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同heartbeat coding-rna-r5-light-monitor；本轮准备不改当前fit→natural→CPU监控对象，也不同时新建等待循环。

## 实际准备记录

2026-10-06 16:05 UTC前读取，fit13410061的root／batch／extern均实际COMPLETED0:0／1:15:53，private1RTX3090／2CPU／16G／gpu034，batch MaxRSS2213712K。STATUS／JOBID匹配，summary记录原三遍各1536／总4608更新及末004608，training.jsonl实际4608行，三个checkpoint各7389798字节。这里只确认终态与结构，tensor数值、每参数曝光及全部原draw仍待原CPU13410742核验，未从训练loss宣布科研成功。

自然13410739实际RUNNING／5:56／gpu034，CPU13410742 PENDING／Dependency；没有读取部分自然分数、标签或门槛。新DEV合同及两个SBATCH已明确同步canonical，远端两脚本bash语法通过，原生DEV两个输出根仍不存在。当前未提交新原生GPU／评价CPU，未打开新DEV数据、参考或模型；同一固定监控仅继续原三根，终态或新故障后由主解释，不因fit已完成先放行DEV。
