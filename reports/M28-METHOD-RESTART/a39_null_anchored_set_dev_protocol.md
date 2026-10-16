# 绝对空集与兼容集合模型的原生 DEV 注释条件方案

2026-10-07北京时间，在A39原生DEV访问前固定。当前仅准备接口和合成检查；fit13415276正在运行，natural13415309及CPU13415319等待依赖，完整末模型结果尚未产生。只有[原TRAIN协议](a39_null_anchored_set_protocol.md)的三个实际root闭合、五项CPU证明及原五条件全部成立后，才执行本方案。A38、A37、A36的科学FAIL仍关闭。

## 执行前提与模型身份

主核对实际4608更新、001536／003072／004608三个checkpoint、原1491窗／16321自然候选和真实CPU证明，再重算较强A31_A3五条件：拟南芥／水稻正确配对分别严格>1531／478，known_nonmatch分别≤347／251，合并至少5-exon正确严格>700。原2916／1086参考窗口分母不变；只超过较弱A37／A38不够。任一FAIL关闭A39，不运行原生DEV，不改系数、标签、阈值、checkpoint、曝光或种子。

infer_null_anchored_set.py绑定fit13415276／natural13415309／CPU13415319，readout准备13415166、preflight13415165和末step004608；在模型、GEN、FASTA或新目录访问前重新要求源身份、完整数值证明与五条件。CPU除原draw、三checkpoint／分支、自然FP32分量及labels／DP／分母四项证明外，必须具有actual_saved_fixed_composite_loss_and_importance_arithmetic_verified=true。配置必须是固定等权mean-known BCE／兼容known-set NLL，两项都按原known数归一化、unknown不进入监督；A38逐候选BCE或A37单集合目标不能替代。

生成器仍为A3末outputs/M28-A3-FIRST-FIT/step_004608.pt、649457参数；评分器为outputs/M28-A39-NULL-ANCHORED-SET-FIT/step_004608.pt、612324参数。两者均原4608更新、eval及无梯度／新增optimizer。推理只用一个A39评分模型，不平均A37／A38模型或分数，也不在推理时计算复合训练loss。

## 原始输入与完整注释

执行仅限原--dev：拟南芥NC_003074.8共23459830bp、水稻NC_089041.1共29936421bp，原data/m1_screen物种genome.fa。24576bp宽、12288bp stride、固定尾窗及双方向，共4scope／8690窗。GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、原四块6144bp、BF16 materialization转FP32、非ACGT转N、方向化和valid-length保持不变；不访问独立物种、Setaria／sealed或A33例外。

每方向窗只编码一次GEN，原A3提议和A39评分复用同一特征与DNA，各有自己的上下文。自然提议仍greedy＋beam4／completed4／cap128及原ORF／内含子／相位规则。A39 gain仍exon additive＋junction additive＋nonadditive−null，FP32、操作点0；旧query／segmentation分数不进入gain。没有reference、标签、特征cache、参考注入、候选补正或标签筛选；全部自然链及来源、非正gain、四分量及null均保存，owner／几何去重／span DP不变。

唯一outputs/M28-A39-NULL-ANCHORED-SET-DEV-INFER含4方向JSONL、完整summary和A39.gff3。自然候选必须逐窗匹配已完成A31_A3同网格几何、顺序和来源；原93195候选是待核对的对照，不预填新结果。全部窗口及GFF完整、模型权重前后相同，才允许COMPLETED。接口支持显式--fasta，但本方案只执行原--dev，不自动换数据。

## 完整评价与比较口径

evaluate_null_anchored_set.py先核对原生作业终态、两模型／TRAIN身份和8690窗完整网格，再配对实际A31_A3记录，核对保存FP32分量并重建装配／GFF。全部预测、空池窗口、源记录与最终CDS坐标／相位闭合后，才打开原参考；不按候选可生成性或缺失预测缩小评价。

沿用全部7728条primary CDS-chain参考、原物种背景bp及alternative／unsupported层；报告逐物种、pooled exact-chain precision／recall／F1、macro及worst-species F1、CDS碱基／背景错误、长链／结构诊断和free→owner→selected损失。十项比较为A39、A31_A3、A31_A23、A3、A6、A、R5_C0、ANNEVO、Tiberius、Helixer。A37／A38没有完整原生结果，不用其TRAIN配对数填充DEV表。

A31_A3同提议，但A39增加上下文／头容量及监督曝光；A31_A23改变生成器。A37／A38／A39的匹配结构TRAIN比较与这种DEV对照不同，不据此作普遍因果主张。原native基线的上游监督／版本来源限制保留，不和cached方法作公平原生耗时排名。这仍为单seed、重复两物种DEV，不是独立泛化、发表成功或复合loss首创证据；不放行新物种、Setaria、claim或公开分享。

## 资源与停止条件

当前准备检查为generanno／private0GPU／2CPU／8G／5分钟／no-requeue，独立M28-A39-NULL-ANCHORED-SET-RAW-CHECKS，只用人工future-PASS收据及伪输入。五项新接线检查覆盖A39源模型与复合目标身份、不能以缓存true或较弱臂绕过原五条件、缺少新实际CPU算术证明时在原生访问前拒绝、评价源身份及新record／score kind。复用已有坐标、搜索、GFF和完整分母实现，不重复未变的人工测试。人工PASS不代表实际TRAIN或科学PASS。

原生GPU暂定private1RTX3090 24GiB／2CPU／16G／8小时／no-requeue。原A31同网格实际3:13:47；A39 preflight去首窗后七个实际load／forward／backward计时共12.480454秒，均值1.782922秒，外推8690窗约4:18:14，简单相加约7:32:01，再留约28分钟操作余量。该preflight包括缓存I/O、BCE及集合目标各自导数与复合反向核验，而原生只做一次前向，样本和候选负荷也不同；这不是原生benchmark或A39速度退化证据。实际自然读出闭合后，可据纯前向计时和活调度调整资源请求并记录，不改变科学命令、模型或阈值。

CPU完整评价为private0GPU／2CPU／8G／5分钟／no-requeue，沿用原A31评价43／69秒及新保存分量核对余量。新JSONL／GFF预计<64MiB、不缓存GEN。实际执行前检查活调度／预约、必要存储及唯一输出；当前只提交合成检查，不提交原生GPU或评价CPU，不打开实际DEV模型、DNA、参考或新监督权重。

五条件全部PASS且当前工程检查实际闭合后，主才提交DEV-INFER；不能预先把TRAIN CPU的afterok当科学放行。评价CPU随后绑定同一真实原生GPU ID的afterok和A39_DEV_INFER_JOB。SBATCH不预建推理入口要求不存在的目录，只在实际成功后写JOBID／COMPLETED。任何工程失败保留状态和日志，不自动重提、延长、续接、取消或覆盖；科学负结果不是工程错误。

本次短准备由主消费，不改变当前13415276→13415309→13415319的固定heartbeat对象。以后若实际原生阶段启动，只复用同一监控对话和同heartbeat；没有新对话或重复等待循环。论文与独立验证总目标仍未完成。

## 本次工程检查与交付记录

2026-10-07北京时间，短检查13415841实际root／batch／extern均COMPLETED、0:0，private／gpu034／0GPU／2CPU／8G，实际18秒、batch MaxRSS752908K，Restarts0；输出STATUS=COMPLETED、JOBID=13415841，五项人工接线测试全部通过（unittest 0.464秒）。只测试人工收据和合成输入，没有真实模型、GEN、DEV或参考访问，也不是五条件科学PASS。提交命令为sbatch --parsable sbatch/M28-A39-NULL-ANCHORED-SET-RAW-CHECKS.sbatch。

源码核对另发现DEV_PROTOCOL误写为不存在的a38_null_anchored_set_dev_protocol.md；在该短作业终态后修正为本A39文档。上述人工测试使用入口常量构造收据，未检出这一文件名错误；保留原检查日志，不将其声称为改正后重跑。改正仅一个协议路径字面量，另核对Python语法、实际常量与目标文档存在，以及本机／集群改动文本一致，不新增重提或科学运行。原生DEV-INFER／DEV-EVAL仍未提交，当前TRAIN三作业和冻结规则未改。
