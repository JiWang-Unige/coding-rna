# A37 自然候选读出与训练曝光核验

2026-10-06，在新性能读出前固定。[原联合上下文协议](a37_context_chain_protocol.md)的全部科学选择、五项继续条件和停止条件保持不变。训练13410061运行中；本流程只使用实际完成的末004608，不读中间权重效果，不修改训练或开放DEV、独立物种、Setaria、A33例外或公开发布。

## 原同池比较与评分

原1491唯一TRAIN窗／1536 draw，完整16321个A3自然候选，几何、顺序和greedy／beam4来源不变，无参考注入或新搜索。源基线为已核验A31_A3读出13398521／CPU13398524，标签1／0／-1和2916 eligible reference/window分母及1086至少5-exon分母保持原口径。预测保留未知竞争，0不是生物学真负，-1不并入错误或新基因；数量是候选／窗口配对，不是去重基因或完整染色体F1。

FP32新整链gain为exon additive＋junction additive＋整链nonadditive−窗口null。保持同一零gain操作点和span非重叠DP，不扫描阈值、重标或修改分数单位。本模型没有15类角色发射，不把局部log-probability差核验移植成新证据。

GPU在没有读取saved label数组／catalog的情况下，逐窗使用原GEN缓存及DNA和末模型，保存全部自然链、来源、gain、四个分量数组、scalar null和选择；全部natural_scores.jsonl关闭后才连接原A31已经核验的标签。保存第二个windows.jsonl和summary，原新模型前后逐tensor一致、参数全部冻结／无gradient／optimizer0。没有新GEN前向、缓存、logits矩阵、R4或参考序列访问。

五条件分别是两物种已选正确严格高于A31_A3、两物种known_nonmatch不增加、pooled至少5-exon已选正确严格增加；基线正确1531／478、known_nonmatch347／251、长链正确548＋152＝700，不退回A6或重判旧FAIL。任一FAIL关闭A37，不改loss、阈值、权重、checkpoint或曝光；PASS只支持另行真实raw-FASTA协议，不等于SOTA、因果机制、独立泛化或投稿成功。

## 准备检查

短CPU先运行四项针对新实现的检查：FP32相消和零gain、完整分量及空池；未知仍竞争与原候选身份／DP保护；逐物种较强基线与严格五条件；实际人工optimizer state与候选调用推算的分支曝光一致。最后一项使用小型随机ContextChain的4个人工case×2遍＝8步，仅验证optimizer行为，不是原科研模型训练，未读真实权重、GEN或特征。

原人口metadata准备核对全部1491窗、16321自然链、原身份／标签／DP／较强基线、2916／1086分母及feature-index覆盖。该metadata阶段真实模型／特征tensor／optimizer0；人工fixture更新与这一原数据阶段分开，不把合成检查当作新模型性能。准备入口readout_context_chain_train.py --mode checks、唯一输出M28-A37-CONTEXT-CHAIN-READOUT-CHECKS，只有实际完成才提交下游。

## CPU 完整核验

读真实4608训练行，逐条对齐三遍原draw、importance、物种、候选标签counts和有限update记录；zero-known draw保留且loss／partition／gold／gradient norm均为零。核对三checkpoint的完整有限FP32模型tensor、固定优化器参数及实际每参数step。SharedContext及null每个draw都调用；exon／链GRU／additive／nonadditive只在非空池调用；junction两分支只在池中实际存在多exon候选时调用，包括未知多exon链，不能按known连接数量或一律4608推定。没有被调用的参数应无optimizer state；被调用而零gradient的参数仍可能执行AdamW step和weight decay。

全部1491保存记录与原自然几何／来源逐条相同，保存数必须确为FP32可表示值。CPU用FP32分别计算exon＋junction＝additive及additive＋nonadditive−null＝gain，逐值严格相等，不放宽容差或用FP64改变原相消。再从原saved labels重新运行DP、计算逐物种counts、原分母／较强基线和五条件，并与GPU记录逐项一致。

该核验没有新模型／backbone前向或optimizer，不是独立神经前向复现，也不是独立第三方评阅。四个actual核验标记分别覆盖原draw／4608更新、三模型tensor及分支optimizer曝光、自然身份及FP32分量、原labels／DP／基线／分母及门槛。工程终态、数值可信和科学五条件分别判断，TRAIN_preparation_supported=false不是工程失败。

## 调度与停止

全部generanno／Slurm／no-requeue、原排除节点及独立%x_%j日志。准备0GPU／2CPU／8GiB／5分钟；自然读出1RTX3090／2CPU／16GiB／40分钟；CPU核验0GPU／2CPU／8GiB／5分钟。原A36单SharedContext实际自然读出21:15，A37增加有序头但不保存发射矩阵；零更新热窗口平均1.630秒包含反向，故40分钟是操作余量，不是A37实测速度或承诺。CPU只处理原metadata、三个小checkpoint和保存分量，原A21核验实际23秒；本轮5分钟为运行余量。新增JSON记录预计低于32MiB，不建特征／隐藏矩阵缓存。

准备实际闭合后，GPU afterok绑定原fit13410061及准备真实ID，显式传A37_CHECKS_JOB；CPU afterok绑定GPU真实ID并显式传A37_READOUT_JOB。源码同时要求实际sacct终态和匹配STATUS／JOBID／summary，不能读取部分训练或分数。

工程失败保留全部原产物，不自动重提、取消、延长或覆盖；主按具体错误及既有失败政策处理，不把科学FAIL当修复。长期等待只复用固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3和同heartbeat coding-rna-r5-light-monitor，提交下游后更新为fit→自然读出→CPU整链，正常RUNNING／有效Dependency静默，整链闭合或新可行动故障一次交接。主负责结果验证及已授权后续，论文总目标未完成。

## 实际执行记录

准备13410681的root／batch／extern已实际COMPLETED／0:0／0:31，STATUS／JOBID／summary匹配；实际private0GPU／2CPU／8GiB、gpu034，batch MaxRSS740904K。四项检查实际PASS／0.359秒，包括小型随机模型人工8步的分支曝光；原数据metadata阶段1491窗／1536 draw／16321自然链／1491 feature-index及身份／标签／DP／较强基线／2916与1086分母对齐，未打开新fit权重或特征tensor、未前向或更新原科研模型。

提交准备前及下游前五项实际调度／reservation检查无相关维护或资源阻塞，实际可用存储156TiB；未解读部分训练性能。

GPU13410739已提交一次，真实SubmitLine为sbatch --parsable --dependency=afterok:13410061:13410681 --export=ALL,A37_CHECKS_JOB=13410681 sbatch/M28-A37-CONTEXT-CHAIN-TRAIN-READOUT.sbatch。准备依赖已满足，Slurm保留afterok:13410061；真实PENDING／Dependency／请求1RTX3090＋2CPU＋16GiB／40分钟。CPU13410742已提交一次，真实SubmitLine为sbatch --parsable --dependency=afterok:13410739 --export=ALL,A37_READOUT_JOB=13410739 sbatch/M28-A37-CONTEXT-CHAIN-RESULT-CHECK.sbatch；真实PENDING／Dependency／请求0GPU＋2CPU＋8GiB／5分钟。两项均Requeue0／Restarts0及canonical command匹配；这些是请求与依赖，不是实际下游分配或完成。

整链读取时原fit13410061仍RUNNING／0:0／30:02／gpu034，下游均PENDING／0:0／正常Dependency。同一固定heartbeat已原生更新为整链cohort M28-A37-CONTEXT-CHAIN-NATURAL-20261006，完整保存prompt／kind heartbeat／ACTIVE／原target／15分钟／通知字段及created_at均实际读回一致。保留gpt-5.6-luna／low，不新建监控。

本cohort手动启动核对已实际完成：turn01a111c7-e5d9-75b2-9955-c6837e26135b／error-null／26.509秒，真实状态证据2026-10-06T15:14:53Z；首次定时核对也已实际完成：同一原对话turn01a111c8-7b06-7191-985c-c747b05bd510／error-null／32.310秒，真实证据2026-10-06T15:15:28Z，两次均fit RUNNING／gpu034及下游正常Dependency，无健康状态向主交接。实际最新turn_context 2026-10-06T15:15:13.843Z确认gpt-5.6-luna／low／原cwd；未将旧fit-only验证算作本cohort验证。

配置读回、启动核对、定时核对均已分别取证，但本cohort终态自动回传尚未发生。工程读出、CPU数值核验和五项科学条件均尚未完成，论文总目标仍在推进。
