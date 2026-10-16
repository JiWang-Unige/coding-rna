# 遗漏正确链的固定发射分数诊断

2026-10-06，实际计算前固定。原A27已证明536个正确链／TRAIN窗配对没有进入原自然池，但“被生成”未必意味着当前局部分数会选择它。本诊断检验遗漏链的原固定gain与既有候选竞争，决定后续是否仍应优先改善搜索，或同时需要改进角色表征／目标。它不是新caller、自由预测或模型性能gate，也不恢复A33失败作业。

## 固定来源和总体

使用已闭合A27 13392411、A24 13391716及原A31 TRAIN读出13398521／独立CPU13398524。保持全部536遗漏（拟南芥421、水稻115；全部query局部排除315、至少一个局部通过221）和324个原metadata顺序命中对照，共860配对、396个原六条TRAIN染色体上的24576bp定向窗口。身份、顺序、参考几何、自然池和旧标签须与真实源记录匹配，不按分数、链长或当前A33结果选择样例。

仅使用原A3末004608／649457参数、原GEN特征cache、现有参考几何与已保存R4标签。参考用于评分诊断，不重新读取catalog或产生新R4标签；不访问DEV、Setaria、sealed、独立物种DNA／GFF、A33中间输出或新的GLM前向。没有训练、参数、梯度、候选搜索、权重选择或阈值扫描。

## 固定计算和解释

先复现每个所用窗口全部原A31_A3自然候选分数，使用原15类角色和单位权重、角色logit减intergenic logit、FP32模型／FP64求和及null0。再对全部原参考链用完全相同公式评分；原命中对照的分数须与自然池同几何成员匹配。

对每个遗漏配对，单独将该参考几何与固定分数加入隔离的诊断池，按原canonical几何顺序和原同strand span非重叠DP选择。参考插入只存在于新诊断输出，不写原自然池、GFF或实际caller。固定分为四种互斥结果：gain≤0；gain>0但仍未选；选中且未移除原已选正确链；选中但移除至少一条原已选正确链。记录移除的原正确／known_nonmatch／unknown，不能把primary竞争者默认为负例。

另在每个窗口一次性加入其全部遗漏参考几何，计算共同选择和原已选三类的变动。单链各自能选的数量不能相加成共同找回数；共同加入使用参考本身，仍不是模型自然生成或真实基因性能。unknown继续参加原DP，不删掉以改善精度；本396窗口子集不产生F1／recall／precision，也不与原7728基因完整DEV比较。

按两物种、原315／221局部机制和one／two_to_four／five_or_more类别保留分母与四状态。324对照另列，不混入遗漏总体。参考／窗重复不是独立gene样本；角色分数不能解释成校准概率或完整HMM likelihood。此诊断不证明新搜索可以自然找回、不会新增假阳性、独立泛化或机制因果。

GPU保存396份有效FP32 logits、860条references和396条windows。CPU仅从保存logits独立计算完整窗口log-probability差，核对旧自然分数和新参考分数、原身份／分类、单独及共同DP选择与实际保存汇总。沿用rtol1e-12／atol1e-8；不同代数求和若改变零点或选择判断则停止解释，不放宽容差或偷偷改分类。CPU不执行模型、backbone或重新打开原feature cache。GPU最终核对全部权重逐tensor未变且无梯度。

## 实际执行范围

入口scripts/experiments/M28-METHOD-RESTART/diagnose_missing_path_emission_train.py，score和check两阶段；四项小型检查覆盖原零点、全局区间竞争、正确／未知位移、单链与共同选择不可相加及命中对照身份。原角色／DP实现不变。

GPU输出outputs/M28-A34-MISSING-PATH-EMISSION-DIAGNOSTIC，private1RTX3090／2CPU／16GiB／20分钟、generanno；CPU输出outputs/M28-A34-MISSING-PATH-EMISSION-CHECK，private0GPU／2CPU／8GiB／10分钟、afterok:实际GPU根。日志位于outputs/M28-PILOT-R5/logs/M28A34EMISSION_JOBID.*及M28A34EMISSIONCPU_JOBID.*。396份logits预计约584MiB，文本低于16MiB，均使用新目录，不覆盖旧结果。

没有新query搜索，预期耗时低于原373窗／约24分钟的A33搜索分片；20分钟只是评分和cache IO的操作上限，不是速度实验或公平计算预算。运行故障保留实际输出；至多一次有证据、同科学范围和原资源时限的最小工程修复采用独立R1目录，再失败则停止，不修改旧A33停止条件。完成须两根实际成功、匹配STATUS／JOBID及上述全部真实产物和复算条件；提交不算结果。只剩外部等待时更新同一固定监控对话／同一心跳对象，不创建第二监控。

若多数遗漏gain≤0或正分仍输给原竞争池，搜索改善本身不足以解决这些已观察配对；若大量正分链在共同诊断池获选且原正确损失少，则支持继续检验自然搜索恢复，但不能替代A33五项固定条件。原失败、未知曝光、sealed限制和完整raw-FASTA竞争要求不变。无Git提交／推送、取消作业、公开发布或goal状态变更。

## 已提交的实际链

源13391716／13392411／13398521／13398524现均Slurm root COMPLETED／0:0；实际A27 references860，missing536／control324，missing机制315／221、同窗同几何无重复，原两输出新目录均不存在后提交。两Python文件本地和远程AST、两SBATCH语法通过；新GPU入口四项实际检查PASS／0.003秒。

| 阶段 | 实际Job | 已核实资源与依赖 | 当前已观察状态 |
| --- | --- | --- | --- |
| 冻结评分与参考辅助选择 | 13403916 | private1RTX3090／2CPU／16GiB／20min，gpu034 | 11:14:35 UTC RUNNING，用时3:38 |
| 保存logits独立CPU复算 | 13403920 | private0GPU／2CPU／8GiB／10min，afterok:13403916 | 同次取证PENDING，有效Dependency |

实际提交命令：

```bash
sbatch --parsable sbatch/M28-A34-MISSING-PATH-EMISSION.sbatch
sbatch --parsable --dependency=afterok:13403916 --export=ALL,A34_SCORE_JOB=13403916 sbatch/M28-A34-MISSING-PATH-EMISSION-CHECK.sbatch
```

实时五项Slurm检查无本任务并发／array／VRAM／维护冲突，private034实际4／8GPU占用、035空闲；BeeGFS实际156TiB可用。选private，不干预其他TE／用户作业，也不取消原A33和A30无效依赖。实际scontrol证明新CPU无GPU请求，依赖和时限匹配；GPU实际分配与请求相同。

仅更新原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor到本两根新cohort，原A33终态与native短预检全部排除；真实保存kind heartbeat／固定target／15分钟／failed_runs_only／ACTIVE和完整prompt逐字匹配。新manual turn01a110eb-3f2f-7e10-8754-9522bbf47601 completed／error-null／38.432秒，11:13:49.037Z取证GPU RUNNING、CPU有效Dependency；实际turn_context11:13:35.075Z仍gpt-5.6-luna／low及同项目。本cohort定时运行和整链终态交接尚未发生，不继承旧cohort验证；健康中间态静默。

本链现已完整闭合：GPU13403916／CPU13403920均root／batch／extern COMPLETED／0:0，分别6:35及2:55；实际860配对、396窗及保存logits齐全，独立复算通过。结果与科学解释见[a34_missing_path_emission_result.md](a34_missing_path_emission_result.md)。主对话直接消费后同一心跳PAUSED；不把manual启动冒充定时终态自动回传。A33原FAILED及其停止条件例外待答复不变，本诊断不部署其替换checker或重提它。
