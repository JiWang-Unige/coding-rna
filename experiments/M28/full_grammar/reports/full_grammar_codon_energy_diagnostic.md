# 最终权重 TRAIN 的神经角色与密码子评分诊断

2026-10-12。方案在fresh×joint、uniform-coarse与codon三条fit仍运行、尚无这三条新分支最终准确率时确定。诊断比较固定final9216的fresh×joint与codon能量模型，回答参考链和实际所选竞争链的分数由哪些项构成；不新增训练、重跑decoder、选择checkpoint或改动冻结DEV方案。

## 模型与输入

对照为fresh×joint fit13494895、TRAIN评价13494898；候选为codon fit13497914、TRAIN评价13497946。须先完成并消费这两条原六TRAIN流程及codon配对13497947，核实全部实际Slurm终态和产物后才提交真实诊断。两个head最终年龄均9216，codon新参数年龄4608；不读取6144／7680快照。仅用同一原六染色体raw输入、冻结GFF、两套严格末权重和原22,688条primary参考profile，保留不支持参考。

参考不参与编码、ownership、模型评分参数或路径选择。重放沿用joint24576、stride12288、nearest-center、两条定向序列、15角色、cap128和同一字面DNA语法。不保存整染色体势函数、大特征、DNA副本或新GFF，不访问DEV、held-out、sealed、Setaria或新物种。

## 分数与解释

对每个支持链，按原角色路径保存15角色gain与base-row数。独立拼接其实际CDS，统计64类完整密码子；每个完整triplet计一次，包含起始ATG和末端stop。将同一最终去均值系数加在第三个编码碱基处，以FP64逐base次序累计role-plus-codon总分。对照也保存相同codon count，但其codon系数固定为零。计数是链曝光，不是unique genomic coverage或估计的codon usage。

每条冻结所选链的重放须落在原八位有效数字GFF舍入界加1e-5内；每strand所选总分与原raw收据须relative1e-8／absolute1e-5一致。参考的删除集合是同strand、span正重叠的全部所选链；相接不冲突，反strand不冲突。竞争者按整个链计分，不截取重叠部分。参考总分减竞争集合总分为合法单链替换margin，保存其15角色和64密码子差值。分解须重构总分、竞争分及margin。

额外保存同一头、同一参考与同一冻结竞争集合的role-only margin。这是去掉密码子项的固定集合评分算术，不是重新求解MAP，也不是zero-codon准确率消融。full margin不超过1e-5且role-only margin大于1e-5，说明密码子项在这一固定替换的数值不等式中起作用；不等于该参考将在其他MAP中被恢复，不证明生物机制因果。仅full-score合法替换能严格改善冻结集合时，停止科学解释并核对评分／坐标／decoder一致性。1e-5只用于数值与近零报告，不是新模型操作点。

两模型按物种、seqid、strand及全部CDS坐标配对，保存双方命中、仅对照命中、仅候选命中、双方遗漏四组，报告全部染色体、逐物种及pooled的类别迁移、角色margin、codon margin和曝光。净新命中必须重构原评价TP差，不预设其符号或大小。不同head的绝对gain不作概率或confidence比较；hit组属于事后描述，不产生subset precision／F1或独立样本。

## 实施与验证

新增独立诊断脚本、七组人工DNA检查及单一分模式SBATCH，复用已完成诊断的真实语法支持、角色路径、冲突与共同尺函数，不改活动训练或生产caller。人工检查覆盖跨多内含子和短外显子的codon闭合、三种stop与生产评分一致、零参数复现role-only与负链几何、整个多竞争集合和空集合的分解、错误和／codon计数／margin拒绝、四hit组及unsupported保留，以及固定集合去项不得改变命中计数。检查不读取研究DNA、参考、checkpoint或native模型。

本轮只提交checks：private-teodoro-gpu，无GPU，2CPU／4GiB／5分钟、Requeue0。真实诊断尚不提交。后续六TRAIN各一个GPU任务，顺序重放两模型，单24GiB RTX3090／2CPU／16GiB／90分钟；CPU汇总2CPU／8GiB／10分钟，afterok依赖整个数组。资源依据既有同更新数两头重放最长约58分钟及汇总48秒，密码子计数会额外增加轻量算术；提交前仍核对当时队列、维护窗口、存储和输出碰撞。预计产物不超过512MiB，不覆盖历史结果。

产物独立保存于`outputs/M28-FULL-GRAMMAR-CODON-ENERGY-CHECKS`和`outputs/M28-FULL-GRAMMAR-CODON-ENERGY-DIAGNOSTIC`；每真实task保存两臂逐参考记录、评分收据、started、STATUS和JOBID，CPU保存完整配对summary。人工检查通过仅说明新增算术可执行；真实末权重评分一致性、准确率收益及论文goal均另待证据。后续长等待只复用原监控对话与唯一心跳，不创建重复监控。

## 已完成的人工检查

实际提交命令为`sbatch --parsable sbatch/M28-FULL-GRAMMAR-CODON-ENERGY.sbatch checks`，返回13499288。提交前扩展后的用户作业数49，低于配置80上限，检查输出与真实诊断目录均不存在，脚本语法有效；无本次维护冲突。实际private-teodoro-gpu／gpu034，0GPU、2CPU／4GiB／5分钟、Requeue0。

UTC2026-10-11 18:12:59–18:13:22执行，root／batch／extern均COMPLETED、ExitCode0:0；batch MaxRSS949724KiB。七组实际unittest全部通过，0error／0failure，unittest耗时7.400秒、脚本7.464秒。STATUS、JOBID和summary均匹配；stderr为正常的七项测试名及OK，不是空文件，也没有异常。research_data_read、checkpoint_read、native_forward、DEV_used、sealed_data_used均false，更新0。

主已核对终态与记录，真实诊断输出目录尚不存在，本轮没有提交GPU重放或CPU全范围汇总，也没有改动三条活动训练或旧评价结果。完整程序的新算术已验证；真实末权重接入、全染色体重放及解释仍未执行，不将人工fixture通过当作正式科研结果或论文成功。
