# A31 评分与候选池 TRAIN 诊断结果

2026-10-06。同一个冻结A3和1491 TRAIN窗下，局部评分主要减少greedy池中的已知非匹配，并不增加正确链；在固定局部分数下扩展为greedy＋beam联合池，才增加正确与长链选择，同时也重新增加非匹配。两项都是给定另一项条件下的计算对照，不是独立因果贡献份额或完整二因素消融。原五项门槛及两个正在进行的DEV variant不变。

## 三格同模型比较

沿用原R4标签及拟南芥2172、水稻744个eligible窗口参考配对。P=正确/(正确+known_nonmatch)，R=正确/eligible配对，F1=2正确/(eligible配对+正确+known_nonmatch)。未知候选参加原零点非重叠选择，但不计已知假阳性，不能据此称真实gene precision或F1。表格用于核对相同几何子集和计数，所有格均为同一原末004608，参数和训练曝光未改。

| 物种 | 池与分数 | 正确 | known_nonmatch | 未知 | 至少5exon正确 | proxy P | proxy R | proxy F1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 拟南芥 | greedy＋query | 1135 | 490 | 147 | 299 | 69.85% | 52.26% | 59.78% |
| 拟南芥 | greedy＋local | 1125 | 251 | 141 | 301 | 81.76% | 51.80% | 63.42% |
| 拟南芥 | union＋local | 1531 | 347 | 315 | 548 | 81.52% | 70.49% | 75.60% |
| 水稻 | greedy＋query | 425 | 265 | 68 | 111 | 61.59% | 57.12% | 59.27% |
| 水稻 | greedy＋local | 425 | 170 | 60 | 111 | 71.43% | 57.12% | 63.48% |
| 水稻 | union＋local | 478 | 251 | 132 | 152 | 65.57% | 64.25% | 64.90% |

固定greedy池后，局部评分使拟南芥正确减少10、known_nonmatch减少239、未知减少6，长链正确增加2；proxy P增加11.91百分点而R减少0.46，F1增加3.63。水稻正确和长链正确均不变，known_nonmatch减少95、未知减少8；P增加9.83百分点、R不变、F1增加4.21。它不能补回greedy未生成的链，不能将这一格描述为提高正确链覆盖。

固定局部分数后扩池，拟南芥正确增加406、known_nonmatch增加96、未知增加174，长链正确增加247；P减少0.24百分点、R增加18.69、F1增加12.19。水稻正确增加53、known_nonmatch增加81、未知增加72，长链正确增加41；P减少5.86百分点、R增加7.12、F1只增加1.42。水稻的精度代价不能被合并或拟南芥收益掩盖，完整DEV需按物种保留全部参考与背景分母。

## 身份与结论边界

greedy几何逐链映射到已有自然联合池，直接取该几何的原局部分数，R4标签不变；unknown从未按标签删池。三格只改变明确的池或分数，零阈值、同一DP选择和同一训练模型不变。两个原格逐字段计数和全部2916参考分母精确复现已核验A31源，新中间格保存1491行。缺少union＋query格，原保存池未保留逐query的beam归属和存在分数，不能从greedy或其他评分猜补；不能估计交互，也不对A23作这项同生成器归因。

这是同模型、单seed、in-sample窗口配对的探索性分解，不证明生物学因果、独立泛化、新颖性或完整caller竞争性。它支持继续检查局部评分与扩大池共同组成的已冻结系统，而不是救阈值、改beam或改写原gate。两个完整raw-FASTA DEV作业已独立运行；本诊断未读其结果或用于模型选择。

## 第四格需要另定义候选级评分

缺格不只是补存日志的问题。原greedy基准在同几何的greedy生成query之间取最大存在logit；原beam则允许每query最多四条完成链，各链共用该query的存在logit。联合池只保留几何和greedy／beam4来源，丢弃query归属。将存在logit直接复制给beam链，会使同一query的存在支持在互不重叠备选上重复进入区间目标，不等于候选链质量。若改取全部greedy和beam生成query的最大logit，又可能重打已有greedy几何的分数，因而不再是只扩池的对照。以上是代码允许的评分语义差异，不是已测得的真实碰撞频率或新性能结果。

原beam累计path_log_score也不能直接补格：几何pointer及EOS的log增量均非正，而原非重叠选择仅接受相对null0的正gain；直接使用这些累计值会选空集。要获得可解释的第四格，须预先固定候选级gain或query与链的联合选择约束，并明确原greedy分数是否保持；这已超出“恢复原保存分数”，不能偷偷加入校准、长度归一化或改阈值。本轮保留三格和条件差，不提交简单复制query分数的联合池消融，不改变既有门槛或正在运行的DEV。

实现依据为run_local_path_caller.greedy_record的greedy-only最大logit、query_beam._beam_query的多完成链与共享query_logit、query_beam_pool.completed_union的几何去重及来源标签，以及core.select_nonoverlapping的严格正gain选择。这些代码事实约束归因设计；没有重新执行模型、搜索、标注或统计，本报告原计数和已完成作业证据均不变。

## 实际执行证据

命令sbatch --parsable sbatch/M28-A31-SCORE-POOL-TRAIN-DIAGNOSTIC.sbatch返回13399561，private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，generanno。两Python AST和SBATCH语法通过，四文件窄范围上传成功；提交前五项实时Slurm检查通过，无输出碰撞、相关维护或并发限制。root／batch／extern均COMPLETED0:0，实际47秒，batch MaxRSS718172K，STATUS／JOBID相符。

三项合成检查全部PASS，内部0.002秒：不同候选次序的几何映射、未知参与选择且不计已知假阳性、缺失几何或同几何标签冲突拒绝统计。真实计算仅读取已闭合的natural_scores.jsonl／windows.jsonl和原draw身份，没有模型或backbone forward、新R4标注、候选搜索、feature cache、optimizer更新、DEV、sealed或Setaria。实际summary和1491行中间格位于outputs/M28-A31-SCORE-POOL-TRAIN-DIAGNOSTIC；源为readout13398521／CPU13398524，命令diagnose_local_path_pool_score_train.py。该短作业由主对话消费，不加入原七根心跳。

## 原自然 TRAIN 池仍有正确链选择空间

2026-10-07。[标签辅助上限协议](a31_score_pool_train_diagnostic_protocol.md)只读已闭合A31_A3原记录。全部1491窗、16321自然链及原2380／12566／1375标签、2916／1086参考分母精确复现；不读取正在运行的A40或任何DEV结果。两个目标在该实际人口上给出相同长链上限，所有known1正确配对都能同时保留，正确链内部的span冲突损失为0。

| 物种 | A31实际正确 | 可共存正确上限 | 未利用正确配对 | A31实际至少5exon正确 | 长链上限 | 未利用长链配对 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | 1531 | 1751 | 220 | 548 | 604 | 56 |
| 水稻 | 478 | 629 | 151 | 152 | 184 | 32 |
| 合计 | 2009 | 2380 | 371 | 700 | 788 | 88 |

首个标签辅助解实际构造出两物种正确严格增加、known_nonmatch均0且不超过原347／251、pooled长链788严格超过700的五条件见证。因此，固定自然池与原非重叠选择器并未在数学上阻止该TRAIN筛选目标；不能把此前同池失败笼统解释成已知正确链彼此冲突。评分和非精确／unknown候选竞争仍有改善空间，但不由此判定哪个新表示或损失能学到该解，也不要求单改评分而排除改进候选生成。

oracle的非匹配与unknown为0是标签辅助赋分的结果：原数组保留全部候选，known0和protected−1在该诊断中赋0，不是原生模型删除unknown或已测得生物学precision1。此处371／88是重复窗口中的自然链配对差额，不是去重gene数；不将eligible参考数减候选数当作新漏检统计，不报告oracle为实际F1、独立泛化或A40 PASS。A40原五条件、失败停止规则及DEV限制均不改。

短CPU作业13419849的root／batch／extern均COMPLETED0:0，实际2026-10-06 22:01:02–22:01:29 UTC、27秒，private gpu034／0GPU／2CPU／8GiB／5分钟，batch MaxRSS725268K、Requeue0／Restarts0；Command／SubmitLine／WorkDir和实际分配与固定SBATCH相符，提交前五项活调度与独立输出检查通过。两个新增人工检查PASS0.004秒，stderr仅通过日志；STATUS COMPLETED／JOBID13419849、1491行windows和4827-byte summary相符，windows文件1573235 bytes。两个Python AST及SBATCH语法通过，四文件与canonical逐字一致；原计数、全部正确链可共存、两个长链上限一致、371／88差额及实际五条件由真实summary整数另行复算一致。

输出独立保存于outputs/M28-A31-NATURAL-TRAIN-SELECTOR-BOUND，不覆盖原分数／标签或任何冻结结果；没有新神经前向、checkpoint／feature访问、优化更新、新搜索或R4、参考catalog、Setaria或A33。短诊断由主消费，不加入原A40心跳或新建监控。以上支持探索性瓶颈解释，不是论文性能或科学放行。
