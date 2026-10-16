# A31 固定评分与候选池 DEV 诊断结果

2026-10-06。同一冻结A3、全部8690窗的保存greedy几何匹配成功。沿事前固定三格路径，局部评分在greedy池上主要减少非匹配primary预测，扩为联合池则补回大量完整链和多CDS片段正确链。合并完整链F1依次52.3200%→56.4037%→66.1671%；这是保存输出的有序条件对照，不是独立因果贡献份额，也不以反复使用DEV选择新的模型或checkpoint。

## 实际闭合

按[a31_score_pool_dev_diagnostic_protocol.md](a31_score_pool_dev_diagnostic_protocol.md)执行。命令sbatch --parsable sbatch/M28-A31-SCORE-POOL-DEV-DIAGNOSTIC.sbatch返回13402574，private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，generanno。实际2026-10-06 09:31:19至09:32:09 UTC，root／batch／extern均COMPLETED／0:0，50秒，batch MaxRSS898120K。STATUS／JOBID一致，三项合成检查全部PASS，内部0.047秒。

两Python AST及SBATCH语法通过，四文件窄范围上传完成；提交前五项实际调度检查通过，private gpu034为mix、gpu035为idle，无相关维护或并发／存储／输出冲突。五个原source root全部COMPLETED／0:0：旧A3推理13357568／评价13357569、新A31 A3推理13399434／评价13399454、保存共同参考13360116。原两个格完整链整数、每scope全部stage及装配计数、旧query计数复现；全部8690窗的旧greedy去重几何恰等于联合池的greedy来源集合。原重复query记录保留，由原owner及重复取最大规则处理；没有按标签删候选。

先完成完整候选匹配和参考无关装配，再打开既有7728条完整primary参考。greedy_local_scores.jsonl实际8690行，diagnostic_selections.jsonl实际四scope行。无模型／GEN前向、新搜索、feature cache、R4标注、训练更新、阈值扫描、Setaria或生产GFF；新增中间格重组选择是诊断结果，不冒称新的原生推理。主对话已独立复算每物种与合并P／R／F1、stage边际整数、至少五片段正确计数及全部条件差，一致。

## 三格完整链测量

分母拟南芥5437、水稻2291，合并7728，unsupported及难例未删。TP为同strand完整有序CDS坐标链精确匹配；“非匹配primary预测”只表示未匹配此共同primary集合，不自动等于生物学错误。至少五片段按CDS坐标区间数分层，不当作基因bp长度。

| 范围 | 固定格 | TP／预测／参考 | Precision % | Recall % | F1 % | ≥5 CDS片段正确链 |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 拟南芥 | greedy＋query | 2519／4297／5437 | 58.6223 | 46.3307 | 51.7567 | 613 |
| 拟南芥 | greedy＋local | 2488／3539／5437 | 70.3023 | 45.7605 | 55.4367 | 620 |
| 拟南芥 | union＋local | 3524／4998／5437 | 70.5082 | 64.8152 | 67.5419 | 1203 |
| 水稻 | greedy＋query | 1219／2264／2291 | 53.8428 | 53.2082 | 53.5236 | 359 |
| 水稻 | greedy＋local | 1240／1952／2291 | 63.5246 | 54.1248 | 58.4492 | 369 |
| 水稻 | union＋local | 1465／2354／2291 | 62.2345 | 63.9459 | 63.0786 | 491 |
| 合并 | greedy＋query | 3738／6561／7728 | 56.9730 | 48.3696 | 52.3200 | 972 |
| 合并 | greedy＋local | 3728／5491／7728 | 67.8929 | 48.2402 | 56.4037 | 989 |
| 合并 | union＋local | 4989／7352／7728 | 67.8591 | 64.5575 | 66.1671 | 1694 |

固定greedy池改评分后，拟南芥TP减少31、预测减少758，非匹配primary预测减少727，≥5片段正确链增加7；P提高11.6801个百分点、R下降0.5702、F1提高3.6800。水稻TP增加21、预测减少312，非匹配预测减少333，≥5片段正确链增加10；P提高9.6818、R提高0.9166、F1提高4.9256个百分点。合并TP净减少10，预测减少1070，非匹配预测2823→1763、减少1060，P提高10.9199、R下降0.1294、F1提高4.0837个百分点。不能把这一步描述为两个物种正确链覆盖都提高。

固定局部分数扩池后，拟南芥TP增加1036、预测增加1459，非匹配预测增加423，≥5片段正确链增加583；F1提高12.1052个百分点。水稻TP增加225、预测增加402，非匹配预测增加177，≥5片段正确链增加122；F1提高4.6294个百分点。合并TP增加1261、预测增加1861、非匹配预测增加600，≥5片段正确链989→1694、增加705，R提高16.3173、P略下降0.0338、F1提高9.7634个百分点。扩池并非只增加正确链或所有指标均更好。

## 参考覆盖变化

以下均为同一7728条参考上的命中数。评分变化不改变free／owner成员关系；联合池扩展还可能补回此前仅见于非owner窗的链，不把owner改善全部称为新free链生成。

| 固定格 | Free参考命中 | Owner参考命中 | Selected参考命中 |
| --- | ---: | ---: | ---: |
| greedy＋query | 4443 | 3898 | 3738 |
| greedy＋local | 4443 | 3898 | 3728 |
| union＋local | 5967 | 5792 | 4989 |

greedy→union的free命中增加1524，owner命中增加1894，selected较greedy-local增加1261。结合[完整DEV及五状态诊断](a31_local_path_caller_dev_result.md)，下一机制准备优先针对完整多片段链生成／保留和正分候选竞争，而不是重复旧质量头、DEV调阈值或救A30失败门槛。本轮未实现或提交新的生成器训练。

## 解释边界与产物

固定A3末outputs/M28-A3-FIRST-FIT/step_004608.pt、null0、局部分数配方、owner／重复取最大／DP均不变。原A3是缓存推理，A31是原生raw-FASTA；完整几何相同不证明编码hidden逐数相同，所以评分对照不是纯同hidden的因果效应。三格缺少union＋query格，不估计交互或贡献百分比，不将有序差相加包装为独立机制份额。不复制query存在logit给beam备选或猜补缺失分数，理由沿用[TRAIN三格诊断](a31_score_pool_train_diagnostic_result.md)。

只解释一个A3生成器的重复使用DEV，没有对A23进行同生成器归因、统计显著性、匹配容量／训练曝光或独立泛化实验。全模型竞争性仍以[完整DEV结果](a31_local_path_caller_dev_result.md)为准：A31 A23虽高于Helixer，仍低于Tiberius／ANNEVO；本诊断没有改变原TRAIN／pilot门槛或论文通过判断。

canonical产物：/home/users/j/jwang/coding-rna/outputs/M28-A31-SCORE-POOL-DEV-DIAGNOSTIC/{summary.json,greedy_local_scores.jsonl,diagnostic_selections.jsonl,JOBID,STATUS}；实现diagnose_local_path_pool_score_dev.py，日志outputs/M28-PILOT-R5/logs/M28A31SCOREDEV_13402574.{out,err}。短CPU作业由主对话直接消费，不启用已消费A31的PAUSED心跳，不新增监控对话。未Git提交、推送或公开发布；保存Markdown未作渲染QA。

