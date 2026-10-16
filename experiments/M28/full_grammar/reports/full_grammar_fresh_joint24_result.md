# 原 TRAIN fresh 抽样与 joint 编码的四格结果

2026-10-12 UTC。B的完整拟合、最终权重核验、六条raw-FASTA推理、共同评价与四格汇总均已完成并由主核验。新格pooled整链F1为70.302624%，较原joint提高3.339861个百分点，两个物种均改善；仍低于Tiberius、ANNEVO及默认原生Carbon-A。F1尺度交互量为−3.987713个百分点，结果不能称超加性互补、因果模块必要性或独立泛化，Nature Communications及以上的发表目标尚未达到。

## 完整执行与相同输入

遵循[先定方案](full_grammar_fresh_joint24_plan.md)：从原13434603末4608及AdamW状态独立出发，追加固定4608更新，只最终9216读出；fresh seed1的4608 draws／4219 unique、原fine/coarse/unknown观察、importance1、冻结b048 backbone和366191参数头不变，无native CDS先验。原两物种六TRAIN、22,688 primary参考及unsupported层全部保留。DEV、sealed、Setaria或新面板均未访问。

fit13494895的root/batch/extern于03:33:24UTC完成，CPU核验13494896于03:34:49完成；此前已核验完整4608行抽样／目标日志和三个checkpoint。六raw实际03:34:50同时启动，numeric roots为13508070／13508071／13508072／13508073／13508074／13494897（任务0–5）；全部root/batch/extern COMPLETED0:0，root结束03:59:08–04:13:11UTC。评价13494898于04:13:11–04:14:21完成，四格汇总13494899于04:14:22–04:15:13完成。整条10个实际roots及batch/extern共30条记录正常终结，不凭数组parent或squeue消失推定成功。

六外层收据、六内层预测收据、STATUS／JOBID、input FASTA及非空GFF齐全；全部指向同一step_009216.pt、source_fit_job13494895。原生整个定向contig双链skyline、FP64评分及原身份保持，caller不读参考／标签／feature cache、不更新参数。完成的共同评价核对GFF计数、DNA／phase／cap及同strand非重叠；四格CPU汇总额外直接比较全部四份caller FASTA与原DNA的大写／非ACGT→N规范化序列，六染色体逐碱基一致，且各格保留相同参考分母。主没有再在登录节点加载权重或运行重数据检查。

四个评价源固定为F00原draw分块13451625、F10 fresh分块13453583、F01原draw joint13453780、F11 fresh joint13494898；旧三格不重训、不重算预测。四格汇总同时保留背景与CDS-base端点，不择优筛掉负交互。

## 完整链计数与物种结果

| 条件 | TP | 唯一预测 | 参考 | Precision % | Recall % | F1 % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 原draw／分块编码 | 12696 | 22862 | 22688 | 55.533199 | 55.959097 | 55.745335 |
| fresh draw／分块编码 | 14140 | 22149 | 22688 | 63.840354 | 62.323695 | 63.072909 |
| 原draw／joint编码 | 14926 | 21892 | 22688 | 68.180157 | 65.788082 | 66.962764 |
| fresh draw／joint编码 | 15739 | 22087 | 22688 | 71.259112 | 69.371474 | 70.302624 |

| 条件 | Arabidopsis F1 % | Rice F1 % | Species macro F1 % | Worst species F1 % |
| --- | ---: | ---: | ---: | ---: |
| 原draw／分块编码 | 61.835611 | 35.412896 | 48.624254 | 35.412896 |
| fresh draw／分块编码 | 67.287503 | 48.869396 | 58.078449 | 48.869396 |
| 原draw／joint编码 | 70.613934 | 53.984471 | 62.299203 | 53.984471 |
| fresh draw／joint编码 | 73.367643 | 59.635781 | 66.501712 | 59.635781 |

F11较F01的primary TP净增813（Arabidopsis+475／rice+338），唯一预测净增195（Arabidopsis−11／rice+206）；pooled precision提高3.078955、recall提高3.583392个百分点。这里是计数净差，不是813条全部新恢复且没有旧命中丢失；四格程序没有运行逐参考交集重放。

| 染色体 | 新格TP | 唯一预测 | 参考 | 新格F1 % | 新格减原joint F1 pp |
| --- | ---: | ---: | ---: | ---: | ---: |
| NC_003070.9 | 5137 | 6943 | 7135 | 72.979116 | 2.793971 |
| NC_003071.7 | 3012 | 4023 | 4293 | 72.438672 | 2.837903 |
| NC_003076.8 | 4610 | 6076 | 6311 | 74.432873 | 2.650015 |
| NC_089043.1 | 1010 | 1649 | 1638 | 61.454214 | 5.878120 |
| NC_089044.1 | 988 | 1657 | 1601 | 60.650706 | 5.460991 |
| NC_089046.1 | 982 | 1739 | 1710 | 56.944042 | 5.646696 |

## 四格差值与交互

定义F00原draw分块、F10 fresh分块、F01原draw joint、F11 fresh joint。下表交互为100×(F11−F10−F01+F00)，F均为count-derived整链F1。

| 范围 | fresh效应分块 pp | fresh效应joint pp | 编码效应原draw pp | 编码效应fresh pp | 交互 pp |
| --- | ---: | ---: | ---: | ---: | ---: |
| arabidopsis_thaliana | 5.451892 | 2.753709 | 8.778323 | 6.080140 | -2.698183 |
| oryza_sativa | 13.456499 | 5.651311 | 18.571574 | 10.766386 | -7.805189 |
| pooled | 7.327574 | 3.339861 | 11.217429 | 7.229716 | -3.987713 |

两个单独改动在该面板均提高F1，联合格亦高于各单改格，但增幅小于F1尺度上的简单相加；species macro交互−5.251686、worst交互−7.805189个百分点。采样同时改变有限监督暴露与梯度序列，F1本身非线性且只有单次拟合；负交互不是统计显著的生物学拮抗、因果冗余或模块必要性证据。协议在既有DEV已读之后固定，不能称DEV前预注册。

## 强基线与不同端点

| 固定原生方法 | Arabidopsis F1 % | Rice F1 % | Pooled F1 % | F11减该方法 pooled pp |
| --- | ---: | ---: | ---: | ---: |
| ANNEVO | 76.830754 | 71.178600 | 75.541796 | -5.239171 |
| Tiberius | 75.417268 | 63.576543 | 72.655456 | -2.352832 |
| Helixer | 68.122369 | 50.294337 | 63.591926 | 6.710698 |

F11在两个物种均超过固定Helixer，但均低于Tiberius和ANNEVO。默认Carbon-A过滤／未过滤为78.171555%／76.063450%，F11分别低7.868931／5.760825个百分点；短窗口两点为77.562301%／75.585179%，差距仍为7.259677／5.282555个百分点。原生模型与辅助网络的标注暴露、head、评分、grammar和过滤不同，不能将共享主权重或名义窗口相同当作decoder-only因果控制。

F11的unstranded CDS-base F1为95.236061%，略低于F01的95.389874%；reference-relative background span覆盖从0.701594%升至1.009303%，wholly-background唯一链从320升至407。整链改善不等于全部端点或背景负担均改善。

新格另有1,316条唯一预测精确匹配其他可评价isoform（Arabidopsis1,202／rice114），不加算primary TP；剩余5,032条不精确匹配primary或其他可评价isoform，不据此判为生物学假基因。

## 下一阶段边界与来源

B已main-consumed，原心跳关闭其监控对象并保留C/D/E；未观察到B的监控终态发送，记为NOT_SENT_BY_WATCHER而非YES。不更换checkpoint、抽样、参数、阈值或进行新增拟合。冻结三候选DEV仍包含B/C/D全体，不因TRAIN排名增删；须待C/D各自完整fit／verify／六raw／评价／配对结果均由主消费后才接续。E仍为单独TRAIN探索，不自动加入DEV。原八／十方法历史结果不改写，封存及对外发布边界仍保持。

- [实际评价13494898](../../outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-TRAIN-DIAGNOSTIC/summary.json)
- [四格输入与算术收据13494899](../../outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FOUR-CELL/summary.json)
- [原始四格表](../../outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FOUR-CELL/comparison.md)
- [原生共同基线13461307](../../outputs/M28-NATIVE-TRAIN-BASELINE-EVAL/summary.json)
- [Carbon-A固定默认结果](carbon_a_original_train_result.md)
- [Carbon-A短窗口结果](carbon_a_context24_train_result.md)
- [冻结B/C/D DEV方案](full_grammar_joint_controls_dev_plan.md)

本机只保存本次收据、结果JSON和四格表；大体积FASTA、GFF、完整训练日志及权重继续保留在Baobab原根。
