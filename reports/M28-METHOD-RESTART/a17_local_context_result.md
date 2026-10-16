# A17 局部注意力生成器：原 gate 全部失败

2026-10-06。固定A17 fresh拟合13387910已 COMPLETED、ExitCode 0:0、Elapsed 1:22:02，完成原4608更新。自然生成13387946已保存完整1491窗自由池，但JSON链标签归一化错误使任务 FAILED；该工程失败不是性能结果。按 [标签恢复协议](a17_label_recovery_r2.md)，只重用原自然池恢复标签，不拟合、不前向、不搜索新候选。恢复13391371 root/batch/extern均 COMPLETED 0:0、19:18；只读CPU核验13391372三项同为 COMPLETED 0:0、14:51。两个R2阶段STATUS及JOBID与实际作业匹配。

## 实际配对覆盖

| 物种／生成器 | 全部自由候选 | correct | known_nonmatch | unknown | 至少5exon correct |
|---|---:|---:|---:|---:|---:|
| Arabidopsis A3 | 8909 | 1751 | 6305 | 853 | 604 |
| Arabidopsis A17 | 8464 | 1690 | 5890 | 884 | 580 |
| Rice A3 | 7412 | 629 | 6261 | 522 | 184 |
| Rice A17 | 7617 | 611 | 6364 | 642 | 171 |

Arabidopsis correct减少61，Rice减少18；pooled至少5exon correct从788减到751（减少37）。原三项严格增加gate全部FAIL，TRAIN_generator_supported=false。

CPU核验确认1491有序窗、16081个A17自由候选、原4608 draw／weight／A3 supervision身份、末模型649457参数及有限tensor、新旧自然池DNA标签、保存池身份及原三gate复算一致。没有reference注入、恢复阶段优化／新头或GLM前向、DEV或Setaria。原失败R1下游13387978仍未执行，未取消；它不是R2核验的成功证据。

## 判断

A17工程恢复完成、标签与统计可信，但生成机制筛选失败。不能把减少候选／Arabidopsis known_nonmatch减少视为整体成功，Rice known_nonmatch还增加103。关闭该候选，不调整seed、曝光、beam、loss、阈值或评价口径，不执行DEV、caller放行或发布。

这些是TRAIN内候选／窗口配对，不是选中基因FP、独立基因F1、泛化或局部context有害的因果证明。A25在本结果读取前已固定独立的initiation监督及attention加项，仍按自身原A3覆盖gate执行；它不重训或覆盖A17冻结结果。
