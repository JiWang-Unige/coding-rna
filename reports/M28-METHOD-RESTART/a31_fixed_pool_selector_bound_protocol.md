# A31 固定池与非重叠选择器上限诊断

2026-10-06。用已闭合的 A6/A7 成员表量化 span 非重叠约束本身的损失，区分“继续改分数有理论空间”和“正确链即使理想评分仍被结构约束排除”。这是标签辅助的探索性上限分析，不是新模型注释、实际性能、独立测试或论文 PASS。

## 固定来源和范围

仅使用 `outputs/M28-A6-A7-SCORE-ATTRITION` 的7728条完整主参考成员表及summary，来源CPU13360116；同尺ANNEVO来自 `outputs/M28-NATIVE-BASELINE-EVAL`，来源13358432。两个实际 root 必须 COMPLETED/0:0、STATUS/JOBID对应。复算A6/A7两物种五状态计数，必须与已接受summary一致，保留全部unsupported参考和5437/2291主分母。

不读取正在运行的 A31输出、FASTA/catalog、DEV feature cache、真实权重或封存Setaria；不前向模型、搜索候选、训练、扫描阈值或写新GFF。原A6/A7结果和gate不变。

## 三种标签辅助上限

分别对完整primary参考、实际A6自由池中的正确参考、实际A6 owner池中的正确参考运行原 `select_nonoverlapping`。给每条正确链单位正gain，在species/seqid/strand内最大化可保留正确链数U；参考链完整span重叠视为冲突，CDS片段互不相交并不自动兼容，首尾相接可以共存，不同strand不互相竞争。

这是任意完美候选级标签评分下的几何上限，未证明实际15类局部分数能实现该赋分。实际自由/owner成员仅由已保存A6_state决定，不新增R4 label或资格筛选，不把所有参考注入自然候选池。完整参考那一行只用于量化选择器自身限制，不代表当前生成器可产生所有参考。

理想precision=1的F1上限为2U/(R+U)，R仍为完整5437/2291。与同尺ANNEVO实际TP、预测数及完整R用整数交叉乘法比较，并计算严格超过它的precision=1最少TP，以及仅在TP达到U时可容许的最大非精确预测数；这是全gene exact-chain FP，不是TRAIN known-label proxy FP。不计背景约束、耗时或曝光公平性。

A6 owner上限仅在A31_A3完成后的实际同池几何核验成立时约束该variant，不适用于不同池的A31_A23或A30。若上限仍低于native目标，单改该固定池分数不足以达到目标；若高于目标，仅排除这个数学阻碍，仍须等待实际完整F1和背景代价。不能以oracle超过ANNEVO宣称模型胜出或泛化。

## 验证和执行

四项纯合成检查覆盖：span而非CDS-overlap、半开边界及最大基数；strand/species隔离；严格native超过及整数FP预算；完整unsupported参考保留及唯一参考key。实际来源若出现同几何不同primary身份，停止并核对原评价语义，不静默缩小分母。

一次private CPU作业，0GPU/2CPU/8G/5min，沿用现有generanno环境；读取已有7728成员，运行小规模DP，预计远低于5分钟，输出仅summary和逐参考oracle成员表。提交前五项调度检查及唯一新目录检查。主对话直接消费这项短诊断，不加入当前七root固定监控，不创建新对话或自动化。结果记录另存，保留原实验所有产物。

实施回执：CPU13400401实际root/batch/extern COMPLETED/0:0/13秒，四项合成检查PASS/0.003秒，STATUS/JOBID匹配。完整7728成员和两物种原状态计数复现；[上限与研发判断](a31_fixed_pool_selector_bound_result.md)另存。诊断没有改变任何模型、门槛、训练曝光或监控对象。
