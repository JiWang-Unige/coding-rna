# A31 固定评分与候选池 TRAIN 诊断

2026-10-06。A31主对照同时改变池和评分，不能直接归因于局部评分。两variant的原五项TRAIN门槛已通过、完整DEV正在独立执行；本辅助诊断只读已闭合TRAIN保存记录，不再拟合模型，也不依据结果选择variant、重判门槛或修改DEV协议。

固定同一个A3末004608和全部1491 TRAIN窗，复算三格：greedy池＋原query分数、同一greedy池＋已有局部分数、完整联合池＋同一局部分数。greedy几何通过显式坐标映射到原自然联合池，保持同一R4标签和零阈值非重叠选择。未知照常参与选择，独立计数，不改known-label proxy的假阳性口径。两原格计数及全部2916 eligible参考分母必须精确复现已核验源，再解释新增中间格。

第一项差是固定greedy池下的评分变化，第二项是固定局部分数下的池扩展。按物种报告正确／known_nonmatch／未知／至少5exon正确变化以及proxy P／R／F1百分点差；这些沿路径的差可在数值上相接，但不是独立因果贡献份额。缺少联合池＋query分数格，因为原保存池没有逐query的beam归属和query存在分数，不能猜补；因此不估计交互或宣称完整二因素消融。

源为TRAIN readout13398521及实际CPU核验13398524，配置／协议和保存顺序不变。只读natural_scores.jsonl、windows.jsonl及原draw身份，不载入模型、DNA缓存或参考catalog，不做新的R4标注、搜索、GLM编码或optimizer更新，不读DEV、sealed或Setaria。本诊断无新的科学继续门槛，原科学PASS／FAIL均不因它改写。

先三项合成检查：不同候选次序仍正确映射局部分数和标签并区分评分与池变化；未知高分链仍参与选择但不变成已知假阳性；缺失几何或同几何标签冲突拒绝统计。随后同一private／0GPU／2CPU／8GiB／5分钟generanno作业只读真实已核验记录，唯一outputs/M28-A31-SCORE-POOL-TRAIN-DIAGNOSTIC，命令sbatch --parsable sbatch/M28-A31-SCORE-POOL-TRAIN-DIAGNOSTIC.sbatch。主对话直接消费短作业，不加入现有七根监控或新建对话。失败只允许一次最小工程修复，保留原失败证据，不救指标。

## 原自然 TRAIN 池的标签辅助选择上限

2026-10-07。在等待A40闭合期间，单独量化已关闭A31_A3自然池的评分空间；原DEV A6／owner池上限不是该1491窗人口的上限。沿用同一readout13398521／CPU13398524、16321自然链、原2380／12566／1375标签及2916／1086参考分母；只读已闭合的natural_scores.jsonl、windows.jsonl和原draw身份，不读取A40输出、模型／checkpoint／features、DNA、catalog、DEV、Setaria或A33。

每窗仍使用原同链向span非重叠选择器及完整候选数组，已知非匹配和protected−1仅在本标签辅助诊断中赋0。第一个目标给known1赋N+1，再给至少5-exon的known1加1，N为该窗总候选数；这保证先最大化正确链数量，再在同数量解中最大化长链数。第二个目标只给至少5-exon的known1单位分，用于其单独计数上限。两个上限各自足够不证明能联合满足；首个目标构造的实际解满足原五条件时，只能说明任意标签辅助评分下有可行见证，不是实际A40／模型PASS。任何单项计数上限不够则表明此固定池下原五条件不可能同时成立；上限足够但未构造见证时保留joint feasibility未证实，不作失败或放行判断。

两个新增人工检查覆盖小池完整子集枚举对最大正确／最大长链的核对、首尾相接、unknown与空池，以及单项上限充分但联合见证未成立的解释；缺失任一原CPU证明必须在标签访问前拒绝。它们与实际统计放在同一次generanno／private0GPU／2CPU／8GiB／5分钟／no-requeue作业，依据原三格诊断实际47秒留余量，不新增独立人工检查作业。独立outputs/M28-A31-NATURAL-TRAIN-SELECTOR-BOUND仅保存summary、1491窗统计及选中索引、STATUS和JOBID，预计低于3MiB；主直接消费短作业，不改变A40链或固定heartbeat对象。诊断工程失败保留并停止，不自动重提或修改科学条件。

这是in-sample窗口／自然链配对的探索性几何上限，不是可学习分数的证明、独立gene F1或泛化。它不改变已关闭实验、A40训练配方或停止条件；不将eligible参考数减正确候选数解释成去重gene漏检，也不将oracle的known_nonmatch0描述成生物学precision1。

实施回执：CPU13419849实际root／batch／extern均COMPLETED0:0，27秒；两项人工检查PASS0.004秒，STATUS／JOBID及全部1491窗、16321链和原计数／分母相符。原自然池全部2380正确配对可共存，长链上限788，标签辅助解满足原五条件；[选择空间与解释边界](a31_score_pool_train_diagnostic_result.md)另记，不将此见证作为A40或论文PASS。
