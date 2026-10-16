# A31 固定局部状态路径 caller

2026-10-06，首次运行前固定。检验现有局部15类预测能否直接支持完整CDS链选择，不再拟合质量头或调整其零点。A6的后池遗漏主要是非正gain，A26虽增加正确链也增加非匹配，本方法不改变这些已关闭实验或其门槛。A23生成覆盖通过后，同一新caller分别作用于冻结A3和A23，区分caller收益与生成器差异。

## 固定分数和选择

每条原自然链按本项目ascending相位渲染CDS、intron、donor、acceptor、start、stop角色，累计coding phase跨外显子保留。沿用candidate_roles及原fine label优先级，start／stop连续3base覆盖端点。分数为候选span内逐base求和：logP(该角色)−logP(intergenic0)。等价于原logit差之和，所有角色权重1、null0，无偏置、温度、长度归一化、长度先验或阈值扫描。模型FP32不变，只用FP64累积分数。

这是固定局部发射势，不是校准的基因存在概率、完整HMM概率或新颖性主张。它不能补回生成器未产生的完整链。候选池仍是已关闭的greedy＋beam4／completed4／cap128几何联合池；不改搜索、ORF、5bp内含子、边界、owner或全局装配。窗内读出继续原零阈值span非重叠最优选择；所有未知候选仍参与选择，不按R4标签筛池。

## 冻结模型与三个主要读出

不更新任何参数，不写模型checkpoint，不继承或继续训练评分头。

- A3 原始greedy／query-logit读出为主基准：冻结A3末004608，仅保留terminated canonical chain，重复几何取最大query logit，沿用原零点与非重叠选择。
- A31 A3：同一A3末模型的局部logits，重用原16321自然联合候选。
- A31 A23：冻结A23末004608及其局部logits，重用已核验18950自然联合候选。

逐variant独立判定，不按新结果改分数配方或挑另一个checkpoint。A3两读出同时改变完成池与分数，不能单独归因于其中一个组件。A31 A3与A31 A23的caller配方相同；两模型参数相差160且共同参数训练结果不同，不称容量严格匹配或多seed因果证据。A6末质量头在同一原自然池的既有选出计数仅作辅助对照，不替代较强主基准，不重跑其训练。

## TRAIN 范围与首次门槛

原1536draw对应的1491唯一TRAIN窗、六条原TRAIN染色体、GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345及R4全部eligible reference／unknown保护不变。cached特征仅用于本次TRAIN准备；无新增GLM编码、DEV、sealed或Setaria。自然分数文件全部关闭后，再打开TRAIN catalog重算候选标签、选择及基准；不注入缺失正确链。原始greedy重放仅核对现存来源并取query logit，不产生新候选池。

此阶段不是完整染色体GFF评价。按候选／窗口配对报告selected correct／known_nonmatch／unknown及所有eligible参考；未知标签不作为假阳性计入known-label proxy，但保持原选择和独立计数。proxy P=correct/(correct+known_nonmatch)，R=correct/eligible references，F1=2correct/(eligible references+correct+known_nonmatch)。这些不能称独立gene precision／recall／F1。

每个A31 variant须同时满足五项预声明条件：拟南芥和水稻known-label proxy F1各自严格高于原A3 greedy；两物种proxy precision各自不低于同一主基准；合并至少5exon selected correct严格增加。用整数交叉乘法判断，不按浮点四舍五入判gate。采用P/R平衡准备条件，不修改A26或A30的既有known_nonmatch数量门槛。全部原计数、未知选择和长链分母完整报告。

任一variant通过，只支持为该明确模型另行固定完整raw-FASTA caller／DEV协议；不自动打开DEV或发布。失败的variant关闭，不救阈值、改角色权重或添加质量头曝光。最终仍须完整参考／背景span、原生ANNEVO／Tiberius／Helixer同输入比较、独立迁移、同源隔离和相称消融。

## 验证与资源

先四项数值测试，再按原metadata固定八项真实TRAIN预检：原R4渲染及正反方向相位对应；logit差与直接log-probability差恒等；逐base平移、padding和span外不变；零分选择null、非重叠目标及空池／非有限拒绝。真实案例核对两个冻结末模型／自然池身份，独立log-probability求和与分数一致、原greedy来源集合相同、所有参数未变，无梯度、optimizer0。八窗不判效果门槛，不按分数换案例。

preflight为private／1RTX3090 24GiB／2CPU／16GiB／10分钟，generanno，唯一M28-A31-LOCAL-PATH-CALLER-PREFLIGHT。通过后全1491窗双模型读出与CPU复算使用独立目录，实际时限依据真实预检与同范围旧读出来确定。若长等待，原固定watcher／heartbeat只扩展当前对象，不新建对话。A30原运行和冻结五gate保持不变。
