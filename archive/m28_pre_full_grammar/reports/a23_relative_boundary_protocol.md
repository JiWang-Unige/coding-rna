# A23 相对边界距离生成器协议

2026-10-06，首次预检及拟合前固定。A22显示本批次536个原已知正链／窗口对未进入自然池，单改评分无法找回这些完整链。本实验检验一个独立生成结构：在原A3边界指针中增加相对单调下界的距离项，能否改善自然完整链覆盖。它不修复首起点排名，也不把距离提示本身称为论文创新；A17局部初始context、A21链评分和本实验是不同机制。

## 唯一结构改动

保留A3全局初始context、query通信、状态／GRU／EOS、起点指针、DNA／ORF语法和greedy＋beam4代码。end与next_start的原geometric logits分别加一个相对距离项。有效长L、原调用下界lower、坐标x，d=max(0,x−lower)，输入为d/L和log1p(d)/log1p(L)，共享Linear2→32、GELU、无bias的Linear32→2；两列分别对应end和next_start。原mask和几何normalizer不变，start不加项，padding不进入L。

这是相对“现有函数下界”的距离，不冒充精确内含子长度：teacher next_start用previous_end+1，raw grammar用previous_end+5，原4bp支持差异保留。不给它添加新入口mask或事后重归一化。新160个参数在全部A3原参数之后初始化；末层权重固定初值0，原共有权重仍fresh seed0逐tensor相同，初始模型分数／监督／匹配不变，不继承训练checkpoint。总参数649617，比A3的649457多160，不是容量严格匹配消融。

## 固定数据目标和继续条件

原六条TRAIN染色体、GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、四块D2048特征、定向24576bp窗口、R4正／alternative／partial／unknown政策不变。沿用1536原draw的顺序、importance权重和三遍曝光，共4608更新；H128/Q16、FP32/batch1、四loss系数1、known-query mean BCE、原Hungarian精确成本、AdamW lr3e-4/wd.01/betas(.9,.999)/eps1e-8/clip1。只有末004608为主读出，不加epoch、seed、beam或类别权重，不改unknown语义。

预检通过后一次固定fresh拟合，每256步保存工程checkpoint。原1491窗先完整关闭自由greedy＋beam4几何联合池，再打开原TRAIN参考标注；width4/completed4/max_exons128不变，无参考注入。配对原A3自然池重算标签；逐物种完整报告correct／known_nonmatch／unknown／全部候选及至少5exon correct。继续条件与A17的生成机制筛选相同：两个物种correct候选／窗口对各自严格增加，pooled至少5exon correct严格增加。FAIL关闭本候选，不救参；PASS只支持另固定完整caller／整链评分协议，不自动执行DEV或发布。候选数不是选中基因FP或独立基因F1。

## 必要验证和执行

四项合成检查覆盖共有参数与零delta分数；真正随lower改变的相对项且start不变；padding、teacher／free一致、末层非零梯度及零末层下输入层初始零梯度；原DNA-only greedy／beam确实使用新指针。随后沿用原metadata选择的八个真实TRAIN case，optimizer0，对照fresh A3的实际四loss与离散分配，核对共有分支有限梯度、新输出事件按实际正链／多外显子支持有无梯度、两物种覆盖两事件及全部参数不变。零末层导致第一层初始零梯度是预期计算，不错误要求其首步非零；训练后其梯度仍待实际更新核实。

预检generanno、1RTX3090至少24GiB／2CPU／16GiB／10分钟。fit上限3小时、自然读出2小时、只读CPU核验30分钟，均独立M28-A23-RELATIVE-BOUNDARY阶段目录；正式提交依实际预检计时与实时队列／维护决定合理运行余量，不改变4608终点。A17实际拟合1时22分、生成池54分38秒提供参考；预检含两模型前向，不当成单步训练计时。完整fit／readout／CPU runner已实施：每256步保存工程checkpoint，记录新增两层的实际梯度，要求整次曝光内两层均接到非零梯度；CPU直接复算末层／输入层权重相对fresh的变化、全部原draw／weight／监督身份与新旧自然池DNA标签及原三项gate。没有额外更新或训练干预。预检不等于拟合或候选质量结果。没有新GLM前向、DEV、sealed、Setaria、作业取消或Git／公开发布。等待复用固定监控对话和同一heartbeat，不新建监控。

## 已有方法和论文边界

[Tiberius原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/)已经把深度网络与可微HMM用于整基因结构预测；结构约束不是新概念。当前[作者仓库](https://github.com/Gaius-Augustus/Tiberius)已有2026年的非哺乳类模型，相关[多clade论文](https://www.biorxiv.org/content/10.64898/2026.04.24.720536v1)为预印本。[ANNEVO作者仓库](https://github.com/xjtu-omics/ANNEVO)当前标为v2.3.3，已提供训练代码。它们的作者指标不能替代本项目同输入、同完整参考的原生比较；最终论文须明确实际版本、模型权重和适用clade，不能把旧冻结结果自动称为最新方法结论。既有结果与模型不自动升级或重跑。相对项是否有效及是否能支持更完整的新方法，目前均未证实。
