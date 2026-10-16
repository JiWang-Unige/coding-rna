# A25 首起点提示的局部注意力生成器

2026-10-06，首次预检和拟合前固定。A24中缺相同首起点为536个遗漏配对的最大单类246个，但不是多数；至少5外显子遗漏仍以后续连接为主。本试验只检验显式首起点表征能否改善自然完整链覆盖，不把A24形态当作模型失败的因果证明。A17和A23原协议、权重、终点及gate不变。

## 结构和监督

基线是A17的query依赖局部注意力，不是均值池化。保留有效长度的Q16等长连续context bins（短contig至少1base）、query通信、全局start/end/next_start指针、GRU/EOS和原DNA/ORF语法。新增Linear H128→1、无bias、权重初值0，对共享hidden的每个有效base产生initiation logit，作为公共加项加入每个query的原局部context注意力logits。它只重加权context，不裁剪链指针、不硬选首起点、不加NMS/topK/参考锚定。推理只读有效hidden，padding不参加context。teacher和原greedy/beam使用同一新window接口。

新增128参数，总649585。全部原参数在新层前按原构造顺序fresh seed0初始化；共有参数与fresh A17逐tensor相同。零initiation权重令初始attention、状态、四项旧loss和匹配与A17相同；不是与全局context的A3初始输出相同。不继承任何训练checkpoint。

新增条件首起点CE：目标为原R4 eligible positive_chains的去重第一CDS base。归一化支持为原fine标签已知且不是13的有效base，另纳入eligible首起点；非eligible class13和start triplet的其他base不作负例。fine IGNORE、noncallable、padding和已mask的alternative/partial位置不进入支持，除非原eligible正链明确给出这个首起点。没有eligible正链时该项为带图0，不训练背景absence。不把这一条件排序解释为生物学无基因概率。unknown位置仍可参加无参考推理attention；loss mask只保证新增项不向未纳入支持的独立logit直接施加负例梯度，不保证共享特征没有其他梯度。

新CE系数固定1，和原四项系数1一起乘原importance weight；原known-query mean BCE、Hungarian成本及unmatched unknown政策不变，新增CE不进入离散匹配。本候选同时新增表征和辅助目标，不是二者分离的因果消融，也不是严格容量匹配消融。

## 固定执行和继续条件

原六条TRAIN染色体、feature revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、D2048×4特征、24576bp定向窗、1536原draw顺序和importance权重不变；三遍共4608更新。H128/Q16、FP32/batch1、AdamW lr3e-4/wd.01/betas(.9,.999)/eps1e-8/clip1不变。每256步工程checkpoint，只有末004608主读出，不额外epoch/seed/beam，不调整unknown政策。

预检通过后一次fresh fit。1491原TRAIN窗先关闭原greedy+beam4自由池，再打开参考建标签，max_exons128/beam_width4/completed4及去重不变。配对原A3自然池复算correct、known_nonmatch、unknown、全部候选和至少5exon correct，逐物种报告。冻结gate与A17/A23相同：两物种correct候选／窗配对各自严格增加，pooled至少5exon correct严格增加。FAIL关闭候选不救参；PASS仅允许准备另行冻结的完整caller/整链评分实验，不自动放行DEV或发布。A17对照为辅助机制比较，不改变该gate。候选覆盖不是选中基因F1。

五项合成测试覆盖零初始化A17一致、新项实际改变attention但不截断指针、辅助支持及unknown直接梯度、teacher/free/padding/旧loss梯度、原raw greedy/beam。随后八个原metadata固定真实TRAIN case做optimizer0预检，核对共有参数、初始四loss与匹配、两物种正链辅助项有限且新层非零梯度、无正链辅助项0、所有参数不变。预检1RTX3090至少24GiB/2CPU/16GiB/10分钟、generanno；正式fit上限3小时、自然读出2小时、CPU核验30分钟，依据实际预检和队列执行。没有新GLM前向、DEV/sealed/Setaria、作业取消、Git或公开发布；长等待复用原固定watcher和heartbeat。

## 先例与论文边界

[Deformable DETR原论文](https://arxiv.org/html/2010.04159)的两阶段变体已从encoder预测选择高分proposal作为decoder query；[DAB-DETR原论文](https://arxiv.org/html/2201.12329)已用显式坐标和尺度调制attention。这些图像研究不证明基因注释收益，但足以否定“数据依赖或位置提示query本身全新”的说法。A25当前只是基因结构生成假设；论文创新和性能优势须由完整raw-FASTA caller、真实同输入原生比较、独立测试/适用范围及相称消融支持。
