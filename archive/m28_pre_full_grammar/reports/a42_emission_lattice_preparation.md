# 角色势语法解码原型

2026-10-07北京时间。A42已完成不依赖query或beam的oriented窗口解码原型，以及另定义上下文的整条contig联合解码接口。最初十二项窗口检查及随后八项联合上下文检查均已通过；这两阶段只用人工DNA和人工势。真实TRAIN八窗工程预检及完整1491窗机制比较另见a42_emission_lattice_train_protocol.md，不能把人工接线通过当科学PASS。以下按阶段保留语法、资源、owner反例及联合接线的实际证据。

## 研究依据与边界

原TRAIN的A27有536条missing参考链，其中315条被所有query的局部选择排除，221条至少有一个query通过局部选择但仍未返回。A28在后一组1708个query路径上记录首次global entry、active extension、completion和bound剪枝分别753、468、279、208次，不能把query事件当独立基因。

A34参考辅助诊断把这536条missing链分为190条非正、179条正分未选及167条正分选中且不挤掉原正确链。167不是自然生成的改善量；只去除搜索剪枝也不能消除其余评分及竞争问题。因此本原型只回答既定加性角色势能否不经过query几何剪枝而完成语法选择，不宣称解决模型评分、SOTA或泛化。

精确语法动态规划本身不是论文创新。角色势仍是外部输入，未增加或训练表示；这里的最优性仅针对单个oriented窗口、下述受限语法和目标，不是全基因组最优、生物学正确或校准概率。

## 语法与目标

保持原contiguous ATG和末STOP、GT或GC至AG、至少5bp内含子、无模糊coding base、无internal STOP及显式max_exons；人工检查含原128上限以内的较小cap，不自动制定真实数据实验。跨外显子保留实际未完成密码子，而非仅相位。内含子允许N；split start或stop、非标准遗传码不在支持范围。

目标是十五类局部logit减intergenic的FP64角色累计，不加query、长度奖励、角色系数、温度或阈值。donor覆盖前外显子末base、acceptor覆盖后外显子首base，随后start和stop三base覆盖，匹配原candidate_roles优先级。一个窗口内按整gene span选择非重叠完整CDS，允许相触；intergenic为0，全零返回null。

对同一donor位置、外显子数和实际codon carry，只保留最高总势历史。donor满足最小内含子长度后，减去该相位的intron前缀势，再合并相同外显子数和carry的成熟历史；任一未来acceptor补同一个前缀势，因而其后续可行边与目标一致。末STOP更新gene外部的最优前缀，新的ATG从此前缀进入。没有site top-k、pointer、query slot、beam或completed数限制，非正partial path仍可延伸。

所有合法edge均进入计算，cap外的未终止路径不被修成STOP。累计采用FP64前缀和及差值，仍有浮点舍入，不能称为逐bit精确实数最优；非零平分保留确定性首次到达，未承诺与原有限池DP的tie身份相同。

## 人工检查与执行条件

src/m28/emission_lattice.py和tests/test_m28_emission_lattice.py不加载任何模型。独立oracle用DNA拼接穷举小序列的所有合法链，不复用ORFIndex，再在完整链集合上用独立interval recurrence求最优gene set。十二项检查覆盖逐角色算术、竞争及相触、多于四个完成基因、零或负势、负partial而正完成、1bp中间和3bp末外显子、同相位不同实际carry、GT/GC和N、cap、softmask和公共logit偏移、25个人工随机小序列以及实际可达的对齐输入错误。人工随机seed固定，只用于代码核验，不是模型seed选择。

唯一命令sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-CHECKS.sbatch。沿用实际成功的generanno环境及Mamba初始化；private-teodoro-gpu、0GPU、2CPU、8GiB、2分钟、no-requeue，唯一根outputs/M28-A42-EMISSION-LATTICE-CHECKS及%x_%j日志。产物仅STATUS、JOBID和小型日志，不产生权重或缓存。提交前核对五项活调度、相关硬约束、源文本和输出冲突；这里没有材料级新增存储。

检查失败先保留并报告；若仅有明确人工fixture或实现错误，可在目标、语法、资源及无研究数据条件不变下修正一次，使用独立R1根，不覆盖原根。再次失败停此准备版本，不修改科学阈值或任何实际实验。人工检查通过也不自动接入A41、读取DEV或独立物种；真实运行还需另固定语料、资源、所有新链的标签保护、owner及整尺评价协议。

## 实际执行

仅一次提交13421286，2026-10-06 23:29:57 UTC提交，23:29:58至23:30:04 UTC运行，6秒；root、batch、extern均COMPLETED0:0，STATUS为COMPLETED且JOBID匹配。实际private-teodoro-gpu、gpu034、0GPU、2CPU、8GiB、2分钟上限，Requeue0、Restarts0，batch MaxRSS468244K。canonical Command、SubmitLine、WorkDir及唯一日志相符，十二项PASS共0.115秒；没有失败或修正重跑。

提交前四个新增文件与canonical逐字相同，两Python AST及SBATCH语法通过，五项活调度及节点、partition、reservation元数据没有相关硬约束，输出根和同名日志无冲突。共享节点的空余不是独占资源保证。一次SSH banner超时后，一次连接恢复；不据此判断任何科研作业失败。Slurm默认CEST加0200；上述提交及起止时间用TZ为UTC的sacct核对，不把默认显示的01:29:58误标UTC。

所有检查仅用人工序列、人工局部势和独立小序列oracle。通过支持在既定受限语法及目标下继续准备，不证明真实注释性能、源模型公平性或论文创新。十二项人工时序不能作为速度benchmark；全长度窗口、128外显子cap的耗时及内存仍未测量。下一步先准备人工长窗资源探测，决定真实TRAIN实验是否可执行；任何真实新搜索、标签配对或完整raw评价须另固定科学范围，不能把本次人工通过自动当作A41或DEV放行。

## 人工长窗资源探测

2026-10-07北京时间，在任何真实新搜索前固定两个24576bp人工窗口，A42解码源码不变。case0在C背景512bp后放一条128-exon链：首ATGAAA，126个CCC中间外显子，末CCCTAA，127个GTTAG内含子。目标span每base的最大角色边际为2，其他角色为负4、intergenic为0；独立pointwise上界为2050，cap128必须返回这条链并达到上界。随后同一人工输入用cap127，必须排除该链且严格低于上界；不把这一cap检查当新科学参数。

case1用torch CPU generator seed0产生24576bp均匀ACGT，十五列潜势全0，cap128必须返回null。它用于观察密集随机motif上的状态和边数，不代表真实基因组、物种分布或最坏复杂度。两case都核对所选链的实际DNA拼接ORF、内含子motif及长度、外显子cap、span非重叠、直接FP64角色势和总势；case0的已知最优来自pointwise上界而非全长穷举。计时只包含decode，另记本CPU进程的peak RSS。

命令sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-LONG-CHECKS.sbatch，固定array0-1且最多并行2，private-teodoro-gpu，每task为0GPU、2CPU、8GiB、3分钟、generanno、no-requeue。根outputs/M28-A42-EMISSION-LATTICE-LONG-CHECKS，各case独立STATUS、array task JOBID、started.json、records.jsonl及完成时summary.json；日志用%x_%A_%a。仅KB级文字产物，无权重、真实DNA、潜势tensor缓存或研究标签访问。提交前读实际调度、源文本和冲突，不修改A41或其唯一监控。

原十二项检查已闭合，不重复运行。这里的具体未决问题是全窗口状态扩展成本及真实128-cap接线：通过才继续准备固定TRAIN协议；超时或OOM先按保存的开始阶段及Slurm资源定位实现瓶颈，不延长、重跑、改语法或宣布机制失败。纯明确fixture或实现错误可在科学范围不变时独立R1修正一次，保留原产物；任何真实TRAIN、DEV或A33仍不执行。

array13421438于2026-10-06 23:41:53 UTC仅提交一次，两task于23:41:54 UTC在gpu034同时开始。13421438_0的root、batch、extern均COMPLETED0:0，23:42:01 UTC结束、7秒，batch MaxRSS373336K；13421438_1同样全部COMPLETED0:0，23:44:47 UTC结束、2:53，batch MaxRSS401352K。实际每task为0GPU、2CPU、8GiB、3分钟，Command、SubmitLine、WorkDir、日志、array0/1、no-requeue和Restarts0相符，各case的STATUS、JOBID和完成summary齐全。没有失败、超时、重跑或资源延长；原源文本及三个新文件提交前逐字同步，语法通过，实际调度、预约、资源和输出无硬冲突。原manifest的WIDTH24576及原A31的max_exons128已从代码核对。

case0的cap128实际返回单条128-exon链，gain2050达到上界，decode0.864621秒、进程peak RSS438664KiB；cap127返回单条127-exon链、gain2032，decode0.862805秒、peak443504KiB。case1实际417个ATG、3098个GT/GC donor、1494个AG acceptor和1121个STOP，返回null、gain0，decode167.549558秒、peak469140KiB。两根仅小型JSON及marker：case0 summary17253字节、records4340字节；case1 summary1443字节、records516字节。Slurm采样MaxRSS与进程峰值分开报告。

三个decode的所选DNA拼接ORF、角色算术、span和cap检查均通过。已知链检查没有全长穷举；随机零势检查也不是物种性能或最坏复杂度证明。该实测说明当前未复用实现能完成这两个人工全窗，但密集随机窗明显慢，不能用稀疏已知链约0.86秒宣称全基因组速度。

## 同目标计算复用

case1实际legal_exon_edges276003，而exon_extensions32010632、donor_states_created28297444、donor_states_dominated27285470、completed_extensions3506456，maximum_active_states2667。源码在每次history扩展重复同一局部边的角色增量，并在比较是否被支配前构造链。改为在既有entry、actual carry、first标志的合法边缓存中只计算一次增量；只在原严格score比较胜出时构造donor历史或完成链。语法、所有边、外显子上限、gain、FP64累计顺序和first-arrival平分规则不变，没有零势快捷返回、score beam或阈值。

先保留原13421438两根和三个原人工答案，不重写旧检查结论。新版本核验命令为sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-REUSE-CHECKS.sbatch：private、0GPU、2CPU、8GiB、3分钟、generanno、no-requeue，独立根outputs/M28-A42-EMISSION-LATTICE-REUSE-CHECKS。新改动使原十二项独立小序列检查需要再运行一次；另用十三项中的一项，对全部三个已保存人工全窗的链身份、各链gain、总gain及每个遍历计数逐值对照，并重核DNA和角色算术。只读这两个人工结果根，不读研究数据，不据时间差宣称公平速度对照。

新核验13421592于2026-10-06 23:53:46 UTC仅提交一次，23:53:47至23:55:54 UTC运行，2:07；root、batch、extern均COMPLETED0:0，STATUS与JOBID相符。实际private、gpu034、0GPU、2CPU、8GiB、3分钟，batch MaxRSS498296K；Command、SubmitLine、WorkDir、唯一日志、Requeue0及Restarts0相符。提交前四个相关文件与canonical逐字一致，两Python AST及SBATCH语法通过，实际五项调度检查及独立输出无硬冲突。无失败、修正重跑或资源延长。

原十二项PASS0.105秒，完整三个已保存人工输入的逐值对照检查PASS116.928秒，共十三项。新cap128和cap127的decode分别0.607728与0.607723秒，进程peak RSS分别443460与448568KiB；随机窗115.377607秒、peak491560KiB。三个输入的所选链、各gain、总gain以及所有遍历计数与原13421438逐值相等，包括随机窗32010632条扩展及2667个最大active状态；没有通过减少边或零势快捷分支制造更快的结果。随机窗本次测量较原167.549558秒减少约52.17秒，单次不同时间运行及热状态、共用节点负载未控制，不用作公平速度benchmark或因果提速比例。

本次代码与人工资源核验支持继续准备固定TRAIN窗口实验，不能替代真实DNAs的耗时、数值闭合、标签保护或性能结果。A41的模型、固定拟合、自然读出、CPU核验、五条件及唯一监控均未改；未读其loss、gain或中间性能，未执行A42真实搜索、新GEN、R4、DEV、独立物种、Setaria、A33或对外发布。

## 全基因组接线边界

现有assembly先按chain midpoint筛nearest owner，再在拥有的候选上作span DP。A42目前输出窗口内最优gene set，不能直接把它当原候选池接入：例如窗口[0,100]与[50,150]中，A链[50,80]的gain5且owner0，B链[60,99]的gain10且owner1。全窗先选会保留重叠的B、丢A，之后owner过滤又丢B；原顺序则丢B后保留A。两链可由C背景上的50及60位置ATG、77及96位置TAA形成完整同向合法ORF，位置55的coding状态3给5、位置88的coding状态2给10即可形成相应加性gain。

上述100bp小例最初仅从源码推导；下节实际运行了原尺寸的正负链及跨窗口反例。它们不是已经存在的A42 raw适配器故障。正式raw接线不能只在已选gene set上补过滤；仅将owner提前、然后输出每个窗口的最优gene set，也不足以支持跨窗口竞争。需要保留必要备选，或另固定联合选择的分数及上下文定义。owner反例闭合时尚无A42真实TRAIN或raw接口执行，不以人工资源通过放行DEV、独立物种或论文claim。

## Owner接线人工反例运行协议

使用原WIDTH24576、STRIDE12288构造两个窗口[0,24576)及[12288,36864)，固定max_exons128，只有C背景人工DNA和人工15-role局部势。新检查调用实际A42解码、现有nearest_owner、assembly及span DP，不读取模型、真实DNA、缓存或参考。检查的失败意味着反例fixture或对源码顺序的理解不成立，先保留实际错误，不能直接推断真实性能。

正链case：A为[15000,20001)、gain5、owner0；B为[15010,23002)、gain10、owner1。窗口0的A42先选B，窗口1同一B的合法几何在全零局部势下分数0并返回null；selected-only assembly应为空，而保留两候选并先owner筛选应保留A5。负链case：窗口0的oriented A为[10000,16000)、B为[0,11001)，gain分别5与10；全局负链A为[8576,14576)、owner0，B为[13575,24576)、owner1。窗口1用实际反向互补及正确偏移核验B的完整ORF和零势，检查同样的损失。两case均独立拼接DNA检查ATG、末STOP、无internal STOP及直接FP64角色势。

第三case固定全局A[8000,14000)、B[12000,21000)、C[18001,24001)，gain5、8、10，A/B为owner0、C为owner1。A与C相容，B分别与A/C冲突。完整owner候选池的全局span DP应返回A+C、gain15；模拟“owner过滤后每窗先作局部DP”只返回B和C，再作全局DP应只保留C、gain10。这是对拟议快捷修法的反例，不冒充已经实现的A42 owner-aware接口。真实A42在窗口0的完整人工势应返回A+C，窗口1返回C；不根据这个case改真实潜势或阈值。

唯一命令sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-OWNER-CHECKS.sbatch：private-teodoro-gpu、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue。独立根outputs/M28-A42-EMISSION-LATTICE-OWNER-CHECKS，仅STATUS、JOBID和%x_%j小日志，不重复原13项已闭合检查，不改变A42解码或原assembly源码。提交前核对五项活调度、资源约束和精确目标冲突。失败保留；只有明确人工fixture或实现错误可独立R1修正一次，资源和无研究数据范围不变；再次失败或超时停本准备版本，不延时或救援阈值。通过也只确定接线要求，不能当科研PASS。

## Owner反例实际结果

13421848于2026-10-07 00:10:41 UTC仅提交一次，00:10:42至00:10:49 UTC运行，7秒。root、batch、extern均COMPLETED0:0；STATUS与JOBID匹配，实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS462256K。Command、SubmitLine、WorkDir及%x_%j日志相符，Requeue0、Restarts0；三相关文件提交前逐字同步canonical，Python AST及SBATCH语法通过，五项活调度、节点与预约无硬冲突或输出覆盖。没有失败、R1、资源延长或研究数据访问。

三项实际检查PASS，0.405秒。正负链均确认：完整owner候选选择返回A5，先A42局部选择再owner过滤返回空。第三项确认：完整owner池全局返回A+C15；模拟owner过滤后局部最优再拼装，只剩C10。所有参与链均通过独立DNA拼接ORF和直接角色势，负链窗口1通过实际反向互补及坐标换算；A42实际完整人工窗口结果与协议一致。

结论是当前selected-only输出不能直接充当原assembly候选池，owner-before-local-DP也不是充分修法。这是接口语义及拟议方案的反例，不是已实现raw接口的失败、真实数据性能或生物学结论。保留原assembly及A42窗口解码源码不变，后续可以固定TRAIN窗口机制比较，但完整raw接线须单独解决备选及跨窗口竞争；不得把窗口最优性或资源核验当作整染色体最优、DEV放行或论文结果。

## 真实TRAIN预检接续

随后已在任何A42真实搜索前另固定a42_emission_lattice_train_protocol.md，原1491窗支持、A3冻结势、单位角色、gain0及max128不变。首次只取每物种每方向原序前两窗，共八窗，不取missing或按新势挑窗。GPU发射预检13421953及CPU无标签解码array13421967已全部闭合：八窗原86候选分数最大误差0、选择索引相同，冻结权重不变、optimizer0；八真实decode为97.31至127.36秒、进程peak RSS不高于720072KiB，完整ORF及直接角色势检查通过。没有标签连接或性能结果。详细资源及实际结果见该协议；上述无真实DNA的边界限定人工准备阶段，不把新增真实TRAIN工程预检仍称为source-only，也不把它当完整1491窗比较或raw性能。

## 整条contig联合上下文人工准备

原owner反例不能靠先输出每窗最优gene set解决。新的准备版本在解码前为每个基因组碱基p固定上下文：只在包含该base的窗口中，选择中心最接近p加二分之一的窗口；距离相同取较小genomic start，正负链使用同一genomic规则。原WIDTH24576、STRIDE12288的grid仍由manifest.starts给出，最后一窗按原定义平移至contig末端；长度不超过WIDTH时只有一窗。这里只用人工窗口，helper接收已给定完整grid，不改manifest或原assembly。

该base的十五列FP32局部势完整复制自选中窗口，不平均、不加权、不重排角色列，也不纳入padding。负链局部势已按窗口反向互补排列，复制到整条反向互补contig的正确位置；全局链坐标转换保持与各gain配对，负链结果顺序为转录方向而非genomic升序。所有窗口按genomic顺序各出现一次，DNA片段、有效长度、FP32及有限性相符才进入联合解码；缺窗、错序、错DNA或padding对应的检查失败时停止，不猜测补齐。

复制后在整条oriented contig上调用未改变的emission_lattice语法DP，不经过局部gene set选择或chain owner过滤。单位角色、FP64累计、gain0、null0、actual codon carry、最少5bp内含子、GT或GC至AG、完整ATG至STOP及max_exons128均不变。链可跨上下文接缝，span可大于原窗口。因此上下文、链分数及允许span均是新定义，不是原A42窗口目标的等价算法，不以人工恢复A5或A+C15宣称在旧分数下实现同样改善。

人工检查分八项：逐base独立枚举owner，含平分、负链十五列、非整stride末窗和单一短contig；完整窗口流及DNA/长度错误；原正负链selected-only owner反例；跨窗口竞争A+C；正负链各一条30000bp单exon完整ORF跨接缝且无法在任一24576bp窗完成；带N内含子的actual carry跨多个接缝及显式cap；24个小序列上下文随机势用原独立DNA穷举oracle核对最优目标，并检查全零返回null。人工seed固定23，仅生成核验势，不作模型选择。所有输入均为人工DNA和人工势，不读研究DNA、缓存、权重、参考、GEN、标签、DEV、独立物种、Setaria或A33，不使用A41修正。

唯一命令sbatch --parsable sbatch/M28-A42-CONTIG-LATTICE-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue，独立根outputs/M28-A42-CONTIG-LATTICE-CHECKS及%x_%j日志，仅STATUS、JOBID和短检查日志。提交前核对五项活调度、实际节点/预约、源文本及明确输出冲突。纯明确fixture或实现错误可在同一目标、语法、资源及人工数据范围内独立R1修正一次，不覆盖原根；再次失败、超时或OOM停本准备版本，不延时、改cap或跳过条件。

人工通过只支持该新上下文接口及受限语法的接线继续准备；不证明染色体级运行成本、真实评分、注释性能、机制创新或泛化。完整1491窗A42及原五条件、独立CPU核验和唯一固定监控保持不变。实际raw仍须在相应科研门槛闭合后另固定科学/资源范围，不能凭本次人工检查放行。

## 联合上下文人工实际结果

13422687于2026-10-07 01:15:26 UTC仅提交一次，01:15:27至01:15:33 UTC运行，6秒；root、batch、extern均COMPLETED0:0，STATUS及JOBID匹配。实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS495040K；canonical Command、SubmitLine、WorkDir、日志、Requeue0及Restarts0相符。四个新增或相关文件提交前逐字同步canonical，两Python AST及SBATCH语法通过，五项活调度、节点/预约和输出根/同名日志无硬冲突。没有失败、R1、资源延长或研究数据访问。

八项检查PASS，共1.615秒。正负链的独立逐base owner枚举与复制结果逐值相等，包含十五列、不变原输入、平分、40000bp非整stride末窗及7bp短contig。错序、缺窗、错DNA、padding、非FP32及非有限势按协议拒绝。原正负链owner反例在旧selected-only assembly中为空，新base上下文目标均返回A5；原非owner竞争链在新上下文中的势为0，明确不是旧chain-owner分数下的等价恢复。跨窗口竞争返回A+C、gain15。

正负链均返回完整30000bp单exon链、gain300000，任一24576bp局部窗均无法独立完成它；这是已知人工链的接缝和坐标结果，不是发现真实长基因。GTNAG内含子及actual codon carry跨多接缝时，cap2返回指定完整链，cap1返回null；新上下文的24个人工随机小序列最优目标与独立DNA穷举一致，全零均返回null。八项总时长和稀疏长ORF不能作为密集染色体速度或内存benchmark。

因此联合上下文接口可继续做不读研究数据的规模准备，原A42窗口比较不变。当前实现仍在整条contig构造局部势、FP64前缀及合法边缓存；染色体级状态、边数和内存成本未闭合，不能因36864bp人工稀疏输入成功就接入真实全染色体。下一步先固定人工规模探测，按实际成本决定是否需要保持同一目标的计算改进，不改A42科学门槛、数据范围或原owner方案。

## 密集人工contig规模探测

固定49152bp均匀ACGT人工DNA，CPU generator seed0，原grid为[0,24576)、[12288,36864)、[24576,49152)，只探测正链；负链接线已有八项人工检查，不把两方向当重复性能样本。case0的三个有效窗15-role FP32势全0，已知全局目标为null0。case1在同一DNA上按genomic窗口顺序用独立CPU generator seed41生成十五列包含intergenic的整数势，取值负3至3；每窗独立，base上下文按新接口固定，不选有利样本或重设seed。

两个case均先复制原base-owned上下文，再调用未改变的DP，单位角色、FP64、gain0、null0及max_exons128相同；没有zero快捷返回、beam或site top-k。记录复制准备时间、完整decode时间、进程peak RSS、motif数、合法边、实际扩展及状态计数，区分全零历史与正分历史/回溯的成本。只对最终完整链独立拼接DNA检查ORF、motif、cap、span及直接角色算术，再核对总势与逐base最大角色边际的上界；case0另核null。case1的长序列没有全池独立穷举，不把合法性和上界检查称为独立最优性证明。

唯一命令sbatch --parsable sbatch/M28-A42-CONTIG-LATTICE-SCALE.sbatch，private、array0至1且最多并行2，每task为0GPU、2CPU、8GiB、12分钟、generanno及no-requeue；独立根outputs/M28-A42-CONTIG-LATTICE-SCALE/case_0或case_1，%x_%A_%a日志。12分钟是参照已有24576bp密集窗约115秒、长度翻倍可能增加候选边配对的运行预估加余量，不是预先认定二次复杂度或保证完成。仅小型JSON/marker和日志，不保存tensor、DNA、模型或缓存，整个任务不读取研究数据、参考/标签、DEV、独立物种、Setaria、A33及A41修正；不修改运行中的A42或固定监控。

提交前核对五项活调度、源文本、节点/预约及输出冲突。语法或角色算术错误保留并停止；仅明确fixture/实现错误可在相同科学/资源范围内独立R1修正一次，不覆盖原根。再次失败、TIMEOUT或OOM关闭此准备版本，不延长或改cap救援；根据保留的开始阶段和实际遍历量决定是否另准备同目标计算改进。两个case都完成，也仅支持这两个49152bp人工输入的可执行性，不据此放行真实染色体或外推每Gb速度。

## 必定非法终点的边界准备

原版legal_ends对每个entry、actual carry和first组合扫描整条contig的donor/STOP事件，即使早已有必定包含的internal STOP或模糊coding base。独立helper exon_end_bound.py以原ORFIndex的bad前缀和分frame STOP索引计算必要上界；不修改原ORFIndex、DP、评分、合法边或正在运行的A42/规模探测。下节将它接入独立解码版本，原版保持不变。

第一模糊base位于q时，合法外显子至多结束于q，仍可在它之前结束或转入内含子。actual carry与本外显子开头完成STOP时，结束必须早于该密码子完成，可保留1bp中间外显子及未完成STOP前缀；不按phase把不同carry合并。其余情况下第一同frame完整STOP位于k时，上界为k加3，保留它作为terminal和在STOP完成前的donor。上界以内仍逐事件调用原ORFIndex.extend，并保留原排序及同终点两种事件；必要上界不代替合法性判断，也没有长度、beam、site top-k或窗口宽度截断。

五项人工检查用独立carry加DNA拼接检查每个原事件，与未裁切ORFIndex及裁切后的合法事件、remainder和顺序逐值比较，覆盖25个seed29小序列、GT/GC、模糊coding/intron、split carry STOP、1bp中间外显子、terminal边界，以及30000bp完整exon不被窗口宽度限制。仅用人工DNA，不读研究数据、潜势、模型或参考。

唯一命令sbatch --parsable sbatch/M28-A42-EXON-END-BOUND-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue；独立根outputs/M28-A42-EXON-END-BOUND-CHECKS及%x_%j日志。提交前核对五项活调度和明确冲突。失败保留；仅明确实现/fixture错误可同范围独立R1修正一次，再次失败或超时停本准备版本，不覆盖、延时或省条件。通过只能支持此必要终点界，不能宣称整个DP提速、染色体复杂度已解决或真实科学PASS。

13422830于2026-10-07 01:27:14 UTC仅提交一次，01:27:15至01:27:20 UTC运行，5秒；root、batch、extern均COMPLETED0:0，STATUS/JOBID相符。实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS451444K，Command、SubmitLine、WorkDir、日志、Requeue0和Restarts0相符。五项PASS共0.081秒；独立DNA拼接、完整原事件及必要界保留的合法事件、remainder与顺序一致，30000bp exon终点33000未截到窗口宽度。四相关文件逐字同步、两Python AST及SBATCH语法通过，五项活调度和输出冲突核对通过；无失败、R1、延时或研究数据访问。

## 索引终点的独立解码接线

在终点界通过后，将原DP过程保留在独立emission_lattice_indexed.py中，复用原角色算术和状态类型；仅legal_ends按必要终点范围索引原事件，之后仍调用ORFIndex.extend。原emission_lattice.py、contig接口、assembly、完整1491窗A42及运行中的49152bp规模探测均不改。原合法边顺序、FP64操作、严格胜出比较、首次到达tie、所有遍历计数和cap保持，不增加新的评分定义；这里与base上下文变更不同，是拟议同目标计算改进。

三个检查分别核对七个人工小序列上的63组整数随机势和七次null，三种cap、actual carry、GT/GC及模糊base，用独立DNA穷举与原DP逐值对照；原24576bp人工128-exon上界及cap127；正负链30000bp完整ORF经过人工上下文后不按宽度截断。比较整个result字典，包括链、各gain、总gain及每项合法边/扩展/状态计数，不只比较最优标量。长序列仍只用已知人工答案及直接DNA/角色核对，不作全长穷举。

唯一命令sbatch --parsable sbatch/M28-A42-INDEXED-LATTICE-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue；独立根outputs/M28-A42-INDEXED-LATTICE-CHECKS及%x_%j日志。提交前核对五项活调度、源文本和精确冲突。错误保留；仅明确实现/fixture错误可同范围独立R1修正一次，再次失败或超时停本准备版本，不改变语法、评分、cap、门槛或正在执行的作业。通过后仍需规模输入的实际身份/分数/遍历闭合，不能由小序列等价宣称真实染色体可用或因果提速。

13422866于2026-10-07 01:31:39 UTC仅提交一次，同秒开始、01:31:47 UTC结束，8秒；root、batch、extern均COMPLETED0:0，STATUS/JOBID相符。实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS497868K，Command、SubmitLine、WorkDir、日志、Requeue0及Restarts0相符。三个检查PASS共2.768秒，完整result字典与原版逐值相等，包括每条链、每个gain及所有遍历计数；小序列另与独立DNA穷举相符，128-exon上界、cap127及两方向30000bp人工链均保持。五个相关文件逐字同步，两Python AST及SBATCH语法通过；代码定向差异仅为必要终点范围及说明，五项活调度和输出冲突核对通过。没有失败、R1、资源延长或研究数据访问。

## 密集规模探测实际结果

第一次提交连接在SSH banner阶段超时，未返回ID；恢复后从squeue、同名sacct、独立输出根和同名日志确认没有创建作业，再重复活调度核对后提交，未盲目重复。唯一实际array13422815于2026-10-07 01:25:39 UTC提交，两task于01:25:40 UTC在gpu034并行开始。task0于01:32:23 UTC完成，6分43秒；task1于01:32:39 UTC完成，6分59秒，各root、batch、extern均COMPLETED0:0，STATUS/JOBID一致。实际每task为private、0GPU、2CPU、8GiB、12分钟，Command、SubmitLine、WorkDir、array及日志、Requeue0和Restarts0相符。task0的JobIDRaw为13422816，task1为13422815，不将raw scheduler ID误当array ordinal。

| 人工势 | 完整decode秒 | 进程peak RSS KiB | 所选完整链 | 总gain |
| --- | ---: | ---: | ---: | ---: |
| 三窗全零 | 397.817506 | 505596 | 0 | 0 |
| 三窗独立整数势 | 413.921085 | 504196 | 88 | 21526 |

两个输入均为同一49152bp人工DNA，实际801个ATG、6102个GT/GC donor、2964个AG acceptor及2251个STOP，base上下文分别覆盖[0,18432)、[18432,30720)、[30720,49152)。复制准备分别0.046501与0.056127秒；Slurm batch采样MaxRSS分别412128K与419724K，不能替代上表进程峰值。两个case的summary和单行records均闭合、stderr均为空。源码及配置逐值核对，未改变原解码器，也没有通过zero快捷返回制造null。

两个case均有530939条缓存合法exon边、63995392次实际exon扩展、56314131次donor状态提交、54247318次已有同key donor历史的比较、2064146次intron历史比较、7252196次完成扩展、429065次cap边界排除及2667个最大active状态；每个计数逐值相等。源码中的donor_states_dominated名称包含已有槽位被比较、可能胜出替换的情况，不能把54247318都解释为被丢弃的新路径。非零case的88条链、各gain和总gain配对一致，最大58-exon；实际DNA、ORF、motif、span、cap和直接FP64整数势检查通过，总gain21526不超过逐base上界142551。全零case为预期null0，case1没有独立全长最优性oracle，也没有真实注释标签或生物学准确率结论。

这确认两个密集49kb输入能在原12分钟/8GiB范围闭合，但不是染色体可行性或公平速度benchmark。两次decode的时间差不能归因于潜势或回溯；与先前24576bp在不同时间/负载下的115秒也不能拟合或宣称二次运行增长。合法边只有约53万、历史扩展却约6400万，说明仅排除必定非法终点并不减少原DP中这些合法历史扩展。下一步先在同一已保存人工配置上核对独立索引版本的完整输出和计数，再据实测决定是否需要保持同一语法/目标的历史合并或批处理；不能因为三个小检查等价就认定整个计算瓶颈已解决。

本轮完整1491窗A42、原五条件、R4保护及唯一固定监控均未改；没有新训练、真实raw、DEV/独立物种/Setaria/A33执行或对外发布。本文记录的88条人工链与21526分不是模型性能或新基因证据，论文目标仍未达到。

## 同一人工配置的索引解码完整回放

原array13422815两项49152bp人工输入均已COMPLETED0:0，且STATUS、JOBID、summary及自然records闭合。新回放只替换为已通过13422866的独立索引解码实现：DNA seed0、三个原genomic窗口、case0全零及case1逐窗CPU seed41整数FP32势、base上下文、单位角色、FP64、gain0、null0和max128全部相同。保留原算法及两项旧产物，人工随机输入不重新挑选，也不加载真实缓存或权重。

新搜索先独立完成并将natural_records.jsonl关闭，随后才打开同case的原summary，核对完整人工配置、原job身份及不读研究数据/模型/标签的范围。比较全部链身份及顺序、各gain、总gain、selected数量、全部八项遍历计数、base上下文及pointwise上界，不用时长或RSS判定算法等价。新输出另以实际拼接DNA检查ORF、motif、cap、span和直接FP64整数角色势。case1仍不是独立全长穷举证明；相同人工答案也不是模型、真实注释或论文创新证据。

唯一命令sbatch --parsable sbatch/M28-A42-INDEXED-CONTIG-REPLAY.sbatch，private、array0至1且最多并行2，每task为0GPU、2CPU、8GiB、12分钟、generanno和no-requeue。根outputs/M28-A42-INDEXED-CONTIG-REPLAY/case_0或case_1，日志%x_%A_%a，只有短JSON及marker，不保存tensor或DNA。12分钟沿用已闭合原case的资源余量，不预先保证或宣称索引版更快；同一节点名不能控制共享负载、热状态或环境时序，单次时间差不作公平因果提速benchmark。

提交前核对原case及三个索引检查终态、相关源码、五项活调度、节点/预约和精确输出冲突。不重跑已闭合的小检查。结果身份、链、分数或计数不一致时保留并停止，不跳过条件；仅明确实现/fixture错误可同语法/评分/资源范围独立R1修正一次，再次失败、超时或OOM关闭本回放版本，不延长或改cap。A42完整1491窗协议、R4保护及固定监控不变；新任务只有人工数据，不打开DEV、独立物种、Setaria或A33。

实际array13422982仅提交一次，两task于2026-10-07 01:44:33 UTC提交/开始；task0于01:46:09结束、96秒，task1于01:46:25结束、112秒，root、batch、extern全部COMPLETED0:0。task0的JobIDRaw13422983、task1为13422982；各STATUS/JOBID正确，每个natural_records.jsonl恰一行、summary非空，两个stderr均0字节。实际private/gpu034、每task0GPU、2CPU、8GiB、12分钟，Command、SubmitLine、WorkDir、array/log、Requeue0及Restarts0相符。新脚本、SBATCH、记录及五个已存在依赖逐字与canonical一致，没有重写依赖；AST/语法与提交前五项活调度、节点/预约和冲突核对通过。

新搜索落盘关闭后才打开原答案，两case的配置/旧job范围身份、全部链与顺序、每个gain、总gain及八项计数均逐值相等，相关三个summary proof均true；实际DNA/ORF/motif/角色/cap/span和直接FP64整数分数也通过。case0仍null0，case1仍88条链、总gain21526、逐base上界142551。两case均保留530939条合法边、63995392次exon扩展等全部原计数，因此不是减少了合法历史路径或使用零势快捷返回。

| 相同人工case | 索引decode秒 | 原decode秒 | 索引进程peak RSS KiB | 索引batch采样MaxRSS |
| --- | ---: | ---: | ---: | ---: |
| 全零 | 89.921090 | 397.817506 | 575040 | 484448K |
| 逐窗整数势 | 105.846233 | 413.921085 | 580360 | 490916K |

这次完整49kb回放的同答案回归通过，且观察到用时更短、进程峰值更高；它不是受控公平benchmark，不能将单次时间比视为因果提速或推广到染色体。当前没有失败、R1、资源延长、模型forward/更新、真实研究输入/标签访问或准确率结论。原1491窗A42及固定监控未改，论文目标未完成。

## 长exon相位等价与历史合并的例外

这一步只有数学/人工准备，不修改或接线任何解码器。对同一acceptor、相同已完成exon数、相同非空carry长度p的中间/末exon，若长度至少3，前3-p个新base会完成incoming codon。完成codon必须为ACGT且非STOP；符合者之后的完整codon合法性、终点remainder及原角色增量只依赖p和新exon DNA，不依赖incoming carry的具体字母。失败者没有合法长度至少3的出口。相位0另成一类；first exon不参与此合并，不同exon数也不合并。

这一等价不能泛化到1–2bp中间exon：它可能保留尚未完成的actual carry，某条历史长出口非法时短出口仍合法，后续codon也可能与同相位另一carry不同。另一个具体限制是原FP64及严格首次到达tie：入口分数1与nextafter(1,+inf)严格不同，加同一个2后均舍入为3；原逐出口比较保留先到的较低入口历史，预先选较高入口会换链身份。因此“先每相位取最大入口”在实数上保持标量目标，不等于原FP64首次到达身份等价。没有证明及实现能处理该问题前，不将这种历史合并接入索引版或冻结A42。

long_exon_phase.py仅返回上述长出口相位/无合法长出口，未被任何生产解码器引用。四项检查覆盖全部64种首三base与21种actual carry、25种带N的随机人工序列、所有长度至少3的终点和terminal/nonterminal，以独立直接拼接DNA同时核对ORFIndex；合法同相位增量另按原角色算术逐值比较。还保留短中间exon反例及上述FP64首次到达反例。有限人工核对支撑这项局部等价和限制，不是全图最优性、合并后的算法证明或模型准确率。

唯一命令sbatch --parsable sbatch/M28-A42-LONG-EXON-PHASE-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue；独立根outputs/M28-A42-LONG-EXON-PHASE-CHECKS及%x_%j日志。提交前相关源码、AST、SBATCH语法、五项活调度与输出冲突核对；错误保留，明确实现/fixture错误仅可同范围资源独立R1一次，再次失败/超时停止，不修改原DP、原评分/阈值/cap/语法或1491窗A42。

实际13423011于2026-10-07 01:46:55 UTC提交/开始、01:47:01 UTC结束，6秒，root、batch、extern全部COMPLETED0:0，STATUS/JOBID一致。实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS451592K，Command/SubmitLine/WorkDir/log及Requeue0/Restarts0相符。四项检查PASS0.577秒，长exon局部相位等价、短exon必须保留actual carry及FP64预先选入口会改变首次到达身份均按预期核对。反例是检查的预期结果，不是作业失败；没有R1、资源延长、解码器接线或真实研究数据访问。三个新源码/SBATCH文件逐字与canonical一致，两个Python AST和SBATCH语法通过，五项活调度及输出冲突核对完成。

下一步仍优先保持原FP64操作及tie身份：可以研究按每个出口delta重新比较或批处理同相位历史，而不能直接使用上述已证实不等价的入口赢家。任何实际DP接线须先独立人工小图对照、cap/长exon/short-carry检查，再复用未挑选的两项49kb人工输入作完整回放；不改变已运行的冻结A42。本文只有计算准备进展，没有科学五条件、真实注释准确率或染色体可行性的新结论。

## 保留出口舍入与短carry的批处理解码

独立emission_lattice_batched.py只用于人工计算准备，原emission_lattice、indexed版本、contig/research调用及1491窗A42均不改。输入明确沿用CPU FP32的15列势，原_Margins、FP64角色前缀/差分和成熟intron归一化顺序不变；不增加beam、motif top-k、长度界、cap变化或零势快捷返回。first exon保持独立；1–2bp中间exon逐actual carry处理；合法长exon只在相同已完成exon数和相位内共享合法终点及角色delta，仍调用原ORFIndex验证。

长组保留原arrival行序，对每个出口将每条入口分数与同delta作numpy逐元素FP64加法，随后argmax保留首次同分行。并非先选最大入口，也没有减少这些FP64加法。减少的是被同key出口赢家替代的Python状态提交与chain构造。固定终点的terminal只有一个相位满足完整codon；非terminal的不同相位具有不同outgoing carry长度及DP key；短/长终点不重合。各长组按原active排序的首次出现顺序建立，完成exon数递增；因此不同组不会颠倒一个实际竞争同key的同分arrival。原跨position严格胜出、null0及gene-set更新保持。

新计数显式分开logical_exon_extensions等逻辑路径数、state_update_proposals实际提交数及short_actual_carry_extensions/long_exit_proposals，不把新物理计数冒充原donor_states_dominated等计数。链身份、顺序、每个gain和总gain必须与原DP逐值相等；另核对逻辑扩展/terminal/cap及最大active与原计数对应。此次只减少等价历史的状态提交，不是声称全长独立最优性、真实注释准确率或染色体资源可行性。

六项人工检查覆盖每出口FP64首次到达与逐Python标量的数值/行身份及实际DNA两carry汇合舍入tie、八个小序列的72组整数势及null独立DNA穷举、short actual carry/GT/GC/intron N/负前缀/邻接六gene、12个固定随机768bp中段上的72个整数或非整数FP32对照、原128-exon上界/cap127、两方向30000bp完整exon。完整identity及逻辑计数不一致即停止；不只比较最优标量，也不调用模型或打开研究数据。

唯一命令sbatch --parsable sbatch/M28-A42-BATCHED-LATTICE-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue，独立根outputs/M28-A42-BATCHED-LATTICE-CHECKS及%x_%j日志。提交前源码/AST/SBATCH、五项活调度、节点/预约和精确输出冲突核对。明确实现/fixture错误仅可同语法/目标/资源独立R1一次；再次失败、超时或OOM关闭本准备版本，不延长、改cap、跳过条件或接线原A42。六项闭合后才进行原两项49kb配置的完整回放，回放仍先完成新自然输出再打开原答案，沿用既有12分钟/8GiB范围及独立新根。

实际13423096于2026-10-07 01:54:50 UTC提交、01:54:53开始、01:55:10结束，17秒；root、batch、extern全部COMPLETED0:0，STATUS/JOBID正确。实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS476856K，Command/SubmitLine/WorkDir/log、Requeue0/Restarts0相符。六项PASS10.687秒，完整链身份/顺序/gain与原版逐值相等，四个逻辑计数对应相等；实际DNA上的较低入口首次到达舍入tie保持。小图独立DNA穷举、非整数FP32、short actual carry、cap127/128及两方向30000bp完整exon均通过。五个相关文件逐字与canonical一致，新增两Python AST及SBATCH语法通过；五项活调度、节点/预约及冲突核对通过，没有R1/资源延长/真实数据访问。long_exon_phase只更新说明以注明该人工调用，语义不变；冻结A42与原DP不改。

## 批处理版同人工配置完整回放

源仍为原array13422815的case0/1，49152bp DNA seed0、原三窗、plus、case0全零或case1CPU seed41逐窗FP32整数势、base上下文、FP64、单位角色、max128/threshold0/null0全部保持。没有重抽输入或改变科学评分；旧产物与原DP保留。replay_batched_contig_synthetic.py仅在六项检查13423096及原case的COMPLETED/JOBID相符后运行；新自然搜索先落盘关闭，之后才打开原summary逐值比较链身份、顺序、各gain、总gain、selected数量、上下文和逐base上界。

逻辑扩展/terminal/cap与最大active必须对应原计数；物理state_update_proposals等按新定义记录，不要求不同算法的物理计数与旧donor比较数相等。另以实际DNA/拼接CDS/角色直接核对每条完整链、motif/ORF/cap/span和FP64整数势。减少状态提交只是一项计算事实，不等于减少FP64加法、目标更好、生物学新gene或模型准确率。无独立全长穷举，人工回归也不支持真实染色体可行性。

唯一命令sbatch --parsable sbatch/M28-A42-BATCHED-CONTIG-REPLAY.sbatch，private、array0至1%2，每task0GPU、2CPU、8GiB、12分钟、generanno、no-requeue。独立根outputs/M28-A42-BATCHED-CONTIG-REPLAY/case_0或case_1，日志%x_%A_%a，短JSON/marker，不保存DNA/tensor。新实现墙时尚未测定；12分钟保持原人工探测资源边界，不承诺更快，也不从一次时间差作公平因果benchmark。

提交前六项检查的root/batch/extern及marker、源文本/AST/SBATCH、五项活调度、节点/预约及精确输出冲突核对。链/分数/逻辑计数不一致停止保留；明确实现/fixture错误同语法/目标/资源独立R1一次，再次失败/超时/OOM关闭，不延长/改cap/跳过case。冻结1491窗A42与唯一固定监控不变，没有新训练、真实raw/DEV/独立物种/Setaria/A33访问、新label或对外发布。

实际array13423130于2026-10-07 01:58:02 UTC仅提交一次，两task于01:58:03在gpu034开始；task0于01:59:34结束、91秒，task1于01:59:50结束、107秒。各root/batch/extern全部COMPLETED0:0，STATUS/JOBID正确，task0 JobIDRaw13423132、task1为13423130。实际每task为private、0GPU、2CPU、8GiB、12分钟，Command/SubmitLine/WorkDir/array/log、Requeue0和Restarts0相符。没有失败、R1、资源延长或接线原A42。提交前五项活调度/预约/输出冲突、replay AST/SBATCH及三项新增文件和两项原人工配置依赖逐字同步核对通过。

两个case均在新自然输出关闭后才读取旧答案，四项summary proof全部true：配置及原job范围、自然搜索闭合顺序、全部链身份/顺序/各FP64 gain/总gain，以及逻辑扩展/terminal/cap/active与原版相符。case0仍null0，case1仍88条完整链、总gain21526、pointwise上界142551；实际DNA/ORF/motif/角色/cap/span与直接FP64整数势检查通过。人工链数及gain不代表模型准确率或生物学新gene。

| 同人工case | batch decode秒 | indexed decode秒 | batch进程peak RSS KiB | batch采样MaxRSS |
| --- | ---: | ---: | ---: | ---: |
| 全零 | 85.691570 | 89.921090 | 525504 | 441776K |
| 逐窗整数势 | 101.870976 | 105.846233 | 527160 | 437532K |

两项均保留63995392次逻辑exon扩展、7252196次逻辑terminal扩展、429065次cap边界排除及2667个最大active。新实际state_update_proposals为11147780，约为原逐历史处理次数的17.42%；这包括出口候选提交，不是成功写入状态数。donor候选9983501、已有donor槽位比较7916688、成熟intron key比较2064146；first缓存边7253、short actual-carry缓存边14147、共享long缓存边77072，不与原逐carry合法边530939混为同一定义。short实际扩展1735773、long逻辑历史扩展62252366与first7253之和恰为63995392；long出口候选9404754与short1735773及first7253之和恰为11147780。全部两case计数逐值一致。

新实现仍完成全部逐历史FP64加法，没有靠剪枝降低逻辑路径数。虽然Python候选提交减少约82.58%，此次完整decode墙时只比上一次indexed短约4秒；两次运行负载/环境时序未受控，不能据此宣称公平提速，更不能把候选数比当作速度比或外推染色体。进程峰值与Slurm采样分别如表，采样不是进程峰值。

实际记录还显示1088111个long_history_groups，其中只有存在合法长出口的组才逐个创建小numpy数组并argmax；旧计数没有单独记录非空组或helper调用次数，不能把总组数直接当调用数。按同acceptor相位合并多exon-count组的数组调用仍是可研究的方向，是否为主要墙时瓶颈尚未通过计时剖析证明。此前整gene列表复制也仍在，当前不是其成本归因实验。下一版若改变分组粒度，继续保留actual carry、每出口FP64加法/首次到达、exon count及全部回放条件，不改冻结1491窗A42或打开研究数据。

## 跨exon count的相位数组块

独立emission_lattice_phase_block.py沿用已闭合batch的语法、CPU FP32输入、原FP64评分/成熟intron归一化、first与short actual carry、max cap和严格胜出。只把同acceptor同相位的多个exon-count组放到一组数组调用中，各count独立为行、actual carry按原字母arrival为列；缺失槽为-inf占位，不成为路径或zero/null候选。每个出口仍先对每条实际入口作FP64加法，再逐carry严格比较，保留首次同分。使用二维count×endpoint工作数组逐carry更新，不生成count×carry×endpoint三维中间数组。

不同count仍逐组写回并按原count递增处理；不同相位在固定terminal终点只有一个可完成，在donor终点有不同carry长度和key，因此把数组调用按相位归块不重排实际竞争同key的同分arrival。整gene-set回溯复制保持，不与这次数组粒度变化混改。原DP、indexed、batch、contig/research调用及1491窗A42均不改。

六项检查复用上一版人工fixture，三方比较原DP、逐组batch和phase block的所有链身份/顺序/gain；上一版全部物理和逻辑计数必须逐值相等，只另记录long_array_blocks、long_carry_column_batches与long_nonempty_history_groups，区分实际块调用、carry列batch和非空历史组。新的数值检查包含多count稀疏carry列、首列缺失、FP64舍入tie及实际DNA两carry汇合。其余仍为72组小图整数势及null独立DNA穷举、short actual carry/GT/GC/N/负前缀/邻接gene、72组768bp中段整数或非整数FP32、128-exon上界/cap127和两方向30000bp完整exon；不单独重跑未改的旧检查作业。

唯一命令sbatch --parsable sbatch/M28-A42-PHASE-BLOCK-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue，根outputs/M28-A42-PHASE-BLOCK-CHECKS及%x_%j日志。提交前源文本/AST/SBATCH、活调度五项、节点/预约和精确冲突核对。明确实现/fixture错误仅可同语法/目标/资源独立R1一次，再次失败/超时/OOM关闭，不延长/改cap/跳过身份或计数条件。通过后才做原两项49kb人工配置完整回放，保持原12分钟/8GiB资源，先新自然输出闭合再读旧答案；没有模型forward/训练、研究输入/标签/DEV/独立物种/Setaria/A33访问或对外发布。

实际13423246于2026-10-07 02:09:53 UTC提交、02:09:56开始、02:10:17结束，21秒。root、batch、extern全部COMPLETED0:0，STATUS/JOBID相符；实际private、gpu034、0GPU、2CPU、8GiB、2分钟，batch MaxRSS482736K，Command/SubmitLine/WorkDir/log及Requeue0/Restarts0相符。六项PASS14.639秒，三方全部链身份/顺序/gain一致，旧batch全部物理/逻辑计数相等；多count稀疏carry列、首列缺失、FP64首次到达及实际DNA舍入tie保持，小图独立DNA穷举、整数/非整数FP32、short carry、cap和两方向长exon通过。原源码不改，定向差异仅为数组块/helper及新增计数，AST/SBATCH、五项活调度/节点/预约/冲突和四个目标文件及原batch测试逐字同步核对通过；计数说明修正后又逐字同步三个改动文件。没有R1、资源延长、研究数据访问或接线冻结A42。

## 相位块版同人工配置完整回放

两项仍是原array13422815的49152bp人工DNA seed0、三窗、plus、全零或逐窗CPU seed41整数FP32势；base上下文、单位角色、FP64、cap128/threshold0/null0不变。replay_phase_block_contig_synthetic.py要求六项检查13423246、原case和上一版batch13423130的COMPLETED/JOBID身份相符，只读这些轻量marker。新自然搜索落盘关闭之后才打开原DP和batch的summary，逐值核对配置/人工范围、全部链/顺序/各gain/总gain、selected数、上下文和pointwise上界。

原四项逻辑计数对应、上一版batch的全部物理和逻辑计数必须相等；新增三个计数区分总历史组、非空历史组、数组块调用和carry列batch，不把所有group视为NumPy调用。占位槽也参与数组运算但不成为实际历史，逻辑路径数不是所有padding加法数。独立直接DNA/角色核对完整CDS、motif、cap/span及FP64整数势。人工同答案不证明全长穷举最优性、真实染色体可行性或准确率。

唯一命令sbatch --parsable sbatch/M28-A42-PHASE-BLOCK-REPLAY.sbatch，private、array0至1%2，每task0GPU、2CPU、8GiB、12分钟、generanno、no-requeue；根outputs/M28-A42-PHASE-BLOCK-REPLAY/case_0或case_1及%x_%A_%a日志，只短JSON/marker，不保存DNA/tensor。运行时间待实际观察，不以一次时间差宣称公平提速。提交前六项检查与两组源case终态、源文本/AST/SBATCH、五项活调度/预约/节点和精确输出冲突核对。明确实现/fixture错误可同语法/目标/资源独立R1一次，再次失败/超时/OOM关闭，不延长/改cap/跳过答案或计数比较。冻结1491窗A42及固定监控不变，研究数据/标签/DEV/独立物种/Setaria/A33均不访问，没有新训练或发布。

实际array13423292于2026-10-07 02:13:49 UTC仅提交一次，两task于02:13:50在gpu034开始；task0于02:15:21结束、91秒，task1于02:15:34结束、104秒。各root/batch/extern全部COMPLETED0:0，STATUS/JOBID一致；task0 JobIDRaw13423293，task1为13423292。实际每task为private、0GPU、2CPU、8GiB、12分钟，Command/SubmitLine/WorkDir/array/log及Requeue0/Restarts0相符。两个summary非空，每个natural_records.jsonl恰一行，stderr均0字节。replay脚本/SBATCH/记录及相关原依赖逐字同步，AST/SBATCH、提交前五项活调度和节点/预约/精确冲突核对通过；没有R1、资源延长或改写原产物。

两case的四项summary proof均true：原配置及两组旧job身份、新自然输出先闭合再读取两组旧答案、全部链身份/顺序/各FP64 gain/总gain，以及原四项逻辑计数和上一版batch全部物理/逻辑计数均相等。实际完整CDS的DNA/ORF/motif/角色/cap/span及直接FP64整数势检查通过。全零仍null0；整数势仍88条完整链、gain21526、pointwise上界142551。计数仍含63995392次逻辑exon扩展、11147780次状态候选提交和2667个最大active，没有缩减搜索空间或改变评分。

| 同人工case | phase-block decode秒 | 前版batch decode秒 | phase-block进程peak RSS KiB | phase-block采样MaxRSS |
| --- | ---: | ---: | ---: | ---: |
| 全零 | 85.328378 | 85.691570 | 529732 | 439248K |
| 逐窗整数势 | 98.748422 | 101.870976 | 530000 | 440308K |

两case新增计数均为8888个long_array_blocks、59209个long_carry_column_batches及1088111个long_nonempty_history_groups。此次全部long_history_groups恰为非空；这说明相同组结构下逐组实现会分别处理1088111个非空组，但不是对旧运行NumPy调用的独立计时测量。新块数比非空组数少约99.18%，carry列batch不是总NumPy调用数或运算数；缺失carry槽的-inf占位也参与数组算术但不成为实际路径。不能把组数/块数比当作提速倍数。

本次decode只比前版短0.363192及3.122553秒，环境/共享负载未受控，未显示可归因的显著墙时改善；进程峰值略高，Slurm采样仍不替代进程峰值。该结果支持同答案计算接线，不证明染色体可行性、准确率、生物学新gene或论文创新。下一轮应先对同一人工配置作有界成本剖析，区分Python状态/回溯、合法边/角色delta、数组构造和数组算术；有具体主要成本后才继续实现，避免仅凭调用计数再叠加计算层。原1491窗A42与固定监控保持不变，研究数据/标签/DEV/独立物种/Setaria/A33未访问，论文目标仍未完成。

## 同人工配置的有界函数成本剖析

profile_phase_block_contig_synthetic.py用标准库cProfile运行未改动的replay_phase_block_contig_synthetic.main；解码器、原回放、人工配置、FP64操作及语法均不改，也不重复六项已闭合小检查。两个case仍为49152bp DNA seed0、原三窗、plus及全零或CPU seed41逐窗整数FP32势，cap128/threshold0/null0/单位角色与原base上下文不变。原13423292两case的COMPLETED/JOBID相符后开始；原回放保持新自然输出先落盘关闭、再打开原DP和batch答案的顺序，随后才打开旧phase-block summary，全部配置/范围、链身份/顺序/各gain/总gain和全部phase-block计数须逐值相等。

剖析覆盖输入重建、上下文复制、完整decode、独立DNA/角色核验和JSON，不包含模块导入。保存小型replay.pstats、profile_summary.json和原回放JSON/marker，记录函数自身与累计时间、调用数、前25项及相关decoder函数；累计时间互相包含，不能相加为总成本。用decode函数累计时间区分外围工作，以函数自身时间定位update、choices、logical_count、long_exon_phase及数组helper等候选成本。单函数内部的tuple复制/状态比较或数组构造/算术未分别计时，不能从函数总时间直接给其中一行作归因；若此项影响下一步，另做有界定向剖析，不据猜测改回溯或数组实现。

唯一命令sbatch --parsable sbatch/M28-A42-PHASE-BLOCK-PROFILE.sbatch，private、array0至1%2，每task0GPU、2CPU、8GiB、12分钟、generanno、no-requeue；根outputs/M28-A42-PHASE-BLOCK-PROFILE/case_0或case_1，日志%x_%A_%a。原49kb同配置未剖析decode约85至99秒，剖析会增加且偏置函数时间；12分钟沿用既有人工资源上限，不保证完成，也不将剖析时间视为未插桩速度或公平benchmark。只有成本排序与具体高调用热点会改变下一项实现选择。

提交前核对原phase-block两项终态/身份、相关源文件、AST/SBATCH、五项活调度、节点/预约及精确输出冲突。身份/分数/计数变化保留并停止；仅明确实现/fixture错误可同目标/语法/资源独立R1一次，再次失败、超时或OOM关闭此诊断，不延长/改cap/跳过case。没有新训练、模型forward、研究输入/标签/DEV/独立物种/Setaria/A33访问、新科学门槛或对外发布；冻结1491窗A42及唯一固定监控不变。

首版array13423384于2026-10-07 02:23:14 UTC提交、两task于02:23:15在gpu034开始；task0于02:25:36结束、141秒，task1于02:25:52结束、157秒，root/batch/extern全部COMPLETED0:0，STATUS/JOBID正确。task0 JobIDRaw13423385，task1为13423384；每task实际0GPU、2CPU、8GiB、12分钟，batch MaxRSS448572K和431132K，Command/SubmitLine/WorkDir/array/log及Requeue0/Restarts0相符。两套summary/profile_summary、单行natural_records及pstats非空，stderr均0；三项新profile proof及原回放proof通过，全部phase-block链/顺序/FP64 gain/计数相等，没有改变自然答案。

但首版成本总表存在实际统计缺口：标准pstats按filename、line及co_name索引；多个dataclass生成的__init__同为<string>第2行，cProfile.snapshot_stats将相同标签覆盖。首版自身时间合计124.624278和136.986227秒，小于decode累计135.255979和149.301935秒，不能当作完整总成本或据此归一化占比；首版pstats和profile_summary原样保留，不用不完整排行决定状态对象优化。该缺口不影响已逐值核对的自然答案，也不将工程COMPLETED改称科学FAIL。

唯一R1只修剖析统计：直接读取profiler.getstats的独立code entry及原self/cumulative/call计数，将_Live、_Solution及Budget构造函数明确命名，保存所有raw函数条目，并记录标准标签碰撞的实际余项数。不再将pstats转换后的总表作为依据或输出新的有损pstats。原回放、DP、人工配置及资源不变；R1唯一命令sbatch --parsable sbatch/M28-A42-PHASE-BLOCK-PROFILE-R1.sbatch，private、array0至1%2、每task0GPU/2CPU/8GiB/12分钟/generanno/no-requeue，独立根outputs/M28-A42-PHASE-BLOCK-PROFILE-R1/case_0或case_1及%x_%A_%a日志。再次统计/身份/分数/计数失败、超时或OOM关闭，不做R2或改科学范围。

R1实际array13423455于2026-10-07 02:29:56 UTC仅提交一次，两task于02:29:57在gpu034开始；task0于02:32:19结束、142秒，task1于02:32:35结束、158秒，各root/batch/extern全部COMPLETED0:0。STATUS/JOBID正确，task0 JobIDRaw13423456，task1为13423455；实际每task0GPU/2CPU/8GiB/12分钟，batch MaxRSS432480K和453164K，Command/SubmitLine/WorkDir/array/log及Requeue0/Restarts0相符。各summary与保存完整raw条目的profile_summary非空，每份natural_records恰一行，stderr均0。三项profile身份/闭合顺序/全部答案与计数proof及raw条目proof均true，原回放四项proof仍true；全零仍null0，整数势仍88条链/gain21526，没有新自然答案或计数变化。

两case的standard_label_collision_surplus均为3，除生成构造方法外，同一行的两个list comprehension也共用标准标签。完整raw自身时间合计136.781367与153.011688秒，decode累计134.948434与150.526413秒，不再出现完整总表反小于其子过程的缺口。插桩decode墙时为135.004966与150.578099秒；与未插桩85.328378/98.748422秒不是可比速度实验，R1也不是解码器提速。

| 函数 | 全零调用数 | 整数势调用数 | 全零自身秒 | 整数势自身秒 |
| --- | ---: | ---: | ---: | ---: |
| decode_phase_block_emission_lattice | 1 | 1 | 33.416455 | 35.365016 |
| logical_count | 8690937 | 8690937 | 27.049399 | 30.969877 |
| update | 11147780 | 11147780 | 18.880669 | 21.410628 |
| _Live.__init__ | 11736452 | 15469140 | 11.549962 | 15.338220 |
| long_exon_phase | 7602025 | 7602025 | 9.136707 | 9.771222 |
| _phase_block_winners | 8888 | 8888 | 2.639574 | 2.752995 |

logical_count的累计时间为34.578750/38.473137秒，long_exon_phase为10.794444/11.411077秒，而数组helper为5.146497/4.900510秒。该插桩结果将下一步从继续调数组块改为消除重复整数计数和同acceptor同actual-carry的DNA合法性判断；函数调用较多会受到剖析偏置，数值不是未插桩成本占比或可期待的提速秒数。update累计24.799715/30.752369秒与_Live构造成本存在嵌套，不能相加；解码主函数自身时间也不能全部归因于某一行或回溯复制。

源语义支持一项直接改进：边表固定后，逻辑exon/terminal/cap计数可按实际历史数乘以对应边数汇总，仍保留每个exon-count状态；相同acceptor及actual carry的long_exon_phase只依赖固定DNA，可在该位置内复用，包括None结果。只有整数记账与纯DNA元数据复用改变，first/short的实际出口处理、每出口FP64加法、严格首次同分、状态提交/回溯和全部既有计数必须保持。下一轮先在独立人工实现做现有六类对照，再在原两项49kb配置中核对全部链/分数/计数；不接线冻结A42或改变科学门槛。此轮没有R2、资源延长、模型更新、研究输入/标签/DEV/独立物种/Setaria/A33访问或对外发布，论文目标未完成。

## 纯DNA元数据复用与整数计数汇总

独立emission_lattice_cached_metadata.py仅将不依赖历史分数的工作复用：固定边表保存edge数量和terminal数量，同acceptor同actual carry的short边表及long_exon_phase结果只取一次，包括None；该位置内累计实际历史数和到cap的历史数，按边数乘历史数记账。long计数在相同phase块内汇总，但各exon-count组仍各自保留，全部实际first/short/long出口按原顺序处理。FP64角色前缀/差分、成熟intron归一化、每出口入口加delta、严格首次同分、状态提交与整gene-set回溯均不改，直接复用原_phase_block_winners。无beam、top-k、长度界、zero快捷返回或科学评分变化。

旧phase-block全部18项物理/逻辑计数必须相等；另记录long_phase_metadata_evaluations实际纯DNA判断次数和logical_counter_batches整数汇总次数，不能当作逻辑路径数、全部运算数或速度倍数。acceptor内元数据表最多覆盖语法的21种actual carry，位置结束后不跨acceptor复用；边表仍按原索引/ORF/角色检查构建。原DP、indexed、batch、phase-block、contig和研究调用、已关闭产物及运行中的1491窗A42保持不变。

人工检查继承现有六类，对照原DP、batch、phase-block与cached版本的完整链身份/顺序/每个gain/总gain，以及全部旧phase-block计数：FP64舍入和稀疏carry列、小图独立DNA穷举、short actual carry/GT/GC/N/负前缀/邻接gene、768bp整数及非整数FP32、cap127/128、两方向30000bp完整exon。另加实际三exon人工链ATGT加GTTAG加A加GTTAG加CAAATAA：carry T在acceptor9的长出口为None，但1bp short形成TA，随后完整CDS合法；cap3必须保留指定链，cap2必须排除。每次纯DNA判断次数不超过该序列acceptor数乘21，检查确实复用而非逐历史重算。

唯一命令sbatch --parsable sbatch/M28-A42-CACHED-METADATA-CHECKS.sbatch，private、0GPU、2CPU、8GiB、2分钟、generanno、no-requeue；根outputs/M28-A42-CACHED-METADATA-CHECKS及%x_%j日志。提交前相关源文本/定向差异、AST/SBATCH、五项活调度、节点/预约及精确输出冲突核对。身份/分数/计数不一致保留并停止；仅明确实现/fixture错误可同目标/语法/资源独立R1一次，再次失败/超时/OOM关闭，不延长/改cap/跳过条件。七项闭合后才重建未挑选的原两项49kb人工配置完整回放，先新自然输出关闭再读旧答案，沿用12分钟/8GiB和独立新根，不打开真实研究输入、标签、DEV、独立物种、Setaria或A33。

实际13423604于2026-10-07 02:43:18 UTC提交、02:43:19开始、02:43:42结束，23秒；root/batch/extern全部COMPLETED0:0，STATUS/JOBID相符。实际private/gpu034、0GPU/2CPU/8GiB/2分钟，batch MaxRSS501948K，Command/SubmitLine/WorkDir/log及Requeue0/Restarts0相符。七项PASS18.073秒，四个版本的完整链/顺序/每个FP64 gain/总gain及旧18项计数逐值相等，独立DNA小图、非整数FP32、cap、两方向长exon通过；明确None长相位下指定1bp short链cap3保留、cap2排除。四个目标文件逐字同步canonical，两Python AST/SBATCH与定向差异、五项活调度/节点/预约/冲突核对通过；原依赖源码不改，没有R1、资源延长或研究数据访问。

## 元数据缓存版同人工配置完整回放

replay_cached_metadata_contig_synthetic.py要求七项检查13423604、原array13422815和phase-block13423292对应case的COMPLETED/JOBID相符，只读这些marker后重建原49152bp DNA seed0、原三窗、plus、全零或CPU seed41逐窗十五列整数FP32势。base上下文、单位角色、FP64、cap128/threshold0/null0不变，不重抽有利输入。新自然搜索落盘关闭后才打开原DP与phase-block summary，逐值比较完整配置/人工范围、全部链身份/顺序/各gain/总gain/selected数/上下文/逐base上界，以及原四项逻辑计数对应和全部旧18项phase-block物理/逻辑计数。两项新增计数仅反映元数据判断与整数汇总次数，不用它们替代搜索路径数。

每条新完整链另做独立拼接DNA/ORF/motif/角色/cap/span与直接FP64整数势核验，全零必须null0，总目标不得超过逐base上界。长输入没有独立全长穷举，有限同答案回归不证明模型准确率、染色体可行性或论文创新。实际墙时/进程peak RSS/Slurm采样分别记录，不以与既往单次时间差作公平因果提速。

唯一命令sbatch --parsable sbatch/M28-A42-CACHED-METADATA-REPLAY.sbatch，private、array0至1%2，每task0GPU/2CPU/8GiB/12分钟/generanno/no-requeue，根outputs/M28-A42-CACHED-METADATA-REPLAY/case_0或case_1及%x_%A_%a日志，只短JSON/marker。提交前检查与原两套case终态/身份、相关源文件/AST/SBATCH、五项活调度、节点/预约和精确输出冲突核对。链/分数/计数/身份变化保留并停止；仅明确实现/fixture错误可同目标/语法/资源独立R1一次，再次失败/超时/OOM关闭，不延长/改cap/跳过case。1491窗冻结A42与固定监控不变；无模型forward/更新、研究数据/标签/DEV/独立物种/Setaria/A33访问或发布。

实际array13423648于2026-10-07 02:48:06 UTC仅提交一次，两task于02:48:07在gpu034开始；task0于02:49:15结束、68秒，task1于02:49:25结束、78秒，root/batch/extern均COMPLETED0:0。STATUS/JOBID正确，task0 JobIDRaw13423649，task1为13423648；实际0GPU/每task2CPU/8GiB/12分钟，batch MaxRSS442744K和457816K，Command/SubmitLine/WorkDir/array/log及Requeue0/Restarts0相符。两份summary非空、natural_records各一行、stderr均空；四项身份/先关闭新搜索/全部答案/旧计数proof均true，独立DNA、角色直接FP64、cap/span与上界检查通过。

| 人工case | 缓存版decode秒 | 进程peak RSS KiB | 旧phase-block秒 | 已选链数 | 总gain |
| --- | ---: | ---: | ---: | ---: | ---: |
| 全零势 | 61.483640 | 534340 | 85.328378 | 0 | 0 |
| seed41整数势 | 71.883514 | 530440 | 98.748422 | 88 | 21526 |

两case的全部链身份/顺序/每个FP64 gain/总gain与原DP及phase-block答案相等，原四项逻辑对应及全部旧18项物理/逻辑计数相等。新增long_phase_metadata_evaluations为62055（原剖析相同函数调用7602025），logical_counter_batches为71744（原logical_count调用8690937）；原逻辑exon_extensions63995392、state_update_proposals11147780和maximum_active_states2667不减。减少的是重复纯DNA判断与整数记账，不是删除搜索路径、状态合并或所有操作；历次单次不同负载墙时不作公平速度比。保留独立人工缓存实现，暂不接入真实A42或宣称注释准确率。人工回放没有R1、资源延长或研究输入访问；随后主对话收到真实A42汇总故障并按其独立协议处理，不能将人工范围声明套给该真实结果处理。
