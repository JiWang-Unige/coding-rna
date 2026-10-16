# A26 候选 CDS 相位状态兼容性证据

2026-10-06，首次预检及拟合前固定。A21在原自然池改善正确链选择，但两物种known_nonmatch增加，原五gate失败。本次是独立表征假设：显式加入候选完整CDS路径对冻结局部状态预测的兼容性，能否改善整链选择组合。它不重训A21，不改变A23/A25生成器实验，也不声称旧hidden完全不含这些信息。

## 原始先例与研究边界

[Tiberius原论文](https://academic.oup.com/bioinformatics/article/40/12/btae685/7903281)第2.4节已有局部15类概率与读框／边界HMM状态的一一对应，以及结构化合法路径解码。因此“局部状态概率、相位一致性和路径证据”不是新概念。A26只检验在原自由候选和unknown保护的集合选择上加入该证据；没有端到端HMM fine-tuning，也不将本文作者的哺乳类准确率移植到当前植物或称为本项目SOTA比较。

## 唯一新表征

保留fresh A21 JunctionChainHead全部七分支和构造顺序。新增无bias Linear28→1，权重初值0，共28参数，总248096。原共有fresh seed0参数逐tensor相同；初始分数和原集合目标与fresh A21相同，不继承任何已训练head。新项加入nonadditive，记录原recurrent nonadditive与phase evidence term，gain仍为additive+nonadditive−原window null。不改变null、集合loss、生成器、候选或操作点。

冻结原A3末004608及其共享上下文／15类segmentation分支。候选坐标在原定向半开窗口中，按本项目内部ascending coding phase生成角色：CDS1/2/3、按上一编码相位的intron4/5/6、donor7/8/9、acceptor10/11/12、start13、stop14。累计相位跨外显子而不在新外显子重置；first/last连续3base覆盖端点角色，遵循原fine label优先级。相触外显子沿用旧头几何并按相同覆盖规则处理；当前自然候选语法不支持split start/stop。角色只由候选几何计算，推理特征不读reference、标签、catalog或原A21分数。

对state1至14分别计算候选角色base上的mean(logP(state)−logP(intergenic0))及log1p(base count)，合计28维；无该角色时mean0/count0。log-odds直接等于原logit差，不加epsilon、温度或二次概率归一化。候选span外及padding不贡献证据。局部logits与全部特征冻结；该项只是学权重的显式兼容性输入，不是已校准基因存在概率，也不是完整HMM路径概率。

## 固定数据曝光及门槛

原A15自然投影／A12核验、1491唯一TRAIN窗、16321自由候选配对，1/0/-1支持2380/12566/1375；原1536 draw顺序和importance权重三遍4608更新不变。冻结原A3 generator及GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345，四块D2048缓存、H128、FP32、fresh seed0、AdamW lr3e-4/wd.01/betas(.9,.999)/eps1e-8/clip1不变。原集合NLL/known候选数不变；-1不进gold/partition，其独立gain直接梯度0，不转为难负类。未知链仍与自然1/0共同参加阈值0的原span非重叠选择。

只有末004608读出，逐物种配对A6及已闭合A21，不扫描阈值、margin、epoch、seed或beam。事前五项gate仍为：两物种selected correct分别严格高于A6；两物种selected known_nonmatch分别不高于A6；pooled至少5exon selected correct严格高于A6。FAIL关闭候选不救参；PASS只支持另固定完整raw-FASTA caller/DEV协议，不自动运行DEV或发布。A21辅助配对不替代原A6门槛。TRAIN候选／窗计数不是独立基因F1或泛化。

## 必要验证与执行

五项合成测试核对跨exon phase、donor/acceptor/start/stop优先级及正反方向R4渲染；逐state精确log-odds/count、padding/span外排除及逐base常数平移不变；fresh共有参数/零项分数一致；局部状态变化实际进入完整链score且局部分支无梯度；空池、相触和unknown直接梯度。随后沿用原metadata八TRAIN窗，optimizer0，核对八窗中全部自然候选的有限28维证据、同fresh A21分数/loss、七共有分支按真实连接支持的梯度、新投影两物种有限非零梯度、generator无梯度以及head参数不变。零更新只验证工程接线，不判性能。

预检1RTX3090至少24GiB/2CPU/16GiB/10分钟、generanno。固定拟合/末读出上限60分钟，CPU核验10分钟，独立M28-A26-PHASE-EVIDENCE各阶段目录。拟合计划复用原A21不可变1491个FP32 hidden GPU缓存18,761,121,792字节，并在首次materialization按相同原候选顺序计算28维证据，额外只保留候选向量而非全量15类logits。新head的live/cached八窗score/loss/梯度等价须实际验证，不能沿用A21收据。每256工程checkpoint，不增加曝光；完整runner和新缓存峰值须实现/核验后才提交正式fit。

没有新GLM前向、候选搜索或reference注入，没有DEV/sealed/Setaria、作业取消、Git或公开发布。若等待，复用原固定watcher/heartbeat，不新建监控。完整raw-FASTA、同输入原生竞争性比较、独立迁移和相称消融仍属于未完成论文目标，不能以本组件的工程完成替代。
