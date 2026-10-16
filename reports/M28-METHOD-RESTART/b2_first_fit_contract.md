# M28 B2 首个拟合与评分位置对照合同

记录日期2026-10-01，在B2任何optimizer更新及新DEV性能之前固定。本轮尚未执行正式拟合。R6已验证两物种均存在合法真实前缀被实际beam或final截断的见证，B2准备已完成；这些结果支持测试一个新机制，不证明完整方法性能，也不自行授予执行权限。

科学问题是：保持同一冻结局部来源、端点、连接、语法及候选预算，新增训练的绝对前缀兼容性用于剪枝前，是否比仅用于固定候选池后的完整链评分，改善真实自由覆盖和完整基因链输出。本轮是单seed探索性位置机制对照，不是两个独立训练复现，也不是R5的延时、K扩展或阈值补救。

## 唯一新学习机制

冻结outputs/M28-PILOT-R5/B1/step_004608.pt的全部参数，包括shared表示、segmentation、endpoint和link头。基座eval、requires_grad=false，只有新增PrefixValue128可训练，共184,369参数；不重新拟合GLM或B1，不继承旧B1 chain-vs-null读出为新增可学习初值。

新头用exon编码、GRUCell前缀状态、entry/donor/final事件特征及几何形成绝对兼容性logit。没有旧proposal累计分的残差叠加；旧proposal分仅保留作控制与来源记录。null固定0，类均衡BCE的logit不能解释为自然候选正确概率。

## 数据与更新规则

只使用现有拟南芥和水稻TRAIN及DEV allowlist。冻结GENERanno revision b0483c23b6b63787b61a6d3a204a9b517d6ba345和原缓存，不复制159GiB特征。train与DEV窗口仍为24,576bp、stride12,288、terminal anchored及原定向规则。test和Setaria继续封存。

原1,536 TRAIN paired draws、原顺序及importance weight重复3遍，batch1，恰好4,608次optimizer更新。有放回及三次重复不称独立样本。新头在冻结基座构造后以torch.manual_seed(0)初始化；唯一正式checkpoint为第4,608次更新后的状态，不按DEV择epoch。只有新头交给AdamW，lr3e-4、weight_decay0.01、betas(0.9,0.999)、eps1e-8、global norm clip1.0，FP32与固定学习率，不扫seed或学习率。

标签政策完整继承R4，primary/支持资格、其他isoform及未知结构保护不变。训练元数据按TRAIN split筛选，再解析该窗口标签；不读DEV特征/参考用于损失。缓存内旧标签不用作监督，按window ID接入R4标签。

## 训练监督与搜索隔离

训练自由池由冻结的原加性DAG生成，FirstCutObserver只观察实际自由搜索组，不改变排序、预算或候选。该首个拟合为固定控制策略上的off-policy前缀学习，不同时引入on-policy搜索改变。

每窗损失为完整链balanced BCE与首次真实截断的配对softplus排序损失之和，系数各1；分别按每窗既有规则归一，再乘原importance weight，不除以batch权重和。完整正链可注入完整链损失池，来源单独标注，不能伪装为自由生成。前缀配对仅来自参考前缀第一次真实被beam或final截断的事件；端点、exon、link或ORF先阻断时不注入虚构路径。

beam比较实际存活的两名竞争者，final比较前两名及cutoff处存活者，重复身份去重。其他isoform精确前缀、未知或非callable跨度标签为-1，不作为负例。没有可用排序对时该项为可微零；自然背景空池不制造正参考。冻结基座用no_grad而非inference_mode产生可进入新头反传的特征，每步不允许基座梯度。

原自由池/首次截断观察可按window ID仅在本次拟合内复用，因为输入、基座、Budget和观察标签全部固定；不缓存新增头的hidden状态、logit或计算图，不改变抽样顺序与更新数，不另复制完整特征。

## 两种解码共享同一新增权重

early与late读取同一最终PrefixValue checkpoint及相同冻结B1基座；解码不读取参考。early只在事件beam和final固定链上限之前使用新增兼容性优先级；late生成原加性固定池，再使用相同新头作完整链读出。除排序位置外，端点、条件化exon、link、ORF及预算完全一致。

Budget保持原设置：1,024bp/类型8端点，每起点/carry/终点类型8条合法exon边，每donor8连接，每事件/carry2条beam，最终8链/kb，最多128exon；GT/GC–AG及物理最短intron5bp不变。proposal_scores仍为旧加性分，不能误作新头gain。

每个窗口将两池候选身份组成排序去重的union，使用一个全新的PrefixWindow一次计算完整链logit，再映射回两池。因此两臂共同候选的最终gain数值完全相同；early搜索中的no-grad cache不复用于最终共同读出。这样控制最终批次/cache数值差异，比较的改变保持为候选形成时的评分位置。

late自由链、顺序、retained links、counts、Budget必须逐窗重放原R5 B1已保存结果；旧proposal数值最大绝对误差≤1e-5。任一控制不匹配则配对对照未完成，不删除问题窗口、不放宽容差或把结果包装成公平对照。

owner按跨度中点归属最近窗口中心、平局较小start，去重取最大同臂gain，最大正增益非重叠selector及gain>0均不变。输出primary protein-coding CDS链，不声称全UTR或所有isoform。

## 新执行资源硬限

| 阶段 | 新增分配硬限 | 估算依据 |
|---|---|---|
| 一次完整拟合 | 1 RTX3090 2CPU 16GiB，3小时 | 12窗准备折算约1.87小时，单次加性控制生成；留启动与波动余量 |
| 全DEV配对推理 | 1 RTX3090 2CPU 16GiB，8小时 | 小样本未训练前移搜索折算约6.2小时，加共同读出与运行余量 |
| 合并及统一评价 | 无GPU，2CPU 32GiB，2小时 | 原全范围评价与两种新输出 |

本新机制的正式阶段最多11 allocated GPUh和4独立allocated CPUh，新增磁盘产物≤2GiB。GPU阶段I/O及CPU等待计入GPUh，不把小样本计时当作完成保证。两GPU阶段依赖串行，推理afterok拟合，最终评价afterok完整推理；不自动启动其他新实验。

准备已用0.027778 GPUh及0.012222 CPUh；R6另用0.040556 GPUh。与R4、R5、已关闭GRU优化分支的历史成本分别保留，不重新开始计费、隐去旧超时或声称本轮数值是整个研究总成本。

输出分别为outputs/M28-B2-FIRST-FIT、outputs/M28-B2-PAIRED-INFER及outputs/M28-B2-EVALUATION，开始前必须不存在。新目录独立，不改R5结果或任何冻结权重。按实时Slurm状态选择合格24GiB GPU分区，墙时不超过上述上限，环境generanno。

任何已启动阶段超时、OOM、非有限、冻结基座出现梯度、输入/标签错配、控制重放不匹配、参考进入自由推理或产物超限均停止该阶段并保留证据；不自动重提、续半成品、延时、改阈值/预算或用不完整DEV排名。未启动的提交拒绝可在命令与科学范围不变时修正调度错误，必须先排除重复提交。

## 全范围结果与下一步判断

只有完整4,608更新、两种解码均完成8,690个DEV方向窗且控制重放通过，才形成正式位置对照表。统一主参考分母7,728（拟南芥5,437、水稻2,291），旧6,450副表单列，不能换分母。

按物种、macro和micro-pooled报告exact-chain P/R/F1、CDS-base、自由候选覆盖、owner及selector逐阶段损失、背景span/CDS覆盖、全背景链/Mb、GFF3结构支持及真实分配成本。CDS-base是最终输出union的覆盖指标，不是原始局部head质量；参考不一致的背景负担不自动称生物学假阳性。

机制继续条件预先设为：early相对同新增头的late在两物种均改善自由完整链数、owner完整链数和最终exact-chain F1，同时两物种背景span覆盖均不增加。部分收益、单物种改善或注入训练池表现均不足以通过该条件；结果仍全部保留并解释，不扫K、阈值、seed或追加训练救结论。

通过条件也只说明该DEV上的位置机制值得继续，不等于优于C0、旧B1或SOTA。C0/B1及ANNEVO、Helixer、Tiberius均以同口径历史缓存列出，注明未重新运行；不宣称独立zero-shot或通用优势。期刊核心主张仍须完整方法收益、独立数据及额外冻结协议支持。
