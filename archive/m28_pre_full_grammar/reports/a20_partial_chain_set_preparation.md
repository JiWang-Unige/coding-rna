# A20 部分观测链集合似然准备

2026-10-05，在准备计算前固定。本项判断：在完整自然候选域中把 unknown 作为潜变量边缘化，是否数学上可计算，并与现有标签保护语义相容。它只使用已经核验的 A15 TRAIN 自然候选及保存 gain，不做训练、新模型前向、自然搜索或 DEV。A15 的原 FAIL 不被本项替代，A17 保持原协议。

## 原始方法先例

[Tsuboi 等的 COLING 2008 原论文第4节、式2–4](https://aclanthology.org/C08-1113.pdf)以与部分观测一致的所有结构概率之和作为似然；其对数是受约束与完整输出空间的两个 log partition 之差，梯度是两种分布的特征期望之差。原文也指出该目标通常不再凹。这里只借用部分观测的统计原则，不将其当作新颖损失首创，也不引用 NLP 实验作为 gene 模型性能证据。

以下链集合形式是本项目的推导，不是原论文中的基因方法。已知正链 G、已知非匹配 N、unknown U 保持原标签。完整分母覆盖全部自然候选的 span 非重叠集合及空集；受观测事件要求全部 G 出现、全部 N 不出现，U 自由。由于 G 必须彼此兼容，受约束 partition 等于 `sum(gain_G) + logZ(U_compatible_with_G)`，loss 为两个 log partition 的差除以原已知候选数。半开 span、相触兼容、原 gain 与 null 含义不改。没有已知候选时事件为全部空间，loss 和梯度严格为零；没有 unknown 时还原 A15 似然。

## 必须检验的语义风险

“不把 unknown 标成0”不等于“unknown 分数没有梯度”。若 unknown 与必选正链冲突，它不能进入受约束分子，却进入完整分母，梯度下降可能压低该 unknown。unknown 与已知负链或其他候选的结构关联也可能产生梯度。与已知部分完全分离的 unknown 组件才会从两个 partition 中抵消。

因此本准备不能把 latent 处理冒充原严格零梯度 ignore 政策。尤其对完整 protected alternative 链，这种间接监督需要明确新的可观测性定义，而不是默默改变既有保护。先用解析反例证明差异，再检查所有实际自然 unknown 中与 gold 冲突的数量及其中 exact protected 几何的数量；报告全部 saved-gain 导数，不挑效果好的窗口。发现实际 protected 冲突就不按“原标签保护不变”的理由开展首训，保留代码作为方法准备，继续设计标签可观测性或候选输出语义。

## 固定计算与决定

新增 `src/m28/partial_chain_set.py` 复用原 interval partition；不修改 A14／A15 runner、旧 loss 或结果。五项检查分别检测完整／受约束分子及梯度是否符合穷举；无 unknown 是否等价 A15；相触独立 unknown 是否抵消；与 gold 冲突的 unknown 是否确有非零梯度；全未知／空集和冲突 gold 的语义是否正确。解析检查在 Slurm CPU 内执行，不建通用训练框架。

随后 `prepare_partial_chain_set.py` 仅读取已完成的 13387977 核验记录、A15 末 fit 配置、1,491 个自然窗口、16,321 个保存候选与 gain，以及原六条 TRAIN 染色体的 catalog 元数据。逐项对齐几何／来源／标签；全部保存 gain 用 float64 重组两个目标及对 gain 的导数。该导数不是评分头参数梯度，也不是重新前向，不能据此推出实际训练效果。原1,536 draw及物种曝光另表累计，不把重复 draw 当独立链。

输出 `outputs/M28-A20-PARTIAL-CHAIN-SET-PREPARATION/{windows.jsonl,summary.json,STATUS,JOBID}`，原结果不覆盖。generanno、private-teodoro-gpu、0GPU／2CPU／8GiB／10分钟上限；预计数十秒至数分钟，新增文本数MiB，不加载 checkpoint、DNA 或 feature tensor。不放行 DEV／test／Setaria，不改任何运行作业、Git 或公开资料。此短项由主对话闭合，不增加监控；长 A17 链仍由同一固定监控等待。

即使算术和所有真实记录通过，也只完成新方法的准备验证，不是模型选择 PASS、准确率、背景、泛化或论文结果。发现当前保护语义不兼容时不拟合该定义，不修改标签规避；完整 raw-FASTA、全部参考、逐物种／worst-species 与 ANNEVO／Tiberius／Helixer 的原目标仍须推进。

## 实际准备结果与决定

作业 13389502 经 `sbatch --parsable sbatch/M28-A20-PARTIAL-CHAIN-SET-PREPARATION.sbatch` 提交，private／0GPU／2CPU／8GiB／10分钟上限。Slurm root／batch／extern 均 COMPLETED、ExitCode 0:0、Elapsed 33秒；STATUS=COMPLETED、JOBID=13389502，windows.jsonl 为 1,491 行，stderr 为空。五项解析检查全部通过，0.228秒；真实候选计算 12.306秒。六份新源码／测试／SBATCH／研究记录均在提交前与 canonical 远程逐字节核对，原模块和训练流程不改。

全部自然候选几何、来源和原标签与闭合 A15／A12 记录一致；保存 gain 重组和所有导数有限，没有 unknown 的窗口还原 A15 loss 与 gain 梯度。下表按唯一窗口／候选配对计数，不能当作独立基因数。

|物种|唯一窗|自然候选|unknown|与 gold 冲突的 unknown|其中 exact protected 链|有冲突窗口|
|---|---:|---:|---:|---:|---:|---:|
|拟南芥|745|8909|853|139|127|120|
|水稻|746|7412|522|35|21|23|
|合计|1491|16321|1375|174|148|143|

148 个 exact protected 链来自 134 个唯一窗口；它们是 catalog 中完整受保护的同几何链，不根据这个计数把它们全称为 alternative 或已确认真阳性。每遍原 1,536 draw 中相应计数为 1,398 unknown、176 个 gold 冲突、150 个 exact protected gold 冲突，涉及 145／136 个 draw；重复窗口按实际曝光保留。两个物种各 768 draw，原自然正／已知非匹配／unknown 曝光分别为 1,806／6,508／869 和 667／6,427／529。

原 A15 loss 对所有 unknown gain 的直接导数确为 0；新 partial loss 的 unknown gain 导数范围为 −0.1140572 至 0.9217435，既可增加也可压低 latent 候选分数。首个按原窗口顺序记录的实际 witness 是水稻 `NC_089043.1|6193152|-` 中一条六段 exact protected 链：原保存 gain=6.8741784，旧 loss 直接导数=0，新 loss 直接导数=0.00632208。这不是近零等价；梯度下降会压低该候选的独立 gain 变量。完整每窗 unknown 索引和导数保留于 windows.jsonl，没有按梯度或预期收益筛窗口。

这里“零梯度保护”只指损失对候选 gain 变量的直接导数，不表示共享评分头更新时 unknown 预测始终不变；本项没有计算参数梯度或进行训练。原 unknown 标签未改为0，但它进入结构竞争后确会受到观测约束的间接监督。这证明当前 partial 定义与原 ignore 语义不同，不证明边缘似然本身数学错误，也不证明它导致 A15 FAIL 或任何准确率变化。

据事前规则，本定义不进入拟合；不改标签或删 protected 链来获得可训练结果。五项软件／算术检查通过与标签语义不相容是不同结论。本项 optimizer_steps=0、模型／backbone forward=false、checkpoint_loaded=false、DNA／feature tensor=false、新搜索=false、DEV=false、test／Setaria=false。后续模型需要先定义哪些观测约束允许多个完整链或只约束可判定区域，而非声称 unknown 边缘化自动保持保护；该设计仍需新的明确协议与真实自然输出验证。

A17 仍沿原局部上下文协议训练；同一固定 heartbeat 已仅保留其三个 roots，配置逐字节读回和同监控对话人工检查已验证，当前新范围的自然定时运行及终态交接尚未验证。A15 已由主对话消费，不重复等待或交接；A20 已闭合，不加入监控。论文目标未完成，本次没有新的性能胜出、新颖性、泛化或发布主张。
