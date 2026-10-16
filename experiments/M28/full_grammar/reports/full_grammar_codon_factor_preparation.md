# 原 TRAIN 的完整 CDS 密码子评分接入

2026-10-12 UTC。最终头9216／codon4608、六raw及共同评价／配对已完成并由主核验；[完整负结果](full_grammar_codon_result.md)的pooled完整链F1为65.468852%，比B低4.833772 pp，两物种及六染色体均下降。工程接入完成不等于性能改善；保留原权重和全部结果，不增训、不改参数，D仍加入原冻结B/C/D DEV。以下保留先定科学设计与历史执行记录，不再代表仍待最终权重。

## 问题与评分变化

现有模型按十五角色逐碱基加分；实际codon carry用于判断结构是否合法，但没有给完整的三碱基密码子独立评分。本次在同一gene-set能量上增加完整拼接CDS的密码子项，每个完整密码子计一次，包括起始ATG和末端TAA／TAG／TGA。CDS拼接来自候选路径本身，不能由参考链、注释或预先选择的外显子生成。

跨内含子时保留已有的实际未完成编码碱基，读到第三个编码碱基才加分。内含子和background不计密码子分数；不允许把stop拆开跨内含子，也不放行内部stop、模糊编码碱基或未完成基因。单碱基中间外显子、三碱基首／末外显子的原角色优先级、非重叠gene-set、空集、cap128及全部原DNA合法性规则不变。

新增64个FP32原始参数，按ACGT字典序记录。每次计算先分别扣除61个sense和3个stop参数的均值，得到62个去均值方向；这不意味着全部方向与神经势函数均可独立辨识。零初始化使原模型嵌套于新模型。参数是学得的能量，不是从目标染色体注释估计的固定codon-usage先验；在训练分区、观察分区、raw解码和输出链评分中使用同一项。

这是性能候选，不把密码子评分本身当作创新性证据，也不把逐碱基内含子负margin解释成已定位的因果缺陷。方案在既有TRAIN及历史DEV诊断已知后确定，不称DEV前预注册，不通过DEV挑选此分支、训练时长、参数或checkpoint。

## 数学与独立实现

完整gene-set的能量为原role-minus-background逐碱基和，加上所有选中CDS拼接字符串的完整密码子评分和。新增转移索引只标记原状态表的codon闭合边，不新增状态或合法路径。分区保持所有可达状态及原一阶梯度；不使用MAP剪枝。解码沿用同stage、实际carry、较低已用exon预算的精确max-only dominance，因为这些状态后续的role和codon评分仍相同。

实现位于独立full_grammar_codon.py、full_grammar_codon_skyline.py及annotate_full_grammar_codon.py；原boundary、skyline、coarse-channel和两个活动分支的文件不修改。fine仍硬观察，coarse仍用固定epsilon0.1的三类通道，unknown仍边缘化；保留绝对background观察校正和原fine+coarse分母，不新增辅助loss或改变源标签。

## 已完成的人工 DNA 验证

13497345的root／batch／extern均COMPLETED、ExitCode0:0，17:01:55至17:02:33UTC，2CPU／4GiB／无GPU／5分钟上限；测试实耗22.394秒，batch MaxRSS704.50MiB。独立oracle枚举CDS拼接合法链及非重叠子集，不读取研究DNA、权重或参考，直接计算拼接密码子的能量。

七组验证覆盖完整及部分观察分区、MAP及输出链评分、FP32／FP64 logit和64参数梯度、零权重对原分区／logit梯度／输出链的严格复现、cap1／2／4／128不增加状态、跨两个内含子的拆分密码子、短外显子及三种stop、全部内部切点和checkpoint跨度、fine／coarse／unknown通道的绝对概率与epsilon梯度、去均值／逐行平移不变性，以及未完成基因和跨拼接内部stop的拒绝。人工零势函数例中，CTT密码子项可以改变MAP；这仅证明表达能力和实现，不是研究数据上的收益。

首轮13497336为FAILED／1:0：6组通过，全未知观察的参数梯度残差8.3266726847×10⁻¹⁷触发逐位零断言。修复只把该人工测试改为FP64绝对容差10⁻¹²，损失本身仍严格为零；模型代码、训练选择和评价规则不变。两轮输出完整保留，不覆盖失败记录。

## 原生预检和固定追加学习方案

预检13497400经afterok13497345提交，实际private-teodoro-gpu／单24GiB RTX3090／2CPU／16GiB／15分钟／Requeue0。沿用B预检的每物种第一个fresh draw，不根据新loss换样本；从原fit13434603的最终4608头与完整AdamW状态各恢复一次。新64参数从零开始，使用同一AdamW超参数；原头和新参数一起做norm1裁剪。每样本只做一次可丢弃更新，原头step4608→4609，新参数step0→1；不保存checkpoint。

预检先核对真实native零参数NLL与原模型完全相同、role VJP一致、新参数导数有限。随后把同一更新后的头和64参数装入严格raw adapter，核对实际joint24576训练势与raw势、参数、完整context MAP和链评分，并运行14bp短原生raw caller的两个方向。final9216 metadata仅为未保存的路由fixture，不能当作已训练9216权重。checkpoint身份必须拒绝错codon顺序、去均值规则、采样seed、编码跨度或参数形状。encoder冻结，无缓存特征、DEV或sealed读取。

已启动的候选固定从同一原4608起点做fresh seed1／4608有序draw／4219 unique／三段各1536的追加学习，沿用B的joint24576、原有限view和ownership、所有十五角色、cap128、FP64分区、epsilon0.1、importance1、原分母、lr0.0003、weight_decay0.01、betas0.9／0.999、optimizer eps10⁻⁸及clip1。原头AdamW年龄从4608至9216，新参数从0至4608；不能把它们写成同一moment年龄。相对于B的唯一机制扩展是零初始化的密码子能量及其64个参数，联合裁剪和随后的优化轨迹是这一干预的组成部分。

唯一正式读出是最后累计9216／epoch6；6144和7680只保留训练快照，不依据TRAIN或DEV筛选。完成权重、draw、观察计数与AdamW核验后，在原六条TRAIN上严格raw-FASTA推理，参考只用于评价，仍保留全部22688参考（含不受支持者）。与B的相同raw输入及评价口径配对，保留逐染色体、逐物种、pooled、macro和worst的完整链及CDS-base计数与P／R／F1；不得用base收益替代完整链结果，负结果和物种差异同样保留。单个同起点实验检验这一追加阶段的表现，不证明泛化、独立重复或论文达标。

现有B+C及新增D流程均接入唯一原固定心跳；不新建监控对话或重复自动化。封存数据、新物种／新面板、DEV反馈调参和对外发布仍不在本阶段范围内。论文goal尚未完成。

## 原生预检实际结果

13497400的root／batch／extern均COMPLETED、ExitCode0:0，17:05:45至17:08:47UTC；脚本实耗150.4185164秒，batch MaxRSS4533448KiB，stderr为空。两条固定实际view仍是B预检的additional step1／2：Arabidopsis NC_003070.9的14389248负链窗口及Rice NC_089046.1的14462976负链窗口，各有49152bp真实context，原已知观察分母均24576。

| 固定 TRAIN view | 零参数 raw NLL | role VJP 与原模型最大绝对差 | 未归一化 codon VJP L2 | 更新后去均值 codon L2 |
| --- | ---: | ---: | ---: | ---: |
| Arabidopsis step1 | 0.00149551027688 | 0 | 0.00158301272388 | 0.000452692112717 |
| Rice step2 | 30.7035429458 | 0 | 26.3764731152 | 0.00235572274797 |

主对话另按输出值复算两条原分母下的weighted NLL（6.085246894855437×10⁻⁸及0.0012493303607496638），核对固定样本、实际context、原头／新参数的不同optimizer年龄、非零有限参数导数与更新，以及所有零权重和raw一致性标记。两次独立可丢弃更新均为原头4608→4609、新参数0→1；源权重和moment未修改。

更新后训练native势与严格raw适配器势逐元素相同，64个原参数相同，真实context的完整MAP链和链评分相同。各次14bp原生raw caller均保留全部DNA行与两个方向；均未预测CDS，因此短输入只证明有限执行和行／方向接入，不证明非空基因的坐标恢复或准确率。人工DNA枚举和真实native接入是不同层次的证据，均不等同于新的正式final9216模型或原六条TRAIN性能。

## 严格 raw 接入与评价端点检查

新增infer_full_grammar_codon_train.py核对最终头9216／密码子4608、严格checkpoint加载、全部64参数及去均值读出，并要求新密码子backend和完整拼接CDS能量身份。原评价入口只接受旧role评分身份，因此本分支独立接入身份核验，不修改旧入口或B／C流程；参考选择、坐标与phase核对、计数评价和聚合仍使用原共同实现。

六条完整TRAIN染色体只从原FASTA构造推理输入，保留全部碱基与两个方向，不给caller参考链或观察模型。评价读取参考前，逐碱基直接比较新分支、B和原DNA经相同大小写／非ACGT转N规范化后的输入。compare_full_grammar_codon_train.py再核对两臂参考总体、支持分层、逐染色体及物种／pooled计数和P／R／F1；macro分别按两物种计算，worst分别取各臂最差物种后作差，不能用最小逐物种差替代worst差。

端点检查13497803的root／batch／extern均COMPLETED、ExitCode0:0，17:19:29至17:19:45UTC，2CPU／4GiB／无GPU／5分钟上限；四组人工fixture检查耗时0.007秒，batch MaxRSS805992KiB。覆盖旧backend／错误checkpoint与完整覆盖身份的拒绝、64参数及密码子顺序错配的拒绝、count-derived P／R／F1与macro／worst计算，以及参考总体、聚合计数或P／R／F1不一致的拒绝。fixture不读取研究DNA、参考或模型权重；通过不等于真实全染色体评价已经完成。

## 正式作业链与等待状态

| 阶段 | 实际作业 | 必需成功依赖 | 资源与时限 |
| --- | --- | --- | --- |
| 追加fit | 13497914 | 已完成并核验的预检13497400，由原参数传入启动核对 | 单24GiB RTX3090，2CPU／16GiB，18小时 |
| 权重与记录核验 | 13497943 | afterok13497914 | 无GPU，2CPU／4GiB，5分钟 |
| 六TRAIN严格raw | 13497945_0–5 | afterok13497943，最多3项并行 | 每项单RTX3090，2CPU／16GiB，2小时 |
| 全范围评价 | 13497946 | 整个raw数组afterok13497945且B评价afterok13494898 | 无GPU，2CPU／4GiB，20分钟 |
| 与B配对 | 13497947 | afterok13497946 | 无GPU，2CPU／4GiB，10分钟 |

五阶段实际提交时间为17:23:06至17:23:38UTC。Slurm控制器显示fit在gpu034实际RUNNING且Requeue0；其余阶段为有效Dependency，评价明确依赖整个raw数组而非单项。实际started.json确认预检13497400、源fit13434603、support13453408及固定配置。最终唯一读出路径为outputs/M28-FULL-GRAMMAR-CODON-FIT/step_009216.pt，输出前缀分别为CODON-FIT、CODON-VERIFICATION、CODON-TRAIN-RAW、CODON-TRAIN-DIAGNOSTIC和CODON-PAIRED；尚未分配的raw numeric roots必须以后从实际Slurm任务映射取得，不能按数组ID猜测。

首次提交要求afterok13497400及13497803时，被Slurm以Job dependency problem拒绝，没有创建作业。取证确认13497400仍为sacct COMPLETED0:0，但已超过MinJobAge300秒并从活动控制器清除；只移除这两个已完成准备作业的提交依赖，保留相同预检作业参数、启动身份核对和科学命令后提交成功。新链各阶段的实时afterok依赖完整保留，未重跑预检或更改训练方案。

唯一心跳coding-rna-r5-light-monitor已通过应用工具扩展到B+C+D，仍绑定原监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3，每15分钟、gpt-5.6-luna／low、正常等待安静。实际保存核验通过：原目标、心跳类型、间隔和通知策略不变，原B+C白名单逐字保留；D仅增加调度、身份和必需产物的结构核对。新增D范围的真实定时追加现已观察：原监控turn01a12c0f-59cf-7352-9d27-e05a6a43fb19于17:42:46.018UTC启动，heartbeat输入包含同一自动化和B+C+D范围，实际model／effort仍为gpt-5.6-luna／low；17:42:58UTC检查三个fit均RUNNING、下游有效Dependency，连接失败0、未读取科研内容或交接。配置与新范围定时运行分别已验证，D整链终态回传仍待实际发生。整条流程完成或真实actionable问题才向原主对话交接一次；不把fit或单项raw完成当作全流程完成。等待中的作业与已通过工程检查不构成准确率、泛化或论文成功结论。
