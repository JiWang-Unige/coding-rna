# 原 TRAIN 三种目标续训对照

2026-10-09，用户解除“不改训练参数、不增训”的阶段限制，授权在原 TRAIN 内自主实施并行实验。本轮检验继续结构似然优化与改变目标家族能否改善完整 CDS 链恢复；不以训练损失下降代替模型竞争力或新颖性证据。

## 固定比较

三个分支均从原作业 13434603 的最终 step_004608.pt 独立恢复 C0 权重和完整 AdamW 状态，各增加 3 遍原顺序 TRAIN，即 4608 步、累计 9216 步。不采用中间 checkpoint 选模。

| 分支 | 原始观察 NLL |
| --- | --- |
| STRUCTURED_CONTINUE | 原完整基因集合结构观察 NLL |
| LOCAL_CHANNEL | 各核心碱基独立的同通道观察 NLL |
| HALF_STRUCTURED_HALF_LOCAL | 0.5 结构 NLL 加 0.5 局部 NLL |

局部项以 15 角色 softmax 定义碱基概率：fine 行使用真实角色概率；fine 缺失且 coarse 已知时，对真实区域内角色赋 0.9 观察概率，区域外角色赋 0.05，再对潜在角色求和；未知行观察概率为 1。fine 优先于 coarse，区域角色集合复用原 REGION_COLUMNS。原 importance / max(1, fine+coarse 行数) 乘数不变，无类权重、边界权重、额外长度或 null 权重。

保持原两物种六 TRAIN 染色体、1536 draws / 1491 unique、seed 0、24576 bp 核心及左右各一个 12288 bp 步长邻窗、固定 owner 拼接。GENERanno 编码器及其 revision 冻结，每次从原 TRAIN DNA 进行原生 BF16 编码并在推理模式外克隆 FP32 特征；只训练原 366191 参数 C0。优化器维持 AdamW lr 3e-4、weight decay 0.01、betas 0.9/0.999、eps 1e-8、clip norm 1。结构项继续使用 cap 128、checkpoint 256 bp、epsilon 0.1、单位角色权重及 FP64 累加。三个目标梯度尺度不同，因此不能声称有效更新幅度已匹配。

## 实施验证与执行

先用人工序列独立枚举局部观察概率，核对局部和混合损失值及一阶导数、fine 优先级、未知行零梯度和行平移不变性。再对各物种原顺序第一个 draw，在每个分支独立恢复原状态后做一个一次性更新；验证原生输入、有限损失和梯度以及 AdamW 从 4608 到 4609。六个一次性更新不保存科学 checkpoint，不进入正式续训，也不用于分支选择。

验证通过后并行提交三个独立续训作业，私有 RTX3090 每分支一张 GPU、2 CPU、16 GiB RAM；结构和混合分支设 18 小时上限，局部分支设 3 小时上限。每遍 checkpoint 仅用于保留工作，正式比较只使用累计 9216 步末权重。原 checkpoint、原训练与冻结结果不覆盖。不自动恢复或取消作业；真实实现失败按已授权范围克制修复，科学负结果如实记录，不为改善数字改变规则。

实际数学检查 13451059 与原生预检 13451086 的 Slurm root/batch/extern 均 COMPLETED 0:0；数学检查 3 项通过，原生预检 6 次独立一次性更新通过且未保存 checkpoint。预定 AT 样本局部原始 NLL 48407.26，结构 NLL 22.87；Rice 样本分别为 48429.37、0.00166。局部项初始尺度明显更大，等权系数不等于等梯度幅度，不事后调混合系数。去除初次编译的结构单样本约 11.3 秒、局部约 1.36 秒；结合原正式末遍 17323 秒 / 1536 draws，为上述 walltime 提供粗略依据，不将两样本推算包装为完整队列实测速率。

2026-10-09 07:50:23 UTC 三个分支开始运行，均为 private-teodoro-gpu / gpu034，1 RTX3090、2 CPU、16 GiB。此时只是正式训练启动，不是训练完成或性能结果。

| 分支 | 实际 Slurm task | walltime |
| --- | --- | --- |
| STRUCTURED_CONTINUE | 13451136_0 | 18 小时 |
| LOCAL_CHANNEL | 13451138_1 | 3 小时 |
| HALF_STRUCTURED_HALF_LOCAL | 13451137_2 | 18 小时 |

## 后续结果判断和边界

每个末权重先完成严格 raw-FASTA 的原六 TRAIN 整染色体双链推理，caller 不读取参考、标签或 feature cache，复用原 exact exon-budget skyline、owner 策略、unit role gain、null 0 和 cap 128。预测完成后用原共同评价尺读参考，报告各物种与 pooled exact-chain precision/recall/F1、CDS base、背景假阳性和重复记录，保留不可支持参考层。必要的有限视图角色混淆分析属于训练机制诊断，不是新参考选择 caller。

新 raw 入口与评价路径绑定检查 13451240 的 Slurm root/batch/extern 均 COMPLETED 0:0，三项人工检查通过；这次检查未运行原生编码器或研究 DNA。旧入口仍严格限定 final4608，新入口仅接受对应分支累计 final9216；新输出路径显式传入原评价函数，原默认输入、参考选择、指标及过滤规则不变。

各分支正式 fit 完成后依赖 CPU 核验全部 4608 条新更新、原 draw/标签数/权重顺序、目标项与归一化算术以及三遍 checkpoint / AdamW 状态。通过后才启动六条 TRAIN raw 推理，最后整个数组 afterok 启动 CPU 评价。每个 raw 数组最多并发 3 个任务，沿用每任务 1 RTX3090、2 CPU、16 GiB、12 小时；提交时末权重推理和评价未完成，当前部分结果见末节。

下列 afterok 链已提交，推理与评价在依赖完成前不运行，不将 PENDING / Dependency 视为失败。

| 分支 | fit | CPU 核验 | 六任务 raw 数组 | CPU 评价 |
| --- | --- | --- | --- | --- |
| STRUCTURED_CONTINUE | 13451136_0 | 13451244 | 13451470 | 13451625 |
| LOCAL_CHANNEL | 13451138_1 | 13451243 | 13451469 | 13451623 |
| HALF_STRUCTURED_HALF_LOCAL | 13451137_2 | 13451245 | 13451471 | 13451624 |

输出分别置于 outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST、OBJECTIVE-CONTRAST-VERIFICATION、OBJECTIVE-CONTRAST-TRAIN-RAW、OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC 的对应分支目录；后三者同样使用 M28-FULL-GRAMMAR 前缀。未覆盖原 final4608 或其冻结预测和评价。

原 final4608 是不变锚点。若只降低各自 NLL、完整链却不改善，按负结果处理；若训练改善，也不能由 TRAIN 宣称泛化或超越现有模型。冻结 DEV 的口径、数据隔离和既有基线不改，不将 DEV 用于训练目标权重、学习率、checkpoint 或阈值选择；后续候选进入 DEV 前记录确定的末权重与方案。封存 / held-out / Setaria 不放行。未授权公开发布、Git 提交或推送，也不恢复旧自动研究框架。

长等待仅复用监控对话 01a0f71f-8845-7432-bd13-bf54c1808fb3 与心跳 coding-rna-r5-light-monitor，gpt-5.6-luna / low、15 分钟；正常运行静默，终态一次交接回原主对话。

2026-10-09 08:03:17 UTC 已更新并启用原心跳；保存类型为 heartbeat、目标仍为原监控对话。08:03:58 UTC 在原对话完成主触发独立首检：三 fit RUNNING、核验 / raw / 评价均 PENDING 有效依赖，未读指标或交接。实际首检 turn 的模型仍为 gpt-5.6-luna / low。这证明配置和一次独立检查成功；本轮真正定时追加及最终自动回传尚未观察。定时检查需要本机开机、应用运行和 SSH 可用，不能保证休眠或断线时即时执行。[定时任务官方说明](https://learn.chatgpt.com/docs/automations?surface=app)。

## LOCAL_CHANNEL 已完成的部分结果

2026-10-09 UTC11:42主对话核对该分支fit13451138_1、CPU核验13451243、全部六raw任务13451469及评价13451623的root/batch/extern均COMPLETED 0:0。fit实际07:50:23–09:33:41，核验09:33:42–09:34:26，六raw最后11:16:41结束，评价11:16:42–11:17:43。核验确认全部新增4608日志原顺序、观测及importance一致、局部目标算术一致，三checkpoint严格加载、权重与AdamW矩有限，最终累计9216。raw及评价沿用原参考和双链整染色体尺，未读DEV或sealed。STRUCTURED_CONTINUE和混合分支尚未结束，整个三分支cohort仍不算完成，也未据此变更冻结候选或训练方案。

与既有final4608评价13443329比较，AT17,739／Rice4,949／合计22,688参考分母不变。表中的末9216是额外局部目标训练后的净变化，不是已匹配更新次数的纯目标因果对照；后者须等STRUCTURED_CONTINUE结束。

| 原TRAIN范围 | final4608链F1 | LOCAL_CHANNEL9216链F1 | 完整链TP／预测／参考 | CDS碱基F1：4608→9216 | 纯背景链数：4608→9216 |
| --- | ---: | ---: | --- | ---: | ---: |
| AT | 50.4871% | 44.1602% | 7821／17682／17739 | 93.0848%→95.8184% | 15→97 |
| Rice／worst species | 28.2780% | 33.3818% | 1722／5368／4949 | 80.8511%→88.4989% | 147→394 |
| pooled | 45.5901% | 41.7290% | 9543／23050／22688 | —→94.0966% | 162→491 |

9216各染色体链F1依次为NC_003070.9 42.9338%、NC_003071.7 43.5955%、NC_003076.8 45.9320%、NC_089043.1 36.0165%、NC_089044.1 33.8489%、NC_089046.1 30.4140%。该分支预测重复链记录均为0。背景span FPR在AT由0.3989%升至0.4871%、Rice由0.7322%升至0.9420%；未用背景阈值重解码。CDS碱基是unstranded union，与完整有向CDS链不同，不将前者提升称为完整gene恢复改善。

当前判断：局部目标续训提高Rice／worst-species链F1及两物种CDS碱基覆盖，但AT链F1下降6.3268个百分点、pooled下降3.8611个百分点且背景链增加，不能称为整体成功。模型family胜负及同更新次数目标影响尚未判定；不以这个部分结果改混合系数、选择checkpoint或退掉其余分支，不由TRAIN宣称竞争力、泛化或论文门槛通过。

来源：`outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC/LOCAL_CHANNEL/summary.json`、该分支VERIFICATION及原final4608诊断。固定监控只保留整个未完成三目标cohort的调度／结构检查，科学部分由主读取；已经读过的LOCAL部分不重复解释或触发新的终态交接。

## 局部目标续训后的完整链错误变化

2026-10-09，对已经完成的LOCAL_CHANNEL末9216和原final4608作逐参考配对。全部22,688条固定主参考保留，没有按命中、语法支持或观测覆盖筛掉记录。新入口`diagnose_full_grammar_local_chain_transitions.py`复用原错误分类和既有参考／观测记录，仅解析保存的预测，不读取DNA、权重或特征，不调用编码器、不重解码或更新模型。作业13453930实际UTC12:40:48–12:41:26、38秒，private-teodoro-gpu／gpu034／2CPU／4GiB／无GPU，root／batch／extern均COMPLETED 0:0；四种配对结果、分类方向及染色体／链方向完整join的人工检查通过。六染色体、两物种及pooled计数均与原共同尺的已完成评价一致。

下表的“丢失／新恢复”指同一主参考完整CDS几何链在两个末权重间的状态变化，不是不同训练样本之间的比较。

| 原TRAIN范围 | 参考 | 正确链保留 | 原正确链丢失 | 新恢复 | 两者均未恢复 | 正确链净变化 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| AT | 17739 | 6347 | 2256 | 1474 | 7662 | −782 |
| Rice | 4949 | 1077 | 286 | 645 | 2941 | +359 |
| pooled | 22688 | 7424 | 2542 | 2119 | 10603 | −423 |

AT多CDS段主参考的正确链由5837变4739，丢失2176、新恢复1078、净−1098；单段由2766变3082，丢失80、新恢复396、净+316。Rice多段由690变825、净+135，单段由673变897、净+224。AT三染色体正确链净变化依次−375／−152／−255，Rice依次+150／+88／+121，没有用物种平均掩盖一个相反方向的染色体。

AT丢失的2256条中，1816条保留主参考CDS两端但内部结构不同，占80.5%；351条同时改变两端及内部结构，89条仅CDS起端不同。Rice丢失的286条中，217条两端不变但内部不同、54条两端及内部均不同、14条仅起端不同、1条无同链CDS重叠。这里的同链最大CDS重叠配对是描述性关联，不是生物实例的一对一匹配或模块因果定位。

| 唯一有向几何事件F1 | AT 4608→LOCAL9216 | Rice 4608→LOCAL9216 |
| --- | ---: | ---: |
| CDS起端 | 79.4107%→82.3616% | 65.2282%→71.5324% |
| CDS止端 | 78.5196%→83.4981% | 58.4647%→73.4322% |
| CDS段边界对 | 77.7122%→70.9395% | 58.2060%→59.6915% |
| 剪接间隔边界对 | 80.5294%→73.9895% | 61.4471%→61.0688% |

AT剪接间隔预测数由70,782增至82,653，而命中由58,325变57,980，precision由82.4009%变70.1487%。CDS两端与碱基覆盖变好不能代替内部链恢复。事件表按染色体和链方向去重，分母不是主参考gene数，也不是phase评价或完整gene准确率。Rice间隔召回上升但precision下降，因此其完整链增益也不代表所有结构事件同步改善。

原观测覆盖分层同样保留。AT全参考跨度fine匹配且无fine冲突组1465条，正确链净−45；无直接core跨度观测组15196条，净−687。Rice对应649条／净+48及3983条／净+298。pooled全跨度fine组2114条保留842、丢失231、新恢复234、仍未恢复807，净仅+3；另209条不能作单链角色join的主参考保留在分母，均仍未命中。无直接core观测不等于未训练、无encoder曝光或独立测试；这些分层不作因果归因或确认性统计。

当前结果说明LOCAL续训后的AT净损失集中在多段内部链，而不是只差CDS端点。两次训练复用同一批draw，未新增直接fine覆盖；全跨度fine组也发生链丢失与恢复，无直接core组不能视为独立泛化。原锚点少4608次更新，尚不能把变化单独归因于目标家族；保持预先声明的STRUCTURED_CONTINUE与混合分支，等待其末权重和全染色体结果后作匹配更新数比较，不提前调混合系数、训练参数或冻结DEV候选。内部变化进一步按相对主参考的间隔新增／缺失分解，不能直接称为生物学假剪接或漏剪接。

完整结果为`outputs/M28-FULL-GRAMMAR-LOCAL-CHAIN-TRANSITIONS/summary.json`和`reference_transitions.jsonl`。新输出独立保存，原错误记录、原预测和运行中fit源代码未改；DEV与sealed均未读取。本轮短诊断由主对话处理，不加入原固定长流程心跳，也不新建监控。

## 两端不变的内部剪接间隔变化

随后用`diagnose_full_grammar_local_intron_edits.py`仅对上述两端不变的丢失／新恢复参考作精确间隔集合差。作业13453974实际UTC12:47:19–12:47:47、28秒，private-teodoro-gpu／gpu034／2CPU／4GiB／无GPU，root／batch／extern均COMPLETED 0:0。新增、缺失、同时变化及无变化拒绝的人工用例通过；每物种与pooled子组分母和13453930完全一致。只读已保存主参考、两组GFF3及配对记录，没有DNA、权重、特征或原生forward。

“额外”指预测中的间隔边界对不在固定主参考中，“缺失”指主参考边界对不在预测中；间隔差分本身不判断其他assessable isoform归属，也不称为生物学假／真剪接。下表仅列原正确链丢失且两端不变的子组，比较LOCAL9216预测与主参考。

| 原TRAIN范围 | 丢失链子组 | 同时缺失及额外 | 仅额外 | 仅缺失 |
| --- | ---: | ---: | ---: | ---: |
| AT | 1816 | 1778 | 37 | 1 |
| Rice | 217 | 213 | 4 | 0 |
| pooled | 2033 | 1991 | 41 | 1 |

AT子组有主参考间隔10,802对、预测12,942对，共有8350对、缺失2452对、额外4592对；Rice对应953／1239／675／278／564。这里对每个参考—预测配对求和，允许同一几何间隔在不同参考中重复计数，不能与上节全染色体唯一事件表互换。混合变化不说明具体间隔之间一一替换、边界位移或某个模块的因果失败；单纯删除若干额外间隔并不足以保证恢复大多数主参考，因为原间隔也已缺失。

新恢复且原4608保持同一CDS两端的674条中，原错误预测668条同时缺失及额外、3条仅缺失、3条仅额外；AT为488／2／1，Rice180／1／2。两个方向都主要是混合变化，不把它包装成LOCAL独有的错误机制。后续匹配更新数结果应同时考察主参考间隔的恢复与额外间隔的抑制，不能仅按间隔数量、CDS覆盖或端点改善判定结构恢复；当前不据此更改运行中的训练或冻结DEV。

结果为`outputs/M28-FULL-GRAMMAR-LOCAL-INTRON-EDITS/summary.json`及`reference_intron_edits.jsonl`，与首次配对结果各自保留，不覆盖或重跑已完成诊断。两个短流程均已由主消费，原固定监控对象仍仅为已提交的四条TRAIN长cohort。

## 主链变化中的同基因其他可评价 CDS isoform

在已知注释范围内，其他可评价isoform不能解释大部分LOCAL主链丢失。全部2542条丢失主参考中，LOCAL9216输出仅对42条发出同一参考gene的其他可评价CDS几何链；2500条没有这种输出。这个判断不要求它是最大CDS重叠预测，也不把其他isoform加算为主链TP。

作业13454143通过`sbatch --parsable sbatch/M28-FULL-GRAMMAR-LOCAL-ISOFORM-TRANSITIONS.sbatch`提交，用`diagnose_full_grammar_local_isoform_transitions.py`执行，沿用generanno环境及5分钟walltime。实际UTC12:58:39–12:59:34、55秒，private-teodoro-gpu／gpu034／2CPU／4GiB／无GPU，root／batch／extern均COMPLETED 0:0。同gene、同两端／不同两端、染色体／链方向及主链排除的人工检查通过。复用原共同尺，在原六TRAIN的已保存input FASTA和参考GFF3上重建可评价CDS：全部22688条主参考的ID、gene ID和CDS几何与保存记录一致；38152条可评价isoform记录的逐染色体数量与既有评价一致。原4608与LOCAL9216的主链计数和全输出其他可评价isoform计数都与已完成共同评价一致，没有改变参考选择、分母或过滤。

| 原TRAIN范围 | 丢失主链参考 | LOCAL输出含同gene其他可评价CDS链 | 占丢失参考 | 没有这种输出 |
| --- | ---: | ---: | ---: | ---: |
| AT | 2256 | 38 | 1.68% | 2218 |
| Rice | 286 | 4 | 1.40% | 282 |
| pooled | 2542 | 42 | 1.65% | 2500 |

两端相同但内部改变的丢失子组中，AT仅18／1816条、Rice4／217条发出同gene且同两端的其他可评价isoform；合计22／2033，剩余2011条没有同gene其他可评价链输出。新恢复的1474条AT主链中，原4608输出已有同gene其他可评价链者31条；Rice为2／645条。两个方向均没有把“主链错误”普遍变成已注释其他isoform命中。

全输出原共同尺的非primary其他可评价链TP在AT为4608的718条、LOCAL9216的614条；Rice为30→42条。它们与上述逐参考计数不是同一分母，不能把全输出isoform TP直接加到主链丢失参考上。完整CDS几何匹配不证明表达、UTR或所有isoform恢复，未匹配当前可评价注释也不证明结构在生物学上虚假。

这项后验解释使用TRAIN参考DNA和注释，但没有向caller提供参考、重跑推理、加载checkpoint／特征、执行原生forward或更新模型；DEV与sealed均未读。结果独立保留在`outputs/M28-FULL-GRAMMAR-LOCAL-ISOFORM-TRANSITIONS/summary.json`及`reference_isoform_membership.jsonl`，不改旧诊断。当前重点仍是匹配更新数下的内部主链恢复与间隔选择，不能将主要链损失归因为已知isoform选择口径，也不据此改正在运行的实验、冻结DEV或训练候选。

本短流程已由主对话核实并消费，不加入原固定长流程心跳，不创建监控对话或自动化。五条已启动训练分支的末权重和后续完整raw结果仍待原固定监控交接；全目标未完成，当前结果不是竞争力、独立泛化或论文成功证据。

## 固定末权重下的目标梯度诊断

混合目标的两个系数均为0.5，但这不保证参数梯度等幅。已有原生预检13451086中，各物种第一个原draw在独立还原final4608后，AT结构／局部／混合梯度范数为0.257126／9.953336／4.987146，Rice为0.0000147689／13.783139／6.891576。这些是采用原importance与known分母后、clip前的参数梯度，不以原始NLL比例代替梯度比较，也不以两例代表整个TRAIN。

2026-10-09，在读取新结果之前固定诊断样本为原draw顺序中每物种前12个、共24个，保留原观测和原context。新入口`diagnose_full_grammar_objective_gradients.py`仅加载不变final4608，在每draw同一组原生potential上分别求结构项和局部项的C0参数梯度；报告归一化后两项范数、局部／结构范数比、两项夹角余弦及混合向量与各项的余弦。混合向量严格为0.5倍两项之和。各物种首draw的三个范数和loss须与既有独立原生预检一致，否则不解释新结果。对每项FP32参数梯度采用FP64向量统计；零范数对应的余弦或比值保留未定义，不强制造出方向。

这项诊断不构造optimizer、不执行更新、不写checkpoint，不改变当前五条fit的输入、源代码、系数、抽样或阈值；只读原两物种六TRAIN及其观测，native编码器保持冻结，没有DEV或sealed。它是原始梯度的固定权重切片，不是已含AdamW旧矩、weight decay及训练轨迹的实际更新贡献，也不是同更新数完整链比较。样本按原顺序选取而非随机总体估计，不据此事后调混合系数或选择候选。

沿用原生预检的generanno环境、1 RTX3090／2CPU／16GiB，24个draw按原结构单样本约11.3秒粗估约4.5分钟，再为原生加载、首次编译及两次参数反传留余量，walltime20分钟。通过`sbatch --parsable sbatch/M28-FULL-GRAMMAR-OBJECTIVE-GRADIENTS.sbatch`提交作业13454272，实际UTC13:16:10–13:21:21、5分11秒、private-teodoro-gpu／gpu034，root／batch／extern均COMPLETED 0:0。输出独立置于`outputs/M28-FULL-GRAMMAR-OBJECTIVE-GRADIENTS`，只保留`gradients.jsonl`逐draw统计和`summary.json`，不保存特征或梯度向量。

固定24个ordinal为0–22及24，各物种12个。首draw原生预检的六组loss／梯度范数全部匹配，所有参数梯度有限，无零范数项，C0权重与载入状态逐tensor相等，native编码器冻结且无gradient；未构造optimizer、未更新或保存checkpoint，DEV与sealed未读。

| 固定 TRAIN 抽样 | 局部项梯度更大 | 局部／结构范数比 min／median／max | 两项余弦小于0 | 混合与局部余弦 min／median |
| --- | ---: | --- | ---: | --- |
| AT | 12／12 | 8.2142／148.0677／837873.9427 | 4／12 | 0.996573／0.999996 |
| Rice | 12／12 | 2.2356／3152.7167／1130385.1051 | 4／12 | 0.931809／1.000000 |
| pooled | 24／24 | 2.2356／267.4525／1130385.1051 | 8／24 | 0.931809／0.999999 |

表中Rice的median按六位小数显示为1，实际0.9999995336331255，并非精确同方向。两项余弦范围AT为−0.941208至0.952763，Rice−0.903522至0.943027；“反向”仅指向量夹角大于90度，不等于所有参数分别反号，也不证明一种目标错误。

这项固定checkpoint切片显示，等系数混合的原始梯度方向明显更接近局部项；它不能称为结构与局部学习贡献平衡的消融。比例、方向及反向例子都保留，不能仅报告大尺度比值或把混合与LOCAL视为等价模型；AdamW旧矩、clip及后续权重变化仍会影响实际更新。也不能将当前AT内部链损失直接归因于这些梯度。原预声明混合分支、匹配更新数STRUCTURED_CONTINUE及其冻结候选保持不变，等全部末权重和raw完整链评价闭合后判断，不据本诊断救援调参或挑checkpoint。

本短流程已由主对话核实并消费；没有加入长流程心跳、改自动化或新增监控对话。接下来复用原固定监控等待已提交五条fit及其核验／完整TRAIN推理／评价，不为未变化的目标继续叠加梯度检查；只有新结果或具体失败再推进下一项研究判断。

## 七个固定末权重的原 TRAIN 汇总

各分支独立评价完成后，用`aggregate_full_grammar_continuation_train.py`一次汇总全部七个预定末权重：原final4608、STRUCTURED_CONTINUE／LOCAL_CHANNEL／HALF_STRUCTURED_HALF_LOCAL末9216，以及双native CDS先验、新抽样、联合24576编码末9216。不按部分结果、F1或训练loss筛掉分支，不新增训练或重跑推理。主控制为同累计9216步的STRUCTURED_CONTINUE；final4608仅是不变锚点，少4608次更新，不能当作纯目标因果对照。

汇总只读原共同尺已经完成的七份measurement summary，沿用原`aggregate`和P/R/F1公式。所有六TRAIN染色体、17739条AT／4949条Rice／22688条pooled主参考都保留；逐染色体核对checkpoint及实际评价job绑定、原参考分母、CDS base union参考分母、assessable isoform记录数、原背景分母与定义和unsupported strata元数据，并核对每层计数、micro、macro及worst-species算术。这里只比较保存的参考元数据，不重读DNA／GFF、重建参考ID或几何，也不把相同分母自动称为新的参考逐记录验证。

完整JSON保留逐染色体及逐物种的exact-chain、unstranded CDS base、背景负担、非primary isoform、重复预测和phase检查数量；摘要列出各方法AT／Rice／pooled／macro／worst F1以及纯背景链数。phase检查数量不是phase准确率，CDS base不是完整链性能，raw计时不是公平速度基准。相对STRUCTURED_CONTINUE的差值统一用百分点，并明确累计步数是否匹配；双先验额外28参数、新抽样提议／曝光／梯度方差及联合编码输入差异继续保留，不因步数匹配宣称全因素匹配。

新入口和CPU脚本独立保存，没有改运行中的fit、caller或共同评价尺。正式输出为`outputs/M28-FULL-GRAMMAR-CONTINUATION-TRAIN-FAMILY/summary.json`和`comparison.md`，须等七份原评价全部完成后由主对话运行；当前不提交这个汇总作业，也不改固定监控对象。人工检查只使用构造的summary，覆盖缺分支／染色体、checkpoint／job错配、参考／背景／strata差异、各层算术与数据边界，不读任何研究数据或权重。准备完成不等于七方法科学结果闭合，也不放行DEV或sealed。

通过`sbatch --parsable sbatch/M28-FULL-GRAMMAR-CONTINUATION-TRAIN-CHECKS.sbatch`提交检查13454458，实际UTC13:31:56–13:32:22、26秒，private-teodoro-gpu／gpu034／2CPU／4GiB／无GPU，generanno环境、5分钟walltime；root／batch／extern均COMPLETED 0:0，六项人工检查全部通过、零failure／error。没有研究数据、权重、原生forward、训练更新、DEV或sealed读取；只产生`outputs/M28-FULL-GRAMMAR-CONTINUATION-TRAIN-CHECKS/summary.json`等检查收据，正式七方法汇总目录仍不存在。该短检查已由主对话核实并消费，不加入长流程心跳，不重复验证未变化的汇总代码。下一步等待已提交五条fit及完整TRAIN评价闭合，按预定全部七权重运行汇总，再作研究判断。

### 2026-10-11七权重正式汇总完成

七份原评价均终态并由主处理后，正式运行预定汇总13493626，private-teodoro-gpu／2CPU／4GiB／5分钟、无GPU、原generanno；UTC13:02:06–13:03:29，83秒，内部测量汇总0.288秒。root／batch／extern全部COMPLETED 0:0，stderr为空，STATUS／JOBID正确；summary107331字节、comparison2122字节，均非空。源代码未改，人工检查13454458未重复执行；提交前正式输出不存在，不覆盖历史结果。

正式输出完整保留原七方法和22688参考，逐染色体参考元数据、原端点与micro／macro／worst算术验证通过。主另核七个实际source job、末权重与步数、各物种／pooled计数直接等于已处理评价及逐参考恢复报告，六种对控制的差值重算相同。未重新读取DNA／注释／预测／权重，不新增baseline或模型推理、训练更新、DEV或sealed访问；只汇总已完成测量，不将相同参考元数据称为新的参考几何验证。

pooled链F1依预定顺序为原4608 45.5901%、STRUCT 55.7453%、LOCAL 41.7290%、HALF 48.6284%、双头先验64.3390%、FRESH63.0729%、JOINT66.9628%；精确计数和逐物种值以实际summary为准。[结果工作稿](../../../../paper/manuscript/M28_full_grammar_results_working.md)已保留全部分支、成熟基线及补充诊断。七权重、单列组合／完整核心和冻结八／十方法DEV家族仍分开；不得把TRAIN最优、汇总完成或描述性消融当作独立泛化、机制因果或论文goal成功。本短汇总已主直接消费，不增添长期监控对象。

## 混合目标末权重拟合与核验已完成

2026-10-09 UTC20:41主对话核对：HALF_STRUCTURED_HALF_LOCAL拟合13451137_2实际UTC07:50:23–20:35:53、45930秒，CPU核验13451245实际UTC20:35:54–20:36:15、21秒；两者root／batch／extern均COMPLETED 0:0，STATUS及核验summary齐全。正式末权重仍为该分支`step_009216.pt`，未选择中间checkpoint。

核验确认全部新增4608次更新、累计9216步、原draw顺序／fine及coarse观测／importance一致，混合目标及归一化算术一致；6144、7680、9216三个checkpoint的权重和AdamW矩有限，末checkpoint严格加载。核验本身没有optimizer更新、native forward、DEV或sealed读取。它证明训练记录和权重可加载，不证明完整链准确率或结构／局部贡献平衡。

既有raw数组13451471已在UTC20:36:16通过依赖自动启动。UTC20:40核查时，任务0–2运行，3–5按预定最多并发3等待；评价13451624保持有效`afterok:13451471_*`依赖。没有重新提交、取消或改变资源／参数。STRUCTURED_CONTINUE当次仍运行、日志到8704步，故整个三目标cohort尚未闭合。混合分支的六染色体raw和完整共同尺评价也尚未完成，不报告新F1或由工程完成宣称研究成功；继续复用原固定监控与全部冻结候选。

来源：Slurm终态及`outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-VERIFICATION/HALF_STRUCTURED_HALF_LOCAL/summary.json`。下一步消费该分支全部raw与评价终态，再与同更新数STRUCTURED_CONTINUE及其他预定末权重作完整比较，不提前调混合系数、选模型或进入DEV。

## 混合目标原 TRAIN 推理的并发路由调整

2026-10-09 UTC20:50现场调度显示私有gpu034／035的16张RTX3090均已分配；共享gpu021具有8张RTX3090、每卡24GB，UTC20:59再核时为IDLE，126个有效CPU及503000 MiB内存未分配。三个尚未启动任务各需1 RTX3090、2CPU、16GiB、12小时，共享分区上限恰为12小时，没有相交GPU维护预约，原六个目标输出路径互不相同。依据用户已明确的并行和优先加速授权，仅调整这一批已提交任务的路由与并发，不增加实验或训练更新。

UTC20:59:39执行`scontrol update JobId=13451471_3,13451471_4,13451471_5 Partition=shared-gpu`及`scontrol update JobId=13451471 ArrayTaskThrottle=6`，均成功；[Slurm官方接口](https://slurm.schedmd.com/scontrol.html)支持以数组子任务身份更新分区及数组并发上限。既有0–2保持在私有节点运行，3–5于UTC20:59:40在shared-gpu／gpu021首次启动；UTC21:02核对六任务均RUNNING，新三项实际分配仍各为1 RTX3090／2CPU／16GiB、12小时。其STATUS和started.json绑定预定三个Rice染色体、HALF_STRUCTURED_HALF_LOCAL末`step_009216.pt`、source fit13451137，reference_used和DEV_used均false。

数组身份13451471、原六任务、科学命令、完整双链raw范围、caller参数、checkpoint、环境、每任务资源及停止时限均不变，没有取消、重提、重跑已启动任务或修改SBATCH源码；其他数组仍沿用原并发设置。评价13451624的原整数组afterok依赖保留，未提前读取参考、DEV或sealed。跨节点／分区的调度与计时只作为实际执行记录，不构成公平速度基准或准确率收益。后续仍由原固定监控跟踪同一cohort，完整链结论须等待全部推理和评价闭合。

## 结构化对照原 TRAIN 推理的并发路由调整

2026-10-09 UTC21:16现场调度显示共享gpu017尚有4张、gpu021尚有5张RTX3090可分配，私有gpu035尚有2张；三个节点每卡24GB，有足够的CPU和主机内存，没有相交GPU维护预约。结构化对照raw数组13451470的六项任务全部尚未启动，仍等待CPU权重核验13451244；六个目标输出目录均不存在。只调整这批已有任务的排程，不新建实验、增加训练或绕过核验。

UTC21:17:59执行`scontrol update JobId=13451470 Partition=private-teodoro-gpu,shared-gpu ArrayTaskThrottle=6`，成功回读数组0–5的并发上限为6、可用分区为上述两个分区。每项请求仍为1 RTX3090／2CPU／16GiB、12小时，保持原排除节点和generanno环境。Slurm支持[多个候选分区择可启动分区](https://slurm.schedmd.com/sbatch.html)及[调整已有数组的并发上限](https://slurm.schedmd.com/scontrol.html)；未来实际分区、节点和同时运行数以启动后的分配为准，不由当前空闲数保证。

原`afterok:13451244`依赖和评价13451625的`afterok:13451470_*`依赖均保留，正式权重仍为STRUCTURED_CONTINUE的`step_009216.pt`、来源fit13451136。数组身份、六TRAIN整染色体双链范围、科学命令、caller、每任务资源和时限未改，没有取消、重提、重跑或修改SBATCH源码；其余尚未启动分支的raw数组保持原并发设置。该调整证明排程已配置，不证明结构化拟合、推理或评价已完成，也不构成公平速度或准确率结果。继续复用原固定监控，不新增监控对话、自动化或DEV／sealed读取。

## 混合目标已完成的原 TRAIN 部分结果

2026-10-09，HALF_STRUCTURED_HALF_LOCAL末9216步的六项整染色体双链raw推理13451471_0–5及CPU评价13451624全部root／batch／extern均COMPLETED 0:0。六项raw最后于UTC21:34:29结束，评价UTC21:34:29–21:35:36运行67秒。各项STATUS、JOBID、外层summary和非空预测齐全，均绑定对应原TRAIN染色体、该分支正式末权重和来源fit13451137；CPU评价已核对完整双链窗口覆盖、原caller身份和GFF往返计数、DNA／phase／cap及同链不重叠约束。推理未读参考、标签或feature cache，没有训练更新、DEV或sealed读取。

评价仍保留AT17739／Rice4949／pooled22688条主参考及原unsupported层。相对原final4608和已完成LOCAL_CHANNEL评价，六染色体参考数、eligible isoform数、CDS参考碱基、背景分母和unsupported元数据一致；主对话按计数重算逐染色体、逐物种、pooled完整链PRF，逐物种及pooled CDS-base PRF和macro／worst，与保存值一致。表中碱基F1是不区分strand的CDS union，背景链是完全位于双链所有注释gene span之外的唯一预测链。

| 原TRAIN范围 | 原4608链F1 | 混合9216链F1 | 混合完整链TP／预测／参考 | CDS碱基F1 原4608→混合9216 | 纯背景链 原4608→混合9216 |
| --- | ---: | ---: | --- | ---: | ---: |
| AT | 50.4871% | 51.6731% | 9196／17854／17739 | 93.0848%→95.9291% | 15→175 |
| Rice／worst species | 28.2780% | 38.3270% | 2016／5571／4949 | 80.8511%→88.4729% | 147→537 |
| pooled | 45.5901% | 48.6284% | 11212／23425／22688 | 90.3601%→94.1577% | 162→712 |

混合目标的pooled完整链F1提高3.0383个百分点，Rice提高10.0490个百分点，macro达到45.0000%；但背景链增加550条，pooled背景span FPR达到1.1435%。pooled TP净增加1246是两个汇总计数的差，不能当作逐参考“恢复1246条”的配对结果；原已提交恢复分析13461885负责区分新增命中与丢失命中。重复链记录均为0，747条非primary eligible isoform精确命中不计入primary TP。

相对同为末9216步的LOCAL_CHANNEL，混合目标AT／Rice／pooled完整链F1分别高7.5128／4.9452／6.8994个百分点，但pooled CDS-base F1仅高0.0611个百分点，背景链再增加221条。这支持混合训练在当前TRAIN完整链终点上的部分改善，不证明目标贡献已平衡、结构模块必要性或泛化。本次状态核查时，同更新数STRUCTURED_CONTINUE仍在拟合，其raw及评价为有效依赖等待；三个目标的完整比较尚未闭合，不能用少4608次更新的原模型作纯目标因果控制，也不据此调混合系数或选择权重。

同一原TRAIN共同尺的固定原生比较器ANNEVO／Tiberius／Helixer pooled完整链F1分别为75.5418%／72.6555%／63.5919%，仍比本混合分支高26.9134／24.0271／14.9636个百分点。比较器的上游训练暴露并未匹配，故这是当前固定配置的相对性能，不是清洁独立测试；混合目标的改善仍未达到论文goal所需竞争力。

本部分结果不改变冻结八／十方法DEV家族、阈值、数据隔离或其他在运行分支，封存／held-out／Setaria不放行。继续等待已提交完整分支及原配对恢复分析，再按既定全部七末权重汇总，不重复局部推理或追加控制。原固定监控保留原三目标cohort未交接标记；本分支已由主对话消费，不等于整组流程或论文goal完成。

结果来源：`outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC/HALF_STRUCTURED_HALF_LOCAL/summary.json`（评价13451624）、原final4608评价13443329、LOCAL_CHANNEL评价13451623及`outputs/M28-NATIVE-TRAIN-BASELINE-EVAL/summary.json`（评价13461307）。完整原生基线说明见[native_train_baseline_result.md](../../../../benchmarks/native_baselines/reports/native_train_baseline_result.md)。以上历史拟合与排程段落保留其当时状态，本节记录后续已完成的推理评价。

## 结构化对照拟合与权重核验完成

2026-10-09，STRUCTURED_CONTINUE拟合13451136_0于UTC22:14:22完成，CPU核验13451244于UTC22:14:24–22:14:41运行17秒；两项root／batch／extern均COMPLETED 0:0，STATUS和JOBID一致。拟合新增4608次更新，累计9216步，正式权重仍仅取该分支末`step_009216.pt`。核验确认原draw顺序、fine／coarse观测、importance及目标算术一致，6144／7680／9216三个checkpoint的权重和AdamW矩有限，末权重严格加载；核验没有训练更新或native forward，拟合与核验均未使用DEV或sealed。

既有raw数组13451470的六项任务已通过核验依赖启动。UTC22:18调度及启动收据确认全部RUNNING，逐项绑定原六TRAIN染色体及其完整长度、STRUCTURED_CONTINUE末9216步权重和source fit13451136，reference_used和DEV_used均false。实际numeric root按任务0–5为13465140／13465141／13465142／13465143／13465144／13451470，与各started.json匹配。评价13451625仍等待整数组有效依赖。此次完成的是拟合记录与权重核验；完整链准确率以及混合目标相对同更新数结构化对照的收益，仍须等待全部raw及共同尺评价闭合。原固定监控、其他分支和冻结DEV候选不变，不增补控制或据训练损失选权重。

来源：上述Slurm终态、`outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST/STRUCTURED_CONTINUE/summary.json`、`outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-VERIFICATION/STRUCTURED_CONTINUE/summary.json`及六项raw启动收据。三目标cohort仍未完成，论文goal仍未达成。

## 原 TRAIN 三目标完整比较

2026-10-09，STRUCTURED_CONTINUE六项raw13451470_0–5全部root／batch／extern为COMPLETED 0:0，最后一项于UTC22:52:17结束；评价13451625于UTC22:52:19–22:53:22运行63秒，同样全部COMPLETED 0:0。六项STATUS／JOBID／外层summary逐项对应原染色体、完整长度、正式末9216权重及来源fit13451136，推理没有参考、标签／feature cache、训练更新、DEV或sealed使用。CPU评价在全部推理完成后核对双链完整网格、caller身份、GFF往返计数、DNA／phase／cap及同链不重叠，再读取参考；三目标TRAIN流程至此全部闭合。

主对话按原计数重算新分支的逐染色体、逐物种和pooled链及CDS-base PRF、物种汇总、背景汇总、macro及worst，均与保存值一致。与原4608、LOCAL_CHANNEL及混合目标相比，六染色体的主参考数、assessable isoform数、参考CDS碱基、背景分母和unsupported元数据一致，仍保留AT17739／Rice4949／pooled22688条主参考。

| 正式末权重 | AT链F1 | Rice链F1 | pooled链F1 | pooled CDS-base F1 | 纯背景链 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原结构目标4608 | 50.4871% | 28.2780% | 45.5901% | 90.3601% | 162 |
| STRUCTURED_CONTINUE 9216 | 61.8356% | 35.4129% | 55.7453% | 93.2001% | 197 |
| LOCAL_CHANNEL 9216 | 44.1602% | 33.3818% | 41.7290% | 94.0966% | 491 |
| HALF_STRUCTURED_HALF_LOCAL 9216 | 51.6731% | 38.3270% | 48.6284% | 94.1577% | 712 |

结构化续训命中12696条主参考，预测22862条唯一链；pooled P／R为55.5332%／55.9591%，macro为48.6243%，worst为Rice35.4129%。相对原4608，AT／Rice／pooled链F1提高11.3485／7.1349／10.1552个百分点，TP净计数增加2730，不等同于逐参考新增恢复2730条。1128条其他eligible isoform命中不计primary TP，重复链为0；pooled背景span FPR为0.2892%，纯背景链比原模型多35条。

同为9216步，结构目标比LOCAL的AT／Rice／pooled链F1高17.6754／2.0311／14.0164个百分点，虽然pooled CDS-base F1低0.8965个百分点，纯背景链少294条。混合目标相对结构目标的Rice链F1高2.9141个百分点，但AT低10.1625、pooled低7.1170个百分点；CDS-base F1高0.9576个百分点，纯背景链多515条。因此，混合目标不是整体链准确率改善：三目标中结构目标的AT与pooled链F1最高，混合目标的Rice及worst-species读数最高，须保留这种物种差异与背景代价，不据单一读数挑选DEV候选。

三者的步数、起点、draw顺序、观测、importance和caller相同，但目标及梯度尺度不同，固定0.5系数不代表有效更新幅度匹配。该完整比较支持当前TRAIN上的目标相关差异，不单独证明结构模块必要性、梯度平衡、因果机制或独立泛化。三个续训分支在两物种均仍低于固定原生基线；其中结构目标pooled仍比ANNEVO／Tiberius／Helixer低19.7965／16.9101／7.8466个百分点，论文goal未达成。

原固定监控已交接本cohort并暂停，主对话已完成三目标结果处理。下一步依冻结三目标DEV方案接入全部三个末9216候选，不等待无关TRAIN支路、不按TRAIN分数预筛，也不改八／十方法家族、阈值或数据隔离。其余七个cohort继续复用原监控；DEV实际提交与执行状态另记，不由本次结果处理推定已运行。封存／held-out／Setaria仍不放行，无公开发布或Git提交。

来源：`outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC/STRUCTURED_CONTINUE/summary.json`（评价13451625）、本记录中的已完成LOCAL／混合结果、原final4608评价13443329及固定原生评价13461307。上述较早段落保留当时状态，本节给出后续整组完成的比较。
