# 原 TRAIN 的编码器低秩适配器对照

2026-10-12 UTC。最终9216拟合、核验、六raw和配对整链已完成并由主核验；[完整结果](full_grammar_joint_lora_result.md)的pooled整链F1为72.958322%，比B提高2.655698 pp，两个物种和全部六染色体均提高，但背景负担上升。Tiberius pooled对照被略超出，不代表拟南芥或ANNEVO统一优胜。原TRAIN范围和冻结B/C/D DEV均不变。以下固定协议与历史提交事实保留。

## 固定科学对照

基线为B的fresh seed1／joint24576分支（fit13494895、原六TRAIN评价13494898）。新分支从相同源fit13434603的最终4608 C0头与完整AdamW状态开始，使用B已经冻结的三段各1536、共4608个fresh draw及4219个unique窗口；每一步的窗口、物种、方向、核心fine／coarse／unknown观察、有限context、ownership、importance1和fine+coarse分母都相同。只读最终累计9216／epoch6权重；6144、7680仅为训练快照，不做TRAIN或DEV模型选择。

相对于B的机制变化只有已预检的编码器低秩残差：固定GENERanno revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345`，26层中的零基22至25层，每层q_proj及v_proj，rank4／alpha4／scale1。八个投影有16个FP32张量，共102400个参数，A按CPU seed0均匀初始化，B从零开始。原366191参数头继续学习，总可训练参数468591；底座、embedding、RoPE、其他投影和MLP不训练。

使用独立可微joint路径，保留原24576bp、4096六聚体token、双向原生SDPA Flash Attention、BF16隐藏表示转FP32、方向后右侧N padding、非ACGT转N及原最近窗口中心ownership。编码器保持eval，无dropout、无gradient checkpointing。投影残差在FP32计算并舍入BF16后加到原BF16输出，不合并到原发布权重。此条件是固定实现，不把它称为独立的新颖性贡献。

完整gene-set合法域、空集、非重叠、cap128、起止密码子和实际跨内含子carry、FP64分区及epsilon0.1都不变，不加入codon能量或其他辅助loss。AdamW仍为lr0.0003、weight_decay0.01、betas0.9／0.999、eps10⁻⁸；头部moments及年龄4608延续至9216，适配器年龄从0至4608，所有训练参数共同norm1裁剪。新增容量、适配器moment年龄及其导致的联合裁剪／优化轨迹都是干预的一部分，不能把结果归为单一生物学因果效应。

## 权重与 raw 接入

每段完成保存一次C0头、16个具名适配器张量和全部AdamW状态，不保存1.2b底座副本。checkpoint记录原revision、训练配置、头与适配器不同年龄，以及适配器参数名到optimizer组的实际注册顺序。CPU核验4608行draw、观察计数、loss算术、三段summary、每份参数和moment的形状／dtype／有限性与年龄，最后用严格raw头加载器读取最终9216。

推理仅接收FASTA、最终checkpoint和输出目录，不读取训练catalog、参考或观察模型。它锁定发布revision，重新加载BF16底座，装入全部16个FP32适配器，关闭所有梯度，并复用原joint全contig caller和十五角色role-minus-background的精确skyline解码。旧冻结编码器、缺适配器、中间权重、错seed、层数、rank或年龄必须拒绝；不能把缺失适配器静默当作冻结模型运行。

## 原六条 TRAIN 评价

| 物种 | 完整染色体 | 长度 bp |
| --- | --- | ---: |
| Arabidopsis | NC_003070.9 | 30427671 |
| Arabidopsis | NC_003071.7 | 19698289 |
| Arabidopsis | NC_003076.8 | 26975502 |
| Rice | NC_089043.1 | 27474823 |
| Rice | NC_089044.1 | 23931887 |
| Rice | NC_089046.1 | 28271460 |

六个raw任务保留完整DNA、全部窗口和两个方向。全部完成后，先逐碱基比较两臂输入与原FASTA经同一大小写／非ACGT转N处理后的DNA，再允许新分支评价读取TRAIN参考。复用原M26 CDS-assessable最长合格CDS／transcript ID tie规则、phase与输出roundtrip核验以及count-derived评价，保留全部22688 primary参考和不受支持分层；训练资格不筛除评价参考。

配对报告保留逐染色体、逐物种、pooled、macro和worst的完整链与CDS-base-union计数及P／R／F1，负结果和物种差异同样保留。macro分别按两物种计算；worst先在各臂内取最差物种再作差，不能以最小逐物种差代替。base收益不得替代完整链收益。单次同起点追加训练只回答原TRAIN条件下的容量对照，不证明独立重复、泛化或论文达标。

## 验证与作业范围

新增六组人工检查覆盖最终身份与16张量形状／dtype、实际CPU原生接口fixture的参数恢复和底座冻结、具名AdamW注册与两种年龄、raw适配器身份，以及配对分母／聚合计数和macro／worst算术。加上原五组低秩与可微joint路径回归，共11组，经无GPU Slurm检查；不读取研究DNA、真实权重、参考或原生模型文件。

| 阶段 | 固定资源 | 时限及依赖 |
| --- | --- | --- |
| 人工检查 | 无GPU，2CPU／4GiB | 5分钟 |
| 追加fit | 单24GiB RTX3090，2CPU／16GiB | 20小时；启动核对原生预检13499640与人工检查 |
| 权重核验 | 无GPU，2CPU／4GiB | 5分钟；afterok fit |
| 六TRAIN raw | 每项单RTX3090，2CPU／16GiB；最多3项并行 | 每项2小时；afterok权重核验 |
| 全范围评价 | 无GPU，2CPU／4GiB | 20分钟；afterok整个raw数组及B评价13494898 |
| 配对 | 无GPU，2CPU／4GiB | 10分钟；afterok新评价 |

全部no-requeue，GPU经private-teodoro-gpu执行。B的近期更新速率约10秒／step，适配器长程训练速度尚未实测；20小时为新增反传留执行余量，不是预计完成保证。正式fit首两步及每64步记录实际更新耗时和GPU峰值，供时限风险判断，不能据NLL调整科学参数或训练长度。fit输出上限96MiB；每条raw预测输出上限128MiB，六条输入FASTA合计约150MiB，整链新增存储预计小于2GiB。提交前核实实际队列、VRAM、并发和存储，不取消已有作业。

实现与输出前缀为JOINT-LORA-FIT、JOINT-LORA-VERIFICATION、JOINT-LORA-TRAIN-RAW、JOINT-LORA-TRAIN-DIAGNOSTIC及JOINT-LORA-PAIRED，各阶段不能覆盖既有结果。短暂连接失败只重查同一handle；真实失败先保留并定位，至多一次不改变科学配置的实现／调度修复，其后的重试仍保留失败输出。科学负结果不属于可重试的工程失败，不自动增训、换seed、层数、rank、epsilon或样本。

整链复用原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3和唯一coding-rna-r5-light-monitor心跳；正常运行及依赖等待保持安静，终态或真实actionable问题向原主对话交接。该分支当前仅限原TRAIN，不向已准备的B／C／D六方法DEV家族添加候选，不解析、保留或使用DEV或封存DNA／注释，不新增物种或面板，不对外发布。

## 实际实施与提交记录

本节以下为历史启动记录，不再代表当前仍运行。完整终态：fit13499817于2026-10-12 08:08:22UTC完成，verify13499818于08:09:39完成，六raw全部08:47:26前完成，eval13499823于08:49:24、paired13499824于08:49:49完成；全30条actual root/batch/extern COMPLETED0:0。实际raw throttle已为6，原%3历史脚本未重提。原watcher08:58:12UTC真实终态交接成功，主核验并记录完整科学结果；不再监控E。新的只读恢复交集诊断在E结果已知后规定，不更改训练或DEV家族。

以下时间均为UTC。人工检查13499807于2026-10-11 18:49:50至18:50:06完成，root／batch／extern均COMPLETED、ExitCode0:0；11项全部通过，零失败、零错误，检查函数用时0.519916972秒，batch峰值RSS954196KiB。源文件AST和六份作业脚本的Bash语法检查亦通过。这些是实现准备证据，不是最终9216权重或完整染色体准确率结果；原生预检13499640的两个更新仍是可丢弃检查，不能代替正式训练。

| 实际阶段 | Slurm handle | 实际依赖 |
| --- | --- | --- |
| fit | 13499817 | afterok人工检查13499807；脚本核对原生预检13499640 |
| 权重核验 | 13499818 | afterok 13499817 |
| 六TRAIN raw | 13499821_0–5，最多3项并行 | afterok 13499818 |
| 全范围评价 | 13499823 | afterok整个13499821数组及B评价13494898 |
| 配对 | 13499824 | afterok 13499823 |

fit于2026-10-11 18:51:13实际启动。18:56UTC核对为RUNNING／gpu034；scontrol确认单RTX3090、2CPU／16GiB、20小时、Requeue0。其余四阶段为PENDING，依赖关系与上表一致；尚未完成任何最终权重、raw输出或评价。singleton FIT／VERIFICATION／TRAIN-DIAGNOSTIC／PAIRED目录不含job-ID子目录；六条raw各用seqid子目录，人工检查按job-ID保存，不能混淆。

实际started.json核对了相同源fit13434603／support13453408、原生预检13499640／人工检查13499807、原seed1的4608 draw／4219 unique、joint24576、epsilon0.1以及新增16张量／102400参数配置。首两个正常fit更新已经记录：头部step4609／4610对应适配器step1／2，每步用时16.821499／10.616146秒；GPU峰值allocated6826383872字节、reserved7310671872字节，读取时stderr为空。两步耗时仅用于执行容量判断，不外推为全程完成保证，也不据此改变科学配置。

原唯一心跳的保存配置已核对：仍为ACTIVE、绑定原监控对话、每15分钟；完整原B+C+D prompt保留，只追加E对象，总长16776字符，其余字段除更新时间外不变。新增E的实际定时追加亦已核验：原监控对话turn01a12c54-0502-7c32-a9de-2b5bef96cdff于18:57:46.291启动，仍为gpt-5.6-luna／low；18:57:46.295的心跳输入含与保存配置完全相同的四对象prompt，调度时间18:57:46.241。其18:57:55实际检查四个fit均RUNNING、全部下游为有效依赖等待，SSH失败0，无科研读取、交接或暂停；18:58:02静默完成。配置成功与一次真实定时运行分别有证据，最终回传仍未发生。B、C、D、E均未终态交接，A已关闭且不再监控。
