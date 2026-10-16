# M28 A10 匹配整窗表示实验准备

2026-10-03。整窗工程预检 13359555 已完成，八窗四块 BF16 缓存精确复现，同 token IDs 的 4,096-token 整窗编码可运行。接下来检验完整窗口的 encoder 表示能否改善同一 query 模型，而不是把旧头直接接入不同特征。用户已授权继续迭代；本轮推进 TRAIN 特征准备、匹配协议和一次固定拟合，不改原 A3 负结果、FULL Q/C0、A6/A7 或 A9。

## 匹配条件

直接对照为已完成 A3，不是 pending FULL 权重。保持 A3 的 CollaborativeGeneQueryModel、649,457 参数、D2048/H128/Q16、四头 pre-norm 残差通信、dropout0/无 FFN 和 known-query-mean 目标；边界、EOS、local、Hungarian、R4 正/保护/未知政策及四损失系数 1 不变。新头应从同一 seed0 初始化完整从头训练，不继承 A3/A6/A9 的已监督参数。

原 1,536 TRAIN draws、1,491 唯一窗、两物种各 768 draws、顺序及 importance_weight_uniform_within_stratum 原样保留。三遍共 4,608 更新，batch1、FP32、AdamW lr3e-4/wd.01/betas(.9,.999)/eps1e-8、clip1；只用最后 step_004608，不选中间 checkpoint、seed 或额外 epoch。六条 TRAIN 染色体和原同源/封存隔离不变。本轮仅更换冻结编码的上下文模式；跨块注意力与位置处理同时变化，不能称单因素因果实验或新的 query 架构。

模型仍为同一 GENERanno preview/revision b0483c23b6b63787b61a6d3a204a9b517d6ba345，SDPA/BF16/eval/frozen/local_files_only。新 feature_spec.encoding_mode 为 one_4096_token_forward、block_bp=24576，特征形状不变。形状及 revision 相同不足以识别新旧表示，后续训练入口必须核对完整 feature_spec，不能从旧路径混入四块特征。

## 此轮实际执行：TRAIN cache

cache_draw_context.py 从原 draws 和已经完成的 full-TRAIN 几何计划定位上述 1,491 个旧 TRAIN 缓存，只读取原 oriented DNA/几何/shape；以 mmap 避免无用地顺读约 24 GiB 的旧 tensor payload。不读取参考 catalog、query 输出、DEV cache 或 test/Setaria，也不执行头/optimizer。每窗保持旧 TRAIN-cache 的非 ACGT→N 规范化和 orientation 后右侧 N padding，要求整窗与四块 token ID 完全一致，再保存有限 4,096×2,048 BF16 新特征。

严格按原 draw 首次出现顺序为唯一窗口编号，draws.jsonl 保留完整原列表、顺序及权重，index.jsonl 保留来源/窗口/新 artifact/实际大小。与工程预检重合的所有 draw 窗口必须精确复现其整窗 BF16 特征；它用于拦截新入口的表示不一致，不检查准确率。不覆盖旧缓存，全部写 outputs/M28-A10-DRAW-CONTEXT-CACHE，预计约 23.4 GiB 特征，准备不少于 26 GiB 可用存储。

Slurm 脚本 M28-A10-DRAW-CONTEXT-CACHE.sbatch：generanno，1 RTX3090、2 CPU/16 GiB、15 分钟，提交前按实际队列选择 eligible partition。八窗整窗同步执行约 0.234–0.239 秒/窗，1,491 次前向粗略约 6 分钟，另需原文件元数据和新 tensor 写入；时限含操作余量，不是完成保证或全 DEV 时长预测。非有限、token/几何/表示错配、OOM 或超时立即停止并保存失败，不自动续跑、延时、换 dtype/encoder 或混用原 feature cache。

## 拟合和后续判断边界

只有上述缓存完整、draw/标签 join 与新头零更新工程检查通过后，才进入该固定配置的一次拟合。run_query_context.py 复用 A3 的 seed0 构造器、TRAIN-only 标签 loader 与原 supervised_readout，独立加载整窗 index，并严格核对原 draws/完整 feature_spec/几何/实际 artifact。工程模式从原 draw 中取每个十二 TRAIN stratum 的首个既有窗口，进行原目标的前向和反向，要求梯度有限、至少一个窗口有非零头梯度、输入无梯度，整个头 state_dict 未变化；合法零监督窗保留并记录，不因零梯度删样本。不创建 optimizer 或读取监督 checkpoint。还以相同有效 tensor/窗口/revision、仅更换为旧 context metadata 的输入检验实际路由会拒绝新旧特征混用。

工程作业 M28-A10-CONTEXT-ENGINEERING.sbatch 为 1 RTX3090、2 CPU/16 GiB、5 分钟，输出 outputs/M28-A10-CONTEXT-ENGINEERING，必须 afterok 缓存。拟合 M28-A10-CONTEXT-FIT.sbatch 使用同资源、原 A3 的 90 分钟操作时窗（旧实际约 31 分钟），新输出 outputs/M28-A10-CONTEXT-FIT，保存三份 epoch 权重且只允许最后 step004608；已在下面两项工程终态核实后提交 13359631。完成末权重后才准备完整 raw-FASTA 推理，不将旧四块 DEV 特征送入新头。自然网格、8690 方向窗、tau0、max_exons128、DNA 语法、owner/去重/正 gain 非重叠选择及 GFF3 phase 沿用原规则；7728 主参考（5437+2291）与全部 unsupported/长链/重叠条目保留，生成完毕后才读取参考评价。

正式结果同时报告与 A3 的逐物种及 pooled/macro 变化、背景、CDS-base union、free→owner→selected 损失及原生 ANNEVO/Tiberius/Helixer 对照。原 A3 的主要继续条件不降低：两物种 exact-chain F1 均严格高于原 A，且两物种背景 span 均不高于原 C0，以整数计数判断；相对较弱 A3 的增量仅是机制比较，不能替代该条件或原生性能排名。反复使用的两物种 DEV 仍是模型开发集，未形成同源隔离独立测试或广泛适用范围证据。

当前拟合尚未完成 4,608 更新，无 raw-FASTA 注释结果或科学 PASS。后续失败保留负结果，不扫阈值/seed、追加曝光或更改分母救本轮；即使开发继续条件通过，论文目标也需要完整原生比较、独立测试与创新证据。新方案是表示机制对照，不因使用长上下文就构成可发表创新；GEN 预训练/监督曝光未知边界继续保留。无 Git 提交、对外发布或封存数据开放。

## 执行记录

缓存 13359590、工程检查 13359607 的 root/batch/extern 均 COMPLETED/0:0，分别 9 分 36 秒和 32 秒，STATUS/JOBID 对应一致；private-teodoro-gpu/gpu035、各 1 RTX3090/2 CPU/16 GiB。实际命令分别 `sbatch --parsable sbatch/M28-A10-DRAW-CONTEXT-CACHE.sbatch`、`sbatch --parsable --dependency=afterok:13359590 sbatch/M28-A10-CONTEXT-ENGINEERING.sbatch`。提交前五项 live 调度检查、真实存储及目录不碰撞核实，缓存前 BeeGFS Available 174,147,048,046,592 bytes。新增 source/sbatch/protocol remote 精确读回与本机一致，AST/bash 语法通过；旧预检、原 A3 代码和所有运行中旧实验未修改。

缓存实际 1,491 份非空 BF16 artifacts，总 25,054,280,916 bytes（23.3336 GiB），物种唯一窗为 745+746。主对话核对 index 1,491 行、draws 1,536 行及真实 feature 文件数量/大小加总，与 summary 一致；后续工程 loader 逐一核对完整新 index、原 draw 列表/顺序/权重、R4 TRAIN 标签连接、窗口几何及所有 artifact 的实际大小。唯一与原八窗预检重合的 draw 窗口 arabidopsis_thaliana|NC_003070.9|15212544|- 整窗特征精确复现；不以一个 witness 声称所有数据重复独立验证。原模型 checkpoint/标签 catalog 无编码参与。

缓存循环 550.2860 秒，逐窗编码/六-mer 对齐/材料化计时 0.2627–2.9084 秒，峰值 allocated 2,522,416,640 bytes。这个计时包含 token 化，和八窗预检已先准备 token 的 whole-window 区间不同，不进行速度排名。十二 stratum 真实原目标梯度检查循环 13.7642 秒，12/12 梯度有限且非零（norm 12.4700–63.9840），零梯度窗 0，头 state_dict 未变，旧 context 正确拒绝，matched 参数 649457。两项均 optimizer0、无新 DEV/test/Setaria 使用；工程模式无监督 checkpoint 或新 backbone 前向。

上述终态由主对话消费后，提交固定拟合 13359631，实际命令 `sbatch --parsable --dependency=afterok:13359607 sbatch/M28-A10-CONTEXT-FIT.sbatch`。提交前重新做五项 live 调度与输出存储/不碰撞核实；私有 gpu035 再有可用 3090，选择 private，Available 174,115,691,429,888 bytes。实际已 RUNNING gpu035，1 RTX3090、2 CPU/16 GiB、90 分钟，STATUS RUNNING/JOBID 正确；启动不是完成拟合或新性能结果。后续需末 step004608、完整 raw-FASTA 8690 窗及 7728 主参考评价，当前没有推理/评价作业，不以中间训练 loss 选择下一配置。

固定监控对话仍为 01a0f71f-8845-7432-bd13-bf54c1808fb3、同一 heartbeat coding-rna-r5-light-monitor，当前 cohort M28-FULL-A6-A7-A9-A10-20261003，仅在旧 11 根作业之外增加 A10fit13359631。整份 prompt/ID/kind/ACTIVE/target/15 分钟/通知策略/created_at 保存读回精确匹配，无 model/thinking override 或新 thread/cron。A10 固定拟合阶段终态是准备后续 raw 推理的接续事件，不是整条注释流程完成；已关闭 cache/工程不重复监控。当前扩展 cohort 尚未观察定时运行，未发送人工启动消息，不以旧 13:42 定时检查替代新对象验证。普通运行静默，真实终态或故障交接一次，主对话负责结果判断。定时监控依赖本机及应用在线，见[官方任务说明](https://learn.chatgpt.com/docs/automations?surface=app)。

## 完整 raw-FASTA 入口与后续作业（14:07 UTC）

在固定拟合正常运行时，完成不依赖其性能的推理／评价准备并提交下游；下面替代上文“没有推理／评价作业”和 fit-only 交接状态，不改变训练配置或科学继续条件。

infer_query_context.py 只允许对应 completed fit 的 step_004608，核对完整 configuration、revision、三遍更新、649457 参数及末权重政策。每个自然方向窗口复用缓存已验证的 encode_full_window：同六-mer IDs、一次4096-token frozen native forward、BF16材料化后转原FP32头，直接从 raw FASTA 大写／非ACGT→N规范化及正反方向 DNA 编码，不读参考、catalog 或任何 DEV 特征 artifact。同一 Q16 greedy、tau0/max_exons128、owner／去重／正gain非重叠选择及既有 exporter 不变；另有显式 --fasta 入口，但它不能进入固定 DEV 评价。长度／窗口／四scope／逐query／assembly收据及原GFF保存到新输出目录。

evaluate_query_context.py 首先核对 terminal fit/infer 与匹配完整自然grid，然后流式核对全部实际 raw record 身份、Q16／tau0、completed／positive数量，重新组装并比较原 assembly counts、最终预测数与 GFF3 identity／坐标／phase；完成这些后才调用 reference_context。没有新增参考过滤，仍保留7728主链和全部unsupported／长链／重叠。记录逐物种及 pooled/macro 的 exact链、CDS-base union、背景、free→owner→selected损失和query重复。主要gate直接复用原A3的整数判定，对原A和R5_C0；另记录匹配A3及本次完成的原生ANNEVO/Tiberius/Helixer差值，不拿旧 cached comparator 代替新 native。六份已完成 comparator 的实际summary字段及A3匹配曝光已读回核实。

CPU 工程检查13359719的root/batch/extern均COMPLETED/0:0，17秒；六项测试PASS、0.653秒。检验 BF16→FP32 入口、raw规范化和两方向单窗前向、带模糊flank的minus半开坐标／phase真实GFF往返、非正gain不输出基因、旧context／中间权重／显式非DEV／残缺grid拒绝，以及实际记录残缺时参考尚未打开。无 native 编码、DEV/test/Setaria 或optimizer。既有共享parser产生ResourceWarning但不影响通过，未扩大修改该parser。六个新增source／test／sbatch本机AST／bash检查通过，部署后整份local/remote diff精确匹配，不覆盖旧运行代码。

| 阶段 | Slurm ID／当前状态 | 固定资源／依赖 | 输出 |
|---|---|---|---|
| 原固定fit | 13359631／RUNNING gpu035 | 1RTX3090、2CPU／16GiB、90分钟 | outputs/M28-A10-CONTEXT-FIT |
| raw全8690窗infer | 13359725／PENDING有效Dependency | 同GPU／CPU／memory、3小时；afterok:13359631:13359719 | outputs/M28-A10-CONTEXT-INFER |
| 全7728主链eval | 13359727／PENDING有效Dependency | 0GPU、2CPU／8GiB、20分钟；afterok:13359725 | outputs/M28-A10-CONTEXT-EVAL |

实际命令为 `sbatch --parsable sbatch/M28-A10-CONTEXT-RAW-CHECKS.sbatch`、`sbatch --parsable --dependency=afterok:13359631:13359719 sbatch/M28-A10-CONTEXT-INFER.sbatch`、`sbatch --parsable --dependency=afterok:13359725 sbatch/M28-A10-CONTEXT-EVAL.sbatch`。提交前按Skill做五项live调度、目录／log不碰撞及存储核实，infer/eval前Available 174086494355456 bytes。私有GPU当时已满，但原fit成功释放3090即满足本轮infer的一张GPU需求，故沿用private；下游开始时间未由Slurm给出，不保证具体完成时间。3小时沿用原greedy full-Q caller时限，不拿A6的beam8小时或整窗encoder单步计时推定实际完成。没有取消、延长或修改其他作业。

同一心跳与固定监控对话更新为 M28-FULL-A6-A7-A9-A10-RAW-20261003，保留旧11根，加A10 fit/infer/eval共14根；消耗完的checks/cache/engineering不监控。现在A10必须整条fit→infer→eval完成才交接，单fit成功且下游正常等待／运行保持静默。ENTIRE prompt、ID／name／heartbeat／ACTIVE／target／15分钟／通知策略／created_at精确读回，未改模型或推理强度。先前fit-only cohort已在原对话自然运行：13:57 UTC、turn01a1020e-5f36-7d80-bae1-41eb2ff5d9ad，completed／error null／cursor49；不能用它替代新下游对象的定时运行验证。本次没有发人工startup消息，新完整raw cohort的定时执行与终态回传仍待自然触发。配置成功、一次检查成功与完整结果回传分开记录。

当前只有工程PASS与真实下游提交，尚无A10全DEV预测或性能结果，论文目标未达到；新结果由主对话按原gate、原生排名及独立测试／创新证据边界解释。

## 完成后的逐参考诊断准备（14:16 UTC）

新增 diagnose_query_context.py，在 A10 fit/infer/eval 全部 terminal COMPLETED 后只读原 raw候选，复用已接受的 outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC/references.tsv 中完整7728主参考几何、A3阶段及三种原生工具 exact成员身份。新诊断不重跑模型，不读新FASTA／DEV feature／监督checkpoint，不新增训练、参考筛选或test/Setaria访问，也不改变原操作点／gate；重放只用于获得 A10 逐链阶段和owner最大gain，而不是重复审计未变化的原A3代码。

按每条参考链记录 A3→A10 的 free、owner、selected gross新增／丢失及五类状态：未生成、owner丢失、owner gain非正、正gain未选、最终selected。与旧表 join 保留species／seqid／strand／ordered半开CDS坐标，phase不属于exact-key；重放出的四scope stage counts、assembly counts及全部selected TP必须与已经完成的A10正式评价一致。针对CSV `False` 字符串、重复query的owner最大gain和nonowner-only丢失做必要检查，不增加哈希或新验收框架。

预先固定描述性分层：参考CDS片段数1／2–4／5–9／10+，基因组链span≤6144／6145–24576／>24576 bp，以及已记录unsupported有／无；每种分层均互斥覆盖全7728，不能删掉长链或难例。6144 bp是原编码块宽度，跨度分层不证明跨块注意力的单因素作用，因为位置处理也改变且基因长度／复杂度相互关联。只报告参考命中率与成员交集，不提供子组precision或F1，不用分层好结果替代完整主评价；A10独有精确命中仍是既有注释参考，不是新基因发现。

工程检查13359759 root/batch/extern均COMPLETED/0:0，16秒；三项CPU合成检查PASS、0.022秒，未读取真实数据或做native前向。四份新增source／test／sbatch的AST／bash与整份local/remote diff精确匹配通过；旧诊断与正在运行的fit/infer/eval源未改。正式诊断脚本 M28-A10-CONTEXT-DIAGNOSTIC.sbatch 使用generanno、private／0GPU／2CPU／8GiB／10分钟，单独新目录outputs/M28-A10-CONTEXT-DIAGNOSTIC，afterok原A10eval及检查13359759；现阶段只有诊断准备PASS，没有新增性能结果。

原14根完整raw cohort已观察到自然定时检查：14:12 UTC、同监控对话turn01a1021c-1b11-7ee0-8684-f0476b5894ad，completed／error null／cursor50，未交接。此前“尚未观察新raw cohort定时运行”已被该现场证据更新；将来加入诊断后的对象仍需另等自然检查，不能沿用原cohort成功证据。

14:18 UTC已实际提交诊断13359766，命令 `sbatch --parsable --dependency=afterok:13359727:13359759 sbatch/M28-A10-CONTEXT-DIAGNOSTIC.sbatch`；前置三项检查实际命令 `sbatch --parsable sbatch/M28-A10-CONTEXT-DIAGNOSTIC-CHECKS.sbatch`，返回13359759。提交前重做五项live调度，核实新目录与按job隔离日志；CPU诊断不申请GPU，沿用private，输出为少量表格／summary，不准备新大缓存。scontrol记录CPU2／8GiB／10分钟、原提交依赖，检查依赖已满足后剩eval有效Dependency；原fit仍RUNNING gpu035，infer／eval有效等待。没有作业取消、延长或补训练。

固定心跳／固定监控对话更新为 M28-FULL-A6-A7-A9-A10-RAW-DIAG-20261003，旧14根不变，加CPU诊断共15根。A10 fit→infer→eval→诊断全部完成或出现真实新故障才交接，单fit／eval成功且下游正常等待时静默；完整诊断不是科学或论文PASS。ENTIRE prompt／ID／name／kind／ACTIVE／15分钟／target／通知策略／created_at精确保存读回，未改模型／推理强度，未创建对话／cron；新对象未发送人工startup，其自然定时运行／完整回传尚未发生。执行记录遵循同对话任务方式的[官方说明](https://learn.chatgpt.com/docs/automations?surface=app)。

## 固定拟合已完成，完整推理已启动（14:20 UTC）

A10 fit13359631的root／batch／extern均COMPLETED/0:0，31分20秒；STATUS COMPLETED和JOBID正确。主对话读回末权重／三遍／4608步／649457参数／同GEN revision／完整context metadata，configuration去除encoder_context后与完成A3的原配置相同；DEV_used_for_selection及test_setaria_used均false。training.jsonl实际4608行，step_001536／003072／004608三份checkpoint各7846736 bytes，均非空，固定拟合loop1865.420338秒。未按中间loss选权重或变更曝光。

原已提交raw推理13359725实际RUNNING gpu035，STATUS RUNNING／JOBID13359725一致；eval13359727与诊断13359766保持有效Dependency。启动并不证明全8690窗完成、7728主参考性能或科学PASS；仍按原规则等待完整输出及正式评价。

## 全流程终态和冻结判断（2026-10-03 15:50 UTC）

fit13359631、raw infer13359725、eval13359727、diagnostic13359766的root／batch／extern全部COMPLETED/0:0；分别31分20秒、1小时15分22秒、29秒、14秒。全部STATUS／JOBID与输出一致，完整8690方向窗、7728主参考和末004608协议完成；上文RUNNING／PENDING仅为当时历史快照。

实际GFF／全7728成员及数值receipt13360847亦COMPLETED/0:0／13秒：5743唯一预测无重复、phase一致，全参考身份／半开坐标／unsupported／既有A3及native字段保留，所有selected成员与实际GFF一致；逐物种／pooled P/R/F1及原gate重算一致。新receipt不重跑推理、优化或新选择。

A10拟南芥／水稻／micro exact-chain F1为56.66094%／58.27528%／57.15982%，均高于matched A3和原固定A；但拟南芥background span84489 bp超过原C0的54809 bp，故冻结gate FAIL。水稻TP较A3少6，不能把F1改善称两物种recall均增加。完整原生ANNEVO／Tiberius／Helixer micro F1分别75.70250%／71.80468%／62.85117%，A10仍低于三者，不具备独立泛化或论文PASS。保留原负结果，不扫阈值、追加epoch、改背景定义或更换主分母救本轮；详见[a10_context_result.md](a10_context_result.md)。

同一监控对话在15:42 UTC自然检查后于15:44:29.300 UTC成功交接A10＋A12，原十根心跳实际PAUSED已读回；主对话消费后仅更新同一心跳为FULL＋A13七根，不重复监控本闭合A10。A13事前固定目标来自A12 TRAIN诊断，不以本轮A10指标后验修改它。
