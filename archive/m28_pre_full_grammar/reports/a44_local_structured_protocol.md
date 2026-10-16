# A44 原局部监督与结构竞争的匹配准备

2026-10-07。[A43 完整 TRAIN 结果](a43_result.md)为 FAIL0/5，模型已关闭。本次原TRAIN零更新预检在一次明确检查接口修复后已通过：全1491窗来源／缓存／gold一致，八窗损失接口与共享参数梯度闭合，模型权重不变。显式1／1 probe下，八窗结构／局部参数梯度范数比为25.27–1391.28，不能将等系数理解成监督作用相同。下一步先固定损失尺度及同容量／同初始化／同draw的局部单独训练对照，再检验额外结构项是否改善完整自然链；本次没有新fit或性能结论。

## 对照问题

两个后续拟议臂均使用 fresh seed0 的原 M28Core C0 表示与十五类发射，366191 参数，不加载 A43 或其他研究 checkpoint，也不增加 query／teacher 任务。保留原六条 TRAIN 染色体、1491 唯一窗口、1536 draws 的顺序与 importance，拟议每臂仍三遍 4608 步及固定末004608，原 AdamW／clip、角色势、null0／threshold0和自然语法不变。局部单独训练只回传原局部项；局部＋结构臂再加入关闭的 A43 有限19189监督配对，不新增 gold、标签或竞争几何。

两臂使用相同原 R4 局部监督，控制表示容量、局部标签信息、draw 和初始化；结构项带来的额外监督及梯度尺度是实验变化，不宣称所有信息曝光或梯度规模都相同。不得把对 A3／A23 的差异解释成纯局部项因果，也不预设普通局部 CE 加结构解码构成新颖性。

未来完整自然 TRAIN 读出仍列较强 A31_A3 的原五条件及2916／1086完整分母、逐物种 known_nonmatch 和单列 unknown；局部＋结构臂还需与局部单独训练对照比较。只有额外结构项在匹配对照上带来实际完整链收益，才支持本机制的进一步研究。混合研究系数及具体继续门槛须在拟合前固定，不能根据训练后效果、DEV 或阈值搜索反推；本次没有选择拟合系数。

## 原监督接线

复用 canonical run_query.py::load_train_inputs，对混合 metadata 先筛 TRAIN，再连接原 R4 窗口和全部 overlapping_transcript_ids。原字段包括裁切的 CDS、完整 coding offset、isoform／gene 分歧、fine_label_unknown_intervals、不可调用 DNA 与 padding；调用原 build_targets，而不是读取缓存旧 R2／R3 标签或将有限完整链以外的位置一律当背景。

在全部1491窗核对局部和关闭有限支持的几何、实际 feature artifact 来源，以及完整正链 exon 几何和2916／1086分母。任一错窗、错缓存或 gold 矛盾即停止；不能修改冻结候选或 R4 记录来通过。链−1仍只从结构分区排除，不自动抹除已解局部区域，也不制造逐碱基负标签。

原局部损失按本窗实际 fine／grouped-region 可监督位置均值；原结构损失按已知候选数归一，两分母不合并。原 window importance 只乘总损失一次。新三项人工来源连接检查仅覆盖本次新增的几何／缓存／完整 gold 连接；原七项局部＋结构数学检查13425410已闭合，不重复运行。

## 实际零更新预检

复用 A43 稳定的八个原 TRAIN 证据窗：每物种最多新增已知、最多新增未知、最多正链 exon、最多正负 span 冲突，按首次 ordinal 解 tie。预检仍先核对全部1491窗及1536draw，只为梯度接线打开这八份原缓存 GEN＋定向 DNA；不挑训练子集、不重新提取 encoder 特征。

两份同型 fresh seed0 模型先逐 tensor 确认相同初始化。对同一次 FP32 前向分别计算原局部、原结构和显式1／1合并项的梯度，核对五个参数模块有限非零，以及局部未知／N／padding位置的局部发射梯度严格为零。根据下方实际人工复现，joint与两分量的相加检查改在新增损失的发射接口，原FP32 rtol1e−4／atol1e−6不变；三次上游FP32内核重放的参数加和误差只记录，不作为通过标准。局部及结构项的范数比与夹角也只是初始梯度规模诊断，不是自动权重算法或训练收益证据；1／1仅为本次接线 probe，不是已选研究权重。

全程无 optimizer、更新或 checkpoint；最终模型逐 tensor 不变，参数 .grad 不积累。记录实际局部位置分母、结构已知数、unknown排除及真实 mask，不重配已有链标签。任何真实来源矛盾、非有限梯度、mask梯度泄漏或权重改变停止，不能按预期数改数据。

先运行 `sbatch --parsable sbatch/M28-A44-LOCAL-STRUCTURED-SOURCE-CHECKS.sbatch`，沿实际 M28 generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立根 `outputs/M28-A44-LOCAL-STRUCTURED-SOURCE-CHECKS`。闭合后运行 `sbatch --parsable sbatch/M28-A44-LOCAL-STRUCTURED-PREFLIGHT.sbatch 13424089 ACTUAL_SOURCE_CHECKS_JOB`，generanno／private／1×24GiB RTX3090／2CPU16GiB10分钟／no-requeue，独立根 `outputs/M28-A44-LOCAL-STRUCTURED-PREFLIGHT`；唯一 M28A44SOURCECHECKS／M28A44PREFLIGHT_%j 日志。时限沿已有同型八窗预检，增加分量梯度重放仍不自动延长。提交前核对新 AST／SBATCH、明确文件传输、五项活调度／节点／预约及精确输出冲突。

仅明确实现或人工 fixture 错误允许同科学范围同资源另根 R1 一次；真实数据／数值矛盾、第二次失败、TIMEOUT／OOM 均停止，不救参或延时。短准备由主处理；仅外部等待需要监控时复用原同一固定轻量监控和心跳，不创建新的监控对话或 cron。

## 当前状态与边界

三项来源连接检查13425615已于07:22:04至07:22:08 UTC实际闭合，root／batch／extern COMPLETED0:0，三项PASS0.003秒，MaxRSS454584K；STATUS／JOBID、canonical命令／工作目录／日志一致。实际generanno／private／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0。新三AST、两SBATCH及明确六文件local／canonical一致；提交前五项活调度／节点／预约和精确输出冲突已核对。

实际预检13425634在第一证据窗的上游参数梯度加和断言中止，未完成八窗或生成summary，不是模型科学PASS／FAIL。原rtol1e−4／atol1e−6下，首个262144元素参数tensor有60014元素不一致，最大绝对差0.0020580291748046875。运行已到这一步说明此前的全1491窗来源／缓存／gold连接及该窗局部mask检查没有抛错，但这是执行路径推断，不是完整预检证明。没有optimizer、参数更新、checkpoint或新模型拟合；原失败记录保留。

先用独立人工CUDA诊断回答：同型C0网络在原运行浮点路径上，即使合并目标的发射梯度正确，上游FP32参数梯度的三次重放是否仍可能不能满足此逐元素断言？命令 `sbatch --parsable sbatch/M28-A44-GRADIENT-ARITHMETIC.sbatch`，generanno／private／1RTX3090／2CPU16GiB2分钟／no-requeue，独立 `outputs/M28-A44-GRADIENT-ARITHMETIC` 和 M28A44GRADARITH_%j日志。只人工24576位DNA one-hot／4096×2048特征／两条几何／fine背景目标及fresh C0，比较原运行设置与仅该人工进程内关闭TF32的诊断；最终恢复原设置。不打开真实输入或权重、训练、改变研究计算设置或放宽容差。一次诊断，无自动重提或延时；仅在证据明确是检查接口的实现问题时，才可按原一次R1规则修正检查位置，不能将真实数据或损失矛盾解释成实现错误。

人工诊断13425702于07:33:37至07:34:14 UTC实际闭合37秒，root／batch／extern COMPLETED0:0，STATUS／JOBID／命令／工作目录及唯一日志一致，stderr空，MaxRSS1758996K。实际generanno／private／gpu034／1RTX3090／2CPU16GiB2分钟／no-requeue／Restarts0。PyTorch2.5.1原matmulTF32 false、cuDNN TF32 true；原设置和仅人工关闭TF32的两格发射梯度相加误差均0，但上游参数均不能通过原逐元素断言，整体相对L2误差分别8.076982783141805e−5和2.963270565409212e−7。因此原断言会在损失组合正确的人工输入上误报，不适合作为新增损失接线判据；关闭TF32也不是修复本研究的方法。浮点运算不保证数学等价计算逐位相同，见[PyTorch数值说明](https://docs.pytorch.org/docs/2.14/notes/numerical_accuracy.html)。这份一般说明不替代本次实际复现或证明模型科学有效。

唯一R1修正检查接口而非放宽容差：在同一实际logits上用原rtol1e−4／atol1e−6核对局部＋结构损失发射梯度；上游三个目标仍逐模块检查有限非零，并保留参数重放加和误差为诊断。原数据、模型、1／1 probe、loss、mask、分母、importance、FP32／原TF32设置及零更新不变。命令 `sbatch --parsable sbatch/M28-A44-LOCAL-STRUCTURED-PREFLIGHT-R1.sbatch 13424089 13425615`已提交为13425772，原1RTX3090／2CPU16GiB10分钟／no-requeue，独立 `outputs/M28-A44-LOCAL-STRUCTURED-PREFLIGHT-R1`、M28A44PREFLIGHTR1_%j日志。保持原失败根不变，不做R2、延时、调参或fit。

## 实际预检闭合与尺度发现

13425772于07:42:38至07:43:18 UTC运行40秒，root／batch／extern均COMPLETED0:0，STATUS／JOBID／canonical Command／SubmitLine／WorkDir和唯一日志相符，stderr空，batch MaxRSS1893448K。实际generanno／private／gpu034／1RTX3090／2CPU16GiB10分钟／no-requeue／Restarts0；修复源AST、R1 SBATCH语法、明确三文件local／canonical及07:41:44 UTC提交前五项活调度／节点／预约／精确输出冲突核对通过。未重复原已闭合的三项来源或七项数学检查。初次调度核对已见终态，全部短准备由主消费，没有启动新的等待监控；原同一固定心跳仍暂停。

实际summary中的六项证明均true：全1491窗原R4／有限支持几何、缓存及完整gold一致；实际局部mask和位置分母；局部、结构、合并项各自的五模块有限非零梯度；损失发射接口梯度相加；两份同容量fresh seed0初始化相同；最终权重不变。原1536draw／importance和19189支持保留，完整gold／长链分母仍拟南芥2172／854、水稻744／232。三遍replay只是后续拟议曝光，本预检实际optimizer0、没有optimizer构建、checkpoint或研究权重读取。

八窗损失发射梯度相加最大误差均0，局部未知／N／padding发射梯度严格0；原cuDNN TF32 true、matmul TF32 false、matmul_precision highest不变。上游参数重放相加最大误差仍有0.001140594482421875至0.02490234375，明确只作诊断，不能宣称上游逐元素严格线性或独立全网络数值等价。

下表来自实际相同前向的两分量反向，范数比为整个366191参数向量的结构／局部L2范数，各项均乘原同一window importance。可监督位置含fine与已解grouped-region位置。链unknown与局部mask是不同粒度，例如水稻“新增未知最多”窗排除10条未知链，但24576个局部位置均可监督。

| 物种及证据类别 | 局部可监督位置 | grouped-region位置 | 结构已知数 | 排除unknown链 | 结构／局部梯度范数比 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 新增已知最多 | 21941 | 173 | 26 | 3 | 182.63 |
| 拟南芥 新增未知最多 | 23143 | 3 | 17 | 4 | 240.35 |
| 拟南芥 正链exon最多 | 24159 | 92 | 24 | 3 | 142.78 |
| 拟南芥 span冲突最多 | 22625 | 4 | 36 | 0 | 331.21 |
| 水稻 新增已知最多 | 24498 | 2 | 22 | 1 | 112.36 |
| 水稻 新增未知最多 | 24576 | 0 | 7 | 10 | 1391.28 |
| 水稻 正链exon最多 | 24573 | 1 | 19 | 0 | 363.87 |
| 水稻 span冲突最多 | 24573 | 1 | 24 | 1 | 25.27 |

这是八个事前选定接线证据窗的初始梯度诊断，不是全部TRAIN分布、训练后作用或独立性能证据。两分量整体梯度cosine为0.2868–0.8161，也不能据此断言所有窗口／模块不存在冲突。结构角色势按链内位置累加，而局部项按可监督位置均值；仅按已知候选数归一结构项并不会使两项具有相同尺度。直接用1／1研究系数并经总梯度clip，可能使局部贡献被结构项淹没，这是需要处理的尺度问题，不是A43退化原因已获证明。

下一拟合协议优先明确两项的监督单位及位置／候选归一化，固定一次方案及局部单独训练对照，不开展多权重搜索或根据DEV选系数。1／1仍仅是本预检probe；本轮没有选择新的研究normalizer或系数。权威实际产物为canonical `/home/users/j/jwang/coding-rna/outputs/M28-A44-LOCAL-STRUCTURED-PREFLIGHT-R1/summary.json`，初始失败及独立人工诊断根均保留。

原 A36／A42／A43 失败保持关闭。本零更新准备闭合时，新的拟合合同、研究系数和实际模型收益尚未确定；raw／DEV、独立物种、Setaria、A33、对外发布或论文主张不在本准备执行范围。

后续另行固定的[匹配拟合协议](a44_paired_fit_protocol.md)将结构项改为NLL／valid_bases，保留局部位置均值及local-only对照；新两项单位检查13425929闭合后，LOCAL13425946与JOINT13425945已于08:05 UTC开始并行训练。此为下一实验，不追改本页原候选均值1／1 probe、数值结论或初始失败记录，实际模型收益仍待完整自然TRAIN读出。
