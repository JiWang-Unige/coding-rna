# A44 局部监督与结构竞争匹配拟合

2026-10-07。固定一次 local-only 与 local＋structure 对照，检验原局部监督之外的有限结构项能否改善完整自然链。上一轮真实零更新预检已闭合，但八个证据窗在候选均值 1／1 probe 下的结构梯度范数为局部项的25.27–1391.28倍。因此本次先明确损失单位，不采用该 probe 的候选均值直接拟合，也不按八窗梯度选系数。模型收益尚未获得证据。

## 固定损失单位

局部项 Llocal 保持原 R4 fine／grouped-region 损失及实际可监督位置分母 M；M为0时保留原可微零。有限集合结构项 Lstruct＝(logZknown − gold_gain)／N，N为本窗实际 valid_bases，不是 width、已知候选数K、可监督局部位置M或预测gene长度。单位角色 gain 仍对链内角色与intergenic logit差累加，FP64累计不变；只改变训练损失的分母，不改变预测链gain、null0、threshold0或DNA语法。

候选集合仍为闭合13424089的19189原TRAIN窗口—链配对；所有−1链仍从logZ和gold排除，全未知集合为可微零。N不要求局部所有位置均有标签，两类unknown的粒度和原mask保持分开。按有效碱基归一是一次预先固定的单位选择，不是已证实的梯度平衡、最优系数或新颖性结论。

| 拟合臂 | 总损失 | 完整链监督参与反向 |
| --- | --- | --- |
| A44_LOCAL | 原 importance × Llocal | 不参与，结构输入为空且系数0 |
| A44_JOINT | 原 importance × (Llocal＋Lstruct) | 原关闭的19189有限支持，结构系数1 |

两臂局部系数1，原window importance只乘总损失一次；没有温度、长度重评分、梯度动态调权、系数扫描或按DEV选择。原A43损失和已闭合A44 probe不改变，以便保留原结果的含义。

## 数据与匹配条件

两臂各自fresh seed0／366191参数的M28Core C0共享表示和十五类发射，不加载任何研究checkpoint、不添加query／teacher任务。相同原六条TRAIN染色体、1491唯一窗、1536原draw顺序及importance，三遍4608步，原AdamW lr3e−4／weight_decay0.01／betas0.9,0.999／eps1e−8／clip_norm1；FP32网络与原运行浮点设置不变。每遍保存001536／003072／004608，主读出只用末004608，不按loss选epoch。

每臂启动先复核真实闭合13425772的原六项零更新证明、新损失单位检查，以及全1491窗局部R4／有限支持几何、cache来源、完整gold、1536draw和2916／1086完整分母。原load_train_inputs在混合元数据预筛TRAIN，只连接原R4窗口及所需全部overlapping_transcript_ids；原cached定向DNA调用build_targets。旧cache标签不使用，原fine未知、isoform／不同gene分歧、N／padding和coding offset不变。只复用GEN缓存，不执行backbone前向。

LOCAL也复核支持身份，但不将其几何或整链标签用于模型损失。JOINT增加完整链监督，故两臂并非总标签信息或梯度大小相同；本对照控制容量、局部标签、初始状态、优化器及draw曝光，不归因A43退化或A3／A23差异为单一因果。

## 固定自然读出与继续条件

两臂拟合后使用相同末004608发射和原完整语法自然TRAIN读出：预测只读原1491窗口的缓存GEN＋DNA，不接gold、有限训练支持、注入链、旧selected候选或R4标签。原单位角色势、FP64、null0／threshold0、max128、ORF／motif／phase／intron／span规则和稳定tie均不变。完整预测文件闭合后再连接原R4标签，逐物种完整分母拟南芥2172／854、水稻744／232；known_nonmatch和unknown单列。

JOINT必须同时满足两组五条件才继续本结构机制：相对较强冻结A31_A3，以及相对本轮LOCAL，均要求两个物种完整正确链分别严格增加、两个物种已知非匹配分别不增加、pooled至少5-exon完整正确链严格增加。任一条件不通过则本次机制FAIL关闭，不追加seed／epoch／系数／阈值搜索。LOCAL相对A31_A3的同五项另列，区分恢复局部监督与额外结构收益；它通过不等于结构机制通过。

上述均为原TRAIN窗口—链配对的机制继续条件，不是独立测试、去重gene、raw assembly优势或论文claim。raw／DEV、独立物种、Setaria、A33及对外发布仍按独立实际协议放行，不由本拟合PASS自动开放。

## 执行与失败处理

新增两个人工检查仅覆盖本次分母／公式／梯度和零已知结构项：独立枚举empty／gold／重叠nonmatch集合，确认N而非K／width，原importance与局部fine／grouped-region梯度；空及全未知池保持精确local控制。命令 `sbatch --parsable sbatch/M28-A44-LOCAL-STRUCTURED-POSITION-CHECKS.sbatch`，generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立 `outputs/M28-A44-LOCAL-STRUCTURED-POSITION-CHECKS`，唯一M28A44POSITION_%j日志。不重跑未改的旧七项loss、三项来源或八窗实际预检。

新检查真实闭合后，分别提交 `sbatch --parsable sbatch/M28-A44-PAIRED-FIT.sbatch A44_LOCAL ACTUAL_CHECKS_JOB` 和同脚本 `A44_JOINT ACTUAL_CHECKS_JOB`。每臂generanno／private／1×24GiB RTX3090／2CPU16GiB2小时／no-requeue，独立LOCAL-FIT／JOINT-FIT根和唯一M28A44FIT_%j日志，空根才创建。原同型A43实际39:42、原A36实际1:15:07；两小时为原局部目标重建及有限支持的操作余量，不保证完成，不自动延长。每臂三个模型／AdamW checkpoint和短日志预计小于64MiB，不产生新features／logits缓存。

提交前核对实际源、AST／SBATCH、明确文件canonical一致、存储、五项活调度／节点／预约及精确根／日志冲突；两臂在资源可用时并行，各只一次seed0。明确实现错误且尚无optimizer更新可同范围同资源另根R1一次；真实数据／数值矛盾、已有更新、第二次失败或TIMEOUT／OOM停止，不自动重提／续训／救参。保留失败和所有checkpoint，不覆盖原结果、不取消其他作业。

短检查由主消费。拟合或后续完整自然读出仅需等待时复用同一固定轻量监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同heartbeat coding-rna-r5-light-monitor，gpt-5.6-luna／low／15分钟。正常RUNNING／PENDING／依赖静默，当前整批终态或可行动故障一次回传原主并暂停；不建立重复goal轮询或新对话。监控只读结构状态及精确完成产物，不解释loss／性能，不提交后续或改变协议。

## 实际检查与拟合提交

新单位检查13425929于2026-10-07 08:03:55至08:04:25 UTC闭合30秒，root／batch／extern COMPLETED0:0，两项PASS0.677秒，训练入口CLI导入也通过；STATUS／JOBID／canonical Command／SubmitLine／WorkDir及唯一日志匹配。实际generanno／private／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS753500K。新三AST、两SBATCH、明确六文件canonical逐字，以及08:02:14 UTC五项活调度／节点／预约／空根检查完成；没有新数据、模型fit或optimizer参加该检查。

两臂首次各提交一次：A44_JOINT为13425945，08:05:19 UTC提交、08:05:20开始；A44_LOCAL为13425946，08:05:20提交、08:05:22开始。命令分别为同一固定脚本的A44_JOINT／A44_LOCAL参数和实际检查job13425929；提交前08:04:49 UTC五项调度、预约、两个精确空根与存储核对通过，BeeGFS可用156T，无实际硬约束阻塞。

08:06:20 UTC实际squeue／sacct确认两臂并行RUNNING／gpu034，各root／batch／extern RUNNING0:0、STATUS／JOBID／Command／SubmitLine／WorkDir及唯一日志相符；实际每臂1RTX3090／2CPU16GiB2小时／no-requeue／Restarts0。启动日志随后均写到64／4608更新，内部计时分别LOCAL136.32907054305542、JOINT136.34986024396494秒，stderr均0字节。该段只是实际训练已进入更新的证据，不外推全程速度、不选模型或改变已固定时限；原预期运行时间并非保证。

原同一固定heartbeat已原生更新为本批两个fit，实际读回ACTIVE／kindheartbeat／固定target／15分钟／原failed_runs_only及完整3992字符prompt逐字一致，created_at1790981930511保持、updated_at1791360487275。manual turn01a11568-5355-7ea2-9084-18376c3afcce于08:08:41.098 UTC的实际turn_context为gpt-5.6-luna／low／原cwd；08:09:31 UTC最小squeue与sacct均见两臂RUNNING／gpu034／0:0，未读日志／产物或交接。配置与一次实际手动检查已确认，本批定时执行及完整自动终态回传尚未验证。主不形成重复等待循环，应用关闭／休眠／连接失败可能延迟回传。

两臂最终4608步、三个checkpoint及完整自然TRAIN读出／数值核验仍未闭合。提交、RUNNING或初段训练loss都不替代完整链门槛，A36／A42／A43 FAIL保持关闭。

## 配对末模型读出执行合同

本节只实现前述固定读出，不改变正在拟合的源、损失、参数、时限或停止条件。先执行 `M28-A44-PAIRED-READOUT-CHECKS.sbatch`：generanno／private／0GPU2CPU8GiB2分钟／no-requeue，六个人工用例检查两臂身份、末步限制、64任务覆盖、R4 unknown／完整分母、三组五条件及损失单位／importance／mask记录；随后只核对原1491窗、1536draw、保存的较强baseline和2916／1086分母元数据，不读取实际拟合checkpoint、feature tensor或TRAIN catalog。它不是拟合结果或科研PASS。

该检查闭合后，`M28-A44-PAIRED-TRAIN-POTENTIALS.sbatch` 以两个原fit的afterok依赖提交，运行入口仍检查两臂实际终态、STATUS／JOBID／summary和末004608身份。generanno／private／1×24GiB RTX3090／2CPU16GiB30分钟／no-requeue；每原窗只加载一次原GEN＋定向DNA，两臂各作冻结发射前向，输出两套1491个有效碱基×15的FP32势，不打开候选池、标签或catalog，不执行backbone。2982份势预计约4.4GB，预留至少5GiB并在提交前核对实际存储；不产生新feature。权重不变、参数冻结及无梯度按实际运行检查。

`M28-A44-PAIRED-TRAIN-DECODE.sbatch` afterok势job，固定array0–63%64：LOCAL任务0–31、JOINT任务32–63，每臂按原ordinal modulo32覆盖原1491窗。每任务generanno／private／0GPU2CPU8GiB2小时15分钟／no-requeue；只用对应势与DNA执行原完整语法，不走beam或旧编译/JIT路线。各精确task根保留独立STATUS／JOBID、自然预测和唯一array日志；提交前累计作业不得超过80、array不得超过100，预约窗口需满足运行。

`M28-A44-PAIRED-TRAIN-READOUT.sbatch` afterok整个64任务array，generanno／private／0GPU2CPU8GiB15分钟。先闭合两臂各1491条完整预测文件，再打开原池／R4标签／catalog作实际标签连接及三组五条件统计。`M28-A44-PAIRED-RESULT-CHECK.sbatch` afterok读出，同环境／0GPU2CPU8GiB15分钟，核对两臂原draw／importance／N与M／4608步、六份有限checkpoint与AdamW曝光、全部2982势及64任务身份，重新计算DNA语法、独立全窗FP64 log-probability角色分数、原R4 unknown与完整分母、逐物种统计及固定五条件。相同模型对原自然池重评分只用于完整语法目标的数值下界，不是预测候选输入、性能比较基线或改阈值。该CPU核验不宣称独立神经前向或完整DP重跑。

上述新stage只一次执行，所有根须空且日志唯一；已闭合原A43单臂势9:55、decode约1:28–1:44、标签读出5:52及CPU核验4:18仅作为操作时限参考，不保证配对运行速度。明确读出实现错误且没有新训练更新时，允许同科学范围、同资源、独立R1根一次最小修复；真实数据／数值矛盾、第二次失败或TIMEOUT／OOM停止，不延长、不改变冻结参数或以科学失败作工程重试。未成功依赖不执行下游；保留失败根和日志，不自动取消任何作业。

完成范围始终是同一原TRAIN机制检验。只有两个五条件同时通过且实际数值核验闭合，才达到本机制的继续条件；不等于严格raw-FASTA新数据性能、独立测试、较SOTA优势、创新性或论文完成。仍复用原唯一固定监控，实际提交下游后更新同一批对象；不另建监控或让监控解释数值。

## 实际配对读出准备与提交

新读出检查13426144于2026-10-07 08:35:25至08:35:59 UTC真实闭合34秒，root／batch／extern均COMPLETED0:0；六项人工检查PASS0.061秒，原预测输入／较强baseline／2916与1086完整分母元数据核对通过。STATUS／JOBID、canonical Command／SubmitLine／WorkDir和唯一M28A44READCHECKS_13426144日志相符；实际generanno／private／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0、batch MaxRSS767276K。没有打开拟合模型、feature tensor或TRAIN catalog，也没有新fit、搜索或标签计算。三份新Python AST、五份SBATCH bash语法、明确八份实现及本协议的canonical逐字一致已核对；既有四份依赖未改。

08:36:56 UTC五项活调度／预约、四个精确空根与日志、存储检查通过，BeeGFS可用156T，实际当前8个账户作业加64任务及另外三个stage为75，未超过80。四下游于08:37:30 UTC一次提交，均原generanno／private／no-requeue；GPU势1RTX3090／2CPU16GiB30分钟，64个CPU解码每个2CPU8GiB2小时15分钟，读出和核验各0GPU2CPU8GiB15分钟。

| 阶段 | 实际job | 实际提交参数 |
| --- | --- | --- |
| 两臂势 | 13426150 | `--dependency=afterok:13425946:13425945 --export=ALL,CHECKS_JOB=13426144 sbatch/M28-A44-PAIRED-TRAIN-POTENTIALS.sbatch` |
| 两臂完整解码 | 13426151_0至63 | `--dependency=afterok:13426150 --export=ALL,POTENTIAL_JOB=13426150 sbatch/M28-A44-PAIRED-TRAIN-DECODE.sbatch` |
| 原R4标签读出 | 13426152 | `--dependency=afterok:13426151 --export=ALL,POTENTIAL_JOB=13426150,DECODE_JOB=13426151 sbatch/M28-A44-PAIRED-TRAIN-READOUT.sbatch` |
| 实际数值核验 | 13426153 | `--dependency=afterok:13426152 --export=ALL,POTENTIAL_JOB=13426150,DECODE_JOB=13426151,READOUT_JOB=13426152 sbatch/M28-A44-PAIRED-RESULT-CHECK.sbatch` |

以上命令前缀均为 `sbatch --parsable`。08:38:28 UTC主实际确认两fit RUNNING／gpu034／0:0，四下游均PENDING普通Dependency；实际Command／SubmitLine／WorkDir／exports／资源／no-requeue／Restarts0和唯一日志匹配。array0–63／throttle64，读出依赖实际为afterok:13426151_*，即整个array而非单任务；账户展开实际75个作业。该状态只证明依赖链已提交，不证明任何下游已运行或取得结果。

原两fit监控同一批扩为69个root（两fit＋势＋64decode＋读出＋核验），仅在整链实际终态或可行动故障回传。原fit-only的08:23:19.807 UTC检查实际为系统heartbeat turn01a11575-bbbe-7571-9223-52b9d8f65257，trigger current_time_iso08:23:19.740Z，gpt-5.6-luna／low，completed／errornull／28332ms／cursor139；journal原标manual已发现需更正。这是fit-only定时检查证据，不是新完整链定时检查或终态回传证据。

同一heartbeat原生更新实际读回ACTIVE／kindheartbeat／原固定target／15分钟／原failed_runs_only，完整9637字符prompt逐字相符；created_at1790981930511保持、updated_at1791362434191。第一条本机读回命令遇到stdin编码SyntaxError，未读取状态；改为ASCII编码标准库命令后exit0核对成功，不属于研究作业或自动化失败。新完整链一次手动启动检查已发给同一监控，不覆盖模型或推理设置；其执行、定时运行及完整自动回传分别待证。

新完整链manual turn01a11588-13a6-7d83-a85e-1b5bb7437699已实际completed／errornull／38812ms／cursor141；实际08:43:21.943 UTC turn_context仍为gpt-5.6-luna／low／原cwd。08:43:46 UTC最小squeue／sacct见两fit仍RUNNING／gpu034／0:0，四下游均PENDING普通Dependency／0:0；未读取产物、日志或科学值，未交接／重提。同一完整链配置和一次手动检查已验证，定时执行及完整自动终态回传尚未验证；仍无完整基因性能结论，论文总目标未达成。

## 配对 raw FASTA 工程准备

现有raw接口写死A43来源、原单组五条件及输出身份，不能用于本轮两臂。新增 `infer_local_structured_emission.py` 在模型、原生GEN、FASTA和输出创建之前，要求实际69根终态及本轮两fit、完整读出和CPU七证明，重新核对原较强基线／逐物种完整分母及两组五条件。它仍要求另行固定 `a44_paired_raw_protocol.md` 和显式FASTA；本节不是该真实执行合同，没有创建它或放行DEV、独立物种、Setaria、A33及公开发布。

推理准备按原双方向窗口逐窗仅编码一次，两臂共享同一GEN和DNA对象，各自产生FP32十五角色势。使用原nearest-containing-window-center逐碱基上下文和tie、原完整语法／FP64／null0／max128；两臂单独解码、保存身份与GFF3，不混合为ensemble，不载入A3／A43生成器、query、候选池、推理标签或feature cache。逐碱基上下文与整contig可允许的span不同于TRAIN窗口目标，不能用TRAIN门槛证明这种原生域的准确率、资源可用性或公平运行时间。

六项人工检查覆盖两臂真实schema／CPU证明／两个五条件、失败条件或缺合同在读取前停止、两份独立末checkpoint身份、非恒定角色势在双方向多窗口拼接与原实现逐值一致且编码对象共享、负链非空CDS／phase／GFF与臂身份、空或重复流不写成功以及真实零预测保留两臂收据。检查只用人工模型／DNA／receipt，不读实际checkpoint、GEN、catalog、参考或DEV，不进行optimizer。原A43接口及正在运行的拟合／读出源码不改。

唯一短入口 `sbatch --parsable sbatch/M28-A44-PAIRED-RAW-CHECKS.sbatch`，generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立outputs/M28-A44-PAIRED-RAW-CHECKS、唯一M28A44RAWPAIR_%j日志。源／AST／SBATCH及canonical明确文件核对、五项活调度／预约／精确根日志检查后提交；人工输出极小，不产生研究缓存。明确实现或fixture错误允许同资源、同人工范围、独立R1一次；第二次失败、数值／坐标矛盾、TIMEOUT／OOM停止。不重测未改旧语法，不改当前69根监控对象，不以此工程PASS替代完整TRAIN条件或真实原生结果。

该唯一人工检查13426251实际于08:59:39至08:59:55 UTC运行16秒，root／batch／extern均COMPLETED0:0、STATUS／JOBID及Command／SubmitLine／WorkDir和唯一日志匹配；实际generanno／private／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS745328K。六项PASS0.229秒，包括双方向非恒定十五角色拼接逐值相同、GEN／DNA对象共享、两臂负链非空phase和GFF。08:57:00 UTC五项活调度／预约与精确空根检查通过，展开账户75作业加本短检查为76<80；三份明确新文件及本协议已canonical逐字核对。新fixture的transcript身份按实际export字段id填写，修正发生于首次提交前，不存在失败或R1消耗。没有实际GEN／模型／FASTA／参考／DEV执行。

## 配对完整评价工程准备

`evaluate_local_structured_emission.py` 只提供helpers，没有真实raw／DEV执行入口。新schema保留两臂模型、来源、末checkpoint、各方向原窗口网格和同一contig顺序；两套完整保存预测及各自GFF坐标／phase／count先核对完，随后才打开显式FASTA和一次完整参考context，最后读取比较收据。读取顺序防止用单臂完整、另一臂残缺或错模型输出评分。

复用未改的完整common ruler算术：两臂均保持原7728完整primary参考及相同背景bp，unsupported／长链／非标准splice等不从分母删除；两臂分别报告原生ANNEVO／Tiberius／Helixer及A31_A3／A31_A23对照、micro／macro／worst-species和背景覆盖，另列JOINT对LOCAL。旧A43来源不接受，两臂不ensemble，不把原池owner覆盖指标强加给全contig解码。旧通用算术产生的GFF路径metadata改为已实际核对的对应arm文件，不重标模型身份或改评价结果。

五项新人工检查覆盖两臂完整receipt／来源／checkpoint／网格、负链phase与gain对应和跨窗长链、任一臂错误在参考入口前停止、两臂完整7728人工分母／背景／逐物种worst值和七方法表、显式FASTA完整性与DNA不一致不以过滤真值救回。复用的7728人工参考fixture不是原DEV数据。新短命令 `sbatch --parsable sbatch/M28-A44-PAIRED-EVAL-CHECKS.sbatch`，generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立M28-A44-PAIRED-EVAL-CHECKS及唯一M28A44EVALPAIR_%j日志；同样一次实施／fixtureR1边界，不改69根监控、不重跑未改旧九项评价检查。

这两步不执行实际raw模型或科学评价，不产生SOTA优势、独立泛化、新颖性、因果匹配总监督、独立神经／发射评分／完整DP回放或公平速度证据。后续真实执行仍须先消费完整TRAIN条件、固定实际FASTA／reference access／原生资源和评价合同；任一门槛不通过按原规则关闭，而不是用本准备放行。

当前完整69根的系统定时执行已另行验证：turn01a11591-33c3-7463-bd15-8f46b739779e实际completed／errornull／22501ms／cursor142，08:53:19.992 UTC turn_context仍为gpt-5.6-luna／low；系统heartbeat trigger current_time_iso08:53:19.937Z属于本完整cohort。最小08:53:19 UTC状态仍是两fit RUNNING／gpu034／0:0、四下游普通Dependency／0:0，没有读日志或产物、没有交接。配置、manual和本轮scheduled均有证据，完整自动终态回传仍未发生；两个人工准备检查由主消费，不加入69根对象。

配对评价人工检查13426288实际于09:07:37至09:07:53 UTC运行16秒，root／batch／extern均COMPLETED0:0，五项PASS0.364秒；STATUS／JOBID、实际Command／SubmitLine／WorkDir和唯一M28A44EVALPAIR_13426288日志匹配。实际generanno／private／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS726460K。09:05:54 UTC五项活调度／预约／精确空根检查完成，展开账户74个作业加此短检查为75<80；新两Python AST／SBATCH语法、明确三份实现及本协议canonical逐字相符。原已消费13426251不重跑、未用作真实raw放行；当前拟合／读出源码和69根监控对象均未改。

原 `scripts/eval_structure_diagnostic.py:51` 的unclosed-file ResourceWarning在人工GFF读取中再次出现，stderr保留；五项新检查、两臂负链phase／GFF、完整人工7728分母／背景和七方法表均通过。该既有句柄警告未被隐藏或当成模型FAIL，当前有限两个GFF路径也不据此扩大旧解析器改动。本批没有实际研究模型／GEN／DNA／catalog／参考／DEV读取或optimizer、新label、真实SOTA评价和公开发布，完整科研目标未完成。

## 实际完整读出超时

2026-10-07原固定监控11:23:20 UTC报告两fit、势与全部64个decode任务已COMPLETED0:0，读出运行，数值核验及机制诊断等待依赖；这仅是调度阶段完成，尚未消费完整七项数值证明或科研门槛。

主在新人工检查提交所需的11:32:07 UTC调度核验中见读出已到时限，11:33:23 UTC针对该新疑点确认13426152于11:17:00–11:32:13 UTC终态TIMEOUT、scontrol Reason TimeLimit，root elapsed15分13秒含终止阶段，固定时限15分钟未延长。sacct root ExitCode0:0、scontrol0:15，batch CANCELLED0:15／elapsed15分14秒／MaxRSS845328K，extern COMPLETED0:0；root的0:0不是成功。实际private／generanno／gpu034／0GPU2CPU8GiB／no-requeue／Restarts0、原Command／SubmitLine／WorkDir及M28A44READ_13426152唯一日志相符，stderr明确时限取消，没有OOM或Python异常证据。

最后stdout进度为1280／1491；两臂A44_LOCAL_natural_scores.jsonl、A44_JOINT_natural_scores.jsonl各1491行已写入。11:36 UTC文本核验时windows.jsonl暂存1224行、summary.json不存在、STATUS仍为旧RUNNING、JOBID13426152相符。进度日志与当时落盘行数分别记录，不凭任一部分记录外推完整分母或补写缺失记录；调度终态优先。源代码在完成所有1491原窗的标签连接后才聚合／写summary，本次没有最终三组五条件，不能把未执行的门槛记成FAIL或PASS，也不能从部分TRAIN记录选择模型。

13426153实际PENDING／DependencyNeverSatisfied，13426529实际PENDING／Dependency；后者的上游已无法成功，不能把尚未运行的数值核验或机制诊断当作已失败科学检验。原checkpoint、完整解码、暂存读出、失败日志和依赖作业保留；没有取消、延长、重提、重新拟合或改损失／阈值／R4未知规则。按本协议TIMEOUT停止该执行链，不用人工数学准备放行raw、DEV或独立物种。

固定70根heartbeat于11:38:21.074 UTC实际触发，监控turn01a11628-4814-7bc0-9eab-a13739b164be用最小调度与stderr证据，将同一故障一次回传原主，send_message_to_thread实际完成、主已收到。该监控turn于11:43:25 UTC完成、errornull，11:46 UTC主核验同一automation已保存PAUSED，kind仍heartbeat、固定target不变、15分钟间隔与failed_runs_only不变、当前cohort journal Handoff sent YES。这验证了本cohort定时可行动故障回传及交接后暂停，不是全部70根／210行或结构白名单完成；监控没有读取科研值、修改作业或作PASS／FAIL解释，仍保留同一监控对话供下一轮复用。平台研究goal未暂停或改变。
