# 联合编码三个固定候选的 DEV 与基线评价准备

2026-10-12 UTC。B/C/D完整TRAIN及冻结六方法DEV R1整链已完成并由主处理，不做性能预筛。Fresh×joint／uniform／codon pooled完整链F1为70.106738%／69.794721%／64.578529%；三者两物种均高于Helixer、低于Tiberius和ANNEVO，D相对B仍为负结果。[完整结果](full_grammar_joint_controls_dev_result.md)保留全部候选、计数、背景和差值。初始前置字段错误失败及一次schema-only恢复均保留，科学配置和旧八／十方法家族不改。这是历史复用DEV，不是独立测试或论文goal完成。以下保留冻结设计及实际执行历史。

## 三个候选与匹配控制

| 固定分支 | 源fit与最终权重 | 主要比较 |
| --- | --- | --- |
| FRESH_JOINT24576 | 13494895，outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FIT/step_009216.pt | 三个原生基线；另有已固定TRAIN四格机制比较 |
| UNIFORM_COARSE | 13496011，outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-FIT/step_009216.pt | fresh×joint；三原生基线为描述性参照 |
| CODON_ENERGY | 13497914，outputs/M28-FULL-GRAMMAR-CODON-FIT/step_009216.pt | fresh×joint；三原生基线为描述性参照 |

三者从同一原fit13434603的final4608头及保存AdamW出发，同fresh seed1／4608有序draw／4219 unique、joint24576、原有限观察与ownership、十五角色、cap128、FP64、最后9216。uniform只在追加阶段把coarse信息置为均匀通道；codon保留epsilon0.1并加入64个零初始化完整拼接CDS密码子参数。codon的新参数AdamW年龄为4608、原头为9216，联合裁剪及优化轨迹属于该干预，不能宣称与增加参数完全分离的单因素因果效应。

6144／7680只是训练快照，不进入DEV选择。DEV不用于checkpoint、分支、阈值、head、上下文、抽样、loss、epsilon、训练时长或下一轮训练选择。负结果、物种差异和未完成分支完整保留；真实执行失败不是零准确率，不通过删除失败分支补成成功家族。

## 执行前置

新研究DEV只在三条原六TRAIN流程均完成末权重核验、六染色体严格raw推理、共同尺评价和四格／配对比较，并由主对话处理结果后执行。完成条件只判断身份与完整性，不设置性能门槛。入口核对下列真实来源及全部30个numeric Slurm roots的root／batch／extern COMPLETED0:0；未来数组任务的numeric roots从各真实JOBID收据及Slurm取得，不猜测父数组加task。

| 分支 | fit | verification | 六任务raw数组 | evaluation | 完整TRAIN比较 |
| --- | --- | --- | --- | --- | --- |
| fresh×joint | 13494895 | 13494896 | 13494897 | 13494898 | 四格13494899 |
| uniform | 13496011 | 13496012 | 13496013 | 13496014 | 配对13496015 |
| codon | 13497914 | 13497943 | 13497945 | 13497946 | 配对13497947 |

入口复用各分支已实现的source身份核验，再核对全部六条raw外收据、最后权重、评价范围、配对源作业与输入／计数一致性标记。某条未完成或失败就停止在前置核验，不打开DEV FASTA、加载新checkpoint或读取参考；常规执行修复须保留失败输出，不能改变科学配置或放宽身份条件。没有用新签名文件或活动goal推定人类授权。

## 固定输入与共同评价尺

范围仅既有两条DEV染色体：拟南芥NC_003074.8／23459830bp、水稻NC_089041.1／29936421bp。保留原大写IUPAC FASTA记录供直接比较，caller仍把非ACGT转N，不改长度或坐标。三候选都直接读取完整FASTA、保留双链和原网格／nearest-center ownership，不给RNA、蛋白、softmask、参考链、训练观察或缓存特征。

共同尺与[原冻结DEV方案](full_grammar_dev_baseline_plan.md)相同：M26 CDS-assessable longest-eligible-CDS／transcript-ID tie，AT5437／Rice2291／总7728 primary参考链。unsupported、长span和重叠参考保留，不以训练eligible或模型语法重选分母。预测按原生gene取最长CDS、transcript-ID tie、complete_only=False；不筛除基线的非标准、partial或不符合M28语法的预测。primary exact-chain包含seqid、strand及全部半开CDS边界，phase和其他assessable isoform另报，后者不算primary TP。

共同评价在读取参考前核对候选与三基线运行身份、全部成功Slurm步骤、原生版本／操作点、完整GFF和候选／原生／原FASTA直接相等。codon分支先严格加载其真实末权重，再比较raw报告的64原参数和去均值读出、完整拼接CDS评分身份、全部窗口与两个方向。原DEV评价器增加显式candidate_validator接入；不传该参数时原身份核验、默认调用形状与评价行为不变。仅新codon接线使用不同评分身份，共同参考选择和计数评价不改。

ANNEVO、Tiberius和Helixer只复用原六份成功预测及固定资产，不更新、不重新训练或重跑，也不宣称未经查询的最新版本。ANNEVO固定42c920f／2.3.3／Magnoliopsida／batch32；Tiberius固定871d497的实际2.0.8隔离包与原2.0.7基础SIF／angiosperms／batch4／seq200070／无softmask或ClaMSA；Helixer固定d17bb49／0.3.7／land_plant_v0.3_a_0080.h5／batch16／subsequence64152及原后处理。实际完整版本、权重路径、命令和收据按原评价器核对，不用缓存成绩表替代真实输出。

## 完整报告与结论范围

每候选报告两物种及pooled完整链TP／预测数／参考数与P／R／F1、CDS-base union、背景span／CDS覆盖与纯背景链／Mb、phase、其他isoform、重复链及描述性结构层。三份比较复用同一基线测量，六方法表要求这些基线逐物种及pooled完整字段直接相等，再从计数复算两端点、背景、micro、macro与worst。worst先分别取各臂最差物种再作差，不能用最小物种差替代。

三候选均与三固定基线给出物种／pooled／macro／worst完整链F1百分点差；uniform和codon的主要匹配控制均为fresh×joint。描述性全表不自动产生赢家、科研PASS或新训练决定。不将base收益代替完整CDS链收益，不把raw FASTA输入称作各方法未见物种，也不把完整CDS链评价扩展为UTR、所有isoform或全部gene类型能力。

两个历史复用开发染色体、单fit、方法相关或未知预训练annotation暴露，不能支持独立泛化或确认性gene级置信区间。新M28推理计时与旧原生工具收据的起止点未统一，不作公平速度优势主张。sealed／held-out／Setaria、新物种／新面板和公开发布不在本次范围内，论文goal保持未完成。

## 实现与固定资源

新增infer_full_grammar_joint_controls_dev.py、evaluate_full_grammar_joint_controls_dev.py和aggregate_full_grammar_joint_controls_dev.py。初始研究输出根为outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-RAW/<branch>/<seqid>、DEV-COMPARISON/<branch>和DEV-FAMILY；R1对应RAW-R1、COMPARISON-R1和FAMILY-R1，均使用完整M28-FULL-GRAMMAR-JOINT-CONTROLS前缀，不覆写原失败结果。每份raw保留STATUS／JOBID、started、原input.fasta、外summary和prediction/{summary.json,full_grammar.gff3}；每比较及家族保留STATUS／JOBID与summary，家族另有comparison.md。

三个研究SBATCH使用以下固定资源：RAW每分支两任务、最多并行2，每项单24GiB RTX3090／2CPU／16GiB／2小时；比较无GPU／2CPU／8GiB／20分钟，必须依赖该分支整个数组；六方法家族无GPU／2CPU／4GiB／5分钟，必须等待全部三份比较。原generanno环境和离线模型资产不变，两次实际提交前均核对调度、维护和存储。每条分支的实际推理耗时仍待观测，2小时是资源上限而非速度保证。

人工检查独立输出在outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-CHECKS，无GPU／2CPU／4GiB／5分钟。六组新人工检查覆盖三TRAIN关闭和真实步骤失败阻止DEV、三入口保持IUPAC与固定权重、codon参数／评分／覆盖错配、评价hook在DNA／参考前拒绝错误、六方法完整计数及worst交换、写出新表而不覆盖旧表；同时重跑原sampling-context五组默认接线检查，验证这次显式hook不改变已有调用。测试只用人工临时文件、短序列和mock，不读研究DNA／参考／历史DEV收据或预测，不加载研究checkpoint或运行native GLM。人工通过仅证明接线，不证明新DEV推理或准确率。

当前唯一原心跳只等待独立M TRAIN整链；已由主消费的DEV R1、B/C/D、短检查和初始失败DEV均关闭、不重复交接。复用原监控对话和同一自动化，不新建对话／cron／重复心跳。常规前置满足后的实施与验证自主推进，不重复请求已授权范围的批准。

## 人工检查实际执行

本节以下为2026-10-11的原11项人工检查历史；它们通过不代表2026-10-12真实研究DEV成功。此次接入修复的新CPU检查只重跑六组受影响接口，并新增实际已完成TRAIN收据／90个root/batch/extern终态关闭核验；不重跑未修改的五组默认接线，也不打开DEV/DNA/reference/checkpoint或native模型。

## 首次研究提交与前置拒绝

2026-10-12 09:49:34UTC，own expanded作业35，新增六raw+三比较+family后45，低于80；private gpu034／035共16张3090、已分配1+2，匹配GPU即时可用，scratch145TiB可用，只有11月6日CPU维护预留。按固定资源提交FRESH raw13519832→compare13519834、UNIFORM raw13519836→compare13519837、CODON raw13519838→compare13519839，family13519840依赖全部三比较；原生基线未重跑。

六raw全部于09:50:24UTC失败：实际roots13519835/13519832/13519841/13519836/13519842/13519838的root／batchFAILED1:0、externCOMPLETED0:0，分别48–49秒。统一stderr为Complete fixed TRAIN comparison required before DEV；每失败目录仅有JOBID／FAILED STATUS，没有started.json、input.fasta、prediction，拒绝在DEV读取和checkpoint加载之前。三比较及family被失败依赖阻止；原作业和失败记录全部保留，不取消、不覆写，不把它们记为零准确率。

原因是infer_full_grammar_joint_controls_dev.py将B四格source_jobs严格对照为F00/F10/F01/F11，但真实producer aggregate_full_grammar_fresh_joint24.py与13494899结果键名为F00_original_split/F10_fresh_split/F01_original_joint/F11_fresh_joint。修复只按这四个真实完整名称及相同jobID匹配；不增加短名兼容／fallback，不放宽检查或改TRAIN结果。人工fixture同步真实schema，并新增错名称／错F11 job在DNA／模型前拒绝的回归。模型参数、候选、最终9216、输入／参考、评分、baseline、阈值和对外范围全部不变。

为保留失败产物，恢复输出为RAW-R1、COMPARISON-R1、FAMILY-R1三个独立M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV根；检查为CHECKS-R1。这只是独立失败恢复命名，不重建科研框架。CPU2／4GiB／5分钟确认六组接口和实际TRAIN closure后，才在同资源／科学命令下做唯一一次研究重试。旧11项检查及原失败六任务不改写为成功；实际R1 handles与结果分别后续记录。

### 2026-10-11原人工检查历史

CPU作业13498761由sbatch --parsable sbatch/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-CHECKS.sbatch提交。实际UTC17:43:29至17:43:50，21秒，private-teodoro-gpu／2CPU／4GiB／无GPU／5分钟上限、Requeue0；root／batch／extern均COMPLETED、ExitCode0:0，batch MaxRSS808948KiB。

全部11组检查通过、0 failure／error，包括六组新接线与五组原sampling-context默认回归。unittest报告0.290秒，summary记录整个检查脚本0.293509931秒；两者不是模型运行速度。STATUS／JOBID／summary均完整一致，summary明确research_DNA_or_labels_opened、historical_DEV_receipts_or_predictions_opened、checkpoint_loaded、actual_native_GLM_forward、DEV_used及sealed_data_used全部false，optimizer_steps_this_run=0。

五份Python AST和四个SBATCH语法通过。该检查完成后，三个研究输出根RAW／COMPARISON／FAMILY仍不存在，三个研究SBATCH未提交；名称中的DEV不代表研究DEV已被读取。原B+C+D TRAIN长链没有变更参数、代码、资源或依赖，仍等待实际末权重与全范围读数。原监控新增D范围的真实定时追加已于UTC17:42:46观察，17:42:58检查三个fit均RUNNING、下游有效Dependency，正常等待没有主对话回传；原模型、low和同一目标对话均经实际turn_context核对。

## R1资格与唯一研究重试实际执行

CPU13519990于2026-10-12 10:00:13UTC FAILED1:0、零秒，尚无CHECKS-R1产物。实际scontrol write batch_script表明提交早于rsync完成，Slurm捕获旧SBATCH、仍指向已存在的原CHECKS目录；test ! -e在运行测试前拒绝碰撞。原13498761结果未覆写，这不是R1测试失败或DEV结果。之后等待精确代码／SBATCH传输完成并核对真实源文件，再提交新资格13520076。

13520076于10:06:54–10:07:36UTC，42秒，root／batch／extern全部COMPLETED0:0。六组受影响接口全部通过、0 failure／error，unittest0.177秒、整个检查4.847795604秒；另外读取实际B/C/D TRAIN收据并确认30roots／90条终态关闭，不打开研究DNA／标签、历史DEV收据／预测、checkpoint或native模型，optimizer_steps=0、DEV／sealed=false。原未改的五组默认接线未重复执行。本次资格由主消费，不作为监控等待对象。

10:14UTC提交前own expanded=35，新增10roots后45<80；private gpu034／035共16张24GiB3090、已分配1+2，无可运行GPU排队阻塞，private可立即启动；scratch145TiB可用，只有11月6日CPU教学预留。同资源、同科学命令、Requeue0，按每分支sbatch --parsable提交RAW.sbatch，然后--dependency=afterok:<整个新raw数组>提交COMPARISON.sbatch，最后--dependency=afterok:13520476:13520478:13520480提交FAMILY.sbatch。没有重提已获受理作业或重跑原生基线。

| 固定分支 | 新raw数组 | task0／task1真实numeric roots | 新比较与成功依赖 |
| --- | --- | --- | --- |
| FRESH_JOINT24576 | 13520473_0–1%2 | 13520474／13520473 | 13520476 afterok13520473 |
| UNIFORM_COARSE | 13520477_0–1%2 | 13520482／13520477 | 13520478 afterok13520477 |
| CODON_ENERGY | 13520479_0–1%2 | 13520483／13520479 | 13520480 afterok13520479 |

家族13520481等待全部三新比较。FRESH于10:14:52UTC、另两分支10:14:56UTC开始；10:15:29UTC六raw全部在gpu034 RUNNING，六started.json均为正确最终9216、源fit／task／species／contig长度、all_three_TRAIN_flows_closed=true、DEV_used=true／sealed=false／reference=false／labels_or_feature_cache=false／optimizer_steps=0，三比较和family依赖均有效。started证明入口错误已修复，不证明完整推理、预测正确或新DEV准确率；必须等六raw、三比较及family完整终态及实际结果。

10:19:11UTC六个真实input.fasta均已生成且STATUS仍RUNNING（AT23753091字节、rice30310640字节，包含FASTA格式开销，不是新的contig长度）；FRESH的AT实际stdout已进入joint24576_encoding正链第512／1909窗口。因此新研究raw确已读取固定DEV并执行原生encoder，而不仅是创建启动收据。完整双链caller和后续评价仍未完成；阶段窗口数不用于准确率或速度优势判断。

同一心跳coding-rna-r5-light-monitor已更新为M/P_R1，原watcher01a0f71f-8845-7432-bd13-bf54c1808fb3/local、main01a0332e-e0a9-7aa2-90f9-159405701833/local、Luna/low、15分钟和通知偏好保留；实际保存的heartbeat类型／目标／ACTIVE／间隔及完整33268字符提示已核对。只监控新R1独立根及原M，初始失败P、13519990／13520076和D等均消费关闭；原失败作业不取消、不覆写，不补发或伪造历史交接。原watcher定时turn01a12fa0-30dc-7d30-a1cb-82ab77b74a94于10:19:49–10:20:07UTC在同一对话追加R1范围，模型／推理经真实turn_context核对仍Luna/low；实际SSH只检查M和新R1，fit与六raw RUNNING、下游有效Dependency，DONT_NOTIFY且未发送主回传。新范围定时追加／单次检查已验证，完整自动终态回传仍待真实事件，不能把配置或单次检查等同于自动流程完成。冻结科学边界和论文未完成判断不变。

## R1完整终态与结果处理

2026-10-12 10:53:22UTC家族13520481完成，六actualraw roots13520474／13520473／13520482／13520477／13520483／13520479、三比较13520476／13520478／13520480及family的30条root／batch／extern全部COMPLETED0:0。三候选原固定最终9216、两条DEV、7728 primary参考及三原生输出均保持；直接输入比较、计数／base／background／macro／worst与全部差值由实际程序和主复算核对。[完整结果](full_grammar_joint_controls_dev_result.md)保留六方法全表、uniform物种权衡及codon负差，不用DEV选择权重或下一轮训练，旧八／十方法家族不改。

原watcher于11:05UTC真实向原main交接P_R1并将同一完整提示心跳暂停，P_R1 handoff=YES；主已消费完整科学结果，更新同一心跳为ACTIVE/M-only并核对原target／15分钟／通知偏好／created_at及完整33966字符提示。终态自动回传已经真实验证，不再作为待验证事项。M原fit仍RUNNING、下游有效依赖，所有数据／科学设置不变；不重复等待／交接P_R1或初始失败P，不创建新对话或自动化。论文goal未完成。

