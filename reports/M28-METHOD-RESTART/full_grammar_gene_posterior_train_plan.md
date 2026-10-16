# 原 TRAIN 整染色体 CDS 后验与排序诊断

2026-10-13。本轮在原六条 TRAIN 染色体、固定 E rank4 与 M dense 最终权重上，执行整染色体条件下的完整 CDS 事件后验，并检验其相对于内部 gain 的排序效用。已有八个 49152bp 原生视图通过分区和事件概率验证，但局部后验不等于整染色体后验。目标是决定这一读出是否值得用于后续模型研究，而非以工程通过代替竞争性注释模型或论文结论。

## 固定输入与概率事件

沿用原 TRAIN 的 AT NC_003070.9、NC_003071.7、NC_003076.8 和 rice NC_089043.1、NC_089044.1、NC_089046.1，长度分别为30427671、19698289、26975502、27474823、23931887、28271460。E 使用 fit13499817，M 使用 fit13510324，均 head9216、native4608；不使用当前新 LR 实验的中间权重或冻结 DEV 选择模型。

GPU 阶段重新从各自已完成 strict raw FASTA 和原 source 权重计算十五角色 FP32 势，使用实际 `manifest.starts`、原24576／12288窗口、相同 native padding 和逐base ownership，包括原尾端窗口。不是三窗口局部训练视图，不拼接局部 posterior。保存双方向完整角色数组；不读取预测或参考 GFF，不训练、不更新权重。

CPU 阶段读取原预测 GFF 的全部唯一完整链作为查询，按双strand恢复 oriented 坐标，使用已验证的 `complete_gene_set_posteriors`。分母仍覆盖该完整 oriented 染色体上的所有合法 complete nonoverlapping gene sets 和空集，实际 carry、GT／GC-AG、min-intron5、cap128不变；两个方向独立。查询链不限制分母候选空间，保存 logP 防止概率下溢。重新计算的每gene gain 必须等于原 GFF 八有效位导出得分允许的舍入误差（relative6e−8／absolute1e−5）；否则停止，不能继续解释成原模型后验。原 GFF、DNA、checkpoint 和冻结结果只读，不输出新预测。

## 排序与参考相对诊断

仅在两模型乘六染色体的十二个后验任务全部完成、source／原预测／原 DNA 收据成立后，评价作业才打开原参考 GFF。原 M26 CDS-assessable longest-eligible 主转录本政策和所有 unsupported strata 保留；逐染色体完整链 predicted／TP／reference 必须精确复现既有 E/M TRAIN 评价。

分别评价 exact-primary 与 exact-any-assessable-isoform 两种参考相对事件。每种事件比较 logP 与内部 gain 的 prediction-pool AUROC 和 average precision，严格合并同分组；它们只评价已输出预测之间的排序，不包含未预测的参考基因，也不是整个注释任务的 detection AP。提供完整同分组 PR 曲线与预先固定保留比例0.25、0.5、0.75、0.9、1.0，不挑最佳阈值，不据本结果改训练或筛写 GFF。主转录本 recall 使用全部保留参考链分母，另列预测内 positive recall；any-isoform不借用主转录本分母。

概率诊断只报告 reference-relative Brier 和固定十个等宽 bin 的 mean structural probability／exact-reference fraction，不拟合温度或校准器。结构包含事件与参考正确事件并不相同；非匹配预测不是已证伪的gene，Brier／bin不能证明生物学校准。按六染色体、两物种及 pooled 报告两模型所有结果，不用 pooled 排名掩盖物种差异。这是同训练范围探索，不能用于独立泛化或新增冻结 DEV 方法。

## 执行资源与验证

新增四组人工检查：minus坐标roundtrip与短exon、同分排序／完整参考分母、log下溢保持排序及单类undefined、固定概率bins／Brier。CPU资格同时读取原 E/M 完整 raw source JSON、十二份原预测 GFF 和已有八案例资格，验证原链数量／唯一性及每mRNA有限得分；不读取研究 DNA、真实参考 GFF、checkpoint 或 native。原通用评价parser有意丢弃score，本读出只恢复原导出mRNA行的得分，不改该parser。只检查新接线与口径，不重复既有四组后验数学测试。

新 GPU 势任务12个，每个单RTX3090／2CPU／16GiB／1h；原最大染色体双strand native编码实测约27分钟，保留原硬件与余量。各自 afterok 一个 CPU后验任务，2CPU／32GiB／4h／无GPU，两方向顺序计算；时限含依据局部实测的规划余量，不是整染色体实测完成保证。全部后验 afterok 一个2CPU／8GiB／30min评价任务。使用 generanno／private-teodoro-gpu／no-requeue，按实际调度核对可用资源；维护和 own expanded80 上限继续有效。

完整势数组共约35.04GiB，metadata／后验／曲线预计小于1GiB，整批声明44GiB上限；新数组单独放入 `outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN/<model>/<seqid>/emissions`，后验放同级 `posterior`，最终 `evaluation`。不改旧 raw 128MiB限额、不重复DNA文件、不保存新weights，现阶段不删除数组。每阶段 STATUS／JOBID／started和summary构成具体完成证据。真实 scientific negative 原样保留，不通过改cap／阈值／范围寻求漂亮成绩。确定性接线错误可至多一次同范围修复，使用独立新输出保留失败；不取消作业、不自动延长时限。

原LR四格批次继续原链。新批次只复用 watcher01a0f71f-8845-7432-bd13-bf54c1808fb3和heartbeat coding-rna-r5-light-monitor，两个完整等待事件分别交接，不以部分encode／posterior成功代表整批完成。正常运行／有效依赖静默；预期末端无squeue行立即查sacct／markers，不能按正常等待处理。配置成功、当前检查成功与实际终态回传分别记录。无新chat／automation／cron，无DEV／sealed／held-out／新panel／对外发布或goal状态变更。

## 实际接线记录

新增三个Python文件AST、四个SBATCH语法和四组本机人工算术检查通过；八份新源／协议按精确叶同步，远程回读全部一致。CPU13548771于03:37:04–03:38:19UTC、75秒，root／batch／extern全部COMPLETED0:0、STATUS／JOBID与summary完整。四组新检查通过，原E/M十二份whole-contig raw JSON／GFF计数、唯一链、有限导出score及旧八案例资格接线通过；程序实耗6.599576秒，没有研究DNA、真实参考、checkpoint或native读取。主已消费，watcher handoff=NOT_SENT_BY_WATCHER；不把资格当作整染色体后验结果。

03:26–03:39UTC调度取证：own expanded约30，新增25正式任务后仍小于80；BeeGFS约161TiB可用，覆盖44GiB声明。private与shared同型3090当前GRES已分配，没有可信更早的shared入口，故保持private；gpu025实际6张而非旧配置8张。private CPU可用内存足以分配32GiB任务，维护预留仅11月6日CPU。本轮不根据nvidia-smi或`sinfo %g`（该字段是允许组，不是GRES使用量）判断空闲。

正式25个作业于03:42:09UTC（北京时间11:42）接受。每对命令为 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-EMISSIONS.sbatch <model> <task> 13548771` 及 `sbatch --parsable --dependency=afterok:<emission-job> sbatch/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-READOUT.sbatch <model> <task> 13548771 <emission-job>`。末端为 `sbatch --parsable --dependency=afterok:13548915:13548917:13548919:13548921:13548923:13548925:13548927:13548929:13548932:13548934:13548937:13548939 sbatch/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-EVAL.sbatch 13548771`，实际evaluation13548940。

| 模型 | task | TRAIN染色体 | GPU势job | CPU后验job |
| --- | ---: | --- | --- | --- |
| E_RANK4 | 0 | NC_003070.9 | 13548914 | 13548915 |
| E_RANK4 | 1 | NC_003071.7 | 13548916 | 13548917 |
| E_RANK4 | 2 | NC_003076.8 | 13548918 | 13548919 |
| E_RANK4 | 3 | NC_089043.1 | 13548920 | 13548921 |
| E_RANK4 | 4 | NC_089044.1 | 13548922 | 13548923 |
| E_RANK4 | 5 | NC_089046.1 | 13548924 | 13548925 |
| M_DENSE | 0 | NC_003070.9 | 13548926 | 13548927 |
| M_DENSE | 1 | NC_003071.7 | 13548928 | 13548929 |
| M_DENSE | 2 | NC_003076.8 | 13548930 | 13548932 |
| M_DENSE | 3 | NC_089043.1 | 13548933 | 13548934 |
| M_DENSE | 4 | NC_089044.1 | 13548935 | 13548937 |
| M_DENSE | 5 | NC_089046.1 | 13548938 | 13548939 |

03:42:39UTC所有25个实际controller handles存在：12GPU势任务PENDING Resources／Priority，12CPU后验及末端PENDING有效Dependency；每项资源、Requeue0、独立log、paired afterok和全部12posterior依赖已核对。M任务2和4的ID有实际跳号，表格取自接受收据，不推算。尚无新整染色体后验或排序结果。

同一原heartbeat已更新为两个独立完整事件 LR13545908／POSTTRAIN13548940，两个handoff均NO，11,446字符完整prompt保留原LR21roots／63记录合同、新25roots／75记录映射和严格结构白名单。先完成的完整事件或需处理失败即回传原main一次并暂停同心跳，另一NO供主消费后恢复；不要求两批同时结束，不让部分encode／posterior成功触发科研结论。完整prompt／ID／heartbeat／ACTIVE／原target／15分钟／通知偏好和created_at精确回读通过。

实际定时turn01a1335d-e779-73c3-ac2a-b31ce0a87cdb于03:45:54–03:46:09UTC追加到原watcher，完整11,446字符instructions一致，实际model gpt-5.6-luna／effort low由当前turn context核对。只查当前两批squeue，两个末端都有返回：原LR两个fit运行、首个新encode13548914已RUNNING，其余新阶段正常资源等待／有效依赖；静默无回传。主随后核验该首encode的started、STATUS／JOBID和真实GPU分配，确为原E13499817 final9216、CPU资格13548771、task0／NC_003070.9／30427671，与实际命令一致。开始运行不是完整势、后验或排序完成。终态消失修正规则及完整自动回传仍待真实终态，不以一次正常状态检查冒充验证。完整论文目标仍未达成。

04:08UTC另追加零GPU的原TRAIN漏检评分诊断13549479，afterok本评价13548940，使用同批数组、不改变本25个作业、权重、预测、后验或排序定义；见[双模型漏检协议](full_grammar_e_m_missed_reference_plan.md)。同一POSTTRAIN监控事件末端现为13549479，共26roots／78终态记录，13548940为中间阶段；原LR事件不变，两个handoff均NO。完整13,434字符扩展心跳配置精确回读通过。04:15:54–04:16:09UTC实际定时turn01a13379-5f59-7573-aa8a-d85329c1bdc5在原watcher核对同一完整扩展prompt／Luna-low，两个末端有效依赖，部分编码和六个CPU后验已运行、无回传；只证明接线和当次正常状态，不证明新诊断完成、终态消失处置或完整结果回传。

## 整染色体数值失败与定位

05:29UTC主调度核验发现评价13548940已为DependencyNeverSatisfied。原CPU后验13548915／17／19／21／23／25／27／29／32共九项root和batch已FAILED1:0，extern0:0；三个其他M任务仍运行。两份代表性stderr一致停在 `complete_gene_set_posteriors` 的前向／反向总质量或空集守卫，未产生完整后验summary。不是OOM、科学性能负结果或已完成概率读出；GPU角色数组独立完成，仍可用于其他原TRAIN对照。旧失败、尚未运行的评价和漏检诊断记录保留，不取消或用成功口径覆盖。主已消费该actionable失败，watcher不补发、不再把旧POSTTRAIN当正常未完成事件。

原错误消息没有保存logZ／reverseZ差值，因此尚不能把原因定性为舍入误差。本轮另做一次2CPU／32GiB／1h纯CPU只读定位，原E task0／NC_003070.9／plus，source emissions13548914与failed posterior13548915。按相同原DNA、全角色数组、query boundaries、cap128和未修改的前后向例程计算两端总质量，记录有限性、实际差值及原rel1e−12／abs1e−8容差是否通过。该诊断不改算法或容差、不输出gene后验／新GFF、不访问真实参考、DEV或sealed，新增产物小于1MiB。若差值不能由明确数值机制解释，不放行结果；数值修复及实际重跑要依据这份证据，而非直接加大容差。定位成功只意味着失败证据完整，不代表后验成立。输出独立 `outputs/M28-FULL-GRAMMAR-POSTERIOR-MASS-DIAGNOSTIC`，同一watcher中作为独立MASSDIAG事件等待，不新建监控。

实际定位job13551206于05:38UTC接受，命令 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-POSTERIOR-MASS-DIAGNOSTIC.sbatch`；AST／SBATCH语法与两个源文件精确同步回读先通过。原算法代码、原容差和旧失败输出均未改写。旧监控05:15UTC只看末端仍有Dependency而漏掉已失败上游，主已纠正同一心跳的最小检查：明确查询全部尚未消费上游与末端，每个缺行立即批量sacct及markers，不能以末端仍存在或merged exit0掩盖父阶段失败。新心跳三个独立事件LR13545908／POOLTRAIN13551150／MASSDIAG13551206，旧POSTTRAIN失败closed；未发生真实自动回传前不声称回传修复已经验证。

## 实际质量差值与稳定候选

定位13551206于05:38:13–06:29:45UTC运行3152秒，root／batch／extern全部COMPLETED0:0，STATUS／JOBID／started／summary完整。原E plus全30,427,671bp、19,078状态，前向logZ844406.7357152087、反向844406.7357164063，均有限且空集守卫通过。绝对差1.1975644156336784e−6，原允许差8.444067357164063e−7，约1.418232倍；相对差1.418232e−12，仅超过原rel1e−12。因此失败定位为质量一致性守卫，不是非有限值或负logZ；长序列累计舍入是优先验证的解释，不因差值小就放行或增大容差。原算法／容差保持，程序3077.459555秒，MaxRSS10,111,536K。

同一watcher实际turn01a133f4-fac2-7552-bb47-c05cb39ec820于06:30:55–06:32:58UTC核完整root／batch／extern及markers／结构白名单，并一次真实回传原main后暂停同heartbeat。主实际收到并独立消费，handoffYES真实；这验证本次数值定位的完整自动交接，不推广成任意未来失败都已测试。已消费POOL不重监／不补发，原LR继续。

新稳定候选每256bp将全部forward或backward状态log mass统一减去当前background分量，用补偿加法累计常数，查询与最终边界完整加回。它只改变数值表示，既不关闭内部边界、不删状态／路径、不改变DNA／角色势／能量／语法／cap／null，也不改gene事件或set概率公式。默认256是固定数值块长，不由正确率选择。保持原rel1e−12／abs1e−8质量容差和P≤1的原守卫，不拟合温度或阈值。

CPU资格13552225于06:42:34–06:42:56UTC运行22秒，root／batch／extern全部COMPLETED0:0；八测试通过：新四组涵盖独立全部set／gene事件、跨query与缩放块不闭合路径、百万bp分离ORF的解析partition／event概率，以及原API／源不变；另在候选数值kernel上跑原四组概率、真实carry、GT／GC／N、cap128、空集、gauge与log下溢检查。原probability／mass容差不变，程序5.515652秒，MaxRSS989984K，没有研究DNA、参考、checkpoint或native读取。候选在单独测试进程接入，原生产源码未替换。主已消费资格，watcher NOT_SENT，不补发或把人工通过当成整染色体通过。

实际真实资格13552583于06:57:42UTC接受，06:57:43开始，命令 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-SCALED-MASS-NATIVE.sbatch`。仅原E task0／NC_003070.9／plus，用原FP32角色数组及全部原MAP查询链，在该资格进程临时替换两条数值kernel；磁盘上的生产算法及旧诊断仍不变。要求原质量／probability守卫通过，全部查询literal-DNA合法、原GFF八有效位gain按原relative6e−8／absolute1e−5保留，并输出logZ差值和全部查询logP。无真实参考、预测／weights更新、新GFF、nativeforward、训练或DEV／sealed。完整范围是2模型×6染色体×2方向，共24方向；本次只覆盖E0 plus一个，其他23方向尚未资格验证。已运行脚本的布尔字段名 `other_eleven_orientations_not_qualified` 命名不准确，不作为覆盖数量；以实际model／task／strand判定，不为这项记录命名问题取消或重启数值任务。

资源2CPU／32GiB／1h／无GPU／private／generanno／Requeue0，基于原51分32秒加小余量；输出独立 `outputs/M28-FULL-GRAMMAR-SCALED-MASS-NATIVE`，声明小于4MiB、不存新数组。三个候选Python／SBATCH和该真实资格源AST／Bash及六份叶精确同步先通过，live own expanded32加新job低于80、private内存充足。所有原失败记录与数据保留，不取消作业或自动延时。结果前不替换生产、不放行完整十二方向或重新解释成accuracy。

主消费MASS后，复用同heartbeat恢复LR13545908与SCALEDNATIVE13552583两个独立NO事件，原LR合同保留。SCALED只有全三终态、STATUS／JOBID／started／summary及固定结构身份才交接；其科学数值只由main解释。没有新对话／自动化／cron、范围／冻结DEV／goal变化。

## R1读出准备

真实资格13552583仍运行时，已准备独立R1脚本与两个SBATCH，尚未提交正式重跑。输出为 `outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-R1`，保留旧失败目录。运行入口必须先同时核验旧坐标／ranking资格13548771、新八测试13552225，以及真实E0 plus资格13552583的完成身份、原质量守卫、literal-DNA合法性与原GFF gain保留；缺任一项即停止，不据人工通过或正在运行状态放行。

R1复用原十二份completed角色数组、固定E/M weights／原GFF全部唯一查询、24方向的共同语法和原数值容差。仅在读出进程接入已资格验证的两条数值kernel，磁盘生产源码不替换。E0 plus复用真实资格中的全部已完成查询，要求原ID顺序完全相同，避免再计算该方向；同task minus与另外十一双向任务仍实际计算。复用时间单列，不当成全pipeline速度优势。十二CPU后验仍各2CPU／32GiB／4h，完整评价2CPU／8GiB／30min，private／generanno／no-requeue，无新GPU／forward／训练／数组／GFF。

评价入口先按实际十二个新job IDs核完整R1数值身份，再复用未修改的原评价器；只在全部十二双向后验完成后读取原参考。原prediction-pool排序、固定保留比例／十bin、full-reference分母和old exact-chain回放全部不改，不加权重／阈值／校准选择或放行DEV／sealed。已通过新Python AST及两个SBATCH语法；代码准备不是实际接线通过或完整后验结果，真实资格与运行收据仍须核验后才提交。

后接漏检诊断的R1接线亦已准备、未提交。原诊断曾将13548940写死，改为核对传入的实际评价job与completed parent记录；原默认评价路径和输出不变，新增两个明确路径参数供R1使用。原角色数组仍来自旧completed emissions，reference profile／评分／grammar／replacement margin／阈值／完整22,688分母不改。

新 `M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE-R1.sbatch` 接收实际R1评价ID，读取R1 evaluation，写独立 `outputs/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE-R1`，源资格仍13549411。2CPU／16GiB／1h／无GPU与原诊断相同，必须afterok实际R1评价且运行时完整parent守卫通过。该更改只解除旧job绑定，不提前放行参考profile或跳过后验完整性。新源AST与SBATCH语法已通过；原失败输出和所有已有作业保持，不取消或重提。

## 原生资格通过与正式R1启动

真实资格13552583于06:57:43–07:50:13UTC运行3150秒，root／batch0:0、extern07:50:14结束0:0，STATUS／JOBID／started／449,912-byte summary齐全。原E0 plus全30,427,671bp、19,078状态，稳定前向和反向logZ在当前FP64表示下均为844406.735716842，差0；原未缩放差1.1975644156336784e−6。原rel1e−12／abs1e−8和gene／set概率守卫均通过，没有放宽容差。该数值行为加上独立八组人工验证，支持累计舍入解释及同数学模型的稳定表示修复，不将一个方向一致当成24方向或生物学正确性。

3467个原MAP查询均literal-DNA合法、原八有效位GFF gain保留；主独立复核ID唯一、全部gain／logP有限、logP≤0、gain总和843912.1072089081与exact-set logP−494.62850793392863及empty-set logP−844406.735716842一致。源权重、DNA、原GFF与磁盘生产代码未改，未读参考或DEV／sealed。程序3124.913857秒，MaxRSS10,184,944K。原watcher实际07:50:25–07:51:42UTC turn01a1343d-c555-7d42-9c3b-dfd7bfe13158核完整结构后一次真实回传并暂停；主收到、取证消费，handoffYES真实。本次完整单方向资格交接通过，不表示未来整R1终态已验证。

R1接线CPU13554449于07:58:07–07:59:36UTC运行89秒，root／batch／extern均COMPLETED0:0，STATUS／JOBID／summary完整。十二份completed source身份与实际E0plus查询ID顺序／有限logP／原GFF gain复用验证通过；没有研究DNA、参考、角色数组、checkpoint或native读取。程序7.657586秒，MaxRSS893988K。主已消费，watcher NOT_SENT，不补发短资格。绑定检查检测的是新R1接线与真实复用，而非重复数值数学测试。

正式十四项于08:01UTC（北京时间16:01）接受，均纯CPU。各posterior命令 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-R1-READOUT.sbatch <model> <task>`；eval使用下面十二个实际ID的afterok，并按E0–5／M0–5顺序传入全部ID；diag afterok实际eval13554510，并以该ID作为唯一参数。

| model | task | TRAIN染色体 | R1 posterior job |
| --- | ---: | --- | --- |
| E_RANK4 | 0 | NC_003070.9 | 13554497 |
| E_RANK4 | 1 | NC_003071.7 | 13554498 |
| E_RANK4 | 2 | NC_003076.8 | 13554499 |
| E_RANK4 | 3 | NC_089043.1 | 13554500 |
| E_RANK4 | 4 | NC_089044.1 | 13554501 |
| E_RANK4 | 5 | NC_089046.1 | 13554502 |
| M_DENSE | 0 | NC_003070.9 | 13554503 |
| M_DENSE | 1 | NC_003071.7 | 13554504 |
| M_DENSE | 2 | NC_003076.8 | 13554506 |
| M_DENSE | 3 | NC_089043.1 | 13554507 |
| M_DENSE | 4 | NC_089044.1 | 13554508 |
| M_DENSE | 5 | NC_089046.1 | 13554509 |

eval13554510 afterok全部十二posteriors；漏检13554511 afterok13554510，共14roots／42root-batch-extern。所有十四controller资源、Requeue0和依赖已核验：十二posterior均RUNNING（十个gpu034、两个gpu035），eval／diag有效Dependency。E0plus已经实际打印native_qualification_reused=true，不再重算已合格方向；E0minus与另外十一双向任务共23方向继续。复用约0.002秒仅是读取／组装成本，不包含原3150秒资格，不能称整体计算提速。

提交前own expanded33加十四仍低于80，private CPU／内存与161TiB剩余storage足够，维护仍仅11月CPU预留。R1含曲线与漏检声明不超过1GiB，旧35GiB数组只读复用，不保存新数组／DNA／weights／GFF。新接线Python AST／相关SBATCH语法及正式五份source精确叶同步回读先通过，没有依赖controller已移除的old source IDs。

同一heartbeat从资格完成暂停中复用恢复，当前两独立事件为LR13545908／POSTTRAINR1末端13554511，两handoff均NO；原LR完整合同保留。新prompt12702字符保存完整14stage／42record合同、严格结构白名单及closed资格事实。正常运行／有效依赖静默，任一完整事件或需处理失败才一次交接并暂停同心跳；没有新对话／重复自动化／cron。完整R1概率排序与漏检结果尚未完成，不把submit／部分复用或单方向资格当成竞争性模型、独立泛化或论文成功。

R1新范围配置精确回读通过：同ID／heartbeat／ACTIVE／原target／15min／failed_runs_only／created_at及完整12702字符prompt匹配。08:20:56–08:21:19UTC实际定时turn01a13459-b264-7d40-971c-f2be17915d27追加到原watcher，完整instructions与实际gpt-5.6-luna／low一致；其最小调度查询包括原LR全部上游／arrays／末端和十四R1 actual IDs，十二posteriorRUNNING、eval／diag有效Dependency，无科学读取或main回传。主08:20快照独立核同十四句柄。这验证实际定时绑定、当前范围和当次正常检查；全十四阶段终态及完整自动回传仍待实际发生。


## R1 后验排序评价完成

2026-10-13，评价13554510已完成。后验比内部gain更能排序参考匹配，两个物种和两个模型方向一致；这是原TRAIN的selected-prediction-pool效用，不是annotation detection AP、独立泛化或生物学校准。漏检诊断13554511仍运行，不能将这个已完成评价阶段等同14-stage整批完成。

十二个后验及评价共13roots／39root-batch-extern均COMPLETED0:0。十二后验08:01:17UTC开始，末项M0于09:47:34UTC完成；评价09:47:36–09:54:24UTC运行408秒，MaxRSS1690992K，程序319.467538秒，实际summary198329291bytes。24个source方向的整染色体messages均通过原mass守卫，其中E0plus来自已核验数值资格，另外23方向实际重新计算；原概率／mass容差、权重、GFF、reference规则和冻结DEV均不变。评价保留全部22,688主参考，原E22324预测／16420主TP、M21927／16106精确重放。

主通过Slurm srun任务13559688读取198MB结果，独立重算所有六chr、两物种、pooled两种标签／两种score的完整同分曲线AP／AUROC、precision／完整主参考recall，重建chr→species→pooled数量与positive分区，全部一致；ten-bin数量和matches也与原预测pool对应。只读源结果，无新forward、训练或预测。紧凑结果副本保存在`outputs/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-R1-RESULT-READOUT/compact.json`，完整曲线仍在原remote evaluation/summary.json；不是新计算的独立实验。

### 主参考链匹配排序

| 模型 | 范围 | 主TP／预测 | Gain AP | Posterior AP | Gain AUROC | Posterior AUROC | Reference-relative Brier |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| E_RANK4 | AT | 13111／17077 | 0.772325 | 0.913445 | 0.558478 | 0.795183 | 0.154309 |
| E_RANK4 | Rice | 3309／5247 | 0.681592 | 0.883152 | 0.610011 | 0.844101 | 0.193207 |
| E_RANK4 | Pooled | 16420／22324 | 0.751671 | 0.907840 | 0.575052 | 0.813012 | 0.163452 |
| M_DENSE | AT | 12873／16859 | 0.754228 | 0.914697 | 0.533935 | 0.800102 | 0.164082 |
| M_DENSE | Rice | 3233／5068 | 0.674407 | 0.890623 | 0.595922 | 0.846322 | 0.197589 |
| M_DENSE | Pooled | 16106／21927 | 0.734968 | 0.910230 | 0.552137 | 0.816459 | 0.171826 |

Pooled主链AP：E0.907840对gain0.751671，M0.910230对0.734968；AUROC：E0.813012对0.575052，M0.816459对0.552137。此结果支持保留完整gene-set后验读出作为后续研究信号，不将AP或AUROC当作整个注释任务的召回／F1，也不由此宣布超过ANNEVO。

### 任意可评价 isoform 匹配

| 模型 | 范围 | 匹配／预测 | Gain AP | Posterior AP | Gain AUROC | Posterior AUROC | Reference-relative Brier |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| E_RANK4 | AT | 14432／17077 | 0.885479 | 0.969857 | 0.659605 | 0.869541 | 0.091831 |
| E_RANK4 | Rice | 3437／5247 | 0.713210 | 0.903840 | 0.627148 | 0.857266 | 0.173439 |
| E_RANK4 | Pooled | 17869／22324 | 0.847113 | 0.958910 | 0.647295 | 0.868730 | 0.111012 |
| M_DENSE | AT | 14200／16859 | 0.864526 | 0.970653 | 0.618046 | 0.874581 | 0.097754 |
| M_DENSE | Rice | 3367／5068 | 0.708887 | 0.912732 | 0.613957 | 0.861765 | 0.175583 |
| M_DENSE | Pooled | 17567／21927 | 0.826724 | 0.960754 | 0.614535 | 0.873659 | 0.115743 |

任意isoform事件不借用主参考总数计算all-reference recall；非主链匹配不等于假基因。这一标签下后验排序同样优于gain，不能仅用primary-negative数量解释生物学错误。

### 全部固定保留比例

以下为pooled主链标签，未搜索最佳阈值、未筛写GFF。两score在这些固定比例恰好保留相同预测数量；同分组保留与ceiling会令M实际比例略大于请求比例。R_all分母始终为全部22,688主参考，而不是预测内positive数量。

| 模型 | 请求保留比例 | 预测数 | Posterior TP | Posterior precision % | Posterior R_all % | Gain TP | Gain precision % | Gain R_all % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| E_RANK4 | 25% | 5581 | 5215 | 93.442035 | 22.985719 | 4187 | 75.022397 | 18.454690 |
| E_RANK4 | 50% | 11162 | 10098 | 90.467658 | 44.508110 | 8593 | 76.984411 | 37.874647 |
| E_RANK4 | 75% | 16743 | 14337 | 85.629815 | 63.191996 | 13019 | 77.757869 | 57.382757 |
| E_RANK4 | 90% | 20092 | 16046 | 79.862632 | 70.724612 | 15357 | 76.433406 | 67.687764 |
| E_RANK4 | 100% | 22324 | 16420 | 73.553127 | 72.373061 | 16420 | 73.553127 | 72.373061 |
| M_DENSE | 25% | 5482 | 5142 | 93.797884 | 22.663963 | 3936 | 71.798614 | 17.348378 |
| M_DENSE | 50% | 10964 | 9965 | 90.888362 | 43.921897 | 8279 | 75.510762 | 36.490656 |
| M_DENSE | 75% | 16446 | 14043 | 85.388544 | 61.896157 | 12699 | 77.216344 | 55.972320 |
| M_DENSE | 90% | 19735 | 15705 | 79.579427 | 69.221615 | 15071 | 76.366861 | 66.427186 |
| M_DENSE | 100% | 21927 | 16106 | 73.452821 | 70.989069 | 16106 | 73.452821 | 70.989069 |

预定义90%预算下，后验比gain多保留E689、M634个主TP，但相对原未过滤输出仍丢E374、M401主TP；分别只保留原输出positive的97.722290%和97.510245%，全参考recall为70.724612%／69.221615%，低于原72.373061%／70.989069%。这说明排序有用且有召回代价，不是新增基因或全面性能胜出。没有选择90%作为部署／冻结DEV operating point，也未重算过滤后的base／背景。

### 高概率区间的参考相对差异

| 模型 | 范围 | P结构 ≥0.9 的预测数 | Mean structural P % | Exact-primary fraction % | Any-assessable-isoform fraction % |
| --- | --- | ---: | ---: | ---: | ---: |
| E_RANK4 | AT | 13016 | 99.002137 | 86.800860 | 94.184081 |
| E_RANK4 | Rice | 3364 | 98.672730 | 82.312723 | 85.077289 |
| E_RANK4 | Pooled | 16380 | 98.934486 | 85.879121 | 92.313797 |
| M_DENSE | AT | 13297 | 99.197304 | 85.756186 | 93.494773 |
| M_DENSE | Rice | 3381 | 98.955475 | 82.194617 | 85.122745 |
| M_DENSE | Pooled | 16678 | 99.148280 | 85.034177 | 91.797578 |

Pooled高概率区间mean structural P为E98.934486%、M99.148280%，但主参考匹配仅85.879121%／85.034177%，任意可评价isoform为92.313797%／91.797578%。结构包含与参考正确事件不同，差异可能包括替代isoform、参考不完整及模型错误；不把这里的P当作基因正确率，不称已完成生物学校准。未拟合温度／校准器。

### 对下一项研究的影响

数据验证结论：可用于带上述局限的原TRAIN机制决策，不能作为独立论文性能／泛化／校准结论。后验排序值得继续，等权组合的同口径后验已并行运行；其完成后检验能否维持新增TP并更好区分非主链匹配输出。先结合漏检13554511的score／grammar证据和LR四格结果，决定下一项原TRAIN objective／表示改进，不把排序改善替代模型本身的竞争性。冻结DEV family、sealed／held-out、新panel和公开claim范围不变，完整论文目标未完成。

