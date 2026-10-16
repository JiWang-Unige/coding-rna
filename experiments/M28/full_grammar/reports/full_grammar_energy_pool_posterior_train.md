# 原 TRAIN 固定组合的完整基因后验诊断

2026-10-13。固定0.5／0.5组合的完整链F1为73.857085%，仍低于ANNEVO的75.541796%。同一22,688主参考分母下，组合16527／22066、ANNEVO16470／20917：组合净多57个主TP，却多1092个非主链匹配输出。因此读出组合模型完整gene-set分布中的gene inclusion概率，检验其排序是否比原gain更能区分参考匹配；不是寻找更好看的阈值、拟合校准或修改原预测。

完成结果：全部十二方向／七正式root／21条Slurm记录成功；组合内posterior AP0.910167高于gain0.761573，但posterior AUROC0.805353低于E／M各自后验，不能称为全面概率提升。固定90%保留预算较gain多603个主TP，却较未过滤输出少471个TP；没有选择该预算为operating point或生成新GFF。完整结果如下。

## 比较者差距与研究问题

AT组合TP13080／预测16938，ANNEVO12931／15922，组合多149TP、867个非主链匹配输出；rice组合3447／5128、ANNEVO3539／4995，组合少92TP、多225个非主链匹配输出。总体precision负担是主要差距，但rice还存在召回短板，不能以pooled数值代替每物种判断。非主链匹配含其他可评价isoform，不能统一称为生物学假基因。

当前组合权重、E fit13499817／M fit13510324末head9216／native4608、原六TRAIN输入和已完成GFF全部固定。E/M势源为原十二个completed双方向emission任务，组合raw13551143／45／46／47／48／49、评价13551150，实际raw入口原生匹配13556298。只重新读取这些势和既有GFF，不训练、不做native forward、不择checkpoint、不生成新预测或全长势文件。

## 概率与数值资格

组合逐角色先promote至FP64再平均，背景差与完整set能量使用同一mean模型。对既有合法非重叠链查询`log P(g included)=alpha_B(start)+E(g)+beta_B(end)−logZ`，分母仍是整条oriented染色体上的全部合法gene sets，含能量0的空集，而不是已选链池。原gene posterior输入API已支持FP32／FP64；只在隔离读出进程使用已核验的分块中心化＋补偿求和messages，原源文件、语法、概率和前后向mass容差不改。

既有数学资格13552225含FP64独立oracle和原API，whole E0正向数值资格13552583保留。新增两个合成资格组专门验证组合后FP64全部set／gene事件与独立oracle一致，以及promote-before-mean保留半单位能量；两支源数组须保持不变。十二个source与六个pool完成产物先核验接线。资格不使用研究DNA／参考／数组／checkpoint；完整染色体适用性仍由正式十二个方向的原mass guard、合法链、原GFF gain逐一验证，短资格不能替代它。

## 固定诊断与评价边界

六个posterior全部完成并核验后才读取原TRAIN参考。复用原M26 CDS-assessable主链及所有assessable isoform规则，保留22,688主参考和unsupported strata；逐染色体重放完整预测数／主TP／参考数，并与E/M R1评价13554510的分母一致。

沿用同一`gene_posterior_ranking`和`P.diagnostics`：posterior versus内部gain的selected-prediction-pool AP／AUROC、固定0.25／0.5／0.75／0.9／1 retained fractions、tie整组保留，十个固定等宽probability bins和reference-relative Brier。AP／AUROC不含未预测参考，不是annotation detection AP。结构inclusion不是参考正确率，概率bin偏差不能称为生物学校准。完整tie曲线仍用于计算，但不在三种聚合层反复保存；全部原query/log保留在六份posterior中。报告六chr、两物种、pooled，以及相同E/M诊断的AP／AUROC／Brier差；不挑最佳cutoff、不拟合温度、也不重新评测过滤后base／背景或声称全面改善。

## 作业范围与完成

一个2CPU／4GiB／5min新资格；通过后六个2CPU／32GiB／3h CPU后验任务并行，每项覆盖一条完整原TRAIN contig的正负两个方向。已测原E0单方向约52.5min，双方向约1.75h；3h包含组合数组读入、不同chain数和调度负载的余量，不自动延长。最后2CPU／8GiB／30min评价，afterok六个实际新后验ID及当前R1评价13554510。源complete旧IDs只作产物身份，不写离开controller的源依赖。全部generanno／private-teodoro-gpu／no-requeue／无GPU；固定新输出`outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-TRAIN`，新增日志／gene metadata／summary声明64MiB，不复制35GiB源数组。

完整七项正式作业及全部root／batch／extern终态、各STATUS／JOBID／started／非空summary、每posterior双方向全contig和完整query列表、评价全原分母重放才构成完成。确定性接线失败最多同范围一次修复，新独立输出保留首次失败；数值或科学负结果不隐藏、不抬容差或挑阈值，不取消当前作业。等待只复用原watcher和同一heartbeat追加POOLPOST独立事件，LR与POSTTRAINR1未交接合同保持不变；无需新chat、cron或重复automation。当前TRAIN探索不能扩展冻结DEV family、封存数据、新物种或论文claim，Nature Communications及以上完整目标仍未达到。

## 资格与实际运行

两份Python AST、三个SBATCH逐项语法和六份新叶文件精确同步回读通过。新资格13556704于09:30:24–09:31:43UTC运行79秒，root／batch／extern全部COMPLETED0:0，MaxRSS895796K。两个新独立数学测试通过，tests runtime2.503秒、程序14.694759秒；十二个emission源和六个组合输出全部接线核验通过，原数值容差未改，未读取研究DNA／参考／数组／checkpoint或执行native。资格由主直接消费，CLOSED／NOT_SENT_BY_WATCHER，不让watcher补发资格交接。

六条正式后验于09:33:36UTC接受，均afterok实际新资格13556704；评价13556802于09:34:31UTC接受，afterok六个实际新后验ID及原R1评价13554510。实际ID按task传参，不按数字相邻推断。

| task | TRAIN染色体 | 后验job |
| --- | --- | --- |
| 0 | NC_003070.9 | 13556785 |
| 1 | NC_003071.7 | 13556786 |
| 2 | NC_003076.8 | 13556787 |
| 3 | NC_089043.1 | 13556788 |
| 4 | NC_089044.1 | 13556789 |
| 5 | NC_089046.1 | 13556790 |

提交前own expanded42，新增七项仍小于80；CPU-only不消耗GPU，目标新输出不存在，3h维护窗口允许。首项13556785于09:34:42UTC在gpu034运行，controller核对实际2CPU／32GiB／3h／noGPU／Requeue0、命令和source参数；其余五项正常Priority／Resources等待，评价七个afterok依赖正确。以上只是实际执行接线，不是后验完成或排序结论。

同一heartbeat `coding-rna-r5-light-monitor` 追加POOLPOST末端13556802与六upstream、7roots／21records产物／结构白名单，原LR和R1正文保留。16,673字符完整prompt精确回读一致，原ID、heartbeat类型、ACTIVE、watcher01a0f71f-8845-7432-bd13-bf54c1808fb3、15min、failed_runs_only和created_at不变，没有新chat或automation。09:50:56–09:53:29UTC实际定时turn01a134ac-19ec-76f0-a973-b8795909734d追加到同一watcher，完整16,673字符instructions精确匹配，实际gpt-5.6-luna／low由该turn_context核对；其已完成exit0的调度命令明确包括六个新posterior／评价13556802以及原LR／R1上游和末端。新六项均正常运行、评价有效Dependency，R1十二posteriors全部结构终态通过而整链未结束，无科学字段读取或主回传。首次新范围定时检查已验证，但七项完整终态及实际回传仍未发生。记录按现有格式区分新资格、实际提交、完整概率数值有效性与科学解释；没有公开发布。


## 完整读出结果与独立主参考核验

六项posterior及评价共7roots／21条root-batch-extern全部COMPLETED0:0。后验最长6299s，最大MaxRSS15129708K；全部十二方向的原mass／event guard和原GFF gain重放通过。评价13556802于11:23:19–11:26:14UTC／175s，MaxRSS1751228K，程序105.371238s，summary166042bytes。STATUS／JOBID／started／双方向查询与完整原分母齐全，原权重／预测／数值容差未改。

独立读出核验13566913（srun，2CPU2GiB5min/noGPU/generanno）实际11:41:36–11:41:42UTC／6s，root、extern及step0全部COMPLETED0:0、MaxRSS98440K。以已保存13449803主参考全几何key重新标注全部22,066查询，独立grouped-score算法复算六chr／两species／pooled主链AP／AUROC、五固定比例／整组tie、Brier和十概率bins；完整22,688主参考与16,527TP重放，全部与正式评价一致。不读DNA／参考GFF／势／权重或做forward。任意isoform数值来自已完成评价及既有同算法数学资格，本额外核验未重新读取isoform参考或独立重算该标签；不要把核验范围扩大成全部指标的额外独立复现。可复现脚本`verify_full_grammar_pool_posterior_readout.py`，小输出`outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-READOUT-VERIFIED/13566913/compact.json`。

| 标签／范围 | 匹配／预测 | Gain AP | Posterior AP | Gain AUROC | Posterior AUROC | Reference-relative Brier |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 主链／AT | 13080／16938 | 0.769036 | 0.914574 | 0.545438 | 0.792764 | 0.158146 |
| 任意assessable isoform／AT | 14417／16938 | 0.882222 | 0.971646 | 0.642549 | 0.871595 | 0.091891 |
| 主链／Rice | 3447／5128 | 0.738809 | 0.892689 | 0.637693 | 0.830905 | 0.194129 |
| 任意assessable isoform／Rice | 3584／5128 | 0.774480 | 0.914621 | 0.661412 | 0.847284 | 0.171132 |
| 主链／Pooled | 16527／22066 | 0.761573 | 0.910167 | 0.572642 | 0.805353 | 0.166508 |
| 任意assessable isoform／Pooled | 18001／22066 | 0.857335 | 0.961563 | 0.648571 | 0.865427 | 0.110306 |

| 固定请求比例 | 预测数 | Posterior TP | Posterior precision % | Posterior 全主参考recall % | Gain TP |
| --- | ---: | ---: | ---: | ---: | ---: |
| 25% | 5517 | 5170 | 93.710350 | 22.787377 | 4177 |
| 50% | 11033 | 10042 | 91.017856 | 44.261283 | 8624 |
| 75% | 16550 | 14251 | 86.108761 | 62.812941 | 13078 |
| 90% | 19860 | 16056 | 80.845921 | 70.768688 | 15453 |
| 100% | 22066 | 16527 | 74.898033 | 72.844676 | 16527 |

固定90%下posterior16056／19860，比gain15453多603TP；相对未过滤16527仍少471TP，全部主参考recall从72.844676%降至70.768688%。没有搜索最优比例或生成新GFF，也没有新base／背景评价。任意isoform匹配的分母只在预测池内使用，不借22,688当作isoform总数。

Pooled primary posterior AP0.910166621，比E0.907840293高0.002326328、比M0.910230273低0.000063652；AUROC0.805353324比E低0.007658778、比M低0.011105912。Brier0.166507918高于E0.163451632、低于M0.171826486。各模型预测集合、positive prevalence及结构分布不同，不能把这些差值当作同一population上的配对性能／校准胜出。对组合自身，两物种的posterior AP／AUROC均高于gain；未进行统计显著性检验，亦无整体概率优势claim。

P≥0.9的16988查询平均结构概率99.095888%，主参考匹配85.501530%，任意isoform92.182717%；结构事件与参考正确事件仍不等价，不称生物学校准。评价问的是预测池排序，不是全参考annotation detection AP；其不能恢复未预测参考。数据验证结论为附这些限制可用于原TRAIN下一项机制选择，不可用于独立泛化或论文优越性主张。

同一原watcher实际turn01a13514-05bd-7e63-9cd7-f9fcc3c9baa1于11:44:27–11:46:20UTC，gpt-5.6-luna/low，已按新scope核验七root／21记录与结构产物，并真实一次发送原main后暂停；main已收到。R2同轮结构核验完成但未第二次交接，主已独立消费R2，不补发假YES。下一项固定局部能量分歧惩罚计划只在原TRAIN，LR对照继续，冻结DEV/封存/原论文goal不变。

