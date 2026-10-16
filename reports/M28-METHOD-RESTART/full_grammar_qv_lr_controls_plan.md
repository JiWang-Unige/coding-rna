# 原 TRAIN 的 q v 参数化与学习率固定四格对照

2026-10-14完成更新：两新全链及四格全部21root／63终态记录0:0，数组numeric身份疑点已消解并独立复算。新增rank4低lr与全秩高lr的pooled完整链F1分别65.928535%／58.982365%，低于各自原设置；native lr效应方向相反，不是纯rank因果或论文目标达成。见 [完整四格结果和身份映射](full_grammar_qv_lr_controls_result.md)；下方原方案／提交／真实失败记录保留。没有重提、补训、改参数、扩充DEV或放行封存数据。

2026-10-13（北京时间）。当前E低秩q/v在原TRAIN pooled完整链F1为72.958322%，M全秩为72.199933%；二者同时改变native学习率（3e−4与1e−5）、参数化和forward舍入。为检验学习率混杂并继续改善完整链性能，只补齐另外两格，不改原E/M、不延长其训练或添加DEV候选。不是从DEV选择学习率，也不是纯rank因果消融或论文已达标。

## 四格与固定读出

| q/v参数化 | Native lr1e−5 | Native lr3e−4 |
| --- | --- | --- |
| rank4／alpha4／16 FP32残差张量／102400参数 | 新独立JOINT-LORA-LR1E5 | 已完成E fit13499817／eval13499823 |
| 八个全秩FP32 master／18874368参数 | 已完成M fit13510324／eval13510351 | 新独立DENSE-QV-LR3E4 |

两个新fit均从原13434603末4608 C0头及保存AdamW和原b048发布encoder出发，不接续E／M最终9216，也不继承可丢弃验证更新。头lr始终3e−4／age4608→9216，native组age0→4608；只改各参数化内部的native lr，追加仍4608。6144／7680只为保存完整性，唯一读出最终累计9216。每格一固定fit／seed，不换seed或重复直到指标变好。

仅原两物种六TRAIN、同seed1 fresh4608有序draw／4219 unique、原fine／coarse ε0.1／unknown、已知行分母／importance1、joint24576／nearest-center ownership、十五角色／FP64分区／空集／cap128及完整gene-set合法性不变。原头366191参数，低秩总468591、全秩总19240559；两者其余encoder权重冻结，整体encoder均因q/v学习而适配。AdamW weight decay0.01／betas0.9,0.999／eps1e−8及共同clip1不变，不加loss／codon／native-CDS先验或缓存特征。

低秩A仍CPU seed0初始化、B为零、scale1，残差FP32再舍入BF16并相加；全秩仍初始化为原q/v的FP32副本、每次前向cast至原BF16。即使学习率匹配，这些数值／参数化／联合clip轨迹仍不同，不能把相应差异孤立归因给rank。发现阶段的四格方案在历史DEV已知后设计，不能称DEV前预注册；此设置依据TRAIN和实际代码的LR混杂，而非DEV成绩。

## 独立实现与必要验证

保留原全部E/M文件、输出、checkpoint和已完成结果。新入口直接复用原编码／语法／目标／raw与共同评价函数，在独立脚本和输出前缀中固定各格的学习率、模型身份与protocol；源发布资产和原4608权重／moments只读。新头／native组实际lr逐步写入日志，CPU核验各step和三checkpoint的真实groups／名称／moments／age。旧E/M checkpoint、缺native_lr、错head_lr和错age不得被新loader当作当前格使用。

每支路两组新增CPU人工检查：真实fit所用native-group构造／AdamW登记及错误lr／age拒绝；新final checkpoint身份及旧权重／缺lr拒绝。使用人工张量与C0随机fixture，不读研究DNA／标签／reference／checkpoint／native模型，不重复已不变的底层数学或原6/11组全套检查。AST及每个SBATCH语法另外核对，检测新接线／路径错误而非科学效果。

CPU通过后，各分支用B/E/M相同的每物种首个fixed fresh TRAIN draw做两次独立可丢弃更新。从头4608/native0到4609/1，检查初始实际势／NLL／角色VJP与冻结parent逐张量相同、选中native参数有限非零梯度和更新、其他encoder冻结、实际head/native组lr分别正确、训练／raw势与完整context MAP一致，以及14bp人工raw两方向。无checkpoint写出、无全染色体准确率或DEV／sealed读取。主核验成功后才提交正式fit，资格验证不等于科学成功。

## 完整原 TRAIN 端点

六raw只读原FASTA、全contig双链、固定末权重与正确16残差或8master恢复；不输入参考／标签／观察模型／feature cache，不作更新。所有raw完成后，CPU先逐碱基比较新格、B和规范化原TRAIN DNA，再访问参考；共同M26 CDS-assessable最长eligible CDS／transcript ID tie、22688 primary（AT17739／rice4949）、全部unsupported／长span／重叠层保留，不按训练资格筛分母。原生ANNEVO／Tiberius／Helixer仍复用已完成13461307固定版本和原输出。

每新格保留对B的预定完整配对：六chr／两物种／pooled／macro／worst的TP／预测／参考与PRF、CDS-base union、背景、phase／其他isoform／duplicates。四格主比较是各参数化内native高lr减低lr及差中差：(dense_high−dense_low)−(lora_high−lora_low)，分别对完整链PRF百分点、base和背景列出；F1非线性、单fit和优化轨迹限制因果解释。Worst分别在各格内取最差物种，再作差，不以最小物种差代替；negative和物种权衡完整保留。

两新格均须完整报告，不能凭TRAIN高低删格；不根据任何DEV读数调lr、rank、层数、head、epsilon、更新数、checkpoint或下一训练。工程失败／非有限更新／OOM／TIMEOUT不记零准确率，保留FAILED，不自动改科学设置重跑。明确接线／调度错误至多一次同科学命令修复；科学负结果不重训。新四格汇总独立写出，不覆写任何旧TRAIN或8/10/6方法DEV表。

## 资源与停止边界

| 阶段每分支 | 固定资源 | 时限 |
| --- | --- | --- |
| 新CPU检查 | 2CPU／4GiB／无GPU | 5分钟 |
| 原生可丢弃预检 | 单24GiB RTX3090／2CPU／16GiB | 20分钟 |
| 正式fit | 同单3090／2CPU／16GiB | 20小时 |
| 最终权重／draw核验 | 2CPU／4GiB／无GPU | 5分钟 |
| 六TRAIN raw数组 | 每项单3090／2CPU／16GiB，max6 | 每项2小时 |
| 共同评价 | 2CPU／4GiB／无GPU | 20分钟 |
| 对B配对 | 2CPU／4GiB／无GPU | 10分钟 |
| 四格完整汇总 | 2CPU／4GiB／无GPU | 10分钟 |

均private-teodoro-gpu、generanno和no-requeue，不安装依赖。E/M实测fit约13.3小时，本次20小时为同类反传留余量，不是完成保证；不得自动延长或从中间snapshot补训。低秩fit96MiB、全秩fit768MiB上限，raw各128MiB，全部预计新增小于4GiB。提交前核对当前配置／Slurm并发／VRAM／maintenance／storage和输出碰撞，不取消其他项目或旧作业。两条fit可并行；后续两组六raw最多12项，仍受当前80 expanded任务与实际GPU容量限制。

独立输出分别outputs/M28-FULL-GRAMMAR-JOINT-LORA-LR1E5-{TRAIN-CHECKS,PREFLIGHT,FIT,VERIFICATION,TRAIN-RAW,TRAIN-DIAGNOSTIC,PAIRED}与DENSE-QV-LR3E4同后缀，四格汇总为M28-FULL-GRAMMAR-QV-LR-FOUR-CELL。实取得jobIDs后记录，不猜数组numeric roots。正常长等待只更新原固定watcher01a0f71f-8845-7432-bd13-bf54c1808fb3及原coding-rna-r5-light-monitor心跳，原main01a0332e-e0a9-7aa2-90f9-159405701833、Luna/low、15分钟不变；不新建对话或重复自动化。已消费M/E及旧批关闭、不再交接；短检查由主核验。

范围仍仅原TRAIN；E/M及两个新格均不进入冻结B/C/D DEV，不访问sealed／held-out／Setaria／新物种／新面板，不对外发布。研究目标保留Nature Communications及以上的竞争性性能或独立更广适用性、机制与论文；本轮提交／验证／负结果都不把它缩减为工程完成或宣布达标。

## 实际接入检查与原生等待

2026-10-13 00:30:53UTC（北京时间08:30:53）提交前own expanded=9，新增四个资格jobs后13<80；private gpu034/035的16张3090当时均已有GRES分配，shared同类也无更早空闲可用，CPU可运行。scratch161TiB可用，只有11月6日CPU维护预留；输出根均不存在。11份Python AST与14份SBATCH语法已通过，精确新脚本传输完成并回读关键七份源／接口后才提交，旧E/M文件未改。

| 分支 | CPU checks实际job | afterok原生preflight实际job |
| --- | --- | --- |
| rank4 native lr1e−5 | 13544938 | 13544940 |
| dense native lr3e−4 | 13544939 | 13544941 |

两CPU均00:30:55–00:31:37UTC、42秒、root/batch/externCOMPLETED0:0；各两组新检查通过、0error/failure，脚本实耗0.686055/0.726476秒，不是native模型速度。配置核对head lr3e−4、native1e−5/3e−4、seed1/joint24576及真实new identities；不读研究DNA/reference/checkpoint/native模型/DEV/sealed。人工toy AdamW step用于组登记验证，不是研究模型更新，研究optimizer steps=0。两CPU由主消费，不加入等待对象。

00:33:57UTC实际13544940 PENDING/Resources、13544941 PENDING/Priority，两CPU依赖已满足；这是真实正常GPU容量等待，不是失败。Controller给出的预测启动时间只为排队估计，不作保证。原生模型未开始加载／更新，没有最终研究权重或准确率；正式fit与下游仅已准备，尚未提交。主先完成原生结果核验，再实施同已固定两条完整TRAIN流程及四格结果处理，常规步骤不重复请求批准。

同一原heartbeat更新为只监控上述两原生资格job，绑定原watcher/main、Luna/low、15分钟和通知偏好不变；完成两支或需行动失败一次回传，主处理后再沿用同一心跳注册实际新fit链。不新建对话或自动化，不检查旧M/E/DEV批或已消费CPU。正式研究结果未产生，论文goal继续active且未完成。

## 原生资格完成与正式两链提交

主对话00:50UTC核验两项真实原生结果：LOW13544940于00:38:49–00:41:11UTC完成（extern至00:41:12），HIGH13544941于00:41:13–00:43:14完成；六条root/batch/extern均COMPLETED0:0。脚本实耗124.925199／105.346944秒。每支两项独立固定首draw实际head/native组lr分别为3e−4/1e−5及3e−4/3e−4，head4608→4609、native0→1；不是继承到正式fit的更新。

两支初始实际native十五角色势、完整gene-set NLL及角色VJP均与原parent逐张量相同；两物种raw NLL分别0.0014955102768796722和30.703542945783738，除原已知行分母24576后的weighted total分别6.085246894855437e−8和0.0012493303607496638，主独立复算一致。LOW zero-B时A梯度0、八个B梯度有限非零且更新非零；HIGH八个q/v矩阵梯度有限非零且更新非零。更新后的训练／raw实际势和完整context MAP均相同。14bp人工raw双链接入完成但没有CDS预测，只证明短输入接线，不能证明非空基因预测或准确率。其他encoder冻结、原weights／moments未改写、未保存checkpoint，无DEV／sealed；最大allocated约6.36／5.91GiB，小于24GiB。主直接消费，监控handoff为NOT_SENT_BY_WATCHER，不伪记自动回传。

00:58:45–46UTC（北京时间08:58）在own expanded9、scratch可用161TiB、无相关maintenance且十个正式输出根不存在时，实际提交以下两条完整链；总新增20 expanded任务，低于80上限。

| 分支 | fit | verify afterok fit | 六TRAIN raw afterok verify | eval afterok整raw | paired afterok eval |
| --- | --- | --- | --- | --- | --- |
| rank4 native lr1e−5 | 13545578 | 13545579 | 13545580_0–5%6 | 13545581 | 13545582 |
| dense native lr3e−4 | 13545584 | 13545585 | 13545586_0–5%6 | 13545587 | 13545588 |

提交命令为`sbatch --parsable sbatch/M28-FULL-GRAMMAR-{JOINT-LORA-LR1E5,DENSE-QV-LR3E4}-FIT.sbatch <本支native13544940或13544941> <CPU13544938或13544939>`；verify/raw/eval/paired依次分别调用对应`-VERIFY.sbatch <fit>`、`-TRAIN-RAW.sbatch <fit>`、`-TRAIN-EVAL.sbatch <fit> 13494895`和`-PAIRED.sbatch <eval> 13494898`，只依赖上述实际新job。旧成功资格已在fit程序核验，不将可能已被controller回收的旧IDs挂作新afterok依赖；B作为只读来源而非旧Slurm依赖。十条scontrol回执逐一核对实际命令／参数、private资源、Requeue0、独立logs及全部新afterok关系，均一致；首次队列短暂raw Reason=None已在正式controller回执中变为正确Dependency。

00:59:05UTC LOW13545578 PENDING/Resources、HIGH13545584 PENDING/Priority，下游均有效Dependency。这是正常GPU排队，不是失败；尚无正式训练开始或新完整链性能结果。四格汇总入口／作业尚未实施或提交，不宣称已执行。原固定心跳已改为仅监控这两条完整链，CPU／native资格及M/E/旧DEV全部关闭，保留真实交接标记；不新建对话或自动化。下一步由主完成四格固定汇总接线，并在真实启动后核验源和实际LR，最终处理全部TRAIN端点。工程提交不代表论文goal完成。

## 实际训练启动及四格只读接入

01:09:26UTC取证：LOW13545578于01:00:59UTC在gpu034、HIGH13545584于01:04:39在gpu035开始，两者root/batch/extern均RUNNING。实际started.json匹配本支job、原13434603、支持13453408及正确已成功CPU／native资格IDs；两支真实首两行更新是head4609／4610、native1／2，headlr3e−4，LOWnative1e−5／HIGHnative3e−4。不以这些初始更新推算最终性能或保证完成时间。

四格只读入口`compare_full_grammar_qv_lr_four_cell.py`及3组人工计数检查现已实施，不加载DNA／GFF正文／checkpoint／native模型；四份新Python／tests／SBATCH精确同步。固定比较各参数化高减低和差中差，输出所有六chr／两species／pooled／macro／worst的完整链与base PRF、背景、phase、nonprimary／其他isoform／duplicates，保留负数。检查共同reference、全部unsupported、22688分母、四份已完成对B／原DNA的逐碱基相等收据；复算PRF和聚合。Worst对每格各PRF取两物种较低值后作差；对有害背景率取各格较高值后作差，macro为两物种平均。未匹配primary／其他assessable isoform不解释成新基因或已证明的假阳性。历史E implicit native lr3e−4绑定其已核验同hyperparameters AdamW来源，另外三格需明确配置native_lr匹配。

首次四格传输因SSH banner超时未落地，01:08:51UTC检查提交被`Unable to open file`拒绝，未创建job；不是作业执行失败。随后精确同步并回读四份源后于01:10:14UTC提交CPU13545776，资源2CPU／4GiB／无GPU／5min。人工3组在0.070秒通过，但旧E source protocol被新入口误写为`full_grammar_joint_lora_plan.md`，实际为`full_grammar_joint_lora_fit.md`，接入guard正确拒绝。此job01:10:17–01:10:19UTC root/batchFAILED1:0、externCOMPLETED0:0，STATUS保留FAILED，尚无最终summary；不记零准确率。只修新汇总绑定的这一处实际名称并作一次同科学范围R1接入检查，不修改E/M／新fits／数据／权重／阈值。

原固定心跳5279字符新scope及全部字段已精确核对，实际在原watcher turn01a132c6-d4cc-7f90-92df-8a5cbf3655b9于01:00:53–01:01:06UTC追加运行，模型gpt-5.6-luna／low。它只查LOW/HIGH新链正常调度并DONT_NOTIFY，未向主发进度或终态交接；这是定时追加与单次检查验证，不宣称完整自动终态回传已验证。短CPU13545776由主消费关闭，不交给watcher。

R1接入CPU13545862于01:15:19–01:15:21UTC、2秒、root/batch/extern全部COMPLETED0:0，summary中3tests／0errors／0failures、旧E/M精确结构接入均成功，无DNA／真实reference GFF／checkpoint／native forward／更新／DEV／sealed。独立比较实际旧fit和新started metadata后，四格所有声明的共同科学配置字段相同。仅新汇总一处旧protocol名字纠正，原FAILED13545776和NOJOB提交记录保留；资格由主消费、watcher未发送，不能伪记自动交接。

四格最终只读作业实际13545908于01:17:27UTC（北京时间09:17）提交：`sbatch --parsable --dependency=afterok:13545582:13545588 sbatch/M28-FULL-GRAMMAR-QV-LR-FOUR-CELL.sbatch 13545581 13545587 13545582 13545588 13545862`。scontrol回执核对2CPU／4GiB／10min／无GPU／private／Requeue0、完整实际参数和两个有效新paired依赖，当前PENDING/Dependency。检查已由主核验且程序会再读资格收据，不为旧完成资格或E/M/B加controller依赖。没有新四格性能值，输出还未生成。

原同一心跳现在按一个完整批次QV-LR-FOUR-CELL-13545908等待两条新全链及汇总共21真实roots／63root-batch-extern记录，单fit／单front chain／单paired不是整批完成；完整成功或需处理失败才向原main一次真实交接后暂停同心跳。全部已消费资格关闭，不重复检查或回传；未创建新的监控对话或自动化。新6643字符完整批次prompt与ACTIVE／heartbeat类型／原ID／target／15min／通知偏好／created_at均已精确保存回读核对；新增汇总对象的定时追加及完整末端回传尚未观察，保留区别于此前已验证的原watcher定时单次检查。最终科学核验、完整四格结果解释和论文更新仍未完成，研究目标不缩减为提交成功。
