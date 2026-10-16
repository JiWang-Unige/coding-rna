# 原 TRAIN 局部能量分歧惩罚试验

完整结果已复算：MINROLE pooled完整链F1 72.772784%、AT74.701405%／rice66.039857%，均低于固定mean-pool；净少537个primary匹配。工程执行完成，固定规则为TRAIN科学负结果，不是独立泛化或论文目标达成。原生raw-FASTA有限接线13572020也已完成并主消费，结果见[接线记录](full_grammar_minimum_centered_raw.md)；下方保留事前方案与完整执行结果。

2026-10-13。等权组合已提高原TRAIN完整链F1但仍低于ANNEVO／Carbon-A；后验排序在组合自身预测池内有效，却未解决全部漏链或证明较两parent更好的概率模型。已完成R2诊断中，E／M漏链约62.16%／61.23%为正能量参考被已选冲突链压过，34.51%／35.60%低于背景。下一项固定TRAIN探索改变角色竞争能量，不扫描后验阈值或改冻结DEV。

## 固定模型和科学问题

每base／role先用FP64计算各branch相对各自BG的能量`a=E_role−E_BG`、`b=M_role−M_BG`，再用`min(a,b)=(a+b)/2−abs(a−b)/2`。背景严格为0。这是固定系数1／2的局部分歧惩罚，不是min完整gene energy、两个MAP输出交集、gene probability平均或生物学错误概率。两支一致时恢复同一模型的背景差能量，独立row gauge不影响结果；完整路径能量不超过同一路径在任一parent的能量。它可能降低误发奖励，也可能过度惩罚长链、缩短或切碎基因；全部负结果和召回损失必须保留。

只用原两物种／六TRAIN，E13499817与M13510324末head9216／native4608及既有已核验35GiB势只读。语法、role集合、null0、cap128、actual codon carry、GT／GC-AG、min-intron、ownership、FP64 whole-contig MAP保持不变；没有训练、native forward、参考驱动选链、温度／阈值／系数搜索，也不读DEV／sealed／新panel。无第一创算法或泛化claim。

## 资格和正式试验顺序

先用三个合成测试组核验独立逐项最小值、全合法set枚举MAP与parent能量上界、source不变／独立gauge／同支退化、whole-path minimum区别及实际carry／exon cap。资格仅2CPU4GiB5min、private/generanno/noGPU/no-requeue；不读取研究DNA／势／参考／weights。通过后再完成原十二emission身份接线和六个CPU重解码，最后沿用完整原评价ruler及全部22,688主参考、unsupported strata、phase/base/background、nonprimary isoform、六chr／两species／pooled／macro／worst和既有E／M／equal-pool／hard-consensus／native基线；不以TRAIN表现预筛参考或只报告获胜物种。每个重新解码2CPU32GiB1h，评价2CPU8GiB30min，资源与规则在正式提交前固定，实际句柄见执行记录。既有完整mean-pool解码均在相同1h资源内完成。

独立正式输出拟定`outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-TRAIN`，资格`outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-CHECKS/<actualjob>`；不得覆盖任何原数组、checkpoint、冻结结果或失败目录。新增输出总声明256MiB，不保存新全长势。零预测是有效科学负结果，保留全部参考，不借空输出换阈值或增训。完成仍需全部实际root/batch/extern、STATUS/JOBID/started/summary/GFF和评价完整；cached重新解码不能称为新的whole-TRAIN raw-FASTA编码接入，后续原生入口资格另行区分。

本方案是在已授权原TRAIN探索内的下一项固定能量机制，不改变原论文goal或封存/DEV协议。两组正在运行的学习率对照不取消、不改参数或延时。外部等待只更新同一watcher和原heartbeat，不新增监控对话／自动化。资格及完整评价之前不宣称研究成功。

## 完成数学资格

13567792于11:55:17–11:56:03UTC运行46s，root／batch／extern全部COMPLETED0:0，MaxRSS935388K，程序15.902969s。三个合成测试组通过、failures/errors为0，原十二个emission来源均已按实际原fit／输入receipt／完整角色数组metadata核验接线。没有读取研究DNA／参考／势数组／checkpoint、native forward或optimizer update；资格summary明确whole-gene minimum未使用。原`equal_role_energy_pool`18行完整源在添加新独立函数前与远程逐字相等，添加没有改变旧mean函数的计算。

主已直接核验资格和全部Slurm终态，CLOSED/main-consumedYES/NOT_SENT_BY_WATCHER。随后接入该固定规则的独立输出driver和同ruler评价，沿用以上范围和资源；实际接线与正式作业记录如下，不把资格成功写成准确率改善。

## 完整运行入口接入

独立driver沿用已完成equal-pool实验的源数组身份、直接双branch原DNA相等、whole-contig skyline、literal DNA／phase／cap检查、GFF roundtrip、same-strand span非重叠、M26完整主参考／isoform选择、base／background评分和逐chr→species→pooled聚合。只把组合kernel换为已资格的局部minimum-centered规则，使用独立输出、明确新identity和GFF ID，不修改旧equal-pool脚本或结果。

比较者固定为原E／M、硬共识、equal-energy pool、ANNEVO／Tiberius／Helixer及原默认98,304bp Carbon-A的两个已完成固定口径；Carbon-A来源13503419，原mean评价13551150。没有重新运行native或选择Carbon-A阈值。除pooled／per-species外，保留macro／worst的同口径差值；所有比较者各chr的参考数／isoform数／unsupported strata／background和base分母必须一致。空新预测仍保留完整参考。

正式研究前新增一次2CPU4GiB5min/noGPU接口接线核验：用实际数学资格13567792和现有比较者JSON检查新driver能正确读取原方法名／完成身份／六scope及同参考人口，不重复三个数学测试，也不读DNA、真实参考、势数组或weights。通过后六个2CPU32GiB1h/noGPU raw任务并行，评价2CPU8GiB30min afterok全部六个actual ID，并记录新binding来源。新summary和source_comparator_jobs独立保存；新输出大小在原256MiB声明内核验。数学资格、接口资格、正式完整评价和科学PASS分别记录，不互相替代。

## 接线资格与实际提交

新driver AST、三个SBATCH语法及四叶精确同步回读一致。SSH曾多次banner timeout，未把连接失败当成作业失败、没有重提／取消现有训练或修改SSH全局配置；随后连接恢复并完成最后版本回读，才进行提交。

接口资格13569064实际12:21:34–12:22:19UTC／45s，root／batch／extern全部COMPLETED0:0，MaxRSS833056K，程序0.724116s。九个固定比较者的方法名、完成来源、六chr／两species／全22,688参考人口及base／background／unsupported元数据一致；未重复数学测试或读研究DNA／真实参考／角色数组／weights。资格由主核验消费CLOSED/main-consumedYES/NOT_SENT_BY_WATCHER。

六个raw实际提交12:25:55–12:25:56UTC，均afterok:13569064，命令参数为task／数学13567792／binding13569064。ID来自逐条提交返回值，不凭相邻数字推算。

| task | 原TRAIN | actual raw job |
| --- | --- | --- |
| 0 | NC_003070.9 | 13569206 |
| 1 | NC_003071.7 | 13569207 |
| 2 | NC_003076.8 | 13569208 |
| 3 | NC_089043.1 | 13569209 |
| 4 | NC_089044.1 | 13569210 |
| 5 | NC_089046.1 | 13569211 |

评价actual13569505于12:27:09UTC提交，afterok上面六actual raw，参数顺序为数学13567792／binding13569064／task0–5六ID。controller核对七项资源、命令及依赖均与计划一致：raw2CPU32GiB1h、eval2CPU8GiB30min、全部private/generanno/noGPU/Requeue0。12:28UTC前五raw在gpu034 RUNNING，task5正常Resources等待，eval有效Dependency。启动前expanded own31，新增七项仍低于80，目标输出不存在、存储约160TiB可用、时限内无相关维护reservation。运行不是结果。

原LR两fit13545578／13545584继续运行，未改变参数或冻结DEV。等待复用同一heartbeat和原watcher，新增MINROLE独立末端13569505与六上游的7root／21条终态及产物／结构白名单合同。数学与binding资格已主消费，不补发交接；七正式阶段成功并齐全后才解释完整TRAIN结果，失败或科学负结果均保留。

新11,696字符MINROLE+LR心跳的ID／类型／原目标／15min／failed_runs_only／created_at及完整prompt精确回读一致，原6,642字符LR合同原文保留。实际新scope定时turn01a13547-4bdd-7963-b4a5-18fcbf05dc93于12:40:27–12:40:50UTC追加同一原watcher，完整instructions精确匹配，实际gpt-5.6-luna／low由该turn_context核对；exit0调度命令明确检查新六raw／13569505及原LR完整句柄。它没有读取科学值或发送正常进度交接。12:41UTC controller确认task5也已在gpu035运行，前五仍运行，评价有效Dependency；前五实际started身份／资格／binding／无参考／无DEV字段由主核对，task5初始加载时尚无started不误判失败。当时仅配置、首次定时检查已验证；完整结果及真实回传随后完成，见下。

## 完整结果和独立复算

实际六raw13569206／13569207／13569208／13569209／13569210／13569211及评价13569505的7root／21条root-batch-extern全部COMPLETED0:0；各raw实际elapsed1261／1132／1261／1218／940／616秒。评价12:51:09–12:55:29UTC／260秒，MaxRSS1683336K，程序230.489714秒。STATUS／JOBID／started／summary、六GFF及完整结构白名单已核验。评价保存前新产物32,125,754bytes，final JSON22,336bytes，在256MiB内。这里是已完成势的CPU重解码与评价，不含native fresh编码或两parent训练成本。

主在完整评价JSON上独立使用原始计数复算precision=TP/P、recall=TP/R、F1=2TP/(P+R)，不调用producer aggregate；核对六chr→两species→pooled的完整链、CDS-base、background全部可加计数和rates、macro／worst，逐chr预测／phase／duplicate与novelty互斥分组、全部九比较者同参考／base／background分母及全部pooled／species／macro／worst差值，全部相等。比较者来自各实际completed eval JSON（13499823／13510351／13461307／13550067／13551150／13503419），不是从新delta倒推旧F1。没有再次读取真实DNA、参考GFF、势、weights或执行新评分；旧数学、接口和literal DNA／phase资格不重复跑。主要JSON为 `outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-TRAIN/evaluation/summary.json`，原全部参考人口与unsupported strata保留。

| TRAIN scope | exact TP | unique predictions | all primary refs | chain F1 % |
| --- | ---: | ---: | ---: | ---: |
| NC_003070.9 | 5131 | 6710 | 7135 | 74.120621 |
| NC_003071.7 | 3019 | 3839 | 4293 | 74.249877 |
| NC_003076.8 | 4609 | 5872 | 6311 | 75.662809 |
| NC_089043.1 | 1107 | 1590 | 1638 | 68.587361 |
| NC_089044.1 | 1067 | 1586 | 1601 | 66.959523 |
| NC_089046.1 | 1057 | 1660 | 1710 | 62.729970 |
| Arabidopsis | 12759 | 16421 | 17739 | 74.701405 |
| Rice | 3231 | 4836 | 4949 | 66.039857 |
| Pooled | 15990 | 21257 | 22688 | 72.772784 |

Pooled precision75.222280%、recall70.477786%、CDS-base F195.481537%；macro70.370631%、worst/rice66.039857%。两species CDS-base F1分别96.942481%／90.788961%。全部21,257预测phase检查，duplicate0；其他assessable isoform非primary匹配1539（AT1405／rice134）。Background的81219532bp分母不变：predicted span545183bp／FPR0.671246%，predicted CDS291285bp，wholly-background264（AT43／rice221）。这些是benchmark reference-relative量，不直接等同生物学错基因率。

| unchanged comparator | AT chain F1 % | rice chain F1 % | pooled chain F1 % | MINROLE minus comparator pp |
| --- | ---: | ---: | ---: | ---: |
| E_RANK4 | 75.315947 | 64.907807 | 72.958322 | -0.185538 |
| M_DENSE | 74.414706 | 64.550265 | 72.199933 | 0.572851 |
| ANNEVO | 76.830754 | 71.178600 | 75.541796 | -2.769012 |
| Tiberius | 75.417268 | 63.576543 | 72.655456 | 0.117328 |
| Helixer | 68.122369 | 50.294337 | 63.591926 | 9.180858 |
| CONSENSUS | 76.240856 | 68.539708 | 74.618817 | -1.846033 |
| EQUAL_ENERGY_POOL | 75.439052 | 68.413218 | 73.857085 | -1.084301 |
| CarbonA_confidence_ge_0.1_or_NA | 78.957727 | 75.258517 | 78.171555 | -5.398771 |
| CarbonA_unfiltered | 76.745109 | 73.530032 | 76.063450 | -3.290665 |

相对mean-pool，预测净少809、primary TP净少537、nonprimary记录净少272；precision仅+0.324247pp，recall−2.366890pp，pooled F1−1.084301pp。AT／rice F1分别−0.737647／−2.373361pp；background span少204920bp、wholly-background少115不能补偿完整链召回损失。MINROLE低于ANNEVO／两Carbon-A口径／hard consensus，pooled还低于E；虽然高于M、Tiberius和Helixer的pooled值，不能只挑这些比较者宣布获胜。该固定局部分歧规则在两species均未改善mean-pool，保留为TRAIN机制负结果，不调系数、筛阈值或改变冻结DEV以挽救读数。

分组复算：matches-both18592／primary15367／other-iso1311／neither1914；E-only420／198／46／176；M-only1486／333／100／1053；neither-parent759／92／82／585。原hard consensus的15,636 primary中有269未被新MAP保留；新增的759只是相对两parent已发出的CDS geometry，不是发现新基因。92 primary、82 other iso、585 neither-exact都报告，不把较多新geometry包装为质量改善或新颖性claim。

## 真实监控交接和后续

原watcher实际turn01a13555-07ae-7241-b307-d8e9bd9ca1b8于12:55:27–12:57:26UTC在同一原对话运行（gpt-5.6-luna／low），完成7root／21终态及产物白名单后一次真实发送原main并暂停同heartbeat。主已消费科学值，MINROLE CLOSED/main-consumedYES/watcher-real-handoffYES；不重复发送，未把watcher结构核验当成科学PASS。

独立严格raw-FASTA入口的固定有限native资格13572020已于2026-10-13 13:15:18–13:17:26UTC完成，三条Slurm记录0:0及七项产物齐全；两个原TRAIN视图和14bp短输入共六方向独立fresh重编码／MAP／GFF／DNA／phase接线通过，详见[接线结果](full_grammar_minimum_centered_raw.md)。它不替代原六TRAIN fresh重编码或独立准确率，也不改变MINROLE的负TRAIN结论。主直接消费MINRAW，真实watcher交接为NOT_SENT，不补发；MINROLE真实交接YES保留。

后续完整LR四格13545908亦已真实交接并主消费，见[四格结果](full_grammar_qv_lr_controls_result.md)，不再等待两fit或重复检查这些关闭事件。同一原watcher／heartbeat当前PAUSED，模型、15min及原通知策略保留。原论文goal仍ACTIVE且未完成，冻结DEV families8／10／6、sealed／held-out／新panel和对外发布边界不变。本负结果已同步到现有Results／Methods和组合稿；没有增训、阈值更改或新科学成功claim。

