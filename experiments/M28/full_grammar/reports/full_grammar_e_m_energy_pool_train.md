# 原 TRAIN 双模型等权结构能量对照

2026-10-13。固定E/M等权组合的原TRAIN pooled完整链F1为73.857085%、CDS-base F1为95.513061%。它高于两parent及Tiberius，但仍低于ANNEVO和硬共识的完整链F1；相对硬共识恢复891个主TP和base覆盖，同时损失precision。整批实际完成，不是全面改善、独立泛化或论文目标完成。

两模型在同一DNA、方向、窗口ownership和完整基因语法上，逐角色先转FP64再取平均，重新求整条oriented染色体合法gene-set MAP。新结果前固定0.5／0.5，不搜索权重、不择优checkpoint或拟合温度／阈值。E/M及硬共识结果已知，因此是post-outcome TRAIN探索，不是独立验证或新增冻结DEV方法。

## 组合规则与解释

设两模型共享语法中的 gene-set 能量为 E_E(S;x) 和 E_M(S;x)，各自减去背景角色势。平均十五角色 logits 得到 E_pool(S;x) = 0.5 E_E(S;x) + 0.5 E_M(S;x)。空集能量仍为0。若各模型 P_j(S|x) 正比于 exp(E_j)，则 P_pool(S|x) 正比于 sqrt(P_E(S|x) P_M(S|x))；各自 logZ 对 S 为常数，求 MAP 不必先计算分区。这里平均的是能量，不是 decoded chain、局部 softmax、gene marginal 或真值正确率。

这是固定 logarithmic opinion pool 的应用，不将组合公式本身称作新颖贡献；一般方法已有[Heskes的先行工作](https://proceedings.neurips.cc/paper_files/paper/1997/hash/59f51fd6937412b7e56ded1ea2470c25-Abstract.html)。原能量尺度未经过跨模型正确率校准，等权并不保证性能更优或生物学校准。新的最优链可以由两支共享的局部证据组合产生，并不局限于旧输出的交集或并集。

原 whole-contig skyline 解码、FP64 per-base forward accumulation、null0、unit角色权重、实际 codon carry、GT／GC-AG、min-intron5 和 cap128 不变。两个方向独立，允许 opposite-strand overlap；不使用观察似然、标签、参考候选或 score gate。

## 输入范围与完整评价

只用原六条 TRAIN：AT NC_003070.9／NC_003071.7／NC_003076.8，rice NC_089043.1／NC_089044.1／NC_089046.1；长度分别30427671／19698289／26975502／27474823／23931887／28271460。E fit13499817、M fit13510324 均末head9216／native4608。输入为已完成整染色体后验流程中的双方向 FP32 势文件，E emissions13548914／16／18／20／22／24，M13548926／28／30／33／35／38；不使用新 LR 实验中间权重、不运行新 forward 或训练。

先确认十二个 completed 源身份、原 raw 收据和逐 base ownership 相同；每条原 FASTA 直接相等后才组合。实际数组必须为 N×15 FP32，源只读；平均结果临时 FP64，不写新完整数组。新 GFF 独立输出，原 DNA／weights／GFF／冻结结果不改写。资格用独立穷举完整 sets 验证平均能量、normalized logarithmic pool、MAP、空集、source 不变、row gauge、FP64 promotion、实际 carry 和 cap；另构造新的最佳 CDS 链，证明并非硬共识或旧链 union 选择。

所有六条新 raw 输出完成后才读取参考。沿用 M26 CDS-assessable longest eligible／tie transcript-ID 主转录本规则，保留全部22,688 primary references和 unsupported strata。导出重读核对唯一链数量、DNA／phase／cap及同 strand span 非重叠。重建六染色体、两物种、pooled PRF、macro／worst-species F1、CDS-base union、背景负担和其他 assessable isoform；参考数量、strata、base 分母和背景必须与既有 E/M、共识和 native 比较者相同。

比较固定 E、M、严格共识、ANNEVO、Tiberius、Helixer，不改任何比较者 operating point。按预测几何分为 matches both、E only、M only、neither parent emitted，再报告各类 exact-primary、exact-other-isoform和 neither-exact。neither-exact 不是生物学假基因。主要判断包含完整链 F1、召回和 base 覆盖是否共同改善，并报告每物种代价；不根据结果追选权重、阈值、预测池或 DEV 方法。

## 执行与完成边界

一个2CPU／4GiB／5min资格作业；六个独立2CPU／32GiB／1h raw CPU作业，依赖新资格；一个2CPU／8GiB／30min评价作业，afterok六个实际新 raw IDs。沿用 generanno／private-teodoro-gpu／no-requeue／无GPU。单条临时 FP64 势最多约3.40GiB，逐 strand 处理；源 mmap 不改写、不额外保存35GiB数组，六份新GFF与元数据声明512MiB总上限。无需已离开controller的旧ID依赖，只核验 completed 产物收据。

输出 `outputs/M28-FULL-GRAMMAR-E-M-ENERGY-POOL-TRAIN/<seqid>` 与 `evaluation`，资格输出独立 actual-job 子目录。STATUS／JOBID／started／非空summary和新GFF，以及 sacct root／batch／extern 完成0:0，才构成完整运行证据。科学负结果保留；确定性接线错误至多一次同范围修复、独立新输出保留失败，不取消作业或自动延长。该实验使用缓存 raw 势，不能将 CPU redecode 的成本当成包括两模型训练／native 编码的完整 pipeline 成本，亦不等于可交付的任意 raw-FASTA 调用器。

若需等待，只在原 watcher 与同一 heartbeat 中追加独立 POOLTRAIN 事件，LR 和 POSTTRAIN 原合同保持不变；完整新七项流程或需处理失败才交接，正常运行／有效依赖静默。不创建新 chat、automation 或 cron，不访问 DEV／sealed／held-out／新物种，不提交或公开发布。完整竞争性模型与 Nature Communications 及以上论文目标仍未达到。

## 资格接线修复

首资格13551061于05:29:39–05:30:09UTC运行30秒，root／batch FAILED1:0、extern0:0。四组人工检查中的平均完整set能量／log pool／MAP、source不变／FP64／gauge和实际carry／cap通过；新最佳链案例失败。原因是人工数据将内部CDS角色1–3错写成donor角色7–9，使两parent都选择空集。仅修该测试的角色列，不改组合公式、解码器或科研参数；首次失败产物保留，未读取研究DNA、参考、weights或源数组。修复后同资源重新运行一次资格，成功前不提交正式raw。

修复资格13551105于05:32:47–05:33:17UTC运行30秒，root／batch／extern全部COMPLETED0:0，STATUS／JOBID／started／summary齐全。四组检查全部通过，十二份角色数组source、原raw收据与完整双方向ownership验证通过；程序14.416444秒，MaxRSS978264K，未读取研究DNA、参考、checkpoint或native。主已核验消费，watcher NOT_SENT_BY_WATCHER，不补发资格交接。

## 实际运行接线

七个正式CPU作业于05:34:44–45UTC（北京时间13:34）接受。raw命令 `sbatch --parsable --dependency=afterok:13551105 sbatch/M28-FULL-GRAMMAR-E-M-ENERGY-POOL-RAW.sbatch <task> 13551105`；eval命令使用六个实际raw IDs的afterok，并按task传入同六ID与checks13551105。没有依赖已关闭旧source IDs。

| task | TRAIN染色体 | raw job |
| --- | --- | --- |
| 0 | NC_003070.9 | 13551143 |
| 1 | NC_003071.7 | 13551145 |
| 2 | NC_003076.8 | 13551146 |
| 3 | NC_089043.1 | 13551147 |
| 4 | NC_089044.1 | 13551148 |
| 5 | NC_089046.1 | 13551149 |

实际evaluation13551150，afterok13551143／45／46／47／48／49；task1前存在实际跳号，不按相邻数字推算。首检查六raw全部RUNNING且对应started／JOBID存在，evaluation有效Dependency；首raw与eval的controller资源、Requeue0、命令和依赖已核对。own expanded35加本批七项仍小于80，BeeGFS161TiB可用覆盖512MiB输出声明，private CPU内存与11月维护窗口允许本次短CPU作业。三个Python AST、三个SBATCH语法及正确相对路径的七份叶文件精确同步回读后才提交。尚无新性能结果，提交与运行不能替代完整评价。

另一个独立数值定位13551206读取原E全染色体角色数组，定位旧POSTTRAIN前后向mass失败；不改变本组合公式或解码，见[后验失败定位](full_grammar_gene_posterior_train_plan.md)。仅复用原heartbeat加入POOLTRAIN与MASSDIAG，原LR完整合同保留，旧后验actionable失败已由主消费，不重监旧失败批次。

同一heartbeat的实际保存11,265字符prompt／ID／类型／原target／15min／failed_runs_only／created_at已回读匹配，应用仅去掉末端换行；原LR合同正文不变。05:45:55–05:46:14UTC实际定时turn01a133cb-c6fb-7001-a292-edcccd6abb97追加到原watcher，完整新instructions与gpt-5.6-luna／low一致；明确查询原LR全部upstream／arrays／末端、新六raw／eval与massdiag，每个当前句柄都有行，运行或有效依赖静默。这只验证绑定、实际定时接线与当次正常检查，不证明未来上游缺失处置或完整终态回传。05:46UTC task1的plus方向实际完成解码，另一方向及其他任务仍运行；没有完整新评价或性能结论。

## 完成证据与完整链结果

七项正式作业及全部root／batch／extern共21条均COMPLETED0:0，全部STATUS／JOBID／started／summary、六份新GFF齐全；六raw完成05:50:19–05:57:53UTC，评价13551150于05:57:55–05:59:16UTC运行81秒。最慢raw1389秒，整批缓存势接线至评价结束1472秒。主核验所有raw源／.5-.5固定身份／双方向全contig／无参考、标签、更新及mutation字段、导出数量与phase数量，并独立复算六chr→species→pooled exact/base PRF、背景、macro／worst、全部六比较者差值和四类预测几何分区，与实际summary一致。全部22,688主参考及unsupported strata保留；没有访问DEV或sealed、挑权重、调阈值或训练。

| 方法 | AT F1 % | Rice F1 % | Pooled precision % | Recall % | F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| E_RANK4 | 75.315947 | 64.907807 | 73.553127 | 72.373061 | 72.958322 |
| M_DENSE | 74.414706 | 64.550265 | 73.452821 | 70.989069 | 72.199933 |
| ENERGY_POOL | 75.439052 | 68.413218 | 74.898033 | 72.844676 | 73.857085 |
| CONSENSUS | 76.240856 | 68.539708 | 81.348525 | 68.917489 | 74.618817 |
| ANNEVO | 76.830754 | 71.178600 | 78.739781 | 72.593441 | 75.541796 |
| Tiberius | 75.417268 | 63.576543 | 74.455290 | 70.940585 | 72.655456 |
| Helixer | 68.122369 | 50.294337 | 60.914669 | 66.515339 | 63.591926 |

组合为16527／22066／22688，pooled F1 73.857085%，比E高0.898763 pp、比M高1.657153 pp。两物种均高于两parent；AT75.439052%、rice68.413218%，macro71.926135%、worst68.413218%。相对Tiberius，AT仅数值高0.021784 pp、rice高4.836675 pp、pooled高1.201629 pp，不把微小AT差或TRAIN结果称为显著／独立胜出。相对ANNEVO仍低AT1.391702／rice2.765382／pooled1.684710 pp；相对硬共识低AT0.801804／rice0.126489／pooled0.761732 pp。因此完整链F1方面并未取代共识或最强既有比较者。

| TRAIN染色体 | TP / 预测 / 参考 | Precision % | Recall % | F1 % |
| --- | --- | ---: | ---: | ---: |
| NC_003070.9 | 5243/6902/7135 | 75.963489 | 73.482831 | 74.702572 |
| NC_003071.7 | 3102/3994/4293 | 77.666500 | 72.257163 | 74.864245 |
| NC_003076.8 | 4735/6042/6311 | 78.368090 | 75.027729 | 76.661540 |
| NC_089043.1 | 1180/1677/1638 | 70.363745 | 72.039072 | 71.191554 |
| NC_089044.1 | 1127/1657/1601 | 68.014484 | 70.393504 | 69.183548 |
| NC_089046.1 | 1140/1794/1710 | 63.545151 | 66.666667 | 65.068493 |

## 恢复召回与base覆盖的代价

相对硬共识，组合多输出2845条、主TP多891，recall提高3.927186 pp，precision降低6.450492 pp。CDS-base precision95.780994%、recall95.246622%、F1 95.513061%，比硬共识高5.123990 pp，恢复了大部分整链交集的base损失；同时完整链F1较低，并不是一个端点上的全面改善。相对E／M，主TP分别多107／421；base F1分别数值高0.077612／0.031177 pp，须与其precision／recall变化一并解释。

| 方法 | CDS-base F1 % | Wholly-background链 | Background span % |
| --- | ---: | ---: | ---: |
| E_RANK4 | 95.435449 | 454 | 1.090450 |
| M_DENSE | 95.481884 | 359 | 0.885714 |
| ENERGY_POOL | 95.513061 | 379 | 0.923550 |
| CONSENSUS | 90.389071 | 245 | 0.400726 |
| ANNEVO | 93.942191 | 383 | 1.190478 |
| Tiberius | 94.428270 | 454 | 0.751776 |
| Helixer | 92.439894 | 1614 | 3.613215 |

组合背景链379、span0.923550%。较E更低但较M更高，span较Tiberius和硬共识更高；背景变化没有统一支配所有比较者。背景及neither-exact都是当前参考定义的类别，不证明对应预测的生物学真假。

## 相对parent的预测几何

| 组合预测分区 | unique预测 | exact-primary | exact-other-assessable-isoform | neither-exact |
| --- | ---: | ---: | ---: | ---: |
| matches both parents | 19220 | 15636 | 1325 | 2259 |
| matches E only | 1209 | 518 | 72 | 619 |
| matches M only | 1175 | 342 | 72 | 761 |
| neither parent emitted | 462 | 31 | 5 | 426 |

四类精确重建22066预测、16527主TP及1474其他可评估isoform。组合保留硬共识全部15636个主TP，但少保留1条共同的非主参考／非可评估isoform预测。因此它并不是把交集直接加回旧union。相对共识的891额外主TP由E-only518、M-only342、两parent均未输出的完整几何31构成；主要收益来自保留单parent的已有正确链，而非大量生成新的正确几何。462条新几何仅31主参考／5其他可评估isoform精确匹配，426 neither-exact不得据此判为假基因，也不得在结果后剔除它们或选择新预测规则。这里的“新”只指相对旧输出的CDS坐标，不是生物学新基因发现。

## 成本与下一步

六GFF共30,422,068字节，最大batch MaxRSS15,118,988K，小于32GiB；六raw加评价按root allocation仅计一次共14,526 CPU-seconds，batch／extern不叠加。这些是缓存角色势重解码成本，不包含两模型训练、native重编码或完整独立raw-FASTA调用器；不能作为速度优势。源完整数组、DNA、weights和原GFF未改变，没有新整数组或checkpoint。

该固定探索达到了恢复召回和base覆盖的目的，但仍不满足领先ANNEVO或独立泛化的论文证据要求，也没有消除precision取舍。保留原权重、.5/.5规则和完整负比较，不追加结果驱动的权重／阈值选择或扩大冻结DEV。原LR四格继续；完整后验数值失败单独定位，不能当成组合性能失败或成功。研究目标仍是竞争性strict raw-FASTA模型及Nature Communications及以上论文，而非这一TRAIN组合实验完成。

主已完整消费POOLTRAIN；本轮watcher未完成自动交接，handoff=NOT_SENT_BY_WATCHER，不补发或伪记YES。06:00:55–06:02:19UTC实际turn01a133d9-8302-7511-adbb-3330723c0b1c虽查出六个剩余root完成，但只用sacct-X、不查batch／extern及完整产物，随后把它称作部分完成。主已补齐全21终态与产物，关闭该事件，并在同一heartbeat明确去掉-X、root全完成时同轮继续step与marker核验；真实回传修复仍待未关闭事件实际终态验证。原LR13545908及MASSDIAG13551206保持独立NO，同watcher／heartbeat／模型／15min，不新建对话或自动化。

