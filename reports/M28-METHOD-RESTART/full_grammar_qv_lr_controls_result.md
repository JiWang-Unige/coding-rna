# 原 TRAIN q v 参数化与学习率四格结果

2026-10-14（北京时间）。完整四格已完成并独立复算：rank4的native lr从1e−5升至3e−4，pooled完整链F1由65.928535%升至72.958322%；全秩则由72.199933%降至58.982365%。两项新设置均低于各自原设置，未产生更好的候选。原监控报告的产物ID疑点已由同一sacct行的JobID／JobIDRaw映射消解，不需改原始记录或重跑。工程成功、数值可信和科学竞争性分别判断。

## 固定范围和实际完成

沿用 [四格固定协议](full_grammar_qv_lr_controls_plan.md)，仅原两物种／六TRAIN／22688主参考，每格一fit、seed1，head4608→9216、native0→4608，head lr3e−4和其余科学设置不变。低lr rank4来源fit13545578／verify13545579／raw13545580_0–5／eval13545581／paired13545582；高lr全秩来源13545584／13545585／13545586_0–5／13545587／13545588。另两格原E高lr fit13499817／eval13499823／paired13499824、M低lr fit13510324／eval13510351／paired13510352均只读。最终汇总13545908，实际checks13545862。

主核对新两链和汇总全部21实际root／63条root-batch-extern COMPLETED0:0。新fit实际elapsed52260秒（rank4低lr）／50739秒（全秩高lr），无自动延长、补训或重提。十二raw的STATUS／JOBID／started／外层summary／内层raw caller逐一对应实际task、完整DNA长度、末step009216和sourcefit；均全contig双链、梯度禁用、正确16 adapter或8 q/v master恢复、无reference／label／cache／观察模型及0更新。两FIT／VERIFY／PAIRED及原DNA相等收据的真实身份、optimizer age、fixed fresh4608／4219、native/head lr、labels/context未改和全部检查完成字段均核对。不存在科学失败记0分、删除较差格或新DEV候选。

## Slurm 数组 ID 疑点

`sacct -j 13545580,13545586 --format=JobIDRaw,JobID,State,ExitCode --parsable2` 的同一记录提供以下映射；所有对应root/batch/extern均0:0。SBATCH的JOBID标记与Python started／summary均保存运行时 `SLURM_JOB_ID`，即numeric JobIDRaw，而不是数组展示名。

| task／TRAIN | LOW JobID | LOW JobIDRaw | HIGH JobID | HIGH JobIDRaw |
| --- | --- | --- | --- | --- |
| 0／NC_003070.9 | 13545580_0 | 13580469 | 13545586_0 | 13578916 |
| 1／NC_003071.7 | 13545580_1 | 13580852 | 13545586_1 | 13579084 |
| 2／NC_003076.8 | 13545580_2 | 13580859 | 13545586_2 | 13579273 |
| 3／NC_089043.1 | 13545580_3 | 13581088 | 13545586_3 | 13579838 |
| 4／NC_089044.1 | 13545580_4 | 13581125 | 13545586_4 | 13579908 |
| 5／NC_089046.1 | 13545580_5 | 13545580 | 13545586_5 | 13545586 |

十二项JOBID／started.job_id／summary.job_id与对应JobIDRaw全部逐字相等。不能通过单独查询numeric子ID返回空，判定这些产物没有调度记录；应从实际array master返回的JobID与JobIDRaw同一行核验。Slurm [sacct字段定义](https://slurm.schedmd.com/sacct.html) 和 [数组环境变量说明](https://slurm.schedmd.com/job_array.html) 区分这两种身份。保留原monitor疑点交接及所有产物，新增本解释，不篡改历史报告或真实失败13545776。

## 全部四格完整链结果

原低秩／全秩与低／高native lr对应 lora_low／lora_high／dense_low／dense_high；低1e−5、高3e−4。PRF为百分数，差值为百分点。Worst先在每格取较低物种F1再作差；不把最小物种差冒充worst差。负值完整保留。

| scope | rank4低lr F1 % | rank4高lr F1 % | 全秩低lr F1 % | 全秩高lr F1 % | rank4高−低 pp | 全秩高−低 pp | 差中差 pp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis | 70.415030 | 75.315947 | 74.414706 | 65.131747 | 4.900917 | -9.282959 | -14.183875 |
| Rice | 51.449340 | 64.907807 | 64.550265 | 37.972225 | 13.458467 | -26.578040 | -40.036507 |
| Pooled | 65.928535 | 72.958322 | 72.199933 | 58.982365 | 7.029787 | -13.217568 | -20.247354 |
| Macro | 60.932185 | 70.111877 | 69.482485 | 51.551986 | 9.179692 | -17.930499 | -27.110191 |
| Worst | 51.449340 | 64.907807 | 64.550265 | 37.972225 | 13.458467 | -26.578040 | -40.036507 |

| cell | primary TP | unique predictions | precision % | recall % | CDS-base F1 % | background span % | wholly-background chains |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| lora_low | 15720 | 25000 | 62.880000 | 69.287729 | 93.806742 | 1.882813 | 1092 |
| lora_high | 16420 | 22324 | 73.553127 | 72.373061 | 95.435449 | 1.090450 | 454 |
| dense_low | 16106 | 21927 | 73.452821 | 70.989069 | 95.481884 | 0.885714 | 359 |
| dense_high | 13412 | 22790 | 58.850373 | 59.114951 | 95.361038 | 0.844658 | 428 |

每格仍含全部22688主参考，所有unsupported／long-span／same-strand overlap层、isoform与background定义不变。rank4低lr和全秩高lr的CDS-base F1分别93.806742%／95.361038%，不能用高base分数掩盖65.928535%／58.982365%的完整链表现。所有新增设置和原E/M都完整保留，非主参考匹配不等于新基因或已证实的假阳性。

## 独立数值复算与解释

主读取四个实际completed评价JSON的原始计数，与最终 [96KB汇总](../../outputs/M28-FULL-GRAMMAR-QV-LR-FOUR-CELL/summary.json) 逐项比较；不调用producer summarize。独立重算36个原chr／species／pooled行的完整链及base PRF、六chr加总、phase／duplicate／其他isoform／nonprimary、全部background计数／分母／rates、每格macro／worst及全部11scope各measure的高减低和差中差，全部一致。原参考人口／unsupported strata、权重来源及四份对B／规范原DNA的访问前相等收据一致。没有重读DNA／参考GFF／checkpoint、执行native forward、再次评分或训练。零预测合法，未用结果筛样本。

学习率敏感性在两物种中方向相反：rank4高lr比低lr的AT／rice F1高4.900917／13.458467pp，全秩高lr却低9.282959／26.578040pp；pooled差中差−20.247354pp。在相同低lr下全秩pooled高于rank4 6.271397pp，而高lr下rank4高于全秩13.975957pp。因此原E／M对比不能孤立解释成rank4普遍更好，或学习率普遍提高／降低性能。不同参数化、BF16 forward舍入、共同clip及optimizer轨迹仍不等同，每格仅单fit／seed；这是固定TRAIN优化敏感性读出，不是纯rank因果、全范围最优lr或独立泛化证据。

原E高lr与M低lr仍为本次各参数化两种已测lr中较好的格，但未改变原竞争者：E pooled72.958322%、M72.199933%，仍低于ANNEVO75.541796%及固定Carbon-A两口径。此次不选择新checkpoint／阈值、扩充8／10／6冻结DEV家族或放行sealed／held-out／新panel，也不发布未发表工作。下一阶段模型改进应针对原完整链评分与竞争结构问题，而非把单fit学习率敏感性包装为论文目标已经达成。

## 真实回传核验

同一原watcher实际turn01a1362c-a6ee-75a3-84bc-3754814450e8在16:50:58UTC启动，原生模型仍为gpt-5.6-luna／low。其真实send_message_to_thread使用裸threadId `01a0332e-e0a9-7aa2-90f9-159405701833`、hostId `local`，省略model/thinking；消息实际到达本研究对话，随后原heartbeat保持全部字段并暂停。此为真实actionable交接，不是测试消息或成功性能通知。主消解数组映射误报、独立消费四格后关闭该事件，不补发第二次回传。暂停监控不改变原平台goal；原完整论文goal仍ACTIVE且未达成。
