# 原 TRAIN 完整链比较结果

2026-10-09。原六条 TRAIN 的18组固定原生预测和共同尺评价已完成。原4608步模型完整链F1为45.590%，ANNEVO为75.542%、Tiberius为72.655%、Helixer为63.592%。差距在两个物种均存在，水稻更大；目前不能声称模型已达到成熟比较器的完整结构竞争力，也不能把问题仅归于跨数据范围泛化。

## 同一输入和共同尺

三个原生工具直接复用原4608步推理的六份大写ACGTN FASTA，共156779632bp。操作点、源码与权重沿用已固定的原生DEV比较，不重新选模型或参数。推理无真实参考、标签、蛋白或RNA输入。全部预测完成后才在CPU评价打开原TRAIN参考，按M26 CDS-assessable、最长eligible CDS、transcript ID tie-break选primary；分母拟南芥17739、水稻4949、合计22688，困难和unsupported参考保留。

逐染色体参考数、eligible isoform数、unsupported分层、参考CDS union碱基和背景碱基均与原4608步控制匹配。背景是双链所有注释gene feature span以外的区域；坐标exact-key包含seqid、strand和有序CDS区间，不含phase，phase诊断另列。CDS-base指标为不区分strand的union，不能代替完整结构指标。输入、参数和评价实现见[native_train_baseline_plan.md](native_train_baseline_plan.md)。

## 完整链准确率

表内F1为百分数。pooled按TP、唯一预测链数和参考数合并后计算，不是染色体F1均值；macro是两物种F1的等权平均，worst为其中较低者。

| 方法 | 拟南芥F1 | 水稻F1 | pooled F1 | macro F1 | worst F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原4608步 | 50.487 | 28.278 | 45.590 | 39.383 | 28.278 |
| ANNEVO | 76.831 | 71.179 | 75.542 | 74.005 | 71.179 |
| Tiberius | 75.417 | 63.577 | 72.655 | 69.497 | 63.577 |
| Helixer | 68.122 | 50.294 | 63.592 | 59.208 | 50.294 |

| 方法 | pooled TP | 唯一预测链 | pooled precision % | pooled recall % |
| --- | ---: | ---: | ---: | ---: |
| 原4608步 | 9966 | 21032 | 47.385 | 43.926 |
| ANNEVO | 16470 | 20917 | 78.740 | 72.593 |
| Tiberius | 16095 | 21617 | 74.455 | 70.941 |
| Helixer | 15091 | 24774 | 60.915 | 66.515 |

三工具相对原模型的pooled完整链F1差值分别为29.952、27.065和18.002个百分点；水稻差值为42.901、35.299和22.016个百分点。全部方法的重复链记录为0。其他eligible isoform但非primary的精确命中分别为原模型748、ANNEVO1498、Tiberius1552、Helixer1339；这些命中未充作primary TP，非primary命中也不直接等于生物学假基因。

同一已完成比较的[分物种及pooled静态图](figures/initial_train_benchmark/initial_train_complete_chain.png)采用统一0–100%刻度，12个F1值直接由TP／预测／参考计数重算，参考元数据及聚合已核对。图仅包含初始末4608模型与三个固定工具，不加入仍在运行的续训或抽样支路，不将未完成结果画为失败。[矢量图及图注](figures/initial_train_benchmark/caption.md)保留精确计数、两份原始summary和无GPU绘图任务13464155的完成证据；它仍是TRAIN诊断，不是独立测试或论文成功。

## 碱基覆盖与背景

| 方法 | pooled CDS-base F1 % | CDS-base precision % | CDS-base recall % | 背景span FPR % | 完全位于背景的唯一链 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原4608步 | 90.360 | 97.208 | 84.414 | 0.618 | 162 |
| ANNEVO | 93.942 | 95.523 | 92.413 | 1.190 | 383 |
| Tiberius | 94.428 | 95.968 | 92.937 | 0.752 | 454 |
| Helixer | 92.440 | 89.818 | 95.219 | 3.613 | 1614 |

原模型CDS碱基precision较高、背景span较少，但完整链precision和recall均落后。该组合说明碱基级保守性或低背景不能代替准确的整条CDS结构；不能仅凭这些汇总定位到目标函数、上下文或某个结构模块。水稻原模型CDS-base F1为80.851%，三个基线分别为89.822%、88.721%、82.156%；拟南芥分别为原模型93.085%、ANNEVO95.240%、Tiberius96.216%、Helixer95.996%。

## 完成证据和后续

数组13461306的0–17全部实际root／batch／extern均COMPLETED0:0。ANNEVO六格15:10:49 UTC开始，最后结束15:13:42；Tiberius最后结束15:23:48；Helixer最后结束15:31:19。CPU评价13461307在15:31:19–15:32:40 UTC运行81秒，全部Slurm记录COMPLETED0:0，STATUS／JOBID／非空summary齐全。原固定监控实际完成结构核验并回传主对话，保留其他五条流程；这不是提交或GPU探测代替结果。

原生评价保存在 `outputs/M28-NATIVE-TRAIN-BASELINE-EVAL/summary.json`，原模型控制为 `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC/summary.json`。主对话独立用TP／预测／参考计数重算全部逐物种和pooled PRF、macro及worst，与保存值一致。逐参考诊断另以Slurm13461637复现原生唯一链命中，比较基线能命中而模型遗漏的链及共同遗漏，方案见[native_train_reference_diagnostic_plan.md](native_train_reference_diagnostic_plan.md)。

本比较是原TRAIN描述性结果，不是独立测试或目标函数因果实验。比较器上游预训练／监督暴露依赖具体工具，不能统一称为未见或清洁zero-shot。原有六条训练分支及四格汇总继续按原依赖执行；不改正在运行的参数、不据此在DEV挑参数、不改变冻结八／十方法评价家族、不放行封存数据。原论文goal仍未达到。

## 完成运行的资源占用

2026-10-09追加原4608模型和三个固定原生工具的24个完成推理任务的执行记录。输入仍为同六份FASTA、156779632bp；不重新推理、读参考／预测正文或改变任何参数。该记录只描述本次配置的端到端资源分配，不是公平速度benchmark，也不是资源费用限制。

原模型为13443323_0–5，ANNEVO／Tiberius／Helixer分别为13461306_0–5／6–11／12–17。当前Slurm accounting的全部24个root及对应batch／extern均COMPLETED0:0，每个root实际分配1 RTX3090。按root的ElapsedRaw合计单GPU秒数，只计一次，不能将batch／extern重复相加；GPU分配小时除以3600，CPU分配小时另乘实际AllocCPUS。前者不是GPU活跃利用率，后者不是TotalCPU或实际CPU使用时间。Slurm将ElapsedRaw定义为运行秒数，MaxRSS是step内任务的最大内存报告值，不是全部任务同时使用的内存之和。[Slurm accounting字段说明](https://slurm.schedmd.com/sacct.html)

| 完成方法 | 每任务CPU／RAM分配 | 六任务单GPU秒数合计 | GPU分配小时 | CPU分配小时 | 单任务wall范围秒 | 最高batch MaxRSS GiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 原4608步 | 2／16GiB | 12006 | 3.335 | 6.670 | 1469–2294 | 7.601 |
| ANNEVO | 8／32GiB | 890 | 0.247 | 1.978 | 138–173 | 3.967 |
| Tiberius | 8／32GiB | 3629 | 1.008 | 8.064 | 578–639 | 9.414 |
| Helixer | 8／32GiB | 2604 | 0.723 | 5.787 | 376–454 | 2.251 |

MaxRSS只取六个batch报告值中的最大值，K值按KiB换算GiB；它不是VRAM、整套工作流内存之和或统一的全进程峰值测量。原模型运行于private-teodoro-gpu的gpu034／035，基线运行于shared-gpu的gpu021／035；GPU型号相同，但CPU分配、节点并发、软件／容器、线程与batch配置及计时边界不同。基线wrapper包含运行资产核对、GPU probe、原生工具和GFF完成检查，不能把Slurm wall解释为纯神经forward。数字也不含训练、共同尺CPU评价、排队或监控等待。

以下为相同染色体的root ElapsedRaw原始秒数，用于重算上表，不按序列长度外推，也不计算或宣称公平加速比。

| 原TRAIN seqid | 原4608步 | ANNEVO | Tiberius | Helixer |
| --- | ---: | ---: | ---: | ---: |
| NC_003070.9 | 2294 | 173 | 639 | 452 |
| NC_003071.7 | 1469 | 157 | 606 | 376 |
| NC_003076.8 | 2034 | 139 | 610 | 425 |
| NC_089043.1 | 2102 | 141 | 608 | 454 |
| NC_089044.1 | 1929 | 138 | 588 | 446 |
| NC_089046.1 | 2178 | 142 | 578 | 451 |

来源命令为 `TZ=UTC sacct -j 13443323,13461306 --format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS,AllocTRES%130,NodeList,MaxRSS --parsable2`；task映射依据实际 `run_native_train_baseline.task_scope`，六份原模型输入对应关系沿用固定TRAIN合同。解析时要求完整且唯一的6＋18个root、每个root的三个成功终态及实际GPU／CPU分配，再按原seqid顺序聚合；原模型12006秒与已保存raw收据一致。原模型内部encoding、joint decoding及外层入口的不同计时边界继续见[原TRAIN模型结果](full_grammar_raw_train_result.md)，不与本表或基线wrapper的不同计时范围混算。

本次配置下原模型累计GPU分配时长更高，现有证据不能支持计算效率优势；这不会改变或延迟已授权的并行研究。若后续完整链性能或独立适用性达到论文要求，仍需随相应最终模型报告实测计算代价，而不能把当前表替换成未经测量的速度主张。当前表没有新增DEV候选、科学PASS阈值、训练／数据放行或公开发布。
