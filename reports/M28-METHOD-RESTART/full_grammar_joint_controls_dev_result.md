# 三个固定联合编码候选 DEV 与基线结果

2026-10-12 UTC。R1六个严格raw任务、三比较及家族13520481整链成功，30条实际root／batch／extern均COMPLETED0:0。Fresh×joint pooled完整CDS链F1为70.106738%，uniform为69.794721%，codon为64.578529%；三者两物种F1均高于固定Helixer、低于Tiberius和ANNEVO。完整评价成功不是论文goal达成，也不是独立泛化。三个候选和所有负差完整保留，DEV不用于后续调参或候选选择。

## 固定科学范围与实际终态

遵循[冻结方案](full_grammar_joint_controls_dev_plan.md)。候选源fit B13494895／C13496011／D13497914，各用事先固定的最终step_009216.pt；6144／7680未选择。原两条历史复用DEV仅拟南芥NC_003074.8／23459830bp和水稻NC_089041.1／29936421bp，caller只读取raw FASTA、双链全contig，不输入参考、标签、feature cache或观察模型，不作推理期更新。B/C/D所有完整TRAIN闭合先于DEV访问。E低秩适配和M全秩Q/V保持TRAIN-only，不加入此家族。

| 阶段 | 实际numeric作业 | UTC开始至结束 | Slurm耗时 |
| --- | --- | --- | --- |
| B AT raw | 13520474 | 10:14:52–10:44:36 | 29:44 |
| B rice raw | 13520473 | 10:14:52–10:51:54 | 37:02 |
| C AT raw | 13520482 | 10:14:56–10:44:36 | 29:40 |
| C rice raw | 13520477 | 10:14:56–10:51:30 | 36:34 |
| D AT raw | 13520483 | 10:14:56–10:44:36 | 29:40 |
| D rice raw | 13520479 | 10:14:56–10:51:53 | 36:57 |
| B comparison | 13520476 | 10:51:55–10:53:05 | 1:10 |
| C comparison | 13520478 | 10:51:31–10:53:05 | 1:34 |
| D comparison | 13520480 | 10:51:53–10:53:05 | 1:12 |
| 六方法family | 13520481 | 10:53:06–10:53:22 | 0:16 |

全部十roots及各batch／extern三步骤成功。每raw的STATUS／JOBID、started、外summary、input.fasta、内prediction/summary.json和full_grammar.gff3完整；末checkpoint／源fit／contig／两方向、b048 encoder、无参考／标签／更新及相应role-only或codon skyline身份核对。D的真实末64参数严格加载和raw读出相等由实际CPU评价完成，不在登录节点重新加载模型或解析大GFF。工程失败初始P、同步碰撞13519990及成功资格13520076均保留原记录，不改写为R1研究结果。

每候选输入、三原生输入和原FASTA在参考访问前直接相等；共同M26 CDS-assessable longest-eligible-CDS／transcript-ID tie保留AT5437、rice2291、总7728 primary。Unsupported、重叠和三个span>24576参考均保留。结构层不是额外筛选分母，nonstandard／partial基线预测不剔除。三份实际比较中的原生逐物种／pooled完整测量直接相等，ANNEVO／Tiberius／Helixer仍为原固定版本及六份预测，不重跑或更新。主独立复算所有方法的完整链及CDS-base计数P/R/F1、pooled加和、背景面积／链率、macro／worst及全部15组候选差值，结果一致。

## 完整链端点

AT／rice的primary分母分别5437／2291；TP／预测／参考列保留全部数，不把非primary链视为生物学错误。

| 方法 | AT TP／预测／参考 | Rice TP／预测／参考 | Pooled TP／预测／参考 | Pooled precision % | Pooled recall % | Pooled F1 % |
| --- | --- | --- | --- | ---: | ---: | ---: |
| fresh×joint（信息coarse） | 3887/5233/5437 | 1466/2310/2291 | 5353/7543/7728 | 70.966459 | 69.267598 | 70.106738 |
| uniform coarse | 3928/5278/5437 | 1427/2339/2291 | 5355/7617/7728 | 70.303269 | 69.293478 | 69.794721 |
| codon energy | 3755/5441/5437 | 1286/2443/2291 | 5041/7884/7728 | 63.939625 | 65.230331 | 64.578529 |
| ANNEVO | 3941/4898/5437 | 1703/2285/2291 | 5644/7183/7728 | 78.574412 | 73.033126 | 75.702502 |
| Tiberius | 3865/4984/5437 | 1576/2443/2291 | 5441/7427/7728 | 73.259728 | 70.406315 | 71.804685 |
| Helixer | 3673/5453/5437 | 1439/3086/2291 | 5112/8539/7728 | 59.866495 | 66.149068 | 62.851171 |

| 方法 | AT F1 % | Rice F1 % | Species macro F1 % | Worst-species F1 % |
| --- | ---: | ---: | ---: | ---: |
| fresh×joint（信息coarse） | 72.858482 | 63.725277 | 68.291879 | 63.725277 |
| uniform coarse | 73.317779 | 61.641469 | 67.479624 | 61.641469 |
| codon energy | 69.038426 | 54.330376 | 61.684401 | 54.330376 |
| ANNEVO | 76.265119 | 74.431818 | 75.348468 | 74.431818 |
| Tiberius | 74.177142 | 66.582172 | 70.379657 | 66.582172 |
| Helixer | 67.456382 | 53.524270 | 60.490326 | 53.524270 |

Fresh×joint减Helixer为AT+5.402100、rice+10.201007、pooled+7.255567 pp，TP净多214／27／241、pooled预测少996；本候选在两物种完整链precision、recall和F1均高于Helixer。不能将这一结果扩展为三候选均提高recall：uniform在rice比Helixer少12个TP，codon在rice少153、pooled少71。

Fresh×joint仍比Tiberius低AT1.318661、rice2.856894、pooled1.697947 pp（pooled少88 TP／多116预测），比ANNEVO低3.406637／10.706541／5.595763 pp（pooled少291 TP／多360预测）。三候选均未超过这两个较强比较者，不能只选择Helixer宣称竞争性目标完成。

## 两个预先固定的匹配控制

| 分支减fresh×joint | AT F1差 pp | Rice F1差 pp | Pooled micro差 pp | Macro差 pp | Worst差 pp |
| --- | ---: | ---: | ---: | ---: | ---: |
| uniform coarse | 0.459297 | -2.083808 | -0.312017 | -0.812256 | -2.083808 |
| codon energy | -3.820056 | -9.394901 | -5.528209 | -6.607478 | -9.394901 |

信息coarse减uniform的pooled F1为+0.312017 pp、AT−0.459297、rice+2.083808 pp；信息分支TP净少2、预测少74，pooled precision+0.663190 pp／recall−0.025880 pp。这个单fit追加阶段没有统一物种收益；TRAIN的pooled差为−0.342257 pp，DEV符号不同不能写成通道必要性或普遍增益。两分支都已先经历4608步信息coarse，不是从零全程移除通道的实验。

Codon减fresh×joint的pooled F1为−5.528209 pp、AT−3.820056／rice−9.394901，pooled TP少312而预测多341（AT−132／+208；rice−180／+133）。这是有效负结果，与原TRAIN负差方向一致，不作为换epoch、重训、删候选或按DEV选择后续科学方案的依据。内容项、增加参数、头更新及联合裁剪轨迹组成此干预；不从单个系数符号推出因果解释或否定coding-content建模一般有效。

## Base与背景不是完整链的替代端点

| 方法 | Pooled CDS-base F1 % | Wholly-background unique链 | Background span % | Background链／Mb |
| --- | ---: | ---: | ---: | ---: |
| fresh×joint（信息coarse） | 95.064575 | 145 | 0.921272 | 5.186640 |
| uniform coarse | 94.960091 | 154 | 0.956828 | 5.508570 |
| codon energy | 94.859006 | 196 | 0.962873 | 7.010907 |
| ANNEVO | 93.793505 | 146 | 1.156313 | 5.222410 |
| Tiberius | 94.098066 | 165 | 0.720517 | 5.902039 |
| Helixer | 92.031173 | 563 | 3.492437 | 20.138472 |

Fresh×joint的base F1为95.064575%，高于Tiberius和ANNEVO的94.098066%／93.793505%，但完整链F1仍更低，不能以base替代主端点。Fresh／uniform／codon背景链145／154／196；codon比fresh多51且背景span上升。Fresh虽比Tiberius少20条纯背景链，背景span却更高（0.921272%对0.720517%），面积与链数是不同端点。背景定义只相对于全部注释gene feature span之外区域，非独立生物学验证。

三个候选均无重复链，所有输出记录phase一致。非primary几何分别2190／2262／2843，其中精确匹配其他可评价isoform为435／416／360，余1755／1846／2483匹配neither；这些层不改变primary TP，不证明未匹配者为错误新基因。三个候选都未恢复保留的三条span>24576参考；该极小层不支持长gene能力或独立负泛化结论。

## 结论与接续

本轮完整六方法家族成功关闭，旧八／十方法表及原生测量不改。Fresh×joint的描述性优势仅相对于固定Helixer；uniform的物种权衡和codon负结果都支持“学习信号／合法结构接入不保证主端点改善”，不支持独立泛化、纯因果机制或论文已达标。历史复用两物种DEV、单fit、方法相关的上游注释暴露仍限制结论；不作clean zero-shot或公平速度优势主张，sealed／held-out／Setaria／新面板和公开发布仍未使用。

源为outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-COMPARISON-R1/{FRESH_JOINT24576,UNIFORM_COARSE,CODON_ENERGY}/summary.json与DEV-FAMILY-R1/{summary.json,comparison.md}，六raw保留独立目录。原watcher于11:05UTC真实交接P_R1（YES），主现已消费；同一心跳后续只等待原M整链，不再检查／交接P_R1。M本已按原TRAIN方案运行，DEV结果不改变其训练、冻结读出或家族成员。论文goal尚未完成。

