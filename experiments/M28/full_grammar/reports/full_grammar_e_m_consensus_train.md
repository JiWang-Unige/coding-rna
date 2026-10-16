# 原 TRAIN 双模型完整 CDS 共识对照

2026-10-13。固定 E/M 完整CDS交集的原TRAIN pooled F1为74.618817%，两个物种均高于两parent及Tiberius，但仍低于ANNEVO；召回和CDS-base覆盖下降，不能称全面改善。规则只看seqid、strand和全部有序CDS区间，不看参考、gain、posterior或已知正确率，不平均跨模型未校准的分数。由于 E/M 单模型结果已经知道，这是post-outcome探索，不是预DEV注册、独立验证或一个已交付的新raw-FASTA调用器。

## 固定范围与评价

两source均末head9216／native4608：E fit13499817／eval13499823，M fit13510324／eval13510351，只使用原六条 TRAIN 的已完成strict raw输出。取预测key的精确交集，使用E的原record保留相同DNA／phase／geometry；它是每支原合法同strand非重叠set的子集。不同方向仍可重叠。近似坐标、局部边界、得分阈值、择优参考命中和union冲突选择均不纳入这个固定规则。

source完成身份／final权重配置／原raw收据先核对；每染色体两支FASTA直接相等、原预测数量／唯一性／DNA／phase成立，并在读取参考前完成交集。原M26 CDS-assessable primary规则保留全部22,688参考及unsupported strata，原E/M的完整链、CDS-base、背景、重复和其他可评估isoform指标全部重建。用固定native评价13461307及相同参考、背景与strata比较ANNEVO、Tiberius、Helixer，不跑新基线或改operating point。

报告六染色体、两物种、pooled count-derived PRF、macro／worst-species F1、CDS-base union、背景负担及phase数量。另分common、E-only、M-only预测的exact-primary、exact-other-assessable-isoform和neither-exact counts；neither不是已证伪的gene。每物种都保留精度与召回损益，不把precision上涨直接称作完整链改善，不用pooled掩盖物种退化。真正支持继续研究共识的证据应包含F1而不只是更少输出或更高precision；负结果原样保留，不据结果补阈值、扩大DEV方法或改规则。

## 执行与结果边界

一个2CPU／8GiB／20min纯CPU Slurm作业，private-teodoro-gpu／generanno／no-requeue，输出独立 `outputs/M28-FULL-GRAMMAR-E-M-CONSENSUS-TRAIN` started／summary／STATUS／JOBID，声明小于2MiB。复用已完成old sources，不添加已经离开controller的old-ID afterok。不读取checkpoint、native features或新角色数组；不改旧GFF、不输出新GFF、不追加训练或访问DEV／sealed。直接在真实source上验证指标重建，不增加只复述set交集的人工测试或第二套评价框架。确定性接线错误至多一次同范围修复并保留失败；科学负结果不重选规则。

本次短作业由main核验并消费，不重复新建监控。若确需长等待，只复用原watcher与同一heartbeat，不改变原LR及POSTTRAIN事件。完整研究目标仍是竞争性raw-FASTA模型、独立适用性证据和论文，而不是交集算术通过或TRAIN单点胜出。

## 实际运行记录

实际13550067于04:46:08UTC（北京时间12:46）提交，命令 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-E-M-CONSENSUS-TRAIN.sbatch`。04:46:09–04:47:35UTC、86秒，root／batch／extern全部COMPLETED0:0，STATUS／JOBID／started／summary和六染色体完成标记齐全。2CPU／8GiB／20min上限／无GPU／private／Requeue0；程序实耗69.599487秒。AST／SBATCH语法及三份新源精确同步回读后才运行，执行前own expanded44加本作业仍低于80，维护仅11月CPU预留。主核验真实终态、结构身份、全部source／两输入DNA先相等／原预测phase与计数／完整parent指标重建，并独立复算六chr→species→pooled counts、PRF、base和背景、macro／worst、全部相对E/M及native差，以及三类agreement分区，全部一致。主已消费，watcher handoff=NOT_SENT_BY_WATCHER，不补发资格／短job交接、不新建自动化。

实际summary18,476字节、batch MaxRSS1,724,020K，小于声明；172 CPU allocation秒只计root一次，batch／extern不叠加。选择器耗时不包含两parent的训练和raw编码／解码，不能作为完整ensemble pipeline的速度或效率优势。

## 完整链结果与召回代价

| 方法 | AT F1 % | Rice F1 % | Pooled precision % | Recall % | F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| E_RANK4 | 75.315947 | 64.907807 | 73.553127 | 72.373061 | 72.958322 |
| M_DENSE | 74.414706 | 64.550265 | 73.452821 | 70.989069 | 72.199933 |
| CONSENSUS | 76.240856 | 68.539708 | 81.348525 | 68.917489 | 74.618817 |
| ANNEVO | 76.830754 | 71.178600 | 78.739781 | 72.593441 | 75.541796 |
| Tiberius | 75.417268 | 63.576543 | 74.455290 | 70.940585 | 72.655456 |
| Helixer | 68.122369 | 50.294337 | 60.914669 | 66.515339 | 63.591926 |

共识pooled为15636／19221／22688，完整链F1 74.618817%，比E高1.660495 pp、比M高2.418884 pp；两个物种都高于两个parent，也高于固定Tiberius（AT+0.823588／rice+4.963164／pooled+1.963360 pp）。但比ANNEVO低AT0.589898／rice2.638892／pooled0.922979 pp。Species macro为72.390282%、worst68.539708%；不因对Tiberius的TRAIN改善而宣称优于最强比较者或泛化成立。

相对E，预测减少3103、主TP减少784，precision升7.795398 pp而recall降3.455571 pp；相对M，预测减少2706、主TP减少470，precision升7.895704 pp而recall降2.071580 pp。完整链F1改善与召回损失同时存在，不把该结果写成无代价的错误修复。

| TRAIN染色体 | TP / 预测 / 参考 | Precision % | Recall % | F1 % |
| --- | --- | ---: | ---: | ---: |
| NC_003070.9 | 5063/6223/7135 | 81.359473 | 70.960056 | 75.804761 |
| NC_003071.7 | 2982/3605/4293 | 82.718447 | 69.461915 | 75.512788 |
| NC_003076.8 | 4566/5515/6311 | 82.792384 | 72.349865 | 77.219685 |
| NC_089043.1 | 1037/1280/1638 | 81.015625 | 63.308913 | 71.076080 |
| NC_089044.1 | 989/1276/1601 | 77.507837 | 61.773891 | 68.752172 |
| NC_089046.1 | 999/1322/1710 | 75.567322 | 58.421053 | 65.897098 |

## Base与背景结果

| 方法 | CDS-base F1 % | Wholly-background链 | Background span % |
| --- | ---: | ---: | ---: |
| E_RANK4 | 95.435449 | 454 | 1.090450 |
| M_DENSE | 95.481884 | 359 | 0.885714 |
| CONSENSUS | 90.389071 | 245 | 0.400726 |
| ANNEVO | 93.942191 | 383 | 1.190478 |
| Tiberius | 94.428270 | 454 | 0.751776 |
| Helixer | 92.439894 | 1614 | 3.613215 |

共识base precision97.285395%、recall84.405752%、F1 90.389071%，比E低5.046378 pp、比M低5.092813 pp。整链完全一致规则会同时排除含正确coding片段的不一致预测；这与base覆盖下降相符，但不证明每条被删记录的生物学真假。背景span降低和245条纯背景预测也只是同参考定义的负担变化，不单独作为全任务改善。

## 全预测 agreement 分区

| 分区 | unique预测 | exact-primary | exact-other-assessable-isoform | neither-exact |
| --- | ---: | ---: | ---: | ---: |
| common | 19221 | 15636 | 1325 | 2260 |
| E-only | 3103 | 784 | 124 | 2195 |
| M-only | 2706 | 470 | 136 | 2100 |

每行三类匹配数精确重建该分区预测数；common与E-only／M-only分别重建原E／M的预测和TP。没有重复common链，19221条原DNA／phase合法记录保留，原GFF和checkpoint均未改写，没有新GFF、forward、拟合、阈值、DEV或sealed访问。这里的neither是参考相对精确几何类别，不是已证伪基因。

## 对下一阶段的影响

严格全链共识确有原TRAIN F1收益且背景更低，但仍输ANNEVO，同时损失召回和CDS-base覆盖。因此保留它为需要独立验证的探索性对照，不将其直接放入已冻结DEV或取代单模型主结果。下一轮值得检验能否通过共享语法下的固定结构能量组合保留更多不一致链的coding证据；这只是由当前代价提出的假设，尚未实施或验证。原LR四格和完整后验／漏检批次继续其原协议，不据此改参数、范围或停止条件。完整raw-FASTA竞争性模型、独立适用性和Nature Communications及以上论文目标仍未达到。
