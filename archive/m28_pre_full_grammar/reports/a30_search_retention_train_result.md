# A30 冻结搜索保留损失未通过候选门槛

2026-10-06。A30已完成原4608逻辑更新、唯一末004608的自然TRAIN读出及实际CPU核验，但原五项候选支持门槛为0/5 PASS。因此关闭本次固定seed0、margin1／系数1方案；条件组合A32未执行，不记为A32失败，也不换checkpoint、系数、阈值或seed救回结果。

## 实际候选结果

1491个原TRAIN窗口，原标签和A3自然池均重新核验。以下为候选/window链计数，不是最终gene TP、F1或独立泛化。

| 物种 | A3正确链 | A30正确链 | A3 known nonmatch | A30 known nonmatch | A3正确≥5 exon | A30正确≥5 exon |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | 1751 | 1750 | 6305 | 7894 | 604 | 605 |
| 水稻 | 629 | 612 | 6261 | 7029 | 184 | 177 |
| 合计 | 2380 | 2362 | 12566 | 14923 | 788 | 782 |

两个物种的正确链均未严格增加，合计正确长链未严格增加，两个物种的known nonmatch均增加，故五项全部FAIL。主对话用这些整数独立复算，与readout及CPU的原五gate完全一致。known nonmatch是相对原参考规则的0标签，不等于生物学真负；unknown仍单列且不改标。

自然候选总数16321→18991，unknown1375→1706。扩大实际池没有改善原门槛：正确链减少18、known nonmatch增加2357、正确长链减少6。此结果否定当前固定配方的TRAIN支持，不证明所有search-aware训练无效，也不能单凭这些计数确定归一化、margin或cutoff资格是唯一原因。

## 工程核验和来源

原fit13398266的root仍为TIMEOUT/0:0；其完整4608行、18个checkpoint至004608、有限649457参数模型及原AdamW、draw／importance weight／teacher监督由reconcile13399750核验。没有恢复训练或额外TRAIN更新，不能把原root改称COMPLETED。

恢复readout13399752与CPU13399755的root／batch／extern均COMPLETED/0:0，分别用时1:00:30／27:35，CPU于07:38:06 UTC结束。STATUS／JOBID相符，free_pool及windows各1491行。CPU重算实际新旧自然池标签、首次丢失cutoff资格与保存margin、全部原曝光和五gate；末模型结构／tensor有限，没有新forward、搜索或优化器更新。

训练记录中的已分配正query disposition合计9084：2506次known nonmatch margin、5530次already returned、1048次primary或unresolved零额外项。这是重复训练曝光中的query记录，不是9084个独立基因，也不与18991个末模型候选/window对混用。资格核验通过不等于目标有效。

权威结果：[readout summary](../../outputs/M28-A30-SEARCH-RETENTION-TRAIN-READOUT-RECOVERY-R1/summary.json)、[CPU summary](../../outputs/M28-A30-SEARCH-RETENTION-RESULT-CHECK-RECOVERY-R1/summary.json)。恢复保留了原模型／seed／曝光／搜索预算和停止规则，没有DEV、sealed或Setaria读取。

## 后续研究边界

[A32固定组合](a32_search_local_caller_protocol.md)因源五gate失败不启动；实现和合成检查保留，不产生真实GPU／CPU结果。独立A31两个原生DEV流程继续，不因A30失败取消或换模型。完整结果闭合后，使用同尺完整链评价和[五状态诊断](a31_dev_score_attrition_protocol.md)区分生成、owner、评分零点和正gain选择损失，再决定下一独立机制；不把本次探索FAIL改写为竞争性或泛化结论。

## 保存监督位置的只读诊断范围

原三epoch额外margin均值0.828587／0.539160／0.465836下降，但原候选门槛仍全部FAIL。进一步诊断只统计已闭合8.9MB训练日志中的首丢失stage、目标节点已闭合exon数、完整参考exon数，按epoch及物种分别报告，并与CPU核验的9084条正query记录／2506个安全竞争对和原epoch均值对应。不读取weights、缓存、FASTA、catalog、R4或DEV，不重训、搜索、重标或改原FAIL。

实际汇总作业13401111的root／batch／extern均COMPLETED/0:0，08:04:24–08:04:27 UTC，用时3秒、batch MaxRSS18484K；private0GPU/2CPU/8G/5分钟，STATUS／JOBID相符。原4608步、9084条正query记录、2506次合格known nonmatch竞争对及原三个epoch均值均对应。只读原保存日志，没有forward、新搜索、R4重标、优化器更新、DEV或Setaria读取，也没有改gate／停止条件。该短作业由主对话消费，不加入A31固定监控。权威输出：[saved loss diagnostic summary](../../outputs/M28-A30-SAVED-SEARCH-LOSS-DIAGNOSTIC/summary.json)。

## 保存监督覆盖了后续搜索阶段，但未改善原门槛

下表只计已核验可用的known nonmatch竞争对；primary／unresolved的零额外项和already returned不混入。每条记录是一次训练曝光中的query观察，不是独立基因。

| 首丢失阶段 | 拟南芥 | 水稻 | 合计 |
| --- | ---: | ---: | ---: |
| local_start_topk | 558 | 221 | 779 |
| global_entry_pool_pruned | 652 | 155 | 807 |
| retained_beam_score_bound | 123 | 47 | 170 |
| continuation_pool_pruned | 178 | 83 | 261 |
| local_next_start_topk | 208 | 81 | 289 |
| completion_pool_pruned | 148 | 52 | 200 |

初始START占779/2506（31.09%）；另1727条记录的目标节点已有非零闭合exon历史，其中667条已有≥4个闭合exon，最大拟南芥21／水稻20。该历史长度是保存目标节点的结构，不是首丢失的搜索迭代编号、失败exon编号或已完成有效gene数；全局entry／bound等阶段也可能涉及不同长度的历史。

| epoch | 物种 | 合格竞争对 | already returned | 目标历史≥4 exon |
| --- | --- | ---: | ---: | ---: |
| 1 | 拟南芥 | 869 | 988 | 210 |
| 1 | 水稻 | 292 | 444 | 63 |
| 2 | 拟南芥 | 524 | 1413 | 176 |
| 2 | 水稻 | 192 | 557 | 51 |
| 3 | 拟南芥 | 474 | 1539 | 137 |
| 3 | 水稻 | 155 | 589 | 30 |

每个epoch的正query观察均为3028条；合格竞争对1161→716→629，already returned1432→1970→2128，primary／unresolved零额外项435→342→271。原margin均值按全部optimizer步计算，包含没有有效额外项的步；其下降不能直接解释为每个活跃竞争对变好，更不等于末模型自由候选支持改善。

据此可以排除“额外监督只实现或触发了最初START，后续阶段完全未参与”的解释。不能由覆盖分布推断监督充分、确定归一化等某一个唯一缺陷，或把降低训练损失作为救回原五gate的理由。当前固定A30仍0/5 PASS，A32仍未执行；下一机制应由闭合A31同尺结果及五状态诊断支持，不对本次参数补试或改写失败。
