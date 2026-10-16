# M28 A9：TRAIN 机制筛选失败，关闭本头的 DEV 扩展

2026-10-03 14:33 UTC。拟合和诊断正常完成，但冻结继续条件失败。新 query-conditioned prefix critic 在原8个 metadata-selected TRAIN pressure windows上减少了自由整链可达性；不能以训练 loss 下降或候选多样性增加宣布改进。本头不进入完整 DEV，不追加训练、换 seed、扫 beam 或改阈值挽救。其他已提交的 FULL、A6/A7、A10 分支继续原协议。

## 实际执行与结果身份

| 阶段 | Slurm job | root/batch/extern | 实际时长 | 资源 |
| --- | --- | --- | --- | --- |
| 固定拟合 | 13359395 | COMPLETED / 0:0 | 1h40m21s | private，1 RTX3090、2CPU、16GiB |
| 固定 TRAIN probe | 13359396 | COMPLETED / 0:0 | 51s | private，同资源 |

唯一正式 critic 为 `outputs/M28-A9-QUERY-PREFIX-FIT/step_004608.pt`，822,786 bytes；生成器仍为冻结的 A3 step_004608，新头67,745参数。原1,536 draws／1,491 unique TRAIN windows按原顺序重复3遍，4,608更新、三轮各1,536步；1491份自由 trace 先关闭，再读取 TRAIN catalog。未用 DEV 选 epoch，未访问 test／Setaria；probe 自身 optimizer 更新0。协议见[a9_query_prefix_fit_protocol.md](a9_query_prefix_fit_protocol.md)。

三轮 mean weighted loss 为0.1274281293、0.1142786260、0.1114135387。每轮已知正／负事件27,674／6,881,928，unknown 5,192,615；68窗没有known事件、grammar-unavailable teacher events为0。这里保留原known-mean目标，loss下降不证明搜索目标或完整链排序学到正确行为；极稀疏正事件是预先记录的风险，不是本结果已证实的因果解释。

## 原27条 reference/window pairs，完整分母

early 是新critic-only beam4；late 是原A3 greedy+beam4自由池。两臂共享同一个新头的完整 STOP readout；不是早／晚完整值不同造成的配对差异。late多一个greedy来源，搜索过程及预算并不完全相等。

| 分组 | early 自由 exact pairs | late 自由 exact pairs |
| --- | ---: | ---: |
| 拟南芥 | 0 | 7 |
| 水稻 | 4 | 8 |
| 合计，分母27 | 4 | 15 |
| 1 CDS segment | 4 | 8 |
| 2–4 CDS segments | 0 | 4 |
| ≥5 CDS segments | 0 | 3 |

完整成员转移：both=4、early_only=0、late_only=11、neither=12。新增free exact pairs=0，丢失11；新增≥5-segment pair=0。原27条身份、链几何、exon分组、draw曝光及顺序与A5记录逐项一致；原8窗无重复、均为TRAIN，预测记录明确未读参考。主对话从实际27条 `references.jsonl` 重算合计、物种和转移计数，从实际8条 `predictions.jsonl` 重算共同(query,chain)值的一致性，全部与summary一致。执行代码在预测文件写出前，对两臂union中的每条返回链调用原DNA合法性检查；本次没有重新编码或进行模型forward。

冻结最低机制条件是 early>late 且至少新增1条≥5-exon reference/window pair，并保持共同值、完整分母与DNA合法性。本次4不大于15、无新增多exon链；`preparation_gate_for_full_DEV=false`与直接重算一致。工程完成和科学门槛失败分别记录，不把它写成运行失败或未执行。

## 候选负担与零操作点视图

| 原8窗池的求和统计 | early | late |
| --- | ---: | ---: |
| completed记录数 | 512 | 616 |
| 去重(query,geometry) | 512 | 537 |
| 每窗去重geometry | 281 | 104 |
| reference-relative标签1 / 0 / -1 | 4 / 277 / 0 | 15 / 88 / 1 |
| complete value>0 的geometry | 0 | 0 |
| value>0 的reference/window hits | 0 | 0 |

early虽然有更少completed记录，却产生更多不同geometry，并减少正确整链。geometry统计是在各窗内去重后求和，不是最终全基因组gene数；标签0是参考不一致，不是已证实生物学假基因。两臂的全部完整值均未越过固定零操作点，不能把4或15条free可达链当作最终注释TP。该附加读出不进入继续gate，不据此调阈值；没有最终GFF、F1或背景负担结论。

## 判断与后续

本结果否定的是当前冻结A3生成器、known-event BCE critic、固定搜索和最后checkpoint在这8个TRAIN压力窗上的最低机制条件；不是对所有学习搜索、cost-to-go目标或新模型路线的普遍否定。当前目标是参考兼容性分类，而非已验证的cost-to-go／DAgger／LOLS／SaIL实现，见[a9_search_objective_audit.md](a9_search_objective_audit.md)。不因已经准备了raw paired helper就绕过FAIL提交DEV。

原始证据保留于 `outputs/M28-A9-QUERY-PREFIX-{FIT,TRAIN-PROBE}`：4608训练行、末头、summary、186,833-byte predictions.jsonl及27条references.jsonl；probe stderr为空。主对话核对终态、保存配置和数值，不冒充独立评阅。没有修改checkpoint、预测、历史负结果或冻结gate，没有Git提交／对外发布。

原固定监控对话在14:27 UTC的自然定时运行成功向主对话交接；用户无须再次批准接续。主对话消费A9结果后，从同一heartbeat移除此两根ID，保留FULL五根、A6三根+A7、A10四根共13根，模型与15分钟间隔不变，完整prompt／身份回读一致。watcher曾声称暂停旧cohort，但当时实际配置仍ACTIVE；本记录不把该声称当作已验证暂停。新13根cohort的下一次定时运行尚待观察。
