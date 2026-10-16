# A31 局部路径 caller TRAIN 结果

2026-10-06。A31 A3和A31 A23均通过首次固定五项TRAIN准备条件：两物种known-label窗口配对proxy F1严格增加、proxy precision不下降，合并至少5exon选出正确链增加。读出和CPU整链闭合，两个variant都支持另行固定的完整raw-FASTA DEV测量；这些数值不是完整基因F1，也不能与原生基线的gene F1直接比较。

## 两物种选择结果

原1491唯一TRAIN窗及自然池不变。拟南芥2172、水稻744个eligible窗口参考配对，总2916；其中至少5exon为854和232，总1086。P=正确/(正确+known_nonmatch)，R=正确/eligible参考，F1=2正确/(eligible参考+正确+known_nonmatch)。未知候选照常参加选择，但不当作已知假阳性计入P／F1，必须单独报告。以下表格用于逐计数核对门槛，而非聚合排名。

| 物种 | 读出 | 选出正确 | known_nonmatch | 未知 | 至少5exon正确 | proxy P | proxy R | proxy F1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 拟南芥 | A3 greedy主基准 | 1135 | 490 | 147 | 299 | 69.85% | 52.26% | 59.78% |
| 拟南芥 | A31 A3 | 1531 | 347 | 315 | 548 | 81.52% | 70.49% | 75.60% |
| 拟南芥 | A31 A23 | 1563 | 310 | 338 | 580 | 83.45% | 71.96% | 77.28% |
| 拟南芥 | A6辅助 | 809 | 98 | 73 | 46 | 89.20% | 37.25% | 52.55% |
| 水稻 | A3 greedy主基准 | 425 | 265 | 68 | 111 | 61.59% | 57.12% | 59.27% |
| 水稻 | A31 A3 | 478 | 251 | 132 | 152 | 65.57% | 64.25% | 64.90% |
| 水稻 | A31 A23 | 492 | 209 | 128 | 156 | 70.19% | 66.13% | 68.10% |
| 水稻 | A6辅助 | 262 | 60 | 23 | 8 | 81.37% | 35.22% | 49.16% |

A31 A3相对主基准正确配对分别增加396／53，known_nonmatch减少143／14；A31 A23增加428／67，known_nonmatch减少180／56。合并长链正确从410到700／736，分别增加290／326。未知选择也从147／68增加到315／132或338／128，不能忽略或解释为全部正确。水稻仍是较弱物种，known-label proxy F1为64.90%／68.10%，并非由拟南芥较高值覆盖其不足。

## 自然池与比较边界

A3自然联合池为拟南芥8909／水稻7412链，共16321，其中正确1751／629、known_nonmatch6305／6261、未知853／522。A23为10459／8491，共18950，正确1783／632、known_nonmatch7796／7231、未知880／628。A3 greedy子池为2370／1818，正确1181／466、known_nonmatch1014／1245、未知175／107。全部未生成的eligible参考仍留在recall分母，不注入缺失正确链。

较大自然池的增加与局部分数变化同时发生，不能将上述主对照差全部归因于评分。随后完成的[同模型三格诊断](a31_score_pool_train_diagnostic_result.md)显示：固定greedy池下，局部评分减少非匹配但正确不增加；固定局部分数扩池才增加正确，同时水稻precision下降5.86百分点。缺少第四格，不能推定完整交互或因果贡献。A6在同一A3联合池，但有额外head训练曝光，不能把更低F1解释成匹配训练的因果证明。A31 A23与A31 A3具有同caller配方，却不是同一生成器或匹配容量实验。局部发射势、相位和零点选择本身也不是新颖性主张。

第四格也不能仅靠补存beam归属自动得到：原存在logit属于query，而每query可返回多条完成链，复制分数可能重复计入存在支持；对全部生成query取最大值又可能改变原greedy分数。候选级gain或query联合选择需要另行明确，详见同模型诊断的评分语义边界。本轮不提交这个失真的替代对照，也不以累计非正beam logscore和零null的空选择证明局部caller更优。

## 完整性与执行证据

GPU读出13398521与CPU13398524的root／batch／extern均COMPLETED0:0，分别1:00:41／19:28，batch MaxRSS1752016K／841680K；两阶段STATUS／JOBID相符，两个保存JSONL各1491行，stderr没有失败摘要。CPU重算全部实际R4标签、DP选择、2916参考分母、proxy和五项整数交叉乘法，两个variant逐项PASS，而非按表中四舍五入值判定。

原A3／A23末004608结构及有限值、自然池和greedy子集身份均核对，GPU参数未变且无梯度；全部自然分数文件先关闭再打开TRAIN catalog。原greedy只重放核对来源，没有生成新自然联合池；没有GLM编码、optimizer更新、DEV、sealed或Setaria使用。GPU内部评分3561.714秒／全阶段3569.866秒，peak allocated157948416 bytes，仅是已有缓存池的TRAIN读出耗时，不是原生完整caller速度。

原始计数及字面验证见outputs/M28-A31-LOCAL-PATH-CALLER-TRAIN-READOUT/summary.json和outputs/M28-A31-LOCAL-PATH-CALLER-RESULT-CHECK/summary.json，复算命令verify_local_path_caller_train.py --source-job 13398521。下一步按[a31_local_path_caller_dev_protocol.md](a31_local_path_caller_dev_protocol.md)分别执行两variant的完整8690窗raw-FASTA注释及7728参考评价。该DEV反复使用，仍不能替代独立迁移、同源隔离、相称消融或论文验证。
