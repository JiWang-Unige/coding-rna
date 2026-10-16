# A8 完成池整链质量固定 TRAIN 原型

2026-10-03，运行前固定。用户已授权继续迭代，原型依据[A5正确前缀诊断](a5_true_prefix_diagnostic.md)：12个遗漏配对中8个局部前四不通过，4个有局部可通过query；其中5个在所有query已有4条更高累计分合法链，原得分的精确top4也不能返回它们。这里只检验完成池保留机制，不把beam/GRU本身称为创新，不改FULL、A6、A7、历史FAIL或封存数据。

## 固定方法与配对对照

同一A3 step004608生成器、同一已完成A6评分头step004608，均eval且全部冻结。TRAIN八个原压力哨兵、24576bp cache/oriented DNA、Q16、活跃beam4、每query完成4、max_exons128不变，无新拟合、GLM前向、DEV、test/Setaria、权重选择或阈值扫描。

独立quality_complete.py不修改原query_beam或运行入口。活跃起点/acceptor及donor历史仍按原累计pointer/EOS分剪枝，DNA/ORF/codon carry规则不变；一直推进至无活跃历史或128 exon上限，不使用旧累计分完成池的提前停止界。非加性整链质量没有该界，不能将旧界直接移植。所有实际探索的合法STOP链在质量头看到之前不按路径分截掉；同几何只评分一次。

一次共同自然搜索，事后给出两个固定完成池读出：continued_path按原path_log_score保留4/query，quality按固定A6 additive+GRU-window-null质量分保留4/query。两者共享完全相同的被探索路径，不由质量更改活跃历史。故本原型是“未提前截断完成池的整链质量读出”，不是质量学习的prefix搜索，也不是全局最优。分数tie按整链坐标，不复制query existence，不过滤非正质量分。

同时保留原greedy/beam基线，比较greedy+continued_path、greedy+quality、原基线+quality联合池；另记录greedy+全部实际探索完成链作为候选阶段诊断，而不是最终注释。expanded_quality明确保留原基线，不能把人为保证的零候选丢失当作科学改善。移除有效路径分停止界应保持原path top4；逐query几何/分数与原A5记录对照，差异须先定位，不把额外计算冒充质量效益。

先将8窗全部自然候选、质量分与运行开销写入predictions.jsonl，再读取原TRAIN-only R4标签和27个参考/窗配对。无参考条件化推进或注入。保持身份/分母，按物种/exon类记录各池参考命中、新增/丢失、候选计数及1/0/-1；0只表示当前主参考非匹配、-1保护，均不是生物学真负或基因precision。不做GFF、背景、F1或主要gate。

## 检查、执行与下一判断

六项合成检查针对非单调质量不能沿用旧分数上界、同路径配对读出、几何去重评分/非正分保留、cap/真实跨exon STOP、padding/softmask/existence、空池及非有限/错配分数。先通过检查再执行8窗。新增analyze_quality_complete_train.py与M28-A8-QUALITY-COMPLETE-TRAIN.sbatch；输出outputs/M28-A8-QUALITY-COMPLETE-TRAIN，generanno/private、1 RTX3090/2CPU/16GiB/30分钟。原21秒八窗beam与57秒432正确前缀探针仅用于操作余量；延长遍历和整链评分预计更慢，30分钟不是稳态基准或全DEV预测。失败保留产物，只允许不改变科学配置的有证据最小工程修复，不扫描宽度或延长科学训练。

只有合法输出、原path控制一致，并且expanded_quality新增至少一个完整参考/窗配对、quality相对continued_path也新增、原基线未丢失，才具备准备完整DEV固定协议的机制依据；不是准确率PASS或论文成功。若全部探索池仍不生成链，优先定位活跃历史/局部评分；若生成但质量不保留，定位整链评分，不在同小样本反复调参救本轮。完整模型效益仍须A6/A7、逐物种完整参考与背景/原生对标及独立测试证据。不提交Git或公开。

## 终态：工程通过，机制继续条件未通过

实际命令 `sbatch --parsable sbatch/M28-A8-QUALITY-COMPLETE-TRAIN.sbatch`，作业13359250，gpu035/private、1 RTX3090/2CPU/16GiB，root/batch/extern均COMPLETED/0:0，1分13秒，batch MaxRSS1599656K。六项检查通过0.482秒；STATUS/JOBID一致。实际数据阶段54.787秒、GPU峰值allocated148842496 bytes。冻结A3生成器及A6末评分头均通过已有严格身份加载，无新optimizer、checkpoint、GLM/DEV前向或封存数据。

8个窗口和27个参考/窗口身份全部保留。主对话从原A5和本轮raw记录独立复算逐query path控制几何/分数（最大差0）、所有池的候选集合/命中/分层/交集及质量分映射；全部实际探索的完成链独立拼接CDS并检查坐标、GT或GC–AG/5bp内含子、真实跨exon阅读框、无编码N/内部STOP和末端STOP，均一致。predictions.jsonl为2398900 bytes/8行，references.jsonl为11262 bytes/27行，summary.json为13494 bytes；完整自然记录关闭后才读取TRAIN标签。

| 固定读出 | 参考/窗命中（27） | 窗内去重候选加总 | 当前主参考相对1/0/-1 |
|---|---:|---:|---:|
| 原greedy+beam基线 | 15 | 104 | 15/88/1 |
| greedy+continued_path | 15 | 104 | 15/88/1 |
| greedy+quality | 11 | 243 | 11/231/1 |
| 基线+quality联合池 | 15 | 307 | 15/291/1 |
| greedy+全部实际探索完成链 | 16 | 4522 | 16/4505/1 |

这些都是窗内候选/参考配对，不是染色体去重基因数、precision或F1。类别0不等于生物学假基因。continued_path逐query与原beam完全一致，质量读出相对该控制新增0、丢失4个配对（单exon2、2–4 exon1、≥5 exon1）。expanded_quality因为保留基线而不丢失，但新增0；mechanism_preparation_supported=false。不把“输出更多合法链”当作模型改善。

拟南芥17配对：基线7、quality5、联合7、全部探索8；水稻10配对：基线8、quality6、联合8、全部探索8。唯一新增进入全部探索集合的是rna-NM_101226.5，18 exon；它未进入quality前四，当前窗A6质量gain为-1.7499676943，即使保留，在该窗冻结零操作点下也不是可选择的正gain。不能将新增进入未截断集合单独归因于去掉停止界：原A5没有保存所有中间完成链，本原型同时不预截断完成池。

quality池候选exon数均值30.05、中位33、最大72；原基线为3.75、2、23。复杂候选增多却少命中4个主参考，是当前输入下的分布/排序观察，不证明这些长链都是生物学错误，也不单独证明additive项或长度是根因，不能据此直接加长度惩罚或降低阈值。尚有11个主参考配对连全部实际探索集合都未生成，完成池重排本身无法恢复它们。

决定：本次负结果保留，不提交A8完整DEV推理，不扫描beam、质量阈值、长度归一或换评分权重救本轮；不改A6/A7或全TRAIN任务。后续优先定位/改进生成器的活跃历史与局部query评分，并结合正在进行的全TRAIN曝光及A6/A7完整结果，区分生成缺失和评分失败。A8仅证明原型工程正确、当前机制假设未获支持，不是论文方法成功、独立泛化或对A6完整DEV的否定。
