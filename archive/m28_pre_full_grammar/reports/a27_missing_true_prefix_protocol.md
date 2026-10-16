# A27 全量遗漏链正确前缀诊断

2026-10-06，在本次模型计算前固定。A24在原1491 TRAIN窗的完成池中发现536个已知正链／窗口配对没有自然生成，但完成池形态不能判定局部指针排名还是活跃历史／完成池竞争。本次用冻结A3和原正确前缀探针补齐这个判断，不提出新caller、不重训失败模型，也不把条件化路径变成自由预测。

## 固定样本和验证对照

来源为已闭合A24（13391716）、A12原正链metadata（13360650）与A6保存的原A3自然池。保留原2916正链身份／几何／来源中的全部536遗漏（拟南芥421、水稻115），涉及396个原TRAIN窗；不按分数、排名或链长挑选遗漏样例。

在每个遗漏窗中，若原自然正链含beam4来源，取原metadata顺序第一条作为必要条件验证对照。实际metadata决定324个对照，不是新挑选的研究正例；它们不并入536遗漏分母，也不与没有对照的72窗作性能比较。两类总860配对，全部原Q16，共13760 query探针。遗漏3711外显子、对照1320，共5031外显子；若无语法不支持，最多80496个完整exon事件。超cap或语法不支持仍保留其身份和已进行的部分事件，不删掉来改善排名统计。

全部窗原定向valid/width均24576，只读原六条允许TRAIN染色体的已有特征缓存和定向DNA。原GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345及A3末004608不变；无新GLM前向、DEV/sealed/Setaria、catalog扩展或新数据。

## 探针与解释

直接复用analyze_true_prefix_train.probe_query，不更改其指针、codon carry或EOS语义。按原A5几何支持log-softmax、合法ATG/GT或GC–AG／最短5bp内含子、donor continue与STOP terminate log-sigmoid，强制正确前缀逐步递归。width4、completed4、maxexons128、同分较小坐标优先均不变；不以query存在分数筛选，不用Hungarian挑“最好query”，不做候选搜索。

每参考／query保存完整逐exon边界、合法／几何名次、组成分数、carry、prospective EOS和首次局部超过4的事件。支持语法的遗漏分为“所有query都被局部必要条件排除”和“至少一个query满足全部局部必要条件”；语法不支持另列，不计为排名失败。最有利query的最大必要名次仅为描述统计，不建议新beam宽度或选择query进行训练。

所有query局部排除，说明该正确链与现有局部排序及固定width不相容；不能据此把增加width当成科学改进。存在局部通过query仍不能区分其他活跃前缀、累计分、提前停止和完成池竞争，也不证明自然搜索会返回。本次不持有全量原池的query级path scores，因此不推断全局精确top4分数上界；原A5的少数压力样例分数界不能移植到这536配对。自然beam命中对照必须至少一个query满足局部必要条件，否则先停止解释并定位身份／定向／数值接线问题。

按物种、原exon类和A24完成池类别分别统计；query首次排除stage计数与reference／window对分母严格分开。TRAIN重复窗口配对不是独立gene recall，诊断不建立性能gate、因果归因或论文创新。

## 执行与下一决定

入口scripts/experiments/M28-METHOD-RESTART/analyze_missing_true_prefix_train.py，唯一新目录outputs/M28-A27-MISSING-TRUE-PREFIX-DIAGNOSTIC，query_profiles.jsonl、references.jsonl、summary.json、STATUS/JOBID。复用原四项真实前缀测试，新增三项聚合检查：reference/query分母不混合、unsupported不伪装rank失败、自然beam对照不能全部局部失败。保存后重读全部profile，核对860身份／Q16、参考几何、path组成和首次排除，以及实际reference记录，不只相信内存汇总。

generanno；private优先、1RTX3090至少24GiB、2CPU/16GiB/40分钟。旧A5实际3456 exon事件、432 probes、57.356秒；本次最多80496事件及396 cache窗，估计约25–30分钟，加操作余量40分钟。这是执行估算，不是稳态速度benchmark；一次运行，不扫描不同配置。新增文本预计低于128MiB，无模型更新或checkpoint写入。故障保留原产物；只允许同科学范围的有证据最小工程修复，不延伸为新方法或超预算训练。

若排除主要在首起点／next_start等局部排序，下一独立表征／目标优先处理相应边界；若较多遗漏仍有局部通过query，下一步须验证真实路径竞争／保留机制，而不能由本诊断直接指认EOS或长度偏好。这个判断与A23/A25/A26的实际固定结果共同使用；三条旧链不改参数、阈值或gate。完成后主对话负责解释，长等待复用同一个watcher/heartbeat，不创建或fork监控。

完整raw-FASTA同输入原生ANNEVO/Tiberius/Helixer比较、最差物种、独立适用范围和论文级机制消融仍未完成；本诊断只减少下一模型迭代的不确定性。没有Git提交／推送、作业取消、旧框架恢复或公开发布。

