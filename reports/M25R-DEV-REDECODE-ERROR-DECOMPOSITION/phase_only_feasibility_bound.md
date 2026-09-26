# Epoch-1 phase-only：全开发集可行性上界

2026-09-08。结论：固定epoch-1/row601，仅移除phase拒绝，无法达到冻结exact-CDS-chain F1 >= 0.55。无需新增全开发集GPU前向来判断这一门槛是否可达。

## 已核定的全量证据

原R4 epoch_1/diagnostic.json的prediction_errors给出：unique reference chains R=6450，unique predicted chains P=2098，exact chains T=1389，nonexact predicted chains=709。与reproduction中的F1=0.3249883013570426一致。

outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE/epoch_1/summary.json的category_counts_by_stratum.all给出correct_chosen_chain_phase_rejection=1149。这不是面板外推，也不是旧original_stage_counts.phase_check=1145；后者来自旧的一对一归因，不能拿来替代按实际chosen-chain重建的全量归因。

对应group_phase_retrospective.py先按chosen_cds构建全候选索引，将所有参考链与实际选定链比较；primary_category先排除已emitted精确链，再将正确chosen链且terminal为phase_check的参考归入该类。没有mixed-terminal类别。单独绕过phase不会改变更早的region state、block、CDS runs或motif选择，故新增精确链数量x至多1149。

src/foundation_probe/train_generanno_structural_heads.py:_validation_metrics使用同物种、染色体、链及完整CDS坐标链的集合F1；原输出不会因仅取消phase拒绝而删除。新增nonexact unique链数记为y>=0。

## 上界推导

F1(x,y) = 2*(1389+x)/(6450+2098+x+y)，0<=x<=1149，y>=0。

它随x递增、随y递减。因此极端乐观地允许1149条全部通过剩余ORF/阈值检查，并且一个新增nonexact都没有：

F1_max = 2*2538/(6450+3247) = 0.5234608642 < 0.55。

这是可达性上界，不是实际B运行分数。真实剩余过滤或新增nonexact只能使这个上界更低。即使6条面板nonexact有其他注释链支持，也不能改变冻结complete-primary定义下的集合计数。原709条nonexact仍被保留。

由一个必要门槛不可达即可否决这个固定机制方案；不需要推测interval F1、FPR或泛化，也不需要跑完整B得到一个更低的数。

## 对原GPU计划的处理

已核定两个开发染色体长度23459830和29936421 bp，6144 bp非重叠窗口、双链共17384窗口。R4代码只在内存保存logits，现有epoch目录没有全量score缓存；面板仅158窗口。完整11通道float16 logits约2.19 GiB（未含其他产物）。面板158窗口338.5秒对应全量约10.35小时纯前向的粗略线性估计，不是运行保证；完整A/B可共享一次前向。没有提交该GPU计划。

由于上述严格上界已经否决固定epoch-1 phase-only达到现有门槛，取消把完整GPU重放作为下一项默认实验。若另有研究目的需要完整效应量，它仍可以作为新的探索性实验单独讨论，但不能以“验证能否过门槛”为理由消耗资源。

## 科学边界与下一决策

本结论关闭的是固定epoch-1/row601的phase-only修正作为达标方案，而非所有GENERanno结构注释方法。面板仍证明phase硬门控会阻断正确链；上界进一步证明这一修正不足以解决主要结构召回瓶颈。

全量归因另显示1634条参考的CDS支持碎裂于多个block、1075条同run数但真边界候选不可达。它们早于phase门控，phase-only不会修复。改变block组装、候选生成、监督或训练属于新的机制范围，不在本次批准中。

因此当前可交付的是“phase-only机制修正不足”的证据闭合NO-GO，不是整个研究目标已经完成。下一步若继续，应先提出针对上游碎裂的有限机制修正与明确预算/停止条件，再请求用户作研究范围决策；不自动训练、不改变冻结评价、不释放Setaria、不把局部NO-GO冒充全路线终态。

本报告是原phase_gate_bypass_result.md和nonexact_forensics_result.md中“建议完整GPU配对对照”的后续修正，以上界证据取代该默认下一步。原报告作为过程记录保留，不覆盖冻结结果。
