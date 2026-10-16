# 原 TRAIN 基线与模型的逐参考缺口诊断

2026-10-09。本诊断区分成熟基线能够精确命中、原4608步模型仍遗漏的参考链，以及两者共同遗漏的参考链，为后续原 TRAIN 模型改进确定具体结构问题。只复用原六条 TRAIN 的完成产物，不新增模型推理或训练，不读取 DEV、封存数据或新物种。

## 配对来源和分母

参考几何、原模型命中标记、CDS段数、跨度与原R4 positive-chain标记来自完成的 `M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE`，作业13449803。原模型错误类别来自完成的 `M28-FULL-GRAMMAR-ORIGINAL-TRAIN-ERROR-PROFILE-R1`，作业13450002；两份记录均绑定原评价13443329和末步004608权重。原模型全范围精确命中9966／22688。

原生命中来自此次18组预测及共同尺评价13461307，只在全评价成功后重新解析已保存的原生 GFF。配对不重读DNA、实际参考GFF或权重，不修改原生预测、参考、参数或运行中的代码。native exact-hit、唯一预测链数、记录数和重复数须逐染色体复现完成评价，再核对逐物种和pooled命中。

参考身份键是species／seqid／strand／有序0-based half-open CDS坐标，phase不是坐标exact-key。参考剖析与错误剖析按species／seqid／strand／transcript ID一一连接，核对CDS段数及原exact命中类别；参考几何必须唯一。六染色体固定分母7135／4293／6311／1638／1601／1710，拟南芥17739、水稻4949、合计22688。缺失、重复、来源或分母不一致即停止，不用删参考修复计数。

## 输出和判断

逐工具列出双方命中、仅基线命中、仅原模型命中、均未命中的四格计数；同时列出四方法全部16种命中组合、任一基线与原模型的四格计数、三基线同时命中而原模型遗漏的参考数。逐参考命中标记保存为 `references.jsonl`，汇总写入 `summary.json`。

按单／多CDS段、跨度≤6144／>6144、原draw R4 positive-chain存在／不存在、原模型七类错误分别分层，并报告每物种。空分母recall为null；所有分层的计数合计须覆盖完整参考。参考级集合交集只说明recall和结构缺口，不为这些子集定义prediction precision或F1，不把仅模型命中称为新基因。

原错误类别来自已保存的描述性最佳CDS重叠关联，不是唯一生物实例配对或因果模块定位。R4标记缺失不等于未训练或独立测试；跨度>6144也不自动代表whole-contig grammar不可支持。全部参考仍进入共同尺，不把共同遗漏当作错误标注的证明。原 TRAIN 及比较器依赖的上游暴露不能建立独立泛化。

后续决策依据是：仅基线命中的缺口集中在哪些原模型错误、CDS复杂度和物种；共同遗漏是否占主要部分。诊断不会自动选择DEV候选、改变冻结八／十方法家族、构造新训练标签或修改已运行的六条训练分支。

## 执行与完成

先以0GPU／2CPU／4GiB／5分钟在Slurm运行合成配对与算术测试，针对ID碰撞、strand、phase、重复链、漏参考、错分母、错误类别漂移、全16组合与空分母。通过后，正式CPU诊断依赖该检查和原生评价13461307的afterok，使用0GPU／2CPU／8GiB／20分钟。两者均为private-teodoro-gpu；输入和输出处于独立目录，不覆盖旧结果。

诊断若需外部等待，只扩展原第六条流程的诊断末端，复用原监控对话和唯一心跳，其他五条依赖不变；已交接的原生数组和评价事件不重复回传。若短诊断在原主分析期间已完成，由原主直接核验消费，不制造第二次等待。完成要求Slurm root／batch／extern均成功、STATUS／JOBID／非空summary与逐参考产物齐全、22688完整分母及源评价命中重放成功；监控只核结构收据，不读数值或逐参考正文。原主负责科学解释和已授权后续工作。

本地8项合成测试及两份Slurm脚本语法检查通过，六份新源文件传输后逐字回读一致。Slurm合成检查13461631的root／batch／extern均COMPLETED0:0，15:37:14 UTC，8项通过，fixture-only且无DNA／实际参考／预测／权重读取、native前向、更新或DEV／封存访问。正式诊断命令 `sbatch --parsable --dependency=afterok:13461631:13461307 sbatch/M28-NATIVE-TRAIN-REFERENCE-DIAGNOSTIC.sbatch` 返回13461637，提交后RUNNING于gpu034。提交前已核用户展开队列、分区、逐节点CPU／内存、预约和存储；仅新增两个CPU短作业，没有请求GPU、取消任务或扩展数据范围。

正式诊断已COMPLETED0:0，15:37:58–15:38:34 UTC，root／batch／extern均成功，STATUS／JOBID与两份非空产物齐全，22688完整分母和全部源评价命中重放通过。主在原生结果分析期间直接核验消费，无额外监控等待。结果见[native_train_reference_diagnostic_result.md](native_train_reference_diagnostic_result.md)。
