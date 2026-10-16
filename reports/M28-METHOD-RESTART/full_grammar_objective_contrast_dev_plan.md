# 三种续训目标的冻结 DEV 比较方案

2026-10-09 UTC，本轮先完成原 TRAIN 三个分支的训练、严格 raw 推理和端到端诊断，再比较它们在既有两条 DEV 染色体上的完整 CDS 链表现。准备阶段只做接线和人工检查，不提交或运行 DEV 推理，不读取 DEV DNA、参考、旧预测或新指标；23:15 UTC 的实际执行记录见末节。三种目标的末权重已在训练前固定，结果不用于选择中间权重、改目标权重或修改评价尺。

## 候选和锚点

固定纳入 STRUCTURED_CONTINUE、LOCAL_CHANNEL、HALF_STRUCTURED_HALF_LOCAL 三个分支各自的 step_009216.pt，包含原 4608 步起点加 4608 步续训。所有成功完成分支按同一口径报告，不按 TRAIN 得分预筛“赢家”。原 step_004608.pt 的冻结 DEV 比较作为不变锚点，复用原结果，不重训或重选。失败、超时、未完成或记录不完整的分支报告其实际状态，不记作零准确率，也不从候选清单静默移除。

当前既有 TRAIN 链为三个 fit → CPU 核验 → 六染色体 raw 数组 → CPU 评价，详见[续训合同](full_grammar_objective_contrast_contract.md)。进入 DEV 前三条 TRAIN 链均须完整关闭；这个顺序检查只读完成身份，不读取分数或制定 TRAIN 门槛。若任一链失败，先处理已授权的实现或执行问题，保留失败证据，不先评估剩余“赢家”。DEV 不进入训练、目标系数、学习率、阈值、cap、窗口或 decoder 选择。

## 不变输入和评价尺

沿用[原冻结 DEV 方案](full_grammar_dev_baseline_plan.md)：拟南芥 NC_003074.8 / 23459830 bp、水稻 NC_089041.1 / 29936421 bp；完整原 FASTA、双链整染色体、nearest-center owner、exact exon-budget skyline、cap 128、unit role gain、null 0、FP64 和原 tie policy。不加入参考边界、softmask、RNA 或蛋白证据，推理不读取标签或 feature cache。保留大写原 IUPAC 序列供输入直接比较，caller 在模型前将非 ACGT 转 N，不改变长度坐标。

参考分母不变：M26 CDS-assessable longest-eligible-CDS / tie-transcript-ID，拟南芥 5437、水稻 2291、合计 7728 unique primary CDS chains。预测端仍按每个原生 predicted gene 最长 CDS / transcript-ID tie，complete_only=False；不按真值匹配或 M28 语法删基线预测。参考中长链、重叠、motif、partial 与模型可支持性按既有规则保留，不用训练 eligibility 缩小测试分母。精细 phase、其他 assessable isoform 命中和结构分层另计，不替代 primary exact chain。

## 固定基线和输出

复用 2026-10-03 六份成功原生输出，不重跑或更新 ANNEVO、Tiberius、Helixer。源码 / 权重 / lineage / native operation point 与原方案一致：ANNEVO 42c920f / 2.3.3 / Magnoliopsida；Tiberius 871d497 / 实际 2.0.8 隔离包 / angiosperms；Helixer d17bb49 / 0.3.7 / land_plant_v0.3_a_0080.h5。每次比较仍核对实际运行收据、原生操作点、Slurm 成功终态和候选 / 基线 / 原 FASTA 逐碱基相等后才读取参考。

三个候选各自输出原共同尺的逐物种 TP / predicted / reference、exact-chain P/R/F1、CDS-base union、背景 span / CDS 覆盖、纯背景链 / Mb、重复链、其他 isoform、phase 和描述性结构层；pooled 按计数 micro 汇总，同时列 macro 和 worst species。汇总表同时展示原 final4608、三个 final9216 和三条固定基线；复用基线计数与锚点不一致时先核对输入或评价问题，不以不同来源的数字拼表。每个候选相对锚点及每基线的完整链 F1 差值按百分点报告，不把 pooled 单点优势解释成两物种都优越。

候选路径为 outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-RAW/<arm>/<seqid>，比较路径为 outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-COMPARISON/<arm>。原 final4608、原生基线和历史评价不覆盖。后续预计沿用每 DEV 元素 1 RTX3090、2 CPU、16 GiB、2 小时及 CPU 比较 2 CPU、8 GiB、20 分钟；提交前以实时队列和维护状态确认资源。当前不把这些估计视为已提交作业或实际完成时间。

## 科学解释和执行状态

该 DEV 已被多轮开发使用，模型和基线的上游物种 / annotation 暴露有方法相关差异或未知，详见[pretrained_exposure.md](pretrained_exposure.md)。此比较是同物种开发集读数，不能声称独立泛化、clean zero-shot、全基因 / UTR / 所有 isoform 能力或论文竞争优势已成立。目标家族对照只有一个 seed，且梯度尺度不匹配；它不能单独建立因果机制、新颖性或收敛结论。

准备代码显式传入候选路径、末权重、身份、协议和方法名，复用原评价函数；原默认 final4608 入口保持不变。人工检查只验证新权重身份、不同 arm 错配、TRAIN 未关闭时不触碰 DEV、原 IUPAC 输入保留及旧路径回归，不运行编码器或研究数据。结果解释、真实 checkpoint / 输入及全部产物核验仍待当前 TRAIN 流程完成后进行。封存 / held-out / Setaria 不放行，不公开论文或模型，不修改现有监控对象。

## 2026-10-09 准备验证结果

仅提交人工 CPU 检查 `13451833`：`sbatch sbatch/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-CHECKS.sbatch`，private-teodoro-gpu / 2 CPU / 4 GiB / 5 分钟，无 GPU。实际 08:20:00–08:20:21 UTC；root、batch、extern 均 COMPLETED / 0:0，STATUS=COMPLETED，summary completed=true，4 项测试通过、0 failure / error。其中旧 final4608 路径回归包含原 5 项人工接线检查。

结果路径为 `outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-CHECKS/summary.json`。未读取研究 DNA / 标签、历史 DEV 收据 / 预测或实际 checkpoint；未做原生编码器前向，optimizer_steps_this_run=0，DEV_used / sealed_data_used=false。验证完成的是工程接线，不是真实末权重或 DEV 科学评价。DEV RAW 与 COMPARISON 作业脚本已准备，均未提交。

固定原监控对话和原心跳保持三条 TRAIN 流程的监控对象不变。主对话已核验 08:18:27 UTC 的真实定时触发确实追加到原监控对话、实际 gpt-5.6-luna / low；08:18:44 UTC 的调度检查确认三个 fit RUNNING、九个后续阶段 PENDING / 有效 Dependency。该次不是人工首检，完整终态回传仍待流程结束。下一步是当前 TRAIN 全流程结束后的结果核验与三目标比较，再依本冻结方案接入 DEV，不重复检查已消费的人工检查作业。

## 2026-10-09 三目标 DEV 实际执行

原三条 TRAIN 链已全部成功闭合。STRUCTURED_CONTINUE 的最后一项 raw 推理于 22:52:17 UTC 结束，CPU 评价 13451625 于 22:53:22 UTC 结束；主对话于 23:08:42 UTC 完成三目标结果核验和比较，详见[原 TRAIN 完整结果](full_grammar_objective_contrast_contract.md#原-train-三目标完整比较)。本轮不按该比较预筛候选，三个预声明末权重全部进入冻结 DEV。

23:15:49–23:15:50 UTC 提交以下六个作业。每个 raw 作业均为 task 0–1、最多同时 2 项；对应 CPU 评价的 afterok 依赖覆盖完整 raw 数组，不是单一 task。

| 固定末权重 | DEV raw 数组 | CPU 共同尺评价 |
|---|---|---|
| STRUCTURED_CONTINUE / final9216 | 13465374_0–1 | 13465375，afterok:13465374_* |
| LOCAL_CHANNEL / final9216 | 13465376_0–1 | 13465377，afterok:13465376_* |
| HALF_STRUCTURED_HALF_LOCAL / final9216 | 13465378_0–1 | 13465379，afterok:13465378_* |

task 0 为拟南芥 NC_003074.8 / 23459830 bp，task 1 为水稻 NC_089041.1 / 29936421 bp。六项实际在 23:15:50–23:15:51 UTC 启动：STRUCT 两项位于 gpu034，LOCAL 两项和 HALF task 0 位于 gpu035，HALF task 1 位于 shared-gpu / gpu021。其实际数字 root 分别为 13465380、13465374、13465381、13465376、13465382、13465378；raw 的 JOBID 仍保存数组 ID_task，started.json / summary.job_id 保存数字 root，终态须按 sacct 对齐。

每项沿用 1 RTX3090、2 CPU、16 GiB、2 小时时限，提交时仅允许 private-teodoro-gpu,shared-gpu 两分区排程；CPU 沿用 0 GPU、2 CPU、8 GiB、20 分钟，均 no-requeue。提交前核验九个目标输出路径不存在、存储充足且没有相交维护窗口；不覆盖原 final4608、原生基线或冻结结果。实际启动收据的 arm / task / 物种 / seqid / 长度 / 对应 step_009216.pt 均匹配，all_three_TRAIN_branches_closed=true、reference_used=false。23:17 UTC 的调度核验显示六项 raw 均 RUNNING，三个 CPU 评价均 PENDING / 有效整数组 Dependency。这是执行和启动证据，不是推理完成或 DEV 科学结果。

原 fixed final4608 与 2026-10-03 原生基线只复用，不新训或重跑；评价口径、7728 参考分母、推理参数和末权重选择不变。本轮不访问 sealed / held-out / Setaria / 新面板，不据 DEV 回调训练或冻结家族，不公开或提交 Git。继续复用原固定监控对话与同一心跳；三候选的六项 raw 和三个评价全部闭合后，主对话负责科学核验与完整三目标 DEV 比较。其他七条已提交 TRAIN / 补充诊断流程保留各自独立完成条件，不为本轮 DEV 增设门槛。

## 三目标 DEV 终态与结果

2026-10-11 UTC，主对话直接核验实际终态并处理结果。六项 raw 13465374_0–1、13465376_0–1、13465378_0–1 及三个评价 13465375、13465377、13465379 的实际 root／batch／extern 全部 COMPLETED／0:0。STRUCT 两项 raw 最迟 10 月 9 日 23:53:05 UTC 结束，评价 23:53:06–23:53:52；LOCAL 最迟 23:55:36 结束，评价 23:55:37–23:56:23；HALF 最迟 23:52:31 结束，评价 23:52:31–23:53:28，extern 23:53:29。没有重提、修改或取消作业。

全部 STATUS／JOBID、非空输入／预测 GFF／内外 summary、对应最终权重与六项数字 root 匹配。正式末权重累计 9216、原 6144 编码、epsilon0.1／null0／unit gain／cap128 保持不变；双链整染色体 skyline、FP64、原 tie policy 与大写 IUPAC 输入相符，推理 reference／label／cache=false、DEV=true／sealed=false、零更新。成功完成的原共同评价在访问参考前验证了真实候选／原生／原 FASTA 逐碱基相等、完整网格和原生操作点；之后核验实际 DNA／phase／cap、预测 GFF roundtrip 和同链非重叠。主对话重新按原计数计算 chain／base P／R／F1、背景合并、micro／macro／worst；三个候选与旧锚点的全部固定原生测量、参考层、CDS reference bp 和背景分母一致。

| 固定方法 | AT 链 F1 % | Rice 链 F1 % | pooled micro F1 % | macro F1 % | 背景链 |
|---|---:|---:|---:|---:|---:|
| 原 final4608 | 50.8442 | 30.0960 | 44.6085 | 40.4701 | 64 |
| STRUCTURED_CONTINUE final9216 | 62.9235 | 36.9312 | 54.7508 | 49.9273 | 79 |
| LOCAL_CHANNEL final9216 | 44.4851 | 36.8666 | 42.1812 | 40.6758 | 173 |
| HALF_STRUCTURED_HALF_LOCAL final9216 | 51.8012 | 42.4519 | 48.9426 | 47.1265 | 245 |
| ANNEVO | 76.2651 | 74.4318 | 75.7025 | 75.3485 | 146 |
| Tiberius | 74.1771 | 66.5822 | 71.8047 | 70.3797 | 165 |
| Helixer | 67.4564 | 53.5243 | 62.8512 | 60.4903 | 563 |

参考分母仍 AT5437／Rice2291／pooled7728。三个候选的 TP／unique predicted 为 STRUCT：AT3362／5249、Rice905／2610、pooled4267／7859；LOCAL：AT2430／5488、Rice873／2445、pooled3303／7933；HALF：AT2840／5528、Rice1025／2538、pooled3865／8066。STRUCT pooled precision54.2944%、recall55.2148%；HALF47.9172%／50.0129%，LOCAL41.6362%／42.7407%。本表所有方法的 worst species 均为水稻，worst F1 即表内 Rice F1；不省略较差物种。所有候选 duplicate_chain_records=0；其他 eligible isoform 命中另记，不充作 primary TP。

STRUCT 相对旧锚点 pooled 提升10.1423个百分点；相对 LOCAL 提升12.5695个百分点，尽管 pooled CDS-base F1 较低1.0514个百分点。HALF 相对 STRUCT 的 Rice 链 F1 高5.5206个百分点，但 AT低11.1223、pooled低5.8081个百分点，背景链多166。pooled CDS-base F1 为 STRUCT92.7560%、LOCAL93.8074%、HALF93.8948%。三目标对照仍显示完整链与 base／背景之间的权衡，不能称混合目标全面改善。

三个候选在两个物种均低于三个固定原生基线；其中 STRUCT pooled 仍低于 Helixer8.1004、Tiberius17.0539、ANNEVO20.9517个百分点。此处只是已完成的七方法三目标子比较；完整[八方法家族](full_grammar_native_cds_prior_dev_plan.md#完整八方法-dev-终态与结果)和单独[十方法扩展](full_grammar_sampling_context_dev_plan.md#十方法扩展-dev-终态与结果)现已分别汇总，原三目标测量不变。不据此调参、选权重或增删冻结候选。单 seed、历史复用同物种 DEV、上游 annotation 暴露限制不变；这不是独立泛化、clean zero-shot、新颖性或论文成功证据。

数据来源：三个 `outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-COMPARISON/<arm>/summary.json`，六项对应 `DEV-RAW/<arm>/<seqid>`，原 `outputs/M28-FULL-GRAMMAR-FIXED-DEV-COMPARISON/summary.json`。旧心跳一直保持上一终态后的暂停；本组三目标 DEV 没有已验证的独立首检、定时追加或终态自动回传，不能把主对话直接核验记成监控回传。该组已由主消费，不再恢复其旧作业的轮询；原固定监控留供真正的新作业复用。
