# M28 A5 query边界分支解码工程验证

2026-10-03。本轮准备独立的边界分支解码原型，检验真实接口是否能保留贪心路径之外的合法完整CDS候选，并测量固定TRAIN窗口上的开销。原型不是新训练、完整DEV实验或新的主要注释模型；不改A3负结果、运行中的全TRAIN Q/C0及原生基线。

## 依据与固定搜索规则

A3支持范围分解确认，未进入free的3285条主参考中，3191条通过窗口与语法支持检查；这不证明贪心搜索是唯一原因。query_decode.py每个query只选择一个起点和一个后续边界，选择后不回溯。新query_beam.py在同一递归状态、指针和EOS接口上保留有限分支，不复用B1的局部端点图或已失败的跨窗max-gain聚合。

规则在运行前固定：活跃边界历史宽度4、每个query最多返回4条完成链、最多128外显子。每步先按累计分保留4个起点/acceptor历史，再对其所有ORF合法结束事件计算边界与prospective EOS分；donor继续历史保留4条，完成STOP链单独保留4条。没有按DEV表现选择宽度、长度归一或阈值。采用确定性坐标/整链tie规则，不向截断的donor路径补STOP。

累计分是所有几何masked pointer log-softmax及对应continue/stop log-sigmoid之和。它不是旧teacher中按每条链取均值的损失，也不是DNA合法集合上重归一化的概率，更不是基因正确概率；不同链长的分数可能偏向短链。保留历史的未来增量均不大于0，因此已有4条完成链且最佳活跃前缀不超过最差完成链时可以结束；这仅是当前已保留beam的上界，不恢复已剪枝历史，不提供全局最优保证。

DNA规则沿用连续ATG、GT/GC–AG、最短5bp内含子、真实跨exon codon carry、无内部STOP/歧义编码碱基及连续末端STOP。返回候选携带query_index与路径分；暂不把同query的多条链复制原existence gain进入最终selector。完整链质量判别与最终注释协议需要后续明确实验，不能用候选数量当模型成功。

## 合成检查与TRAIN烟测

新增六项合成检查针对实际失败：贪心最高起点无合法结束时保留第二起点、GT/GC及未知内含子跨exon carry、外显子上限/编码N/短内含子不伪造完成、累计分与独立单exon公式一致、padding/softmask/重复执行与query身份、跨剪接产生STOP被拒绝。原模型和贪心入口不改，不重复整套未变化检查。

烟测固定复用已完成full-TRAIN readiness的八条哨兵，每物种四类，选择来自已有TRAIN元数据；不是随机准确率样本。只读取既有全TRAIN计划中的这八条cache及oriented DNA，使用相同seed0未训练CollaborativeGeneQueryModel，两种解码共享一次实际头前向。head为eval且requires_grad false、optimizer更新0、无checkpoint读取/保存、无DEV前向或参考目标解码、无test/Setaria、无新GLM前向。旧cache标签字段不参与监督。

记录贪心与新beam的完成/唯一候选、交集与各自独有、每窗运行时间、显存，以及全部输出链独立拼接CDS/剪接检查；不记录参考命中或F1。缓存工程烟测不冒充新解码器已完成原始FASTA端到端或训练末权重性能。

## 执行边界和下一决定

入口smoke_query_beam.py，输出outputs/M28-A5-QUERY-BEAM-ENGINEERING；sbatch/M28-A5-QUERY-BEAM-ENGINEERING.sbatch使用generanno、1 RTX3090、2CPU/16GiB、10分钟。申请时窗是小规模工程预算，不是可靠全DEV运行预测。仅允许按实际错误作不改变搜索规则/数据/模型的最小工程修复，原始失败产物保留；不把超时或较差候选数量作为扫描宽度理由。

通过只证明候选接口、DNA合法性与八窗工程可运行。原全TRAIN末权重及统一评价完成后，再结合实际模型覆盖、背景和原生对标决定下一次固定结构实验；不在本烟测中调参、选择训练checkpoint或开放独立测试，不提交Git或对外发布。结果与作业证据追加在本文。

## 首次检查与有证据的测试修正

首个作业13358593 FAILED/1:0、6秒，STATUS FAILED exit=1；六项检查中五项通过，一项断言失败，真实TRAIN烟测尚未执行。失败检查将外显子cap下所有完成输出都当作非法，但序列ATGCGTTAGTTTAA的0–9拼接CDS为ATG CGT TAG，确实是合法直接终止链；原贪心偏好的donor4路径与它不同。

CPU-only srun诊断13358598逐一输出四个实际案例，确认只有上述合法0–9替代链被错误排除；含coding N、4bp内含子和无ATG案例均返回空完成池。修正检查允许真实0–9 STOP，仍明确拒绝伪造两exon延续；未修改搜索代码、权重、分数、beam4、completion4、max_exons128或样例数据。

第一次FAILED目录/日志保留。一次有界R2工程重跑使用新目录outputs/M28-A5-QUERY-BEAM-ENGINEERING-R2及新M28A5BEAMR2日志，原sbatch保留；smoke入口仅增加输出路径参数，所有科学/工程输入和规则不变。六项检查必须通过后才执行原八窗烟测。

## 工程终态与推进判断

R2作业13358604已COMPLETED/0:0，32秒，STATUS COMPLETED、JOBID一致；六项合成检查通过（0.070秒）。private-teodoro-gpu/gpu034、1 RTX3090、2CPU/16GiB、generanno，实际命令为 `sbatch --parsable sbatch/M28-A5-QUERY-BEAM-ENGINEERING-R2.sbatch`。前置定位13358598也已COMPLETED/0:0、4秒、0GPU。首次13358593失败产物未覆盖。

八条既有TRAIN哨兵全部完成，每物种四类和原窗口身份一致。贪心每窗16条完成记录、去重后1条；beam每窗64条完成记录，去重后七窗4条、一窗5条。逐窗唯一数加总为贪心8、beam33，其中25为对应窗口中的beam独有；这个加总不是染色体全局去重参考命中数。每窗原贪心唯一链都在本次beam池中，但这只是八窗未训练观察，不保证训练后贪心池始终为beam子集。

所有返回链均通过独立CDS拼接和splice检查。raw records进一步核对：beam512条完成记录和greedy128条均为单exon；合成多exon/codon-carry检查通过，不代表本八窗实际已生成复杂正确基因。未训练Q16在本样本高度重复，记录数增加也不能当作生物学多样性或正确完整链召回。没有参考命中、F1、背景指标或科学gate的新结果。

八窗工程循环8.3465秒，GPU峰值allocated150166528 bytes、reserved192937984 bytes。beam逐窗解码0.1111–0.4583秒；首窗greedy1.0453秒，后七窗0.0395–0.0893秒。执行顺序和首调用初始化可能影响耗时，未做稳态/重复基准，不据总平均宣称加速，也不外推训练后复杂链及全DEV时长。

原始证据在outputs/M28-A5-QUERY-BEAM-ENGINEERING-R2/summary.json（8250 bytes）与records.jsonl（155769 bytes、8行）。主对话复核窗口/物种类别、greedy/beam交集与各自独有守恒、完成预算、记录计数和raw exon直方图，与summary一致；来源是缓存与未训练头，不是原始FASTA端到端或独立注释评价。三份初始Python AST、原/R2 bash语法及新增/修正文件远程精确读回通过；合成验证失败原因和修法保留在上节。

结论：原型具备进入后续固定机制实验的工程条件，尚无证据证明改善A3遗漏的3191条支持链或实际注释质量。下一实验必须在读取新性能前明确权重、输入曝光、候选路径分与最终链质量分的区别，并保留完整参考、逐物种F1和背景条件；不能因本次候选增多直接替换运行中的Q/C0或宣布成功。beam搜索本身不自动构成创新，贡献判断仍需完整模型效益及既有[相关工作边界](../../docs/M28_learned_search_related_work_20261001.md)。当前全TRAIN与原生基线仍按各自冻结协议继续，论文目标未完成。

验证结论仅适用于“工程可运行且记录一致”，不适用于性能排名或论文主张；Markdown结果已保存并读回，未以UI预览作为数值验证证据。

## 固定末权重TRAIN机制诊断

未训练头烟测不回答候选是否正确。本轮另做八窗TRAIN机制诊断，唯一权重固定为已完成A3的step004608，不查看/选择全TRAIN中间权重。两种解码共享同一次A3头前向；greedy保持原实现，beam保持宽度4、每query最多4完成链、max_exons128及原累计分，不调整长度/阈值/预算。

哨兵与上文相同，来自既有full-TRAIN readiness每物种四类；仍是TRAIN元数据选择的压力样例，不是随机总体、独立测试或全DEV。先对八个窗口仅用既有TRAIN cache/oriented DNA生成并保存全部free候选，再读取R4的TRAIN-only窗口与对应catalog做观察。正/保护/未知政策不变，计数分母为这八窗的原eligible positive参考窗配对，不替换7728条DEV主参考；同一参考若跨窗出现按配对计数。报告各自命中、共同命中、新恢复、丢失、均未生成，按物种和exon类分层，并记录原A3 draw曝光与candidate_targets的参考相对1/0/-1类别。

本诊断不做新的optimizer更新、checkpoint写入、DEV前向、gain选择、GFF/背景/F1或test/Setaria。它只能决定当前beam路径是否具备继续准备完整链评分的机制依据；单个压力样本成功不等于性能获益，候选增加却不恢复真实链时不能据此宣称搜索改进。原A3 FAIL、全TRAIN协议及原生基线依赖不改变。

新增analyze_query_beam_train.py、三项分母/迁移/候选类别合成检查和M28-A5-TRAINED-TRAIN-DIAGNOSTIC.sbatch，输出outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC；generanno/private-teodoro-gpu、1 RTX3090、2CPU/16GiB、10分钟。一次固定诊断，不按结果扫描宽度、混合来源或重开A3救援。完整训练后的注释性能及正式后续结构实验仍另按事前协议判断。

## 已训练八窗诊断终态与下一假设

作业13358638已COMPLETED/0:0，root/batch/extern均成功，56秒，STATUS COMPLETED且JOBID一致；三项合成检查通过（0.002秒）。实际命令为 `sbatch --parsable sbatch/M28-A5-TRAINED-TRAIN-DIAGNOSTIC.sbatch`，private-teodoro-gpu/gpu034、1 RTX3090、2CPU/16GiB、generanno。严格读取A3末权重step004608，新增optimizer更新0、可训练参数0、checkpoint写入0；八窗无参考生成并保存完毕后才读取TRAIN标签。没有新GLM/DEV前向、test/Setaria、GFF或性能gate。

计数单位为eligible positive参考窗配对，不是染色体去重基因或DEV分母。27个配对中，共同命中8、仅beam命中7、仅greedy命中0、两者均未命中12，故greedy为8、beam为15。逐物种和exon类如下；这些压力样例不支持总体准确率外推。

| 分层 | 参考窗配对 | greedy命中 | beam命中 | 新恢复 | 丢失 | 均未命中 |
|---|---:|---:|---:|---:|---:|---:|
| Arabidopsis | 17 | 4 | 7 | 3 | 0 | 10 |
| Rice | 10 | 4 | 8 | 4 | 0 | 2 |
| 单exon | 9 | 4 | 8 | 4 | 0 | 1 |
| 2–4 exon | 9 | 2 | 4 | 2 | 0 | 5 |
| ≥5 exon | 9 | 2 | 3 | 1 | 0 | 6 |

逐窗去重候选数加总greedy29、beam94，交集19、beam独有75、greedy独有10，联合池104；这不是全染色体去重。原未训练烟测的greedy子集观察在训练后已不成立：beam丢失10条greedy候选，虽然本样本未丢失原正确参考配对，不能保证其他窗口无正确链损失。后续若保留旧候选可采用同一冻结生成器的greedy与beam联合池，不能用这次局部零丢失假定全局单调改善。

candidate_targets参考相对类别为：greedy精确匹配8、已知非匹配21、保护/未知0；beam精确匹配15、已知非匹配78、保护/未知1。类别0是相对于当前主参考的可计数非匹配，不等于生物学真负或已证明假阳性；类别-1不进入负样本。候选增多同时带来更多非匹配，不能以15/94报告基因precision，更不能取代正式选择、背景及F1评价。≥5 exon仍有6/9配对未生成，复杂链遗漏未解决。

27个正配对所在的具体窗口都未进入原A3 draw；唯一被原draw抽到的哨兵是Arabidopsis零正无保护窗。窗口未抽到不等于基因未见，基因可能在其他训练窗口出现；同属TRAIN且经元数据选择，不是同源隔离泛化或独立测试证据。

八窗预测循环21.3625秒、总22.4967秒，GPU峰值allocated148045312 bytes。逐窗beam0.3398–3.9490秒；首窗greedy受初始化影响，未进行稳态测速，不能按小样本外推全DEV时长。原始产物为outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC/summary.json（14098 bytes）、predictions.jsonl（366680 bytes、8行）、references.jsonl（11085 bytes、27行）。主对话核对27条raw reference状态、逐窗/物种/exon/exposure分母与边际守恒、候选类别/交集及exon直方图计数，均与summary一致。

判断：同一已训练生成器上的分支搜索确实新增7个完整参考链配对，具备准备完整链质量评分的机制依据；不是DEV性能提升或新模型成功。下一假设是冻结生成器，对同权重greedy+beam联合候选学习整链质量，明确区分free与仅TRAIN的参考注入，保留保护标签和原始raw-FASTA/完整参考/逐物种背景gate。先验证梯度、输入身份和冻结边界，再决定正式训练；不把beam、已有链GRU本身包装为创新，不重开A3阈值/权重救援。全TRAIN Q/C0及原生基线保持原协议继续，论文目标仍未达成。
