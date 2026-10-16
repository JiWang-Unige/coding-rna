# A19 全窗 query 保护与局部监督覆盖诊断

2026-10-05，在新计算前固定。问题是原A3的全窗 unknown query 政策是否经常由局部保护跨度触发，以及两物种实际监督曝光如何不同。这是已有TRAIN记录和标签元数据的描述性诊断，不改标签、权重、匹配或运行中的A15／A17，不以此自动批准新的首训。

## 已有事实与方法先例

`unmatched_query_policy`遇到unknown chain span、完整protected非positive链或有效DNA中的非ACGT时，未匹配query均标为-1；原query的后续pointer支持仍为全窗。原A3已保存三遍监督计数：每遍拟南芥query 1／0／-1为2242／2831／7215，水稻786／8069／3433。这些是query／draw计数，不是独立基因数或梯度强度；实际importance权重和known-query mean分母仍影响更新。

局部位置query和稀疏标签处理本身已有先例。[DAB-DETR原论文第4.2–4.4节](https://arxiv.org/html/2201.12329)用显式box坐标形成位置query、逐层更新anchor，并以宽高调制位置attention；这是动态soft ROI pooling，不证明query输出有硬地域约束。[SparseDet原论文第3.4–3.5节](https://arxiv.org/html/2201.04620)用提案置信度和重叠关系区分已标注、潜在未标注与背景区域，对未标注区域施加一致性监督。它们是原始视觉任务研究，不提供gene模型准确率，也不能据此把本项目当前unknown链改成负例。

对本项目的比较推断：A17只限制初始pooling，Hungarian分配、首起点、后续splice和跨度仍为全窗；不能直接按初始bin授予负标签。即便一个起点不在保护跨度中，继续链仍可能跨入保护范围。未来若检验局部监督，必须先定义可输出链的地域／可判定性语义，保留跨区长链、partial／alternative和unknown，而非仅缩小ignore区域。以上先例不证明本项目失败的原因，也不建立新颖性结论。

## 计算范围

独立CPU脚本`audit_local_supervision.py`复用原TRAIN metadata loader，只取六条允许染色体、1491唯一窗和原1536 draw及权重。核对A3实际4608条记录的epoch／ordinal／window／weight／监督计数，逐窗重建metadata可识别的两种保护原因，并与三个epoch保存汇总一致比较。不开checkpoint或特征tensor、不读DNA、不做模型／GLM前向、搜索、optimizer、DEV或test／Setaria。

保护空间是chain_unknown_intervals与完整protected非positive链span的半开区间并集；合并重叠和相触区间，避免重复累计bp。另按A17固定Q16／有效长度划分初始pooling bins，报告不与并集相交的bin数量。区间合并和半开边界各有一项解析检查。非ACGT精确位置未从元数据获得，因此含保存`noncallable_valid_base`原因的窗保留监督计数，但不进入空间覆盖统计，明确报告其数量；不为本诊断读取约50GiB特征cache。

分别报告1491唯一窗与原每遍1536 draw的计数、逐物种、guarded窗和regional分母、并集bp、并集外bp比例及初始bins计数，保留所有窗口。并集外不是生物学真负背景，也不是已经合法回收的query负标签。空间统计未按importance加权，不宣称其等于训练损失或梯度比例。

独立输出`outputs/M28-A19-LOCAL-SUPERVISION-DIAGNOSTIC`，windows.jsonl／summary.json／STATUS／JOBID。generanno，private／0GPU／2CPU／8GiB／10分钟；预期数十秒到数分钟，新增文本小于数MiB，不复制原数据。短诊断由主对话闭合；两条长链继续复用既有固定监控，不新增自动化。

## 解释与下一步

若大量保护窗仍有较多并集外区域，这只表明全窗保护与metadata空间范围有显著粒度差，支持继续设计具有正确语义的局部／候选监督；不证明忽略负监督导致精度下降，也不修改原A3或A17规则。若并集外空间很少，或主要由缺少精确位置的非ACGT造成，不以局部保护为主要新架构理由。完整raw-FASTA性能、逐物种和worst-species、背景负担、原生ANNEVO／Tiberius／Helixer以及独立适用范围仍是论文目标，当前计数不能替代这些证据。

原论文以公开方法名定向检索取得，只查询公开概念，未发送未公开代码或结果。Exa返回传输错误，未获得有效结果；依据项目允许的公开工具回退，用web检索和作者arXiv HTML核对两篇实际方法。OpenReview浏览器验证及CVF首次403不作为论文内容或阴性证据；不绕过登录或权限。只读上述必要方法段，不称穷尽性新颖性检索。写作按项目原有事实、推断与待验证事项分层；写作参考工具当前未暴露，未声称新取得作者样式样本。

## 实际结果与研究判断

CPU作业13389147经`sbatch --parsable sbatch/M28-A19-LOCAL-SUPERVISION-DIAGNOSTIC.sbatch`提交，Slurm root／batch／extern均COMPLETED、ExitCode=0:0、Elapsed=00:00:25；目录STATUS=COMPLETED、JOBID=13389147。两项解析检查和全部源记录核对通过，实际source fit job=13357566、末权重A3 step004608。脚本4.338秒，没有权重、DNA或特征tensor读取，也没有模型／GLM前向、搜索、更新或标签变更。

1491唯一窗中782有保护原因，10301个query为unknown；去掉缺少非ACGT精确位置的7窗后，775个空间可定位保护窗共19,046,400bp，登记保护跨度并集为2,876,072bp，之外为16,170,328bp／84.90%。9463／12400个固定初始bins不与并集相交；它们不是已经合法获得负标签的query。

下表按原每遍1536 draw报告，重复窗口按实际曝光保留；每物种768 draws／12288 queries。空间比例只在可定位的保护draw内计算，不把无保护窗口或noncallable缺少位置的窗口混入分母。

|物种|query正／已知负／unknown|有保护draw|空间可定位保护draw|保护并集外bp／空间分母bp|并集外比例|
|---|---|---:|---:|---|---:|
|拟南芥|2242／2831／7215|572|562|11745902／13811712|85.04%|
|水稻|786／8069／3433|237|237|4934173／5824512|84.71%|
|合计|3028／10900／10648|809|799|16680075／19636224|84.95%|

每遍pooled query unknown占10648／24576＝43.33%；拟南芥7215／12288＝58.72%，水稻3433／12288＝27.94%。保护原因可重叠，508次unknown_chain_span加586次complete_protected_nonpositive_chain不能相加作为保护窗口数；10次noncallable均来自拟南芥，对应7个唯一窗，所有计数保留。9759／12784个空间可定位保护draw的初始bins不与保护并集相交。三遍实际记录的监督计数与原A3 summary每遍一致，且原draw顺序和importance权重逐项一致。

结论是登记保护跨度通常只覆盖保护窗的一小部分，而存在性监督保护整个未匹配query集合，两物种都有这类粒度差。unknown query曝光比例的物种差异也确实存在，但这些计数不识别漏监督的因果效应，不证明它解释拟南芥性能，也不能把84.95%改写为可恢复的负监督比例。未知的真实基因和全窗续链支持仍构成必要限制。

下一步保留“模型可输出域与可判定监督域如何一致”的独立设计方向，优先考虑实际完整候选或有明确链语义的条件监督，不直接把A17的pooling bins当成可标负slots，不硬裁长链来取得更易通过的结果。A15已完成末head而自然读出／五gate未闭合，A17仍按原协议训练；先消费这两条实际机制结果，再决定匹配的完整raw-FASTA方案或独立下一实验，不重训旧失败配置、不事后换阈值。论文目标仍需对标三个原生模型、完整参考分母、背景及独立适用范围，此诊断只提供设计依据。
