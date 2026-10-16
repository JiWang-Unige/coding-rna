# 双头 CDS 支路的冻结 DEV 评价方案

2026-10-09冻结，2026-10-11补执行记录。双头支路只在原TRAIN上续训；其独立DEV评价沿用原固定输入、基线、参考分母及完整CDS链评价尺，不改已冻结的三目标实验。以下准备记录保留当时状态；满足TRAIN前置后已提交DEV，实际作业见末节。

## 固定候选与主要控制

唯一新候选为双头CDS支路末step_009216.pt，原C0参数从4608步增加4608次更新，新28系数从0步增加4608次更新；不使用中间checkpoint。主要机制控制为原STRUCTURED_CONTINUE末9216，保持结构目标、原draw顺序、标签、importance和追加更新次数相同。共同梯度裁剪包含新增28系数，因此不能宣称有效优化幅度完全相同。

完整报告家族包括原final4608锚点、三种目标final9216、双头支路final9216及ANNEVO／Tiberius／Helixer。所有成功完成候选全部并列，不按TRAIN排名预筛。进入新支路DEV前，四条TRAIN流程须完成各自末权重、全部更新核验、六染色体严格raw推理及端到端诊断；此顺序条件只检查完成身份，不读取分数或设TRAIN性能门槛。真实失败或未完成不计作零准确率，不静默移除候选。

## 固定推理与共同评价尺

DEV仍为拟南芥NC_003074.8／23459830 bp和水稻NC_089041.1／29936421 bp。完整原FASTA大写保留IUPAC供直接比较；caller在模型前把非ACGT转N，长度和坐标不变。使用冻结b048编码器与两个原生分类头、6144 bp编码、固定末权重角色系数、原nearest-center owner、双链整染色体exact exon-budget skyline、cap128、unit gain、null0、FP64及原tie policy。不加参考边界、softmask、RNA／蛋白证据，推理不读标签或feature cache。

M26 CDS-assessable longest-eligible-CDS／transcript-ID tie参考分母固定为AT5437、Rice2291、共7728 unique primary CDS链。预测端取每个原生gene的最长CDS／transcript-ID tie，complete_only=False，不按真值或本模型语法删基线预测；长链、重叠、motif及partial限制按[原DEV方案](full_grammar_dev_baseline_plan.md)保留。其他assessable isoform、phase和结构层另计，不替换primary exact chain。

沿用原六份成功ANNEVO／Tiberius／Helixer原生输出，不更新权重或operation point、不重新运行基线。新支路单独调用原共同评价函数，每次在读取参考前核验候选／基线／原FASTA直接相等、运行收据及实际Slurm终态。新结果路径独立，原三目标和final4608默认入口均保持不变。

## 固定结果与解释

每个方法报告逐物种TP／predicted／reference、exact-chain P／R／F1、CDS-base union、背景span／CDS覆盖和纯背景链／Mb、重复链、其他isoform、phase与描述性结构层。pooled按计数micro汇总，同时报告macro和worst species。另汇总全部八个方法；五个M28报告中复用基线的逐物种测量及pooled值须直接相等，参考分母和保留层须相同，否则先查输入／评价差异，不拼成比较表。

双头支路相对主要匹配控制、旧锚点、其他目标及每条原生基线的F1差值均按百分点报告，包含AT、Rice、pooled micro、macro和worst species。不用最佳物种、单点pooled优势或训练NLL代替完整结论，不按DEV选择head、阈值、上下文、loss、学习率、checkpoint或后续训练。单seed且只有两条历史复用开发染色体，不计算把gene当独立重复的确认性置信区间。

GENeranno及基线的上游物种／annotation暴露存在方法相关差异或未知，见[预训练暴露记录](../../../../paper/supplement/pretrained_exposure.md)。本支路显式复用已有CDS监督分类头，位置输入不是本项目的新模型创新。即使该开发比较提升，也不能据此宣称clean zero-shot、独立泛化、所有gene／UTR／isoform能力或Nature Communications论文成功；目标仍需要真正的竞争力、适用性与新颖性证据。

## 准备验证与后续执行

人工检查固定覆盖四TRAIN未关闭时不触碰DEV、真实driver接口的大写IUPAC保留、双头末权重／raw身份、先核输入再读参考，以及完整八方法汇总的分母／micro／macro／worst和基线一致性。只使用人工序列、临时元数据和mock，不读研究DNA／标签／checkpoint或旧DEV预测，不运行原生编码器。

待TRAIN关闭并完成结果处理后，新支路DEV raw为两个任务，每任务1RTX3090／2CPU／16GiB／2小时；CPU共同评价2CPU／8GiB／20分钟。完整家族汇总只读五份已完成比较报告，2CPU／4GiB／5分钟，不生成新预测。实际提交前再核队列、维护窗口和资源；当前这些是准备脚本而非已提交DEV作业。

独立输出为outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-DEV-RAW、NATIVE-CDS-PRIOR-DEV-COMPARISON及M28-FULL-GRAMMAR-CONTINUATION-DEV-FAMILY；第二项同样带M28-FULL-GRAMMAR前缀。汇总保存summary.json和comparison.md。sealed／held-out／Setaria不放行，无公开发布或Git提交。已运行的三条cohort仍由原固定监控对话和同一心跳跟进，本次准备不增加或修改监控对象。

## 人工准备验证实际结果

仅提交CPU准备检查13453287，private-teodoro-gpu、2CPU／4GiB／5分钟，无GPU，沿用实际成功的generanno环境。root／batch／extern均COMPLETED 0:0，UTC10:50:35–10:50:49、14秒。四项测试通过、0 failure／error：旧三目标或新支路TRAIN未关闭时不读DEV；driver保留大写IUPAC并调用固定末权重接口；新双头身份和预测先行／参考后置；完整八方法汇总拒绝遗漏候选、变化基线和错误micro算术，固定差值算术与人工计数一致。

outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-DEV-CHECKS中的JOBID、STATUS和summary.json齐全。全部测试只用了人工元数据、人工序列及mock：research DNA／labels、历史DEV receipt／预测、实际checkpoint和原生GLM均未读取或执行，optimizer_steps_this_run=0、DEV_used／sealed_data_used=false。四份新Python语法及四个SBATCH语法通过。证明的是准备接线与人工算术，不是真实末权重、DEV输入核验或模型准确率。

DEV RAW、COMPARISON及完整FAMILY脚本均只准备、未提交，正式输出目录未创建；不会因人工检查成功提前打开DEV。原新支路训练13453177和既有两条长cohort继续由同一心跳跟进。该同一监控在UTC10:48:27实际定时触发后追加到原对话，覆盖三cohort、gpt-5.6-luna／low，10:48:46完成DONT_NOTIFY；并未创建新监控或添加本次已消费短检查。

## 2026-10-11实际DEV提交与启动

原三个目标与双头CDS支路的末权重、更新核验、六条TRAIN raw和端到端评价均已完成；主对话核对Slurm终态、运行身份、固定22688参考分母和计数算术，并完成结果处理。双头链路13453177→13453180→13453182_0–5→13453184全部root／batch／extern COMPLETED 0:0，最终评价于10-10 UTC01:41:44结束。此项顺序前置按完成身份判断，不因性能增删候选。TRAIN结果已写入[结果工作稿](../../../../paper/manuscript/M28_full_grammar_results_working.md)。

10-11 UTC12:49:39实际提交两个DEV raw任务13493570_0–1，完整数组afterok比较13493571，以及afterok13493571的原八方法汇总13493572。三个原目标DEV比较13465375／13465377／13465379和final4608比较13449926已经完成，汇总仍须读取全部五份固定比较，不只依赖最新候选。每个raw沿用单RTX3090／2CPU／16GiB／2小时，数组最多并行2；仅路由扩至private-teodoro-gpu,shared-gpu，科学命令不变。比较2CPU／8GiB／20分钟，汇总2CPU／4GiB／5分钟，均无GPU。

UTC12:50:47启动核对：两个raw均RUNNING，gpu034；实际numeric Slurm roots分别13493574／13493570，从12:49:40启动。两个下游CPU均PENDING有效Dependency。这只是调度启动，不是完成推理、指标结果或监控定时回传。RAW的JOBID、started及outer summary按actual SLURM_JOB_ID核对，而不是数组ID_task。

完整八方法结果独立写入原预定输出路径，不覆盖三目标单项报告。原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同一心跳coding-rna-r5-light-monitor供此批剩余DEV流程复用，不创建新对话／自动化。新批监控配置、实际首次检查和定时追加分别取证，不能从提交成功推定。未新增训练、修改冻结评价、放行sealed／held-out／Setaria或公开发布。

## 完整八方法 DEV 终态与结果

2026-10-11 UTC，native raw13493570_0–1的实际数字roots13493574／13493570于13:20:27／13:27:46结束；比较13493571于13:27:46–13:29:38、八方法汇总13493572于13:29:38–13:29:58执行。全部root／batch／extern为COMPLETED／0:0，STATUS／JOBID、两项输入和GFF及内外summary、汇总summary／comparison均齐全、非空并匹配实际作业。

末9216、b048冻结编码器和两原生头、28系数、6144编码、双链全染色体、单位15角色gain／null0／cap128／FP64均符合固定接入。推理reference／label／feature cache=false、DEV=true／sealed=false、更新0。已完成共同评价在读参考前直接核对输入和原生操作点，之后核对真实DNA、phase、cap和GFF roundtrip。主对话独立核算八方法逐物种／pooled链与base P／R／F1、背景、macro／worst、全部匹配差值和结构层；五份M28报告的参考层及三个原生方法完整测量直接相等。没有按支持性删参考或基线预测。

| 固定方法 | AT F1 % | Rice F1 % | pooled micro F1 % | macro F1 % | 背景链 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原 final4608 | 50.8442 | 30.0960 | 44.6085 | 40.4701 | 64 |
| STRUCTURED_CONTINUE final9216 | 62.9235 | 36.9312 | 54.7508 | 49.9273 | 79 |
| LOCAL_CHANNEL final9216 | 44.4851 | 36.8666 | 42.1812 | 40.6758 | 173 |
| HALF_STRUCTURED_HALF_LOCAL final9216 | 51.8012 | 42.4519 | 48.9426 | 47.1265 | 245 |
| NATIVE_CDS_PRIOR final9216 | 69.2972 | 51.5164 | 63.7761 | 60.4068 | 121 |
| ANNEVO | 76.2651 | 74.4318 | 75.7025 | 75.3485 | 146 |
| Tiberius | 74.1771 | 66.5822 | 71.8047 | 70.3797 | 165 |
| Helixer | 67.4564 | 53.5243 | 62.8512 | 60.4903 | 563 |

参考仍AT5437／Rice2291／pooled7728。双头TP／unique predicted：AT3653／5106、Rice1223／2457、pooled4876／7563；pooled precision64.4718%、recall63.0952%，CDS-base F1 93.6717%。所有方法worst species为Rice；双头duplicate0、其他assessable isoform命中440另计，两个物种全部prediction records的phase均一致。长span参考AT1／Rice2均保留而未命中。

相对同更新数STRUCT，双头AT／Rice／pooled F1提高6.3737／14.5852／9.0253个百分点；相对Helixer为+1.8408／−2.0078／+0.9249，macro低0.0835个百分点。双头在两个物种仍低于ANNEVO／Tiberius。该结果不支持全物种优越，新增28系数和复用注释预训练CDS头也不自动构成新颖性或有效优化匹配。单seed、历史复用DEV、不依DEV调参／选权重和独立泛化限制不变；论文goal未完成。

来源为`outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-DEV-COMPARISON/summary.json`及`outputs/M28-FULL-GRAMMAR-CONTINUATION-DEV-FAMILY/{summary.json,comparison.md}`。五份M28比较来源job为13449926／13465375／13465377／13465379／13493571。另有预定fresh／joint十方法扩展，但原八方法summary和测量不改变，见[扩展终态](full_grammar_sampling_context_dev_plan.md#十方法扩展-dev-终态与结果)。

原监控13:40:58 UTC完成六raw＋三比较＋两汇总整批结构核验并向原主对话成功交接一次，随后暂停同一心跳。主对话实际观察scheduled turn01a12b31-c7ad-7b41-9e56-604cdb9a0c80于13:40:45–13:41:57运行，原watcher保持idle，保存类型heartbeat、原target／15分钟／通知策略／created_at未变。本节指标是主对话科学处理，不是监控解释；该DEV批次不再轮询。未新增训练、封存数据访问或对外发布。
