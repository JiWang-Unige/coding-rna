# 新抽样与联合编码的冻结 DEV 比较方案

2026-10-09 UTC冻结，2026-10-11补执行记录。原TRAIN上新增的新抽样和24576联合编码各有唯一末9216候选，本方案在两者正式结果可用前固定其后续开发集比较。原[八方法家族](full_grammar_native_cds_prior_dev_plan.md)及其入口、候选、结果不变；本次是另行记录的两个预定对照，不按TRAIN排名增删候选，不回写原家族。以下准备记录保留当时状态；满足六条TRAIN前置后已提交DEV，实际作业见末节。

## 固定候选和机制控制

新抽样候选固定为`outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-FIT/step_009216.pt`，原C0从final4608及保存AdamW继续4608次更新，使用已经固定的seed1均匀抽样、4608 draws／4219 unique、importance1。24576联合编码候选固定为`outputs/M28-FULL-GRAMMAR-JOINT24-FIT/step_009216.pt`，相同起点和追加更新数、原1536 draws重复三遍／1491 unique及原importance，每个原24576神经读取窗口一次原生共同编码。两者都不选择中间checkpoint，不混合这两个改动，也不在DEV后改上下文或抽样。

各自的主要匹配控制是已固定的STRUCTURED_CONTINUE末9216。联合编码相对该控制只改变编码器可见范围，读取窗口及grammar上下文不同概念须分开表述。新抽样相对该控制保持预期均匀目标，但改变抽样提议、有限训练曝光与梯度方差，不能把差异单独归因为覆盖量。原final4608是少4608次更新的不变锚点，不是匹配更新次数的因果控制。

进入这两个DEV推理前，原三个目标、双头CDS、新抽样和联合编码共六条新TRAIN分支须完成末权重、全部更新核验、六染色体严格raw推理、端到端评价及主对话结果处理。顺序检查只判断完成身份，不设性能门槛或按分数排除分支；真实失败、超时或未完成明确报告，不能记作零准确率或隐藏。保存失败证据后按已经声明的执行／实现修复规则处理。

## 相同输入和共同评价尺

DEV范围仍仅拟南芥NC_003074.8／23459830 bp与水稻NC_089041.1／29936421 bp。保留原大写IUPAC序列做输入逐碱基直接比较，caller把非ACGT映射N但不改变坐标。新抽样用原6144编码；联合编码用原b048 GENERanno、原SDPA flash、4096六碱基tokens、BF16隐藏状态后FP32头，短窗口按原定向后右补N至24576。backbone冻结，不使用native CDS头；两者均保持相同raw网格、双链nearest-center ownership、整定向染色体exact exon-budget skyline、cap128、unit gain、null0、FP64和原tie policy。参考、标签、RNA／蛋白证据、softmask和特征缓存不进入caller。

参考仍为M26 CDS-assessable longest-eligible-CDS／transcript-ID tie，AT5437、Rice2291、总7728 primary链；unsupported参考层保留。预测仍按原生gene取最长CDS／transcript-ID tie，complete_only=False，不按真值或M28语法过滤基线。精细phase、其他assessable isoform及结构层另报，不替换primary exact-chain口径。详见[原冻结共同尺](full_grammar_dev_baseline_plan.md)。

ANNEVO／Tiberius／Helixer复用原六份成功原生输出和原operation point，不重跑、不更新权重。每个新比较在读取参考前先核实际checkpoint、运行收据、Slurm终态及候选／基线／原FASTA直接相等。现有原final4608、三目标和双头结果均不覆盖。

## 完整报告和解释

两个新增候选均报告AT、Rice、pooled micro、macro和worst-species exact-chain P／R／F1及原计数，另报CDS-base union、背景span／CDS覆盖、纯背景链／Mb、重复链、其他isoform、phase和描述性结构层。分别与STRUCTURED_CONTINUE、旧4608及三原生基线比较，差值按百分点报告；失败或未完成与负性能结果区分。

原八方法报告保留原义；另一个扩展比较表可并列原八方法和这两个预先指定候选，共十方法。它必须完整保留候选、相同参考分母及复用基线逐物种计数和汇总值，不能只挑两新候选中的赢家，也不能把新增方法当作原八方法冻结时已有的确认性设计。表内任何复用基线或输入差异先查清再汇总，不将不同评价来源拼成一个排名。

这是单seed、两条历史复用同物种开发染色体的描述性比较，不把gene当独立生物重复计算确认性置信区间。上游物种及annotation暴露限制仍适用，见[预训练暴露记录](pretrained_exposure.md)；开发表现不等于clean zero-shot、独立泛化、完整gene／UTR／所有isoform能力或论文成功。DEV结果不用于训练损失、抽样、head、阈值、上下文、学习率、checkpoint或后续训练选择。sealed／held-out／Setaria不放行，无公开发布或Git提交。

## 后续实施

已准备两个候选共用的直接raw入口、共同尺评价入口、独立扩展汇总及人工绑定检查，原默认入口和八方法汇总不变。每候选两个raw任务、最多并行2，每任务单RTX3090／2CPU／16GiB／2小时；共同评价2CPU／8GiB／20分钟，扩展表2CPU／4GiB／5分钟。实际研究DEV raw／评价／扩展表三个SBATCH仅准备、尚未提交，输出目录未创建；仅人工CPU检查已提交，执行研究DEV时须再次核队列和维护。

入口分别为`infer_full_grammar_sampling_context_dev.py`、`evaluate_full_grammar_sampling_context_dev.py`及`aggregate_full_grammar_sampling_context_dev.py`，前两者显式指定FRESH_DRAWS或JOINT24576。raw输出在`outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-RAW/<branch>/<seqid>`，比较在`M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-COMPARISON/<branch>`，十方法扩展在`M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-FAMILY`，后两项同样带outputs前缀。扩展入口复用原八方法的直接验证与测量，但只写新路径，不更新原八方法文件。

CPU人工检查单独使用`M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-CHECKS`，2CPU／4GiB／5分钟、无GPU、沿用实际generanno环境。五项测试只用临时人工元数据、序列及mock：原四分支或任一新TRAIN未关闭时不读DEV；两个raw driver保持大写IUPAC并调用各自固定末权重；新checkpoint身份及输入先于参考；全部十方法的micro／macro／worst、差值、固定基线和参考分母；实际扩展写出summary／comparison但不动原八方法文件。未读研究DNA／参考／历史DEV收据或预测，不加载checkpoint或执行原生GLM。这类人工接线验证不证明真实末权重、DEV输入相等或科学准确率。

当前唯一固定监控只检查已经提交的TRAIN长流程，不能把本方案中的未提交DEV对象加入心跳。新DEV任务将按同一固定监控和同一自动化复用规则处理，不新建监控对话或重复自动化。

## 人工检查执行

CPU检查作业13453843通过`sbatch sbatch/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-CHECKS.sbatch`执行，private-teodoro-gpu／2CPU／4GiB／5分钟、无GPU，原generanno环境。实际UTC12:26:39–12:27:40、61秒、gpu034；root／batch／extern均COMPLETED 0:0，五项测试通过、0 failure／error，人工用例本身约0.03秒。四份新Python和四个SBATCH语法通过，运行中训练源代码与原八方法文件未改。早期RUNNING快照没有完成报告，随后以实际终态和summary核实成功，未因观察超时或导入启动时间重提。

输出`outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-CHECKS`中的JOBID、STATUS、summary.json齐全。实际测试验证了两个新TRAIN完成身份与原四条前置、各自末权重的driver参数、大写IUPAC、参考读取前的输入接线，以及十方法的参考分母、重复基线测量、micro／macro／worst和匹配控制差值算术。新扩展真实写出临时人工summary／comparison、原八方法人工文件未变。测试只用了人工序列、临时元数据及mock：research_DNA_or_labels_opened／historical_DEV_receipts_or_predictions_opened／checkpoint_loaded／actual_native_GLM_forward／DEV_used／sealed_data_used均false，optimizer_steps_this_run=0。它不证明真实研究输入已核对、末权重准确率或原生DEV执行成功。

主核对研究DEV raw／comparison／family三个新输出根均不存在，三个研究SBATCH未提交；不能把这个名称含DEV的人工CPU检查当作放行DEV。原固定心跳继续原四条TRAIN长cohort，不把已消费的短检查或未提交DEV对象加入，不重建或重复监控。

后续实际DEV raw由同一预备SBATCH分别带FRESH_DRAWS、JOINT24576参数提交，两个任务各自完整数组afterok后调用相应比较；扩展表需全部六条续训DEV比较和既有final4608锚点报告齐全，不只依赖两个新增比较。所有实际作业ID将在真正提交时记录，不从这些预备命令推定已执行或数据已放行。

## 2026-10-11实际DEV提交与启动

原三个目标、双头CDS、新抽样和联合编码共六条TRAIN前置已完成末权重、全部更新核验、六染色体raw、同尺端到端评价及主对话结果处理；各Slurm root／batch／extern均COMPLETED 0:0。新抽样13453573→13453581→13453582_0–5→13453583的最终评价于10-10 UTC02:43:59结束；联合编码13453777→13453778→13453779_0–5→13453780于02:01:26结束。主对话核对原22688参考分母、输入身份、逐染色体／逐物种／pooled计数和macro／worst，不以分数筛选候选。结果见[结果工作稿](../../docs/M28_full_grammar_results_working.md)。

10-11 UTC12:49:39–40实际提交FRESH_DRAWS raw13493573_0–1及完整数组afterok比较13493576；JOINT24576 raw13493577_0–1及完整数组afterok比较13493578。两者都用各自固定step_009216.pt，不选择中间权重。独立十方法扩展13493579依赖afterok13493572:13493576:13493578；13493572为原八方法汇总，且扩展程序仍核对全部七份固定M28比较。未改变原八方法，也未加入combo／complete-core为DEV候选。

每个raw单RTX3090／2CPU／16GiB／2小时，每候选最多并行2，实际路由private-teodoro-gpu,shared-gpu，科学命令和推理参数不变。每比较2CPU／8GiB／20分钟；扩展2CPU／4GiB／5分钟，CPU后处理均无GPU。UTC12:50:47四个raw均RUNNING：新抽样gpu035／gpu021、实际numeric roots13493575／13493573，12:49:40启动；联合编码gpu017、roots13493580／13493577，12:49:43启动。两个比较及十方法扩展均PENDING有效Dependency。启动不等于完成推理或获得指标，JOBID／started／outer summary须与actual SLURM_JOB_ID一致。

这两项与双头DEV、八方法汇总合为一批剩余冻结DEV流程：六个raw、三个比较、两个汇总全部闭合后才能整批交接。复用原监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3与同一心跳coding-rna-r5-light-monitor；不新建对话／重复自动化，不轮询已完成旧TRAIN作业。监控配置、实际首次检查、真实定时追加和终态回传分别记录，当前启动记录不证明后三者。无新训练、DEV调参、封存数据或公开发布。

## 十方法扩展 DEV 终态与结果

2026-10-11 UTC，新抽样raw13493573_0–1的实际roots13493575／13493573于13:21:45／13:27:56结束，联合编码13493577_0–1的roots13493580／13493577于13:20:53／13:28:14结束。比较13493576执行13:27:56–13:29:38，13493578执行13:28:25–13:29:38；独立十方法扩展13493579执行13:29:58–13:30:12。所有root／batch／extern均COMPLETED／0:0，STATUS／JOBID和全部非空收据、GFF及summary／comparison匹配。

原八方法汇总13493572保持原文件与完整测量；扩展13493579直接复用并附加两个早先固定final9216，不加入combo／complete-core或删去负支路。原八方法顺序和measurements与十方法前八项逐项相等；七份M28比较的三原生测量、参考层、AT5437／Rice2291／pooled7728分母相同。主独立核算全部十方法链／base PRF、背景／isoform合并、micro／macro／worst、所有差值及60项描述性结构PRF，全部一致。20个物种×方法的预测记录phase均一致，duplicate0。长span参考AT1／Rice2保留，七个M28均未命中，不据此推断长gene能力。

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
| FRESH_DRAWS final9216（独立扩展） | 67.4379 | 51.9087 | 62.6852 | 59.6733 | 188 |
| JOINT24576 final9216（独立扩展） | 70.4990 | 57.5017 | 66.6229 | 64.0003 | 121 |

新抽样TP／unique predicted：AT3585／5195、Rice1217／2398、pooled4802／7593；pooled P63.2425%／R62.1377%，base F1 93.9940%。联合编码为AT3765／5244、Rice1305／2248、pooled5070／7492；P67.6722%／R65.6056%，base F1 95.0958%。两个分支其他assessable isoform命中413／415另计，不充作primary TP。所有方法worst species均为Rice。

相对STRUCT，新抽样AT／Rice／pooled F1提高4.5145／14.9775／7.9344个百分点，但在两个物种仍低于三个原生基线。联合编码相对STRUCT提高7.5756／20.5704／11.8721个百分点；相对Helixer提高3.0426／3.9774／3.7717，但相对Tiberius低3.6781／9.0805／5.1818，相对ANNEVO低5.7661／16.9302／9.0796个百分点。

联合编码比Helixer少42个pooled TP、少1047个预测，precision高7.8057个百分点而recall低0.5435；AT净TP多92、Rice少134。完整链F1超过Helixer不能称两物种recall都更好。纯背景链121比Helixer少442，但比STRUCT多42；共同背景区域27,956,441 bp内span FPR为0.6851%，不能把121条纯背景链当全部非exact错误。基线和候选的其余isoform、结构及背景定义均保留。

来源为各branch `outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-COMPARISON/<branch>/summary.json`及`outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-FAMILY/{summary.json,comparison.md}`。原八项来源见[八方法终态](full_grammar_native_cds_prior_dev_plan.md#完整八方法-dev-终态与结果)，新增比较job13493576／13493578。完整结果同步到[结果工作稿](../../docs/M28_full_grammar_results_working.md#complete-frozen-dev-family-and-predeclared-extension)。

原监控13:40:58 UTC自动完成整批结构检查并向原主交接一次，同一心跳暂停供新任务复用；终态自动回传已实际观察，不再写作未验证。单seed、历史复用DEV、上游annotation暴露限制不变；不据此选模型／checkpoint／上下文或改变下一轮训练，不声称独立泛化、机制因果、新颖性或论文目标完成。新的两模型TRAIN评分回放方案在主读取这些新指标前已经固定，DEV不是其设置来源。无sealed／held-out／Setaria访问、公开发布或Git提交。
