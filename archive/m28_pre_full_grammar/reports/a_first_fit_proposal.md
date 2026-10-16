# M28 A gene query首训提案

2026-10-01起草。2026-10-02用户先明确回复“批准，这种问题不要再询问我，直接开始迭代运行”，随后要求“不要考虑成本问题，尽快发文章给科学界提供新的模型才是更重要的事情”。本轮单一A配置、一次拟合→完整DEV推理→评价及结果分析已授权；经济成本门槛与人为总资源/产物容量预算现已取消，不重复询问同范围常规步骤。模型、TRAIN曝光、唯一末checkpoint、数据隔离和科学评价条件不变。三阶段现已COMPLETED/0:0，缓存权重/完整DEV评价完成；原pilot继续条件因两物种背景span增加而未通过，详见[a_result.md](a_result.md)。这不是论文或独立泛化成功。检验从头训练的gene-query完整方案，跳过B1固定候选池；B2负结果不证明A有效或查询输出形式本身有新颖性。

## 已知事实与当前实现状态

A的纯合成指针、DNA语法、矩形匹配、批量teacher成本及query到公共合并/GFF3的接口已核验；2026-10-01工程job13319189及batch均COMPLETED/0:0，6项新人工DNA跨接口测试通过（0.433秒），仅6CPU秒、0GPU，无真实数据或optimizer。既有GFF解析器未显式关闭读句柄，测试产生ResourceWarning但坐标/phase回读通过，未改冻结评价器。TRAIN元数据盘点见outputs/M28-A-QUERY-TRAIN-CENSUS/summary.json。自然TRAIN正目标每窗最多11条，原paired draws最多10条，最多37个exon。原draw的两物种元数据保护为567/768和237/768；这些只是下界，未检查全窗DNA。不能据此把其余全部query定为生物学负例。

新增GeneQueryModel及joint_query_losses已接通从头初始化的共享模块、局部分支和query目标，外部GLM/DNA张量保持无梯度；balanced_query_components实现下文类均衡存在性和完全一致的匹配成本，unmatched_query_policy保守保护未知整窗且保留正例。旧GeneQueryPointer.for_window仍detach外部输入，matched_query_components仍是普通存在性BCE，不能混用两条目标路线。run_query.py已实现prepare/fit/infer代码入口：TRAIN预筛与R4重建目标、固定计时选择与历史成本报告、从头三遍原draw、末checkpoint读取及完整DEV grid到query记录。query_output.py连接已完成query与公共owner/去重/selector；annotate_query.py与evaluate_query.py已补齐原始FASTA编码和冻结统一评价代码。prepare接线已于2026-10-02在真实TRAIN/cache/CUDA完成固定8窗forward/loss/backward/无参考query，job13331164及batch COMPLETED/0:0、optimizer更新0；原成本门槛当时未通过，故当时未提交后续阶段；用户随后取消成本阻断，正式拟合13344559、完整DEV推理13344568、评价13344571现均完成。实际训练权重与完整DEV指标已核验，但原生GLM/FASTA仍未执行；详情见[a_execution.md](a_execution.md)和[a_result.md](a_result.md)，不把缓存DEV完成称作原生FASTA端到端或科学PASS。

2026-10-01合成工程核验：job13320418中9项裸接口及8项普通匹配通过，新测试7/8通过；一项断言误用places=6，与事前绝对容差1e-6不一致，差值9.5367431640625e-7。使用原一次源码修复机会改为显式delta=1e-6，不改模型/目标/fixture或容差；原失败完整保留。job13320471及batch COMPLETED/0:0、STATUS=COMPLETED，8项新联合测试全部通过（0.352秒），含query-only梯度、独立穷举完整目标/null常数和标签保护。两轮合计35CPU秒、0GPU、4986 bytes运行产物；A准备总148CPU秒=0.041111CPUh、19207 bytes，修复机会已使用。未训练D2048/Q16/H128模型实测583153可训练参数、optimizer更新0；合成PASS不是训练、准确率或完整FASTA能力。

run_query入口合成工程job13321036及batch COMPLETED/0:0，STATUS=COMPLETED，7项新测试通过（1.108秒），覆盖新接线的TRAIN/cache身份和标签隔离、计时类别/预算门槛、末A checkpoint身份、R4重建目标→shared/loss/backward及shared→无参考query输出。仅内存人工metadata/DNA与未训练张量，0GPU、0 optimizer更新，无真实cache/checkpoint或TRAIN/DEV读取；既有未变测试不重跑。Slurm实际45CPU秒、MaxRSS976864KiB、918 bytes运行产物；截至本入口阶段A准备累计193CPU秒=0.053611CPUh、20125 bytes，仍低于300CPU秒/1MiB。生产命令代码是python scripts/experiments/M28-METHOD-RESTART/run_query.py --mode prepare|fit|infer，三种模式该合成工程阶段均未调用。该工程阶段当时prepare/fit/infer目录、原10.166667GPUh/4CPUh/2GiB预算及停止条件未变，首训待批准；这是历史记录，现行取消成本门槛的授权见下文。

2026-10-02评价/FASTA接口合成工程job13321229及batch COMPLETED/0:0、STATUS=COMPLETED，8项新测试通过（0.580秒）。覆盖人工FASTA两链→替身encoder→未训练小头→owner/GFF、原四块six-mer/BF16编码对齐、A完整scope/checkpoint/阈值来源、两物种整数F1/背景gate、模拟参考loader保留unsupported与GFF phase回读。未执行生产CLI或读取真实TRAIN/DEV/cache/checkpoint，0GPU、0 optimizer更新；冻结解析器再次出现已知ResourceWarning，回读通过，未改它。实际17CPU秒、batch MaxRSS828152KiB、1428 bytes运行产物；A合成准备累计210CPU秒=0.058333CPUh、21553 bytes，仍在300CPU秒/1MiB内，未新增修复或重提。

evaluate_query.py主分母保持7728、旧6450副表，不按A支持条件删参考；两物种均比完整C0提高F1且背景span不增加才通过原pilot继续条件，不是论文/独立泛化成功。annotate_query.py的原生FASTA入口仍要求显式数据和匹配末A checkpoint，本轮未执行；原生GLM/FASTA的数据范围不因取消成本门槛而自动扩展。真实TRAIN计时53898.220694080秒的旧拟合投影及budget_gate=false完整保留，不改写为通过；它们不再控制本轮是否启动，也不再追加warmup计时来证明预算。run_query.py保留完成身份校验与实际产物大小报告，取消成本gate和2GiB人为中止；除这两项工程控制外，训练与科学配置不变。两项针对变更的测试在正式FIT作业内通过（0.036秒），未重跑其他未变测试。

## 推荐的单一模型

使用原SharedContext：冻结GLM六聚体2048维输入、hidden128的双向GRU及DNA卷积分支；从头初始化共享上下文和15类局部分支，不继承C0、B1或B2监督权重。共享部分的初始参数按原C0 seed0构造顺序逐张量一致，然后初始化query头。推荐Q=16、hidden128，只做这个固定容量，不作Q扫描。Q16高于已观察的TRAIN下界11，预留5个null位置；这不是DEV或其他物种的容量保证。

query损失必须能反传到本次从头训练的共享上下文。GeneQueryModel入口明确detach外部GLM/DNA输入，内部query窗口保留新h的计算图；合成query-only损失已给新shared非零有限梯度，外部输入无梯度，未知query行不获监督梯度。正式A训练接线仍须使用这一路线，不能误调用旧冻结输入的裸接口。没有额外stop-gradient消融或可选模式。每次重新建立window状态，匹配no-grad成本不用于最终反传。

## 推荐的标签保护

逐窗调用现有R4 build_targets，不改primary选择、训练资格或draw身份。全部positive_chains必须一对一匹配；N>16、链ID/坐标错配或支持合同冲突即停止，不按分数删目标、换更短isoform或增加Q。

未匹配query采用一个可检查的保守窗口政策。如果有unknown_chain_spans、任一完整protected_chain不在positive集合中，或任一有效DNA碱基不callable，则所有未匹配query为-1；否则为0。正query始终按明确positive链监督，不因窗口含未知区就删正例。此政策比B1的逐候选跨度判负更保守，可能损失负监督；这是首轮方案的明确代价，不是已证明最优政策。不存在输出依赖的query判负、参考注入自由解码或把缺失注释称生物学真负。

N=0且unmatched=-1时query目标是可微零，仍有可用局部监督；N=0且unmatched=0时所有query为参考相对负例，不制造边界/EOS目标。训练日志按物种报告保护原因、正/负/未知query数量和实际监督质量。metadata保护804/1536只是当前draw的下界，DNA读取后可以增加保护，但不能改positive资格或重抽样。

## 推荐的损失与完全一致的匹配目标

每窗目标为importance_weight乘以local_loss、存在性BCE、边界CE、EOS BCE之和，四个系数固定1。local_loss沿用fine/grouped likelihood；边界/EOS分别先对每个正gene的事件取均值，再对N个正gene取均值。teacher使用原有有序半开坐标目标，推理另受DNA语法限制，不称DNA条件化似然。

存在性项按窗口内已观察正/负类各占一半质量；只出现一个已知类时该类占全部质量，-1没有贡献。因此不是把Q16直接代入当前普通BCE成本。设p_i=softplus(-z_i)、n_i=softplus(z_i)，匹配矩阵中存在性项应为：

- unmatched=0且0<N<Q：每个gene到query i贡献0.5 p_i/N减0.5 n_i/(Q-N)，其余全负项是分配不变的常数。
- unmatched=-1，或N=Q只有正类：贡献p_i/N。
- N=0：无需分配，按上述空窗目标直接计算。

每个矩阵单元另加该链的boundary_mean/N和EOS_mean/N。Hungarian离散分配后，用当前h和被分配链重新计算可反传的完整目标。新增合成测试已穷举Q3、N1/2/3和unmatched为-1/0的全部分配：所选分配达到完整类均衡目标最低值，每个分配的成本加null常数与完整目标在固定绝对容差1e-6内一致，参考逆序最低损失一致；不能沿用旧delta/Q却改训练loss。query不足或非有限成本必须停止。此验证只支撑合成实现的目标一致性，不证明真实训练或科学收益。

## 数据曝光与比较含义

TRAIN染色体仅为拟南芥NC_003070.9、NC_003071.7、NC_003076.8，以及水稻NC_089043.1、NC_089044.1、NC_089046.1。沿用R4标签与原1536 paired draws、原顺序和importance weight，重复3遍，batch1，4608更新；seed0、FP32头、AdamW lr3e-4、weight_decay0.01、betas(0.9,0.999)、eps1e-8、global clip1.0和固定学习率。所有checkpoint保留，但唯一主评价状态为step004608，不用DEV选epoch。

复用outputs/M28-PILOT-R4/features的已有1491个唯一TRAIN和8690个DEV窗口缓存，只读引用，不复制159.139GiB，不增加数据或预训练暴露。revision仍为b0483c23b6b63787b61a6d3a204a9b517d6ba345。缓存内旧标签不用作监督，按ID连接R4标签；训练catalog先按TRAIN预筛，不能为方便解析DEV参考记录。

原draw三遍提供9084次完整正链曝光但仅2813个不同primary转录本ID；有放回与重复遍数不是独立样本。C0、B1的同曝光完整末checkpoint/GFF3作为已运行对照，不为本轮重训它们。A与C0/B1是完整方案比较，包含query、辅助目标、梯度路线和更保守负监督等差异，不能单独归因于输出格式或宣称等GPU耗时因果消融。B2额外拟合过新头，单列历史诊断，不作A的曝光匹配主对照。

## 无参考输出与固定操作点

推荐τ=0，gain为query存在性logit减τ，不等于自然gene正确概率。每个query使用全部DNA合法事件、单一greedy历史，max_exons=128。沿用连续ATG、GT/GC–AG、intron至少5bp、真实codon carry、无内部STOP/未知编码碱基和连续末端STOP；拆分起止、非规范剪接和其他遗传密码仍不支持。模型未结束的链不能修补为完整gene。不设置B1式endpoint/link/beam限额，不借参考选择更优query。

所有完成链含非正gain先保存为free记录；公共规则按全局跨度中点最近窗口中心归属，平局较小start，完全同链取最大同臂gain，随后最大正增益同链向非重叠选择，0选择null。保留raw query、owner损失及冲突阶段，不改变原评价的重叠参考分母。只输出primary protein-coding CDS链/GFF3，不冒称UTR、全部isoform或全局最优搜索。

## 执行资源与历史计时（不作成本门槛）

| 阶段 | 当前Slurm资源/单作业时窗 | 实际执行 |
|---|---|---|
| 固定TRAIN接线计时 | 1 RTX3090 24GiB、2CPU、16GiB、10分钟 | 13331164 COMPLETED/0:0，110秒、0更新；历史budget_gate=false |
| 唯一正式拟合 | 1 RTX3090 24GiB、2CPU、16GiB、4小时 | 13344559 COMPLETED/0:0，2606分配秒 |
| 完整DEV 8690窗无参考推理 | 1 RTX3090 24GiB、2CPU、16GiB、6小时 | 13344568 COMPLETED/0:0，5033分配秒 |
| 合并、GFF3回读和统一评价 | 0GPU、2CPU、32GiB、2小时 | 13344571 COMPLETED/0:0，84分配秒 |

以上终态于2026-10-03 Asia/Bangkok结果处理时核实；单作业walltime是当前调度设置，不是经济预算或再次请求批准的门槛。取消原10.166667GPUh、4独立CPUh及2GiB人为合计上限，仍记录实际资源和产物大小。真实VRAM、磁盘空间、调度约束及科学一致性错误仍须处理。提交前按coding-rna-baobab技能检查实时Slurm；输出使用新目录，不覆盖冻结资产。

已完成计时的固定样本为每物种四类原TRAIN draw：正目标最多、总正exon最多、零正且无metadata保护、正例且有metadata保护；同分取最小原draw ordinal。八类/八唯一窗口只使用未训练权重，完成forward、匹配、loss/backward和无参考query，optimizer更新0，不读DEV。它们是工作量/标签保护检查，不是准确率选择。

历史公式max(full step)×4608+1800产生14.971728小时投影，超过旧4小时门槛；decode分量投影52.390194分钟。原summary及逐窗计时保持原值。首窗冷启动被重复外推的解释仍是工程推断，没有进行同窗warmup确认；用户取消成本门槛后无须另跑计时或以改公式得到“预算PASS”。这些成本数字不构成A的科学结果。

## 完成条件与下一科学判断

只有4608更新、8690窗完整推理、保存范围身份一致、无参考自由生成和GFF3回读通过，才计算A正式结果。主分母固定7728（拟南芥5437、水稻2291）；长链、训练/语法不支持参考全部保留，旧6450副表单列。每物种与macro/micro报告exact-chain P/R/F1、最终CDS union、实际free/owner/selected链恢复、Q饱和/unterminated、背景span/CDS覆盖、全背景链/Mb和真实资源。背景指标不是生物学真FP，局部head与最终CDS union不混称。

推荐继续条件是在两物种均相对完整C0改善exact-chain F1，且两物种背景span覆盖都不增加。C0当前同口径F1为44.3290%/36.7338%，背景span覆盖约0.6698%/1.0026%；实际gate用原summary计数算出的完整精度，不用这里四舍五入值。它是进入下一独立迁移/方法消融的pilot条件，不是SOTA或期刊成功标准。不通过就保留当前配置负结果，不扫τ、Q、seed或追加epoch救结论；若未完成则记未检验充分，与科学FAIL分开。

ANNEVO、Tiberius、Helixer以同口径已有缓存展示且注明未重新原生运行，不称最新版本独立重现。拟南芥/水稻已用于开发；backbone暴露未知，不称独立泛化或干净zero-shot。test/Setaria继续封存。即使A通过pilot，也仍需原始FASTA端到端交付、独立未用于选择的物种、冻结对标版本与数据、机制消融和稿件/作者决策，不能把这份提案或工程PASS当Nature Communications证据。

## 授权与停止

用户于2026-10-02在原主对话01a0332e-e0a9-7aa2-90f9-159405701833先批准这轮单一A拟合/DEV评价/结果处理，再明确取消成本考量并强调尽快提供科学界可用的新模型。当前直接推进同配置的实施、验证和结果分析；不因进入常规下一步或成本估算再请求批准。成本指令不开放sealed/test/Setaria，不授权改变冻结科学规则、恢复Auto Research、对外发布或取消他人作业；原生FASTA及独立物种的数据协议仍需具体确定。

非有限loss/梯度/分数、错配、N>Q、冻结输入被更新或自由推理读参考时停止受影响阶段并保留证据，不能用半成品排名。超时/OOM/磁盘不足属于实际执行失败，报告事实后在现有授权和不改变科学配置的范围内处理；不把旧成本gate或人为2GiB当停止原因。工程修复不自动授权改seed、曝光、操作点、追加epoch或参数扫描。科研pilot条件仍按原协议判断，不能为加速投稿把探索结果包装为确认性成功。
