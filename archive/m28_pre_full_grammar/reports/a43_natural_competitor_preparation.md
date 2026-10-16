# 自然竞争结构监督的准备方案

2026-10-07。A43的一次拟合、完整1491窗自然读出及CPU核验已闭合；七项数值证明true，但原五项继续条件全部false，两物种完整正确链均为0，科学FAIL0/5，按原规则关闭模型。详见[完整结果](a43_result.md)。A43检验将完整语法自然选出的结构加入有限已知链监督竞争，而非扩大冻结评分的搜索；真实支持保留原gold并加入1863个已知非匹配、469个未知，9个原自然池外正确几何全部属于原注入gold，实际19189个配对。下文保留准备和执行历史，不把旧“尚未执行”阶段当当前状态；A36与A42仍关闭，不用新实验覆盖旧FAIL。

## 证据与待检验假设

[A42已闭合结果](a42_emission_lattice_train_protocol.md)在原1491个TRAIN窗口中损失1039个正确配对，仅增加11个；1034个丢失正确配对的已选重叠竞争者标签仅为known_nonmatch。全部1778个丢失正分配对通过必要局部替换检查，表明没有找到换回单条旧链即可提高保存目标的具体反例；不能据此证明全空间算法无误或训练目标是因果。

[A36](a36_structured_emission_protocol.md)已经用fresh表示和十五类发射接受有限候选的compatible-known-set似然监督，且完整TRAIN五条件FAIL。因此“新增结构loss”不是新机制。A43的假设限于：向同类训练目标暴露在完整语法推理中实际胜出的已知非匹配几何，可能改善对这类结构的区分。A36另有正确链对null的非正问题，新增竞争者也可能加重拒绝正确链；不预设本机制必然成功。

## 拟议监督支持

原池拟沿用A12的16321个自然候选加536个仅训练参考注入、2916个正配对、12566个已知非匹配和1375个未知，保留原1491窗、1536draw顺序及importance。新自然输入拟使用已闭合A42 R1的所有3630个已选窗口—链配对及原R4保存标签，不按错误标签、676个退步窗口、gain大小或长短挑选；不重新生成候选、不重新配标签。A42这些几何在标签连接前自然搜索已闭合。

每窗以完整exon几何连接，不以候选索引或外层span去重。原顺序、全部原正链和0/−1标签不变；共享几何的标签必须相同，新增几何依A42原自然输出顺序附加。新增标签1若不在原正监督域，停止核对原参考支持，不能默默扩大gold或删除该链。新增0进入已知竞争，新增−1保留未知且不进入loss分区及gold；后续自然推理不能过滤0或−1。新支持大小、每物种新增0/−1及所有新增1是否已在原gold中，尚需真实零更新连接取证，不能按3630直接宣称新增3630个负例。

拟复用原structured_emission_loss，不新设margin、温度、长度归一、角色权重或损失权重。令K为合并后所有0/1候选，G仍为原全部正链，loss为(logZ(K)−sum gain(G))/|K|，logZ覆盖K的全部span兼容集合和null。分母随已知支持改变，原空已知draw若新增0则获得背景监督；这改变每窗监督密度，不能声称相对于A36只改变了负例身份而训练梯度尺度完全匹配。空或全未知draw保留可微零，不删draw、不修改importance。

这是冻结A3来源的离线有限竞争池，既不是当前新模型在线挖掘，也不是全合法语法上的分区函数或完整语法训练。536个参考注入仅用于训练，任何自然预测、继续条件和最终raw-FASTA输出都不得使用它们。

## 拟议模型与继续条件

真实支持及零更新梯度接线闭合后，在用户已直接要求的原TRAIN模型迭代范围准备独立fresh seed0的A36同型366191参数表示和15类发射；只用原缓存GEN与DNA，冻结预训练encoder，不从A36、A41或A42研究checkpoint续训。保留1536draw×3的4608步、原AdamW/importance、末004608、单位FP64角色势、null0及gain0；不试多seed、挑中间模型或调阈值。下方实际准备与一次fit合同明确命令和资源；原模型FAIL和停止条件不重开。

最终目标是完整raw-FASTA注释，不把有限池重评分改善作为本机制成功。拟议TRAIN门槛用同1491窗自然完整语法推理与较强A31_A3原计数比较：两个物种正确严格增加、known_nonmatch均不增加、pooled至少5-exon正确严格增加；同时列A42与A36不同监督/解码支持，不能把额外曝光对不同旧模型的比较当作纯因果消融。完整新自然输出闭合后才连接原R4标签，原2916/1086分母不变，unknown另列。FAIL关闭独立模型，不救参；PASS仍只是下一raw协议与DEV准备依据，不是独立泛化、SOTA或论文创新。

以上是新的监督支持假设，沿用户直接要求继续迭代的原TRAIN范围实施，不以本文件或人工检查存在推定批准，也不为同范围常规下一步另设批准。该范围不包括新增gold、额外数据域、重新访问原始参考、新阈值/曝光搜索、raw/DEV、独立物种、Setaria、A33或公开发布；这些原科学或人类边界保持有效。

## 人工接线与停止条件

src/m28/natural_competitor_pool.py只合并调用方传入的几何与标签，没有数据读取、模型或搜索副作用。五项新检查覆盖完整几何/顺序/标签/gold不变且不改输入；新已知竞争者确实进入原loss及梯度、相冲突unknown无负梯度但自然选择仍允许它；真实支持错配与新增gold被拒绝；原空或全未知draw不丢弃、空已知draw新增0可获得背景监督；保存结果相对原自然池的成员计数与缺失/乱序/共同标签漂移被拒绝。只核对本次新连接，不重复原完整DP或角色算术的全部检查。

唯一命令sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-CHECKS.sbatch。沿用generanno/private/0GPU2CPU8GiB2分钟/no-requeue，独立根outputs/M28-A43-NATURAL-COMPETITOR-CHECKS及M28A43COMPCHECKS_%j日志；只人工tensor、STATUS/JOBID及短日志，不读取真实JSON、DNA、参考、缓存或研究权重，不建optimizer或checkpoint。提交前相关源/AST/SBATCH、五项活调度、节点/预约及精确输出冲突核对。明确实现或fixture错误可同范围同资源另根R1一次；再次失败、TIMEOUT/OOM关闭准备，不延长或盲目重提。PASS只支持候选合并的人工接线，不构成真实支持核验、训练许可或准确率证据。

## 原自然池成员诊断

先回答会改变训练选择的问题：A42已选已知非匹配几何中，有多少不在原16321个自然候选里？若没有新增已知非匹配几何，拒绝将相同竞争支持作为新训练机制；若存在，再核对其是否及如何进入原A12完整监督域。只证明成员差异，不证明训练能改善评分。

census_natural_competitor_support.py沿用A42竞争诊断已经打开的三份JSON：原A31 natural_scores、A42 R1 natural_scores及windows标签行，先核对同源实际终态、七CPU证明、冻结配置/原五条件，再逐窗全身份/ordinal/长度/共享标签检查。按物种和0/1/−1分别统计新已选几何已经在原自然池或不在其中；旧自然全标签须为2380/12566/1375，新已选3630及物种标签计数、1491窗和2916/1086完整分母须与闭合源相符。这里原池是全部自然候选而非旧selected；新增正确几何不在自然池不代表不在训练gold，不能提前断言其违反A43监督域。

五项人工检查闭合后可在同保存结果分析范围运行sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-SUPPORT-CENSUS.sbatch，generanno/private/0GPU2CPU8GiB2分钟/no-requeue；根outputs/M28-A43-NATURAL-COMPETITOR-SUPPORT-CENSUS及M28A43COMPCENSUS_%j唯一日志，只有小summary/STATUS/JOBID。它不开A12注入gold、TRAIN catalog/参考、真实DNA、缓存或研究tensor，不重配标签、搜索、前向或optimizer，不放行新训练。相同提交前核对和一次明确实现/fixture修复上限适用；成员数据错配则停止，不能把矛盾当新机制成功。

## 人工检查的实际闭合

13423961于2026-10-07 03:32:45 UTC提交，03:32:48至03:33:03运行，15秒；root/batch/extern全部COMPLETED0:0，五检查PASS0.515秒，STATUS/JOBID相符。实际private/gpu034/0GPU2CPU8GiB2分钟/generanno/no-requeue/Restarts0，batch MaxRSS769100K；canonical Command/SubmitLine/WorkDir和唯一日志匹配。三个Python AST、两个SBATCH语法及六个新增文件local/canonical逐字一致，提交前五项活调度/节点/预约和精确冲突通过。原M20两处未提交改动保持5增3删及3增3删，无Git提交、研究输入打开、新模型拟合或监控创建；原A42固定监控仍暂停。

这一PASS确认人工新增已知竞争者进入有限集合loss及实际梯度、unknown既不作负监督也不被自然选择过滤，不能证明真实训练支持、完整语法推理或准确率改进。随后只运行上述三份保存JSON的成员诊断，依据实际新增已知几何决定是否继续准备A43，结果如下。

## 原自然池成员的实际结果

13423970于2026-10-07 03:34:09 UTC提交，03:34:12至03:34:26运行，14秒；root/batch/extern均COMPLETED0:0，STATUS/JOBID相符、stderr空。private/gpu034/0GPU2CPU8GiB2分钟/generanno/no-requeue/Restarts0，batch MaxRSS755892K，canonical命令/工作目录/唯一日志匹配。原1491窗/16321自然候选/3630新已选、完整分母、标签计数和冻结A42配置/source_jobs均相符，两项实际保存行证明true。只有三份原允许保存JSON，optimizer0，无额外研究输入或标签生成；未使用R1或延时。

下表按全部原自然候选的完整exon几何区分A42已选成员，不是按旧selected区分，亦不是训练gold成员判断。每格为已在原自然池／新增几何配对。

| 物种 | 已知非匹配 | unknown | 正确 |
| --- | ---: | ---: | ---: |
| 拟南芥 | 108／1351 | 80／327 | 739／9 |
| 水稻 | 106／512 | 23／142 | 233／0 |
| pooled | 214／1863 | 103／469 | 972／9 |

实际新增自然几何共2341（1863个已知非匹配、469个未知、9个正确），共享1289；合计3630。原自然标签完整计数仍2380正确/12566已知非匹配/1375未知，拟南芥1751/6305/853、水稻629/6261/522。不同于此前按旧selected得到1874个新增已知非匹配：其中11个早已在旧自然池但原来未选，本成员诊断没有把它们误称新训练几何。

据此继续准备A43的独立监督支持是有实际对象的，不是重复相同有限池；但1863个新已知非匹配仍不能证明增加监督能够改善模型，A36正确链对null的拒绝风险仍存在。9个新正确几何只确认不在原自然池；须在原A12既有536个训练参考注入中逐几何核对其资格，不能新增正例救援。若均属于原gold，拟合并支持为16857+1863+469＝19189；这是条件推导，尚不是已验证真实监督池大小。unknown继续保护，新增gold、当前模型在线挖掘、全语法训练和首训权限都不能由本成员统计推定。

## 完整保存监督支持的连接

依据用户直接要求继续迭代，在原TRAIN监督和保存结果的分析范围，下一步读取原A12 windows保存监督记录及原A36训练第一遍1536条draw收据，连接上述三份自然结果JSON和原draw记录。原A12的536个注入正链是已保存训练几何，不是重新访问原始参考文件或新增正例。本步不开TRAIN catalog、原始参考、DNA、缓存、checkpoint或DEV；不重新标签、搜索、模型前向或更新。记录文件的存在不充当授权依据，也不为既有授权的常规连接新设批准要求。

build_natural_competitor_support.py逐窗核对原R4保存窗口的十字段几何投影、TRAIN六条允许染色体、首次draw顺序、旧自然前缀和0/1/−1、原候选来源和计数；完整原gold必须span兼容且对应既有2916/1086分母。以完整exon几何合并全部A42自然输出，9个原自然池外正确几何须逐一属于原reference_injected正链；不加入新正链、不丢弃矛盾链。first-pass收据须逐draw与原species/顺序/importance一致，原replay计数还须与A36实际CPU证明和A12summary相同。任何真实来源、标签或gold矛盾停止本机制准备，不能按预期数修数据。

只在全部连接成立后，独立根outputs/M28-A43-NATURAL-COMPETITOR-SUPPORT以x模式写1491条最小windows监督记录、1536条未经修改draw副本及summary。候选只保留几何/标签/原来源或新增自然索引，不将旧模型数值gain或参考上下文传给未来模型；报告每物种新增0/−1、9个注入匹配、原gold数及原/合并draw计数。19189仍是待验证预期，实际所有输入一致后才能宣称该有限支持已构建。

先运行sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-SUPPORT-CHECKS.sbatch：三个新增人工adapter检查覆盖注入gold的自然几何连接与重复draw统计、来源/标签/分母/importance矛盾拒绝，以及原空已知draw获得新背景监督但不被删除。只测本次新记录连接，不重复已关闭的五项loss/成员检查。闭合后运行sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-SUPPORT.sbatch ACTUAL_CHECKS_JOB。两者generanno/private/0GPU2CPU8GiB2分钟/no-requeue，独立根SUPPORT-CHECKS和SUPPORT，唯一M28A43SUPCHECKS_%j/M28A43SUPPORT_%j日志；提交前核对新源/AST/SBATCH、五项活调度、节点/预约及精确根/日志冲突。明确实现或fixture错误允许同范围同资源另根R1一次；数据矛盾、第二次失败、TIMEOUT/OOM停止，不延时或修改科学支持。闭合仅表明有限监督artifact可用，不是梯度预检、新拟合、模型改善、raw注释或泛化成立。

## 完整监督支持的实际闭合

13424084于03:52:45 UTC提交，03:52:47至03:53:00运行13秒；三个adapter检查PASS0.022秒。13424089于03:54:11 UTC提交，03:54:12至03:54:34运行22秒，四项实际完整支持/gold/分母/draw证明true。两作业root/batch/extern全部COMPLETED0:0，STATUS/JOBID、canonical Command/SubmitLine/WorkDir及唯一日志相符；均private/gpu034/0GPU2CPU8GiB2分钟/generanno/no-requeue/Restarts0，batch MaxRSS分别750816K、813480K。支持构建stderr空；1491条windows及1536条draws文件已闭合，五个本次文件local/canonical逐字相同，两AST/两SBATCH语法及提交前活调度/节点/预约/精确冲突检查通过。没有R1、延时、Git提交或监控重启。

9个原自然池外正确几何全部在拟南芥原reference_injected监督中匹配，水稻0；原536个注入及2916个正配对全部保留，不加入新gold。有限监督池实际19189＝16857+1863+469，标签为2916正/14429已知非匹配/1844未知。拟南芥2172/7656/1180，水稻744/6773/664。原完整2916/1086分母、原顺序和全部标签保持一致；候选记录没有携带旧数值gain或原参考上下文。

原1536draw replay仍16806自然/555注入、3028正/12935已知非匹配/1398未知，与A12和A36实际记录一致。合并后replay新增2417个自然配对（1931个0、486个−1），标签为3028/14866/1884；zero_known_windows仍3，并非人工空池例子的背景监督变化在真实数据中已经发生。唯一窗口与重复draw配对数分开报告，不称独立gene或新增生物负例。没有原始参考/catalog、DNA/feature/checkpoint读取、标签生成、前向或optimizer。下一步只检验实际新池通过fresh同型表示的零更新梯度，不把artifact闭合当训练收益。

## 实际新支持的零更新梯度预检

依据用户要求持续迭代，已关闭支持允许在原TRAIN缓存输入上执行必要梯度接线，不改变既有数据域。preflight_natural_competitor.py读取13424089的完整新监督artifact、原draw和feature-plan/summary，只打开八个预定证据窗口对应的原缓存GEN+DNA，fresh seed0同型366191参数表示/15类发射；不加载任何研究checkpoint、不建optimizer、不写checkpoint、不更新预训练encoder或重新配标签，不开原始参考/catalog、DEV/独立物种/Setaria/A33。

每物种依首次ordinal稳定选择四个不同的、有正和已知非匹配的窗口：最多新增0、最多新增−1、最多正链exon、最多正负span冲突；前两类必须确有新增对象。这只是梯度失败证据选择，训练仍保留全部1491窗/1536draw，不能据此挑选训练窗。实际role势经原structured_emission_loss和原importance回传，五参数模块梯度须有限非零、初始权重最终逐tensor不变。以同实际势独立leaf接原chain_set_loss，loss精确相同，新增已知0的gain导数不得为负，所有−1导数必须为零；记录正导数数量和最大值，数值零不能误说每个新负例都已产生非零梯度。该导数不是自然预测标签过滤器，也不证明训练后的抑制效果。

命令sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-PREFLIGHT.sbatch 13424089，generanno/private/1×24GiB RTX3090/2CPU/16GiB/10分钟/no-requeue，独立根outputs/M28-A43-NATURAL-COMPETITOR-PREFLIGHT、唯一M28A43PREFLIGHT_%j日志。只有八条梯度receipt和小summary/STATUS/JOBID，提交前新AST/SBATCH语法、实际源、五项活调度/节点/预约与精确输出冲突检查；不重复原已关闭的人工数学检查。真实数据/梯度错配、非有限值或权重改变停止准备；明确实现错配允许同范围同资源另根R1一次，第二次失败或TIMEOUT/OOM停止，不改变loss、阈值、支持、seed或延时。PASS仅作为下一独立fit的工程就绪证据，不是TRAIN五条件、raw模型或发表证据。

### 预检实际闭合

13424115于03:58:48 UTC提交，03:58:51至03:59:21运行30秒，root/batch/extern均COMPLETED0:0，STATUS/JOBID/summary、实际配置/支持13424089/八窗口及Command/SubmitLine/WorkDir/唯一日志相符、stderr空。private/gpu034/1RTX3090/2CPU16GiB10分钟/generanno/no-requeue/Restarts0，batch MaxRSS1832668K。内部13.506秒、第一窗9.523秒，其余七窗0.432至0.854秒；峰GPU allocation163786240字节。这不是全训练计时或速度benchmark。提交前新AST/SBATCH及三文件local/canonical一致，活调度/节点/预约/精确冲突通过。

八窗五模块实际梯度均有限非零，366191参数的state前后相同，fresh seed0/optimizer0/无研究checkpoint读取。新监督与原feature几何、same-potential loss连接相符；22个unknown配对的gain导数全为0。41个新增已知非匹配证据配对中29个FP64 gain导数严格正、12个数值零，部分正导数极小（如7.30e−182），进一步回传FP32网络时可能下溢。因此只能说同实际监督池确实接入有限loss且两个最多新增0窗口存在约0.02的抑制gain导数，不能说所有新竞争者初始都产生有效非零参数梯度，更不能说它们训练后已被抑制。该数值情况保留，不改温度/score/loss救梯度。

## 一次独立拟合合同

原TRAIN新表示拟合沿用户的直接模型迭代授权推进，不要求对已授权的常规下一步再批准；新支持及零更新预检仍须实际闭合。fit_natural_competitor.py只接相同配置/同支持job的真实闭合预检，fresh seed0同型366191参数，仅19189有限保存监督、原1536draw顺序/importance×3＝4608 AdamW更新、原超参数/null0/阈值0。训练范围仍六条原TRAIN染色体，不访问原始catalog/参考、DEV或独立物种；不继承研究checkpoint、不运行GEN前向、不新搜索/配标签。

仅独立根outputs/M28-A43-NATURAL-COMPETITOR-FIT写training.jsonl、001536/003072/004608三个checkpoint及summary，主科学读出只用末004608。每遍真实counts须与已构建支持的原draw replay一致，非有限loss/gradient/state即失败；保存有限FP32模型与完整AdamW状态用于后续实际步数/配置核验，不能据训练loss挑epoch。训练新增监督身份和known分母同时改变，新增竞争池又来自已训练A3，故不是总曝光/梯度尺度匹配的单因素因果消融，不称全语法训练或创新已成立。

命令sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-FIT.sbatch 13424089 ACTUAL_PREFLIGHT_JOB，generanno/private/1×24GiB RTX3090/2CPU/16GiB/2:30/no-requeue，唯一M28A43FIT_%j日志。时限保留同型原A36已经实际1:15:07完成的一次fit上限；新有限支持约增加14%，八个梯度证据窗口不是随机计时总体，不能以其热启动速度承诺全拟合时长。小训练日志/三个checkpoint预计低于64MiB，无新feature/logits缓存；提交前新AST/SBATCH/实际源、实际预检、存储、五项活调度/节点/预约与精确输出冲突核对。只有确实尚无optimizer更新的明确实现故障允许一次同科学/同资源另根R1；已有更新、第二次失败、TIMEOUT/OOM不自动重提、续训、延时或改参，保留原失败证据。后续完整语法自然TRAIN读出及实际数值核验仍未完成；提交/训练终态都不替代原五条件或论文目标。

### 实际拟合启动与固定监控

13424143于2026-10-07 04:04:44 UTC提交，04:04:47开始gpu034/private，唯一命令sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-FIT.sbatch 13424089 13424115。启动实际root/batch/extern RUNNING0:0，STATUS RUNNING/JOBID13424143，1RTX3090/2CPU16GiB/2:30/generanno/no-requeue/Restarts0、Command/SubmitLine/WorkDir/唯一日志相符；04:05:42启动证据日志为step64/4608、内部34.091秒、stderr空。新fit AST/SBATCH、三文件local/canonical一致、当前存储156T可用以及提交前五项活调度/节点/预约/精确冲突检查通过。只提交一次，没有续训/延时、Git提交或修改原研究结果。RUNNING及记录64步不是拟合闭合、科学五条件或准确率证据。

只复用固定轻量对话01a0f71f-8845-7432-bd13-bf54c1808fb3和同heartbeat coding-rna-r5-light-monitor，将旧A42已消费对象替换为唯一fit13424143；native实际ACTIVE/fullprompt/fixedtarget/15分钟/原通知字段逐项读回一致，未新建对话或cron。当前同对话manual turn01a1148b-8a3e-7bc2-a453-7a199de8c3f0完成/errornull/43.970秒，04:07:52UTC最小snapshot为RUNNING/gpu034；04:07:31.682Z实际turn_context核验gpt-5.6-luna/low/同cwd。配置和一次manual执行已验证，本轮scheduled运行与完整终态回传尚未验证，不能沿用旧A42的scheduled证据冒充本轮。监控只在当前fit实际终态/可行动故障时一次回传并暂停；不读loss/权重/性能、不改科研选择或提交下游。主不形成重复等待轮询，独立准备后续完整语法自然TRAIN读出；该读出和实际数值核验尚未提交或完成。应用须运行且连接可用，休眠或断线时不能保证即时回传。

## 末模型完整语法自然 TRAIN 读出合同

readout_natural_competitor_train.py只用实际闭合fit13424143的A43末004608，不用早期checkpoint或旧A36/A42权重。预测输入从原1536draw和feature-plan提取首次出现的1491个十字段窗口几何及原缓存GEN+DNA；GPU发射和CPU搜索不打开训练候选池、注入正链、R4标签或catalog。输入仍是原TRAIN窗口，不是未知组装的raw-FASTA独立泛化测试；本实现将自然推理的接口和监督池分开，但不据此宣称已经完成论文主模型。

先在Slurm运行六个新增人工adapter检查：原预测几何/顺序不接gold、实际A43末模型身份、无标签评分及null保留、新几何R4与unknown/完整分母、原五条件严格不等号、真实训练收据的原权重/union已知分母/三遍顺序及空已知draw。随后只连接原保存metadata，核对全部1491窗口/1536draw/原缓存artifact、16321自然池、较强A31_A3计数及2916/1086分母；此准备不打开真实feature tensor、研究模型或TRAIN catalog。旧完整语法DP、人工数学检查和八窗计时不重复运行，算法与实际A42闭合版本相同。

四阶段依次提交：TRAIN-POTENTIALS在fit实际闭合后冻结366191参数，FP32前向保存1491份有效碱基15类发射和DNA；TRAIN-DECODE按原ordinal modulo32进行32个完整语法搜索；TRAIN-READOUT先证明全部32分片终态与预测文件闭合，恢复原顺序写完自然文件后才加载原R4元数据；RESULT-CHECK在CPU核验实际训练与全部保存预测。角色单位权重、FP64累加、null0/gain0、max128、原ORF/motif/phase/intron/span规则和稳定tie不变，无query/beam/在线监督筛选、GEN前向、优化器更新、新参数、DEV/独立物种/Setaria/A33。

RESULT-CHECK逐步核对4608条真实训练收据和原1536draw×3，重算已知union NLL/importance及每遍counts/mean，读取三份实际有限FP32模型和AdamW状态核验配置/参数/步数。新自然geometry逐条核对DNA ORF/motif/cap/span、有效碱基FP32势和独立whole-window FP64 log-softmax角色差；旧自然候选只在搜索全部闭合后按同一A43发射重新评分、做null0 span DP，完整语法目标不得低于该有限池目标。这是数值下界，不是旧A31准确率基线，也不是独立模型前向或第二次完整DP最优性证明。最后重新连接未改R4的0/1/−1、完整2916/1086分母和较强A31_A3原计数，重算原五条件与按物种/unknown计数。七项实际证明全部闭合后才解释科学结果，不拿训练loss、局部链或下界PASS冒充性能改善。

SBATCH命令均用独立M28-A43-NATURAL-COMPETITOR-{READOUT-CHECKS,TRAIN-POTENTIALS,TRAIN-DECODE,TRAIN-READOUT,RESULT-CHECK}脚本与同名输出根，generanno/private/no-requeue、唯一M28A43READCHECKS/POT/DECODE/READ/CHECK日志。CHECKS为0GPU2CPU8GiB2分钟；POTENTIALS为1RTX3090/2CPU16GiB30分钟，afterok实际fit并要求实际closed CHECKS；DECODE为0GPU2CPU8GiB2:15、array0–31%32，afterok全POTENTIALS；READOUT和RESULT-CHECK均0GPU2CPU8GiB15分钟，依次afterok全数组/READOUT。传入实际CHECKS_JOB/POTENTIAL_JOB/DECODE_JOB/READOUT_JOB，不猜ID。预测发射预计约2.2GiB；资源沿用已实际闭合的同型1491窗发射及32分片语法流程上限，不据此承诺新fit或搜索速度。

每次提交前核对相关源/AST/bash语法/canonical内容、五项活调度/预约/节点及精确目标根/日志和实际存储；不重复旧源审查或主线程等待轮询。明确新adapter实现错误可同科学/同资源另根R1一次；来源/标签/数值矛盾、第二次失败、TIMEOUT/OOM停止，不改科学协议、重训、救参、延时或覆盖结果。当前fit已更新权重，任何其失败不自动重提。原A36/A42 FAIL保持关闭；A43原五条件FAIL关闭该模型，PASS仅支持后续raw/DEV独立协议准备，仍不自动放行新域、公开发布或论文主张。整条下游实际提交后只扩展原固定监控和同heartbeat对象，终态/可行动失败一次回传并暂停；当前新准备与下游尚未执行，不能写为闭合。

### 新读出准备的实际闭合与下游提交

13424381于2026-10-07 04:33:48至04:34:15 UTC运行27秒，root/batch/extern全部COMPLETED0:0；六个新adapter检查PASS0.015秒，原1491窗口/1536draw/feature artifact/16321自然池/较强基线计数/2916及1086完整分母实际连接成立。STATUS/JOBID/summary匹配，generanno/private/gpu034/0GPU2CPU8GiB2分钟/no-requeue/Restarts0，batch MaxRSS777224K，实际Command/SubmitLine/WorkDir及唯一日志一致。只人工tensor和原保存metadata，无真实feature tensor/研究权重/catalog/新labels或forward；stderr中的六项unittest均OK，非空测试日志不误称空stderr。三AST/五bash语法通过，本次九文件local/canonical逐字一致；实际存储156T可用，提交前各项活调度/节点/预约及根/日志冲突核对通过，没有R1或Git提交。两次scontrol只读查询的参数拼写错误已按实际CLI改正，未改变或重提作业，不影响实际终态证据。

随后四阶段分别于04:35:36/04:36:02/04:36:30/04:36:51 UTC各提交一次，实际为potentials13424395 afterokfit13424143/exportCHECKS_JOB13424381；decode13424400 afterok13424395/exportPOTENTIAL_JOB13424395，array0–31%32；readout13424410 afterok完整13424400_*/export同POTENTIAL_JOB及DECODE_JOB13424400；CPU13424413 afterok13424410/export再加READOUT_JOB13424410。实际Command/SubmitLine/WorkDir、原资源/时限/excludes/no-requeue/Restarts0、唯一日志及全部数组依赖逐项核验，提交前每次均完成五项fresh调度/节点/预约/存储与精确冲突检查。04:37:11 UTC启动snapshot为fit RUNNING/gpu034、四下游PENDING0:0普通Dependency；没有研究终态或准确率结果，不把排队当失败。

原固定heartbeat和同固定监控对话扩展为M28-A43-NATURAL-COMPETITOR-FULL-TRAIN-20261007，36实际root纳入同完整流程：fit1/potentials1/decode32/readout1/CPU1，不在fit单独闭合时提早交接。native完整prompt/ACTIVE/kindheartbeat/固定target/15分钟/原failed_runs_only/created_at1790981930511读回全部相符，updated_at1791347967261；CPU白名单false字段与实际producer逐名一致，manual调用前已改正两处监控字段别名，未更改科学代码或协议。现仅结构与七实际CPU证明可被监控读取，模型/行内容/labels/counts/gates由主处理。先前单fit的scheduled turn01a11499-86f2-7400-8a91-5e3c8b10567d于04:22:48.364Z确为gpt-5.6-luna/low，实际heartbeat trigger04:22:48.553Z，04:23:02Z snapshot RUNNING，completed/errornull27.043秒/cursor125；该定时证据属于原单fit，不冒充扩展cohort的运行验证。新完整cohort manual检查已发给同对话，实际执行和定时终态回传仍须后续证据。

新完整cohort manual turn01a114a9-5cbc-7211-ac12-bf6c7208551a已完成/errornull35.137秒/cursor128；04:40:06.119Z实际turn_context保持gpt-5.6-luna/low/原cwd，04:40:22Z最小snapshot为fit RUNNING/gpu034/0:0、四下游普通Dependency/PENDING0:0。同监控已追加本次证据，没有读取研究产物、向主报健康进度或改变native配置；主核验完成记录，不建立并行作业等待轮询。新完整cohort配置与一次manual检查已验证，scheduled执行及全流程自动终态交接仍未验证；不存在完整预测、数值证明或科学五条件结果。论文目标仍未达到，下一步是由同固定监控在实际终态/可行动故障交接，再由主验证并解释保存结果。

## 单模型 raw FASTA 接口准备

infer_natural_competitor.py将A43末004608的366191参数表示接到显式FASTA原生编码与完整contig语法，不加载A3 query生成器、不依赖训练候选池、注入链或参考标签。GEN仍为原冻结六mer/四块编码，保留缓存原BF16 materialization后FP32发射的实际实现；本次人工检查使用替代编码器，不能称为原生GEN数值复现。每个原WIDTH24576/STRIDE12288方向窗只调用一次编码及单表示，传入实际有效DNA长度，十五列FP32势交由已有contig helper，不调用原HMM、greedy/beam或旧selected-only owner装配。

窗口最优gene set不能直接送入旧owner装配：A42人工反例已经证明局部选择后过滤会丢链，owner-before-local-DP也不能保留全部跨窗竞争。新接口因此复用已闭合的逐碱基上下文定义：每个genomic base选包含它且中心最近的窗口，平分取较小genomic start；完整复制该窗口的十五列，不平均、不重排负链角色列、不含padding。随后整条oriented contig一次联合语法DP，最后把完整chain及其gain成对转genomic坐标和GFF3 phase。允许跨24576bp接缝的完整链；这是不同于TRAIN独立窗口及旧whole-chain owner的上下文、评分和span域，不把接口等价性或人工长链当作真实模型收益。

原contig_emission_lattice.py及原emission_lattice.py保持不变，单位角色/FP64/null0/gain0/max128/完整ORF与splice/carry规则不变；不引入未经验证的快解码替代或窗口裁切。已闭合的49152bp密集人工规模仅证明两个人工输入可执行，不能据此外推染色体成本。真实raw科学/资源合同a43_natural_competitor_raw_protocol.md仍不创建，亦无GPU执行SBATCH或新数据运行；实际TRAIN五条件和七CPU证明、较强A31_A3原逐物种计数与2172/854及744/232完整分母都必须真实成立。入口在模型/GEN/FASTA和输出创建前核对实际36根与配对source_jobs，重算五条件，不相信孤立PASS字段；还要求另行固定科学/资源合同与显式已授权FASTA。文件存在是协议前提，不是人类许可证据；不能默认读取DEV、独立物种、Setaria或A33。

七项新增检查只覆盖当前适配层：源身份/七CPU证明/未重基线及完整物种分母；失败五条件或缺合同/显式FASTA时禁止原生模型、输入和输出；人工A43末状态只载入一个冻结表示、早期状态拒绝且旧生成器禁用；规范化/双链单表示路由/负链CDS phase/GFF3；30000bp人工完整链通过原两窗发射接到一次联合选择；零/负gain、无ORF及空contig；人工空或重复FASTA流不得写成功summary/GFF。长链稀疏人工序列与小随机模块不是染色体benchmark或训练模型预测，旧DP全套oracle/规模检查不重复。

唯一短检查命令sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-RAW-CHECKS.sbatch，generanno/private/0GPU2CPU8GiB2分钟/no-requeue，独立根outputs/M28-A43-NATURAL-COMPETITOR-RAW-CHECKS及M28A43RAWCHECKS_%j日志。只人工DNA/tensor/替代编码与临时文件、STATUS/JOBID和短日志，不读取研究权重/GEN/FASTA/缓存/参考/标签或DEV，不建optimizer；受控fixture的模拟合同不创建真实协议。本次新脚本/检查及本节明确传输canonical，提交前AST/bash/canonical逐字、五项fresh调度/节点/预约与精确根/日志冲突核对。明确adapter/fixture错误可同科学/同资源另根R1一次；第二次失败、TIMEOUT/OOM或数值/语法矛盾停止，不延长、改cap或修改运行中的TRAIN流程。该短准备由主处理，不纳入原唯一完整流程监控，不新建监控或改其模型/对象。

### raw 接口检查实际闭合

13424499于2026-10-07 04:54:23 UTC提交，04:54:24至04:54:39运行15秒；root/batch/extern全部COMPLETED0:0，STATUS/JOBID13424499匹配，七个新检查PASS0.502秒。实际private/gpu034/0GPU2CPU8GiB2分钟/generanno/no-requeue/Restarts0，batch MaxRSS585000K，Command/SubmitLine/WorkDir及唯一日志相符。两份Python AST及bash语法、本次四文件和两个未改语法helper的local/canonical逐字核对通过，五项fresh调度/预约/节点及根/日志冲突检查通过；没有R1、资源延长、Git提交、真实FASTA/模型/GEN/参考/缓存或新域读取。stdout为空、stderr含七项成功unittest，保留原日志不误说stderr为空。

人工A43末状态仅构建一个冻结366191参数module，早期状态被拒绝且旧生成器入口没有调用；这是fake权重/接口路径检查，不是研究checkpoint数值证据。人工负链完整CDS导出为[(0,5,phase2),(10,14,phase0)]，GFF坐标分别1–5和11–14，相位保持2/0，joint gain140。30000bp人工链[3000,33000)跨原[0,24576)与[12288,36864)两窗，任何单窗都不能覆盖它，新发射路由的joint输出为完整链、gain300000；每方向每窗一次替代编码和单表示，共四次，未使用局部selected-only装配。零/负势、无ORF及空/重复流均按原约定返回null或禁止成功标记。稀疏人工长链和手工势不构成真实长基因发现、密集染色体速度或独立最优性新证明；原DP及上下文helper未改变。

consumer中的较强基线和分母还逐字段与已闭合A42 CPU保存metadata核对：拟南芥A31_A3 selected2193/正确1531/已知非匹配347/unknown315/至少5exon正确548、eligible2172/long854；水稻861/478/251/132/152、eligible744/long232，完全一致。此次只复核consumer常量，未重新评价旧A42或修改其FAIL；实际A43五条件仍等待完整流程，人工PASS记录不放行raw执行。

同固定监控的完整cohort scheduled turn01a114b4-fec2-7821-ac06-ea642a93ab37已实际完成/errornull24.361秒/cursor129；04:52:48.550Z turn_context保持gpt-5.6-luna/low/同cwd，04:52:48.710Z heartbeat user trigger含current_time_iso04:52:48.449Z。04:52:58Z实际最小snapshot为fit13424143 COMPLETED0:0、potentials13424395 RUNNING/gpu034/0:0、其余三阶段普通Dependency/PENDING0:0，未读取拟合/发射研究产物或提早交接。native配置仍ACTIVE/同固定target/15分钟/原通知字段/updated_at1791347967261，未更新自动化或增加短检查目标；配置、一次manual和本cohort一次scheduled都已验证，全流程自动终态回传仍未发生。训练结束不是数值曝光、自然性能或五条件闭合，后续仍由同唯一监控等待36根完整终态，再由主解释结果。

## 整条 contig 的完整评价接口准备

evaluate_natural_competitor.py仅提供人工可测的适配函数，不增加真实raw或DEV执行入口。显式FASTA的实际producer将species写为user_input；评价必须另行获得固定的contig到物种/长度映射，不能凭名称自动认作原DEV。入口先限制到原完整评价scope，重建每方向整条contig的已选链，核对原网格/逐碱基上下文、oriented与genomic转换、成对gain顺序、正gain/null0总和、max128/最短intron5/同链与链间非重叠、所有预测数量及GFF3坐标和phase。关闭保存预测和GFF后才允许读取合同内显式FASTA，再进入参考；与完整参考上下文逐碱基比较规范化DNA，预测逐条通过原ORF/phase/splice规则，不以过滤错误预测修好指标。

评价复用原完整5437/2291、合计7728条primary CDS chain分母和全注释gene背景，不按训练支持、窗口宽度或预测可达性删参考；零预测物种仍保留。与A31_A3、A31_A23及ANNEVO/Tiberius/Helixer共享原完整参考和背景尺，报告逐物种、pooled、macro与worst-species F1，保留原生基线的曝光限制。联合contig没有旧free/owner候选阶段，不制造候选覆盖率或同候选池因果主张；nonexact-primary也不直接称生物学假阳性。当前JSON未保存发射tensor，gain求和与链坐标核对不是独立神经前向、发射得分或完整DP重放；真实执行合同仍须保证显式FASTA在推理与评价之间不变。

九项新检查只使用人工JSON/GFF/DNA及模拟收据，含负链多gene顺序与phase、跨两原窗的30000bp完整链、空预测scope、相同名称长度但DNA错配、预测语法违规、参考读取前的保存结果/输入完整性，以及7728条人工参考和人工原生比较表的完整分母。人工参考包含超过窗口宽度的记录，用于检测意外删真值，不是实际DEV/原生性能。唯一短命令sbatch --parsable sbatch/M28-A43-NATURAL-COMPETITOR-EVAL-CHECKS.sbatch，generanno/private/0GPU2CPU8GiB2分钟/no-requeue，独立outputs/M28-A43-NATURAL-COMPETITOR-EVAL-CHECKS和M28A43EVALCHECKS_%j日志。只执行新增adapter检查；不打开研究checkpoint/GEN/缓存/真实FASTA/参考、TRAIN行内容、DEV、新域或旧科学结果，不建optimizer，不改运行中的TRAIN流程。提交前核对新AST/bash、本次明确文件与canonical一致、五项fresh调度/节点/预约及精确根/日志冲突；只明确实现/fixture错误允许同科学同资源另根R1一次，数值/来源矛盾、第二次失败、TIMEOUT/OOM停止。该短准备由主处理并消费，不扩展原唯一完整流程监控；真实raw合同和执行仍未放行。

### 完整评价接口检查实际闭合

13424655于2026-10-07 05:11:45 UTC提交并开始，05:12:02结束17秒；实际root/batch/extern均COMPLETED0:0，STATUS/JOBID匹配，九项新增人工检查PASS0.477秒。private/gpu034/0GPU2CPU8GiB2分钟/generanno/no-requeue/Restarts0，batch MaxRSS755796K，实际Command/SubmitLine/WorkDir及唯一日志相符。两份新Python AST及bash、本次四文件local/canonical逐字核对通过，05:11:19 fresh调度/节点/预约和精确根/日志冲突检查通过；没有R1、延时、研究数据读取、实际GEN前向、新模型更新或真实DEV/原生评价。人工完整尺保留全部7728条、超窗口参考、零预测物种及背景外新增链，错误GFF phase/坐标、DNA错配和不合法预测都被拒绝；这些是fixture验证，不是模型性能数字。

stdout为0字节，stderr3753字节含九项unittest成功和既有scripts/eval_structure_diagnostic.py:51的ResourceWarning。实际parser的enumerate(open(path))未显式关闭文件；本轮未扩大修改到公共parser，也未抑制警告，原始日志保留。接口检查已消费并从等待监控排除。提交前主只读调度可见32个A43 decode分片正常RUNNING，完整读出与CPU仍等待依赖；固定监控另有05:07:58 scheduled最小snapshot，fit/potentials已COMPLETED0:0、32 decode RUNNING、readout/CPU普通Dependency，尚无全流程终态交接或科学五条件结论。下一步仍由同固定监控等待真实完整闭合，再由主处理原TRAIN数值和五条件；本次评价准备不放行raw执行。

## CPU 编译版精确递推准备

实际同人工49kb剖析与缓存回放表明，数组分块没有消除大量Python状态扩展及对象构造；缓存版仍需61.483640/71.883514秒。新的emission_lattice_compiled.py仅在独立人工范围内使用已安装generanno的Numba0.62.1，不安装依赖、不接线原DP/contig/raw或运行中的1491窗A43。保留原标量索引递推、全actual carry、exon-count和所有合法出口，21种carry按原Python字符串lexical顺序编码，前驱节点替代整tuple复制。前驱arena保留曾被保留的节点至单次解码结束，染色体规模的内存成本尚须实测，不能承诺当前实现已经适用全染色体。

角色数组保持原Torch FP64 margin/cumsum/差分值，Numba只编译元数据、合法性、逐出口原顺序加法、严格胜出比较和状态/回溯；fastmath关闭，不做相位合并、预选入口、beam/top-k、cap/阈值/角色权重变化或全零快捷返回。成熟intron减prefix与acceptor加prefix仍分步进行。只有必要ORF终点范围排除必定非法事件；必须与原全部链/顺序/每个gain/总gain及八项计数逐值相等，实际nopython signature存在，不允许object-mode回退。热调用耗时与第一次编译加调用耗时分开记录，不能用既往不同负载单次墙时作公平提速结论。

八项新增编译实现检查复用已支持的人工用法：独立小图DNA穷举；全部相位/actual carry、GT/GC/内含子N、1–2bp短exon、负前缀、邻接gene与null；768bp整数/非整数FP32；128-exon与cap127；原两方向跨窗30000bp完整exon；实际FP64首次到达反例；原Torch margin数组逐值相同；None长相位下合法1bp短carry不得消失。只运行新实现与原标量对照，不重复未改的所有优化版本。唯一命令sbatch --parsable sbatch/M28-A43-COMPILED-LATTICE-CHECKS.sbatch，generanno/private/0GPU2CPU8GiB5分钟/no-requeue，独立outputs/M28-A43-COMPILED-LATTICE-CHECKS和M28A43JITCHECKS_%j日志。5分钟包括首次LLVM编译和原人工对照，不保证完成，不自动延长；提交前新AST/bash、指定文件canonical逐字、五项fresh调度/节点/预约及精确根/日志冲突核对。仅明确编译接口/fixture实现错误允许同语法同资源另根R1一次，链/分数/计数矛盾、第二次失败、TIMEOUT/OOM停止；通过后才回放未挑选的原两项49kb配置，沿用原12分钟8GiB。无模型/GEN/训练、研究输入/标签/真实FASTA/DEV/新域或公开发布；该短准备由主消费，原36根唯一监控不变。

首轮13424814于2026-10-07 05:28:27 UTC提交、05:28:28开始、05:28:56结束28秒；root/batch FAILED1:0、extern COMPLETED0:0，STATUS FAILED exit=1/JOBID匹配。实际private/gpu034/0GPU2CPU8GiB5分钟/no-requeue/Restarts0，batch MaxRSS735072K，Command/SubmitLine/WorkDir和日志相符。setUpClass第一次LLVM lowering失败，Ran0tests，没有kernel自然结果或任何研究输入/模型更新，不能记录八项PASS。Numba0.62.1将typed-dict pop返回值推断为OptionalType(Tuple(float64,int64))，即使前面已有membership判断，static_getitem也不能lower；完整stderr保留具体异常，不是ORF/分数/tie不一致。

唯一R1只把已确认存在的completed/pending键从pop改为先getitem、再del，避免Optional类型；Python语义和删除顺序相同，不改递推、评分、语法、cap、FP64/tie、人工输入或检查。命令sbatch --parsable sbatch/M28-A43-COMPILED-LATTICE-CHECKS-R1.sbatch，资源和5分钟上限不变，独立outputs/M28-A43-COMPILED-LATTICE-CHECKS-R1及M28A43JITCHECKSR1_%j日志。提交前重新核对改动源AST/bash/canonical和五项fresh调度/节点/预约及冲突。再出现编译/答案/分数/计数问题或TIMEOUT/OOM即关闭该计算原型，不做R2、延时或接线研究流程。

### 编译原型 R1 实际结果与停止

13424849于2026-10-07 05:33:11 UTC提交并开始，05:33:36结束25秒；root/batch FAILED1:0、extern COMPLETED0:0，STATUS FAILED exit=1/JOBID匹配。实际private/gpu034/0GPU2CPU8GiB5分钟/generanno/no-requeue/Restarts0，batch MaxRSS802708K，Command/SubmitLine/WorkDir及唯一日志相符。首次编译加人工null调用14.797153秒，Numba0.62.1实际nopython signature为1、fastmath=false；不打开研究输入、optimizer0。stdout232字节、stderr1806字节及原失败根完整保留，Ran8tests17.518秒，七项通过、一项ERROR，不是八项PASS。

实际通过项为小图独立DNA穷举、短actual carry及GT/GC/N/负前缀/邻接gene、128-exon及cap127、两个方向30000bp完整exon、FP64首次到达反例、原Torch FP64数组相等，以及None长相位下合法短carry。dense随机矩阵的首个768bp整数/cap2对照已完成整个结果字典（含原八项计数）相等，随后继承的batch版测试读取其专属state_update_proposals字段，触发KeyError；其余整数/非整数、cap4/8及重复未执行，不能报告完整随机矩阵通过。该错误没有显示链或数值不一致，也不能因此推定尚未完成的对照正确。

事后只在新compiled测试类覆盖该方法，保留seed71、原12次序列及caps2/4/8、整数和非整数FP32输入、原逐值完整结果比较，去掉不适用于标量编译实现的batch计数减少断言；不修改旧batch测试、不给编译器伪造计数字段、不改任何递推代码。此fixture修正未重新运行，原R1仍FAILED，完整验证未通过。按本轮一次R1停止条件，编译原型的执行分支关闭：不提交R2、不运行49kb回放、不延时或接线原1491窗/contig/raw流程。七项人工通过不构成速度、染色体内存、模型准确率或泛化证据。

该短准备由主消费并排除监控，不是A43真实TRAIN科学FAIL；A36/A42原FAIL仍关闭。固定监控05:37:58 UTC的最新最小snapshot仍为fit/potentials COMPLETED0:0、32个decode RUNNING、readout/CPU普通Dependency，完整五条件尚未闭合。原唯一36根监控、模型、资源、参数、原数据与科学停止条件不变；下一主科学动作仍是消费实际完整终态交接，不由编译准备放行raw/DEV或新域。

## 原局部监督与有限整链监督的接线差异

2026-10-07。A3／A23的十五类局部发射有直接的局部监督，A43没有该项；但“缺少局部监督导致新模型较差”仍是待检验假设，不是代码比较能够证明的因果。三个保存fit metadata均为STATUS COMPLETED、completed=true、记录optimizer_steps4608，参数数分别649457／649617／366191。A3／A23保存loss_coefficients均为[1,1,1,1]；A43保存local_cross_entropy=false且loss为exact_compatible_known_candidate_set_NLL_divided_by_known_count。前两者没有A43专属的local_cross_entropy字段，不能把字段缺失解释成未使用局部监督。

实际固定入口为M28-A3-FIRST-FIT.sbatch→run_query_known_mean.py::supervised_readout→query_known_mean.py::joint_known_query_losses；M28-A23-RELATIVE-BOUNDARY-FIT.sbatch→run_query_relative_boundary.py::fit也调用同一supervised_readout。总损失按原window importance乘四项之和：local、known-query existence BCE、teacher boundary CE和termination BCE，四系数均1。四项都连接共享表示，原局部评分不能视为只用完整链标签训练的发射头；此处没有验证各项实际非零梯度或贡献大小。A43 fit_natural_competitor.py则只用关闭的有限候选几何和0／1／−1接structured_emission_loss，没有build_targets／local_loss或query teacher项；当前fit和协议不因此改动。

局部项也不是在所有位置强制一个精细类别的普通CE。canonical src/m28/labels.py::build_targets使用原R4窗口及其全部overlapping_transcript_ids，裁切primary CDS到窗口，保留完整转录本的coding offset；isoform分歧、不同gene重叠分歧、fine_label_unknown_intervals及不可调用DNA按原规则屏蔽。local_loss在精细类别已解时取−log p(role)，精细类别未知但区域已解时取−log sum p(role in region)，区域组为intergenic、CDS及intron；fine和region均未知的位置不计入分母。均值按本窗实际可监督位置计数，不按已知候选数计数。

原局部监督并非直接取自feature缓存中的旧标签。canonical run_query.py::load_train_inputs按原1536draw的1491个ID连接outputs/M28-LABEL-POLICY-R4下的windows_train.jsonl及其全部overlapping_transcript_ids对应的transcripts.jsonl；混合transcript目录先用train_records筛出TRAIN行，再解析和核对原六条TRAIN染色体、窗口／转录本身份及R4正例资格。随后supervised_readout以该R4窗口、catalog和缓存的原定向DNA调用build_targets。validate_cached明确不读取或使用缓存中的R2／R3标签字段。

cache_features.py的新写出分支保存features、sequence、window、model、revision，window仅为十字段几何；复用R3分支仍指向原artifact。这里核实的是写出规范和当前训练接线，不是每份历史tensor内容的独立读回。若后续检验原局部监督与结构竞争，必须复用原R4标签来源和mask，而非将旧缓存标签或A43整链正例直接当作逐碱基目标；原GEN特征与DNA仍可按原几何复用，不需要据此重新提取backbone特征。当前没有打开catalog或tensor、重建真实标签、选择混合loss权重或启动新fit，A43当前数据域与协议不变。

因此，复用原局部项会带入有限完整链配对以外的逐碱基信息，例如受原mask约束的裁切CDS／区域及背景监督。精细角色未知与候选链−1是不同粒度的保护对象，不能把后者直接转成逐碱基intergenic标签，也不能仅从现有2916个正链推断所有其余位置都是已知背景。加入原局部项需要明确其原R4标签来源、mask、位置分母及与整链项的相对归一化；相同1491窗、draw和更新数不足以宣称监督信息、容量或梯度尺度已匹配。

若A43完整结果仍不支持自然注释收益，保留原局部监督并与实际结构竞争共同学习可作为下一机制候选；目前不能判断它必要或有效，也不据此选择辅助权重、重开A36／A42、改A43 loss或启动新fit。后续实验必须先区分局部监督、query任务梯度、表示容量与竞争支持的变化，再检验收益是否到达完整raw输出；普通CE加结构解码本身不构成新颖性结论。

本节依据canonical当前源码和三份小fit配置metadata，未打开训练行、tensor、DNA、缓存、R4 catalog／参考或运行中的预测行，没有神经前向、梯度重放或optimizer。源码位置为scripts/experiments/M28-METHOD-RESTART/下的run_query_known_mean.py:49、run_query_relative_boundary.py:65、fit_natural_competitor.py:51、run_query.py:72／93和cache_features.py:43／63，以及src/m28/下的query_known_mean.py:65和labels.py:61／125。本机缺run_query_known_mean.py、run_query.py、cache_features.py及src/m28/query_known_mean.py，不复制镜像或以本机缺失作科学阴性证据。完整A43数值与五条件仍由原唯一监控在终态交接后处理，raw／DEV、新域及发布边界不变。

## 局部与整链监督共存的人工接线准备

local_structured_emission.py仅组合原labels.local_loss与structured_emission_loss，保留原精细／区域mask和位置分母、有限已知候选分母与−1排除，不从整链标签生成逐碱基背景。两项接同一十五类logits，显式返回各损失与各自监督计数；原window importance由未来调用者只乘总损失一次。混合系数无默认值，调用者必须显式提供；人工检查使用的1／1及2／3不是已选研究权重。没有新模型、参数、数据源、teacher任务或推理接线，也不修改A43原loss。

七项新增人工检查检验已知位置3／已知候选2的独立归一化与区域公式、共享FP32发射的梯度对独立闭式公式、链−1但局部区域可监督、局部全未知但链已知、双未知／空池的可微零，以及原build_targets在人工R4元数据上的细标签未知／N／padding规则和错窗／已知padding拒绝。直接解析仅限单exon人工竞争集，不据此证明完整语法partition、真实R4标签一致性、神经网络梯度规模、训练收益或模型创新；可微零也不是AdamW权重不变的证明。本准备不打开研究catalog、tensor、DNA、checkpoint或运行预测，不执行optimizer或真实标签生成。

唯一命令sbatch --parsable sbatch/M28-LOCAL-STRUCTURED-LOSS-CHECKS.sbatch已提交为13425410，06:51:05至06:51:14 UTC运行9秒，root／batch／extern均COMPLETED0:0。七项人工检查PASS，Ran7 in1.644s／OK；stderr包含正常unittest逐项结果，不是空文件。实际generanno／private／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS474072K；独立根outputs/M28-LOCAL-STRUCTURED-LOSS-CHECKS的STATUS／JOBID及唯一M28LOCALSTRUCTCHECKS_13425410日志、canonical Command／WorkDir匹配。新两AST／SBATCH语法、四个明确文件local／canonical逐字及06:50:18 UTC五项fresh调度／节点／预约／精确输出冲突检查通过。没有重提、延时、真实输入、optimizer、新fit或A43接线；该短作业已由主消费，不加入原36根监控。PASS仅支持本人工接线；A43实际FAIL已见结果文件，下一实验仍须先固定科学协议和公平对照，未选择混合权重，也未放行raw／DEV、新域或公开发布。
