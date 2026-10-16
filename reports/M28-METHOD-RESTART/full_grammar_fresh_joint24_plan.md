# 原 TRAIN 新抽样与联合编码的交互实验

2026-10-12 UTC执行更新：完整B流程13494895–13494899已由主核验并消费，全部30条root／batch／extern正常完成，六raw及四份输入DNA一致性通过。Pooled整链F1为70.302624%，比原joint提高3.339861个百分点，四格F1交互−3.987713个百分点；仍低于Tiberius、ANNEVO和Carbon-A。完整[结果记录](full_grammar_fresh_joint24_result.md)保留负交互与不同端点。C/D/E仍待完成，冻结三候选DEV尚未执行，不改变任何固定设置。

2026-10-11。补齐已有抽样×编码设计的第四格：同一固定 fresh4608 抽样序列是否与单窗口联合24576编码互补。依据为已完成原TRAIN对照、新抽样和联合编码三支路；不是依据DEV排名挑选新的权重。本方案在完整DEV结果已读取后固定，不能称为DEV前预注册。新的支路只进入原六TRAIN诊断，不加入原八方法DEV或独立十方法扩展。

## 固定四格与问题

| 单窗口编码 | 原1536 draws重复三遍 | 原seed1 fresh4608 draws |
| --- | --- | --- |
| 独立6144编码 | STRUCTURED_CONTINUE，评价13451625 | FRESH_DRAWS，评价13453583 |
| 联合24576编码 | JOINT24576，评价13453780 | 本支路，尚无结果 |

三个已有控制直接复用，不重复训练。每格均从原13434603末4608 C0及保存的AdamW状态独立出发，新增4608更新，仅末9216为正式读出；不接续任何9216权重。新格沿用原seed1固定4608 draws／4219 unique、importance1和顺序，复用准备13453348及全部fine事件支持13453408，不重新抽样、改变标注或增加窗口。上下文读取保持24576及相邻12288 stride、nearest-center ownership；唯一相对fresh改变是每个原读取窗口由四个6144块改为一次4096-token原生共同编码。

各格报告完整链count-derived precision、recall、F1，按染色体、两物种、pooled、species macro／worst区分，并保留background与CDS-base端点。编码在原抽样上的差值、编码在fresh上的差值和事前固定F11−F10−F01+F00交互量都报告，无论正负。单seed TRAIN交互是描述性、非线性F1尺度上的比较，不是因果互补、统计显著性或独立泛化。新抽样同时改变有限监督暴露及梯度序列，不将其缩写成纯覆盖效应。

## 不变的模型与范围

仅原拟南芥NC_003070.9／NC_003071.7／NC_003076.8和水稻NC_089043.1／NC_089044.1／NC_089046.1，原22,688 primary参考和原native基线。全部209语法不支持参考保留。DEV、held-out、sealed、Setaria和新物种均不读取。已读DEV不用于抽样、学习率、阈值、checkpoint、候选家族或后续训练选择。

复用原b048 GENERanno冻结backbone、原bidirectional SDPA flash与RoPE、六碱基token、BF16隐藏状态后FP32头；短窗口定向后右补N至24576，仅有效碱基进入原势函数。无native CDS先验。C0仍366191参数、原seed0起点；结构观察NLL、fine优先、coarse epsilon0.1、unknown边缘化、known-row分母、unit role gain、null0、cap128和FP64分区／skyline均不变。AdamW仍lr3e-4、weight decay0.01、betas0.9／0.999、eps1e-8、clip1，恢复全部旧moments及step4608。

## 实施与验证

独立入口为`fit_full_grammar_fresh_joint24.py`、`annotate_full_grammar_fresh_joint24.py`与`infer_full_grammar_fresh_joint24_train.py`，复用既有fresh训练和共同评价实现，仅绑定新的明确身份及联合编码生产函数；不修改任何已有支路或活动评分回放代码。

先在两物种各一个既定fresh视图上进行独立可丢弃4609更新，检测模型／AdamW接入、有限损失和冻结backbone。实际更新后严格raw adapter装载的势函数须逐张量等于训练生产接入；原联合编码八见证预检13453758必须完成，新的抽样和编码身份fixture须拒绝漂移。此预检不保存checkpoint，正式训练重新从原末4608及moments开始，不继承可丢弃更新。

预检通过并由主消费完成证据后，提交固定4608更新的正式fit、CPU全日志／三个checkpoint核验、六项原TRAIN raw FASTA和完整数组afterok共同评价；其后CPU四格汇总。只最后9216读出，6144和7680用于运行记录核验；不按中间loss或任何评价选权重。raw caller只接收FASTA与匹配最终权重，不读参考、标签或特征缓存。共用TRAIN评价器本身检查raw身份及导出DNA，但不直接比较不同支路的推理FASTA；新的四格汇总因此在解释前逐碱基比较四份caller输入与原TRAIN DNA，按既定大写、非ACGT→N映射保留长度，拒绝真实输入或参考分母错配。primary政策和固定native输出不变，不把同一路径或相同长度当作DNA相等。

预检单RTX3090／2CPU／16GiB／15分钟；正式fit相同资源／18小时，原fresh及joint实测总时间分别51,698及47,490秒，实际预检速率仅作资源核对，不改变更新数。CPU核验2CPU／4GiB／5分钟，六项raw每项单RTX3090／2CPU／16GiB／2小时、max3，原joint六染色体实测最长2,214秒，留出新评分下decoder复杂度的运行余量；CPU共同评价2CPU／4GiB／20分钟，四格汇总2CPU／4GiB／10分钟。全部使用现有generanno环境和唯一新日志名。fit文件上限96MiB，无全范围特征缓存；原始DNA、checkpoint、预测和冻结结果不覆盖。独立输出为`outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-{PREFLIGHT,FIT,VERIFICATION,TRAIN-RAW,TRAIN-DIAGNOSTIC,FOUR-CELL}`，预期合计小于2GiB。

常规同科学设置实施与明确修复按已有授权自主推进；失败或科学负结果如实报告，不自动增训、改变参数、改冻结DEV或取消旧作业。长等待只复用原固定监控对话和唯一心跳，原六项评分回放仍保留监控对象；不新建监控或重复自动化。完成全部运行与共同评价后再判断完整链竞争差距，当前尚无本支路结果，原论文goal未完成。

## 实际接入与预检

本方案及入口于2026-10-11约14:19 UTC固定，随后仅同步新增文件。主实际读取三个已完成TRAIN评价：13451625／13453583／13453780，pooled完整链F1分别55.7453%／63.0729%／66.9628%，三个参考均22,688，DEV／sealed均false；没有用新DEV指标设定实验参数。实际队列显示private可用24GiB RTX3090、作业数未触及80上限、当前无相关维护、输出无碰撞，BeeGFS可用147TiB；未修改同账户其他TE作业。

新预检13494740实际于14:21:03–14:25:30 UTC运行，root／batch／extern均COMPLETED／0:0，STATUS／JOBID一致、stderr空。主消费两条实际记录：既定fresh步骤1为拟南芥NC_003070.9／14389248／负链，步骤2为水稻NC_089046.1／14462976／负链；各自独立从4608更新至可丢弃4609，训练生产势函数与严格raw adapter视图逐张量相等。抽样seed／unique／epsilon／编码块漂移fixture拒绝；原joint见证13453758及fresh支持13453408均通过既有完成门槛。仅2次可丢弃更新，无checkpoint、正式训练或DEV／sealed读取；原起点文件未改。此结果证明接入，不证明正式性能。

新四格汇总的count-derived PRF、错误F1拒绝及零／负交互算术检查通过，四个新Python入口AST及六个SBATCH语法通过。该局部算术检查不读科研数据、不运行模型；真实四输入DNA比较和新格性能必须等待全部raw、共同评价及四格汇总完成。

14:29:38 UTC实际提交正式fit13494895、afterok核验13494896、afterok六项raw数组13494897_0–5（max3）、完整数组afterok共同评价13494898及其afterok四格汇总13494899。fit立即RUNNING于gpu034；主读取scontrol确认评价依赖为afterok:13494897_*，不能仅凭某一raw task通过进入评价。fit从原末4608开始，不继承预检的可丢弃更新；全部后续路径独立，既有控制和评分回放作业不动。目前只有提交／启动，没有新格末权重或性能结果。

14:31:52.174 UTC更新原唯一心跳coding-rna-r5-light-monitor，绑定同一原监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3，并保留原主01a0332e-e0a9-7aa2-90f9-159405701833。主直接读回ACTIVE／heartbeat／15分钟／failed_runs_only／created_at未变，双事件prompt逐项相等；原评分回放13494264／13494270仍保留。任一独立cohort整批完成或发生actionable问题才向主交接一次并暂停，由主消费后继续另一事件；不把单阶段成功或配好监控称作全流程完成。按[官方定时任务说明](https://learn.chatgpt.com/docs/automations?surface=app)，本机监控仍要求电脑开机、应用运行及连接可用。

双事件人工首检turn01a12b61-6b67-7831-8c91-12ee6a17ecf2实际14:32:47–14:33:24 UTC运行，14:33:01检查；真实定时turn01a12b62-c392-7330-8782-e2565b4c0931于14:34:15–14:34:52运行，14:34:28检查。主读取实际turn_context，两者均gpt-5.6-luna／low；真实heartbeat输入14:34:15.374Z包含原automation_id及两个当前事件。首检和本轮定时追加分别验证，未读取科研内容、交接或暂停。旧评分回放当时只有task0仍RUNNING，汇总等待有效依赖；新fit运行、其余阶段有效Dependency。两事件整批终态回传均尚未观察。

主另读实际fit started确认job13494895、configuration与通过的预检逐项相等，起点仍13434603／原4608及moments，不是预检权重。scontrol逐项确认核验、raw、完整数组评价及四格汇总的afterok链。启动检查时training.jsonl仍为缓冲期0字节、stderr空，不能据此判定停止或宣称已完成更新数；完整日志和末权重核验待终态。

## 原心跳恢复后的实际定时运行

原评分回放已于2026-10-11T14:49:55Z由同一监控对话交接并由主完成科学分析；原唯一心跳于14:56:34.801Z恢复ACTIVE，仅监控本fresh×joint流程，详细结构白名单、原target、15分钟和通知策略不变。恢复后的真实scheduled turn01a12b85-18f4-7162-9f44-9b35c28e3620于15:11:45–15:12:07 UTC追加到原监控对话；实际heartbeat输入15:11:45.394Z包含原automation_id及仅事件B的范围，实际turn_context15:11:45.495Z确认为gpt-5.6-luna／low。15:12:00检查fit13494895 RUNNING于gpu034，核验、六项raw、共同评价及四格汇总均PENDING有效Dependency。没有检查已闭合事件A，没有读取科研内容、回传终态或暂停；连接失败0，DONT_NOTIFY静默结束。配置恢复与恢复后的实际定时追加现均有直接证据，但本支路末权重、完整评价和终态交接仍未完成，不能记为科学成功。

## 最终拟合与权重核验收据

2026-10-12 UTC主核对实际TZ=UTC sacct：fit13494895于2026-10-11 14:29:38至2026-10-12 03:33:24运行，核验13494896于03:33:24–03:34:49运行；各自root／batch／extern共6条记录均COMPLETED／0:0，STATUS、JOBID和非空summary.json一致。拟合实际configuration与started及verification逐项相同，从13434603末4608及原AdamW moments独立出发、追加4608、累计9216；原seed1的4608 draws／4219 unique、每pass每物种768 draws保持。训练参数仍366191，backbone冻结，无native CDS先验、缓存特征、原targets改变、DEV或sealed使用。

CPU核验13494896已经读取全4608条日志并核对抽样／观察／importance顺序、已记录目标算术和三个6144／7680／9216 checkpoint；所有权重和AdamW moments有限、严格装载通过。主依据实际核验器的检查范围及完成收据消费这部分工程证据，没有在登录节点再次加载权重、运行native forward或进行新的optimizer update。正式读出保持[最终9216收据](../../outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FIT/summary.json)与[独立CPU核验收据](../../outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-VERIFICATION/summary.json)，不按中间loss或DEV选择权重。

六raw均于03:34:50UTC开始；实际sacct任务映射与各自started.json一致，均指向同一最终step_009216.pt和source_fit_job13494895，reference_used／DEV_used=false：

| Task／染色体 | 实际numeric root |
| --- | --- |
| 0／NC_003070.9 | 13508070 |
| 1／NC_003071.7 | 13508071 |
| 2／NC_003076.8 | 13508072 |
| 3／NC_089043.1 | 13508073 |
| 4／NC_089044.1 | 13508074 |
| 5／NC_089046.1 | 13494897 |

此时六raw均RUNNING，未据运行中的ExitCode字段判定成功；原完整数组评价13494898及四格汇总13494899尚未消费，四格输入DNA相等性和准确率／交互量仍待这些阶段。B完整流程保持未完成，原监控handoff B=NO；不因此开始冻结B/C/D DEV评价、不新增作业或改变正在运行的B/C/D/E参数。该部分核验以后直接复用，仅在新失败或具体未决疑点下复查。
