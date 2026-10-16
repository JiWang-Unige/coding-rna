# A30 当前搜索竞争损失：执行记录

2026-10-06。按[a30_search_retention_protocol.md](a30_search_retention_protocol.md)首次固定协议执行并闭合。原fit的Slurm root仍为TIMEOUT，但完整4608行、末004608及model／AdamW经13399750核验；恢复读出13399752和CPU13399755均COMPLETED/0:0，无新增训练更新。原五项候选门槛0/5 PASS，关闭本配方，A32不执行；见[完整TRAIN结果](a30_search_retention_train_result.md)。下方运行中记录均为历史。

## 实际预检

sbatch --parsable sbatch/M28-A30-SEARCH-RETENTION-PREFLIGHT.sbatch返回13398207。提交前五项实时Slurm检查通过；gpu034的8卡用4、gpu035的8卡用0，均3090 24GiB；账户无活动作业，待排GPU无可运行资源竞争，没有相关维护截断或输出目录碰撞。沿用实际generanno环境，private／1RTX3090／2CPU／16GiB／10分钟。

root／batch／extern均COMPLETED0:0，实际2:50，batch MaxRSS2054216K。STATUS、JOBID及summary配置匹配。八项原metadata固定案例全部通过；完整fresh A3 649457参数相同、原四loss／分配／保护逐case相同，原自然搜索几何、顺序、分数和首次淘汰匹配。额外项及总损失梯度有限，两物种共享与query参数实际收到非零额外梯度；无正query案例额外项及梯度为零。所有模型参数不变，optimizer0，无监督checkpoint继承、新GLM、DEV、sealed或Setaria。

报告内部54.665159秒；其中包含两个模型比较、额外搜索、原搜索重复核对及反向，不能当正式每步速度。八项额外采集／带梯度回放分别0.919／0.490／0.000126／0.0895／0.532／0.292／0.000128／0.0867秒。实际peak_GPU_allocated_bytes407204864。27条已分配正query全有raw支持，首次淘汰均local_start_topk，cutoff均known_nonmatch；因此本实际预检不覆盖训练后全局或深层首淘汰，相关分数语义由已有组件测试及A28原搜索取证支持。此工程通过不表示候选覆盖或科学gate通过。

## 正式链

拟合保留原4608更新、末004608和每256步工程保存；每步保留全部正query采集记录及安全资格，未改变原四loss或Hungarian成本。readout先关闭1491窗原自然联合池，再加载TRAIN参考；CPU只读已保存结果，重算标签、负例资格、margin平均／importance weight和首次固定五gate，不做模型前向或搜索。

资源依据：A23原四loss拟合约64分钟；本轮metadata案例采集均不超过0.919秒，但未覆盖拟合后更深搜索，不能保证实际训练时长。fit请求2:30／1RTX3090／2CPU／16GiB，readout沿用已完成同范围约44分钟的2小时GPU上限，CPU沿用同范围约14分钟的30分钟／0GPU／2CPU／8GiB。预计新增checkpoint及日志低于2GiB，正式提交前核对实际存储和实时调度。所有新阶段独立目录，不覆盖旧实验。

本机两正式runner AST及三SBATCH语法、远程相同解析与语法检查均通过；窄范围上传成功。每次提交前五项实时调度检查均完成，没有GPU维护或并发上限冲突，新阶段目录不存在；BeeGFS可用139TB，满足新增产物。私有gpu034／035有可用3090，无需共享队列迁移。账户其他TE作业不在本轮监控或操作范围内。

实际提交sbatch --parsable sbatch/M28-A30-SEARCH-RETENTION-FIT.sbatch返回13398266，提交后RUNNING gpu034；sbatch --parsable --dependency=afterok:13398266 sbatch/M28-A30-SEARCH-RETENTION-TRAIN-READOUT.sbatch返回13398270，有效Dependency排队；sbatch --parsable --dependency=afterok:13398270 sbatch/M28-A30-SEARCH-RETENTION-RESULT-CHECK.sbatch 13398266 13398270返回13398273。fit自身核对已完成预检13398207及配置，不附已移出活动控制器的短预检依赖。日志为outputs/M28-PILOT-R5/logs/M28A30RETAINFIT_13398266.*、M28A30RETAINREAD_13398270.*、M28A30RETAINCHECK_13398273.*。

提交后整链快照：fit13398266 RUNNING gpu034／STATUS RUNNING／JOBID相符；readout13398270和CPU13398273均PENDING有效Dependency。实际首步额外margin1.533319、weighted_total16.620918、clip前总梯度65.527878，均有限，reference_injected_candidates0。日志已到64／4608，额外margin1.836537，证明新目标参与实际更新，不是性能结果。

## 运行时限风险

首64步训练内部181.566409秒，128步361.628800秒，第二段64步180.062391秒，约2.825秒／步；简单外推约3小时37分钟，超过原2:30上限。两个早期区段速度一致，仍不能证明最终耗时。主对话重新完成五项实时调度检查，确认当前资源无GPU维护截断，尝试仅延长运行中13398266至4:30：scontrol update JobId=13398266 TimeLimit=04:30:00，实际返回Access/permission denied，exit1；运行中扩时需要特权，见[Slurm scontrol文档](https://slurm.schedmd.com/scontrol.html)。没有重试权限命令、取消、重提、改科学命令或声称扩时成功，原2:30仍有效。

该风险不是科学FAIL。目前继续原作业并保留每256步checkpoint；若实际TIMEOUT，固定监控将回传，主对话先核对末可读的完整本候选model／optimizer和已flush日志，判断同一冻结4608曝光的工程恢复。不增加总更新次数、不选择更好checkpoint、不覆盖原失败记录，不预称恢复实现已完成。当前没有额外训练作业。

固定监控仅复用原对话01a0f71f-8845-7432-bd13-bf54c1808fb3／heartbeat coding-rna-r5-light-monitor，当前对象仅上述三根。同一原生心跳03:35:07.749 UTC保存ACTIVE，整个新prompt和id/name/kind heartbeat/target/15分钟/failed_runs_only/原created_at精确读回通过。原对话manual startup实际turn01a10f47-bf7b-77a0-9c12-df0404447251完成、error null、43.108秒，03:35:55.300 UTC快照为fit RUNNING、两个下游有效Dependency。实际turn_context03:35:22.809 UTC仍gpt-5.6-luna/low／相同项目，无覆盖或新对话。当前cohort定时执行及终态回传尚未发生，不能由配置或manual推定。旧A23／A25误读已纠正，白名单使用字面new_head_or_backbone_forward，不再取证旧作业。预检由主对话直接消费，没有另开监控。本机休眠或断线时不能保证即时唤醒。

后续04:09 UTC记录：上述是03:35时的历史状态；A30-only实际定时快照03:50:20.494Z、04:05:22.437Z均健康且未交接。A31独立caller读出13398521／CPU13398524已另提交，同一心跳04:08:52.926Z仅扩展为五根、两条独立整链，A30三根及字面白名单保持不变。新combined cohort的配置成功不代表其定时运行或终态交接已经发生。A30仍RUNNING，未改模型／曝光／已提交时限，未启动恢复作业。

## checkpoint 只读准备

2026-10-06。实际训练速度仍提示原2:30上限存在不足风险；本轮没有再次尝试扩时、取消、重提或修改运行中的runner。`inspect_search_retention_checkpoint.py`仅在CPU检查早已闭合的step001024和对应前1024条日志，不读取或采用其后的正在写入段。这是工程准备，不是选择恢复checkpoint或实施恢复训练；原4608逻辑更新、seed0、末004608、margin／系数及五gate不变。

本机两个新Python文件AST及SBATCH语法通过，三个文件窄范围上传完成后，重新完成五项实时Slurm检查并确认唯一输出未存在。private有CPU资源、无相关维护截断或并发限制。命令sbatch --parsable sbatch/M28-A30-CHECKPOINT-PREPARATION.sbatch返回13399061，private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，generanno，唯一outputs/M28-A30-CHECKPOINT-PREPARATION。

root／batch／extern全部COMPLETED0:0，实际2:30，batch MaxRSS749064K，STATUS／JOBID相符。四项检查PASS，内部0.717秒：合成线性模型恢复AdamW后下一更新与连续执行逐tensor相同，仅恢复权重却重置优化器则不同；错误超参数或moment拒绝接续；原draw顺序、importance weight及1536步epoch边界保持，checkpoint之后的丢弃残行不进入前缀；只有实际Slurm root TIMEOUT才允许考虑checkpoint选择，RUNNING／PENDING／COMPLETED／FAILED均拒绝。

真实旧checkpoint的649457参数结构及所有模型tensor有限；AdamW原超参数、moment形状和数值均可载入。对应日志1024行逐身份、epoch、draw序号及weight相符，拟南芥506／水稻518。下一原draw为全局step1025／epoch1／ordinal1024，窗口oryza_sativa|NC_089046.1|17461248|-，weight1.333333333333333；这是旧checkpoint的接续身份，不是当前运行位置或未来恢复起点。

46个带AdamW状态的参数张量，其应用计数范围781到1024。续跑必须保留各自的原moment及计数，不能统一重置为全局step1024，也不因计数不同擅自否定checkpoint。原训练路径并非每步接通所有query分支；本检查不将该范围解释为生物学机制。源checkpoint没有保存RNG状态，当前模型训练module没有启用dropout、原draw不重采样；本检查仍不证明跨进程GPU位级等价。

检查时源13398266的实际sacct root仍RUNNING0:0。真实模型检查没有forward、候选搜索或optimizer更新；optimizer短步只发生在合成小线性模型测试中。没有新GLM、DEV、sealed或Setaria读取，没有修改源JOBID／STATUS、checkpoint、日志、资源或协议。恢复训练尚未实现或提交。

若源实际TIMEOUT，主对话须重新核对届时最新已保存文件的完整model／optimizer及同一checkpoint日志前缀，按工程步序而非效果选择可读的完整状态，保留原失败产物和丢弃尾段记录，再决定同一4608逻辑曝光的接续。不能用本旧step001024检查直接宣称最新状态可恢复、从头另训或追加epoch。本短只读作业已由主对话消费，不加入原五根心跳，不创建新监控。

## 条件续跑实现

2026-10-06 05:20 UTC。上述只读准备之后，`resume_query_search_retention.py`已实现同一A30轨迹的条件接续，但没有提交恢复GPU作业。入口先核对原13398266的实际Slurm root为TIMEOUT，再读取届时最高已保存步的model／AdamW和相匹配日志前缀；不完整状态会停止，不按效果回退到较早checkpoint。源训练仍RUNNING，不能由时限风险推定已经失败。

接续保留原draw顺序、importance weight、原四loss／Hungarian保护、当前搜索损失margin1／系数1、clip1和原AdamW的逐参数moment及应用计数。模型移到CUDA后重新装载优化器，核对moment与参数的device相同。只复用已保存checkpoint对应的前缀，从下一原draw继续到逻辑step4608；末004608仍是唯一评价checkpoint，readout和CPU沿用原1491窗及五gate。新产物独立位于FIT／TRAIN-READOUT／RESULT-CHECK各自的RECOVERY-R1目录，不改源checkpoint、日志或失败标记。原checkpoint之后的尾段不进入最终接续轨迹；重算丢失尾段可能增加实际执行的更新与耗时，不能把逻辑4608误称全部物理执行次数。

本机三个Python文件AST和新SBATCH语法通过，四文件窄范围上传完成。提交前五项实时调度检查通过：private有CPU资源，无相关维护或并发限制，唯一输出未存在。`sbatch --parsable sbatch/M28-A30-RESUME-IMPLEMENTATION-CHECKS.sbatch`返回13399170，private／0GPU／2CPU／8GiB／5分钟，generanno；root／batch／extern均COMPLETED0:0，实际48秒，batch MaxRSS752032K，STATUS／JOBID相符。三项unittest全部PASS，内部1.130秒：RUNNING源在权重／TRAIN输入／输出绑定前拒绝接续；合成4608行准确恢复三epoch及统计、拒绝跳步；合成小模型的原四loss加带importance weight额外项及clip1更新与直接计算逐tensor相同。第三项额外项是合成可微函数，不是实际搜索margin资格的重新验证。

此作业没有真实TRAIN forward、搜索或optimizer更新，短更新仅发生于合成CPU模型；没有读取DEV、sealed或Setaria。它证明实现的限定工程行为，不证明实际TIMEOUT恢复成功、最新源checkpoint可读、跨进程GPU位级等价或科学gate通过。源RNG未保存而当前训练module未启用dropout的限制仍保留；自有中断状态接续不等于从其他实验继承监督checkpoint。该短检查已由主对话消费，不加入原五根固定心跳。

## 终态产物核对与只恢复读出

2026-10-06 06:04 UTC。原13398266 root为TIMEOUT0:0／2:30:08，batch及extern为COMPLETED0:0；不能用batch或STATUS将根状态改称成功。实际源training.jsonl恰有4608行、最后身份step4608／epoch3／ordinal1535，18个256步checkpoint至step004608均非空，summary completed=true／optimizer_steps4608／唯一末004608／原配置与revision，STATUS为COMPLETED。训练内部wall8937.342326秒；日志末打印完整summary，未见错误堆栈。旧13398270为DependencyNeverSatisfied、13398273继续等待；原记录和排队作业不修改、不取消。

这些事实使继续训练不再必要，但仍需验证末checkpoint的model／AdamW和完整日志。新reconcile_query_search_retention_fit.py在CPU核对原TIMEOUT、完整末summary、649457有限模型、AdamW超参数／moment／逐参数计数、全部4608行原draw与weight／三epoch统计和原监督曝光。只有实际CPU root COMPLETED0:0及相匹配产物才允许进入readout；部分summary或RUNNING源会停止，不能以batch成功替代根成功。

恢复只沿用原末模型和原readout／CPU评价实现；原FIT保持绑定源目录，新读出与评价使用独立RECOVERY-R1目录。没有TRAIN optimizer step、fresh seed、新checkpoint选择、候选注入、阈值或五gate变化。读出仍先完成全部1491窗reference-free自然池，再加载TRAIN参考；CPU重算原负例资格、margin、完整曝光和原五gate。任一科学FAIL仍关闭A30，不用超时理由救回。条件训练接续driver保留未启动，本次不使用它的fit入口。

核对作业private／0GPU／2CPU／8GiB／5分钟；readout沿用已完成同范围约44分钟的2小时／1RTX3090／2CPU／16GiB，CPU核验30分钟／0GPU／2CPU／8GiB。BeeGFS实际155TB可用，新增读出与短记录不超过原2GiB估计。固定同一心跳已由主对话暂停，并已消费原TIMEOUT事件，待取得实际新job IDs后替换原A30监控对象；A31两个DEV原作业和依赖保持不变。

六文件窄范围上传、本机和远程两Python AST／三SBATCH语法通过。每次提交前五项实时调度检查通过：private gpu034有空余3090、gpu035空闲，无相关维护、并发上限或输出碰撞。sbatch --parsable sbatch/M28-A30-TIMEOUT-FIT-RECONCILIATION.sbatch返回13399750；--dependency=afterok:13399750 sbatch/M28-A30-RECOVERED-TRAIN-READOUT.sbatch 13399750返回13399752；--dependency=afterok:13399752 sbatch/M28-A30-RECOVERED-RESULT-CHECK.sbatch 13399750 13399752返回13399755。

核对13399750 root／batch／extern实际COMPLETED0:0／41秒，batch MaxRSS737340K，STATUS／JOBID相符。三项合成检查PASS，内部0.004秒：RUNNING源在权重／TRAIN metadata前拒绝；未满4608的summary不能走只读出恢复；batch成功不能替代核对根COMPLETED0:0。真实终态核对内部3.695452秒，649457参数模型及全部tensor有限，46个AdamW状态的原超参数／moment／device合法、逐参数应用计数3492至4608，全部4608行身份／importance weight／三epoch统计和原监督曝光相符。没有模型forward、搜索、optimizer更新、源STATUS／JOBID／checkpoint／日志修改、DEV或Setaria读取。此通过不表示自然候选五gate或最终gene F1通过。

06:11 UTC实际13399752 RUNNING gpu034，13399755 PENDING有效Dependency，A31两个DEV推理继续RUNNING，其CPU依赖不变。同一原生心跳06:12:00.279 UTC仅改为本次恢复三根＋原四个DEV根；完整prompt／ID／kind heartbeat／原固定对话／15分钟／failed_runs_only／created_at精确读回通过，A31字面白名单未变。源root TIMEOUT是明确保留的历史，不能再次当新故障；新CPU summary也保留该源状态及核对job ID，不将原fit根改成成功。当前新cohort的manual启动、定时执行和终态回传还须分别取得实际证据；旧cohort实际06:01:50.980 UTC定时快照只能说明当时原fit仍RUNNING，不替代本新cohort验证。

后续实际readout stdout已到64／1491自然TRAIN窗，证明已进入原末模型的自然搜索，不是完整候选或gate结果。scontrol确认13399755提交参数为核对13399750／读出13399752，afterok:13399752有效、0GPU／2CPU／8GiB／30分钟，未误绑定旧失效readout。

同一watcher手工启动turn01a10fd8-82c7-7751-855d-bd591b6cf9a4实际completed／error null／33.503秒，06:13:40.412 UTC取证快照确认核对已完成、恢复readout和两个A31 infer运行、三个CPU等待有效依赖；实际turn_context06:13:29.961 UTC保持gpt-5.6-luna／low、同项目，未覆盖模型。新cohort配置和manual检查已验证，定时运行及终态交接尚未观察；不把旧cohort定时记录或manual当新定时成功，也不能保证休眠或断线期间即时回传。原协议、科学停止条件和数据边界不变；主对话等待整链结果再验证、解释和推进。

## 最终读出和停止判决

2026-10-06 07:44 UTC。readout13399752／CPU13399755实际root、batch、extern全部COMPLETED/0:0，分别1:00:30／27:35。保存自然池及标签各1491行、STATUS／JOBID相符，原曝光、model、标签资格及五gate核验闭合。原fit root TIMEOUT不重分类，恢复TRAIN optimizer更新仍0。

原五项全部FAIL，主对话按整数独立复算一致：正确链拟南芥1751→1750、水稻629→612；合计正确≥5 exon788→782；known nonmatch分别6305→7894、6261→7029。关闭本次A30配方，条件A32未执行；不做参数救回或把未执行记为失败。完整解释与证据集中于结果文件，不新增重复台账。

固定watcher定时turn01a11029-b85c-7b11-b28a-55bff4571ff6已完成/error-null，实际07:42:25.068 UTC向主对话发送一次A30结构终态交接，并暂停同一heartbeat且读回PAUSED。其结构通过不是科学PASS；五gate及停止判决由主对话处理。A31两条独立DEV链继续，由同一固定监控更新对象后等待，不新建聊天或自动化。
