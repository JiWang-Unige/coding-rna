# 冻结上下文上的局部发射修正

2026-10-07北京时间，首次实际梯度检查前固定。A41检验保留A31的冻结单位局部路径势、在原A3上下文上学习零初始化线性角色修正。一次固定拟合13420820、自然读出13420821及CPU核验13420824已全部闭合；六项工程证明通过，但原五科学条件仅2/5通过，本版本FAIL关闭。原机制及五条件未改，不修改A40协议或作业。

## 与旧实验的区别

A26在fresh链头末端加入28维候选mean和count证据，并不保留已训练A31最终分数。A36从头学习上下文及局部发射；A37至A39从头学习上下文及链头；A40把冻结局部证据作为新链头的输入，也没有固定的A31输出项。A41固定原上下文和原发射，将零初始化修正直接加到原路径势，因此初始状态必须逐候选等于A31。它不混合A38或A39的已训练分数，不重新判定任何旧FAIL。

这个修正等价于在固定上下文上学习额外的线性局部角色边际，并非新的非线性表示或已证实的论文创新。从实数表达看，它没有扩大既有仿射角色边际的函数类；实际实现仍分别进行原发射和修正的FP32计算、FP64累计，不能默认为合并FP32权重后逐值相等。初始保留不保证训练后保留原已选正确链，已有监督A3曝光与新增可训练参数均须显式报告，不能称为全fresh或匹配曝光的因果消融。

## 固定评分定义

原A3末004608全649457参数冻结，实际只前向其SharedContext和15类segmentation，不运行query搜索。对原H128上下文h，加一个Linear128至14及bias，共1806参数，全部初值0；intergenic修正固定0。各状态修正为delta_s(h)=W_s h+b_s，s为1至14。source上下文和原局部logits均detach，只有这1806参数可有梯度。

每候选gain为原FP64 A31角色减intergenic求和，加FP32修正按同一原角色路径的FP64累计。两项分别累计后相加，不先合成FP32局部logits造成额外舍入；原相位、start和stop覆盖优先级、半开跨度、相触、valid length和padding规则不变。原gain0、null0与span非重叠DP不变，没有长度归一化、角色系数、温度、margin或阈值扫描。零修正时分数与本次冻结源直接计算值必须严格相等；本次源与已保存A31值仍使用原rtol1e-12和atol1e-8核验，并须精确复现原已选索引。

采用未修改的(mean-known BCE加mean-known兼容集合NLL)/2；两项保持原known数量归一。标签-1不参与两项监督及分区，独立gain直接梯度0，但共享修正参数的更新可能改变其分数，推理继续保留unknown。空池和全unknown保留可微零。zero-start的唯一可训练层直接接收梯度，没有被零输出层阻断的可训练上游；冻结上下文零梯度是设计要求，不是梯度遗漏。

## 原数据与结果边界

复用原六条TRAIN染色体、1491唯一窗、1536draw顺序及importance；自然池16321条，标签1、0、-1分别2380、12566、1375。没有536参考注入、新候选、新GEN前向或新R4。原2916 eligible reference/window配对和1086至少5-exon配对分母保持，不把候选配对当独立基因或完整染色体F1。

一次独立拟合固定seed0、原draw三遍4608 AdamW更新、lr3e-4、wd.01、betas(.9,.999)、eps1e-8、clip1；保存三checkpoint，仅末004608自然读出。零更新预检和下文接线检查实际通过、A40整链闭合后才提交，不根据A40中间分数改当前修正定义。A41修正定义及上述科学选择在A40结果读取前已固定。

后续原五条件不降级：拟南芥正确严格高于1531，水稻正确严格高于478；known_nonmatch分别不高于347和251；pooled至少5-exon正确严格高于700。任一FAIL关闭本版本，不救loss、角色权重、阈值、seed、epoch或checkpoint。PASS仅支持另固定完整raw-FASTA协议，不自动执行DEV、独立物种或Setaria，不代表SOTA或论文成功。此固定池修正不能生成原池缺失的正确链。

## 预检与资源

五项合成检查核对零修正完整保留原分数及DP选择；非零修正实际进入原角色求和、padding排除及直接算术梯度；zero-start有非零可训练梯度且unknown不监督；空池和全unknown为可微零；已正确饱和的冻结gain可给有限零梯度而非接线失败。实际预检复用原population和较强A31完成收据、八个原metadata窗口，不按新gain选窗；逐窗检查原初始分数及已选索引、有限修正梯度、unknown直接梯度零、冻结源无梯度、两组权重前后相同。零梯度case及非零候选gain导数分别记录，八样本中两物种均须实际出现非零修正梯度，不替换饱和样本。optimizer0，无checkpoint；不检查不存在的fresh上下文梯度。

入口preflight_emission_residual.py，唯一输出M28-A41-EMISSION-RESIDUAL-PREFLIGHT；generanno、private、1RTX3090至少24GiB、2CPU、16GiB、5分钟、no-requeue及唯一日志。沿用同系列实际SBATCH环境，不把旧框架通用coding-rna环境或自动决策恢复为研究要求。预计小型JSON低于1MiB，不建hidden缓存或全量logits。任何失败保留原输出并报告，不自动重跑或转入拟合。

A40原13418488、13418506、13418507及固定轻量监控保持不变。此短预检由主对话取证，不创建新监控，不覆盖原结果、不取消作业、不提交Git、不对外发送或发布。

## 实际零更新结果

一次提交sbatch --parsable sbatch/M28-A41-EMISSION-RESIDUAL-PREFLIGHT.sbatch，作业13420270于2026-10-06 22:25:19至22:25:59 UTC运行，root、batch、extern均COMPLETED0:0，40秒。实际private、gpu035、1RTX3090、2CPU、16GiB、5分钟上限，Requeue0、Restarts0；canonical Command、SubmitLine与WorkDir相符，batch MaxRSS1802040K。提交前五项活调度及预约检查没有相关硬约束或输出冲突，五个新增文件与canonical逐字一致，Python AST及SBATCH语法通过。

五项合成检查全部PASS，0.495秒。原八metadata窗口全部完成，实际初始gain及DP已选索引不变，与原保存A31分数最大误差0；八窗修正梯度norm均有限非零，unknown直接gain梯度全部0，冻结A3无梯度，两组state前后逐tensor相同。两物种均有非零修正梯度，实际zero-gradient case为0；两个水稻窗口只有5和8个角色行有非零weight梯度，其余六窗14行，不能把整体非零说成每个角色在每窗都被监督。optimizer0，无checkpoint。八窗内部9.310秒、GPU峰值116326400字节，只是接线工程时序，不是速度或性能benchmark。

完整原人口和较强基线身份核对通过：1491唯一窗、1536draw、16321自然配对及2380、12566、1375标签不变；draw重复支持16806配对及2473、12935、1398，不能混同唯一窗支持。完整参考分母2916和1086保持。主核对八窗保存的等权复合loss算术、完成标记、JOBID和真实allocation，与summary一致。这是零更新保存记录核对，不是独立神经前向复现或科学PASS。

该零更新根仅有STATUS、JOBID、summary.json，summary11748字节；没有新缓存、logits或权重。该阶段A40输出、新GEN、候选搜索、新R4、DEV、Setaria和A33均未访问，旧结果及作业未改。13420270闭合时尚无A41拟合入口或训练结果；后续实施见下节，不把预检loss或协议存在当作拟合结果或授权。

## 固定拟合和完整结果接线

fit_emission_residual.py只训练margin.weight和margin.bias，冻结A3全649457参数；每draw重放同一冻结上下文，不建立隐藏缓存，不解码新候选。三个checkpoint同时保存1806参数的修正、实际AdamW状态及原冻结源；空池和全unknown仍执行原AdamW过程。两参数tensor的实际step应为1536、3072、4608，不意味着每个角色行在每窗均有监督。

正式拟合前的READOUT-CHECKS先运行五项新人工adapter检查，再用旧JSON完整核对1491窗、16321自然配对、较强A31选择及2916和1086分母；A41占位仅为原基线，不是新结果。新检查检测候选缺失、乱序、reference用于评分或陈旧选择；保护unknown和完整分母；FP64分量求和被误转FP32；缺失零更新证明；错误optimizer population及step。不打开新checkpoint、特征tensor、模型前向或optimizer。准备通过后只提交一次fit和其afterok自然读出及CPU核验；任何失败保留，不自动重跑、续训或延长。

readout_emission_residual_train.py仅用末004608及原1491自然池，逐窗保存原FP32十五类发射、FP32十五列修正（intergenic列恒0）及两项FP64链势。全部1491潜势文件及自然分数关闭后才连接原保存标签。输出仅A31_A3和A41两臂、原DP及完整分母，没有原生DEV或新搜索。

verify_emission_residual_train.py在CPU重核4608实际draw和importance、三遍固定loss算术、三个有限FP32修正checkpoint及实际optimizer step、三份冻结源与原A3逐tensor相等。再逐候选从保存的两组FP32局部潜势重算FP64分量、原A31分数和最终和，重算原标签配对、DP、完整分母与五条件。沿用原局部势跨设备rtol1e-12及atol1e-8；保存的两个FP64分量之和与JSON最终gain必须严格相等，不加新的FP32量化容差。此核验不从保存上下文重放线性头，更不是独立神经前向或生物学因果验证。

六项实际CPU证明分别覆盖原draw及4608更新、固定复合loss与importance算术、三修正checkpoint及optimizer曝光、三冻结源状态及optimizer排除、保存局部潜势与FP64 gain算术、原标签DP基线分母及五条件。科学PASS还需全部原五条件，而不是这六项工程证明。

generanno、private、no-requeue及独立A41根不变；CPU接线检查0GPU、2CPU、8GiB、5分钟，fit1RTX3090、2CPU、16GiB、1:45，natural同GPU配置30分钟，CPU核验0GPU、2CPU、8GiB、15分钟。fit时限参考预检七个热窗口平均约1.064秒外推约82分钟再留操作余量；非随机八窗口不作速度benchmark或时间保证。三checkpoint及训练日志预计低于32MiB；自然两套局部潜势按全部24576bp上限约4.4GB或4.1GiB，提交前检查实际可用存储。CPU核验时间留同系列实际7至9分钟核验的操作余量。长等待只复用原固定监控对话和心跳，不创建新监控。

## 实际接线与一次提交

接线13420720于2026-10-06 22:50:51至22:51:25 UTC完成，root、batch、extern均COMPLETED0:0，34秒，private0GPU、2CPU、8GiB，batch MaxRSS724728K。五项新adapter检查PASS，0.062秒。完整原1491窗、16321候选、较强A31配对及2916/1086分母通过；A41占位仅原基线，新checkpoint、特征tensor、模型前向和optimizer均未使用。原13420270零更新证据及A40三作业、七项CPU证明闭合后提交；A40判FAIL3/5，不救旧模型。

正式链仅一次提交：fit13420820于2026-10-06 22:57:16 UTC提交、22:57:17 UTC实际开始，private-teodoro-gpu、gpu035、1RTX3090、2CPU、16GiB、1:45，RUNNING；natural13420821于22:57:22 UTC提交，afterok13420820，export FIT_JOB13420820/CHECKS_JOB13420720，1RTX3090、2CPU、16GiB、30分钟；CPU13420824于22:57:29 UTC提交，afterok13420821，export FIT_JOB13420820/READOUT_JOB13420821，0GPU、2CPU、8GiB、15分钟。两个下游PENDING/Dependency正常，实际Command、SubmitLine、WorkDir、唯一日志及no-requeue/Restarts0均已核对。提交前实际分区、队列、节点和预约没有相关硬约束；三个输出根与日志无冲突，scratch可用157T，不把共享余量当独占资源保证。

根依次outputs/M28-A41-EMISSION-RESIDUAL-FIT、M28-A41-EMISSION-RESIDUAL-TRAIN-READOUT、M28-A41-EMISSION-RESIDUAL-RESULT-CHECK。仅主在完整结果闭合后按原五条件解释，作业启动不等于训练成功、科学PASS或论文结果。复用原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3与coding-rna-r5-light-monitor，排除已消费A40及A41两次准备。完整接线代码与SBATCH提交前逐字同步canonical，Python AST及SBATCH语法通过；无Git提交、旧输出覆盖、新数据访问或对外发布。

## 实际完整结果与关闭

2026-10-07 00:31:12 UTC固定监控交接后，主对话读取实际Slurm、完成标记及两份结果summary。13420820的root、batch、extern均COMPLETED0:0，2026-10-06 22:57:17至23:47:55 UTC，50:38，batch MaxRSS2033688K；13420821同样全部COMPLETED0:0，23:47:55至2026-10-07 00:05:40 UTC，17:45，1528464K；13420824全部COMPLETED0:0，00:05:41至00:17:32 UTC，11:51，827772K。三STATUS为COMPLETED、JOBID逐项匹配，没有重跑或延长。

CPU六项实际证明均true：原draw及4608更新、固定loss与importance算术、三修正checkpoint及optimizer曝光、三冻结源相同且不在optimizer、FP32保存局部势至FP64 gain、原标签/DP/基线/分母及五条件。margin.weight及bias实际step均4608，三checkpoint为001536、003072、004608；source649457权重不变、无梯度和更新。原A31分数最大误差0，结果重算与自然读出相同。这是保存结果工程闭合，不是独立神经前向或科学PASS。

| 原条件所用计数 | A31_A3 | A41 | 条件 |
| --- | ---: | ---: | --- |
| 拟南芥正确 | 1531 | 643 | FAIL |
| 拟南芥known_nonmatch | 347 | 155 | PASS |
| 水稻正确 | 478 | 144 | FAIL |
| 水稻known_nonmatch | 251 | 153 | PASS |
| pooled至少5-exon正确 | 700 | 174 | FAIL |

完整2916 eligible及1086长链reference/window配对分母不变，原自然池16321及2380/12566/1375标签不变。A41拟南芥已选915条（643正确、155known_nonmatch、117unknown、139长正确），水稻331条（144、153、34、35）；合计1246条、787正确、308known_nonmatch、151unknown。原已选3054条、2009正确、598known_nonmatch、447unknown。正确净减1222、长正确净减526；降低known_nonmatch不能补偿正确和长链的大幅丢失。

结论为TRAIN科学FAIL2/5，TRAIN_preparation_supported=false，按原停止条件关闭A41，不更换loss、阈值、角色权重、seed、epoch或checkpoint，不转入raw、DEV、独立物种或Setaria。上述总计数本身不区分非正势与正势竞争，完整保存分数诊断见末节；这次失败不代表所有局部发射模型无效。A42继续使用其原先固定的A3源及搜索定义，不混入A41修正或按本结果改规则。固定监控已交接并暂停，本轮事件只消费一次。

## 历史监控启动记录

A41完整TRAIN运行期间同heartbeat曾更新为ACTIVE；完整prompt、ID/name、kind heartbeat、固定target、15分钟、failed_runs_only及原created_at逐字段读回相同。manual监控turn01a11372-b67a-7031-bd0e-d45ae5526070完成、error null、47.063秒，确认当时fit RUNNING/gpu035及下游正常Dependency，没有主交接；实际turn_context2026-10-06T23:00:47.454Z仍gpt-5.6-luna/low/原cwd，无模型覆盖。该快照只证明配置和一次手动检查，不是定时或终态证据；实际终态回传和关闭见上一节。

A41首次定时触发已验证：2026-10-06 23:15:16.410 UTC的原生heartbeat包含同automation ID，实际追加到固定监控对话的turn01a1137f-f93c-7b81-8d2c-884eb970edf3；该turn完成、error null、38.428秒。对应turn_context为gpt-5.6-luna/low/原cwd，只记录当时fit RUNNING及下游正常Dependency，没有主交接或指标解读。该次定时证据与最终回传分别记录，不把历史运行快照当当前状态；同监控目前已复用为A42整链。

## raw FASTA 接口准备

infer_emission_residual.py复用原GEN四块编码、A3 greedy加beam4、窗口和span-owner规则、同链去重、原非重叠DP及GFF3输出。每oriented window只编码一次GEN、一次冻结A3 shared；该同一上下文用于原自然候选、原十五类发射及1806参数修正。保持两套FP32角色势分别FP64累计后相加，不把合成局部logits或最终gain转FP32。候选、来源、两项gain、null0及非正候选均保存。负链仍用原坐标和CDS phase，不产生UTR或alternative-isoform覆盖主张。

实际TRAIN三作业闭合、配置和六CPU证明相符、原较强A31基线与两物种完整分母相同、五条件重算全部PASS，以及另固定执行范围，仍是任何真实原生运行的必要前提。假收据和文件存在不代表授权。当前独立DEV协议、A41原生infer/eval根均不存在，也未提交原生执行SBATCH；上述入口本身不触发DEV访问。正式fit、natural和CPU代码及作业未改。

evaluate_emission_residual.py先核对完整原8690窗口/四scope、冻结源和末004608身份，再逐窗配对原A31_A3的自然几何、顺序和来源；原冻结FP64分量须匹配同raw输入的A31分数，沿用rtol1e-12/atol1e-8，最终两项FP64和严格相等。原DP、saved counts、assembled prediction及GFF坐标phase关闭后才读原参考。未来评价保持全部7728 primary-CDS参考和各物种/背景分母及结构诊断；不把TRAIN unknown或原生非exact预测直接视为生物学假阳性。同尺对照包含ANNEVO、Tiberius、Helixer、A31_A3、A31_A23及既有内部方法，报告分物种、pooled、macro及worst species。它是已复用DEV的评价接口，不是独立泛化或公平速度证据。

原人工检查13421014于2026-10-06 23:09:00至23:09:18 UTC FAILED1:0，18秒，private0GPU、2CPU、8GiB，batch MaxRSS776188K；extern COMPLETED0:0，STATUS为FAILED exit1。七raw检查中六PASS、一ERROR：missing-separate-contract样例先替换ROOT为临时目录，再构造假checkpoint相对路径，触发pathlib ValueError；评价四项尚未执行，不记为FAIL或PASS。仅把假收据构造移至ROOT替换之前，评分、协议参数及任何科研作业未改。原失败根与日志完整保留，不在原根重跑。

一次独立CPU修正检查13421042于2026-10-06 23:10:58至23:11:25 UTC完成，root、batch、extern均COMPLETED0:0，27秒；private0GPU、2CPU、8GiB、2分钟上限、gpu034，batch MaxRSS775184K，Requeue0/Restarts0，canonical Command、SubmitLine、WorkDir、STATUS/JOBID相符。七raw检查PASS0.183秒，四评价检查PASS0.373秒。覆盖实际五条件/六证明/源身份/分母与加载前阻断、零和非零修正/FP64舍入/空池、每方向一次编码及负链phase、保存池/评分/DP/GFF/全尺对照。均为人工样例、模型、GFF和未来假PASS收据，没有研究checkpoint、缓存、真实FASTA或reference访问，也没有optimizer。原根M28-A41-EMISSION-RESIDUAL-RAW-CHECKS及独立R1根保留；新文件提交前逐字同步canonical，四Python AST及SBATCH语法通过，修正仅影响一项fixture及R1输出日志路径。

四评价检查包含一条既有scripts/eval_structure_diagnostic.py第51行open(path)未显式关闭的ResourceWarning，在人工GFF检查中可达；没有造成十一项检查失败，不据此判断科学结果失效。本轮未扩改共享参考解析器；实际原生执行前如该警告影响文件句柄或运行，再作针对性修复。以上是接口准备，不是完整注释运行、训练PASS、外部基线提升或论文结论。

## 关闭结果的分数诊断

diagnose_emission_residual_attrition.py只读取已闭合A41及较强A31的自然分数、保存标签与summary。先核对原三作业COMPLETED0:0、配置、末004608、六CPU证明及原五条件FAIL，再访问分数；FAIL允许解释失败，不放行raw。逐窗核对原几何、分量严格FP64和、冻结分量与原A31相符，以及两臂DP计数和全部2916/1086分母；未知标签仍独立，不重建TRAIN catalog或标签。

分类复用原transitions：保留、新增、原已选变为非正、原已选仍正但未选，以及两类从未选择。额外只记录仍正却丢失的正确链实际与新已选哪些标签集合重叠，沿原半开gene-span和相触规则；核对移除这些竞争链并替换为该链不能提高保存DP目标。此处不作新的decode、阈值扫描或角色消融，不把未知竞争链称为负例，也不从描述性分类推断训练loss的因果机制。

唯一命令sbatch --parsable sbatch/M28-A41-EMISSION-RESIDUAL-ATTRITION.sbatch；generanno、private、0GPU、2CPU、8GiB、5分钟、no-requeue、独立根outputs/M28-A41-EMISSION-RESIDUAL-ATTRITION和%x_%j日志，先运行两项新增人工接线检查再执行完整保存JSON诊断。参考A40同类纯JSON诊断实际短于1分钟，5分钟留启动和文件读取余量；预期小型summary低于1MiB，不新增权重或潜势。任何失败保留并报告，不在旧根重跑，不自动救科学FAIL。A42作业、评分和监控对象均不改；该短诊断主取证，不创建或改变后台监控。

## 正确链丢失主要来自非正势

唯一诊断13422424于2026-10-07 00:57:46至00:58:15 UTC运行，root、batch、extern均COMPLETED0:0，29秒；private gpu034、0GPU、2CPU、8GiB、5分钟上限，batch MaxRSS743360K，Requeue0/Restarts0及canonical Command/SubmitLine/WorkDir/唯一日志相符。STATUS COMPLETED、JOBID13422424，summary8948字节；两项新人工检查PASS0.008秒。提交前四相关文件逐字local/canonical一致，Python AST与SBATCH语法通过，五项活调度及预约无相关硬冲突，新根和日志无冲突；不覆盖旧结果或用户M20改动。

完整1491窗、16321自然链配对及两臂保存计数、2916/1086分母复现。冻结分量与原A31最大误差0，两项FP64分量严格和、实际已选重叠竞争集与替换目标、两臂和完整分母三项诊断证明均true；主再核对两物种及pooled共12分类组的pool分解、保留/新增/丢失净变化，以及正势丢失与竞争标签集合计数相符。以下是同一TRAIN窗口/链配对，不是独立基因或独立测试。

| 配对组 | 原已选 | 保留 | 非正势丢失 | 仍正但未选 | 新增 | 新已选 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 拟南芥正确 | 1531 | 623 | 897 | 11 | 20 | 643 |
| 水稻正确 | 478 | 131 | 335 | 12 | 13 | 144 |
| pooled正确 | 2009 | 754 | 1232 | 23 | 33 | 787 |
| 拟南芥至少5-exon正确 | 548 | 139 | 407 | 2 | 0 | 139 |
| 水稻至少5-exon正确 | 152 | 35 | 110 | 7 | 0 | 35 |
| pooled至少5-exon正确 | 700 | 174 | 517 | 9 | 0 | 174 |

原已选正确共丢失1255配对，其中1232变为非正、23仍正但未选；新增33后净减1222。长正确丢失526，其中517变为非正、9仍正但未选，新增0。冻结正势未变而最终和非正，直接说明已学修正抵消并压过原正势；主要丢失不是原DP漏选仍有正势的链。保持gain0时，仅改变后处理竞争规则无法恢复这些非正势原正确链，不据此改变原阈值。

23个仍正却丢失的正确配对中，19个实际已选重叠竞争集合只含known_nonmatch，4个只含unknown；长正确对应8和1。这里统计丢失链对应的竞争集合，不是19或4个独立竞争基因；没有label1、混合或无已选重叠组。移除整组竞争链并替换为该链不能改善原保存DP目标，故不能把这23个归因于DP算术错误。unknown仍是受保护未知，known_nonmatch也不等于生物学假阳性。

known_nonmatch原598：保留99、非正丢失472、正势未选27、新增209，最终308；unknown原447：保留121、非正丢失318、正势未选8、新增30，最终151。因此非匹配净减少290不是无代价筛除，而伴随正确和长链的严重丢失。

本诊断确认保存分数层面的直接失败路径，不识别BCE或set NLL哪一项的因果贡献，也不证明固定上下文或全部发射模型不可用。A41继续FAIL关闭；不再对同一仿射修正作阈值、checkpoint、角色权重或loss补救。A42固定原A3势、无A41修正，继续检验完整搜索是否增加正确和长链；它针对搜索支持而非本次新修正造成的负势，不因本诊断改其科学规则。没有权重/特征tensor、神经/GEN前向、optimizer、新搜索、新R4、TRAIN catalog/reference、DEV、独立物种、Setaria或A33访问，也没有raw运行或对外发布。
