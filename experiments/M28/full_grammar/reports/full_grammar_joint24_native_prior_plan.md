# 原 TRAIN 联合上下文与双头 CDS 先验的交互实验

2026-10-09。在读取组合结果之前固定问题：联合24576原生编码与两个冻结native CDS头的角色先验是否互补。组合caller和原生预检已通过，现接入事前固定第四格的正式训练及原六TRAIN染色体评价；不启动额外控制、扫参、扩大数据或宣称创新成立。上下文单独支路已有全部1491固定权重TRAIN视图的表示证据，双头已有六位置语义和真实梯度接线证据；它们仍不是末权重整染色体性能或互补性结果。

## 2026-10-11实际四格结果

组合fit13457841、核验13457846、六raw13457849_0–5及评价13457852已全部COMPLETED 0:0；评价于10-10 UTC05:12:08结束。预定四格汇总13460976于05:12:23完成，root／batch／extern成功，STATUS／JOBID／summary／comparison齐全。主核对全部四格来源、22688参考、主要测量与已完成评价直接相等、micro／macro／worst及六种预定差值算术；不按分数移除任何格。

四格pooled链F1依F00／F01／F10／F11为55.7453%／64.3390%／66.9628%／62.4512%。组合比联合编码单独支路低4.5115百分点、比先验单独支路低1.8877；预定F11−F10−F01+F00为pooled−13.1052、AT−9.9675、Rice−23.5104百分点。组合仍比结构控制高6.7059百分点，且纯背景链190少于两个组件的259／320，不能省略这项端点权衡。当前不支持在完整链F1尺度上声称两个改动互补，但也不把单seed TRAIN描述量升级为因果拮抗、显著性或泛化。完整记录见[结果工作稿](../../docs/M28_full_grammar_results_working.md#reference-recovery-losses-and-additional-fit-supervision)。无新训练、调参、DEV候选变更或封存访问。

## 固定四格比较

| 原生编码 bp | native CDS 角色先验 | 末权重来源 | 实际 fit job |
| --- | --- | --- | --- |
| 6144 | 无 | STRUCTURED_CONTINUE9216 | 13451136_0，运行中 |
| 6144 | 双头 28 系数 | NATIVE_CDS_PRIOR9216 | 13453177，运行中 |
| 24576 | 无 | JOINT24576 末9216 | 13453777，运行中 |
| 24576 | 双头 28 系数 | 新组合末9216 | 13457841，运行中 |

三个已有控制直接复用，不重新训练。四格均从原13434603末4608 C0及原AdamW状态出发，原两物种六TRAIN、1536 draws／1491 unique及三遍顺序／importance保持不变，新增4608更新，正式只末9216。原fine／coarse／unknown、known分母、epsilon0.1完整集合结构NLL、24576读取窗及12288邻窗、nearest-center owner、unit role gain、null0、cap128及FP64分区／skyline规则不变。单seed0，AdamW lr3e-4、weight decay0.01、betas0.9／0.999、eps1e-8、共同clip1，不用局部项救援调权。

新组合保留原366191参数C0，加原双头支路同样的零初始化Linear(2,14,bias=False)，总366219参数。旧C0组恢复AdamW step4608和全部moments，新28系数从step0及无moments起，最终分别9216和4608。原b048编码器和两个六位置BF16分类头均冻结，两个头都保留，不按TRAIN或DEV选择head。原每窗口一次4096六碱基token共同编码，BF16隐藏状态转FP32；短窗口定向后按同规则右补N到24576，只用真实有效碱基。native CDS减NON_CODING的两列FP32分数按原六位置顺序映射到十四非背景角色，背景role0不改。

先报告各物种、worst-species及pooled完整链P/R/F1，再看CDS base、背景、重复预测及unsupported层。若四格末权重和原六染色体raw共同评价齐备，分别报告每个支路的控制差值；可描述F11−F10−F01+F00的百分点交互量，但这是单seed、TRAIN、非线性F1尺度上的描述，不作显著性、机制因果或泛化结论。新增28参数和共同裁剪改变优化过程，匹配样本和累计更新次数不等于有效幅度完全匹配。主链不是UTR／所有isoform，GEN上游监督暴露限制仍有效；预训练分类头本身不是本项目创新。

## 实施与验证

新严格caller只导入原推理模块及固定模型定义，不导入训练计划、参考catalog或观察损失。复用JOINT24576真实窗口编码和全定向染色体skyline，神经接口复用原双native适配器。仅接受组合协议的累计9216／epoch6和原revision，拒绝旧6144先验、无先验JOINT24576及中间权重；输出每染色体128MiB上限，不保存特征或冻结head副本。

人工检查覆盖组合身份拒绝、4096 token产生完整24576碱基的两头六位置排列，以及1／14 bp短窗口零先验与原C0输出精确相同。用构造特征、构造神经权重及构造分类头，不读研究DNA或权重，不运行原生模型。

原生预检固定原序列中每物种第一个draw。每例独立恢复原末权重与AdamW；同一次真实joint原生forward的features上，零先验窗口势必须与原C0逐元素相同，结构NLL须与既有JOINT24576预检13453758的同draw一致。然后做可丢弃一次更新，旧组4608→4609、新组0→1，所有可训练参数及梯度有限，编码器和两个头无梯度。更新后的神经权重仅在内存序列化、strict加载后，使用同一raw适配器及owner拼接重算，势函数须精确相同。另用14 bp人工DNA输入执行原生双链raw推理核心；这不包含FASTA文件入口和末9216身份核验，完整严格FASTA路径仍需真实末权重生成后验证。两次更新不保存科学checkpoint，不带入正式fit，也不将更新后的神经权重冒充末9216权重。

预检使用generanno环境，1 RTX3090／2CPU／16GiB／10分钟；依据已有双头预检52秒和联合预检98秒留加载与接线余量，不由预检loss宣称性能。CPU人工检查2CPU／4GiB／5分钟、无GPU。输出分别为`outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-CHECKS`和`M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-PREFLIGHT`，后者同置outputs，不改活动源代码或原结果。实现失败保留产物，只允许一次同科学设置的明确修复。

预检通过不等于完整链收益。本阶段不改变原七权重TRAIN汇总、原八方法DEV或已事前固定的十方法扩展，不自动将组合加入任何DEV家族，不读取DEV、sealed、held-out或Setaria，不公开、提交或推送。长等待仍只复用原监控对话及唯一心跳，短检查由主处理，不新增监控。

## 正式执行链

独立新增`fit_full_grammar_joint24_native_prior.py`、CPU核验、严格FASTA TRAIN入口及原评价器的明确路径绑定，不修改任何活动控制组代码。正式fit只接受匹配13454690预检和本轮人工绑定检查的完成证据；从原末4608及原AdamW状态重新开始，不继承任何可丢弃更新。每遍1536固定draw，保存6144／7680／9216三次权重和旧／新参数组moments，正式只最后9216；上限96MiB，无自动恢复、延训或指标驱动重试。

CPU核验必须对4608实际训练日志逐行还原draw、species、fine／coarse／unknown、known分母、importance及损失计算，核对三遍均值和三个checkpoint的strict模型加载、有限FP32及old／new AdamW组步数与moment形状。核验成功后，六任务每次只从原TRAIN FASTA保留固定一条染色体，再经新严格末权重FASTA caller完整双链、原owner、全定向染色体skyline推理。评价必须等全部六任务成功并核实实际输入、末权重身份、窗口全覆盖及GFF3后，才读取同范围参考；直接复用原TRAIN评价器，保留22688主参考、unsupported、isoform、背景和原科学口径，不将训练eligible窗口过滤评价分母。

人工绑定检查只用构造元数据、构造4608行日志及mock评价器，检测旧6144先验或无先验JOINT身份误接、未完成预检／核验、错末权重路径、样本／importance／损失计算漂移和误用原评价器默认路径。检查不读研究DNA、参考或权重，不运行native模型，不进行更新；检查失败则不提交正式fit。

资源为private-teodoro-gpu、generanno：fit1 RTX3090／2CPU／16GiB／18小时；CPU核验2CPU／4GiB／5分钟；raw数组0–5、最多3同时，每项1 RTX3090／2CPU／16GiB／4小时；评价2CPU／4GiB／20分钟。fit沿用单独JOINT训练实际约10秒／draw的量级及加载／新先验余量。raw时间依据已完成LOCAL六染色体单项28:41—38:47，并为原生联合注意力差异留保守余量；不是按12小时分区上限分配。预期新增fit加六FASTA和每染色体128MiB预测上限约1GiB，不保存原生特征或head副本。实际150TiB可用空间足够，不替代输出上限。

依赖严格为fit→afterok CPU核验→afterok六raw→整数组afterok CPU评价，任一失败保留产物，停止自动下游；只允许协议内一次明确实现修复，不自动改参数、放行DEV、增训或延长作业。组合四格比较另做事前约定的独立TRAIN读出，不追加到已有七权重汇总或DEV家族。

## 实际执行证据

CPU人工检查13454611于2026-10-09 13:51:32—13:51:53 UTC完成，Slurm root／batch／extern均COMPLETED、0:0。三项测试全部通过，summary及STATUS齐备；未读研究DNA、模型权重或DEV／封存数据，未运行原生编码器或更新权重。该结果只允许继续原生接线预检，不代表正式组合模型通过科学评价。

原生预检命令为`sbatch --parsable sbatch/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-PREFLIGHT.sbatch 13454611`。实际作业13454690使用private-teodoro-gpu、gpu035、1 RTX3090／2CPU／16GiB，于13:58:00—13:59:17 UTC完成，Slurm root／batch／extern均COMPLETED、0:0，STATUS和summary齐备。原TRAIN ordinal0拟南芥和ordinal1水稻各独立更新一次；每例三个真实joint窗口的零先验logits逐元素等于原C0，起点结构NLL与原13453758相同。起点raw NLL分别24.060655879005026和0.0015818514566063868，只用于同输入接线回放，不用来判断组合质量。

两例旧AdamW组均4608→4609，新先验组0→1并产生非零系数，有限参数／梯度及共同clip检查通过。内存strict神经权重回放后，原生raw适配器及owner拼接势逐元素相同；14 bp人工DNA双链均保留全部有效碱基，零完整CDS预测只表示该人工输入执行完成。编码器和两个原生头保持冻结、无梯度；只解析使用原六TRAIN染色体，未使用缓存特征、DEV或封存数据。未写科学checkpoint、未启动正式fit，两个可丢弃更新不计入任何科学末权重。上述证据支持组合工程接线可执行，不支持互补性、完整链收益、泛化或论文创新成立。

新增绑定CPU检查13457836于14:07:36—14:08:14 UTC在gpu034完成，38秒，root／batch／extern均COMPLETED、0:0，五项检查全部通过。只用人工元数据、人工日志和mock评价器，未读DNA、参考、权重，零原生forward、零更新。未改活动控制代码或冻结评价规则。

正式执行实际提交如下；每次使用返回的真实作业号串联，不依据预想ID提交：

```sh
sbatch --parsable sbatch/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-FIT.sbatch 13454690 13457836
# 返回13457841
sbatch --parsable --dependency=afterok:13457841 sbatch/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-VERIFY.sbatch 13457841
# 返回13457846
sbatch --parsable --dependency=afterok:13457846 sbatch/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-TRAIN-RAW.sbatch 13457841
# 返回13457849
sbatch --parsable --dependency=afterok:13457849 sbatch/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-TRAIN-EVAL.sbatch
# 返回13457852
```

fit13457841于14:09:34 UTC在gpu035启动，1 RTX3090／2CPU／16GiB／18小时。主已核对实际RUNNING、STATUS／JOBID／started.json匹配组合配置、13454690及13457836，stderr无错误。实际调度显示核验13457846等待afterok:13457841，raw13457849的六项（max3）等待afterok:13457846，评价13457852等待afterok:13457849_*整数组；均是有效依赖，不是已完成结果。原五个未完fit仍继续，新增组合不重跑已有三控制。现在没有组合正式末权重、完整FASTA回放或科学指标可供解释。

原唯一心跳coding-rna-r5-light-monitor于14:12:44.358 UTC仅更新对象至五条长cohort，原监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3、主对话、gpt-5.6-luna／low、15分钟和安静监控设置均保持不变。完整保存提示22837字元、heartbeat／ACTIVE／原target／created_at回读匹配。实际原对话手动首检turn01a12103-48ec-7853-b7f7-8a2d5c0f65e5于14:13:45.919—14:14:17 UTC完成，确认当前五流程的六个fit均RUNNING，其余阶段有效Dependency；未读科学产物或回传／暂停。实际turn_context仍gpt-5.6-luna／low。该证据是更新配置及手动首检，不代替定时追加或最终自动回传。主仅负责结果科学核验及已授权的比较，不形成另一套定时轮询。

更新后的实际定时追加另已核实：原对话turn01a12105-c3ca-72f2-a322-a1458a7da0c7由同一automation于14:16:28.360 UTC触发，actual turn_context14:16:28.447仍gpt-5.6-luna／low；输入实际含新13457841／13457852和五cohort范围。调度取证包含全部28个当前作业root、exit0，14:16:45检查当前六fitRUNNING、下游有效Dependency，未读科学产物、未回传或暂停；14:16:54完成并DONT_NOTIFY。配置、手动首检、真正定时追加分别成立，最终完成或异常后的自动回传仍待观察。

## 四格独立 TRAIN 汇总接入

在读取任何四格正式末9216指标之前，实现独立`aggregate_full_grammar_context_prior_train.py`。F00固定STRUCTURED_CONTINUE／评价13451625；F01固定6144双头先验／13453184；F10固定24576无先验／13453780；F11固定24576双头先验／13457852。只读这四个已经完成的原同尺`summary.json`和STATUS／JOBID；任一缺失、未完成、作业号或最终权重错接、六染色体范围／22688主分母／参考元数据变化，均不生成比较。直接复用原TP／预测／参考计数计算与物种汇总，核对phase检查数量、CDS base、背景、isoform和unsupported元数据；本步不重新读取研究DNA、参考、预测GFF或checkpoint，也不重检参考ID／几何。

固定输出四格逐染色体及物种／pooled完整测量，分别计算F10−F00、F01−F00、F11−F01、F11−F10、F11−F00，以及F11−F10−F01+F00的完整链F1百分点描述量。pooled先汇总计数再算F1，macro为两个物种的算术平均；worst先在每格取物种最小值，再计算同样对比，不是取各物种差值的最小值。最差物种切换必须保留这一解释。只记录单seed TRAIN非线性F1尺度上的描述，不排序选模、不生成gene级伪独立区间，不作为互补机制因果、显著性、泛化或创新成立证据。

人工CPU检查覆盖六种预定对比的独立计数期望值、物种最小值切换、全部四格及六染色体完整性、错job／错末权重、参考元数据漂移和聚合算术错误；不读取研究数据、权重或运行native模型。检查及正式汇总均2CPU／4GiB／5分钟、无GPU，实际数据仅小型评价JSON，预计低于1MiB。检查通过后才提交正式汇总，afterok同时绑定上述四个评价实际ID，不重提任何fit或raw。输出分别为`outputs/M28-FULL-GRAMMAR-CONTEXT-PRIOR-TRAIN-CHECKS`、`outputs/M28-FULL-GRAMMAR-CONTEXT-PRIOR-TRAIN-COMPARISON`。不改七权重TRAIN、八／十方法DEV或活动训练源文件。正式汇总加入原第五cohort的末端监控，四个控制未齐造成的有效Dependency保持安静，不新建监控对话或自动化。

人工检查实际13460966于2026-10-09 14:27:55—14:28:31 UTC在gpu034完成，36秒，Slurm root／batch／extern均COMPLETED、0:0；STATUS／JOBID／summary齐备，五项测试全部通过。六种差值的独立计数期望、最差物种切换、错末权重／作业和参考／聚合漂移均检查通过；未读取研究DNA、权重或DEV／封存数据，零native forward、零更新。

正式汇总实际提交：`sbatch --parsable --dependency=afterok:13451625:13453184:13453780:13457852 sbatch/M28-FULL-GRAMMAR-CONTEXT-PRIOR-TRAIN-COMPARISON.sbatch`，返回13460976。它只等待四个已有评价，不新增fit／raw，不改现有作业参数。现在只有人工检查完成和正式汇总依赖提交证据，没有四格科学比较结果。

主已核对实际13460976为PENDING／Dependency，四个afterok均尚未满足；private-teodoro-gpu、2CPU／4GiB／5分钟、无GPU、Requeue=0符合约定。原唯一心跳于14:30:05.198 UTC仅增加这一末端汇总，五个cohort及原ID均保留，第五条完成条件必须含13460976；评价13457852完成但汇总仍等其他控制时保持安静。保存的24631字元提示、heartbeat／ACTIVE／原target／created_at／15分钟／通知设置回读完全匹配，未新建对话、自动化或改模型。此前14:16定时证据不代替本次新增末端的运行验证；本次实际定时追加和最终自动回传需分别从真实事件核实。

本次新增末端后的实际定时运行已观察：原监控turn01a12113-7f8e-7650-bac0-2df222de7d12于14:31:28—14:31:52 UTC完成，actual turn_context14:31:28.460仍gpt-5.6-luna／low。14:31:43定时检查包含13460976，汇总为PENDING有效Dependency，当前六fit仍RUNNING，其余未完阶段依赖正常；未读科学产物、未回传或暂停，DONT_NOTIFY。这是同一原对话的真实定时追加，不另发重复手动检查。新配置保存和新增末端定时检查均已验证，正式四格结果及最终自动回传仍未完成。
