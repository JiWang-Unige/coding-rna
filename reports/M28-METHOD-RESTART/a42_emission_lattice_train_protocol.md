# 冻结发射上的完整语法搜索

2026-10-07北京时间，在任何A42真实DNA搜索前固定。比较对象为原A31_A3，唯一机制改变是以A42完整受限语法窗口DP替代query/beam产生的有限自然池选择。原A3的SharedContext与十五类发射冻结；没有A41修正、新训练、GEN编码器前向、阈值或角色权重变化。1491窗发射13422268及32片搜索13422269已闭合；唯一汇总修复R1的13423713/13423721均完成，七项CPU证明和直接FP64算术通过，但五项科学条件全部FAIL，A42机制关闭、不救援。正确配对2009降至981，至少5-exon正确700降至81；这是原TRAIN窗口结果，不是raw注释、独立泛化或论文成功。

## 研究问题和范围

原A34诊断表明缺失参考链中存在评分非正、正分未选及可以选中的不同组，不能从搜索去剪枝推断性能改善。这里检验相同局部势、相同语法和gain0下，放开query/beam几何截断后自然已选完整链的变化。精确语法DP本身不是创新；冻结A3有既有监督曝光，不能称fresh模型或与外部方法的公平曝光因果消融。

仅原六条TRAIN染色体：拟南芥NC_003070.9、NC_003071.7、NC_003076.8；水稻NC_089043.1、NC_089044.1、NC_089046.1，两个方向。窗口支持固定为原1536draw首次出现顺序的1491唯一窗，不扩至25512全TRAIN窗，不选择536missing或按新gain挑窗。所有cached输入只读取原TRAIN index对应的artifact；已有缓存中的R2/R3标签字段不作为目标或选择输入。

原自然池16321条仅用于冻结势及较强A31基线重放，不注入A42、不限制新搜索。原baseline jobs为13398521和13398524，源checkpoint为outputs/M28-A3-FIRST-FIT/step_004608.pt，649457参数全部冻结。实际只调用shared与segmentation，valid长内保持FP32十五类logits；不运行query、pointer、beam、GEN或optimizer，不写checkpoint。缓存模型及revision沿用原验证函数，误配即停止，不另建兼容路径。

## 固定目标与支持

沿用原WIDTH24576、STRIDE12288、valid长度和oriented半开坐标，max_exons128，contiguous ATG/末STOP、GT/GC至AG、min_intron5、实际codon carry、无模糊CDS和internal STOP；内含子N仍允许。十五角色减intergenic以FP64累计，全部单位权重，null0、gain0、span非重叠且相触允许。没有score beam、site top-k、完成链数量上限、长度归一、角色系数、校准、温度或阈值扫描。全零null及确定性首次到达tie规则沿用已检查A42，不强求非零平分与旧有限池身份相同。

每次原A31候选的本次冻结分数须与旧值以原rtol1e-12、atol1e-8相符，并精确重现旧selected索引。原池全部链须在上述A42语法内；若不成立，停止同支持比较，不删链或扩大语法救援。新已选链用实际DNA拼接独立核验完整ORF、motif、cap及span，再用candidate_roles直接FP64求和对照A42保存的每链和总gain，沿用同容差。新最优总势不能低于同窗原池最优（在原数值容差内）；提高总势不等于增加正确链。

A42窗口返回最优gene set，不返回全部合法候选；所有合法边在标签连接前已经参与推理，不能按参考筛掉unknown或known_nonmatch再解码，也不宣称已选链数是总候选空间大小。

## 八窗工程预检

案例选择仅使用原free_pool中的species、strand和原1491唯一窗顺序：species按名称排序，每物种依次正负链，各取首两个，case0至7按该顺序。无需TRAIN参考、标签数组、missing分类或新势；案例ID在GPU产物中固定，CPU直接重放，不替换耗时或失败窗口。

preflight_emission_lattice_train.py的potentials阶段核对完整原人口元数据及原A31完成证据，仅为这八个窗口加载已有TRAIN缓存、前向冻结source并保存FP32发射、实际oriented DNA及原池/分数/索引。全部冻结state前后逐tensor相同、requires_grad及grad检查必须通过。预检只重放原baseline，不读取A41拟合或其读出；每个新文件关闭后才写GPU完成summary。

decode阶段每case在CPU只读已完成GPU的对应payload及原人工准备job完成状态，运行一次A42 max128，记录decode秒数、进程peak RSS、实际motif计数、遍历计数及已选链/势，完成独立ORF、角色算术和同窗baseline重算。没有参考标签、TRAIN catalog或新R4计算；八窗不计算TP、FP、F1、五条件或机制成功。这是资源和数值接线，不是随机抽样benchmark、最坏复杂度或完整TRAIN结果。

唯一GPU命令sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-TRAIN-POTENTIAL-PREFLIGHT.sbatch；generanno、private、1RTX3090至少24GiB、2CPU、16GiB、5分钟、no-requeue。根outputs/M28-A42-EMISSION-LATTICE-TRAIN-POTENTIAL-PREFLIGHT，仅8份局部势/DNA/原pool payload、summary、STATUS与JOBID；上限约12MiB加小型文字，需提交前核对存储。预计前向与写出低于原A41五分钟预检资源，不把旧八窗时序当本轮保证。

CPU命令为sbatch --parsable --dependency=afterok:GPU_JOB --export=ALL,POTENTIAL_JOB=GPU_JOB sbatch/M28-A42-EMISSION-LATTICE-TRAIN-DECODE-PREFLIGHT.sbatch；固定array0-7且最多并行8，generanno、private、0GPU、2CPU、8GiB、每task3分钟、no-requeue。根outputs/M28-A42-EMISSION-LATTICE-TRAIN-DECODE-PREFLIGHT/case_0至case_7，独立STATUS/JOBID/started.json/summary.json；唯一%x_%A_%a日志。人工密集窗实际115.38秒，3分钟留操作余量，不代表真实DNA保证。开始marker只记录已固定案例身份，不读取标签；summary只有正常完整计算后写出。

提交前执行五项活调度、资源及预约检查，核对精确输出/日志冲突和源文本，不改变A41或原结果。任一数值/语法/身份失败保留并报告；只有明确实现或fixture错误可独立R1修复一次，科学支持和资源不变。再次失败、TIMEOUT或OOM关闭本预检版本，不延时、续算、跳窗或改阈值。通过只支持准备完整1491窗资源接线。

## 完整TRAIN比较的冻结规则

完整运行复用原1491窗、同一源模型、上述目标及全部原较强baseline。八task所有实际产物与Slurm核对已完成；仅据实际资源固定32分片、时限、存储及下文提交命令，不改变窗支持或科学参数。不能把本协议或八窗通过当作全量已运行。

全部1491新自然推理文件及局部势关闭后，才连接原R4 TRAIN target与candidate_targets，保持primary、protected、unknown spans、callable全span及既有eligible定义。对新几何须重新调用原冻结标签函数，不能把旧索引标签套到新链；这是同一R4的标签配对，不重建或改变R4政策。原A31池重新连接后须逐值等于旧标签数组和计数，原2916 eligible reference/window配对及1086至少5-exon配对分母保持。不把它们当独立基因；unknown继续在推理竞争中保留，不作失败负例，不因未知数量增长删链。

原五条件不降低：拟南芥selected_correct严格高于1531；水稻严格高于478；selected_known_nonmatch分别不高于347和251；pooled至少5-exon正确严格高于700。任一FAIL关闭本机制版本，不救角色权重、阈值、loss、seed、epoch或checkpoint。完整窗比较需CPU独立重算保存势、原baseline、标签、分母、选择及五条件；工程闭合与科学PASS分开。PASS只支持另固定raw协议，不证明泛化、外部SOTA、因果机制或论文成功。

## raw边界与交付

owner人工13421848三检查已证明：selected-only窗输出再过滤owner不成立；owner过滤后仅保留局部最优也不足以跨窗竞争。这里明确只比较窗口机制，不把它改名为raw annotation。未来完整raw必须单独固定必要备选或联合decode的上下文、评分及owner定义，验证原完整参考分母与独立数据约束；不会在本预检中隐式改变原assembly。

当前交付推进完整1491窗比较；raw-FASTA注释、DEV/独立物种/Setaria、公平外部对照、论文和公开发布均未完成且不由本比较自动放行。A41完整链已闭合并科学FAIL关闭；短接线检查主对话核对，完整搜索等待只复用同一固定监控，不新建或重复轮询。

## 实际GPU发射与CPU提交

GPU13421953于2026-10-07 00:18:34 UTC只提交一次，00:18:35至00:19:03 UTC运行，28秒；root、batch、extern均COMPLETED0:0，STATUS/JOBID匹配。实际private、gpu034、1RTX3090、2CPU、16GiB、5分钟上限，Requeue0/Restarts0，canonical Command、SubmitLine、WorkDir及%x_%j日志相符；batch MaxRSS1468704K。提交前五相关文件逐字同步、Python AST及两SBATCH语法通过，五项活调度、节点与预约无硬冲突，两个根和日志均无冲突，scratch可用157T是共享余量而非独占保证。

八case原ordinal依次3、10、0、8、2、5、1、4，对应原协议species/strand/首两窗顺序，全部valid24576。原候选数7、14、12、9、10、11、11、12，共86；原A31分数最大误差全0、已选索引全部精确相同，86链全部通过实际DNA的固定语法核验。完整人口元数据1491/1536/16321相符。冻结source前后逐tensor相同、无requires_grad及grad，optimizer0、无checkpoint、query/GEN、新候选搜索或标签访问。实际feature revision为b0483c23b6b63787b61a6d3a204a9b517d6ba345。

八份payload共12017120字节，summary7378字节，STATUS/JOBID19字节，GPU stderr空；内部前向及保存共5.537628秒、GPU峰值allocated108832256字节。Slurm采样MaxRSS、GPU allocated和内部时序含义不同，八个非随机窗不作公平速度或全流程时限保证。这里实际读取研究checkpoint、原TRAIN缓存和DNA，不再是source-only人工检查；但未读取参考标签、TRAIN catalog、A41新checkpoint、DEV、独立物种、Setaria或A33。

CPUarray13421967于00:19:51 UTC仅一次提交，afterok13421953、export POTENTIAL_JOB13421953；0至7八task于00:19:52 UTC均在gpu034开始，throttle8。实际每task0GPU、2CPU、8GiB、3分钟，原命令、%x_%A_%a唯一日志（实际job名M28A42DECODEPRE）以及Requeue0/Restarts0相符。提交前再次五项活调度无硬冲突，GPU完整产物已闭合且CPU根不存在。启动时八STATUS均RUNNING，不把开始或started marker记为解码PASS；等待全部task、summary及直接算术核对后解释。

## 八窗实际闭合

13421967_0至_7的root、batch、extern全部COMPLETED0:0，八STATUS为COMPLETED、JOBID逐项为13421967_i；started及summary均非空，八stderr全空。共同起点2026-10-07 00:19:52 UTC，最后case3于00:22:19 UTC结束，没有失败、TIMEOUT、OOM、R1、重跑、续算或资源延长。各task运行1:56、2:14、1:56、2:27、2:04、2:19、2:15、2:06；batch MaxRSS依次438104、455636、447688、633708、463680、474184、489884、605940K，Slurm采样与下面进程峰值分开。

| case / 原ordinal | 物种 / 方向 | decode秒 | 进程peak RSS KiB | 已选链数 |
| --- | --- | ---: | ---: | ---: |
| 0 / 3 | 拟南芥 + | 97.314821 | 504944 | 0 |
| 1 / 10 | 拟南芥 + | 116.552106 | 519932 | 6 |
| 2 / 0 | 拟南芥 - | 97.757588 | 512576 | 2 |
| 3 / 8 | 拟南芥 - | 127.361953 | 720072 | 0 |
| 4 / 2 | 水稻 + | 105.672909 | 527808 | 3 |
| 5 / 5 | 水稻 + | 121.591410 | 540540 | 3 |
| 6 / 1 | 水稻 - | 117.773038 | 544688 | 0 |
| 7 / 4 | 水稻 - | 106.182625 | 708032 | 1 |

八case均与GPU案例、原ordinal、window ID、configuration、protocol、source job13421953及feature revision逐项相符。CPU原A31分数最大误差全0，原selected索引精确相同，全部原池链语法支持。新选择总15条，实际DNA拼接ORF、motif、cap、span非重叠及每链/总FP64角色势核验均通过；三个null case的空选择算术也闭合，但不把空链检查说成非空ORF证据。所有case最优势不低于原池，case1、2、5、7严格增大，另外四个相等；这只是目标函数变化，没有标签，不能说四窗性能提高。

真实窗exon_extensions从19506278至39716819，最大active状态均2667；三个null也实际遍历完整语法，没有零势或空答案快捷退出。case3虽返回null仍解码127.361953秒，说明不能用已选链数预测耗时。八个同节点并行、非随机样本支持本预检可完成，不是公平速度、全基因组耗时或最坏复杂度保证。

本阶段工程结论为PASS：真实DNA上的冻结势重放、语法支持、数值接线及有限资源运行成立。optimizer0，没有checkpoint、模型/GEN CPU前向、参考标签连接、TRAIN catalog、新R4、DEV、独立物种、Setaria或A33。没有计算原五条件，不能写TRAIN科学PASS、泛化或完整annotation；owner备选问题仍未解决，原assembly与A42解码源码未改。

据此继续准备完整1491窗的同势比较，不再重复八窗或已闭合人工检查。资源初步采用一次GPU发射30分钟、32个CPU shard（ordinal模32，最多每片47窗）及一次CPU结果核验；按本预检最大127.36秒外推每片约100分钟、再留约25%和启动/IO余量，候选时限2:15。这仅是操作资源依据，非运行保证；完整入口、精确资源、产物及命令须在提交前写入本协议，不据八窗gain改变科学配置，不把待实施写成已提交。

## 完整1491窗执行合同

run_emission_lattice_train.py分为preparation、potentials、decode和readout四阶段，verify_emission_lattice_train.py做最后独立CPU重算。preparation先执行五项新adapter人工检查：新几何重新配对原R4且unknown保留；原政策/分母漂移或陈旧选择拒绝；分片乱序恢复及缺失/重复拒绝；null与FP64直接角色算术；五项严格条件。随后仅用原JSON核对1491/1536/16321、原A31标签/选择及2916/1086分母，A42占位不作结果；不读取权重、特征tensor或TRAIN catalog。只新增这五项接线检查，不重跑已闭合的13项解码检查或八窗预检。

potentials加载原冻结source，保存window_000000.pt至window_001490.pt各窗FP32十五类发射、实际oriented DNA、原池/分数/选择及身份；每窗重放原分数、索引和全部原链语法。逐tensor冻结state前后相等、grad及requires_grad为零，全部文件关闭后写summary；不读取新标签、TRAIN catalog或A41权重。decode固定array0-31且最多并发32，原ordinal模32，每片46或47窗，顺序不改，每窗一次完整max128搜索和直接ORF/角色/原目标核验；不跳失败或耗时窗、不作局部续算。各片started.json、natural_scores.jsonl及summary分开，任何未闭合片阻止下游。

readout在全部32片COMPLETED0:0、1491行及身份闭合后，先关闭按原ordinal合并的natural_scores.jsonl，才打开原TRAIN catalog及旧保存标签，调用同一build_targets/candidate_targets配对新几何。旧池标签逐值/计数及全部原分母相等；不按旧索引给新链贴标签。windows.jsonl保留原及新标签，summary只报告两臂已选计数，不把已选集合冒充全部合法池。显式记录new_geometry_labels_computed=true和R4_policy_or_catalog_changed=false；不将重新计算新几何说成没有新标签计算。

CPU核验重新读取实际1491潜势及32片完整身份，从全窗candidate_roles角色赋值独立直接求FP64原池和新链势、原DP、新链DNA ORF/motif/cap/span、总势不低于原池。独立计数代码重算冻结R4标签、unknown保护、分母、物种聚合和五条件，核对保存结果。七证明覆盖完整原人口、冻结源保存势及零更新收据、原A31分数/DP/基线、新链语法/角色算术、32片先搜索后标签、原R4标签/分母和五条件。不会重跑1491遍新DP或独立神经前向；冻结state证明引用GPU实际收据，不冒充CPU重新前向验证。

所有阶段generanno、private、no-requeue、2CPU、offline；仅发射阶段1RTX3090至少24GiB/16GiB/30分钟，其余0GPU/8GiB：CHECKS5分钟、每decode片2:15、readout15分钟、result-check15分钟。潜势上限约2.20GB（24576×15×4×1491）加约37MB DNA及原池metadata，预计完整新增低于2.4GB；提交前核对实际存储。2:15按八窗最大127.36秒×47约100分钟再留操作余量，并不保证全部窗复杂度。任一身份/数值/语法失败保留并报告；仅明确实现/fixture错误可一次独立R1、另立输出根及下游作业、不改科学参数和原资源。再次实现失败、TIMEOUT、OOM或科学FAIL关闭本版本，不续算、延长、跳窗、重调或混用旧输出。

固定命令依次为sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-TRAIN-CHECKS.sbatch；准备闭合后sbatch --parsable --export=ALL,CHECKS_JOB=CHECKS_JOB sbatch/M28-A42-EMISSION-LATTICE-TRAIN-POTENTIALS.sbatch；sbatch --parsable --dependency=afterok:POTENTIAL_JOB --export=ALL,POTENTIAL_JOB=POTENTIAL_JOB sbatch/M28-A42-EMISSION-LATTICE-TRAIN-DECODE.sbatch；sbatch --parsable --dependency=afterok:DECODE_JOB --export=ALL,POTENTIAL_JOB=POTENTIAL_JOB,DECODE_JOB=DECODE_JOB sbatch/M28-A42-EMISSION-LATTICE-TRAIN-READOUT.sbatch；sbatch --parsable --dependency=afterok:READOUT_JOB --export=ALL,POTENTIAL_JOB=POTENTIAL_JOB,DECODE_JOB=DECODE_JOB,READOUT_JOB=READOUT_JOB sbatch/M28-A42-EMISSION-LATTICE-RESULT-CHECK.sbatch。大写JOB是实际已取证作业ID占位，不由文件或协议存在推定提交完成。

独立根依次outputs/M28-A42-EMISSION-LATTICE-TRAIN-CHECKS、TRAIN-POTENTIALS、TRAIN-DECODE/shard_00至31、TRAIN-READOUT及M28-A42-EMISSION-LATTICE-RESULT-CHECK（前四省略共同M28-A42-EMISSION-LATTICE前缀）。各根STATUS/JOBID及非空summary必需；GPU另需1491潜势，decode每片started及46/47自然行，readout自然及label文件各1491行，CPU七证明均true。日志%x_%j或array%x_%A_%a唯一，不覆盖原根。每次提交前五项活调度、预约/节点/存储、精确根和日志冲突检查；source逐字同步及AST/SBATCH语法核对。完整链只复用固定监控01a0f71f-8845-7432-bd13-bf54c1808fb3及同一heartbeat，不读取中途gain/标签/科学条件。

## 实际接线与完整链提交

13422248于2026-10-07 00:44:54至00:45:22 UTC完成，root、batch、extern全部COMPLETED0:0，28秒、private gpu034/0GPU/2CPU/8GiB/5分钟、batch MaxRSS787220K。五项新人工adapter检查PASS，0.029秒；原1491窗、1536draw、16321自然池、较强A31及2172/854拟南芥和744/232水稻reference/window分母相符，总2916/1086。STATUS/JOBID及实际Command/SubmitLine/WorkDir、唯一日志、Requeue0/Restarts0匹配；没有特征tensor、模型前向、TRAIN catalog或新几何标签计算。

完整链仅一次提交：potentials13422268与decodearray13422269于00:46:54 UTC提交；readout13422270及CPU13422271于00:46:55 UTC提交。GPU于00:46:55 UTC在gpu034实际开始，1RTX3090/2CPU16GiB/30分钟；array0-31%32、afterok13422268及POTENTIAL_JOB相符，每片0GPU/2CPU8GiB/2:15；readout实际依赖afterok13422269_*，覆盖整array而非单task，0GPU2CPU8GiB/15分钟；CPU依赖afterok13422270，同CPU15分钟。三个下游PENDING/Dependency正常，不能据提交或启动写计算完成或科学PASS。所有canonical Command、SubmitLine、WorkDir、exports、日志和no-requeue/Restarts0已核对；array pending-parent的4294967294日志显示是Slurm占位，实际task采用原%A_%a，不作真实日志文件名。

两次提交前各五项活调度及预约无相关硬冲突，private gpu[034-035]均24GiB RTX3090、有效CPU与内存足够；MaxTime7天覆盖原资源，未改变分区/资源/科学参数。全部新根和精确日志无冲突，scratch可用157T是共享余量，不是独占保证。三个Python AST、五SBATCH语法与全部十个新增或相关记录文件逐字local/canonical一致后运行；原M20未提交改动保留，无Git提交或公开发布。

固定监控沿用原对话01a0f71f-8845-7432-bd13-bf54c1808fb3及coding-rna-r5-light-monitor，原A41终态交接已消费，重新指向本四阶段35个root。原生ACTIVE更新、完整新prompt和ID/name/kindheartbeat/target/15分钟/failed_runs_only/created_at读回一致，未创建/fork监控或用cron替代。一次manual检查于2026-10-07 00:50:00 UTC核对GPU RUNNING及下游正常Dependency，不读产物或科学值、不向主交接；同线程turn01a113d6-3bf3-7522-8488-a985c527baf4 completed/errornull/43.480秒，实际turn_context保留gpt-5.6-luna/low。配置和单次手动检查与以下本轮定时证据分别成立，不继承A41记录为本轮成功；完整终态回传仍未发生。主不重复等待轮询，终态由同一轻量监控交接后再按原五条件解释。

本轮首次定时执行已验证：2026-10-07 01:03:47.271 UTC实际heartbeat user trigger包含同automation ID及current_time_iso01:03:47.057，追加到固定监控对话turn01a113e3-5173-7493-8ec5-5aa826f7142f；该turn completed/errornull/20.796秒，实际turn_context仍gpt-5.6-luna/low/原cwd。01:03:58 UTC最小调度快照为potentials13422268 COMPLETED0:0、32个decode task全部RUNNING、readout及CPU正常等待依赖；不读取潜势、gain、标签或条件，也未向主发送终态交接。此为同一对话的定时检查成功，不是GPU完整产物核验、完整1491搜索闭合或科学PASS。

## 完整搜索闭合及汇总修复

13422268实际于2026-10-07 00:46:55至01:02:53 UTC完成、15分58秒。32个13422269_i均于01:02:54开始，最早task29于02:32:22、最后task16于02:47:43结束，全部root COMPLETED0:0，各片1:29:28至1:44:49，未超原2:15资源上限。原readout13422270于02:47:44开始、02:48:21失败，37秒、root/batch FAILED1:0而extern COMPLETED0:0，batch MaxRSS760600K；STATUS为FAILED exit=1，JOBID相符，合并natural_scores有1491行而windows为0行。原CPU13422271保持PENDING/DependencyNeverSatisfied，未运行、未取消；不能说CPU核验失败或已有完整科学结论。

同一固定监控于02:48:58 UTC取证并向主交接一次，随后暂停同heartbeat，turn01a11443-7560-7fe0-9211-3d2bc4179b47完成/errornull。具体异常为run_emission_lattice_train.py readout的Original TRAIN target and natural window geometry differ；尚未进入首窗build_targets。首窗TRAIN自然记录有原GEOMETRY十字段，R4窗口另有chain_unknown_intervals、fine_label_unknown_intervals、label_policy、motifs、no_annotated_gene_either_strand、overlapping_transcript_ids、positive_ids、primary_complete_ids八字段；所有共同字段逐值相同。直接比较整字典误报，不能据此宣称真实坐标漂移；其余1490个窗口仍须R1逐窗完整验证，不用首窗事实替代全量。

唯一R1仅改readout的目标记录比较：按项目既有GEOMETRY十字段与自然窗口逐值相等，潜势与自然窗口仍严格整字典相等；R4全字段记录原样传入build_targets。独立CPU同样对齐十字段，但保存和比较的label row仍用原自然窗口几何，不能把R4标签元数据改写到自然预测记录。原冻结R4配对、unknown保护、旧标签逐值/旧baseline/2916及1086分母、评分容差、五条件和七证明全部保留，不改科学协议。新增实际形状人工回归证明附加标签字段不造成错配、目标数据不被修改，并分别扰动全部十个几何字段均须拒绝；readout R1开始先运行这项和原五项adapter检查。

仅重新提交两项下游：sbatch --parsable --export=ALL,POTENTIAL_JOB=13422268,DECODE_JOB=13422269 sbatch/M28-A42-EMISSION-LATTICE-TRAIN-READOUT-R1.sbatch；闭合CPU命令为sbatch --parsable --dependency=afterok:READOUT_R1_JOB --export=ALL,POTENTIAL_JOB=13422268,DECODE_JOB=13422269,READOUT_JOB=READOUT_R1_JOB sbatch/M28-A42-EMISSION-LATTICE-RESULT-CHECK-R1.sbatch。原32片已全部COMPLETED0:0，提交前实际终态与STATUS/JOBID/summary齐，R1入口仍以sacct及原闭合证明检查全部源；不依赖已移出调度器的旧array。沿用generanno/private/每项0GPU2CPU8GiB15分钟/no-requeue，独立根outputs/M28-A42-EMISSION-LATTICE-TRAIN-READOUT-R1和outputs/M28-A42-EMISSION-LATTICE-RESULT-CHECK-R1及新job名%x_%j日志。两个SBATCH只显式将现有模块的输出常量指向R1根，不新增搜索、神经前向或原发射/解码重跑，不覆盖旧根。

提交前相关源码定向差异/AST/SBATCH、五项活调度、节点/预约和精确输出/日志冲突须通过；源发射和32片在R1入口再次验证全部原身份/配置/自然闭合证明及1491排序，不混入人工缓存答案。再次实现失败、TIMEOUT或OOM即关闭，无R2、延长、续算、跳窗或门槛救援。若CPU完整闭合后科学条件FAIL，同样关闭机制版本；结果处理范围只含原TRAIN，没有raw/DEV/独立物种/Setaria/A33或对外发布。长等待继续复用唯一固定监控及原模型，切换到新R1两项下游和原33个闭合源，原失效依赖不再监控。

2026-10-07 02:56:57 UTC首次R1提交尝试指定afterok:13422269，被sbatch以Job dependency problem拒绝，没有生成作业、日志或R1根；sacct按新job名未找到作业，旧array的scontrol已不能查询。去掉已结束且无法接受的旧调度依赖，保留实际源终态及程序内所有源检查、原科学命令/资源/输出根；这是提交前调度拒绝纠正，不是第二次运行或R2。

R1实际readout13423713于2026-10-07 02:58:11 UTC提交并在gpu034开始；原自然搜索源依次13422268/13422269，仅一次实际运行。六项adapter检查PASS0.091秒，包括附加R4标签元数据和全部十字段拒绝错配；当前RUNNING不作完整读出或科学结论。CPU13423721于02:58:51提交、afterok13423713，同原源export和READOUT_JOB13423713；尚需完整结果核验。两次提交前各五项活调度/节点/预约和精确根/日志冲突通过，七个相关文件local/canonical逐字一致、Python AST与两SBATCH语法通过；原源markers完整，scratch可用156T为共享余量，新自然JSON约1.8MB加标签和summary，不重新写2.2GB发射。全部原失败根及旧未满足CPU依赖保留，不取消其他作业、Git提交或公开发布。

原生将同heartbeat更新为R1 ACTIVE，完整新prompt、id/name/kind heartbeat、固定target、15分钟、failed_runs_only与created_at逐值读回一致，没有新对话、fork、cron或模型覆盖。两个实际Command/SubmitLine/WorkDir/exports/新日志与Requeue0/Restarts0相符，CPU登记afterok13423713(unfulfilled)正常。同监控03:00:41 UTC手动快照为readout RUNNING/gpu034、CPU PENDING/正常Dependency，无产物或科学值读取；实际turn01a1144e-1fdf-7872-ba3f-d5266020d3f1的03:00:26.763 UTC context仍gpt-5.6-luna/low/原cwd。当前仅证明配置与手动检查成功；本R1定时检查及整链自动终态回传尚未观察，不能继承原失败交接作为R1完成。

## 保存结果的链竞争诊断

diagnose_emission_lattice_competition.py只解释已闭合A42与原A31_A3的窗口选择差异，不重训或救援科学FAIL。先要求实际readout13423713/CPU13423721/原baseline13398521及13398524完成，核对原配置、source_jobs、1491/1536/16321、七证明、原2172/854及744/232分母、两份汇总计数和原五条件算术；PASS或FAIL均可描述，未闭合或缺证明禁止打开预测行。随后仅读取R1的natural_scores/windows及原A31 natural_scores三个JSON，按原ordinal/完整窗口身份对齐，不开局部势tensor、研究checkpoint、DNA、TRAIN catalog或参考，不重新配R4标签。

新池与旧有限池不同，必须以完整exon几何连接而不是复用候选索引。同一物种内分别核对correct、known_nonmatch和unknown的retained/lost/gained，正确至少5-exon另列；每类retained加lost精确等于旧已选数、retained加gained等于新已选数，全部物种/完整分母与闭合汇总相等。共同几何的原FP64分数和标签须一致，保持原数值容差。统计单位是原TRAIN窗口—链配对，不是独立基因；unknown不能视为负例或生物假阳性，新增correct也不等于发现新基因。

另核对已选重叠链的必要替换条件：旧已选正分链已通过A42语法，若该链在新集合丢失，则删除新已选中所有与其完整span有真重叠的链，再加入旧链，仍是允许的非重叠gene set。因此重叠链的gain总和不能低于旧链gain超过原数值容差，否则保存答案存在具体可改进的目标，不把问题直接归给模型。相触不算重叠；同势平分不强求旧有限池身份。此为局部必要条件，不是独立全空间最优证明。正确丢失链另记重叠竞争者的标签集合，每旧链一个组合，多旧链共享竞争者不计作独立竞争者或可相加全局gain贡献。按窗口记录目标增大或数值相等与正确数增/平/减的交叉，不能从目标提高推断性能提高，也不宣称训练loss因果。

四项新的人工保存行检查覆盖不同候选池的几何连接与两物种净计数、unknown/长链/联合标签组合、同分和相触、缺失/乱序/共同几何标签分数漂移，以及总势已高于baseline却有单条正分链可改进的实际失败。逐一缺少七项CPU证明时，真实行读取必须先被阻止；不重复原DP、R4或六项已闭合adapter检查。检查唯一命令sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-COMPETITION-CHECKS.sbatch，generanno/private/0GPU2CPU8GiB2分钟/no-requeue，独立根outputs/M28-A42-EMISSION-LATTICE-COMPETITION-CHECKS及M28A42COMPCHECKS_%j日志。

四项人工检查和真实R1整链闭合之后，才可提交sbatch --parsable sbatch/M28-A42-EMISSION-LATTICE-COMPETITION.sbatch，generanno/private/0GPU2CPU8GiB5分钟/no-requeue，独立根outputs/M28-A42-EMISSION-LATTICE-COMPETITION及M28A42COMP_%j日志，小型summary/STATUS/JOBID；不是放行raw或DEV。每次提交前source/AST/SBATCH、五项活调度、节点/预约及精确冲突核对。真实身份/计数/替换条件失败保留并停止解释；仅明确诊断实现或fixture错误可同输入/资源另根R1一次，重试失败、TIMEOUT/OOM关闭诊断，不修复或重启已关闭A42、不改评分/标签/数据/五条件/资源或覆盖旧根。

13423810于2026-10-07 03:11:26 UTC提交、03:11:27开始、03:11:39结束，12秒、root/batch/extern全部COMPLETED0:0，STATUS/JOBID相符。四项人工检查PASS0.010秒；实际private/gpu034/0GPU2CPU8GiB2分钟/generanno/no-requeue，batch MaxRSS750580K，Command/SubmitLine/WorkDir/log及Requeue0/Restarts0相符。五个新增或相关文件逐字local/canonical一致，Python AST与两SBATCH语法、五项活调度/节点/预约/精确冲突核对通过；无真实行、局部势/模型tensor、标签重配或新搜索。该PASS仅放行上述已闭合保存JSON诊断，不改A42科学结论。

## 完整TRAIN实际结果与关闭决定

R1 readout13423713于2026-10-07 02:58:11至03:04:00 UTC完成，5分49秒、batch MaxRSS846752K；CPU13423721于03:04:01至03:07:47完成，3分46秒、batch MaxRSS846704K。两项root/batch/extern均COMPLETED0:0、STATUS/JOBID正确，原private/gpu034/0GPU2CPU8GiB15分钟/no-requeue/Restarts0及canonical命令、exports、日志不变。两份summary的原配置/源身份/计数/完整分母/五条件一致，两1491行文件非空；全部七证明true，max_original_A31_score_error和max_new_gain_error均0。没有R2、资源延长、跳窗、推理重跑、模型更新或封存数据访问。

| 原TRAIN窗口配对 | A31_A3正确 | A42正确 | A31_A3 known_nonmatch | A42 known_nonmatch | A31_A3 unknown | A42 unknown | A31_A3至少5-exon正确 | A42至少5-exon正确 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | 1531 | 748 | 347 | 1459 | 315 | 407 | 548 | 72 |
| 水稻 | 478 | 233 | 251 | 618 | 132 | 165 | 152 | 9 |
| pooled | 2009 | 981 | 598 | 2077 | 447 | 572 | 700 | 81 |

原eligible配对仍为拟南芥2172/至少5-exon854，水稻744/232，合计2916/1086；不是独立基因分母。原五条件逐项均false（0/5）：两个物种正确未严格增加，known_nonmatch均增加，pooled长链正确未严格增加。A42科学FAIL并关闭，不因工程通过、较高目标、unknown增加或实现提速改门槛，不接入真实raw/DEV/独立物种/Setaria，不重新调角色权重、阈值、loss、seed、epoch或checkpoint救援此版本。known_nonmatch是原R4可评价非匹配，不直接等同生物假阳性；unknown继续保留为未知。

同一固定监控于03:12:59 UTC手动终态核对原33源根及两个R1根的白名单均通过，并只向原主交接一次。原1491个非空潜势、32片全部46/47行合计1491、先搜索闭合后标签与R1七证明齐，主已消费并读取完整科学汇总；该次是手动终态回传，不冒充新的定时触发。已验证数值与标签配对不解释为何正确下降；后续保存JSON的链竞争诊断只区分具体目标疑点和竞争关系，不把探索性归因当作新的确认性结论或继续开放A42。

## 保存结果竞争诊断的实际结果

13423863于2026-10-07 03:18:56 UTC提交、03:18:57开始、03:19:10结束，13秒；root/batch/extern全部COMPLETED0:0，STATUS/JOBID相符，stderr空。实际private/gpu034/0GPU2CPU8GiB5分钟/generanno/no-requeue/Restarts0；batch MaxRSS753200K，Command/SubmitLine/WorkDir及唯一日志匹配。输出根outputs/M28-A42-EMISSION-LATTICE-COMPETITION的配置、四个实际source_jobs、1491窗/16321原候选/3630新已选、两物种原计数/2916及1086分母/原五条件均与已闭合R1结果一致。三个实际保存行证明均true；只读取约定三份JSON，没有新搜索、tensor/研究权重/参考打开、标签重配、optimizer或新数据访问。提交前源文本/AST/SBATCH、五项活调度/节点/预约及精确输出冲突通过，无诊断R1或资源延长。

按完整exon几何而非索引分解如下。每格为保留／丢失／新增的原TRAIN窗口—链配对，长链为正确且至少5-exon。

| 物种 | 正确 | known_nonmatch | unknown | 长链正确 |
| --- | ---: | ---: | ---: | ---: |
| 拟南芥 | 737／794／11 | 99／248／1360 | 80／235／327 | 66／482／6 |
| 水稻 | 233／245／0 | 104／147／514 | 23／109／142 | 9／143／0 |
| pooled | 970／1039／11 | 203／395／1874 | 103／344／469 | 75／625／6 |

正确净变化11−1039＝−1028，长链6−625＝−619，known_nonmatch1874−395＝+1479，unknown469−344＝+125，分别精确解释原闭合汇总，而非新的性能评价。1039个丢失正确配对中1034个的全部已选重叠竞争者标签集合仅为0，5个仅为−1；625个丢失长链正确配对中621个仅为0，4个仅为−1。拟南芥正确丢失794中789仅0/5仅−1，水稻245全仅0；长链分别478/4及143/0。标签集合按每个旧丢失链计一次，共享竞争者不得当独立样本；0仍是原R4已知非匹配，−1仍受未知保护，不改称生物假阳性。

全部1778个旧已选正分丢失配对（1039正确、395已知非匹配、344未知）通过必要局部替换检查；所有对应新重叠链gain总和均严格高于旧链gain，最小差0.12405398488044739，非数值平分。没有发现删除这些重叠链再换回单条旧链即可提高保存目标的具体反例。这不是独立全空间最优证明，也不能排除所有实现问题；旧链、共同几何分数、DNA语法及全量角色算术的证明仍来自已闭合原CPU核验，不冒充本JSON诊断重新计算。

1491窗的成对模式：拟南芥478窗目标增加且正确减少、135窗目标增加且正确不变、132窗两者均不变；水稻对应198/182/366。合计676窗目标增加而正确减少，317窗目标增加而正确不变，498窗两者均不变；没有正确数增加的窗口。11个新增正确配对发生在同时有其他正确链损失的窗口，不能当作11个净改善窗口或独立基因。

当前可支持的研发判断：在这组固定TRAIN、A3发射及单位角色势下，扩大合法结构支持提高了保存目标，却使原正确完整链主要让位于可评价非匹配结构。故继续放宽同一搜索或改变其计算效率不构成准确率改进方案；A42保持科学FAIL关闭，不通过过滤标签0/−1、参考注入推理、改权重/阈值或挑checkpoint救援。下一独立机制应检验评分对实际结构竞争的区分能力，而不是重复同势搜索。

这仍是探索性问题定位，不证明有限候选训练是因果，也不证明所有加性势或结构训练失败。现有A36已经通过局部发射训练compatible-known-set似然，不能把“新增结构loss”包装成未尝试的新机制；其训练支持是闭合有限原候选，非整个合法语法。下一机制若使用完整语法竞争或修改训练曝光，应另定训练目标、unknown保护、监督支持、固定曝光与公平比较方案；本诊断不自动批准首训、改变任何冻结规则、开放raw/DEV/独立面板/Setaria/A33，亦不是论文claim或公开发布。
