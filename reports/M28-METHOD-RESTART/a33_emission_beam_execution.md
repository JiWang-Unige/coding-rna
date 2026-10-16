# A33 发射分数驱动搜索：实际执行

2026-10-06，主对话持续科研目标保持ACTIVE。本轮依据预先冻结的[a33_emission_beam_protocol.md](a33_emission_beam_protocol.md)，没有新参数／优化更新、阈值扫描、DEV／Setaria／sealed访问、Git提交或对外发布。

当前结论：八窗预检通过，但完整TRAIN分片0在指针重放检查FAILED，合并／独立CPU没有执行，五项科学条件均未评价。只读诊断13403620证明失败窗按原beam批量计算路径可逐链精确复现，单路径重放则有小幅浮点差异；原搜索、容差和停止条件未改，不能据此把原FAILED改为完成或报告科学PASS。

## 已闭合的工程证据

原预检13402858 FAILED／1:0，但五项搜索合成检查PASS／0.888秒。首个真实窗在新搜索前被旧坐标tuple与JSON list的直接相等比较拦住；自然输出0行。原记录保留，唯一R1修复只规范容器表示，不放松身份／顺序／来源或数值条件。

R1预检13402927 root／batch／extern全部COMPLETED／0:0，实际10:00:06至10:01:29 UTC，1:23，batch MaxRSS1773576K；八行自然输出和匹配STATUS／JOBID／summary齐全。八窗旧自然池与A31_A3分数复现、新完成链raw grammar和指针／EOS重放通过、增量gain与直接log-probability差求和误差均0。指针重放最大绝对误差2.0796884200535715e-05，满足原FP32容差，不宣称零误差。128个query深度2至24，cap存活0，模型权重不变、无梯度。内部八例43.422秒；新增唯一链5／7／6／12／8／16／7／23，未标注正确性；没有科学gate结果。

完整读出新实现检查13403081 root／batch／extern全部COMPLETED／0:0，27秒，batch MaxRSS748416K，实际三项测试PASS／0.005秒：1491窗四分片完整唯一合并及漏／重窗拒绝、五项整数gate不混淆F1与precision／长链收益、原GT／GC／N-intron及codon carry语法和非法链拒绝。远程新Python AST及四SBATCH语法检查通过。

## 已提交的完整 TRAIN 链

| 阶段 | Job | 资源／依赖 | 10:16 UTC实际状态 |
| --- | --- | --- | --- |
| 自然读出四分片 | 13403217_0至_3 | 各private RTX3090×1／2CPU／16GiB／2h；afterok:13403081 | 四项均RUNNING gpu034 |
| 合并与R4标注 | 13403223 | private／0GPU／2CPU／8GiB／40min；afterok:13403217全部任务 | PENDING，有效Dependency |
| 独立CPU核验 | 13403224 | 同CPU资源；afterok:13403223 | PENDING，有效Dependency |

自然array实际ArrayTaskThrottle=4；scontrol核实每项申请／分配1 RTX3090与原2CPU／16GiB／2h。下游实际无GPU请求，依赖分别为afterok:13403217_*和afterok:13403223。array任务0的独立Slurm数字JobId13403218与alias13403217_0为同一任务，输出JOBID采用alias，不混淆任务身份。

提交命令按实际先后：

```bash
sbatch --parsable sbatch/M28-A33-EMISSION-BEAM-TRAIN-CHECKS.sbatch
sbatch --parsable --dependency=afterok:13403081 sbatch/M28-A33-EMISSION-BEAM-TRAIN-NATURAL.sbatch
sbatch --parsable --dependency=afterok:13403217 --export=ALL,A33_ARRAY_JOB=13403217 sbatch/M28-A33-EMISSION-BEAM-TRAIN-READOUT.sbatch
sbatch --parsable --dependency=afterok:13403223 --export=ALL,A33_ARRAY_JOB=13403217,A33_READOUT_JOB=13403223 sbatch/M28-A33-EMISSION-BEAM-RESULT-CHECK.sbatch
```

原1536 draws首次出现的1491唯一窗按ordinal modulo 4分为373／373／373／372，合并恢复原次序。原A31_A3池完整保留，冻结模型新forward复现旧分数，再加入新搜索链。全部自然文件关闭、四任务成功后才打开原TRAIN R4 catalog；原2916 eligible配对分母和unknown规则不变。独立CPU从原序列／R4及保存FP32 logits重算语法、直接分数、标签、非重叠选择、分母和五项gate。新搜索计算量增加，不能声称matched compute。

产物根outputs/M28-A33-EMISSION-BEAM-TRAIN-NATURAL/shard_0..3、outputs/M28-A33-EMISSION-BEAM-TRAIN-READOUT、outputs/M28-A33-EMISSION-BEAM-RESULT-CHECK；预期自然文件共1491窗，合并natural_scores.jsonl及windows.jsonl各1491行。日志outputs/M28-PILOT-R5/logs/M28A33NATURAL_13403217_i.{out,err}、M28A33READOUT_13403223.{out,err}、M28A33RESULT_13403224.{out,err}。尚未闭合这些产物或解释结果，不由提交／运行推定支持条件通过。

提交前实际五项scheduler检查无并发／array／显存／维护冲突，private034有4／8 GPU空闲、035闲置，选择private；BeeGFS实际155TiB可用，预计新增保存logits约2.2GiB。保留旧A30等待依赖与其他TE作业，未取消、移动或修改它们。

## 固定监控和后续判断

仅复用原监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3与原automation coding-rna-r5-light-monitor；原生更新kind heartbeat／固定target／15分钟／failed_runs_only／ACTIVE，实际读回并逐字匹配新prompt。沿用gpt-5.6-luna／low，无新chat、fork、cron或主模型变更。新cohort M28-A33-EMISSION-BEAM-TRAIN-20261006只监控四自然任务及两CPU根；旧已消费作业不再检查。实际手动启动turn01a110b9-3f9e-7590-b4d4-9b540c36de0c完成／error-null／31.115秒，10:19:11.501Z取证四自然任务RUNNING、两CPU有效Dependency，未读中间产物或回传主对话；实际turn_context10:18:58.402Z仍为gpt-5.6-luna／low及同项目。本cohort定时执行与完整终态回传尚未验证，不能继承旧cohort成功。

常规RUNNING／PENDING／有效Dependency静默；六根整链完成及结构齐全，或新可行动故障，监控只发一次证据接续并暂停同一心跳。主对话负责核实五项原条件：两物种proxy F1各严格升、两物种proxy precision各不降、合并≥5 CDS selected correct严格升。任一FAIL关闭固定配方；全部通过也仅支持另定raw-FASTA协议，不是native gene F1、独立泛化、机制因果证明或论文PASS。当前剩余工作是整链终态取证、独立CPU结果核实和科学解释；研究目标未完成。

## 完整读出的数值检查失败

同一监控实际scheduled turn01a110c5-5e5e-7c20-b98c-afbaf478e77d于10:32:27.558Z取证并交接：13403217_0 root／batch FAILED／1:0、extern COMPLETED／0:0，用时9:15，batch MaxRSS1333372K；其他三个分片仍RUNNING，两个CPU仍等待依赖。主对话已独立核实stderr、原STATUS及实际文件：137个自然记录／137份logits，最后ordinal544，失败窗为下一个ordinal548，无完成summary。不是完整1491窗读出，不是五项科学条件FAIL，也没有新效果标注。

实际错误位于原指针／EOS重放assert_close，原记录期望-13.335641860961914、独立重放-13.335674285888672；绝对差3.24249267578125e-05、相对差2.4314485268783047e-06，超过原atol1e-5／rtol1e-6。没有修改容差、权重、搜索或预算，没有重提／取消数组；保留原失败及其他分片原执行。当前不能宣称完整TRAIN工程通过或科学gate通过。

额外只读诊断在运行前固定：仅ordinal548、原冻结A3／原TRAIN cache及原新搜索、不读catalog或标签；分别按原beam全部合法end事件的批量recurrent与单路径scalar recurrent重放，并区分FP32累加、Python FP64累加及end／EOS的原FP32组合次序。核对原局部分数、原始语法、最终发射分数及冻结tensor，只记录误差来源，不改原validator或输出。独立输出outputs/M28-A33-POINTER-NUMERICAL-DIAGNOSTIC、private1RTX3090／2CPU／16GiB／5分钟、generanno，日志M28A33NUMDIAG_JOBID.*。诊断可分辨“批量／逐条算术差异”与“同计算路径仍不复现”，但不能凭诊断结果放宽既定条件、恢复原失败root或解释未执行的科学效果。

## 数值诊断实际结果与完整分片终态

诊断13403620 root／batch／extern COMPLETED／0:0，用时36秒，batch MaxRSS1665464K；STATUS／JOBID匹配，64行paths.jsonl及summary完整。实际原次序ordinal548为arabidopsis_thaliana|NC_003076.8|3612672|+，仅TRAIN cache／原冻结A3、无GEN forward、catalog、效果标签或优化更新，参数前后逐tensor相等。旧A31局部分数、全部新完成链raw grammar及直接发射gain均复现。

| 数值对象 | 64条链实际最大误差 |
| --- | ---: |
| 原批量GRU／原FP32 end+EOS组合／FP64累计，对照原beam存储分数 | 0 |
| 单路径GRU分数与批量GRU分数，两者均用原分组及FP64累计 | 3.2663403544574976e-05 |
| 单路径FP32累计与同路径独立FP64组件求和 | 9.742550901137292e-06 |
| 每exon结束state的批量／单路径最大元素差 | 2.86102294921875e-06 |

原容差在64条中仍有1条FAIL：query0／rank2，五exon链[21298,21896)、[21997,23187)、[23287,23467)、[23602,23742)、[23781,23803)。原beam与匹配批量重放均-13.335641670972109；单路径FP32为-13.335674285888672，单路径同分组FP64为-13.335674334375653。原assert将存储分数转成FP32后期望-13.335641860961914，故错误日志的绝对差与未转换的FP64差略有不同。单纯换FP64累加不能消除批量／单路径差异，不能以“只有累加精度”概括。

另由保存的逐步start／end／EOS操作数独立重算FP32 end+EOS及两种FP64组件和，64行全部复现；原tol判断与唯一失败行一致，同批量路径与原分数全部精确相等。这个取证支持验证器比较了不同浮点执行路径的解释，没有发现本窗同批量路径、语法或最终发射分数不复现；它不是全1491窗独立通过的证明。

三个其他分片现均root／batch／extern COMPLETED／0:0：_1 24:43／MaxRSS1362416K／373行及373份logits，_2 24:21／1510592K／373，_3 24:17／1524160K／372，完成summary／STATUS／JOBID匹配。加上失败分片137行，实际共1255自然窗，缺236；没有将部分覆盖当完整人口。合并13403223处于DependencyNeverSatisfied，CPU13403224等待该依赖，无新R4效果标注，不取消或覆盖它们。

恢复需要原“失败不自动重提”停止条件的一次明确工程例外，已提交具体选择：只修正与原beam一致的重放执行路径、独立目录重跑分片0，复用三个已完成分片；搜索、权重、每次执行原资源／时限、容差及五项科学条件不变，原FAILED保留。重跑增加实际工程计算量，不声称总计算预算或matched compute不变。尚未收到该例外批准，没有部署替换checker或恢复作业；继续其他已授权路线。原固定监控已实际定时发现故障、成功交接一次后PAUSED，不再轮询已消费的A33终态。

该scheduled turn01a110c5-5e5e-7c20-b98c-afbaf478e77d实际completed／error-null／62.726秒，固定监控对话仍idle，同一心跳PAUSED；配置成功、同对话实际定时检查与故障单次回传均已分别验证，完整TRAIN成功链则没有。独立人工短contig三原生预检13403567_0..2全部COMPLETED／0:0，结果见[native_short_contig_preflight.md](native_short_contig_preflight.md)；三种均真实空GFF，不作准确率／泛化或A33效果证据。

## 一次工程恢复的准备范围（未执行）

若该停止条件例外获准，最小修改只限本A33工程重放和混合源合并：

- 与原搜索相同的全部合法end事件／相同次序批量计算GRU及EOS，再取实际链的分支；严格保持原START先加至FP64累计、end+EOS先FP32相加再加至累计的结合次序。诊断trace的同批量结果在本64链为零，不以该单窗结果推定整个TRAIN零误差。
- 不全局改query_search_score_replay.replay_node_score：该函数参与既有梯度／教师路径，原scalar语义及历史复现保留。新验证仍检查raw grammar，且不放宽原atol／rtol，不改搜索、tie-break、角色分数、local top4、Q16／beam4／completed4／cap128或原五项条件。
- 先在原八个元数据预检窗及ordinal548做同批量重放；通过后只在独立R1目录重跑分片0的全部373窗，不拼接原137行。原失败分片没有最终权重不变summary，故不把其部分前缀冒充已闭合自然源。
- 合并显式核实新分片0实际成功root和原13403217_1／_2／_3三个成功root、各自STATUS／JOBID／summary、373／373／373／372实际完整文件及原1491顺序。新源不能沿用“整个原array成功”的判据；旧array仍FAILED。保持原文件，新增R1合并与独立CPU目录。
- 全部自然文件闭合后才读取原TRAIN R4标签；复用原独立CPU直接log-probability求和、raw grammar、标签／非重叠选择／2916分母及五项整数条件核验，不因数值诊断或恢复通过解释科学效果。
- 新GPU仍1RTX3090／2CPU／16GiB／2h，新CPU仍0GPU／2CPU／8GiB／40min；使用原环境。实际重跑会增加工程计算量，不作matched-compute主张。提交后只更新原固定监控／同一心跳对象；原无效依赖作业不取消、不复用。

以上是依据已保存代码和实际诊断准备的可执行范围，不是部署、作业提交、协议例外批准或科学支持结论。
