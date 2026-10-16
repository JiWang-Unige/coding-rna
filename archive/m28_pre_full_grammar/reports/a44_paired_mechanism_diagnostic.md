# A44 配对正确链分数与选择诊断

2026-10-07，真实结果读取前固定。本诊断解释同一A44 LOCAL／JOINT末模型的原TRAIN结果，区分正确链分数、有限支持中的已知冲突分数差与最终自然选择；它不增加训练、搜索或新的科学门槛，也不救活失败的结构机制。

## 来源与运行条件

仅使用原13424089闭合支持中的19189窗口—链配对、原1491唯一TRAIN窗口及2916完整gold配对（至少5-exon为1086），和当前两臂末004608的保存势、自然预测及闭合读出。两fit为13425946／13425945，势13426150、64解码13426151_0至63、读出13426152、CPU核验13426153。实际源需整链终态、明确STATUS／JOBID及原七项数值证明；不能凭训练作业结束或不完整输出开始诊断。实际完整分母与正确／长链选择必须回到原逐物种读出。

科学FAIL且工程／数值闭合的结果允许诊断，保留原两组五条件和原TRAIN_preparation_supported布尔值；不得因分数改善、平均margin改善或LOCAL恢复而改成PASS。数据或数值核验失败则不执行。本诊断不加载模型、GEN、feature tensor、原始注释catalog、DEV或独立物种，不执行前向、反向或优化器，只读取已有TRAIN支持元数据及两臂保存的势／预测。

## 描述量与解释范围

每条标签1完整链用原FP64单位角色相对intergenic的logit差求和，另分为coding phase、intron phase、donor phase、acceptor phase、start和stop六项。六项和必须等于原完整gain；它不是gene概率或逐碱基准确率。保存JOINT减LOCAL的同链gain及六项变化，零阈值正／非正与最终是否选中分别记录，不用分数符号推定选择结果。

有限冲突margin仅定义为该正确链gain减去固定支持中与其span重叠的最高标签0链gain，两个模型各自取最大值。标签−1不当作负例，不重分类；没有标签0冲突时为null，并单列其分母。该差值不是全语法regret、完整兼容集合之间的最优差，也不代表所有推理时竞争结构；允许各臂的最高冲突身份不同。有限支持含原gold注入和旧机制产生的几何，其覆盖与选择偏倚不能由本诊断消除，分数变化也不是额外监督与梯度匹配的因果效应。

按同一窗口／链几何记录both_selected、newly_selected、lost_selected及neither_selected；逐物种、至少5-exon完整分母与原读出核对，再分别报告各组gain、有限margin及角色变化的样本数／均值／中位数／范围／符号数。窗口—链配对不冒称去重生物gene，标签0仅是原参考口径下的已知非匹配，不是生物真负。正gain但未选中只能表明未达到最终自然选择，不能单凭这些统计断定某个未保存竞争路径是原因。

## 工程检查与执行

新增四个人工用例覆盖实际配对来源与七证明（允许真实科学FAIL、拒绝旧源／假PASS）、独立写明的六角色数值与共同logit平移不变性、极大未知分数及不重叠非匹配不制造负冲突、完整链几何选择转换及无冲突null分母。只使用人工tensor和既有人工收据，不读取实际保存势、模型或结果。检查通过不等于实际诊断或性能收益。

检查脚本M28-A44-PAIRED-MECHANISM-CHECKS.sbatch使用generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立outputs/M28-A44-PAIRED-MECHANISM-CHECKS及M28A44MECHCHECK_%j日志。实际诊断M28-A44-PAIRED-MECHANISM.sbatch使用同环境／0GPU2CPU8GiB15分钟／no-requeue，独立outputs/M28-A44-PAIRED-MECHANISM及M28A44MECH_%j日志；预计输出小于16MiB，只流式加载当前窗口两臂势，不覆盖源。人工检查闭合后可用afterok:13426153依赖提交，入口仍验证全部69源及实际人工检查job，未满足依赖不读取数据。

两阶段均只一次；明确实现错误可保持原数据、源、资源、统计和科学门槛，用独立R1根作一次最小修复，数据／数值矛盾、TIMEOUT／OOM或第二次失败停止，不延长时限。实际诊断提交后才将它加入原固定监控的同一批，复用同一对话及heartbeat，不新增监控；源69根加诊断为70根。监控仍只检查状态和完成产物，不解释margin或改门槛。整批终态或可行动故障一次交接原主，主负责结果解释。原A36／A42／A43 FAIL及A44原停止条件不变，真实raw／独立数据和公开发布仍未由本诊断放行。

## 实际检查与依赖提交

13426518于2026-10-07 09:39:29–09:39:41 UTC闭合，root／batch／extern均COMPLETED0:0，用时12秒、batch MaxRSS758508K；四项人工检查PASS0.180秒。实际generanno／private／gpu035／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，STATUS／JOBID／Command／SubmitLine／WorkDir及唯一M28A44MECHCHECK_13426518日志相符。新增两Python AST、两SBATCH bash语法及明确五文件canonical逐字一致；没有实际势、模型、预测或标签参与检查。

09:41:17 UTC五项活调度／预约、精确空根与日志及存储核对通过，账户展开71项加本诊断为72，小于80；BeeGFS可用156T。随后执行 `sbatch --parsable --dependency=afterok:13426153 --export=ALL,CHECKS_JOB=13426518 sbatch/M28-A44-PAIRED-MECHANISM.sbatch`，实际09:41:18 UTC提交13426529。09:43:04 UTC主核对PENDING普通Dependency／0:0、afterok13426153未满足及准确Command／SubmitLine／WorkDir／exports／0GPU2CPU8GiB15分钟／no-requeue／Restarts0和唯一日志；只证明依赖接线成功，不是实际诊断结果。

原heartbeat coding-rna-r5-light-monitor仅扩展本批为70根，原固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3不变。原生更新读回ACTIVE／kindheartbeat／15分钟／failed_runs_only及完整11035字符prompt逐字相符，created_at1790981930511保持、updated_at1791366246794；新增诊断使用独立protocol及source_configuration结构白名单，不误套原源configuration字段。配置、扩展后的实际手动和定时检查均已验证；完整自动终态回传尚未发生。仍无性能、机制收益或论文完成结论。

扩展后manual turn01a115c1-019c-7690-86b6-94a40f572280实际completed／errornull／30531ms／cursor145；09:45:32.896 UTC的实际turn_context为原gpt-5.6-luna／low／cwd。09:45:52 UTC最小状态检查见两fit及势COMPLETED0:0、64解码均RUNNING、读出／CPU／诊断均正常PENDING Dependency0:0；未读取产物或科学值，未交接。本批配置与一次真实手动检查已验证；扩展后的实际定时证据如下，不用原69根的定时证据替代。

70根扩展后的实际定时 turn01a115c8-23ea-7102-96c8-feeb29e12cab已由主消费，completed／errornull／25255ms／cursor146。原固定监控的实际09:53:20.410 UTC turn_context仍为gpt-5.6-luna／low／原cwd；自动触发current_time_iso09:53:20.360 UTC含准确70根cohort。09:53:20 UTC记录两fit及势COMPLETED0:0、64解码RUNNING、读出／CPU／机制诊断正常PENDING Dependency0:0，无产物读取或交接。配置、手动检查与扩展后的定时检查均已验证，完整自动终态回传仍未发生。
