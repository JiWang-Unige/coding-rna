# A15 自然候选域运行记录

2026-10-04 11:56 UTC的首次提交记录保留在下文。原fit后来TIMEOUT，无可恢复checkpoint，原读出与核验未执行；这不是TRAIN科学gate的FAIL。故障和保持科学配置的工程恢复见 [恢复记录](a15_natural_chain_set_recovery.md)。冻结选择见 [A15协议](a15_natural_chain_set_protocol.md)，A14原FAIL不变。

## 已核实的预检

预检13367298的root／batch／extern均COMPLETED、0:0，实际140秒；Slurm显示13:49:58–13:52:18 CEST，即11:49:58–11:52:18 UTC。STATUS与JOBID匹配。四个新候选域检查PASS（0.451秒）；八个按原metadata规则在自然视图选出的真实窗口，exon／sequence／additive／nonadditive／null均有限非零梯度，head权重不变，生成器／link无梯度，optimizer更新0。实际梯度循环27.419秒，不含准备读取开销。

1491唯一TRAIN窗的自然候选／窗口对16321，标签1／0／-1为2380／12566／1375，移除536注入正链；72个窗口失去全部gold。唯一窗零已知为2。原1536 draw每遍重放自然16806，标签2473／12935／1398，零已知draw为3，注入0。与A14的draw／optimizer曝光相同，但正例曝光减少。

## 已提交的固定链

以下均在远程 `/home/users/j/jwang/coding-rna` 提交。11:53 UTC即时squeue显示fit运行gpu034，两下游PENDING／Dependency；scontrol核对实际依赖、命令和资源相符，不以提交当完成。

|阶段|实际作业|命令|资源／时限|
|---|---|---|---|
|fit|13367306|`sbatch --parsable sbatch/M28-A15-NATURAL-CHAIN-SET-FIT.sbatch`|private／1RTX3090／2CPU／16GiB／45分钟|
|末004608 TRAIN读出|13367308|`sbatch --parsable --dependency=afterok:13367306 sbatch/M28-A15-NATURAL-CHAIN-SET-TRAIN-READOUT.sbatch`|private／1RTX3090／2CPU／16GiB／20分钟|
|只读结果核验|13367310|`sbatch --parsable --dependency=afterok:13367308 sbatch/M28-A15-NATURAL-CHAIN-SET-RESULT-CHECK.sbatch 13367306 13367308`|private／0GPU／2CPU／8GiB／5分钟|

日志按M28A15FIT／M28A15TRAINREAD／M28A15RESULTCHECK及实际job ID保存在原 `outputs/M28-PILOT-R5/logs`。独立输出路径按协议，不覆盖任何旧结果。四段均使用generanno；本轮没有DEV／sealed／Setaria访问、其他作业修改或取消、Git提交或公开发布。

## 文件及等待交接

仅新增自然域投影、四检查、独立runner／只读核验脚本、对应sbatch与协议；未修改A14代码和checkpoint。Python AST、bash语法通过，新传输文件逐个与本地直接字节比较一致，不使用哈希。上次A14 tar生成的七个AppleDouble元数据旁文件已确认类型并可恢复移到 `/home/users/j/jwang/m28-a14-appledouble.FbeIpm`，没有删改科学文件；本次用明确逐文件scp，没有同步整个项目。

固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3、heartbeat coding-rna-r5-light-monitor复用；gpt-5.6-luna／low／15分钟，完整新prompt、kind、ACTIVE、目标对话、间隔及原通知设置已实际保存回读一致。人工启动的11:55:19 UTC快照确认fit RUNNING／两下游有效Dependency，原监控turn01a106c4-7294-7a73-965f-c85568357391完成且error-null，实际turn_context仍gpt-5.6-luna／low；未向主对话发送健康回执。此人工检查不算定时证据，本A15批次定时检查及最终回传尚未观察到。正常运行／依赖静默，全部三根作业结构完成或新可行动故障才一次交接后暂停。本机需保持开机和应用运行，见[官方定时任务说明](https://learn.chatgpt.com/docs/automations?surface=app)。

未完成：4,608更新末权重训练、全16321自然候选读出、CPU实际记录核验及原A6五条件判断。主对话在整链终态后解释A6／A14／A15配对结果；PASS也不等于原生raw-FASTA F1、背景、独立泛化或论文成立，FAIL不加轮次、扫参数或选择其他checkpoint。

## 首次运行的后续终态

2026-10-04 12:38:30 UTC，13367306 root TIMEOUT／0:0、batch CANCELLED／0:15、extern COMPLETED／0:0，实际45分28秒；stderr明确time limit。原日志持久化至step896，没有001536或其他checkpoint，也没有summary；不能把日志当作模型／optimizer状态。可能存在未刷盘的尾部更新，不把896宣称为全部物理更新。13367308进入DependencyNeverSatisfied，13367310等待Dependency；旧链不覆盖、不取消，不再重复监控或交接。

固定监控于12:40:11.910 UTC的实际定时检查发现故障，12:40:30 UTC向本主对话交接，原生PAUSED已回读；对应turn01a106ed-6230-7461-b37d-35a3dc361d35完成且error-null。主对话独立核实并消费后，仅对仍有效的A17三段链重新启用同一固定心跳。A15零更新计时13367598由主对话直接处理，恢复配置须依据实际计时决定，不能把TIMEOUT改写成模型性能结论。
