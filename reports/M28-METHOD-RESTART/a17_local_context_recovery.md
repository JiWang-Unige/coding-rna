# A17 超时取证与唯一工程恢复R1

2026-10-05，恢复训练提交前固定。A17原45分钟时限未覆盖实际4608步训练，未完成模型，没有自然候选读出或科学gate结果。此次仅修正运行时限、保存频率和输出路径，保持模型、seed、输入、全部loss／匹配／标签规则和停止终点不变。

## 实际失败与证据

13367504于2026-10-04 12:28:41–13:14:01 UTC运行，root TIMEOUT／0:0、batch CANCELLED／0:15、extern COMPLETED／0:0，45分21秒；`outputs/M28-PILOT-R5/logs/M28A17LOCALFIT_13367504.err`明确time limit。原fit目录只有JOBID／遗留RUNNING标记及训练日志，没有任何checkpoint或summary。日志最后持久化step576／epoch1／ordinal575，不保证是全部物理更新，未刷盘尾部未知。无可恢复模型／AdamW状态，不能声称精确续训。

原readout13367506和CPU核验13367508实际CANCELLED by 0／0:0，Start None、Elapsed0；未执行，不称为科学失败。这里不推定取消者或原因，也没有由本Agent取消、重置或覆盖旧作业。原目录、日志及所有科学输入保留。

固定监控于2026-10-05 13:45:06.845 UTC实际取得终态，发送一次新故障交接并原生暂停，主对话独立核实。该heartbeat的触发时间仍记录2026-10-04 13:15，实际取证晚于触发，不能声称休眠／断线期间即时回传；同一monitor turn完成且error-null，不创建替代监控。

## R1范围与操作时限

从完整A3相同seed0 fresh649457参数初始化开始；原已通过三项检查和八个真实梯度case的预检仍适用。独立`recover_query_local_context.py`保留原fit计算路径：每draw读取原cache、`window_inputs`、局部context模型、原`supervised_readout`和同FP32 AdamW／clip1。未修改原runner、架构、任何loss、Hungarian匹配、unknown／partial规则、原1491窗／1536 draw的顺序和importance权重。3遍＝4608更新，不从旧A3或失败A17权重继承。

每256步保存同格式model／AdamW状态，包含001536／003072／004608；只有末004608可进入科学读出，不选较早checkpoint。失败尝试的576条持久更新和未知尾部是额外丢弃计算，不宣称所有物理更新合计4608。

独立输出`outputs/M28-A17-LOCAL-CONTEXT-{FIT,TRAIN-READOUT,RESULT-CHECK}-RECOVERY-R1`，三段afterok链，直接复用原readout和CPU核验，仅改变输出路由。generanno／private、两GPU阶段各1RTX3090／2CPU／16GiB，CPU0GPU／2CPU／8GiB不变。操作时限固定为fit8小时、自然TRAIN读出4小时、CPU核验2.5小时；fit根据45分钟仅持久化576步的观察外推约6小时再留余量，不是精确ETA。读出与CPU标签核验都要重新读取全部1491原cache，A15同源读取计时提示不能继续套用45／10分钟；A15的head阶段计时不能代表A17生成／梯度开销。

仍先关闭全部1491自然greedy＋beam4池，再打开TRAIN catalog；无参考注入／DEV／sealed／Setaria／新GLM编码器前向。CPU仍核对原4608 draw／权重／A3监督身份、真实新旧自然链标签、全有限末模型以及原三项generator gate。PASS仅支持准备完整caller/head协议，FAIL不改seed、宽度、阈值或曝光救参。未知／已知非匹配候选数单列，不包装成选中基因F1或泛化。

## 实际提交与等待

2026-10-05 13:54:41 UTC，fit13387910实际启动gpu034，与A15 R1 fit13387908并行。13:56–13:57 UTC即时squeue／scontrol确认两fit RUNNING、四下游均有效Dependency，实际private、GPU类型／CPU／memory／命令／参数／时限与事前配置一致。

|阶段|实际job|实际依赖／状态|操作时限|
|---|---:|---|---|
|A17 R1 fit|13387910|RUNNING gpu034|8小时|
|末generator自然TRAIN池|13387946|afterok:13387910，PENDING Dependency|4小时|
|原CPU标签／监督／gate核验|13387978|afterok:13387946，传入13387910／13387946，PENDING Dependency|2.5小时|

SBATCH为`M28-A17-LOCAL-CONTEXT-RECOVERY-R1-{FIT,TRAIN-READOUT,RESULT-CHECK}.sbatch`，日志前缀M28A17R1FIT／M28A17R1TRAIN／M28A17R1CHECK加真实ID。R1 Python AST／bash语法通过；完整逐draw循环以AST去除唯一新增的256保存分支及显式R模块引用后，与原fit严格相同。原科学runner未改动且remote/local直接字节比较一致，全部R1传输文件也直接比较一致，不引入哈希或重复架构审计。

同一固定监控对话／heartbeat已原生启用，完整新cohort和prompt／固定目标／15分钟／原通知设置／created_at实际回读相同，gpt-5.6-luna／low保持。当前六root scope、交接和限制详见[A15共同等待记录](a15_natural_chain_set_recovery.md)。原A17故障及未执行旧下游不再监控或重复交接；本R1范围尚无定时执行或终态回传证据。

未完成：4608更新末模型、1491自然自由池、真实DNA／候选标签与监督曝光核验、原三项TRAIN gate。当前训练运行不等于模型成功、泛化、独立F1或论文结论。本次不取消或更改其他作业，不提交Git／公开发布。
