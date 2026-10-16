# A17 局部上下文生成器：执行记录

2026-10-04，UTC 12:31／Asia-Bangkok 19:31的首次提交记录保留在下文。该fit后来TIMEOUT、下游未执行，无可恢复checkpoint；不视为TRAIN科学gate FAIL。2026-10-05的实际取证与保持科学配置的工程恢复见[恢复记录](a17_local_context_recovery.md)。原[固定协议](a17_local_context_protocol.md)不变。

## 已完成的实施与预检

新增独立 `query_local_context.py`、三项检查、`run_query_local_context.py`、`verify_query_local_context_train.py` 和四份 SBATCH，不修改原 A3/A15 源码、权重或冻结结果。局部 pooling、全窗 pointer、teacher/free 共用入口和完整 seed0 初始参数身份已落实。local/remote AST 与 bash 语法通过，显式传输的新文件与协议逐字符回读一致，没有同步仓库或提交 Git。

预检作业 13367495（gpu034，private／1RTX3090／2CPU／16GiB）root/batch/extern 均 COMPLETED／0:0，UTC 12:26:05–12:27:41，96秒，STATUS／JOBID正确。三项新检查通过（1.682秒），八个原 metadata-selected TRAIN case 的实际梯度检查通过（36.471秒）；shared／segmentation／query 的梯度有限且可达，权重不变、optimizer0、无监督checkpoint继承。完整649,457参数与A3从头seed0初始化逐tensor相同。

这只确认架构路由与训练工程可执行，不是自然正确链改善、F1或论文证据。预检没有DEV／test／Setaria或GEN编码器前向。

## 已提交的三段链

| 阶段 | job | 提交命令 | 已核实状态 |
|---|---:|---|---|
| 固定4,608步fit | 13367504 | `sbatch --parsable sbatch/M28-A17-LOCAL-CONTEXT-FIT.sbatch` | RUNNING gpu034 |
| 末权重自然TRAIN池 | 13367506 | `sbatch --parsable --dependency=afterok:13367504 sbatch/M28-A17-LOCAL-CONTEXT-TRAIN-READOUT.sbatch` | PENDING，有效afterok:13367504 |
| 只读CPU数值核验 | 13367508 | `sbatch --parsable --dependency=afterok:13367506 sbatch/M28-A17-LOCAL-CONTEXT-RESULT-CHECK.sbatch 13367504 13367506` | PENDING，有效afterok:13367506 |

资源为协议中的 generanno／private，fit/readout各1RTX3090／2CPU／16GiB／45分钟，CPU 0GPU／2CPU／8GiB／10分钟。提交前五项实际Slurm调度检查通过，private存在可用3090 slots，无切断本轮的维护；BeeGFS可用155TiB。其他TE作业和所有已有任务保持不动。

输出分别在 `outputs/M28-A17-LOCAL-CONTEXT-FIT`、`TRAIN-READOUT`、`RESULT-CHECK`。自然生成先关闭全部1,491窗自由池，再打开TRAIN catalog标注；CPU对真实缓存DNA/标签、更新曝光与三项预定筛选作独立重算。尚未读取或解释任何末模型结果。

## 固定等待与下一决定

复用原监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3` 和同一 heartbeat `coding-rna-r5-light-monitor`，gpt-5.6-luna／low／15分钟／failed_runs_only，主模型不改。监控对象为A15原三根13367306/08/10及A17新三根13367504/06/08；两条流程中任一完整终态或新故障才一次交接，然后暂停同一心跳，由主对话处理并为剩余流程重新启用。正常状态静默，不读性能，不创建第二自动化或对话。

完整新prompt及实际TOML的ID/name/kind/ACTIVE/目标对话/间隔/通知策略/原created_at逐项精确核对通过。随后手动启动确实追加到原监控对话，turn `01a106e7-5824-7df3-804d-1a83b3b5c3ee` 在 UTC12:33:07 开始、12:34:27 completed／error-null，实际turn_context仍为gpt-5.6-luna／low。12:33:18.574的scheduler snapshot显示A15/A17两fit均RUNNING gpu034，四个下游均为有效Dependency；未做未完成产物检查、性能解释或主对话健康消息。

这是新六根范围的一次手动取证，不是定时运行或终态交接；新范围定时执行尚未观察，不声称完整自动回传已经成功。原A15单流程在同一对话12:09及12:24的定时运行已核实，但不能充当新范围的验证。新预检由主对话直接消费，不重复加入监控。

接续时先核对三根terminal/exit/STATUS/JOBID和全部产物，再消费CPU实际重算的完整TRAIN读出。三项生成机制筛选全部通过，仅支持另写匹配caller／整链评分协议；失败则保留并关闭本架构，不调宽度/seed/末权重/阈值救本轮。两者都不构成完整raw-FASTA科学PASS，也不自动进入DEV、封存数据或发布。

自动等待需要电脑开机、应用运行及SSH可用，不能保证休眠／断线期间即时回传，见[官方定时任务说明](https://learn.chatgpt.com/docs/automations?surface=app)。长期论文目标仍未完成；本轮没有修改goal状态或恢复历史自动研究框架。
