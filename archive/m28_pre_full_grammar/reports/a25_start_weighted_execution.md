# A25 首起点提示生成器：执行记录

2026-10-06。本轮已按[固定协议](a25_start_weighted_protocol.md)完成，原三gate全FAIL，关闭A25。不更改A17/A23或读取DEV、sealed、Setaria。结果见[a25_start_weighted_result.md](a25_start_weighted_result.md)。下文提交和启动快照为历史取证。

## 终态

13392060／13392061／13392062三根及各batch、extern均COMPLETED0:0，实际59:46／44:16／14:19，STATUS和JOBID匹配。完整原4608曝光、1491窗自然池及CPU核验通过；合并正确候选2380→2315，≥5外显子788→721，known_nonmatch12566→17036。原科学gate全false，不重训补救。监控错误使用inference字段名；实际协议字段new_head_or_backbone_forward=false，主对话已核对并排除误报。固定监控已回传并暂停，旧作业全部已消费。

## 已完成预检

预检13391985 root/batch/extern均 COMPLETED、ExitCode0:0、Elapsed38秒；STATUS=COMPLETED、JOBID一致。五项合成测试全部通过（0.693秒），八个原metadata固定真实TRAIN case零更新检查13.771秒。649585参数、共有fresh A17参数、零prior时实际四旧loss与离散匹配一致，所有参数不变，optimizer0/nocheckpoint继承/no新GLM。两物种eligible正链case新initiation层梯度有限非零；无positive case新增CE=0，但该层仍可通过原context/existence路径接到梯度，这是预期，不是给背景增加辅助absence监督。

## 已提交正式链

实时Slurm五项资源/队列/维护检查支持private-teodoro-gpu：gpu034配置8卡/用7、gpu035配置8/用6，RTX3090每卡24GiB；没有相关GPU维护冲突、并发限额冲突或新阶段目录碰撞。BeeGFS可用136TB。沿用已在实际M28作业验证的generanno，不往base安装依赖。新代码与SBATCH经逐文件直接字节比较确认本机/远程一致；两正式Python runner语法及三个SBATCH语法通过。

提交 `sbatch --parsable sbatch/M28-A25-START-WEIGHTED-FIT.sbatch` 得13392060，提交时 RUNNING gpu034；private/1RTX3090/2CPU/16GiB/3小时。训练本身读取匹配COMPLETED预检产物，不附已完成短预检的Slurm依赖。

提交 `sbatch --parsable --dependency=afterok:13392060 sbatch/M28-A25-START-WEIGHTED-TRAIN-READOUT.sbatch` 得13392061，PENDING有效Dependency；同GPU/CPU/memory，上限2小时。

提交 `sbatch --parsable --dependency=afterok:13392061 sbatch/M28-A25-START-WEIGHTED-RESULT-CHECK.sbatch 13392060 13392061` 得13392062，PENDING有效Dependency；private/0GPU/2CPU/8GiB/30分钟。

输出分别为 `outputs/M28-A25-START-WEIGHTED-FIT`、`...-TRAIN-READOUT`、`...-RESULT-CHECK`，日志为 `outputs/M28-PILOT-R5/logs/M28A25FIT_13392060.*`、`M28A25READOUT_13392061.*`、`M28A25CHECK_13392062.*`。fit按原4608终点、每256工程checkpoint执行；CPU复算新首起点监督支持/原draw-weight-query政策、末新增层权重变化、自然池标签及原三项gate，不新增模型前向或搜索。

## 历史等待与启动记录

正式拟合、完整1491窗自然读出及CPU三gate尚未完成，未作任何覆盖提升或科学PASS判断。A23仍独立运行。等待复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor，更新对象为A23+A25，沿用gpt-5.6-luna/low、15分钟。A17原gate失败已消费，不继续检查或重报。新cohort配置、启动检查、定时检查、完整终态回传分别取证，不将配置成功当作自动回传成功。

本轮启动取证：fit13392060 RUNNING/STATUS RUNNING/JOBID一致，日志已到128/4608，stderr无错误；首步新增initiation层梯度norm10.392739167065962，有限非零。该工程连接证据不是候选覆盖收益，后续仍按末004608及三gate判断。

同一原生心跳于2026-10-05T17:05:14.135Z保存 ACTIVE；整个prompt和id/name/kind heartbeat/target/15分钟/failed_runs_only/原created_at1790981930511直接读回完全一致。新cohort手工启动turn01a10d07-439f-72e3-9d02-79b6778eaf5c在原固定对话完成、无错误、23.849秒，实际调度快照17:05:55.831Z：两个fit运行，各下游有效Dependency。实际turn_context17:05:42.342Z仍gpt-5.6-luna/low，无模型覆盖；没有新对话、fork或cron。新cohort定时运行及完整终态回传尚未观察到，不能由手工检查推定；上一cohort16:54的定时检查记录保留为历史。
