# A21 剪接连接整链头执行记录

2026-10-05。A21 固定拟合、末权重TRAIN读出和CPU核验全部完成，原五条件3 PASS／2 FAIL，候选关闭，不运行DEV或救参。完整选择数与科学边界见[结果](a21_junction_chain_result.md)，科学定义见[冻结协议](a21_junction_chain_protocol.md)。A15保持原FAIL，A17独立流程不因本候选改变。下文首次运行快照保留为历史执行证据。

## 准备结果与一次工程修复

原 preflight13389874 root／batch FAILED1:0／42秒，extern COMPLETED0:0；六项结构测试PASS／0.577秒。故障是准备检查要求所有八窗口都给连接分支非零梯度，但原第五窗口 `oryza_sativa|NC_089043.1|26025984|+` 的9个候选全部单外显子。其他七窗口都有已知多外显子候选。该尝试只做零更新检查，没有 optimizer、拟合或 checkpoint，不是科学 FAIL。

修复不替换窗口、不改模型或 loss：有已知多外显子时连接两分支必须有限非零，没有时必须无梯度或精确零；仍要求两物种都实际覆盖连接分支。原失败输出和日志保留。唯一修复重试13390039 root／batch／extern COMPLETED0:0／34秒，STATUS／JOBID匹配；六项结构测试PASS／0.344秒，实际准备16.3675秒。

实际1491唯一TRAIN窗／16321自然候选的标签1／0／-1为2380／12566／1375，零已知窗2；原1536 draw重放为16806候选／窗口对，标签2473／12935／1398，零已知draw3。共有9328多外显子候选、33946相邻连接token，没有相触连接。八个原metadata窗口全部通过适用的有限梯度检查，七个适用连接分支；两物种有覆盖，所有未知独立gain的直接梯度为零，参数不变、optimizer更新0、生成器冻结。它不证明独立基因F1、改进原因或泛化。

## 已提交流程与验收产物

拟合与读出 root13390188，命令为 `sbatch --parsable --dependency=afterok:13390039 sbatch/M28-A21-JUNCTION-CHAIN-FIT-READOUT.sbatch`。在private-teodoro-gpu请求1 RTX3090／24GiB VRAM／2CPU／16GiB／60分钟，generanno环境；原1536 draw×3＝4608更新，仅末step004608科学读出。先生成不可变FP32隐藏缓存并在同八窗口核对live／cache分数、loss、裁剪梯度及梯度有无状态，再优化；该缓存阶段不是已完成训练。拟合结束后从磁盘严格装载末权重，在同缓存上读取全部1491自然窗口。

CPU核验root13390189，命令为 `sbatch --parsable --dependency=afterok:13390188 sbatch/M28-A21-JUNCTION-CHAIN-RESULT-CHECK.sbatch 13390188`。private-teodoro-gpu／0GPU／2CPU／8GiB／10分钟，初始有效Dependency等待，不改变阈值或进行新模型前向。原preflight依赖已满足，Slurm可清除其依赖字段；保存的SubmitLine仍记录afterok:13390039。

首次运行快照确认13390188 RUNNING／gpu034／elapsed2:44，FIT与TRAIN-READOUT的RUNNING及JOBID均匹配；日志到immutable_hidden128／1491窗，stderr无错误。这是冻结隐藏缓存 materialization，不是128次更新或训练完成。13390189 PENDING／有效Dependency。

拟合和读出目录分别为 `outputs/M28-A21-JUNCTION-CHAIN-FIT`、`outputs/M28-A21-JUNCTION-CHAIN-TRAIN-READOUT`，二者JOBID都为13390188；核验为 `outputs/M28-A21-JUNCTION-CHAIN-RESULT-CHECK`，JOBID13390189。只有两root均COMPLETED0:0、三个目录完成标记和全部约定产物齐全，经真实4608条更新、1491条几何／标签／分量／选择／profiles／原五门槛核验后，主对话才解释科学结果。单fit summary、提交、缓存完成或运行中状态不等于该结果。

提交前实时Slurm显示private可用RTX3090，gpu034／gpu035为24G；shared没有更早可信完成优势，故保留private路由。当前并发远低于80，GPU维护保留未覆盖此执行窗，BeeGFS可用约136TiB。本轮不取消或修改其他作业、不增加seed／epoch、不访问DEV或sealed、不提交Git或公开发布。

## 固定等待交接

同一固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3、同一heartbeat coding-rna-r5-light-monitor已更新为 `M28-A17-A21-TRAIN-20261005`，只包含A17三roots13387910／13387946／13387978和A21两roots13390188／13390189。native更新成功，完整prompt、heartbeat类型、ACTIVE、固定target、15分钟间隔、failed_runs_only和原created_at实际保存读回一致；新手动启动的真实turn_context为gpt-5.6-luna／low，主模型未改变，没有新对话、fork、cron或第二自动化。

新范围的手动启动取证于15:19:35.630 UTC完成：A17 fit13387910 COMPLETED0:0、readout13387946 RUNNING、CPU check13387978有效Dependency；A21合并root13390188 RUNNING、CPU check13390189有效Dependency。实际同一监控turn01a10ca5-d7f3-79a1-8752-2142acd3225f completed／error-null／37.070秒，未发健康状态消息；原两条流程都尚未闭合，不解释指标。

新范围的15:28:23.473 UTC定时快照仍为健康运行／有效依赖；同一监控turn01a10cad-edc1-7073-b845-0420cccf99a8于15:29:57 UTC completed／error-null／109.738秒，未发送健康消息。随后15:43:32.346 UTC定时取证确认A21两root成功且白名单齐全、A17仍未完成，向原主对话成功交接一次。同一监控turn01a10cbb-a96d-72b2-9baf-a26abc60fd94于15:45:54 UTC completed／error-null／166.380秒；原心跳实际TOML为PAUSED，updated_at对应15:45:45.006 UTC。配置、启动、定时追加及本轮终态回传均已验证；主对话已核验实际科学FAIL，后续只恢复同一心跳中的A17。

主对话消费A21后直接恢复原同一心跳为 `M28-A17-TRAIN-REMAINING-AFTER-A21-20261005`，只保留A17三root。15:51:39 UTC实际保存的完整新prompt、ACTIVE、heartbeat类型、固定target、15分钟间隔、failed_runs_only及原created_at均读回一致。手动启动15:52:20.469 UTC确认A17 fit已完成、readout仍运行、CPU check有效Dependency；同一监控turn01a10cc3-f275-7601-a4a9-054a340c1420于15:52:31 UTC completed／error-null／20.510秒，真实turn_context仍gpt-5.6-luna／low，无健康回传。缩小范围的定时运行及A17终态回传尚未观察，不能沿用A21回传代替。

完整流程终态或可行动故障只向原研究主对话交接一次，再暂停同一心跳。正常运行／有效依赖静默，监控不读取性能、解释科学结果或改变任务。长期本机等待需要电脑开机、应用运行和连接可用，不能保证休眠时即时检查，见[官方定时任务文档](https://learn.chatgpt.com/docs/automations?surface=app)。
