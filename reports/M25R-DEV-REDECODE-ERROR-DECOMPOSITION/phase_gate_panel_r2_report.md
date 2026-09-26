# M25R epoch-1 phase 面板 R2：等距抽样 CPU 预检

日期：2026-09-08。**预检完成；158个原方向化窗口超过64上限，GPU消融未执行。** 这是执行预算阻塞，不是phase机制的科学失败。

## 授权与实施

本轮用户批准：仅修订为层内等距抽样、重新冻结128条lineage、执行一次CPU预检；资源上限2CPU、16GB、10分钟，暂不提交GPU。

保留epoch1、row601及64 reference-exact phase failures / 32 reference-nonexact phase failures / 32原emitted controls。每层仍使用原largest-remainder配额；exact组冲突标签仍为**原首次phase失败位点**，未改成any-start冲突。species/seqid/strand/lineage_id排序后按 `floor((j+0.5)*N/n)` 选择，不使用scores、链长或窗口覆盖优化。等距覆盖有序列表减少前缀的位置偏向，但不是概率随机抽样，仍不能外推总体性能。

新增 `phase_panel_systematic.py` 和3项聚焦测试。原 `phase_panel_preflight.py::main` 仅增加显式output/selector/selection参数，默认前缀行为保留；R2复用原候选池、对账及窗口核算，不复制一套算法。另记录CDS数、已存phase-check数、chosen-chain跨度和必要score位置距tile边界≤6bp的描述标记，均不参与选择。

使用coding-rna-baobab Skill，在Baobab源项目通过Slurm执行；当前AGENTS与cluster_config已核实。未恢复旧框架、未改他人工作、未提交Git、未调用新的评阅、未读FASTA或checkpoint内容，更未读取Setaria。

## 作业与验证

- 命令：`sbatch --parsable sbatch/M25R-E1-PHASE-PANEL-R2-PREFLIGHT.sbatch`
- Job `12520106`，private-teodoro-gpu，gpu035，2CPU、16GB、10分钟上限；AllocTRES无GPU，Requeue=0。
- Slurm COMPLETED，exit0:0，2026-09-08 16:28:14–16:28:19 Europe/Zurich，用时5秒，MaxRSS164016K。
- 8项测试通过，耗时0.09秒：原5项核算回归加3项等距选择测试；Bash syntax通过。无作业重投。
- 128条唯一lineage与固定候选池1149/2129/2098对账通过；manifest先以exclusive-create落盘，再核算窗口。不因超限重新选择。

## 结果

| 组别 | lineage数 | 原方向化窗口数 |
|---|---:|---:|
| reference-exact、phase拒绝 | 64 | 81 |
| reference-nonexact、phase拒绝 | 32 | 40 |
| 原输出对照 | 32 | 37 |
| 去重合计 | 128 | **158** |

方向分解：Arabidopsis +=56、-=63；rice +=19、-=20。当前三组窗口恰无跨组重叠，组内已去重。旧前缀面板需137，新等距面板需158；二者不能混用。仅exact组81也超过64，因此不是移除少量对照就能保持原实验解决的问题。

窗口核算仍包括全部motif竞争候选分数、chosen边界和所有chosen CDS starts（含acceptor+2），按原6144bp方向化tile映射。ORF只需真实sequence，不增加推理窗口。region/candidates保持冻结，核算并非完整pipeline重放。

状态 `STOP_WINDOW_BUDGET_EXCEEDED`。GPU作业数=0，forward数=0，训练数=0。原路径重放、phase旁路、ORF/threshold通过计数均为NOT_RUN，不填0或猜测。32/64与8/32继续仅为探索性参考，不能解释为路线否决标准。

## 产物与下一步

独立目录：`outputs/M25R-E1-PHASE-GATE-PANEL-R2/`，包含frozen_panel.json、window_coverage.json、preflight.json、PREFLIGHT_JOBID、STATUS和logs。R1 `outputs/M25R-E1-PHASE-GATE-PANEL/` 及全部原R4/勘误/回溯结果保留。

**建议下一步只调整已查明的资源合同，不再换样本：**若继续同一128条等距面板，需要明确批准将窗口上限从64提高到158（970752bp原窗口输入），其余1GPU×2小时、8CPU、32GB、零训练、原路径先完全重放一致、失败/超时停止不重投保持不变。该预算能否容纳推理仍需提交前保守估算；本CPU预检未实测GPU速度，不能保证2小时。批准前不实现/提交GPU实验。

若64窗口不可改变，则需要用户另行选择更小样本或不同的窗口簇诊断设计；那是新设计，不自动截断当前面板。本轮授权工作完成，到此停止，不自动重训或扩大实验。
