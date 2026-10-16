# M28 A首轮执行记录

2026-10-03，Asia/Bangkok更新。用户已批准本轮单一A拟合/完整DEV评价及结果处理，并取消成本门槛。三个正式阶段均COMPLETED/0:0，4608更新、8690方向窗及评价完成；A在两物种的完整链F1均改善，但背景span都增加，原pilot继续条件未通过。完整结果、代价与下一科学问题见[a_result.md](a_result.md)。原科学配置和数据隔离未改；本轮不是新的可发表模型或独立泛化成功。

## 当前执行与固定科学范围

Q16/H128/seed0、四个loss系数1、原1536 paired TRAIN draws按原顺序三遍、4608更新、R4目标、tau0/max_exons128及唯一step004608主评价checkpoint保持不变。C0/B1和原GLM feature cache只读复用；主参考5437+2291=7728，旧6450仅副表。两物种均比完整C0严格提高exact-chain F1且背景span都不增加，才通过pilot继续条件；不是论文成功或独立泛化证据。

| 阶段 | 作业 | 命令入口 | Slurm资源/单作业时窗 | 终态 |
|---|---|---|---|---|
| 固定TRAIN计时 | 13331164 | run_query.py --mode prepare | 1 RTX3090 24GiB、2CPU、16GiB、10分钟 | COMPLETED/0:0，历史成本gate=false |
| 唯一正式拟合 | 13344559 | run_query.py --mode fit | 同GPU/CPU/memory、4小时 | COMPLETED/0:0；2606分配秒 |
| 完整DEV推理 | 13344568 | run_query.py --mode infer | 同GPU/CPU/memory、6小时 | COMPLETED/0:0；5033分配秒 |
| 统一评价 | 13344571 | evaluate_query.py | 0GPU、2CPU、32GiB、2小时 | COMPLETED/0:0；84分配秒 |

上表为2026-10-03结果处理时fresh sacct及STATUS/JOBID核实的终态。提交时RUNNING/Dependency快照已由终态取代，不以旧状态代替调度事实。当前walltime为调度设置，不是经济停止/审批门槛。原10.166667GPUh、4独立CPUh、2GiB人为合计上限取消，仍记录实际资源。此前合成工程210CPU秒、21553 bytes、0GPU单列保留。本轮不运行原生FASTA，不开放test/Setaria，不重训控制，不扫参数，不发布或提交Git。

## 取消成本阻断的具体改动与验证

run_query.py的validate_preparation保留completed、CONFIG、0 optimizer、revision和八类/八窗口身份，但不再要求历史budget_gate.passed=true；artifact_bytes只报告实际大小，不再因人为2GiB中止。prepare的旧预算报告函数和全部原计时文件保持不变，未改写历史false。没有改模型、损失、数据、draw顺序、更新数、checkpoint选择或科学评价阈值。

tests/test_m28_query_runner.py只新增/调整对应这两项变更的测试：成本gate=false或缺失时仍接受正确准备身份，不接受不完整/错误窗口数/已更新/错revision/配置；模拟3GiB产物不再人为中止。FIT脚本先执行这两项测试，再进入真实fit。远程精确回读、Python AST与FIT bash -n通过；FIT stderr显示两项测试均ok，Ran 2 tests in 0.036s / OK。未重跑未变旧测试。

初期FIT日志的step64（41.304秒）和step128（73.013秒）仅为当时已开始optimizer的证据。现在完整summary确认4608更新，唯一主评价末checkpoint step_004608.pt为7044788 bytes；TRAIN loss下降不作DEV准确率或模型成功判断。

## 本轮调度证据

提交前已读远程规则、cluster_config.yaml及当前执行脚本，沿用实际generanno环境。五项实时Slurm检查：private gpu034的RTX3090配置8、已用7，gpu035满载；用户当时无活动作业，未见影响此GPU cohort的预约。私有可提供适配资源，故三阶段均private-teodoro-gpu。BeeGFS可用176530879676416 bytes；三个后续输出目录当时不存在。保留已有工作树修改，不恢复退役框架。

sbatch --parsable sbatch/M28-A-FIRST-FIT.sbatch返回13344559，立即观察RUNNING；sbatch --parsable --dependency=afterok:13344559 sbatch/M28-A-DEV-INFER.sbatch返回13344568；sbatch --parsable --dependency=afterok:13344568 sbatch/M28-A-EVALUATE.sbatch返回13344571。三ID逐一核实，未重复提交或取消作业。非有限数、错配、N>Q、错误参考读取等仍停止受影响阶段并保留证据；实际硬件/调度失败按具体原因处理，未完成不能算科学FAIL或进行正式排名。

## 完成后的只读query诊断

job13347285完成全部已保存A方向窗的同窗不同链/重复query计数，private-teodoro-gpu、1CPU/2GiB、0GPU、5分钟调度时窗，实际root COMPLETED/0:0、2分配秒、STATUS=COMPLETED、stderr为空。只读取M28-A-DEV-INFER，不读取参考或模型，不重新解码、不修改候选/score/阈值，输出新目录outputs/M28-A-QUERY-DIVERSITY。新脚本AST与sbatch bash -n、精确回读通过；计数对齐原inference完成字段。75.4430%的正gain完成query为同窗重复链，具体统计及科学边界见[a_result.md](a_result.md)。这是原授权的结果诊断，不是新拟合；因短作业当轮即完成，没有重启等待心跳。

[A2固定协议](a2_query_collaboration_proposal.md)据此提出唯一query协同结构变化，用户于2026-10-03明确批准并要求不再重复询问；新增实现和9项工程检查通过，固定拟合→完整DEV推理→评价已启动，见[a2_execution.md](a2_execution.md)。A2科学收益仍待完整结果，不变更本轮A结论或既定背景标准。

## 固定监控接续

固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3（local）的真实定时turn 01a0fdc9-a037-7962-ac66-3094e9be1269已核对三阶段完成字段，并向原主对话成功交接一次。本轮实际追加到原监控对话的定时检查及自动回传已观察；turn_context确认gpt-5.6-luna/low，不改变主对话模型。主对话随后独立核实终态并完成科学结果处理。

完成心跳按应用内置规则再次删除了自动化及状态目录，而非用户偏好的暂停保留。主对话通过官方工具恢复同一自动化ID coding-rna-r5-light-monitor为PAUSED，同target及15分钟间隔，保存字段回读一致；保留同一监控对话、恢复本cohort Handoff sent YES记录，不新建监控对话或cron，不重复交接。未修改应用或绕过其指令，未来自动暂停保留仍有既知产品限制。当前无等待对象，下一批仅更新同一对话/心跳的真实作业对象。

## 历史准备计时与已取消的成本停止

13331164 root及batch均COMPLETED/0:0，STATUS=COMPLETED。实际110秒×1GPU=0.030556 allocated GPUh，110秒×2CPU=0.061111 allocated CPUh；batch MaxRSS1942304KiB，GPU峰值allocated307266048/reserved413138944 bytes。八固定类别/八唯一窗口、optimizer_steps=0，DEV_forward_used=false、test_setaria_used=false；日志及产物合计70072 bytes，stderr为空。这只是生产接线证据，不是训练或原生FASTA验证。

首窗load→backward=11.306037477秒（load1.384039382、shared forward4.288843201、matching/loss2.430062816、backward3.202953228）；其余七窗step为1.207379479、0.423808973、0.541426679、0.705770499、0.671879637、1.002185567、0.460110951秒，shared forward范围0.066187459–0.068093879秒。八窗循环16.875896457秒，环境启动包含于Slurm110秒。

历史公式max(full step)×4608+1800得到53898.220694080秒=14.971728小时，超过旧4小时；max(decode)×8690+1800=3143.411614457秒=52.390194分钟。历史budget_gate=false使当时停止，原summary及逐窗证据保持原值。首shared forward约后续63倍强烈提示一次性冷启动被重复外推，但未同窗warmup确认，不保证实际fit/infer时长。

旧文档提出的PREPARATION-R2 warmup复测与新计时公式从未执行；现被用户取消成本门槛的决定取代，无须通过复测得到“预算PASS”。上述三阶段现已完成，不再存在计时规则批准待办。正式科学判定见[a_result.md](a_result.md)：完整链改善，但背景条件未通过；不是旧成本gate阻断。原计时文件不改写。
