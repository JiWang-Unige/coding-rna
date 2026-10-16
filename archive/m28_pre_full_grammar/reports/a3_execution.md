# M28 A3已知query平均损失执行记录

2026-10-03，Asia/Bangkok。用户在收到A2负结果后明确要求“继续尝试”。已按[a3_known_query_mean_protocol.md](a3_known_query_mean_protocol.md)完成新目标实施与8项合成检查，并启动一次固定拟合、完整DEV推理及评价。截至本轮终态核对，拟合、完整DEV推理、评价与关键计数复算均已完成；预先规定的科学继续条件失败，详见[a3_result.md](a3_result.md)。工程完成不等于科学PASS或论文模型完成。

## 实施和单一变化

新增src/m28/query_known_mean.py，存在性损失按窗内全部已知query平均，不再按正负两类等权平均。Hungarian配对成本与完整窗目标一致：label0替换正链时使用正BCE减负BCE再除以Q；label-1时只用正BCE除以N。匹配后重新计算带梯度的原teacher损失。边界/EOS/local口径及未知政策不变。

新入口run_query_known_mean.py复用原TRAIN身份/目标构建、A2模型类和完整seed0初始化；Q16/H128、649457参数、1536 draws三遍/4608更新、原importance、四loss系数1和AdamW/FP32保持。新增configuration.existence_reduction=mean_over_known_queries。没有新头、新标签、阈值扫描或监督checkpoint继承。

evaluate_query_known_mean.py保留8690自然网格、7728主参考、原GFF3回读、组装和计量；新增A2已完成缓存对照，核对配置除目标标记外相同、参数与更新数相同、每遍TRAIN已知/未知监督计数相同。主要gate仍两物种F1严格高于A、背景span均不高于C0，原始整数比较，不以A2弱结果降低门槛。

所有旧A/A2/B2源码、checkpoint和结果保持不变，未提交Git或推送。新文件经远程git apply --check后应用；4个新Python文件AST和4份sbatch bash语法通过，9个源码/测试/脚本/协议文件逐字符远程回读一致。只新增本轮文件，不改无关未提交工作。

## 新目标相关检查

job13357559及root为COMPLETED/0:0，STATUS COMPLETED，8项测试1.211秒；分配82秒×1CPU、4GiB、0GPU，private-teodoro-gpu/generanno。没有真实TRAIN/DEV前向、optimizer更新或参考注入。

通过不等类数和全正类的Q4/N1、N2、N4穷举，确认Hungarian达到新完整窗目标的最低值；零logit时单正三负的直接梯度分别为-1/8和+1/8，而旧按类平均为-1/4和+1/12。未知直接logit梯度为零，但共享context/协同层可通过已监督query接收间接梯度。空正例的已知负窗与全未知窗保持不同，容量不足不删正例。

还核对真实训练入口调用新损失、R4目标和原importance未绕过；A3完整初始化与A2逐张量完全相同，参数649457；错误A2/旧目标/中间步数和部分DEV网格被拒绝。未重跑不变的旧模型/解码测试；这些检查不证明性能。

## 调度与当前状态

提交前fresh五项Slurm检查显示gpu034/035各配置8张RTX3090、已用3/1，私有合格资源可用，没有影响本轮的GPU预约；选择private-teodoro-gpu，不触碰他人作业。BeeGFS可用175930158350336 bytes，四个A3目录均不存在。没有经济成本、总GPUh或2GiB门槛；真实资源/维护/输出碰撞限制仍遵守。

| 阶段 | Job ID | 入口 | 资源和操作时窗 | 启动取证 |
|---|---:|---|---|---|
| 合成目标检查 | 13357559 | test_m28_query_known_mean.py | 0GPU、1CPU、4GiB、5分钟 | COMPLETED/0:0 |
| 唯一拟合 | 13357566 | run_query_known_mean.py --mode fit | 1 RTX3090 24GiB、2CPU、16GiB、90分钟 | RUNNING gpu034 |
| 完整DEV推理 | 13357568 | run_query_known_mean.py --mode infer | 同GPU/CPU/memory、3小时 | PENDING Dependency，afterok:13357566 |
| 统一评价 | 13357569 | evaluate_query_known_mean.py | 0GPU、2CPU、32GiB、20分钟 | PENDING Dependency，afterok:13357568 |

三次正式sbatch各返回上述ID，逐次立即squeue核对；scontrol实际回读两条afterok依赖正确，未重复提交。正式提交约2026-10-03 06:43 UTC，即13:43 Asia/Bangkok，Slurm提交时间08:43为CEST。

06:45 UTC前后主对话核实fit STATUS RUNNING，日志已到step192，64/128/192循环累计秒为36.3787/60.7903/85.7551，stderr为空。这只是实际optimizer更新和有限训练分量的启动证据，不用TRAIN loss宣称收敛、预测质量或保证结束时间。

新产物独立位于outputs/M28-A3-QUERY-ENGINEERING、M28-A3-FIRST-FIT、M28-A3-DEV-INFER、M28-A3-EVALUATION，log沿用outputs/M28-PILOT-R5/logs/M28A3..._jobID。推理3小时操作余量参考A2原2小时中断，不改变模型、科学停止条件或允许额外epoch。

## 固定轻量等待

复用同一监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3（local）与同一automation ID coding-rna-r5-light-monitor，heartbeat/ACTIVE/15分钟，完整prompt、目标对话、状态和间隔保存回读精确匹配。只监控cohort M28-A3-13357566-13357568-13357569，Handoff sent NO；已关闭A2事件保留但不继承其YES或重复发送。

人工启动turn01a10082-7beb-7cc0-8248-8bcc227db682于2026-10-03 06:45:13–06:45:39 UTC完成，实际turn_context核对gpt-5.6-luna/low。它在同一固定对话核实fit RUNNING gpu034、infer/eval PENDING Dependency并更新自身记录，没有向主对话发进度交接。配置和人工单次检查均成功；本A3 cohort的真实定时追加与完成回传尚未观察，不用旧A2成功代替。

正常RUNNING/PENDING保持静默；空squeue必查sacct，整链完成并且末权重/8690无参考推理/7728评价/非空GFF3完成字段齐全，或失败/超时/依赖不满足/缺产物时一次回传。主对话继续验证、解释及本轮结果记录，监控不能重提/取消/延时/改参数/访问封存数据或判断科学PASS。

等待依赖本机开机、应用运行和SSH连接，不保证休眠/断线即时回传。应用完成心跳可能按内置规则删除而非暂停保留；固定对话仍复用，必要时经官方工具恢复同一ID/target，不改应用或绕过规则。当前没有平台goal循环或第二套轮询。

## 已完成结果与接续

13357566/13357568/13357569均为COMPLETED/0:0，分配1876/3614/32秒；全部STATUS完成，唯一step004608和完整8690窗身份一致，7728共同尺与非空GFF3已核实。评价终态2026-10-03 08:15:48 UTC。最终拟南芥/水稻F1为51.7567%/53.5236%，背景span157687/368489bp；两物种F1>A和背景≤C0四项全部false。结果、机制代价、计数复算和历史工具差距见[a3_result.md](a3_result.md)。

固定监控对话在08:15同一对话进行了真实定时检查，当时评价仍RUNNING，未发送交接。其后主对话直接消费完整终态，不伪称自动回传成功。复用同一对话和同一heartbeat，已移除完成A3与工程smoke，当前仅等待另行固定的full-TRAIN五作业；不开新监控对话或cron。

本轮负结果保持，不改tau、checkpoint或追加训练。已启动的full-TRAIN Q/C0保持其独立协议和新匹配C0 gate；提交或工程检查通过不等于结果完成。不开test/Setaria，不自动独立迁移或发布，论文目标仍未达到。
