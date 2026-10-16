# M28 A2查询协同执行记录

2026-10-03，Asia/Bangkok。用户明确批准[A2固定协议](a2_query_collaboration_proposal.md)并要求同范围不再重复询问。唯一结构已实施，9项新合成接线检查全部通过；单次4608更新拟合已完成，初次DEV推理因2小时调度时限中断。同末权重补跑及原评价作业均已COMPLETED/0:0，8128复用+562新窗口=8690。完整计数复算一致，但两物种F1均低于A、背景span均高于C0，预先规定的四项科学条件全部未通过；不是可发表模型完成。完整结果见[a2_result.md](a2_result.md)。原A结果及冻结资产不变，不提交Git或对外发布。

## 实施与固定范围

新增src/m28/query_collaboration.py中的CollaborativeQueryWindow、CollaborativeQueryPointer与CollaborativeGeneQueryModel。旧A模块按原顺序初始化，然后新增LayerNorm与四头自注意力；q=LayerNorm(initial_A)，initial_A2=initial_A+Attention(q,q,q)，dropout0，无FFN。存在性头及全部teacher/自由指针初始状态共同使用initial_A2。训练入口保留shared梯度，原始外部feature路径仍detach；训练和推理都在所有query上通信，不读取参考决定通信对象。

旧query_model.py、query_pointer.py、query_training.py、query_decode.py、query_output.py、旧A运行/评价器没有改。A2独立入口run_query_collaboration.py复用原TRAIN身份、R4目标重建、损失、输入和无参考输出函数；evaluate_query_collaboration.py复用原参考、组装、GFF3回读与公共计量实现。不存在通用实验注册框架、可选开关或新的标签/阈值策略。新模型649457参数，比A增加66304；不是等参数对照。

保留Q16/H128/seed0、D2048与原revision、六TRAIN染色体、原1536 draws及顺序/importance weights三遍、4608更新、batch1、四loss系数1、原AdamW/clip/FP32。复用A已完成的TRAIN准备身份，不重跑成本计时或控制拟合；实际fit仍按原函数核对1536 draws/1491唯一窗/物种各768等身份。唯一主评价step004608、tau0、max_exons128、完整8690方向窗，原DNA语法、owner、去重、正gain非重叠选择不变。test/Setaria不开放。

## 针对新机制的工程验证

工程job13348629，private-teodoro-gpu，1CPU/4GiB、0GPU、5分钟时窗；root、batch和extern均COMPLETED/0:0，STATUS=COMPLETED、JOBID=13348629。9项测试在1.043秒通过，Slurm分配68秒×1CPU，batch MaxRSS791320KiB，无真实TRAIN/DEV、原checkpoint或optimizer更新。

测试逐张量核对完整D2048/Q16/H128 seed0旧A参数起点；验证残差pre-norm公式与置零新增输出层恢复原状态、跨槽位通信、未知logit无直接监督但表示有间接梯度、新attention/LayerNorm/shared获得有限非零梯度、Hungarian仍达到穷举全窗目标最低值，以及训练/自由DNA解码使用相同协同状态。另检查A2末权重身份拒绝A/错步数、完整DEV网格身份和F1对A/背景对C0的整数gate。未重跑未变旧测试。

4个新Python文件AST通过，4份sbatch分别bash -n通过；远程git apply --check后应用，8个新源码/测试/sbatch逐文件精确回读一致。未修改既有未提交工作。以上只支撑实现接线，不证明准确率、新颖性或泛化。

## 初次作业及当前恢复

| 阶段 | 作业 | 固定入口 | 资源与walltime | 当前已核实状态 |
|---|---|---|---|---|
| 合成接线 | 13348629 | test_m28_query_collaboration.py | 0GPU、1CPU、4GiB、5分钟 | COMPLETED/0:0 |
| 唯一拟合 | 13348639 | run_query_collaboration.py --mode fit | 1 RTX3090 24GiB、2CPU、16GiB、90分钟 | COMPLETED/0:0，4283分配秒 |
| 完整DEV推理 | 13348644 | run_query_collaboration.py --mode infer | 同GPU/CPU/memory、2小时 | TIMEOUT/0:0，7227分配秒 |
| 同末权重补跑 | 13353734 | resume_query_collaboration.py | 1 RTX3090、2CPU、16GiB、1小时 | COMPLETED/0:0，399分配秒 |
| 统一评价 | 13348647 | evaluate_query_collaboration.py | 0GPU、2CPU、32GiB、20分钟 | COMPLETED/0:0，35分配秒 |

三次sbatch --parsable分别返回上述ID，立即squeue核对无重复提交。2026-10-02 19:08 UTC前后取得的启动快照：拟合STATUS=RUNNING、JOBID一致，TRAIN日志已到step192；64/128/192的循环累计秒为47.7358/80.4724/113.4389，stderr当时为空。这些只是已发生optimizer更新的证据，不用于DEV排名、收敛结论或保证完成时间。

新目录固定outputs/M28-A2-QUERY-ENGINEERING、M28-A2-FIRST-FIT、M28-A2-DEV-INFER、M28-A2-EVALUATION；产物身份arm=A2、GFF3为A2.gff3。提交前重新取证五项Slurm调度/预约：gpu034配置8张RTX3090、已用7，无影响本轮的GPU预约；适配私有资源可用，选择private。其他TE项目作业未触碰。BeeGFS可用176113480368128 bytes，三个正式输出目录不存在，无覆盖。没有经济成本gate或人为总GPUh/2GiB停止门槛；walltime仅为调度设置。

## 评价与等待接续

主分母7728=5437+2291、历史6450副表，unsupported/长链/重叠参考保持。A2两物种exact-chain F1均严格高于A、背景span均不高于C0才通过本轮继续条件，按整数计数比较；完整CDS-base union及所有背景指标照实报告。A/C0及其他工具只读缓存，不称最新原生重现。

新评价器在同一遍原自由候选读取中记录每物种同窗重复query，和已完成A只读冗余结果比较；同时比较free→owner→selected精确参考覆盖。未知标签与间接协同梯度区别已记录；机制诊断不替代主gate。DEV是开发数据，不是独立确认性证据。未通过不事后改tau、Q、seed、epoch或checkpoint救本轮；实际运行失败与科学FAIL分开。

以下为初次启动时的历史快照：固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3（local）和同一automation ID coding-rna-r5-light-monitor复用，heartbeat/15分钟/ACTIVE，target与完整prompt保存回读精确一致。仅当前cohort M28-A2-13348639-13348644-13348647，交接标记重置NO，不重查已关闭A/B2。人工启动检查实际turn 01a0fe05-0ab0-7c33-a80d-300069aac956完成，核实fit RUNNING gpu034、infer/eval PENDING Dependency，未向主对话发进度交接；实际turn_context核对model=gpt-5.6-luna、effort=low，不改变主对话。该人工检查不等于本cohort真实定时追加/自动回传已经观察；完成/失败需处理才一次交接，由主对话验证及解释。

完成心跳可能按应用内置规则删除自动化而非用户偏好的暂停保留，这一既知限制没有用应用修改或指令绕过解决。固定监控对话保留，必要时通过官方工具恢复同一自动化ID/target。定时等待依赖本机开机、应用运行及连接可用，不保证休眠/断线时即时唤醒。当时完整结果验证、科学解释及本轮记录尚待完成；最终关闭状态见末节。

## 2026年10月3日同权重恢复

实际定时turn 01a0feb7-acf4-73f3-9b72-51b7dd0b7782于2026-10-02T22:24:05Z至22:24:58Z核实原推理TIMEOUT、评价DependencyNeverSatisfied并一次交接主对话；不是科学gate失败。主对话fresh sacct确认fit13348639 COMPLETED/0:0，optimizer_steps4608，末checkpoint step_004608.pt为7847716 bytes；原推理13348644 TIMEOUT/0:0，stderr明确“2026-10-03T00:17:27 ... DUE TO TIME LIMIT”（Baobab CEST）。原推理STATUS仍残留RUNNING，不能代替调度终态；原文件/日志保留，不改写为完成。

恢复入口resume_query_collaboration.py只针对这个已确认的13348644源目录。逐记录核对有序window ID/index、无参考来源、Q16/raw query、score及counts，并重建完整网格。无换权重、重训、追加epoch、调tau或科学评价选择；新产物写outputs/M28-A2-DEV-INFER-RECOVERY，不覆盖原四scope文件。

每个已有方向scope用同一末checkpoint回放首窗，要求离散链/终止原因完全一致、query logit最大绝对差不超过1e-5、正gain状态不变；不匹配就停止，不把不同输出混成完整推理。只复用生产者已写完换行的记录；可能截断的末记录不算完成，重新生成。2026-10-02 22:35 UTC前后恢复日志已确认8128复用窗口，新窗口至少64个；最终缺失总数562。计数不重复加入4个接线回放窗，具体完成数量及回放差异由最终summary记录。

4项新恢复测试在job13353644通过：完整/截断尾记录、错身份/重复/完整坏JSON拒绝、回放离散/logit/正gain差异拒绝、计数重建。COMPLETED/0:0、STATUS COMPLETED，4 tests/0.008s，Slurm64分配秒×1CPU、0GPU。新sbatch分别bash -n通过，两个新Python源码精确回读一致；未重跑未变模型的旧9项测试。新增恢复不恢复Auto Research或修改goal。

补跑job13353734启动时RUNNING gpu034，最终COMPLETED/0:0，private-teodoro-gpu、1 RTX3090 24GiB、2CPU/16GiB、1小时操作时窗；这不是成本门槛。旧fit保留，仅补缺失方向窗。scontrol update JobId=13348647 Dependency=afterok:13353734成功，回读为afterok:13353734(unfulfilled)，随后PENDING Dependency；不取消旧评价、不创建第二个评价作业。evaluate_query_collaboration.py只将输入路由到恢复目录，并把新增实际产物大小计入记录，冻结计量/参考/gate/输出目录不变。

提交前fresh五项Slurm与存储核实：gpu034/035均配置8 RTX3090、各已用3，没有影响本轮GPU预约；BeeGFS余量约161TiB。补跑/恢复测试和最终评价目录无碰撞。原训练/失败推理的资源与产物单列，恢复时长不伪装为原完整推理耗时。

旧心跳在超时交接后再次按内置规则删除。官方工具恢复同一automation ID coding-rna-r5-light-monitor、同target01a0f71f-8845-7432-bd13-bf54c1808fb3，ACTIVE heartbeat/15分钟，完整保存配置精确回读；新cohort M28-A2-RECOVERY-13353734-13348647交接标记NO，旧TIMEOUT事件已处理不重复发送。人工启动turn01a0fec2-7224-7c22-a104-4b0cafc25ed8核实fit完成、补跑RUNNING、原评价PENDING Dependency，实际turn_context确认gpt-5.6-luna/low。该人工检查不代表定时回传；随后同一监控对话的实际定时检查完成交接，最终状态见末节。

## 最终完成与科学负结果

2026-10-02 22:40:24 UTC（2026-10-03 05:40:24 Asia/Bangkok）原评价job13348647结束。fresh sacct核实fit13348639、补跑13353734、评价13348647均COMPLETED/0:0；对应STATUS及summary完成字段一致，infer/eval stderr为空。原推理TIMEOUT及残留STATUS保留，不伪装成功。

恢复summary最终为8128 reused+562 new=8690，四个首窗回放discrete_equal=true、max_abs_query_logit_difference=0.0，无截断尾记录；源文件未变、模型/操作点未变，step004608仍唯一权重。原评价作业复用成功，无重训或第二评价。最终A2.gff3非空，7507964 bytes。

主对话独立复算整数P/R/F1、背景、宏/微口径和阶段守恒均与summary一致。A2拟南芥/水稻完整链F1为50.8991%/55.0203%，均低于A的53.3975%/56.3767%；背景span覆盖1.6088%/1.7952%，均高于C0的0.6698%/1.0026%。四项主要条件均false。机制诊断改善不能替代此负结果；完整原始计数、CDS-base代价和free→owner→selected分析见[a2_result.md](a2_result.md)。

固定监控对话的真实定时turn01a0fecf-7ac5-70e3-99fa-d2dbe8ef91a2于2026-10-02 22:50:05–22:51:06 UTC完成检查并成功向原主对话一次接续，实际turn_context为gpt-5.6-luna/low。已观察本cohort真实定时追加与自动回传，不仅是配置/人工检查成功。应用完成规则再次删除心跳；主对话通过官方工具恢复同一ID coding-rna-r5-light-monitor、同target01a0f71f-8845-7432-bd13-bf54c1808fb3为PAUSED，15分钟间隔。旧批次不继续轮询，不新建监控对话或重复交接。

本轮授权内的实现、一次固定拟合、完整DEV推理、验证、解释及记录已完成；论文新模型目标尚未达到。保留原A与A2负结果，不事后扫阈值、改checkpoint、追加训练、开放test/Setaria或对外发布。
