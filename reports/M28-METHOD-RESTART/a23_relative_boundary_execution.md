# A23 相对边界结构执行记录

2026-10-06。正式固定TRAIN链已完成，原三项生成覆盖gate全部通过；错误候选负担同时增加，完整caller质量与论文竞争性未证实。科学选择保持首次固定版本，不增加seed、更新次数或阈值搜索。完整结果见[a23_relative_boundary_result.md](a23_relative_boundary_result.md)。下文提交和启动状态为历史取证。

## 终态

13391608／13391609／13391610三根及各batch、extern均COMPLETED0:0，实际1:04:03／43:36／13:58，STATUS和JOBID匹配。完整4608曝光、1491窗自然池与CPU原三gate核验闭合。合并正确候选2380→2415，≥5外显子788→832，known_nonmatch12566→15027。保留A23，不将候选覆盖PASS等同最终精度。固定监控已定时回传，主对话已消费全部终态；同一心跳暂停等待下一实际作业，不再监控这些已完成对象。

## 实际预检

作业13391430，private-teodoro-gpu／1 RTX3090 24GiB／2CPU／16GiB／10分钟上限；root、batch、extern均COMPLETED 0:0，实际38秒。输出M28-A23-RELATIVE-BOUNDARY-PREFLIGHT，STATUS与JOBID相符。四项模型／padding／原始DNA greedy与beam检查全部通过，用时1.030秒。

八项真实TRAIN预检12.869889秒，optimizer0，所有模型参数未变，共有A3 fresh seed0逐tensor相同，初始四loss与离散匹配逐case一致。模型649617参数，新增160。两物种都存在end和next_start新输出列非零梯度；只有单外显子的拟南芥most_positive case next_start不适用，实际梯度为0，未错误要求其非零。无正链窗口两列均为0。末层零初始化导致输入层首步零梯度，符合预声明计算；正式训练将核验其更新后实际梯度，而不是事后改变初始化。

预检共有分支梯度有限且实际接通。八项计时包含两个模型的前向、监督及一次反向，不能据此声称单模型训练速度或生成覆盖改善。无新GLM前向、DEV、sealed或Setaria访问。

## 固定正式执行

run_query_relative_boundary.py沿用原A3 loss与TRAIN loader，fresh一次、1536 draw×3、4608更新，末004608唯一主读出，每256步工程保存。readout先在全部1491窗关闭raw DNA greedy＋beam4自然联合池，再加载原TRAIN参考标注；原A3池与标签直接重算并配对，输出逐物种完整计数与三项原生成器gate。verify_query_relative_boundary_train.py在CPU独立核对完整曝光、实际标签／池身份／gate、末模型结构和新增层实际梯度及权重变化。脚本本机AST解析和三项SBATCH语法检查通过，窄范围远程传输逐文件直接cmp一致；远程SBATCH语法通过。

fit／readout分别1RTX3090／2CPU／16GiB，上限3小时／2小时；核验0GPU／2CPU／8GiB／30分钟。预计拟合参考A17约82分钟，新增结构仍有未知余量，给3小时而非改变曝光终点。全部新阶段独立目录，不覆盖既有候选、数据或checkpoint。科学FAIL关闭A23；PASS只准备另行固定的完整caller和评分协议。等待使用同一个固定轻量监控，不把提交当结果。

首次提交带已完成预检13391430的afterok依赖，被Slurm拒绝Job dependency problem，没有创建job。实际sacct仍确认预检三记录COMPLETED0:0，scontrol show job已Invalid job id，控制器MinJobAge为300秒；与短预检已移出活动控制器的事实一致。账户队列没有A23fit，输出目录不存在。仅移除这个已完成阶段的Slurm排队依赖；训练自身仍严格核对相同预检STATUS、configuration和零更新布尔。资源／科学命令／曝光不改，fit后的readout和CPU依赖保持afterok。该事件是提交拒绝，不是训练失败或科学FAIL。

实际命令为sbatch sbatch/M28-A23-RELATIVE-BOUNDARY-FIT.sbatch，返回13391608；sbatch --dependency=afterok:13391608 sbatch/M28-A23-RELATIVE-BOUNDARY-TRAIN-READOUT.sbatch，返回13391609；sbatch --dependency=afterok:13391609 sbatch/M28-A23-RELATIVE-BOUNDARY-RESULT-CHECK.sbatch 13391608 13391609，返回13391610。每次提交后即时squeue核对，主取证2026-10-05 16:38:08 UTC：fit RUNNING gpu034，readout及CPU均PENDING有效Dependency；三个上限分别3小时／2小时／30分钟。未提交重复fit、未改科学命令、未取消旧作业。实际训练／自然覆盖与gate仍未闭合。

首次运行路径核验：fit STATUS RUNNING／JOBID13391608、stderr无错误尾部，stdout已报128／4608次更新。并发读取时training.jsonl最后已flush记录为step64，拟南芥正链case的新增输入层／输出层梯度范数分别0.011060079157248555／0.3339856755563225，总梯度范数9.010366439819336，均有限。这直接补齐“零末层之后输入层能收到实际梯度”的工程证据，不是最终拟合完成、覆盖改善或科学PASS。不会为该检查增加更新次数或重跑未变检查。

同一heartbeat coding-rna-r5-light-monitor已启用，kind heartbeat／固定target01a0f71f-8845-7432-bd13-bf54c1808fb3／15分钟／failed_runs_only／created_at与完整新提示精确读回通过。新cohort为M28-A17-LABEL-R2-A23-TRAIN-20261006，A17原fit13387910＋标签13391371＋CPU13391372，以及A23上述三根；不监控已消费的旧失败依赖。相同对话manual startup实际2026-10-05 16:39:21.553 UTC取证，A17标签与A23fit RUNNING，三个下游为有效Dependency；20.305秒完成且error null，实际gpt-5.6-luna/low。没有正常态主对话消息。当前仅配置和单次检查验证成功；本新cohort定时执行与终态回传尚未发生，不保证机器休眠或断线时即时唤醒。任一完整链终态或新可行动故障时仅交接一次再暂停同一心跳，主对话继续已授权结果分析。
