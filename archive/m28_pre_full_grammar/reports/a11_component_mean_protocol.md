# A11：冻结 A3 上的零更新 TRAIN 组件均值探针

2026-10-03，在本新探针读取结果之前固定。用户已授权在原模型目标内继续迭代；本轮不新增拟合，也不挽救或覆盖A9 FAIL。只复用原A5的8个metadata-selected TRAIN pressure windows／27个primary reference/window pairs和原A3最后step_004608。FULL、A6/A7、A10、旧评价和test／Setaria边界不变。

## 待证伪问题与三种评分

A3训练中boundary CE按每链2n个边界均值、termination BCE按每链n个事件均值，二者系数均1；原A5 beam却累计两类log分数。A5正确前缀诊断说明某些复杂链已有四条更高累计分竞争者，但不能据此断言长度惩罚是唯一原因。本探针分离EOS相对比例与完成池跨长度竞争，不扫归一指数、宽度或阈值，也不称均值评分本身创新。

设B为该路径几何pointer log-softmax之和，E为donor continue／STOP terminate log-sigmoid之和，n为完成链CDS segment数。

| 新读出 | live／entry排名 | 完成池排名 |
| --- | --- | --- |
| RAW_SUM，无bound | 原累计B+E（沿用原逐事件浮点加法） | B+E |
| COMPONENT_SUM，无bound | B+2E | B+2E |
| COMPONENT_MEAN，无bound | 与COMPONENT_SUM同一次实际搜索 | (B+2E)/(2n)=B/(2n)+E/n |

COMPONENT_SUM和COMPONENT_MEAN共享一份真实搜索轨迹、全部实际合法STOP候选，只维护各自完成池；完成池不反馈live。每层live链都具有相同CDS段数，因此不因“部分链长度归一”另改active搜索。两个weighted读出之差只来自完成池保留顺序；RAW_SUM对weighted还包含EOS相对比例变化，不能混为单一长度作用。全部关闭提前停止，因为均值分没有原累加分的单调上界。RAW_SUM与原保存A5 bounded beam的返回身份及共同分数另作重放核对，不假定移除bound必然在浮点tie下维持身份。

这里的“训练组件比例”不等于完整训练目标复现：仍保留原DNA beam的最短5bp内含子mask，teacher几何normalizer用上一端点+1；query existence、local损失及Hungarian项不进入此完成排名。推理不读取参考，不把pointer/EOS分数叫正确概率或最终gain，不能直接接零阈值宣布annotation。

## 固定输入、读出与继续条件

仅加载 `outputs/M28-A3-FIRST-FIT/step_004608.pt`，全部生成器冻结，optimizer更新0、新头参数0、不写checkpoint、不执行新GLM forward。只读取既有8窗TRAIN特征／定向DNA；A3 seed／权重／原几何normalizer／canonical ATG、GT/GC–AG、实际codon carry、Q16、beam4、每query完成4、cap128不变。原greedy自由池原样复用，三臂都与同一greedy池取union；原greedy+bounded beam作为第四个ORIGINAL参照。

在全部新自由预测关闭后才读取原TRAIN catalog／27条参考；unknown／protected按原R4 candidate标签政策保留。保存所有27身份、每物种及1／2–4／≥5段计数、mean对component sum的新／丢链，以及每窗候选去重和1／0／-1标签负担。候选池标签0只表示参考不一致，不是生物学假基因；这些计数不是最终GFF、precision／F1或背景负担。

唯一最低机制条件：mean的pooled free exact pairs严格高于component sum，至少新增1条≥5段reference/window pair，mean在每个物种的总命中不低于ORIGINAL，全部原27保留，RAW_SUM重放身份一致且共同分数差≤1e-3。全返回链须经过原DNA合法性检查，weighted池共享搜索不变。所有条件事前写定；结果失败就关闭本读出的DEV扩展，不扫系数／归一指数／beam、改阈值或追加训练。通过也只表示值得另行冻结完整raw-FASTA annotation读出及对照，不是独立科学PASS或论文成立。

## 检查、资源与失败处理

新增5项CPU合成检查：组件均值与长短链排序、单exon均值对teacher成分及原score、weighted共享搜索／existence不影响排名、实际carry/GC/unknown/cap语法、softmask/padding/tie/Q身份。单exonteacher数值对应不扩写成多exonmask等价。先通过本checks再提交真实探针，不重复未变化的历史测试。

`M28-A11-COMPONENT-MEAN-CHECKS.sbatch`：generanno、private、0GPU／2CPU／8GiB／5分钟，新目录 `outputs/M28-A11-COMPONENT-MEAN-CHECKS`。真实入口 `probe_query_component_mean.py`，`M28-A11-COMPONENT-MEAN-TRAIN-PROBE.sbatch`：generanno、1 RTX3090／2CPU／16GiB／20分钟，新目录 `outputs/M28-A11-COMPONENT-MEAN-TRAIN-PROBE`，保存8份预测、27参考及summary，预估<10MiB。与原A5 GPU／内存相同；关闭bound可能延长搜索，20分钟是本次硬停止条件而非已测吞吐保证。根据提交时实际队列选择最早可信完成的合格分区，接近时优先private。

非有限、源身份／几何／参考错配、DNA非法或超时／OOM留证停止，不续半成品或自动延时。仅明确源码／测试实现错误允许一次不改变上述科学配置的最小工程修复，原失败保留；没有指标失败后重跑权限。所有状态与结果据Slurm／实际文件记录，主对话作解释。长等候只复用既有固定monitor和heartbeat，不创建新对话或cron。没有Git提交、公开发布、新数据或test／Setaria访问。

## 实际执行记录（2026-10-03 14:49 UTC）

上述协议、独立decoder、探针、测试及两个sbatch共6个文件已定向部署到Baobab，逐文件文本回读一致；Python AST及bash语法通过。未修改旧decoder、A3权重或其他研究分支。

CPU检查 `13359870` 根作业／batch／extern均COMPLETED、ExitCode 0:0，耗时6秒；5项测试全部通过（测试用时0.552秒）。检查产物 `outputs/M28-A11-COMPONENT-MEAN-CHECKS/STATUS` 为COMPLETED。合成检查只验证实现，不是性能或机制条件通过。

随后仅提交一次真实探针 `13359876`，依赖 `afterok:13359870`；private、gpu034、1 RTX3090／2CPU／16GiB／20分钟。在14:48:28 UTC实际调度快照中为RUNNING。产物目录 `outputs/M28-A11-COMPONENT-MEAN-TRAIN-PROBE`；日志 `outputs/M28-PILOT-R5/logs/M28A11MEANPROBE_13359876.out/.err`。尚未读取探针结果，未运行DEV扩展。

复用原监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3` 与heartbeat `coding-rna-r5-light-monitor`（gpt-5.6-luna／low，每15分钟）。当前cohort `M28-FULL-A6-A7-A10-A11-RAW-DIAG-20261003` 包含14个根作业；A9已交接关闭，A11 CPU检查已消费，不重复监控。配置实际回读确认同一心跳类型、固定目标、模型、间隔及ACTIVE状态；此前13根作业范围的14:42 UTC定时检查已在原对话完成，新14根范围尚未观察到一次定时运行。监控只判作业与产物完整性，不解读A11指标或机制条件；A11终态／新故障交接一次后暂停该轮心跳，主对话负责完整27条参考的结果判断。机器休眠／断线时不保证即时唤醒。

## 探针终态（2026-10-03 15:04 UTC）

13359876 root/batch/extern COMPLETED／0:0／1m41s，8预测和27参考完整；14:57:38.974 UTC同一固定对话的新14根范围自然心跳已实际运行并交接A11及A6/A7终态。主对话CPU13359929核对真实预测池、原27身份、共享搜索和raw score重放通过。所有臂命中15/27，MEAN无新增／丢失，也无新增≥5段参考，冻结机制gate FAIL，关闭DEV扩展。详细结果统一记录于[a11_component_mean_result.md](a11_component_mean_result.md)，不将完成池长链增加冒充准确率收益。
