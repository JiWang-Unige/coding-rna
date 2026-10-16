# M28 A9 冻结首拟合与 TRAIN 配对诊断

2026-10-03，在新增头任何 optimizer 更新及 trained-critic 解码前固定。用户已批准继续新模型迭代，不再逐步重复询问。目标及数据范围不扩展，旧 B2/A8 负结果不挽救或覆盖；FULL/A6/A7 保持原冻结运行。

零更新准备13359370已COMPLETED/0:0（55秒、7 tests通过）。8 TRAIN windows、27 primary reference/window pairs、216 exon，432 positive/35,557 negative/33,700 unknown prefix events；299 positive仅在TRAIN注入/teacher-history池，不是free覆盖。START为27/29,712/32,357，NEXT_START为189/731/553，DONOR为189/4,663/533，STOP为27/451/257（正/负/unknown）。所有正事件齐全、grammar-unavailable为0、8窗梯度有限非零、生成器无梯度、头权重未变。标准库独立回读身份/positive前缀/标签计数并从logit重算loss，最大差2.904e-8。START极稀疏是已知风险，保留known mean目标，不在看到训练或DEV后改类权重。

## 唯一拟合

冻结A3 `outputs/M28-A3-FIRST-FIT/step_004608.pt` 的全部生成器和已有GENERanno输入；只训新query-conditioned PrefixValueHead(H128)67,745参数，seed0、FP32。与旧B2共享“学习前缀评分”的概念，差别与限制见[a9_query_prefix_preparation.md](a9_query_prefix_preparation.md)，不是独立创新已成立。

原1,536 draws/1,491 unique TRAIN windows原顺序重复3遍，恰4,608更新，importance weight不变。AdamW lr3e-4、wd.01、betas(.9,.999)、eps1e-8、norm clip1；唯一正式头为最后step_004608，不按TRAIN probe或DEV选择epoch。继承生成器4,608更新再新增头4,608更新，不宣称与原C0/A3总曝光或参数相等。

先用冻结生成器完成1491窗原greedy/beam4 label-free trace，关闭记录后再加载TRAIN catalog。每draw重新计算冻结表示、全部初始ATG事件、原自由路径前缀、冻结A3 matching及teacher-history counterfactuals；不缓存新增头值/图，不在搜索中注入参考。真实prefix的TRAIN注入单独处理，unsupported语法正事件仍保存并计数，不删目标。每事件类型known BCE均值，再对有known样本类型等权平均，乘原importance；-1保护、空池零梯度、无虚构positive与准备口径相同。

拟合新目录 `outputs/M28-A9-QUERY-PREFIX-FIT`。资源1 RTX3090、2CPU、16GiB、3小时；既有A6固定自由池生成约2006秒，加本轮8窗反向平均约1.04秒/窗折算4608 draws约1.89小时，留运行余量，不保证完成。输出自由trace、4608训练行、3个epoch checkpoint/summary，预计<300MiB。超时/OOM/非有限/基座梯度/身份错配停止保存，不续半成品、自动延时或改变科学设置。

## 拟合后 TRAIN 机制筛选

afterok拟合，固定最后头与同一A3生成器；重复既有8 metadata-selected TRAIN pressure windows，不换窗口或读DEV/test/Setaria。early用新critic全部合法entry/end排名，beam4/complete4/cap128，无错误raw-score提前终止；late使用既有同checkpoint原greedy+beam4真实自由池。late多一个greedy来源且early搜索工作量不同，不称唯一变量/预算或运行时间完全相等。

两臂(query index, exact exon chain) union做一次新的完整STOP readout，再映回两臂；共同query-chain gain精确相同。query不同的同geometry是不同slot compatibility，未来窗口读出须显式max，不把它们当多个基因或自然正确概率。本诊断只检验自由完整链可达性，不产生最终annotation/F1/background。

新目录 `outputs/M28-A9-QUERY-PREFIX-TRAIN-PROBE`，1 RTX3090/2CPU/16GiB/30分钟，预计<1MiB。预测全部关闭后才接入原27 reference/window identities，报告early/late、两物种、new/lost与>=5-exon new，不只列改善案例。

进入后续完整DEV准备的最低机制条件：early总free exact pairs严格多于late，且至少新增1个>=5-exon reference/window pair，共同值一致、原分母27完整、全部返回链DNA合法。它是探索性TRAIN筛选，不是准确率/泛化/论文gate；单物种损失和新增candidate负担仍完整报告。失败保留，不扫beam、阈值、seed或追加训练挽救本头；若通过，另在任何DEV推理前固定完整paired raw-FASTA/统一owner/selector/F1/background协议，而非直接宣布成功。

该两作业在现有授权内各提交一次；独立fixed watcher复用同一对话和心跳，正常等待静默，拟合失败或两阶段完成一次唤醒原主对话。本轮不公开发布或Git提交。

## 实际提交与启动

源码、两sbatch和协议部署后逐字符回读一致，AST/bash语法通过；实时Slurm gpu035有2张未分配RTX3090、CPU22/126、RAM208/503GiB，无本轮GPU维护限制，因此选private-teodoro-gpu。以下命令各提交一次：

```bash
sbatch --parsable sbatch/M28-A9-QUERY-PREFIX-FIT.sbatch
sbatch --parsable --dependency=afterok:13359395 sbatch/M28-A9-QUERY-PREFIX-TRAIN-PROBE.sbatch
```

拟合13359395于gpu035 RUNNING，probe13359396 PENDING有效Dependency；拟合作业开始前新增union readout后共8项合成测试通过（1.322秒）。2026-10-03 12:39:31UTC已有32/1491固定自由trace，尚未完成新头4608更新，不把启动、测试或free trace进度称作模型效果。

同一heartbeat `coding-rna-r5-light-monitor` 的完整prompt已加入这两root IDs及最小完成字段，固定target仍01a0f71f-8845-7432-bd13-bf54c1808fb3、模型沿用gpt-5.6-luna/low、15分钟；ID/kind/目标/间隔/状态/全prompt/创建时间已实际回读一致。原FULL/A6/A7仍监控，不创建第二个monitor或cron。该新cohort的定时检查与最终回传尚待实际观察；旧cohort成功回传不替代本轮证据。协议和科研结果不由watcher改变。

## 等待拟合期间补齐诊断记录

2026-10-03，在probe13359396仍为PENDING/Dependency时补齐原协议已要求的candidate负担与new/lost报告。`fit()`、CONFIG、LIMITS、new_head、checkpoint identity的AST与提交快照完全一致；正在运行的fit不改模型、监督、更新、采样或资源。probe的原free继续条件表达式不变，只增加已保存free池的读出记录，不做新搜索、forward、loss或阈值选择。

每窗每臂分别列completed记录数、(query,geometry)去重数、geometry去重数，避免把同一链多query和greedy/beam重复当多个gene；对同geometry显式max over query的complete value，另列>0的候选数及R4 reference-relative 1/0/-1标签。0分不算正值，protected/unknown/noncallable保留-1。这个窗口池诊断不是owner/selector后的annotation，不称生物学precision或F1；common query-chain值相同不保证两臂同geometry最大值也相同，因为访问query集合可变。

每条原27 reference/window identity保存both/early_only/late_only/neither及zero-value视图命中，summary按物种与exon class列全部状态。该负担统计和zero-value视图不进入原继续gate，不用它改变阈值、checkpoint或旧结果。新CPU-only检查入口`M28-A9-QUERY-PREFIX-REPORT-CHECKS.sbatch`，private/0GPU/2CPU/8GiB/5分钟，独立输出目录；测试计数去重、zero值、真实unknown政策、空池与共同值不一致的可见报告，未触及真实DEV、TRAIN数据读取或optimizer。

实际提交一次`sbatch --parsable sbatch/M28-A9-QUERY-PREFIX-REPORT-CHECKS.sbatch`，job13359420及batch/extern均COMPLETED/0:0，15秒、2CPU/8GiB/0GPU（0.008333独立CPUh），batch MaxRSS732804KiB，3 tests PASS/0.003秒。四个修改/新增文件部署后逐字符回读一致；拟合核心及原probe gate表达式的AST均与首提交快照相同。该完成只验证新增汇总口径，没有新模型结果；原fit13359395和probe13359396不重提、不改变依赖/资源。

新cohort已有真实定时检查：`read_thread`返回本轮最外层`<heartbeat>`，automation_id为原ID、current_time_iso=2026-10-03T12:42:38.398Z；同固定对话turn01a101c9-b3c0-7d22-a1a3-eb1a0f9eb02d在12:42:50.666Z实际读取11根Slurm状态并于12:43:03Z正常完成。watcher曾错误按memory里的待启动NO将它记作人工检查；主对话已据真实触发消息修正为scheduled，并在同一heartbeat明确触发类型取证，非据历史NO/YES猜测。运行/依赖状态健康，未交接；本轮完整终态回传仍未发生。心跳自身普通通知设为failed_runs_only， meaningful终态仍一次发到原主对话；仅此心跳设置变化，ID/target/15分钟/模型/主对话设置未改，实际保存全prompt及字段回读一致。

## 最终状态：本头的冻结 TRAIN gate FAIL

2026-10-03 14:33 UTC，fit13359395与probe13359396均COMPLETED/0:0。原27条完整分母和共同query-chain值已从实际记录核对；early自由整链4条，late15条，新增0、丢失11、新增≥5-exon链0，因此原继续gate失败。本头不进入DEV，不增加训练、改变beam／阈值／seed；记录全部物种损失、候选负担和零操作点视图，见[a9_query_prefix_result.md](a9_query_prefix_result.md)。FULL、A6/A7、A10继续原运行与同一固定监控；前文等待／尚无结果为当时历史状态，不是当前终态。
