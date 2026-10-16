# A32：search-retention generator 与固定 local caller 组合

冻结日期：2026-10-06。本协议在主对话读取 A30 最终自然池及其 CPU 门槛结果之前写定。初始仅准备短合成实现检查，实施回执见下文；不预判 A30 科研 PASS。

## 来源与停止规则

A30 reconciliation13399750、readout13399752、CPU13399755 必须实际 COMPLETED/0:0，并完成原标签与五门槛核验。原五门槛为：两个物种各自自然池正确链严格增加、合计正确 ≥5-exon 链严格增加、两个物种各自 known-nonmatch 不增加。任一 FAIL，A32 不执行；不得换 seed、checkpoint、系数、搜索范围或阈值补救。

原 fit13398266 的 Slurm root TIMEOUT 保留，不改称成功。已完成的原始4608步与最后004608 checkpoint 经13399750核验，A32 只读该模型，不恢复训练或增加曝光。A31 readout13398521/CPU13398524 是已闭合的主基线来源。

## 固定组合与人口

主比较 A32_A30 对 A31_A3，而非较弱的 greedy 或 DEV 选出的胜者。两者同649457参数、seed0、3×1536训练步、最后004608、同 GEN 缓存 revision 与固定15类 segmentation/local score。A30 唯一已实施的训练目标差异是原批准的 current-search first-loss competition；A32 本身0新参数、0 optimizer updates。

复用已闭合的 A30 greedy+beam4 union，不再搜索、扩 beam、注入 reference 或调用 GEN backbone。beam4/completed4/max128 exons、原 ORF 与5bp intron规则不变。评分复用 `sum_valid_candidate_span_role_logit_minus_intergenic_logit`，role权重全1、null0、FP32 forward/FP64累加、原 span非重叠DP、不扫描阈值。

沿用原1491个唯一 TRAIN windows、同六条 TRAIN contigs、同顺序与有效长度。原始完整自然 score 文件关闭后，才读上游保存的逐窗口 R4 labels。不得打开 DEV、sealed/test 或 Setaria；不重新生成 R4 标签。上游 A30 与 A31 的真实 CPU 核验分别提供标签有效性依据。

池/标签/基线按完整 geometry 与顺序对应；基线逐窗口 counts 和汇总必须精确复现。保留 Arabidopsis2172、rice744，共2916 eligible window/reference pairs；正确 ≥5-exon 分母854/232。unknown 仍参与 DP 竞争，独立报告，不算已知 FP。这里 P=TP/(TP+knownFP)、R=TP/eligible、F1=2TP/(eligible+TP+knownFP)，均是 TRAIN window/chain proxy，不是全基因组 gene F1。

## A32 预先固定的五项 caller 门槛

两个物种各自 proxy F1 严格增加、各自 proxy precision 不下降、合计 selected correct ≥5-exon 严格增加。用整数交叉乘法判断，均与 A31_A3 比较：Arabidopsis TP1531/knownFP347/正确长链548；rice TP478/knownFP251/正确长链152。

五项全 PASS 仅支持另外明确固定的后续 raw-FASTA annotation 研究；不自动启动 DEV、放行封存数据、改变 claim 或对外发布。任一 FAIL 关闭本组合，不改 role weight/null/选择器救结果。A30 与 A31 原门槛不变。

这是同容量、同曝光、单seed的目标/模型/自然池组合对比；公共表征权重和池均可变化，不能宣称纯 pool 效应、独立生物学因果、模型新颖性或泛化成立。

## 实现检查与资源

四项纯合成 CPU 检查覆盖：原 A30 任一门槛 FAIL 即停；停止先于实际权重/cache/population访问；unknown竞争及池/标签/基线错配；新五门槛的严格性与 precision约束。无真实模型、标签或 optimizer。复用已核验的 local scorer/selector，不修改运行中的 A30/A31 源码。

真实读出保存自然 scores、选择、标签及 counts；CPU 重新计算逐窗口配对、选择、分母、proxy和整数五门槛。CPU使用已经独立核验的上游R4标签，不独立重跑 emission forward、catalog或新R4。

短合成检查：private0GPU/2CPU/8G/5min。未来真实读出：private1RTX3090/2CPU/16G/1h；原 A31 两模型加 greedy replay 用时1:00:41，此处单模型且无新搜索/greedy replay。未来 CPU 保存文件检查：0GPU/2CPU/8G/5min。新独立输出预计小于2GiB，真实提交前核实存储与实时队列。

现在只允许提交短合成检查，由主对话消费、不加入现有七root监控。真实 GPU 不预提交；须先消费 A30 完整链与五门槛。未来需要等待时，仅更新同一固定监控对话及心跳，不新增聊天或重复轮询。工程错误可以依据具体证据最小修正，失败回执保留；科学 FAIL 不作工程重试。

## 实施回执：2026-10-06 06:50 UTC

实现、四项合成测试、三份 SBATCH 和本协议已窄范围上传到 canonical remote。两个 Python 文件本地/远程 AST 解析、三份 shell 脚本语法检查通过；未修改运行中的 A30/A31 代码或已有实验输出。真实 GPU 脚本使用实际 `gpu:nvidia_geforce_rtx_3090:1` 类型。

提交前五项实时调度检查已读取；private gpu034/035 有资源，账号队列远低于80作业限制，已声明的预约无冲突。短 CPU13400182于06:49:06 UTC提交、06:49:07启动、06:49:32结束，实际 root/batch/extern 均 COMPLETED/0:0，用时25秒、MaxRSS718764K；STATUS=COMPLETED、JOBID匹配。四项 unittest 全 PASS/0.004秒，实际0GPU/2CPU/8G/5min，无真实权重/cache/标签访问或 optimizer update。

A32 真实 readout 与结果 CPU 尚未提交，三处 A32 输出在短检查提交前均不存在；只有短检查生成了其自己的输出。A30 recovered readout 和 A31 两条 DEV inference 当时仍运行，下游三个 CPU 均等待有效依赖。现有七root固定轻量监控与15分钟心跳不变，短13400182已由主对话消费、不加入该监控。实现准备通过不等于 A30 gate、A32 proxy、gene F1 或论文结果。

## 来源闭合后的不执行判决

2026-10-06 07:44 UTC。A30三根13399750／13399752／13399755均实际COMPLETED/0:0，但[原候选五门槛](a30_search_retention_train_result.md)全部FAIL。按本协议，A32真实GPU读出与CPU检查不执行；不是A32科学FAIL或工程失败，不产生模型、评分或DEV结果。合成检查和条件实现保留，不换源checkpoint、seed、评分、role权重、null、选择器或门槛补救。独立A31原生DEV与其诊断不受此停止判决影响。
