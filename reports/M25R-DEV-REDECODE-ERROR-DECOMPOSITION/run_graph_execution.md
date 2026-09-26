# M25R-E1-RUN-GRAPH-R1 execution

2026-09-08：用户在当前对话明确回复“批准”，授权固定 R1 合同的实现、一次 CPU 合成验证、一次全量 GPU 前向及 A 重放、A 通过后的 CPU B 解码/评价。原合同中的“待批准”是提案时状态，不要求重复批准。

- 输出：`outputs/M25R-E1-RUN-GRAPH-R1`；保留历史所有结果。
- 算法 CPU 作业：`12522142`，`sbatch sbatch/M25R-E1-RUN-GRAPH-TEST.sbatch`，private-teodoro-gpu，无 GPU，2 CPU / 16 GiB / 1 h。COMPLETED 0:0，19 s，MaxRSS 485660 KiB。8 项测试通过，含 250 个随机小图的正负链穷举对照。
- 唯一 GPU 作业：`12522183`，`sbatch sbatch/M25R-E1-RUN-GRAPH-FORWARD.sbatch`，private-teodoro-gpu，1 RTX3090 24 GiB / 8 CPU / 32 GiB / 16 h。作业内先做两个导出测试，然后前向 17384 个窗口、原 A 重放及独立 A 评价。提交前 Slurm 确认两私有节点各剩一个 GPU，有足够 CPU/内存；无维护 reservation，BeeGFS 空余约 212 TiB。
- CPU 解码：脚本 `sbatch/M25R-E1-RUN-GRAPH-DECODE.sbatch` 已准备，尚未提交。只在 GPU 正常完成且 `A_replay_check.json` 的 `passed=true` 后提交；8 CPU / 32 GiB / 4 h，无 GPU。

实现文件：`run_graph_core.py`（候选/21 种密码子后缀/图 DP）、`run_graph_run.py`（隔离推理、A 重放、B 解码、评价进程），均在原诊断脚本目录。核心使用坐标扫描、donor 激活堆及按 run ID 的前缀最优树，避免所有 donor×acceptor 显式枚举；同时处理 ±6 bp 导致的 motif 坐标与 run ID 顺序不同。

实现细节不改变候选或评分：

- 内部终止匹配时把非 ACGT 字符映为 C（C 在三种 stop 的任意位置均不出现），保留原“不匹配 stop”的行为；原 FASTA/GFF/ATG 和 stop motif 不作替换。
- 21 种普通后缀状态外，未完成首个 ATG 的短首 CDS 有两个初始化状态 A/AT；这是保证原完整 ORF 条件的起始状态，不扩充生物学接受规则。
- 起始/终止密码子跨短 CDS 时，B GFF 输出实际 CDS 上的分段 codon feature。独立审计按拼接 CDS 首末三碱基坐标核对。原 A 输出及坐标评价不改动；不把跨内含子的连续 3 bp 伪作真正 codon feature。
- 完整 GFF 落盘并退出预测进程后，评价进程才读取开发参考。B 候选、评分、选路不读取原真链、oracle 诊断或参考注释。

任何作业/测试失败、超时、A 不一致、输出不完整或独立结构审计失败立即停止；本轮没有修复重投、追加前向、调权/阈值、训练、Setaria 或 ablation 授权。B-A 不作为 ablation gain。只有单模型必要门槛通过也不能宣布整个方法 GO。

跨时段临时静默监控：已询问用户，尚无新授权时不创建。结果处理属于已批准范围；作业提交不是结果交付。
