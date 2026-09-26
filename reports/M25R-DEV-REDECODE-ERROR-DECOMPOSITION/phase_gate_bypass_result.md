# M25R epoch-1 phase-gate bypass R2：冻结面板结果

日期：2026-09-08。状态：COMPLETED_FROZEN_PANEL_DIAGNOSTIC。

## 结论

在冻结的 128 条候选 lineage 上，仅绕过 phase 一致性拒绝、保持候选边界、权重、ORF 检查与边界阈值不变后，64/64 reference-exact phase-failure 候选恢复输出；32 条 reference-nonexact phase-failure 中 9 条也输出，另 17 条被 ORF 检查拒绝、6 条被边界阈值拒绝。32/32 原输出对照保持不变。

因此，phase gate 是此面板 64 条 reference-exact 候选的直接阻断点；绕过它也增加 reference-nonexact 输出，不能据此宣布整体精度、F1 或泛化提升。不能把此结果归因为 phase-head 架构缺陷，也不能据此正式关闭整个研究路线。当前开发集 NO-GO 和 Setaria embargo 不变。

## 执行与可信性

- Slurm 12520168：COMPLETED，ExitCode 0:0；2026-09-08 17:10:36–17:18:05 Europe/Zurich，耗时 7分29秒。
- private-teodoro-gpu；1 RTX3090、8 CPU、32 GiB，2小时硬上限；MaxRSS 3272676 KiB，PyTorch 峰值 allocated GPU memory 2446588416 bytes。
- 作业内前置测试 32 passed。158 个唯一原始窗口前向；无训练、无阈值搜索、无 Setaria 文件读取。
- A 原行为回放 128/128 全部核对通过后才计算 B；32 个 emitted controls 坐标、phase、codon、边界分数等保持不变。
- A/B GFF 分别含 32/105 个 transcript，结构审计分别 32/32、105/105 合格。结构合法不等于生物学正确。
- 本实验是固定候选局部门控对照，不是完整流水线重放。使用 R2 分层系统抽样面板；历史 64-window preflight STOP 保留，实际按用户后续批准的 158-window 上限执行。

## 配对结果

| 冻结组 | n | A 输出 | B 输出 | B ORF拒绝 | B阈值拒绝 |
|---|---:|---:|---:|---:|---:|
| reference-exact phase-failure | 64 | 0 | 64 | 0 | 0 |
| reference-nonexact phase-failure | 32 | 0 | 9 | 17 | 6 |
| emitted control | 32 | 32 | 32 | 0 | 0 |

Exact 恢复：拟南芥 54/54、水稻 10/10；正链 33/33、负链 31/31。首次失败位置 token-conflict 25/25、无 conflict 39/39，二者都恢复，不能单独归因为 token alias。64 条 exact 样本全部为 multi-CDS，不能外推到 single-CDS exact 候选。

Nonexact：拟南芥 22 条中 7 输出/13 ORF拒绝/2阈值拒绝；水稻 10 条中 2/4/4。Multi-CDS 28 条中 9/15/4；single-CDS 4 条中 0/2/2。17 条 ORF 拒绝中，frame-length 失败 11 条、internal-stop 失败 16 条，二者重叠 10 条；这些分项不可相加为互斥终态。

9 条新增 nonexact lineage：拟南芥正链 342、2975、4633、5800，负链 232、772、6740；水稻正链 2217、负链 3380。Reference-nonexact 不自动等同于生物学假阳性。

预先记录的探索性标志：exact恢复至少32条为 true；nonexact输出最多8条为 false（实际9条）。它们不是路线 kill criteria，也不是事后可修改的确认性门槛。分层富集面板不能用于总体 precision、F1、置信区间或全基因组收益估计。

## 下一步建议（未执行，需新的范围授权）

优先设计一次保持 epoch-1 权重和其他规则不变的完整开发集 phase-only 配对对照，检验全局 exact 指标、计数和新增 nonexact 代价；在确认净收益之前，不进入新 head 训练或 Setaria。先用本次保存的 9 条新增 nonexact lineage、分数及参考注释做针对性 CPU 诊断，可帮助解释额外输出，但不能替代完整开发集对照，也不能自行判定参考之外的生物学真伪。完整对照需要另行冻结范围及预算，不属于本次 158-window 批准。

## 证据位置（Baobab 项目根目录下）

- outputs/M25R-E1-PHASE-GATE-PANEL-R2/frozen_panel.json、window_coverage.json
- outputs/M25R-E1-PHASE-GATE-BYPASS-R2/resolved_run.json、window_runs.jsonl、scores/
- outputs/M25R-E1-PHASE-GATE-BYPASS-R2/A_replay.jsonl、replay_check.json、paired_results.jsonl、summary.json、STATUS
- outputs/M25R-E1-PHASE-GATE-BYPASS-R2/A_panel_predictions.gff3、B_panel_predictions.gff3、A_structural_audit.json、B_structural_audit.json
- logs/M25RPHASEBYP2_12520168.out、logs/M25RPHASEBYP2_12520168.err
- scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/phase_gate_bypass.py
- sbatch/M25R-E1-PHASE-BYPASS-R2.sbatch

本次结果分析完成；临时监控在交付时删除。不追加作业、不重跑、不改冻结结果、不提交 Git。
