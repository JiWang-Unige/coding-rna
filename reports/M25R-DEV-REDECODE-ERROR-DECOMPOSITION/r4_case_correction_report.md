# R4 CPU 大小写勘误：完成报告

日期：2026-09-08。用户已明确批准本次 CPU-only 勘误及完整结果交付。

**交付结论：勘误完成，工程与本次科学一致性核验通过；M25R development NO-GO 不变。** 原“约36% emitted structures 不合法”和“大量 reference 不兼容 canonical grammar”的解释撤回。未进行新训练、baseline、SegmentNT、GPU re-decode、第三方评阅或 Setaria 访问；没有 Git commit/push。

## 作业与实现

- 唯一提交：`sbatch --parsable sbatch/M25R-R4-CASE-CORRECTION.sbatch`。
- Job `12496928`，`private-teodoro-gpu`，gpu034，2 CPU、16 GB RAM、30 分钟上限；Slurm AllocTRES 仅 CPU/memory/node，无 GPU。
- Slurm `COMPLETED`，ExitCode `0:0`，elapsed `00:03:06`；MaxRSS `2548800K`。开始时间 2026-09-08 11:16:18 Europe/Zurich。没有重投、重跑或后台监控。
- 提交后将该作业 `Requeue=0`，禁止 scheduler 自动重排；交付 sbatch 也显式写入 `--no-requeue`，未因此重投。
- 环境复用原 R4 的 `generanno`，通过 `/opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh` 激活；未安装依赖。所有测试与真实数据处理均在该 CPU Slurm 作业内完成，登录节点只做文本、文件与调度操作。

实现范围：

1. `scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/redecode_error_decomposition.py`：仅在 `load_species()` 中将读取结果转大写，不修改公共 reader 或冻结 decoder。
2. 同目录新增 `posthoc_r4_case_correction.py`：只读原 GFF3、reference/lineage/diagnostic 和已批准 A/rice development 数据；复用原 parser、validity、coordinate metrics、canonical 和 strata 函数。未加载 checkpoint 或调用模型 forward。
3. `tests/test_m25r_case_correction.py`：5 项聚焦回归；与原 15 项测试一起执行，**20 passed in 24.86s**。
4. `sbatch/M25R-R4-CASE-CORRECTION.sbatch`：CPU-only 一次性执行脚本。

数据仍为 Arabidopsis NC_003074.8 和 rice NC_089041.1 的原 6,450 validation references，primary transcript/split/config 不变。公共训练读取器、原始数据、checkpoints、threshold、radius、motif 字典、decoder grammar 均未改动。工作树中的既有旧框架退役改动被保留，未恢复。

## 独立输出与冻结原始结果

原结果保留于：

`outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4/`

纠正结果单独写入：

`outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4-CASE-CORRECTION/`

输出包含：

- `posthoc_case_correction_summary.json`；
- 每 epoch 的 `correction_summary.json`、`structural_validity_corrected.jsonl`、`reference_case_correction.tsv`；
- `STATUS`、`JOBID` 以及 Slurm stdout/stderr。

最终 STATUS：`POSTHOC_CORRECTION_COMPLETED_REVIEW_REQUIRED`。这里的 REVIEW_REQUIRED 表示停止自动研究推进、交付给用户，不要求额外第三方评阅。

程序运行前后直接比较原 R4 全部文件的 size/mtime，结果一致；所有写入路径均为独立 correction 目录。原 R4 STATUS、scientific_status 和旧 ledger 保留作历史记录，由本报告及新 summary 明确 supersede，不静默覆盖。

## 纠正后的结果

| 指标 | epoch 1 | epoch 2 | epoch 3 |
|---|---:|---:|---:|
| 原错误 validity 有效数 | 1338/2098 | 745/1180 | 429/667 |
| 纠正后独立有效数 | 2098/2098 | 1180/1180 | 667/667 |
| validity audit coverage | 100% | 100% | 100% |
| 所有 validity 组件失败数 | 0 | 0 | 0 |
| transition reachability（不变） | 3727 | 3993 | 4017 |
| 原错误 motif reachability | 1693 | 1758 | 1738 |
| 纠正后 motif reachability | 3616 | 3852 | 3862 |
| motif / 全部 6450 references | 56.06% | 59.72% | 59.88% |
| canonical subset motif count | 3614 | 3849 | 3859 |
| canonical subset motif fraction | 59.59% | 63.46% | 63.63% |
| emitted exact chains（不变） | 1389 | 907 | 499 |

canonical-compatible reference 从 **2666/6450** 纠正为 **6065/6450（94.03%）**，三个 epoch 集合一致；不兼容数从错误的 3784 降为 385。3,399 条 reference 的 canonical 标签发生纠正。这里 canonical-compatible 包含 ATG/stop、GT–AG、完整 ORF/no-internal-stop 等冻结 grammar 条件，不只是 intron motif。

按物种，canonical-compatible 为 Arabidopsis **3948/4151**、rice **2117/2299**。主 stage assignment 不变，仅 canonical/noncanonical 分层重新计数。全部 exact emitted chains 均属于纠正后的 canonical 集合，且逐 reference 满足 transition 与 motif reachability。

有 2/3/3 条 noncanonical reference 仍满足 motif reachability 并不矛盾：motif reachability 不检查完整 ORF/frame/internal-stop；不能把 motif 可达等同于可输出合法完整链。

## 完整性与复现核验

- 20 项聚焦测试通过；新增覆盖诊断入口 mixed-case、正负链 validity、trace/production 大小写不变性、状态重建、索引算法与原 reachability 函数等价及重叠 block 拒绝。Python syntax、sbatch `bash -n` 和 tracked diff whitespace 检查通过。
- 从原 GFF3 重新解析并计算四项 coordinate metrics，逐值等于 R4 原 `reproduction.*.replayed`：

| 指标 | epoch 1 | epoch 2 | epoch 3 |
|---|---:|---:|---:|
| exact CDS interval F1 | 0.120414079 | 0.078302211 | 0.038261788 |
| exact CDS chain F1 | 0.324988301 | 0.237745740 | 0.140227624 |
| pooled intergenic FPR | 0.012468286 | 0.005077354 | 0.003415162 |
| gene-count ratio | 0.325271318 | 0.182945736 | 0.103410853 |

- 每 epoch correction reference ledger 均有 6,450 个唯一 key，与原 ledger 完全相同。交付前另行进行纯文本字段比较：除 canonical、motif 和已废弃 truth-assisted 三列外，所有单元格逐一不变，包含 main stage、transition、CDS-count 和 edge/span 分层字段。
- corrected validity ledger 行数分别 2098/1180/667，与原 GFF3 和 summary 一致，组件失败为零。
- 从已保存 block_span/runs 重建 I/C/G 状态，使用原 transition grammar、event anchors、半径和一对一匹配规则；每条 reference 的 transition flag 均精确复现旧值，而不仅是总数相等。
- exact-chain identity 与源 stage ledger 逐条对账；所有 exact chains 均被 corrected motif 覆盖，各 epoch 均满足 `exact <= motif <= transition <= 6450`。
- canonical/noncanonical strata 已重算，其他 strata 与原 summary 完全一致。

## Truth-assisted：明确不可恢复

最终 summary、各 epoch summary 以及全部 correction reference rows 均明确标记：

`NOT_RECOVERABLE_WITH_SAVED_ARTIFACTS`

原 897/577/320 不再沿用，也未用 exact-chain 或 motif count 代替。原 lineage 只保存 candidate/chosen coordinates、部分 chosen scores 和 phase checks，无法恢复每个 truth coordinate 所需的 boundary/phase logits。未为补齐该项重跑推理。

## 科学解释与下一步建议

1. **撤回错误归因。** 不能再说约36%的 emitted structures 不合法，也不能说大部分 references 不兼容当前 canonical grammar。大小写修正后 94.03% references 兼容，所有 emitted models 通过独立结构规则检查。
2. **M25R NO-GO 不变。** 原 gene-count ratio 仍只有 0.3253/0.1829/0.1034，远低于原 0.80–1.20 合同；结构合法性不是准确率或足够召回。原 0/5625 admissible 起点未因本次勘误变成模型成功。
3. **motif 本身不支持此前的强瓶颈解释。** 在纠正后的 canonical subset，transition 与 motif counts 分别完全相同（3614/3849/3859）；这里可达 transition 周围都能找到对应 truth motif。对全体 references，reachability 仍受冻结 state path/grammar/radius 条件限制，不能称 backbone 上界或 F1 上界。
4. **可达性改善而完整链恢复下降。** motif count 随 epoch 增加，但 exact count 从1389降到499；这与保留的 phase-check attrition 增长和完整 grouping 损失共同表明问题不只是有没有 motif。该差额不是可兑现的 oracle gain，因为完整 truth-assisted 检查不可恢复；不能据此直接授权放松阈值或更换 decoder。
5. **下一步建议：** 若继续研究，优先用现有 corrected ledger 与 candidate traces 区分 complete-block grouping 和 phase-check 损失，再决定是否值得设计单一新对照；不得把本次勘误推成 backbone、LoRA 或监督方式的因果结论。此建议未执行，未自动启动任何研究分支。

本次授权的实施、一次 CPU 作业、结果复核、canonical/motif 勘误和报告交付均完成。无临时监控需要保留。
