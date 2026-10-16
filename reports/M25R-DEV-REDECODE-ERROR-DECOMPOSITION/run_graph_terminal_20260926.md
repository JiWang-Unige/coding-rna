# M25R 固定直接结构注释路线：NO-GO 结案

日期：2026-09-26（Europe/Zurich）。权威计算产物位于 Baobab `/home/users/j/jwang/coding-rna`；本报告及小型评价 JSON 同时交付到本机。

## 结论

**工程恢复完成；冻结 R1 开发评价有效，但科学判定为 NO-GO。按当前 goal 允许的 NO-GO 终态，关闭本项目已测试的 M25R 固定直接结构注释路线，不再追加局部调参。**

联合图解码确实产生了很大的结构恢复收益：完整精确链从 1,389 增至 3,351，CDS interval F1 从 0.1204 升至 0.7345，chain F1 从 0.3250 升至 0.5303，基因数比例恢复到 0.9594。但 interval、chain 和 intergenic FPR 三个必要门槛均未通过；FPR=0.05606 还触发冻结的 `>0.030` 立即停止条件。结构合法、数量接近参考，不能替代准确性和假阳性控制。

关闭范围是当前 GENERanno 1.2B CDS-preview、M25R 已训练 heads/LoRA、固定 checkpoint、候选生成和已测试解码修正的组合，不是“GENERanno 架构普遍不可能成功”。不宣布成功方法、不声称跨物种泛化、不释放 Setaria，也不把这个负结果包装成已达到发表标准的正面方法。

2026-09-26 用户明确要求“这种情况不再需要我的授权，以自动推进完成goal为第一要义”，因此本次工程修复、结果处理和有边界的收尾不再另设批准循环。冻结科学门槛、数据限制与原始结果保留要求未取消；没有新训练、阈值搜索、后台监控、Git commit/push 或公开发布。

## 执行与证据完整性

原 GPU 作业 `12522183` 已完成 17,384 个原窗口前向并保存约 2.19 GiB float16 logits，但 A 跨进程评价因 JSON 内层坐标 list 不可哈希而失败。这是工程错误，不是模型失败；原失败目录和 STATUS 保持不变。

恢复作业 `13225219`：

- `private-teodoro-gpu`，gpu034，0 GPU / 8 CPU / 32 GiB，申请上限 4 小时。
- Slurm `COMPLETED`，ExitCode `0:0`，08:49:56–09:01:55，耗时 **00:11:59**，batch MaxRSS **3,876,300 KiB**。
- 7 项 JSON 往返、只读缓存和 GFF 导出测试全部通过，pytest 28.16 秒；stderr 为空。
- A 的 2,098 条完整 CDS 链及 phase 集合与 R4 完全一致，四项冻结指标绝对误差均为 0；独立结构审计 2,098/2,098 合法。
- B 使用同一份缓存；生成、评分、选路和导出进程不读取参考 GFF。完整预测保存后，独立进程才读取开发参考评价。
- B 模型与 GFF 的链/phase 集合一致，无重复；独立审计覆盖 6,188/6,188，全部合法。
- 开发范围仍为拟南芥 `NC_003074.8`、水稻 `NC_089041.1`，6,450 条 complete-primary 参考链；没有按 canonical 子集缩小主分母。

最小修复在 JSON 输入边界恢复 CDS tuple，并允许独立输出目录读取旧缓存。没有改原指标函数或图核心。此前合成图验证为作业 `12522142` 的 8 项测试；本次没有把恢复测试冒充新的全算法独立评阅。

## 冻结主结果

| 指标 | 原 A / R4 epoch 1 | R1 B | 原必要门槛 | B 判定 |
|---|---:|---:|---:|---|
| exact CDS interval F1 | 0.120414079 | 0.734515969 | ≥0.80 | 失败 |
| exact CDS chain F1 | 0.324988301 | 0.530305428 | ≥0.55 | 失败 |
| exact coding-gene F1¹ | 0.324988301 | 0.530305428 | ≥0.50 | 通过 |
| intergenic FPR | 0.012468286 | 0.056056158 | ≤0.020 | 失败；触发 >0.030 停止线 |
| predicted gene-count ratio | 0.325271318 | 0.959379845 | 0.80–1.20 | 通过 |
| 独立结构合法比例 | 1.000000 | 1.000000 | ≥0.99 | 通过 |
| matched-gene strand accuracy | — | 0.999604587 | ≥0.98 | 通过 |
| exact-matched CDS phase accuracy | — | 0.996203857 | ≥0.90 | 通过 |

¹ 冻结实现中 coding-gene F1 使用同一 complete-primary exact-chain 定义，不是另一项独立证据。A 的破折号仅表示此表不重列该诊断，不能当作未计算或零分。

链计数：R=6,450，P=6,188，TP=3,351，缺失参考链=3,099，reference-nonexact 预测链=2,837；`2×3351/(6450+6188)=0.5303054281`。链 precision=0.541532，recall=0.519535。这里 nonexact 指不匹配冻结参考，不断言每条都是生物学假基因。

相对 A，interval F1 增加 **0.614102**，chain F1 增加 **0.205317**，FPR 增加 **0.043588**。A→B 按完整 CDS 链加 phase 比较：保留 2,003，新增 4,185，删除 95；这是一种联合解码机制的变化，不能单独归因为取消 phase、跨 block、ORF 约束或评分。

FPR 是预测 gene span 落入冻结参考 intergenic 区域的碱基比例，不是 nonexact 链占预测链的比例。strand accuracy 的条件分母是 5,058 个 overlap 匹配基因；phase accuracy 的条件分母是 22,918 个 exact-matched CDS，均不是全参考基因的准确率。

## 每物种诊断与错误边界

| 物种 / 染色体 | 参考链 | 预测链 | 精确链 | interval F1 | chain F1 | FPR | count ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| 拟南芥 NC_003074.8 | 4,151 | 4,317 | 2,348 | 0.752671225 | 0.554558337 | 0.127439262 | 1.039990364 |
| 水稻 NC_089041.1 | 2,299 | 1,871 | 1,003 | 0.693841520 | 0.481055156 | 0.012003732 | 0.813832101 |

原门槛是 pooled 判定，未事后增加逐物种 gate。逐物种结果揭示显著异质性：拟南芥 chain 已接近/越过 pooled 数值门槛，但 intergenic FPR 为 12.74%；水稻 FPR 为 1.20%，完整链恢复仍不足。整体失败已经由原 pooled 门槛成立，不依赖新增约束。

既有 overlap-degree 诊断给出：132 个预测跨度重叠多个参考 transcript，59 个参考 transcript 重叠多个预测跨度；拟南芥为 109/27，水稻为 23/32。这是忽略 strand 的跨度重叠代理，不能直接称为 132 个已证实 fusion 和 59 个已证实 split。完整边界偏移分布保存在 `B_evaluation.json`，未用于调参或选择结果。

## 有限机制修正的证据链

1. **原 M25R 固定方案失败。** 原 5,625 个有限网格点没有 admissible tuple；R4 完整三 epoch 重放保留原坐标指标。后续大小写勘误证实旧“约36%结构非法”是诊断假阴性，应撤回；三 epoch 实际输出全部合法，canonical-compatible reference 为 6,065/6,450。该勘误没有改善原模型召回。
2. **仅取消 phase gate 不足以达标。** 全量回顾性计数最多新增 1,149 条精确链；即使零新增 nonexact，chain F1 上界仍为 `2×(1389+1149)/(6450+2098+1149)=0.5234608642<0.55`。这是固定 phase-only 方案的严格乐观上界，不是富集小面板外推。
3. **跨 block 原 run 图有空间，但 oracle 不等于成绩。** 放宽边界分数与 phase 的可达性诊断有 3,585 条 canonical 完整参考路径，理想零 nonexact 的 chain F1 上界 0.7144993。它只说明值得测试一次固定无参考选路，不能声称实际模型达到 0.7145。
4. **固定无参考 R1 实测仍失败。** B 恢复了 3,351 条精确链且 gene count 达标，但仍有 2,837 条 reference-nonexact 预测，interval/chain 精确度不足，FPR 超过硬停止线。工程与结构合法性已验证，不能再用尚未运行或接口错误解释这一次 NO-GO。

可以支持的判断：原整 block/phase 解码是一个实质性的恢复瓶颈，替换它能释放结构信号；但当前已训练 emissions、候选和未校准的固定评分/选择组合，不能同时达到完整结构与 intergenic 特异性要求。

不能支持的判断：全部剩余错误都由评分函数导致、换权重必定成功、backbone 没有结构信息、所有 GENERanno 方案不可能成功。边界候选缺失、监督、head、上下文、评分与选择仍有耦合；当前研究不通过事后继续扫描来分离所有因素。

## 正式停止范围与未执行项

当前 goal 以 **NO-GO** 终态收尾：不再对同一 M25R 固定路线追加 phase-only 全量推理、R1 阈值/权重/半径/距离扫描、额外 seed 或同方案续训；不重复已有完整 GPU 前向。保留原始数据、checkpoints、失败证据、全部纠正结果、缓存、A/B 预测和代码快照。

Setaria 为 **NOT_ACCESSED / NOT_RUN**，不是失败或零分。原 full-versus-ablation 为 **NOT_EXECUTED**；本次 B−A 不是原 ablation gain，不能冒充 learned structural heads 的独立增益。开发必要门槛已失败，无理由为了完成表格而打开盲测或消耗消融计算。

长上下文、不同 heads/监督、结构化训练或全新 decoder 是尚未检验的另一个研究提案，不属于这个 NO-GO 的普遍否定，也不是当前收尾后自动启动的工作。任何未来重启须有实质新机制、独立预算与预先冻结评价，不能只是为了过门槛调整本次操作点。

本轮通过内置浏览器与 ChatGPT Pro 做了有限方法学商议；其意见用于检查“固定方案失败”与“普遍不可能”的措辞边界，不构成对远程产物的读取、独立代码测试或科学验证。原消融的含义仍以项目合同为准，不采用宽泛的“新解码消融”替代。

## 复现与审查入口

权威恢复目录：`outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY/`。

| 产物 | 用途 |
|---|---|
| `A_replay_check.json` | A 的冻结指标与完整链/phase 复现 |
| `B_prediction_complete.json` | 四方向候选/预测数量、固定目标分数、reference-free 生成记录 |
| `B_evaluation.json` | pooled/逐物种指标、错误诊断、门槛与停止判定 |
| `A_to_B_chain_changes.json` | 完整链加 phase 的增删清单 |
| `A/B_models.json`、`A/B_predictions.gff3` | 实际预测，而非仅有指标 |
| `A/B_predictions_structural_audit.json` | 全量独立合法性审计 |
| `recovery_tests.xml`、`logs/M25RGRAPHRECOVER_13225219.out/.err` | 测试及执行日志 |
| `reproduction_code.tar.gz` | 本次代码、配置、测试、恢复脚本快照 |
| `run_graph_run_before_recovery.py` | 原接口故障版本，便于审查最小修复 |

基础 Git HEAD 为 `01fbbf15544f18ef0fca19688adc1aa0d3284248`，但存在未提交改动，单独 checkout 该 HEAD 不足以复现；应在既有项目中采用本次源代码快照。原缓存位于 `outputs/M25R-E1-RUN-GRAPH-R1/scores/`，checkpoint 和开发数据仍在原位置；快照不包含数据或权重。环境为既有 `generanno`；原前向固定 torch 2.5.1、transformers 4.49.0、peft 0.19.1，6144-bp 原切窗、bf16→float32→float16 保存，图内 float64 累加。

实际提交命令为 `sbatch sbatch/M25R-E1-RUN-GRAPH-RECOVERY.sbatch`。脚本明确先测试/A 复现，再 B 解码、独立评价。**不要向已有结果目录重投**：脚本拒绝已有 A/B 模型，保存函数也拒绝覆写。若以后仅为复现，应由调用者指定新的输出路径，复制原 A 三个产物，并在 CPU Slurm 作业内依次调用：

```bash
python scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py evaluate --arm A --output-dir NEW_OUTPUT
python scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py decode --output-dir NEW_OUTPUT --score-source-run outputs/M25R-E1-RUN-GRAPH-R1
python scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py evaluate --arm B --output-dir NEW_OUTPUT
```

以上是复现说明，不代表另一次运行已提交。历史原 R1 失败记录、phase-only 上界、oracle 报告与勘误报告保留；它们的“尚未执行/等待批准”属于当时状态，由本终态报告和 `run_graph_recovery_execution.md` 更新，不回写历史预测或把旧 STATUS 改成成功。
