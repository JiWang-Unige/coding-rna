# 20 M34-NTV2-SCALEUP 阶段一：中期说明（训练器实现与 CPU 测试状态，2026-10-21）

状态：**只完成了 CPU 部分**（实现 + 单元测试 + CPU 上的 DP 干跑）。**没有启动任何 GPU 任务（含 GPU smoke），没有训练，没有读 DEV**。阶段一（T2、T10）的结果、与 M31 B / Helixer / Tiberius / ANNEVO / GENERanno 的对比、T10 ≥ 70 判据、其余 5 个 arm 的运行建议，都要等 GPU 阶段之后在本文件追加。

## 0. 排队条件的变更（用户 2026-10-21 决定）

GPU 顺序：M33 → **M35（PlantCAD2 + LoRA）** → 按 NT-v2+LoRA 与 PlantCAD2+LoRA 谁更好决定 M34 骨干 → M34。

- 排队脚本 `m34_queue.sh` / 门控 `m34_gate.py` 已改成：**`~/coding-rna/outputs/M35-PLANTCAD2-LORA/M35_DONE` 与 `~/coding-rna/outputs/M34-NTV2-SCALEUP/GO_STAGE1` 都存在**、且本机没有 M32/M33/M35 进程（含占 GPU 的进程）、连续 3 次检查（间隔 5 分钟）都通过，才会开始 smoke → T2 → T10。`GO_STAGE1` 只能由协调者创建，代码里没有创建它；文件第一个词选择骨干（`ntv2`，空文件默认；或 `pc2`）。
- **排队脚本没有在 ont 上启动**（此前也没有启动过），所以没有需要停止的后台排队进程。GPU 阶段需要时在 ont 上执行：
  `cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts && nohup bash m34_queue.sh > ~/coding-rna/logs/m34-queue.out 2>&1 < /dev/null &`
- 其余 5 个 arm（T4、T6、k50U、k25U、k50F）只实现、测试，队列里不会启动它们。

## 1. 实现（ont `~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts/`；本地副本同路径，只新增，未 commit）

| 文件 | 内容 |
|---|---|
| `m34_train_lora.py` | 多物种训练器 = M31 B 配方（NT-v2 冻结 + q/v LoRA 块 13–27，层 21+27，C0 头，配分函数 CRF，并行 DP 12 进程，两阶段「无梯度前向 → DP → 分块重算反传」）。样本来自预先生成的 `schedule_<arm>.jsonl`（从第 0 个抽样开始；M31 的 A 阶段偏移不适用）；3,072 步 cosine，每 64 步存可恢复检查点，每 768 步存模型文件；`training.jsonl` 每步一行（含按物种的 NLL、infeasible、跳过的全未知窗口）。 |
| `m34_train_common.py` | `CatalogLabeler`：直接用 `catalog/<物种>/transcripts.jsonl` 还原 M34 标注器（不再重新解析 GFF3，每个物种省几分钟），targets 与 FULL / U（未选基因位点 unknown，边缘化）/ F（填背景）逻辑全部复用 `m34_labels.Labeler`；`example()` 与 M31/M30 同构造；`register()` 把 8 个新物种的 `__val8M` 伪 contig 登记到冻结的 M29 辅助函数（序列读取、参考构建），不改任何冻结文件。 |
| `m34_backend.py` | **可切换编码器**（环境变量 `M34_BACKEND`）：`ntv2` = NT-v2-500M（12,288 bp 窗、步长 6,144、层 21+27、q/v LoRA，M31 代码原样调用，环境 generanno）；`pc2` = PlantCAD2-Small（M32 几何：8,192 bp 窗、步长 6,144、单次前向同时得到两条链（RC 严格等变）、autocast bf16、六碱基均值池化、头宽 8,196，环境 plantcad2；报告 17 §1）。 |
| `m34_infer.py` / `m34_val.py` | 推理 + 整区域精确 skyline 解码 + 内部验证打分（Ath/Osa 留出染色体与新物种 `__val8M` 伪 contig；复用 M31/M32 的 `c1_infer`/`c1_val` 的解码与打分函数）。 |
| `m34_ckpt.py` | 可恢复检查点（头、LoRA、优化器、RNG、步数、配置；配置变化则拒绝恢复）。 |
| `m34_supervise.sh` / `m34_infer_retry.sh` | OOM 等中断自动从最近检查点恢复并记录到 `interruptions.jsonl`；GPU 空闲显存不足时等待；连续 4 次无进展则 FATAL。 |
| `m34_stage.sh` | 单个 arm 的全流程（单实例锁）：训练到 768/1536/2304/3072，每个停点推理 + 内部验证；**checkpoint 只按 Ath/Osa 两条留出染色体的 pooled F1 在 {1536, 2304, 3072} 里选**（`m34_select.py`，并列取步数大者）；对未包含在该 arm 的新物种验证段做跨物种迁移评估；**DEV 只评一次**（哨兵 `dev_started.json`/`DEV_DONE`，模型 md5 不一致则拒绝）；再跑 `c1_eval --compare-dev` 与 1 Mb 块 bootstrap。 |
| `m34_smoke.sh`、`m34_smoke_schedule.py` | GPU smoke（≤20 步，条件满足后才执行）：T10 压力调度（每个物种外显子最多的基因和 CDS 跨度最长的基因 + 真实抽样），10 步 → 停 → **恢复** → 再 10 步；T10_k50U、T10_k50F 各 3 步；20 步模型对一段新物种验证伪 contig 做端到端推理 + 打分；写 `smoke/summary.json`（显存、单步时间、infeasible、预测训练时长）与 `SMOKE_OK/SMOKE_FAIL`。 |
| `m34_queue.sh`、`m34_gate.py`、`m34_env.sh`、`wait_stage.sh` | 排队 / 门控 / 环境选择 / Mac 侧事件驱动等待器（返回条件：新队列或阶段事件、新指标行、FATAL、≥3 次新中断、队列进程消失、STAGE1_DONE，或 50 分钟心跳）。 |
| `m34_boot_pair.py`、`m34_cpu_dp_check.py`、`m34_select.py` | DEV 上 T10−T2（及对 M31 B、GENERanno、Carbon-A 过滤）的配对 1 Mb 块 bootstrap；CPU 干跑；checkpoint 选择 / DEV 哨兵。 |

**PlantCAD2 + LoRA 的边界**：`pc2` 后端的编码器、几何、池化、训练/推理路径都已实现并测试，**但 PlantCAD2 的 LoRA 目标（Mamba2 的 `in_proj/out_proj`，且 `mamba_fwd/mamba_rev` 权重共享）由 M35 定义，我没有自行设计**。接口约定：`experiments/M35-PLANTCAD2-LORA/scripts/pc2_lora.py` 提供 `add_lora(model, cfg) -> params`、`lora_state(model)`、`load_lora_state(model, state, cfg)`、`CFG`；该文件不存在时 `pc2` 后端只能 `--lora none`（冻结编码器，只训头），队列在 `pc2` 且该文件缺失时会停下（事件 `pc2_lora_module_missing`）。`pc2` 的输出目录为 `outputs/M34-NTV2-SCALEUP/pc2_<arm>/`，头初始化取 M32 runA 的 `best_step` 模型。M35 的头学习率/LoRA 学习率若与 M31 B（1e-4 / 3e-4）不同，用 `--lr-head/--lr-lora` 传入。

## 2. 测试（ont，CPU 核 25–31，nice；全部通过）

`tests/test_m34_train.py`（NT-v2 后端，19 项）与 `tests/test_m34_pc2_backend.py`（PlantCAD2 后端，4 项）。

- **多物种抽样确定性**：T2、T4 的抽样序列由 `m34_schedule.draw_schedule` 重新生成，与磁盘文件**逐行逐字节相同**；7 个 arm 的调度均满足：24,576 个抽样、步号正确、不含 DEV / 留出 / 验证 / 未用 contig、物种计数在等权重二项分布 5σ 内、窗口 ≤ 24,576 bp；物种集合严格嵌套（T2 ⊂ T4 ⊂ T6 ⊂ T10）；U 与 F 共用同一抽样序列；基因优先抽样只落在该 arm 的已标注基因集合内（k25/k50/FULL）。
- **标签掩码**：`CatalogLabeler` 与解析 GFF3 的原 `Labeler` 在 Brachypodium 验证段上**目录逐项相同、三种处理的 targets 逐位相同**；在 Brachypodium 与 Dioscorea（注释最弱的物种）上验证：U = 未标注基因位点整段 unknown（细/粗标签均为 IGNORE）、F = 背景、FULL = 有标签，已标注基因的标签在三种处理中相同；拟南芥/水稻在部分标签 arm 里始终全标注，新物种标注比例精确为 k（±0.2%）。
- **与 M31 单物种路径在 T2 设置下的一致性**（12 个 T2 抽样窗口，拟南芥+水稻，两条链、基因优先/均匀）：基因组上下文、窗口几何、取向后的上下文字符串、核心偏移**逐位相同**；其中 **10/12 个窗口的 targets 完全相同，同一随机势下的 DP 损失与梯度逐位相同（`==`）**。另外 2/12 个窗口 targets 不同：差异位置 3,096/294,912（1.05%；另测 200 个窗口为 1.7%），**没有一处与 M31 的标签相矛盾**——全部是 M31 冻结 R4 标签因「备选转录本不一致」置为 IGNORE 的位置，M34 目录只含主转录本，所以在这些位置给出主转录本的标签。这是方案 19 §2.3 已记录的唯一标签口径差异，需要写进 T2 对 M31 B 的比较说明里（已叠加 sampler 差异）。
- **检查点恢复**：中途「崩溃」后从第 4 步检查点恢复，结束时头 + LoRA 参数与不中断运行**逐位相同**（含优化器状态与 RNG）；配置变化时拒绝恢复。
- **PlantCAD2 后端**（RC 等变伪编码器，CPU）：窗口几何、取向字符串、**训练路径（窗口势 + 冻结 `owned_context_logits`）= 推理路径（区域拼接）**（两条链，容差 2e-5）、梯度到达头。
- 其他：smoke 调度构造正确（每个物种外显子最多 / 跨度最长的基因，含 > 24 kb 跨度）；门控脚本在缺少标志文件时拒绝放行且不会把自己当成 M32/M33 进程；8 个新物种验证伪 contig 的参考数与目录统计完全一致（565–955 条，合计 5,793），序列长度 8,000,000。
- **CPU DP 干跑**（`m34_cpu_dp_check.py`）：10 个物种的压力窗口（外显子最多至 79、CDS 跨度最长至 168 kb）+ 真实 T10 抽样共 32 个，随机势下 **InfeasibleObservation = 0**；DP 单窗口 13–22 s（4 个进程，CPU 被他人任务占满）。目录中 eligible 基因最多 79 个外显子，未触及 128 的上限。

**未验证（必须等 GPU）**：编码器 + LoRA 的真实前向/反传、显存与步时、GPU 上的恢复、推理脚本的端到端运行（这些由 `m34_smoke.sh` 在放行后检查）。

## 3. 时间估算（沿用方案 19 §4，未更新）

smoke 之前不修改：每个 arm ≈ 22.5–24 h（3,072 步 × 23.4 s + 验证 + DEV），阶段一（T2 + T10）≈ 47 h；M31 B 实测步时 18–35 s，CPU 被占满时 DP 是瓶颈。ont 当前 CPU 负载 ≈ 55–60，上述 CPU 干跑里 DP 单窗口 ≈ 17 s（4 进程）。

## 4. 已知风险与注意事项

1. NT-v2 后端的 GPU 路径沿用 M31 B 的写法，但整套脚本尚未在 GPU 上跑过；第一步是 smoke。
2. `pc2` 后端依赖 M35 的 LoRA 模块；若 M35 选择的 LoRA 设计需要逐块前向（类似 NT-v2 的 `forward_layers`），`m34_backend._pc2_encode` 需要相应调整。
3. T2 与 M31 B 的差异混有 sampler（基因优先）与上述标签口径；T10 − T2 不受影响。
4. 本任务没有改动 M31/M32/M33 的文件和进程，没有改冻结代码；`m34_common.py` 只加了一个守卫（`a1_common` 已加载时不再往 `sys.path` 头部插入 M29-A1 目录，避免 PlantCAD2 的 py3.11 环境被 py3.10 包遮蔽）。
5. 测试注意：用 `~/envs/generanno/bin/pytest`（不要 `python -m pytest`，multiprocessing spawn 会重新进入 pytest）；`pc2` 测试需 `M34_BACKEND=pc2 PYTHONPATH=~/envs/m32_extra ~/envs/plantcad2/bin/python -m pytest ...`（该环境下不使用 DP 进程池）。

## 5. GPU 阶段（待办，另行安排）

1. 协调者在用户决定骨干后创建 `GO_STAGE1`（内容 `ntv2` 或 `pc2`），M35_DONE 出现后在 ont 启动 `m34_queue.sh`（见 §0）。
2. smoke → T2 → T10 自动运行；等待器：Mac 上 `wait_stage.sh`（后台，50 分钟心跳）。
3. 完成后在本文件追加：训练曲线、内部验证选择、DEV 结果（按物种/段数分层、CDS 碱基 F1、背景负担、1 Mb 块 bootstrap CI，T10−T2 配对差）、与 M31 B 63.61 / Helixer 62.85 / Tiberius 71.81 / ANNEVO 75.70 / GENERanno fresh-joint 70.11 的对比、T10 ≥ 70 判据、其余 arm 的建议；并写 `M34_STAGE1_DONE`。
