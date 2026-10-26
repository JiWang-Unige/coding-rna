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

## 6. GPU 阶段执行记录（2026-10-22 至 10-23；追加）

骨干：M35 报告 21 的主判据成立（DEV Δ +8.61，CI 下界 +7.57 > 0），协调者 10-22 17:56 创建 `GO_STAGE1` = `pc2` 与 `M35_DONE`，18:09 启动 `m34_queue.sh`。pc2 smoke 共放行四次：

| 时间 | 结果 | 原因与修正 |
|---|---|---|
| 10-22 18:06 | 失败（schedule） | `m34_smoke.sh` 在 `O` 赋值前 `mkdir`，`pc2_smoke/` 不存在；调整 `mkdir` 顺序 |
| 10-22 18:19 | 失败（train_t10_a） | `m34_train_lora.py` 里局部变量 `ctx` 遮蔽了编码器上下文（`'numpy.ndarray' object has no attribute 'enc'`）；改名为 `ctx_logits` |
| 10-22 21:38 | 失败（train_t10_a，OOM） | 分块反传保留约 13 GB Mamba 激活；`m34_backend.load` 对 pc2 启用 M35 同款 LoRA 块层检查点（`enable_layer_checkpoint`） |
| 10-22 22:00–22:32 | **全部计算通过，仅资源检查未过** | T10 20 步含断点恢复、k50U/k50F 各 3 步、切片推理与打分、infeasible 0、峰值显存 3.80 GB；步时中位 **64.4 s**（均值 65.8、最大 93.3）> 脚本写死的 `step_s_median < 60` |

**偏离说明（用户 2026-10-23 决定）**：`step_s_median < 60 s` 是按 NT-v2 后端（约 23 s/步）写的算力检查，不是科学判据。pc2 后端的步时在 M35 已记录为 45–110 s（独占 45–66 s），报告 21 据此估 M34 每 arm 56–94 h。用户决定把 pc2 的门槛放宽为 **< 120 s**（`m34_smoke.sh` 按 `M34_BACKEND` 区分，NT-v2 仍为 60 s），其余检查不变；对应最坏情形约 102 h/arm，按实测 64.4 s 约 55 h/arm，阶段一（T2 + T10）约 5 天。方案 19 的训练配置（3,072 步、抽样、判据）未改。重启队列前不重跑 20 步 smoke：对第四次的既有产物（`pc2_smoke/`）按新门槛重新执行脚本末尾的判定块，`SMOKE_OK` 中记录 `step_s_limit`=120，队列据此直接进入 T2。重启队列的时间与后续事件见 `~/coding-rna/outputs/M34-NTV2-SCALEUP/queue.jsonl` 与看板日志。

## 7. T2（pc2）结果与 PlantCAD2 预训练暴露核实（2026-10-24）

**T2**（拟南芥 + 水稻，3,072 步，M34 抽样器）：内部验证 768 / 1536 / 2304 / 3072 步 = 52.58 / 70.75 / **72.67** / 72.39，按规则选 2304 步。768 步的低分来自多报（预测 8,708 条对参考 5,894），1536 步起消失。DEV（只评一次）**73.88** [71.80, 75.63]（拟南芥 75.38、水稻 70.22；1 / 2–5 / ≥6 段 78.56 / 72.92 / 71.69；背景负担 0.19%）。配对差：对 NT-v2+LoRA +10.27 [9.21, 11.46]，对 GENERanno fresh-joint +3.77 [2.80, 4.88]，对 Carbon-A 过滤 −3.53 [−4.41, −2.76]。

**跨物种迁移**（T2 2304 步在 8 个未参训物种的 8 Mb 内部验证区段上，只评一次）：pooled 68.84（P 69.3 / R 68.4）。按 PlantCAD2 预训练暴露分组（名单：报告 10 附录 A，65 个基因组、每属一个、整科排除茄科）：

| 暴露 | 物种（F1） | 分组 pooled F1 |
|---|---|---|
| 同物种在预训练内 | 毛果杨 74.08、二穗短柄草 69.91、胡萝卜 71.21 | 71.7 |
| 同属不同种（Citrus clementina、Dioscorea alata 在内） | 甜橙 73.11、几内亚薯蓣 54.82 | 64.9 |
| 不在预训练内 | 番茄 74.15（茄科整科排除）、椰枣 62.45（无棕榈目）、菜豆 68.89（同族 Vigna、Glycine 在内） | 68.2 |

DEV 物种拟南芥、水稻都在预训练内。三组之间差别与单子叶 / 双子叶混杂（最弱两种都是单子叶），每组只有 2–3 个物种，不能把差别归因于预训练暴露；论文须写明"8 个新物种中 5 个在 PlantCAD2 预训练里同种或同属"。P2b 用的 B 集（Dryas、Eschscholzia）不在预训练内（同科 / 同目近缘）。名单来自 `plantcad/plantcad2-c4096` 数据集与报告 10 附录 A，未直接读论文补充表 1。

## 8. P2b 执行准备（2026-10-26）

本节只记录 P2b（方案 19 §3.1，T10 − T2 在 M30 评估集 B 上的零样本 exact-chain F1 配对差）的脚本、验证与排队状态；**判据、统计量、阈值均未改动**，结果待 `M34_STAGE1_DONE` 之后追加。

**文件**（`experiments/M34-NTV2-SCALEUP/`，冻结的 M29/M30/M31/M32/M33/M35 文件与 `m34_infer.py` 均未改）：

- `scripts/m34_p2b_infer.py`：导入 `m34_infer`，把 `T.register` 包一层，使每次注册（主进程、spawn 出的解码子进程）之后再执行 `M30.register()`，于是 `m34_infer` 能读 7 条 B contig 的序列并解析区域名（与 M33 `m33_infer.py` 模块级调用 `M30.register()` 是同一思路）。用法与 `m34_infer.py` 相同；另有 `--probe CONTIG,...` 只读不跑模型，用于检查主进程与解码子进程看到的序列长度与校验和。
- `scripts/m34_p2b_score.py`：用 `m30_score.Evaluator` 与 `SC.load_gff_chains`（CDS ≥ 96 bp、不做后验过滤，即 M33 `m33_evalB.py` 的计数）对每条 B contig 打分，再做预注册的配对 1 Mb 块自助法。统计量逐项沿用 DEV 的 `c1_boot.py` / `m34_boot_pair.py` / M35 `c1_boot.py`：块 = (contig, 首个 CDS 起点 // 1 Mb)，块集合取参考链与各被比较运行的预测链之并集（本次共 276 块，仅供 sanity；P2b 的块数由 T10/T2 的预测决定）；`default_rng(0)`，2,000 次 `integers(0, n, n)`；每次重采样的同一组块下，两模型、pooled 与两个物种用同一重采样（物种 F1 取重采样块中属于该物种者，与 M35 一致）；95% CI 取 2.5/97.5 百分位。判据写成代码：pooled 差 ≥ +2.0 且 pooled CI 下界 > 0，且两物种点估计差均 > 0。输出 JSON 含每物种与 pooled 的 tp/pred/ref/F1/P/R、每条 contig 计数、各自 CI、配对差与 CI、判据三项检查、`DONE_*.json` 里的模型路径与 step，以及 `m36_zero_shot_arm_Z`（T10 的 B 集运行目录）。
- `scripts/m34_p2b.sh`：ont 驱动，单实例（`flock`）。等 `M34_STAGE1_DONE` 存在，且本用户没有 `m34_train/m34_infer/m34_stage/m34_queue/m34_supervise/m34_smoke/m36` 进程、也没有属于 coding-rna 的 GPU 计算进程（连续 3 次、间隔 60 s 通过）后才开始。同账号的 ab-initio-TE 项目也用 GPU（实测 `te_benchmark` 环境的 `chrom_span_eval.py`），不计入阻塞。开始后：读取 `pc2_T10/best_step.txt` 与 `pc2_T2/best_step.txt`（T2 必须是已登记的 2304，否则拒绝；T10 要求 `ARM_DONE` + `selection.json` 存在），首次使用时把两个模型的 md5 记入 `p2b/model_<arm>.json`，之后任何重启遇到不同文件就拒绝；先 T10 后 T2 对 7 条 B contig 推理（`m34_infer_retry.sh` 式重试，按 (contig, 链) 续跑，GPU 空闲显存 < 5 GB 时等待），再 CPU 打分并写结果。
- `tests/test_m34_p2b_infer.py`（pc2 环境，3 项）、`tests/test_m34_p2b_score.py`（generanno 环境，5 项）。

**验证（ont，仅 CPU，`nice -n 10`，`taskset -c 25-31`；没有跑任何 GPU 作业）**：

- 注册：未包装的 `T.register()` 不认识 `NC_141367.1`（KeyError，作为反例）；包装后 7 条 B contig 的长度、物种与 M30 `lengths.json` 一致，M34 验证伪 contig 仍可读，`NC_141373.1` 与 `NW_028238483.1` 全长读取并通过 `resolve`/`parse_region`；spawn 出的解码子进程里看到的长度和序列 md5 与主进程相同。不设 `PYTHONPATH` 也能导入（`pc2_adapter` 自己把 `~/envs/m32_extra` 加入 `sys.path`）。
- 打分器复现 M33 NT-v2 B 未过滤（M33 的 `evalB_dryas` + `evalB_esch` 目录）：**Dryas 49.92（tp/pred/ref = 4,665/9,181/9,508）、Eschscholzia 44.09（5,404/12,148/12,367），与报告 18 §3 及 M33 结果文件的 tp/pred/ref 逐项相同**；M30 Gate 臂同样与 M33 结果文件一致（Dryas 71.97、Eschscholzia 71.38）。
- 自助法：M33 NT-v2 B 对自身，三个范围的差值精确为 0、CI = [0, 0]。M33 NT-v2 B 对 M30 Gate 只做报告，不是 P2b：pooled −25.03 [−25.94, −24.13]，Dryas −22.05 [−23.16, −20.93]，Eschscholzia −27.29 [−28.40, −26.21]；NT-v2 B 的 pooled F1 46.61，块自助法 CI [45.33, 47.75]。注意 P2b 按预注册采用 DEV 的 1 Mb 块统计量，所以这些 CI 与报告 18 里 M33 用的基因单位自助法 CI（Dryas [48.95, 50.95]）数值不同，不要混用；块自助法下同一运行的 Dryas CI 为 [47.78, 51.73]。
- 驱动：在假根目录（NT-v2 B 与 Gate 的 GFF3 软链当作 T10/T2 的运行结果）上演练了闸门之后的全部步骤（模型登记、跳过已完成推理、打分、清理 T2 的 `pot`、摘要事件、`P2B_DONE`），并演练了两种拒绝：模型文件被改动（退出码 11）与 T2 步数不是 2304（退出码 1）。演练目录已删除。
- 沙盒结果保存在 `~/coding-rna/outputs/M34-NTV2-SCALEUP/p2b_sanity/sanity_M33_NTv2B_vs_self_and_Gate.json`，与真正的 `p2b/` 目录分开。

**GPU 费用**：B 集共 273.3 Mb（Dryas 102.8 + Eschscholzia 170.5）。方案 19 §4 估计约 2.4 h/模型；T2 的 DEV 实测（53.4 Mb，编码 1,396 s、总墙钟 1,500 s）外推约 2.1 h/模型；若按保守的 ≈15 kbp/s（GPU 被共享、CPU 解码拖尾）则 ≈5 h/模型。两个模型合计约 4.5–10 h，先 T10 后 T2；另需约 5–10 min CPU 打分（ont 负载高时偏慢）。每个模型的推理在显存上与 DEV 相同量级（峰值 < 3 GB），解码用 4 个进程。

**输出路径**：`~/coding-rna/outputs/M34-NTV2-SCALEUP/p2b/`：`T10/`、`T2/`（m34_infer 布局，`<contig>/prediction/full_grammar.gff3`、`chains/`、`DONE_*.json`）、`model_T10.json`、`model_T2.json`、`p2b.jsonl`（事件，每行一条 JSON）、`p2b_result.json`（结果与判据）、`P2B_DONE`。日志：`~/coding-rna/logs/m34-p2b*.{out,err}`。查看进度：`tail -n 5 ~/coding-rna/outputs/M34-NTV2-SCALEUP/p2b/p2b.jsonl`。

**M36 的 Z 臂**：T10 的 B 集运行目录 `~/coding-rna/outputs/M34-NTV2-SCALEUP/p2b/T10` 就是 M36 的零样本（Z）臂：同一模型、同一推理与解码路径、同一 ruler、不做后验过滤；其 `pot/` 势能文件默认保留（约 33 GB，供 M36 以后需要 B 集后验时使用，不必重跑 GPU；`P2B_KEEP_POT_T10=0` 可关闭），T2 的 `pot/` 在打分后删除。M36 复用时仍应按各自协议自行核对基模与判据。

**状态**：驱动于 2026-10-26 03:17 CET 在 ont 后台启动（PID 3174339，`setsid nohup bash m34_p2b.sh`），`p2b.jsonl` 第二行是 `waiting`（原因：无 `M34_STAGE1_DONE`、T10 的训练与阶段进程仍在）；启动后没有任何 GPU 工作；T10 训练仍在进行（2026-10-26 02:30 步 ≈ 2,330/3,072，≈ 58 s/步）。
