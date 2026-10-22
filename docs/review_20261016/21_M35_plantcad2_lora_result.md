# 21 M35-PLANTCAD2-LORA 结果：PlantCAD2-Small + LoRA 对 NT-v2 + LoRA 的同条件比较（ont，2026-10-21 至 10-22）

预声明协议：`experiments/M35-PLANTCAD2-LORA/protocol.md`（写于任何训练输出之前；判据与各项设计未改）。范围：给 M32 的冻结 PlantCAD2-Small 配方加 Mamba2 `in_proj/out_proj` LoRA（上 12 块，rank 8，alpha 16，fwd/rev 共享），头从 M32 runA step 2048 热启动，与 M31 B 同条件（1,536 步、同抽样、同损失、同 DP），决定 M34 多物种放大的骨干。单 seed；DEV 为开发面板且两个物种都在 PlantCAD2 的 65 个预训练基因组内，结论只限于"这个配方、这两个物种"。

## 0. 结论先行

| 项 | 结果 |
|---|---|
| **DEV exact-chain F1（只评一次，step 1536 由内部验证选定）** | **72.219**（95% CI 70.05–73.97）；拟南芥 74.05、水稻 67.84 |
| **主判据 Δ = M35 − M31 B（NT-v2+LoRA 63.607），1 Mb 块配对 bootstrap 2,000 次** | **+8.61 [+7.57, +9.80]**，CI 下界 > 0 ⇒ 协议规则 1：**M34 骨干取 PlantCAD2-Small** |
| 对 PlantCAD2 冻结（M32，59.875） | +12.34 [+10.81, +14.29]：LoRA 增益（M31 中 NT-v2 的 LoRA 增益为 +12.7，量级相同） |
| 对 GENERanno fresh-joint（70.107，注释监督骨干 + 同一语法） | **+2.11 [+1.28, +3.05]**，两物种同向（拟南芥 +1.19 [+0.46, +1.91]，水稻 +4.11 [+1.91, +6.13]）：纯自监督植物骨干 + LoRA 首次超过注释监督骨干 |
| 对 Carbon-A 过滤（77.406） | −5.19 [−6.23, −4.36]；拟南芥 −3.93、水稻 −8.19 |
| 分层 F1（1 / 2–5 / ≥6 段） | 77.85 / 71.61 / 68.74（M31 B：75.71 / 63.71 / 53.98；GENERanno：79.85 / 67.75 / 65.81；Carbon-A 过滤：80.49 / 76.23 / 76.49）：增益几乎全在多外显子连接，≥6 段比 NT-v2+LoRA 高 14.8 点 |
| CDS 碱基 P / R / F1 | 96.48 / 92.32 / 94.36；背景负担 83.4 kb（0.30%），为全部对照中最低 |
| 内部验证曲线（NC_003071.7 + NC_089044.1） | step 512：65.66 → 1024：70.64 → 1536：**71.33**（M31 B 同步数 56.14 / 61.73 / 63.26）；1,536 步时仍在上升，按协议不追加训练 |
| 成本 | 步时 45–110 s（GPU 与他人作业共享时 110 s；独占约 45–66 s），比 NT-v2+LoRA 的 18–35 s 慢 2–3 倍；三段训练合计约 21 h 墙钟（含两次内部验证）；DEV 推理 2,800 s + 打分；峰值显存 3.8 GB |

## 1. 判定与对 M34 的含义

- 协议主判据成立（CI 下界 +7.57 > 0），**M34 阶段一用 `pc2` 后端**：PlantCAD2-Small 8,192 bp 窗 / 6,144 步长、最后一层、六碱基均值池化、头从 M32 runA step 2048 热启动、LoRA 零初始化、头 lr 1e-4 / LoRA lr 3e-4（与 M35 相同，M34 默认值一致）。
- 代价：M34 每个 arm 3,072 步，按 66–110 s/步估算 **56–94 h**（方案 19 按 NT-v2 估 23 h）；T2 + T10 约 5–8 天。放宽步数或换骨干都需要修改预注册方案，本报告不建议改。
- 次判据（描述）：两物种配对差同向；≥6 段 F1 领先 14.8 点；背景负担最低；吞吐是主要劣势。
- 对论文的含义：Fig. 2 的"骨干无关"有了第三个点，且顺序为 PlantCAD2+LoRA > GENERanno 全语法 > NT-v2+LoRA > PlantCAD2 冻结 > NT-v2 冻结；与 Carbon-A 的差距从 13.8 点缩到 5.2 点。零样本跨物种仍未测（M33 只测了 NT-v2 B）。

## 2. 执行记录

| 项 | 内容 |
|---|---|
| 训练 | `c1_stage_lora.sh runB 1536 2`：三段 0→512（至 10-22 01:28）、512→1024（02:00–09:10）、1024→1536（09:27–16:36），每段结束做内部验证；无中断、无 NaN；`interruptions.jsonl` 只有 launch/done 事件 |
| checkpoint 选择 | 仅按内部验证 pooled F1：1536（71.33）> 1024（70.64）> 512（65.66）；`best_step.txt` = 1536 |
| DEV | `c1_dev.sh runB` 一次（哨兵 `dev_1536/DEV_STARTED` 17:02），推理 `DONE_ed03de0c.json`（encode 2,800 s），打分 `dev_eval_step1536.log`，结果 `dev_1536/eval_step1536_DEV.json`；配对 bootstrap `boot_dev_1536.json`（`c1_boot.py`，2,000 次，seed 0，同块同重采样） |
| 共享 GPU | M37 装配 B 推理（10-21 22:55 至 10-22 00:16、09:08–09:10）与另一用户的作业曾与本训练共用 GPU；步时随之在 45–110 s 波动，损失与梯度不受影响 |
| 协议偏离 | 无（§偏离记录为空）。2026-10-22 在 `pc2_lora.py` **末尾追加**了 M34 后端合同的包装函数（`CFG`、`add_lora(model, cfg)`、`load_lora_state(model, state, cfg)`），只是对既有 `install_lora/load_lora_values` 的封装，M35 自身的调用形式保持不变，不影响已完成的训练与 DEV |

## 3. 数据与路径

ont `~/coding-rna/outputs/M35-PLANTCAD2-LORA/runB/`：`training.jsonl`（逐步）、`stage_log.jsonl`、`model_{000512,001024,001536}.pt`、`dev_1536/`、`dev_eval_step1536.log`、`boot_dev_1536.json`。代码 `experiments/M35-PLANTCAD2-LORA/{scripts,tests,protocol.md}`。
