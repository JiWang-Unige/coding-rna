# M35-PLANTCAD2-LORA 预声明协议（写于 2026-10-21，在任何 M35 训练、内部验证或 DEV 输出产生之前）

问题：给 M32 的冻结 PlantCAD2-Small 配方加 LoRA，与 M31 B（NT-v2 冻结 + q/v LoRA，DEV 63.61）做同条件比较，决定 M34 多物种放大训练用哪个骨干。

## 设计（固定）
- 骨干：PlantCAD2-Small（88 M，权重 FP32 + BF16 autocast），窗口 8,192 bp / 步长 6,144，最后一层 `last_hidden_state`（1,536 维），一次前向得到双链（RC 严格等变），6 碱基均值池化后接 M28 C0 头（同 M32）。
- LoRA：Mamba2 的 `in_proj`（768→3,224）与 `out_proj`（1,536→768），**24 块中的上 12 块（0 起第 12–23 块）**，rank 8，alpha 16（缩放 2，同 M31 B），A~kaiming_uniform(a=√5)，B=0（初始为恒等）。参数 604,416（M31 B 为 491,520，1.23 倍）。`mamba_fwd` 与 `mamba_rev` 共享同一对 LoRA（保持预训练的 in/out 投影绑定；同时把 HF 加载打断的冻结基权重重新绑定，并在重绑前断言两份权重逐元素相等）。LoRA 在 RCPS 共享的 mixer 里，对正向半通道与 RC 半通道作用相同，反向互补等变保持（测试验证）。
- 初始化：头从 M32 选定 checkpoint（runA step 2048，内部验证 58.74）热启动，与 M31 B 从 A step 2048 热启动同构。LoRA 零初始化。
- 训练：与 M31 B 相同——fresh-joint 抽样（seed 1，抽样序列接在冻结训练已用的 2048×8 个 draw 之后）、同损失/目标/DP、batch 8、AdamW wd 0.01、clip 1、头 lr 峰值 1e-4、LoRA lr 峰值 3e-4（均 cosine 到 1e-5，5% warm-up）、**步数上限 1,536**（M31 B 同）、每 512 步内部验证、早停 patience 2、每 64 步存可恢复 checkpoint。
- **偏离任务文字的一处**：任务写"步数上限与 M32 一致"（M32 为 4,608）。因为头是从 M32 step 2048 热启动，而且要与 M31 B 做"同条件"比较（M31 B = 从 A step 2048 热启动 + 1,536 步），我取 1,536 步。4,608 步会让 M35 多 3 倍训练预算，比较不再同条件。若 1,536 步时曲线仍在上升（同 M31 B），在报告中如实写明，不追加训练。
- 内部验证 / checkpoint 选择：与 M31/M32 相同的整条留出染色体（拟南芥 NC_003071.7、水稻 NC_089044.1，5,894 条参考）；只按内部验证 pooled exact-chain F1 选 checkpoint。DEV（NC_003074.8、NC_089041.1）选定后**只评一次**（`c1_dev.sh` 一次性哨兵）。不使用 Setaria、独立面板物种或 sealed 数据。
- 评估：同 M31/M32 的项目 ruler（M26 CDS-assessable 主参考、exact CDS 链、CDS≥96 bp 规则）；报告 pooled/按物种、按外显子数（1 / 2–5 / ≥6）、CDS 碱基 P/R/F1、背景负担、1 Mb 块自助法 95% CI（2000 次，seed 0），并与 M31 B、M32 冻结、NT-v2 冻结 A 作**配对**差（同块同重采样）；与 Helixer/Tiberius/ANNEVO 仅比既有 DEV 值；GENERanno fresh-joint 与 Carbon-A 过滤用同函数重算。

## M34 骨干决策规则（先于结果固定）
主判据：DEV pooled exact-chain F1 的配对差 Δ = M35（PlantCAD2+LoRA）− M31 B（NT-v2+LoRA），1 Mb 块自助法 95% CI。
1. CI 下界 > 0 ⇒ 推荐 PlantCAD2-Small 作为 M34 骨干（同时报告成本：步时、推理吞吐、8 kb 上下文对 12 kb 的取舍、预训练数据暴露）。
2. CI 上界 < 0 ⇒ 推荐 NT-v2（沿用 M34 已备好的方案）。
3. CI 跨 0 ⇒ 两者无法区分；**默认推荐 NT-v2**（理由：M34 的数据/调度已按 NT-v2 配好；NT-v2 预训练明确排除植物与注释，来源更干净；PlantCAD2 预训练窗口取自已注释基因 ±5 kb 且含拟南芥/水稻），并写明 PlantCAD2 在哪些指标上占优（按物种、≥6 段、背景负担、吞吐），供用户权衡。
次判据（只用于描述，不改变上面的规则）：按物种配对差（拟南芥 / 水稻）、≥6 段 F1、背景负担、墙钟成本。单 seed、DEV 为开发面板（对 PlantCAD2 是样本内），结论只限于"这个配方、这两个物种"，不外推到独立面板。

## 不动的东西
不改 M31/M32/M33/M34 的文件和进程、不改冻结代码；不提交 Baobab 计算、不用 Modal；GPU 只在 `M33_DONE` 存在且无 M33 进程后使用（`scripts/start_after_m33.sh` 门控）。

## 偏离记录
（执行中若偏离，在此追加：时间 / 原因 / 是否已看过相关结果。）
