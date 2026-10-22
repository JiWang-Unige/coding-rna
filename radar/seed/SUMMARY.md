# coding-rna 项目现状（子代理盘点，2026-10-22；总结由主会话落盘）

1. 在做什么：基因组语言模型编码器 + 完整 CDS 语法的配分函数 CRF，做植物 ab initio 蛋白编码基因注释。
2. 里程碑实际到 M37（M35 在跑，M36 只做了 CPU 准备，M37 已完成）；README 仍写 M1–M28，已过时。
3. 定位变了三次：10-16 超过 Carbon-A → 10-17 放弃拼原生精度 → 10-19 确认 GENERanno preview 即 Carbon-A 主权重，核心改为“骨干无关、能用不完整注释训练、输出校准后验的概率化基因语法层”。
4. 10-21 起“目标基因组自举适配”（M36）升为第一主张。
5. paper/ 下稿件是 Codex 时期旧主张，未重写；论点以 docs/paper_outline_20261021.md（C1–C5）为准，C6 自举适配为归纳。
6. DEV（拟南芥、水稻各 1 条复用染色体，单 seed）F1：M28 GENERanno 70.1；NT-v2+LoRA 63.6；PlantCAD2 冻结 59.9；ANNEVO 75.7；Tiberius 71.8；Helixer 62.9；Carbon-A 过滤 77.4。
7. M35（PlantCAD2+LoRA）训练完成，内部验证 71.33（M31 B 63.26）；ont 上在跑 DEV 推理，DEV 只评一次。
8. 卡点：新物种零样本崩塌（Dryas 49.9、Esch 44.1，比 Carbon-A 过滤低 21.6/26.9），问题在多外显子连接；M34 十物种放大是唯一已排队修复，目标 DEV ≥ 70，<65 止损。
9. 正面证据：边缘化比填背景高 16.5/14.2 点（M30）；后验 τ=1 即校准，ECE 0.031，过滤后 DEV +1.68、新物种 +2.3–2.5（M33）。
10. M37 多装配联合解码按预注册判据 NO-GO（headroom 10.3% < 20%）；副发现势平均 +1.53/+2.42 F1。
11. 算力：全部在 ont 串行；Baobab 12-01 前不提交；不用 Modal；ont 与另一用户及 TE 共享。
12. 排队：M35 DEV → 定骨干 → M34 smoke/T2/T10（~47 GPU h）→ M36（~28–30 h）→ 其余 5 臂（~115 h）→ +2 seed → 面板推理。
13. 独立面板 12 物种数据已下载冻结（331,106 条参考，md5 一致），尚无任何推理与基线。
14. 等决定（5）：M34 骨干并创建 GO_STAGE1（最急，等 M35 DEV）；势级集成是否另立协议；主应用选 isoform 还是泛基因组；面板基线放哪跑（拟）；目录改名时机。
15. 已决定（5）：主线目标；放弃拼精度；不用 Modal 且 12-01 前不提交；GPU 顺序；自举升为第一主张。
16. 来源：本地记忆 4 文件；coding-rna-reorg 的 README、AGENTS.md、docs/review_20261016/00–26、paper_outline、paper/；Baobab 只读；ont 进程/GPU/M35 metrics.jsonl。
17. 没读透：部分 review 文档只读开头；M28 子实验用摘要；experiments/*/results 未逐个复算。
18. 待核对：M35 pc2_lora.py 就绪；M28 checkpoint 备份；Baobab 余额（10-16 为 147,263）；21 号文档不存在；CADB 是否覆盖面板物种；5 个 GenBank 面板注释质量。
19. 竞品数字为 05 号文档转述的作者自报值，与 exact-chain F1 不可直接比。
20. 图号 F6、F7 与 claim_C6 为适配看板归纳，大纲未分配。
