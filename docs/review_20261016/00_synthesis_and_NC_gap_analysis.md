# 00 项目总梳理与 Nature Communications 差距分析（2026-10-16）

本文件汇总同目录 01–05 五份子报告（均由只读盘点生成），并给出发表差距与推进方案。数值出处见各子报告。

| 报告 | 内容 |
|---|---|
| [01_code_inventory.md](01_code_inventory.md) + [01_path_mapping.tsv](01_path_mapping.tsv) | 代码/工程盘点、最终管线依赖闭包、建议目录结构与旧→新映射（251 条规则覆盖 1,874 文件） |
| [02_remote_inventory.md](02_remote_inventory.md) | Baobab 远端结构、checkpoint、环境、git、credits、ont |
| [03_research_history_M1_M27.md](03_research_history_M1_M27.md) | M1–M27 时间线、方案谱系、正/负结论清单、评估体系 |
| [04_M28_current_state.md](04_M28_current_state.md) | M28 子实验谱系、当前方法、最佳结果、草稿 claim 证据表、未完成项 |
| [05_literature_positioning.md](05_literature_positioning.md) | 竞品（Carbon-A、Helixer、Tiberius、ANNEVO、GeneCAD、PlantGeneAnn、GENATATOR 等）与发表证据标准 |

---

## 1. 项目一句话现状

冻结的 GENERanno-eukaryote-1.2b-cds-annotator-preview 编码器 + 366k 参数角色头（BiGRU+卷积，15 个逐碱基角色）+ 完整 CDS 语法下的全集合配分函数 CRF（精细标签硬约束、粗区域标签走观测通道、未知位置边缘化），推理时对整条 contig 做精确 MAP 并导出 GFF3。

| 设置 | 我们 | Helixer | Tiberius | ANNEVO | Carbon-A（同主权重） |
|---|---:|---:|---:|---:|---:|
| DEV（拟南芥+水稻各 1 条染色体，7,728 参考，历史复用）pooled exact-chain F1 | **70.107**（fresh-joint） | 62.851 | 71.805 | 75.703 | 未跑 |
| TRAIN（头部拟合用 6 条染色体）| **72.958**（rank-4 LoRA，未上 DEV） | — | 72.655 | 75.542 | **78.172**（过滤）/ 76.063 |

- 单 seed、无置信区间；无任何独立物种/盲测；两物种均在 GENERanno/Carbon-A 训练快照参考与 Helixer 训练集中。
- Setaria 盲测未消耗，但 Setaria 已在 Tiberius 训练列表、Helixer 作者测试集中，作为"新颖盲测"价值有限。
- 早期正面结论（M9–M19 粗粒度低 FPR 编码区检测，gene-body F1 0.88）与大量负面结论（跨物种固定调用器、CRF/Tversky/阈值校准无增益、M25R 结构头 no-go、codon 能量有害）证据较扎实但均为开发级（无 L3 确认性证据）。

## 2. 工程现状（影响一切后续工作）

1. **代码不完整地分散在三处**：Baobab 工作树（事实来源，HEAD `01fbbf1` + 924 条未提交改动）、GitHub `origin/main`（`d77d225`，含 `src/m28/core.py`、`training.py`、`export.py`、`c0_decode.py`、`candidates.py` 等最终管线必需模块）、本地（落后 origin/main 39 个提交，M28 全部 untracked）。本地无法运行最终管线；3 个脚本两边都找不到，大概率只在 Baobab。
2. **checkpoint 只在 scratch，无备份**：M28 LoRA `outputs/M28-FULL-GRAMMAR-JOINT-LORA-FIT/step_009216.pt`（5.7 MB）、dense-QV（231 MB）等；M10–M24 的 LoRA 均未保存，无法复载。`/home` 已 98%。
3. **代码按 Slurm 作业号门禁组织**（199 个脚本硬编码 job id），核心逻辑散落在 `scripts/` 里互相 import；无 CLI/pyproject/README/LICENSE/环境文件；117/151 测试本地不可 import。
4. 建议结构：`cdsgram/`（暂名）= `src/cdsgram/{data,model,decoder,train,eval,cli}` + `tests/` + `slurm/` + `experiments/<Mxx>/`（证据归档）+ `benchmarks/` + `paper/` + `archive/`（A5–A44 等废弃分支，不删）。详见 01。

## 3. 资源现状

- Baobab PI 组余额 **147,263 credits**（< skill 规定的 200k 阈值 → 任何提交前须征得用户同意），计费年到 2026-11-30；jwang 10 月 1–16 日已用 101,800，按此速度约 3 周耗尽。
- ont（RTX 4090 24 GB，免费）可连通，但当时已被他人占 14 GB，且无本项目副本、无 `generanno` 环境。

## 4. 对 Nature Communications 的差距（按严重程度）

### 致命级（不解决基本无法投 NC）

| # | 差距 | 现状 | 需要补的实验 |
|---|---|---|---|
| F1 | **没有独立、多物种的泛化评估** | 2 物种、各 1 条复用 DEV 染色体；同行标准已到 Helixer 45 种、ANNEVO 566 种、Tiberius 多谱系 33 种 | 构建 ≥12 个植物物种的独立面板（单/双子叶、基部被子植物、裸子/蕨类可选、绿藻），优先选**晚于 GENERanno/Carbon-A 训练快照发布**的 assembly（时间留出），并做蛋白同源隔离审计；所有工具同口径运行 |
| F2 | **同主权重的 Carbon-A 全面领先，新颖性与必要性受质疑** | Carbon-A TRAIN 78.2 vs 我们 73.0；Apache-2.0、带 GFF3 流程和基因置信度 | 必须在独立面板上与 Carbon-A 头对头；若精度不胜，论文主张必须换到 Carbon-A 做不到的能力（见 §5） |
| F3 | **论文核心主张从未被直接检验** | 标题/摘要核心是"未知位置边缘化而非填背景"，但该对照**从未做过**；粗通道两面板方向相反；codon 有害 | 关键消融：边缘化 vs 填背景 vs 仅精细标签，在相同训练预算下多 seed 比较 |
| F4 | **统计稳健性缺失** | 全部单次拟合 | 最终配置 ≥3 seeds；按基因/染色体 bootstrap 95% CI；配对检验 |

### 重要级（审稿人大概率会提）

| # | 差距 | 需要补的内容 |
|---|---|---|
| I1 | 最终模型未冻结 | LoRA 仅 TRAIN；需要用 DEV 一次性在 fresh-joint / LoRA / (fresh-joint+LoRA) 中选定，之后冻结 |
| I2 | 长基因与多段基因 | 跨窗口参考 0 召回；多段基因是与 Carbon-A 差距主因（95% 净差）| 需分析或改进（更长上下文/窗口拼接），至少作为明示局限 |
| I3 | 生物学层面证据 | BUSCO 完整度、蛋白水平（与 UniProt/同源蛋白比对）、RNA-seq/Iso-Seq 支持的"参考未注释基因"案例 |
| I4 | 不确定性/后验 | CRF 给出真实边缘后验，但尚未系统评估校准（可靠性曲线、ECE、跨物种稳定性）与 Carbon-A 置信网络（AUROC 0.876）对比 |
| I5 | 计算成本 | 每 Mb 运行时间、显存，对比各工具 |
| I6 | 软件可用性 | pip/conda 包、CLI、HF 权重、Docker、文档、示例数据、CI |

### 次要级

- 草稿与最新结果不同步（摘要/12_publication_strategy 落后），过度表述需修正（例如把 TRAIN 上"超过 Tiberius"当作缩小差距的证据）。
- 评估分母阶梯（6,450 / 7,728 / 22,688）不可互比，论文中需统一。
- 未结项目 A33/A44/E-M 能量池、qvo-LoRA（7 个 sbatch 已写未提交）：建议除 qvo-LoRA 外全部按"探索性未完成"归档，不再投入。

## 5. 论文定位选项

**A（推荐）：概率化完整基因注释 + 部分监督新物种适配。**
利用本方法相对 Carbon-A/Tiberius/ANNEVO 真正不同的结构：(i) 未知位置边缘化使模型可以从**不完整注释**（仅高置信基因、BUSCO 基因、蛋白同源支持的基因、少量 Iso-Seq）在新物种上训练/适配；(ii) 配分函数给出**校准的基因级后验**。核心实验：在独立面板物种上，用 k% 部分标签适配 → 与填背景基线、零样本模型、Carbon-A、Tiberius 比较完整链 F1 随标签比例的曲线；外加后验校准与基于后验的质量控制。即便零样本精度低于 Carbon-A，"用少量部分证据即可超过/追平"是有实际意义且可检验的主张，也直接检验 F3。风险：适配实验需要 GPU 训练预算；若适配增益小则仍回到 B/C。

**B：正面硬拼 SOTA 工具论文。** 继续堆叠精度（qvo-LoRA、fresh-joint+LoRA、更长上下文），在独立面板上需要 ≥ Tiberius 并接近 Carbon-A。按现有差距（DEV 低 ANNEVO 5.6 点、TRAIN 低 Carbon-A 5.2 点）成功概率低。

**C：降档快速发表。** 以现有证据 + F1/F4 的最小补充，投 Genome Biology / NAR Genomics & Bioinformatics / Bioinformatics / Plant Journal 方向的方法+分析论文。

## 6. 推荐推进顺序（方案 A；credits 与 ont 约束下）

| 阶段 | 内容 | 主要算力 | 估计 Baobab credits |
|---|---|---|---:|
| P0 工程整合 | Baobab 工作树快照提交到独立分支 → 合并 origin/main → 生成 `cdsgram` 新仓库（复制，不移动旧文件）→ Baobab 上 CPU 跑 pytest → checkpoint 备份 | 登录节点 + 1 个短 CPU 作业 | <50 |
| P1 独立面板 | 核对 GENERanno/Carbon-A/Tiberius/Helixer/ANNEVO 训练物种清单，选 12–15 个植物 assembly（时间留出优先），下载与标准化，同源审计 | CPU（ont 或 Baobab 小作业） | ~1–2k |
| P2 冻结最终模型 | DEV 一次性评估 LoRA（及 qvo-LoRA，可选），选定后 3 seeds | ont 优先（需先部署 generanno 环境与基模）| 0–6k |
| P3 全面基准 | 我们 + Carbon-A + Tiberius(植物模型) + Helixer + ANNEVO（+BRAKER3 ab initio 可选）跑面板；BUSCO；运行时间 | ont 为主，Baobab 3090 补充 | ~10–25k |
| P4 核心机制实验 | 边缘化 vs 填背景 vs 仅精细；部分标签比例适配曲线；多 seed | ont + Baobab | ~15–30k |
| P5 写作与发布 | 校准分析、图表、稿件重写、软件包/HF 权重/Docker | 本地 | 0 |

12-01 新计费年开始前，Baobab 只做小规模、必须的作业；重计算尽量放 ont 或推到 12 月。
