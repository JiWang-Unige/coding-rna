# pipeline/m28：M28 full-grammar 最终管线（冻结参考实现）

冻结 GENERanno-eukaryote-1.2b + 15 角色头；训练 = 完整 CDS 语法下的部分观测（精细/粗区域标签 + 未知位置边缘化）sum-product；
推理 = 整条 contig 的 max-only exon-budget skyline 解码并导出 GFF3。最终报告模型 = fresh-joint（joint 24,576 bp 编码 + fresh 4,608 draws 续训）；
DEV 固定家族 = fresh-joint / uniform-coarse(ε=2/3) / codon。

本目录是**冻结的参考实现**：不再在这里加功能。下一代架构的新代码写进 [`src/cdsgram/`](../../src/cdsgram/__init__.py)。

## 布局

```
pipeline/m28/
├── src/                       # 必须在 PYTHONPATH 上：`from src.m28.xxx import ...`
│   ├── __init__.py
│   ├── m28/                   # 40 个模块：labels / core / full_grammar_* / skyline / codon / lattice 遗留依赖 ...
│   ├── screen_anchor/         # M1-M9 数据/解码库（census.py 依赖 data.assign_splits；与 foundation_probe 一并保留）
│   └── foundation_probe/      # M25 结构头训练库（same_scope.py -> train_generanno_structural_heads）
├── scripts/                   # 94 个脚本，平铺；入口脚本 + 被入口 import 的"库脚本"（裸 import 同级文件）
└── sbatch/                    # 按阶段分组的最终作业链（73 个）
    ├── 00_upstream/           # 数据/标签/特征/C0 pilot 等上游产物的生成作业（产物已在 Baobab outputs/ 中）
    ├── 01_data/               # fresh draws 抽样计划 + 支持性检查
    ├── 02_fit_base/           # coarse-channel 首次拟合（4,608 更新）及核验
    ├── 03_fit_continuations/  # fresh-draws / joint24 / fresh-joint24 / uniform-coarse / codon 续训（preflight -> fit -> verify）
    ├── 04_train_inference/    # 6 条 TRAIN 染色体逐 contig 推理（array 0-5）
    ├── 05_train_evaluation/   # TRAIN 评估、配对比较、four-cell
    ├── 06_dev_family/         # DEV 固定家族（fresh-joint / uniform-coarse / codon）检查-推理-比较-聚合
    └── checks/                # 单元/数值核验作业
```

这与旧仓库的相对布局一致（`src/` 与 `scripts/` 平级），所以 `foundation_probe` 等库里“以其上三级目录为根”的写法无需改动；
但**脚本自身**用 `Path(__file__).resolve().parents[3]` 定位的是仓库根（`outputs/`、`data/` 所在处），这一层数在平铺到 `pipeline/m28/scripts/` 后恰好仍为 3。

## PYTHONPATH 约定

```bash
export PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts"   # 在仓库根执行
# 运行 check_full_grammar_* 和测试时再加 $PWD/tests
```

- `pipeline/m28` 提供 `src.m28.*`、`src.screen_anchor.*`、`src.foundation_probe.*`（`pipeline/m28/src/__init__.py` 是普通包，会胜过仓库根的 `src/`（命名空间目录，内含 `cdsgram`））。
- `pipeline/m28/scripts` 提供脚本之间的裸 import（`import fit_full_grammar_objective_contrast as C` 等）；以脚本方式运行时该目录自动在 `sys.path[0]`。
- 所有 sbatch 已改为 `PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts[:$PWD/tests]"` 并调用 `pipeline/m28/scripts/<name>.py`。
- 若 `pip install -e .`（安装 `cdsgram`）后运行本管线，不会冲突：已安装的是顶层包 `cdsgram`，不是 `src`。

## 最终链路（阶段 → 作业 → 参数）

`$1/$2` 为上游作业号（脚本会核对 `outputs/<STAGE>/JOBID` 与内置期望作业号）。数据准备需要的前置输出（`outputs/M28-DATA-MANIFEST-R2`、`M28-LABEL-POLICY-R4`、`M28-CORE-LABEL-SMOKE-R2`、`M28-PILOT-R4`、`M28-FULL-TRAIN-*`、`M26-*` 参考与 `M25R-*` 的 `outputs/`）均已在 Baobab 上。

| 阶段 | 作业（`sbatch/<stage>/M28-…`） | 主要脚本 |
|---|---|---|
| 01 | `FULL-GRAMMAR-FRESH-DRAWS-PREPARATION` → `-SUPPORT-CHECKS` → `-SUPPORT`（array 0-15）→ `-SUPPORT-AGGREGATE <array_job>` | `prepare_full_grammar_fresh_draws.py`、`check_full_grammar_fresh_draws_support.py` |
| 02 | `COARSE-CHANNEL-CHECKS`、`-FIT-PREPARATION`、`-FIRST-FIT`、`-FIT-VERIFICATION` | `fit_full_grammar_coarse_channel.py {prepare,fit}`、`verify_full_grammar_coarse_channel_fit.py` |
| 03 | `FRESH-DRAWS-{PREFLIGHT,FIT <src>,VERIFY <fit>}`；`JOINT24-{PREFLIGHT,FIT,VERIFY}`；**`FRESH-JOINT24-{PREFLIGHT,FIT <preflight_job>,VERIFY <fit_job>}`**；`UNIFORM-COARSE-*`；`CODON-{PREFLIGHT <checks_job>,FIT,VERIFY}` | `fit_full_grammar_{fresh_draws,joint24,fresh_joint24,uniform_coarse,codon}.py`、`preflight_full_grammar_codon.py` |
| 04 | `*-TRAIN-RAW <fit_job>`（array 0-5）；`ORIGINAL-TRAIN-RAW` | `infer_full_grammar_*_train.py raw`、`infer_full_grammar_original_train.py`（模型加载/窗口编码/解码/GFF3：`annotate_full_grammar_*.py`） |
| 05 | `*-TRAIN-EVAL`、`*-PAIRED`、`FRESH-JOINT24-FOUR-CELL`、`ORIGINAL-TRAIN-DIAGNOSTIC`、`CODON-ENERGY` | `infer_full_grammar_*_train.py evaluate`、`compare_full_grammar_*_train.py`、`aggregate_full_grammar_fresh_joint24.py`、`evaluate_full_grammar_original_train.py` |
| 06 | `JOINT-CONTROLS-DEV-CHECKS` → `-RAW`（array 0-1，每个 branch）→ `-COMPARISON` → `-FAMILY` | `check/infer/evaluate/aggregate_full_grammar_joint_controls_dev.py` |

评估核心：`common_ruler.py`、`eval_gene_body_mask.py`、`eval_structure_diagnostic.py`、`reference_pair.py`/`semantic_reference_pair.py`（同尺：CDS-assessable primary、完整 CDS 链精确匹配）。
基线运行器 `run_native_baseline.py` 被 `evaluate_full_grammar_dev_baselines.py` 直接 import，因此留在本目录；其余基线见 [`benchmarks/`](../../benchmarks/README.md)。

## 内容边界

- 除最终入口外，`scripts/` 还包含**被最终入口 import 的库脚本**，其中不少本是"诊断/对照分支"文件（如 `fit_full_grammar_objective_contrast.py` 是训练循环、`fit_full_grammar_joint24.py` 是 joint 编码、`diagnose_full_grammar_training_observations.py` 是观测目标构造），以及 pre-full-grammar 遗留（`annotate_query.py`、`run_query.py`、`query_*`、`emission_lattice*`、`junction_chain` 等）。它们必须与最终管线放在一起，否则 import 链断裂。
- `src/m28` 还收纳了全部 `full_grammar_*` 库模块（含 E/M 能量池、基因后验、scaled-mass、LoRA/dense-QV 的**库**，对应驱动脚本在 `experiments/M28/full_grammar/`）。
- 分支驱动脚本（LoRA、dense-QV、E/M、native-CDS-prior、complete-core 等）在 [`experiments/M28/full_grammar/`](../../experiments/M28/README.md)；A5–A44 与旧 pilot/B2/R5 等在 [`archive/m28_pre_full_grammar/`](../../archive/README.md)。

## 已知限制

1. **作业号门禁**：约 200 个脚本含形如 `job_id != '13434603'` 的硬编码，并要求 `outputs/<STAGE>/summary.json` 的 `completed/protocol/configuration` 与内置常量一致。在已有 Baobab `outputs/` 上可复核；在新环境以新作业号从头重跑会被门禁拒绝。
2. **protocol 标识字符串未改**：脚本里的 `PROTOCOL = 'reports/M28-METHOD-RESTART/....md'` 被写入 checkpoint 并在加载时逐字比较（例如 `state['protocol'] != PROTOCOL`）。为保持与既有 checkpoint 兼容，这些字符串**刻意保持旧值**；对应协议文件已移到 [`experiments/M28/full_grammar/reports/`](../../experiments/M28/full_grammar/reports/)（pre-full-grammar 的在 `archive/m28_pre_full_grammar/reports/`）。
3. `query_training_census.py` 硬编码 `ROOT = Path('/home/users/j/jwang/coding-rna')`（Baobab 路径）。
4. 本机（Python 3.9，无 torch）只做了静态检查；实际运行需在 Baobab `generanno` 环境验证（清单见 [`docs/REORG_20261016.md`](../../docs/REORG_20261016.md)）。
5. 仍引用已移走路径的“非管线”文件（归档的 sbatch/测试/报告）未被改写，仅保留作历史记录。
