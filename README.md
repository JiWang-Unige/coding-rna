# coding-rna：基于基因组语言模型的植物 ab initio 蛋白编码基因注释

本仓库是 GENERanno-eukaryote-1.2b 基因组语言模型 + 完整 CDS 语法 CRF 的 ab initio 基因注释研究项目（里程碑 M1–M28）的研究代码与证据归档。
项目目标、权限与科研边界以 [`AGENTS.md`](AGENTS.md) 为准（`CLAUDE.md` 仅指向它）。

> 2026-10-16 起仓库按"研究仓库重组"整理：只做目录与路径的系统化移动（`git mv`，保留历史），**不是**把 M28 重写成 pip 包。
> 新旧路径对照、代码修改清单与验证结果见 [`docs/REORG_20261016.md`](docs/REORG_20261016.md)。

## 当前状态（截至 2026-10-16 的盘点，详见 [`docs/review_20261016/`](docs/review_20261016/)）

- 当前方法线：**M28 full-grammar**。冻结 GENERanno 编码器 + 15 角色 BiGRU 头（约 366k 参数）+ 完整 CDS 语法下的部分观测配分函数训练 + 整条 contig 的精确 max-only skyline 解码。
- 开发集（DEV，拟南芥+水稻各 1 条染色体，历史复用、**非独立测试**）上最终报告模型 fresh-joint 的 pooled exact-chain F1 为 70.107%；
  同尺下 Tiberius 71.805、ANNEVO 75.703、Helixer 62.851。同一 GENERanno 主权重的 Carbon-A 在 TRAIN 上高于我们所有变体。
- 没有独立物种/盲测评估，没有多 seed 与置信区间；论文主张已在 [`07_revised_claims_and_plan.md`](docs/review_20261016/07_revised_claims_and_plan.md) 中重写，**最终架构还会因超越 Carbon-A 的实验而改变**。
- M28 管线现作为"冻结参考实现"放在 [`pipeline/m28`](pipeline/m28/README.md)；新代码写进 `src/cdsgram/`。

## 目录导航

| 路径 | 内容 |
|---|---|
| [`pipeline/m28/`](pipeline/m28/README.md) | **当前 M28 full-grammar 最终管线（冻结参考实现）**：`src/`（`src.m28` 等模块）、`scripts/`（入口与被 import 的库脚本）、`sbatch/`（按阶段分组的最终作业链）。需 `PYTHONPATH` |
| [`src/cdsgram/`](src/cdsgram/__init__.py) | 新包骨架（占位名 `cdsgram`，src 布局）：`data model decoder train eval cli`，目前为空，不含 M28 代码 |
| [`tests/`](tests/README.md) | 仅保留针对 `pipeline/m28` 内代码的测试（其余测试随各自实验/归档） |
| [`experiments/`](experiments/README.md) | 历史里程碑证据：`M1-M9 / M10-M19 / M20-M24 / M25R / M26 / M27 / M28`（脚本、sbatch、configs、reports、tests、docs） |
| [`benchmarks/`](benchmarks/README.md) | 基线运行器与同尺面板：ANNEVO / Tiberius / Helixer / Carbon-A / native baselines |
| [`paper/`](paper/README.md) | 稿件草稿（`manuscript/`）、图表（`figures/`）、发表策略（`strategy/`）、补充材料 |
| [`docs/`](docs/README.md) | 项目历史文档（台账、路线图、结果日志）、`review_20261016/` 系统梳理、`REORG_20261016.md` |
| [`archive/`](archive/README.md) | 废弃/过时内容（不删除）：M28 pre-full-grammar（A5–A44 等）、备份文件 |
| [`environment/`](environment/README.md) | `requirements-generanno.txt` 与历史环境脚本 |

`outputs/`、`data/`、`logs/`、`refs/` 等运行目录由 `.gitignore` 忽略，仍位于仓库根（Baobab 上 `/home/users/j/jwang/coding-rna/outputs/...`），脚本里对它们的相对路径语义未改变。

## 怎么跑最终管线（Baobab）

训练、评估与重数据处理一律经 Slurm（登录节点只做轻量操作）。环境：conda 环境 `generanno`（见 [`environment/requirements-generanno.txt`](environment/requirements-generanno.txt)）。

```bash
cd /home/users/j/jwang/coding-rna          # 仓库根（outputs/ data/ 在这里）
export PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts"   # sbatch 内已设置；交互/单测时需手动设置
sbatch pipeline/m28/sbatch/01_data/M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION.sbatch   # 以此类推，阶段与参数见 pipeline/m28/README.md
```

单元测试（CPU 即可，需 torch/numpy；本机 Python 3.9 无 torch，无法运行）：

```bash
PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts:$PWD/tests" python -m pytest -q tests      # 或 python -m unittest discover -s tests
```

**重要**：最终链路上的脚本以"Slurm 作业号门禁"串联（脚本内硬编码上游作业号与 `outputs/<STAGE>/{STATUS,JOBID,summary.json}` 检查），
因此它们用于在 Baobab 既有 `outputs/` 上复核/续跑，而不是在全新环境用新作业号从头重跑。详见 `pipeline/m28/README.md` 的"已知限制"。

## 工作约定

- 用简体中文沟通；保留数据划分、sealed/held-out 限制、冻结口径与阈值（Setaria 盲测数据未读取，不得自动放行）。
- 历史文本中的 Auto Research 阶段闸/goal 记录已退役，只保留其中的真实数据限制与结果。
- 旧 `README.md.orig`（2026-09-27 状态）保存在 [`archive/backups/`](archive/backups/README.md.orig)。
