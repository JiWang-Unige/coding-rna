# Artifact Registry / 产物与目录契约

> 由 `/artifact-registry` 初始化和维护。所有脚本、配置、训练结果、第三方软件输出都应按本契约落位，避免“文件散落，不知道该读哪个”。

## 1. Directory contract

| Category | Path pattern | What goes here | Never put here | Retention |
|---|---|---|---|---|
| Project docs | `docs/*.md` | 决策、计划、结果摘要、投稿策略 | 大模型权重/大数据 | 永久 |
| Per-experiment notes | `docs/experiments/<exp_id>.md` | 单实验结构化记录 | 原始日志全集 | 永久 |
| Evaluator contract | `docs/19_evaluator_contract.md` | 指标定义、评估器脚本、split/claim 可比性合约 | 未验证猜测 | 永久 |
| Baseline reproduction ledger | `docs/20_baseline_reproduction.md` | SOTA/基线复现索引、metric/split/rawness 已核实事实 | 大日志/权重 | 永久 |
| Code review log | `docs/21_code_review_log.md` | 训练前代码审查、blocker、waiver | 完整代码 diff 大段粘贴 | 永久 |
| Framework upgrade log | `docs/22_upgrade_log.md` | 框架升级、兼容修复、新 skill/doc/hook 记录 | 研究结果主记录 | 永久 |
| Source archive | `refs/{pdfs,repos,supp,dossiers}/` | 文献、仓库、补充材料、数据/指标 dossier | 未验证指标结论 | 永久 |
| Reusable scripts | `scripts/` | 通用框架脚本、项目可复用脚本 | 单次实验临时代码 | 永久 |
| Per-run generated scripts | `scripts/experiments/<exp_id>/` | 单次实验专用但会影响结果的训练/eval/数据转换 wrapper | checkpoint、日志、未使用草稿 | 永久，必须被 config/sbatch/docs 引用 |
| Pipeline scripts | `pipelines/<pipeline_id>/` | 流程化 DAG 的 stage 脚本/wrapper | 训练 checkpoint | 永久 |
| Configs | `configs/<exp_id>.yaml` or `configs/pipelines/*.yaml` | run/pipeline 参数 | 结果 JSON | 永久 |
| Slurm scripts | `sbatch/<exp_id>.sbatch` | 提交脚本 | 训练日志 | 永久 |
| Full training state | `runs/<exp_id>/` | checkpoint、tensorboard、训练中间产物 | paper/pdf | 至少到投稿完成 |
| Metric reports | `reports/<exp_id>.json` or `outputs/<exp_id>/metrics/metrics.json` | validate_goal 可读指标 JSON；本项目历史与当前训练脚本主要使用 `outputs/<exp_id>/metrics/metrics.json` | 大日志 | 永久 |
| Outputs/status | `outputs/<exp_id>/STATUS` | run 状态(RUNNING→COMPLETED)、job 对账小文件 | checkpoint | 永久小文件 |
| Code-review gate | `outputs/<exp_id>/code_review_gate.json` | `/code-review-gate` 机器产物：verdict+independence+reviewed_files(sha256)；`pre_submit_gate`/`submit_guard` 据此 deny | 大 diff | 永久小文件 |
| Code-review waive | `outputs/<exp_id>/code_review_waived.json` | 机器记录的代码审豁免 `{"reason":...}`（deny 的可审计 override） | 口头确认 | 永久小文件 |
| Per-run env snapshot | `outputs/<exp_id>/env.txt` | conda/pip freeze + GPU/CUDA（可复现） | — | 永久小文件 |
| Logs | `logs/<exp_id>/` | stdout/stderr、软件日志 | 唯一结果摘要 | 保留关键日志 |
| External software outputs | `software_outputs/<tool>/<run_id>/` | 其他软件生成的原始输出 | 手改后的最终表 | 保留原始+hash |
| Data raw | `data/raw/` | 原始数据只读镜像/链接 | 预处理覆盖版 | 永久/链接 |
| Data interim | `data/interim/` | 中间转换 | 最终发布结果 | 可清理但需可重建 |
| Data processed | `data/processed/` | 可复现实验输入 | 临时scratch | 至少到投稿完成 |
| Analysis notebooks | `analysis/notebooks/` | 探索性图表/诊断 | 生产级 pipeline | 可保留 |
| Manuscript assets | `manuscript/` | 图、表、补充材料草稿 | 原始数据 | 投稿期 |
| Optional worktrees | `worktrees/<exp_id>/` | 多方向同时改共享代码时的 git worktree | secrets/raw data/checkpoints/claim结果 | 临时，合并或丢弃后清理 |

## 2. Exp ID naming
- 正式模型实验：`EXP-A-NNN`（screen）/ `EXP-B-NNN`（full/scale）/ `EXP-M<milestone>-NNN`。
- SOTA 随机重训：`SOTA-<model_slug>-SF<frac>-S<seed>`。
- 流程化 pipeline run：`PIPE-<pipeline_slug>-NNN`。
- 临时探路：`SPIKE-<topic>-NNN`。

## 3. Required run bundle
每个真实 run 至少要有：
- `configs/<exp_id>.yaml`
- `sbatch/<exp_id>.sbatch` 或清楚记录的本地/远程命令
- 若本 run 生成了专用脚本：`scripts/experiments/<exp_id>/...`，并在 config/sbatch/实验文档中引用
- `outputs/<exp_id>/STATUS`
- `reports/<exp_id>.json` 或 `outputs/<exp_id>/metrics/metrics.json`
- `docs/experiments/<exp_id>.md`
- `docs/21_code_review_log.md` 中的 `/code-review-gate` 记录（若本 run 改过训练/eval/data/config/job 代码；full/scale 必须有）
- `docs/06_results_log.md` 中的 `## Result: <exp_id>`

## 4. External software bundle
每次调用第三方软件至少要有：
- `software_outputs/<tool>/<run_id>/command.txt`
- `software_outputs/<tool>/<run_id>/version.txt`
- `software_outputs/<tool>/<run_id>/stdout.log` / `stderr.log`
- `software_outputs/<tool>/<run_id>/inputs.sha256`
- `software_outputs/<tool>/<run_id>/outputs_manifest.tsv`
- `docs/13_pipeline_blueprint.md` 的 ledger 行


## 5. Optional worktree bundle

只有共享代码冲突时使用 `worktrees/<exp_id>/`。必须同时记录：
- `docs/17_parallel_workspace.md` 的 Matrix 行（exp_id、branch、orthogonal axis、owner、status）。
- branch/worktree 路径与 merge/drop decision。
- 训练产物仍按 `runs/<exp_id>/`、`reports/<exp_id>.json`、`outputs/<exp_id>/STATUS` 落位，不放在 worktree 内作为唯一副本。
- 不自动 commit/merge；合并前必须有 smoke/test、code-plan-reviewer 或人工审阅。

## 6. Audit 2026-06-14 — v4.1 post-upgrade placement check

- `python3 scripts/artifact_registry.py --init --list-contract` reported no newly created directories; standard v4.1 directories already exist.
- Central ledgers exist: `docs/19_evaluator_contract.md`, `docs/20_baseline_reproduction.md`, `docs/21_code_review_log.md`, `docs/22_upgrade_log.md`.
- Historical Slurm scripts currently live under `scripts/run_*.sbatch` and `scripts/setup_helixer_container_m1.sbatch`. They are retained in place because old docs/results cite those paths. New Slurm scripts should use `sbatch/<exp_id>.sbatch`.
- `reports/` is currently empty; historical metric JSONs live under `outputs/<exp_id>/metrics/`. This is now an accepted project-specific metric-report location, and `scripts/artifact_registry.py --audit-run` accepts either `reports/<exp_id>.json` or `outputs/<exp_id>/metrics/metrics.json`.
- `docs/experiments/` currently contains only selected experiment notes; future completed experiments should run `/exp-log`, but no bulk backfill was performed in this audit.
- `software_outputs/` currently has no external six-file bundles beyond `.gitkeep`; future pipeline/external-tool runs must use the bundle contract in §4.
