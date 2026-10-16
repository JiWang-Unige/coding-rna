# benchmarks：基线运行器与同尺面板

| 目录 | 内容 |
|---|---|
| `native_baselines/` | ANNEVO / Tiberius / Helixer 的原生运行与评估脚本、sbatch、测试、报告（M28 同尺口径：同 FASTA、固定默认操作点） |
| `carbon_a/` | Carbon-A（同一 GENERanno 主权重的外部完整流程）的资产准备、运行时验证、TRAIN 对标脚本、sbatch、报告 |
| `panels/` | 物种面板配置（`m12_publication_panel.yaml`、`m13_distance_generalization_panel.yaml`） |

注意：`run_native_baseline.py` 被最终管线的 `evaluate_full_grammar_dev_baselines.py` 直接 import，因此该文件位于
[`pipeline/m28/scripts/run_native_baseline.py`](../pipeline/m28/scripts/run_native_baseline.py)，而不在本目录。
M12B 同尺基线（M10–M19 时期）的脚本与报告在 [`experiments/M10-M19/`](../experiments/README.md)。
