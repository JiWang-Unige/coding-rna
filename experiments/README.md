# experiments：历史里程碑证据

按里程碑分组归档脚本、sbatch、configs、reports、tests、docs。这里的内容是**证据**，不要求仍能运行（只移动、不改内容；
其中引用旧路径的脚本/sbatch/报告保持原样）。要重跑某次实验，请先对照 [`docs/REORG_20261016.md`](../docs/REORG_20261016.md) 里的旧→新路径映射。

| 分组 | 里程碑 | 子目录 |
|---|---|---|
| `M1-M9/` | screen-anchor 与早期探针（M1–M9）；原 `scripts/` 根目录 75 个文件、`configs/m1_*` | `scripts/ sbatch/ configs/ tests/ docs/` |
| `M10-M19/` | LoRA/校准/跨物种/多支系（M10–M19） | `scripts/<实验名>/ reports/<实验名>/ configs/ docs/` |
| `M20-M24/` | clean 面板冻结、SOTA 错误分析、非 CRF 路线、直接结构诊断 | 同上 + `sbatch/` |
| `M25R/` | 结构头与 M25R 重解码误差分解（含 M25） | `scripts/ reports/ sbatch/ configs/ tests/` |
| `M26/` | 同尺机制比较 | `scripts/ reports/ sbatch/ tests/` |
| `M27/` | 等位基因完整性（paired-native） | `scripts/ reports/ sbatch/ tests/` |
| `M28/` | 见 [`M28/README.md`](M28/README.md)：full-grammar 非最终分支与诊断 | `full_grammar/{scripts,sbatch,tests,reports,src}` |

说明：

- M26 中被最终管线 import 的 6 个模块（`reference_pair`、`semantic_reference_pair`、`partial_semantics`、`qualify_araport`、`same_scope`、`background_attribution`）已随管线移到 [`pipeline/m28/scripts/`](../pipeline/m28/scripts/)，对应的 M26 sbatch 与测试里的旧路径因此失效（见 REORG 文档“已知风险”）。
- 评估/数据库 `src/screen_anchor`、`src/foundation_probe` 同理位于 [`pipeline/m28/src/`](../pipeline/m28/src/)；M1–M25 的脚本若要运行，需要把 `pipeline/m28` 加入 `PYTHONPATH`。
- `.gitignore` 中的 `**/reports/**/predictions/` 继续忽略可重算的预测转储。
