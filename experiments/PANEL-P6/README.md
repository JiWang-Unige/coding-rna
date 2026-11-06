# PANEL-P6：面板基线与 C1 判定

记录与协议：`docs/review_20261016/31_P6_panel_baselines.md`（§2 评分协议写于任何面板分数之前）。

| 文件 | 作用 |
|---|---|
| `scripts/p6_common.py` | 路径、核心 10 物种 |
| `scripts/p6_genomes.py` | 整核基因组 FASTA（软屏蔽版给 Vipsania，大写版给其余工具）；`split` 生成 Carbon-A 逐序列输入 |
| `scripts/p6_lib.sh` | GPU 闸门（GO 哨兵、TE 进程、其他 coding-rna 作业、显存与磁盘）、事件日志、资源守护 |
| `scripts/p6_native.sh TOOL` | ANNEVO / Tiberius / Helixer，命令同 M28 DEV 原生基线 |
| `scripts/p6_carbona.sh` | Carbon-A（M29/M30 管线：GPU 推理 + CPU 后处理） |
| `scripts/p6_vipsania.sh` | Vipsania ft（判据）与 pre（描述），同 M38 |
| `scripts/p6_queue.sh` | ont 串行队列：ANNEVO → Tiberius → Helixer → Carbon-A → Vipsania |
| `scripts/p6_check_gff.py` | 原生 GFF 合法性检查 |
| `scripts/p6_score.py` | `validate`（DEV / set B 复算）、`status`、`convert`、`panel`（全部完成后打分一次，输出 C1 判定与 seed 规则） |
| `results/` | `validation.json`、`genomes_manifest.json`、暴露审计、最终面板结果 |

ont 上的状态：

```bash
ssh baobab 'ssh lab "tail -n 5 ~/coding-rna/outputs/PANEL-P6/events.jsonl"'
```

推送代码：`./push.sh`（经 baobab 跳板）。
