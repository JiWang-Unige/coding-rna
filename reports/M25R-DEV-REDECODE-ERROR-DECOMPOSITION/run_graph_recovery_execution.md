# R1 CPU 恢复执行记录

2026-09-26 用户明确要求：此类工程故障不再需要逐次授权，以自动推进完成原 goal 为先；必要时可在内置浏览器与 ChatGPT Pro 商议。这取代原 R1 合同对这类不改变科学方案的工程修复所设的再次审批要求。冻结科学阈值、候选图、模型权重、评价定义及 Setaria 封存仍保持。

原故障目录 `outputs/M25R-E1-RUN-GRAPH-R1` 完整保留。新恢复目录为 `outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY`。原运行脚本另存于新目录的 `run_graph_run_before_recovery.py`，便于复核失败版本；图核心和原指标函数未修改。

修复仅包含：在 JSON 输入边界将 CDS 区间列表恢复成 tuple；增加独立输出目录及只读旧 score 来源参数；新增 4 个正/负链与空/非空预测 JSON→原指标回归用例，以及 1 个只读 score 来源用例。与已有 2 个 GFF 导出用例一起在 CPU Slurm 作业内先验证。

执行顺序：回归测试 → 复制原 A 模型/GFF/完成记录到新目录 → A 独立审计及 R4 冻结指标/完整链集合对照 → 只在 A 通过后从旧 logits 解码 B → 预测落盘并退出后独立进程评价 B。

提交脚本：`sbatch/M25R-E1-RUN-GRAPH-RECOVERY.sbatch`。资源：private-teodoro-gpu，0 GPU / 8 CPU / 32 GiB / 4 小时。关闭作业内 CUDA 可见设备；不训练、不加载 checkpoint、不重新前向、不改变门槛、不访问 Setaria。此脚本取代原尚未执行的 B CPU 后续步骤。后台监控未创建，Git 未提交。

真实指标失败按原科学停止规则处理，不自动扫描阈值、权重或距离。工程故障先定位，再在新授权范围内最小修复；不能把工程失败当作科学 NO-GO。

## 实际执行

- 作业 `13225219`，2026-09-26 08:49:56 Europe/Zurich 开始，节点 gpu034；Slurm 分配 8 CPU、32 GiB、0 GPU。
- 7 项恢复/导出测试通过，pytest 28.16 秒；没有跳过或失败用例。
- `A_replay_check.json`：原 2,098 条完整 CDS 链及 phase 集合一致；interval F1、chain F1、count ratio、FPR 四项误差均为 0；独立结构审计无非法模型。
- B 在不读取参考注释的解码进程中完成 6,188 条预测并落盘，随后启动独立评价。独立结构审计 6,188/6,188 合法。
- 新目录保存 `reproduction_code.tar.gz`，包含此次源代码、配置、测试和恢复脚本。基础 Git HEAD 为 `01fbbf15544f18ef0fca19688adc1aa0d3284248`，但该 HEAD 不包含所有当前修改；复现应采用代码快照，不能只 checkout HEAD。没有提交或推送 Git。
- Slurm 终态 `COMPLETED`，ExitCode `0:0`，09:01:55 结束，elapsed `00:11:59`，MaxRSS `3876300K`；stderr 为空，最终输出 STATUS 为 `COMPLETED_R1`。
- 科学判定为 `R1_NO_GO`：interval F1=0.7345159688（门槛 0.80）、chain F1=0.5303054281（门槛 0.55）、FPR=0.0560561582（门槛 0.020，且触发 >0.030 立即停止）；count ratio=0.9593798450、结构合法性=1.0。
- 全量主评价、逐物种诊断及 A→B 链增删均已完成。完整结论与有限路线关闭见 `run_graph_terminal_20260926.md`。不再调参，不访问 Setaria，未执行原 full-versus-ablation；B−A 不作为原消融收益。
