# M25R-E1-RUN-GRAPH-R1：终态核实与有限恢复提案

核实日期：2026-09-26（Europe/Zurich）。本记录不授权修复、重跑、B 解码、训练或 Setaria 访问。

## 当前结论

**原作业工程失败，R1 科学结论未定。** 全量前向与 A 预测导出完成，A 独立结构审计通过；评价接口异常导致冻结指标/完整链集合复现未完成。B 未执行。不得将这次异常称为模型指标不达标，也不能据此关闭整个 GENERanno 路线或宣布 GO。

## 本次读取的权威证据

Baobab 项目：`/home/users/j/jwang/coding-rna`。通过用户当前已配置的 `baobab` SSH 别名及 Bamboo ProxyJump 连接；直连 login1 的 SSH 端口拒绝连接。这只是连接路径变化，不是作业失败原因。

| 项目 | 实际证据 |
| --- | --- |
| GPU 作业 | `12522183`，Slurm `FAILED`，ExitCode `1:0` |
| 起止时间 | 2026-09-08 22:00:14 → 2026-09-09 08:46:53，集群显示时间 |
| 资源/用时 | private-teodoro-gpu；1 RTX3090、8 CPU、32 GiB；10:46:39；batch MaxRSS 5793828 KiB |
| 前向完成 | `inference_complete.json`：17384 窗口，38589.4603 秒 |
| 保存的 logits | scores 目录中 2 条染色体 × 双链 × region/boundary/phase，共 12 个 npy 文件，约 2.19 GiB；本次检查目录与完成记录，没有重跑推理或重新扫描全部数组 |
| A 预测 | `A_prediction_complete.json`：2098；拟南芥正/负链 716/673，水稻正/负链 380/329 |
| A 独立结构审计 | emitted=checked=valid=2098，invalid=0，coverage=1.0，全部 11 项失败计数为 0 |
| A 与冻结 R4 的复现 | `A_replay_check.json` 不存在；指标和冻结完整链集合对照均尚未执行完成 |
| B | JOBIDS.json 中 `decode_and_B=null`；不存在 B 预测完成记录；本次未提交作业 |

上述路径均位于 `outputs/M25R-E1-RUN-GRAPH-R1/`。错误日志为 `logs/M25RGRAPHFWD_12522183.err`；阶段日志为对应 `.out`。原始 STATUS 保持 `FAILED_FORWARD_OR_A exit=1`，不覆写成成功。

## 具体故障及覆盖缺口

失败命令是 GPU 脚本最后一个 Python 进程：

```text
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py evaluate --arm A
```

`run_graph_run.py:266–267` 从 `A_models.json` 读回预测。JSON 将每个 CDS 二元 tuple 还原成 list。调用 `train_generanno_structural_heads.py:464–465` 时，`tuple(model["cds"])` 仅转换外层，内层仍是 list；把整条链加入 set 触发：

```text
TypeError: unhashable type: 'list'
```

错误发生在主指标计算入口（`run_graph_run.py:282`），早于读取冻结 R4 指标及完整链集合比较（:285–293）。此前同次运行已完成 A GFF 的独立审计、导出 GFF 与 JSON 模型的一致性检查、开发参考范围检查。结构审计通过不是冻结复现通过。

这是本次新增调用接口的实现遗漏：合成图测试覆盖选路，导出测试覆盖 GFF/codon 与结构审计，却没有覆盖“JSON 保存→读取→原指标函数”的实际跨进程入口。原 10 项测试通过不能支持这一入口已经验证。没有证据表明需要改变模型权重、阈值、评分或候选图。

## 已采取的动作与授权边界

遵守 R1 合同第 6 节：核实日志、产物和故障位置后停止执行。没有修改运行代码，没有补算 A，没有启动 B，没有重投 GPU/CPU，没有创建后台监控，没有访问 Setaria；也没有修改或删除原始产物。整体研究目标尚未完成。

## 提请批准的一次 CPU 恢复（尚未实施）

建议只对工程失败作一次明确例外，不改变冻结科学合同：

1. 在 JSON 输入边界把每个 CDS 区间规范为 `(start, end)` tuple，保持坐标、相位、顺序和评价公式不变；原评价函数不改。增加覆盖真实 JSON 往返后调用原指标函数的合成回归测试，同时覆盖空预测和正/负链。
2. 使用新的独立恢复目录，保留原失败目录和日志。复用原 A GFF/模型与已保存的 logits；不重新前向、不重新训练、不覆盖已有审计文件。避免原 `open('x')` 输出保护造成第二次碰撞。
3. 在一个 Slurm CPU 作业内依次完成回归测试、A 指标与冻结完整链集合对照；只有两者及结构审计全部通过，才按原 R1 算法完成 B 解码、独立评价及 A→B 对照。
4. 资源上限：private-teodoro-gpu、**0 GPU / 8 CPU / 32 GiB / 4 小时**，替代原尚未提交的 CPU 后续步骤，不追加 GPU 前向。任一新失败/超时/不一致立即停止，没有第二次隐含恢复。
5. 原阈值、候选图、评分、门槛、ablation 未执行限制和 Setaria 封存不变。结果用于原 R1 判定；本提案本身不是执行批准，也不授权后台监控或对外发布。

该恢复能回答原 R1 科学问题；反之，仅凭当前工程异常，不能形成“合理机制修正后仍不达标”的证据闭合 no-go。

## 整体研究目标的收尾边界（2026-09-26，既有证据复核）

读取现有报告及原 JSON 小型汇总，不运行新诊断。此次明确区分四类结论：

| 对象 | 可支持的判断 | 不能推出的判断 |
| --- | --- | --- |
| 原冻结 M25R / epoch-1 row601 | 已有开发集 NO-GO：exact CDS interval F1=0.1204141、chain F1=0.3249883、gene-count ratio=0.3252713；大小写勘误不改变这些坐标指标 | 不能用已撤回的“约36%结构非法”解释失败，也不能推出所有 GENERanno 方法无效 |
| 固定 epoch-1 / row601，仅取消 phase 拒绝 | 机制达标 NO-GO：最多增加 1149 条精确链；即使零新增 nonexact，chain F1 上界仍仅 0.5234609，低于 0.55 | 面板 64/64 恢复不代表全开发集可以达标，不值得为这一已被上界否决的问题重新做全量前向 |
| 原 run 支持、允许跨 block 的放宽候选图 | 3585/6450 条参考存在 canonical 完整路径；理想零假阳性 chain F1 上界 0.7144993。可达性尚未排除该机制 | 0.7145 不是实际模型分数，不保证原边界阈值或固定打分能选中这些路径，也不证明其余门槛可达 |
| 固定 R1 无参考图解码器 | 合成测试完成；全量 logits 已保存；真实 B 尚未执行，当前停在 A 评价接口故障 | 不能把工程异常写成 R1 科学 NO-GO，更不能把 oracle 上界写成真实收益或 Setaria 泛化 |

数值来源：原 R4 `epoch_1/diagnostic.json` 的 reproduction/prediction_errors，`M25R-R4-GROUP-PHASE-RETROSPECTIVE/epoch_1/summary.json` 的全量 chosen-chain 分类，以及 `M25R-E1-BLOCK-ASSEMBLY-REACHABILITY/summary.json`。两项上界用已核定计数重算：`2*(1389+1149)/(6450+2098+1149)` 与 `2*3585/(6450+3585)`。不是新增模型评价。

原 R4 的结构勘误已有正式结果：`r4_case_correction_report.md` 已取代 `r4_terminal_review.md` 中“CPU 勘误尚未实施”的历史状态。三 epoch 的输出全部通过纠正后的结构审计；canonical-compatible reference 为 6065/6450。最终文章或结案说明不得重新使用旧的 2666/6450、旧 invalidity 比例及其归因。

完成整个目标仍有下列明确条件，不能靠增写报告消除：

- **若 R1 实际 B 失败：** 根据固定必要门槛、每物种结果、计数/FPR、结构审计和 split/fusion 给出局部 NO-GO；结合原 M25R 与 phase-only 负证据，提交关闭“本项目已定义的直接结构注释路线”的正式决策。范围必须明确，不声称证明所有未来 GENERanno 方案不可能成功。依照 R1 合同，正式路线关闭仍交由用户决定。
- **若 R1 必要门槛全部通过：** 仍没有全模型对消融收益或 Setaria 盲测结果。只能提出原要求的后续验证方案，不能用 B-A 替代 ablation，也不能直接释放 Setaria 或宣布 GO；相关执行/封存放行仍需原有批准。
- **若不恢复此次工程失败：** 可以按用户决定停止投入，但结论只能是“已知固定方案失败、R1 评价未完成”，不能标记为已经证明合理机制修正仍失败的科学 NO-GO。

因此当前不应继续投入 phase-only 达标验证、重复 GPU 前向或同一 R1 的阈值/距离/权重扫描。是否实施上一节的一次 CPU 恢复，是现阶段唯一待用户决定的执行分支；未借本次整理增加训练、试验、监控或结果发布范围。
