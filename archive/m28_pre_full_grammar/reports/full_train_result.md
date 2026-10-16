# FULL 完整 TRAIN 对照：冻结 DEV 门槛 FAIL

2026-10-04，曼谷时间。Q 与 C0 的完整固定训练、raw-FASTA DEV 推理、共同评价及终态数值核验均已完成。Q 在水稻的 exact-chain F1 高于 C0，但拟南芥下降，两物种背景预测 span 均增加，因此原冻结门槛 FAIL。Q 的 pooled F1 也低于已完成的 ANNEVO、Tiberius、Helixer 原生运行；本轮不能作为论文性能胜出或独立泛化证据。

## 实际执行与协议

原协议是 Baobab `reports/M28-METHOD-RESTART/full_train_matched_protocol.md`。两臂均使用完整自然 TRAIN 的 25512 个窗口、六条 TRAIN 染色体、seed0 fresh、三遍固定顺序，共 76536 次更新。保留两物种等权的原 stratum expected-risk、FP32 AdamW 参数及末 `step_076536.pt`，不按 DEV 挑 checkpoint。GEN revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345` 冻结，原四块 6144bp 编码覆盖 24576bp；训练使用 BF16 材料化特征，头为 FP32。

Q 保留原 A3 known-query-mean 目标及 Q16 协作查询，649457 个可训练参数；C0 为 local-only 头与冻结 native HMM 解码，366191 个可训练参数。两臂共享／分割层初始权重配对，训练曝光和优化器相同，但目标、参数量和输出机制不同，不能称单一组件消融。它也不是“旧 A3 只增加数据”这一因素的独立因果证明。

| 阶段 | Slurm ID／终态 | 实际耗时 | 曼谷时间 |
|---|---|---|---|
| Q fit | 13358083／COMPLETED 0:0 | 8时25分43秒 | 10-03 15:06:29–23:32:12 |
| C0 fit | 13358084／COMPLETED 0:0 | 8时25分42秒 | 10-03 15:06:30–23:32:12 |
| Q raw 推理 | 13358085／COMPLETED 0:0 | 1时04分11秒 | 10-03 23:32:14–10-04 00:36:25 |
| C0 raw 推理 | 13358086／COMPLETED 0:0 | 1时06分35秒 | 10-03 23:32:13–10-04 00:38:48 |
| 原共同评价 | 13358090／COMPLETED 0:0 | 35秒 | 10-04 00:38:49–00:39:24 |
| 终态数值／实际文件核验 | 13361286／COMPLETED 0:0 | 35秒 | 10-04，原评价完成后 |

六项 root／batch／extern 均已核对为 COMPLETED 0:0，STATUS／JOBID 对应。原评价 extern 为36秒；其余耗时列为 root。最后一项为 private／0GPU／2CPU／8GiB／5分钟上限的 CPU 核验，不进行优化器更新或新模型推理。原五项资源、时间上限和停止条件未更改。

两臂各完成完整 8690 个有向 DEV 窗口。推理入口使用原生 raw FASTA，不读取注释或训练特征缓存；输出分别为 Q 6772、C0 7715 条唯一 CDS 链。共同主参考固定为 7728 条：拟南芥5437、水稻2291。长链、重叠或其他 unsupported 参考层仍留在分母，不因方法支持范围缩小参考。sealed／Setaria 未访问。

## 原冻结门槛与完整 DEV 结果

原继续条件为：**两物种 Q exact-chain F1 都严格高于匹配的 FULL C0，且两物种 background span 都不增加**。表中 P／R／F1 为百分比；background 是所有已注释 gene feature span 的并集之外，任一链方向的预测 span 并集 bp，不等同于确认的生物学假阳性数。

| 物种／方法 | TP／预测／参考 | P | R | exact-chain F1 | background span bp |
|---|---|---|---|---|---|
| 拟南芥 Q | 2813／4761／5437 | 59.08 | 51.74 | 55.17 | 102778 |
| 拟南芥 C0 | 3127／5331／5437 | 58.66 | 57.51 | 58.08 | 34545 |
| 水稻 Q | 1267／2011／2291 | 63.00 | 55.30 | 58.90 | 167892 |
| 水稻 C0 | 1171／2384／2291 | 49.12 | 51.11 | 50.10 | 139104 |
| pooled Q | 4080／6772／7728 | 60.25 | 52.80 | 56.28 | 270670 |
| pooled C0 | 4298／7715／7728 | 55.71 | 55.62 | 55.66 | 173649 |

Q−C0 的拟南芥 F1 为 −2.9118 个百分点、TP少314；水稻 F1 为 +8.8066 个百分点、TP多96。background span 分别为 C0 的2.9752倍与1.2070倍。四个冻结条件只有水稻 F1 条件 PASS，其余三项 FAIL，整体 `passed=false`。

pooled F1 虽提高0.6131个百分点，但 TP少218、预测少943，不是两物种召回或 F1 一致改善。物种 macro F1 为 Q 57.0353%、C0 54.0879%，与 pooled 的加总计数口径不同；任何聚合改善都不能替代原逐物种门槛。

CDS-base unstranded union F1：拟南芥 Q 81.46%／C0 96.00%，水稻82.68%／91.31%，pooled 81.8555%／94.5200%；pooled CDS-base recall 为72.2837%／93.4356%。因此不能只用链精度的提高描述 Q 的输出质量，碱基覆盖仍明显不足。

Q 的6772条最终输出都通过现有完整 CDS 结构检查；C0 为7714／7715，拟南芥有一条 `ambiguous_CDS`，保留在评价中而未删除。两臂 GFF 相位均与实际装配记录一致。结构合法不等于参考 exact 命中，也不能据此宣称新基因。

## 与已完成原生基线的同尺比较

下表采用[原生基线完成结果](../../../benchmarks/native_baselines/reports/native_baseline_result.md)，不是 FULL 原评价中附带的旧缓存工具数字。终态核验将原生比较单独加入新 receipt，没有覆盖原 FULL summary 或冻结 gate。7728 条参考的物种／ID／染色体／方向／CDS 坐标与原生参考表逐项一致，background 分母也一致。

| 方法 | 拟南芥 exact-chain F1 % | 水稻 exact-chain F1 % | pooled exact-chain F1 % |
|---|---|---|---|
| FULL Q | 55.17 | 58.90 | 56.28 |
| FULL C0 | 58.08 | 50.10 | 55.66 |
| ANNEVO | 76.27 | 74.43 | 75.70 |
| Tiberius | 74.18 | 66.58 | 71.80 |
| Helixer | 67.46 | 53.52 | 62.85 |

Q 在水稻高于 Helixer，但在拟南芥和 pooled 低于 Helixer；两物种均低于 ANNEVO 与 Tiberius。Q pooled 距 ANNEVO 19.4266、Tiberius15.5288、Helixer6.5753个百分点，不能称达到当前原生基线竞争力。

原生版本／权重来源及曝光限制沿用原生完成记录：ANNEVO author-updated Magnoliopsida，Tiberius2.0.8／angiosperms（实际2.0.7依赖 SIF 中固定安装2.0.8），Helixer0.3.7／a0080。物种预训练曝光因方法而异或未知；GEN上游监督曝光未知，Tiberius作者元数据包含水稻，其他曝光声明也不能代替完整审计。本比较是反复使用的 bounded DEV 上的同尺结果，不是统一无污染 zero-shot 测试或独立泛化证据。

## 已定位的损失阶段，不是修复效果预测

从实际 raw 候选重建后，按同一主参考计算不同阶段“至少出现一次的唯一 exact 参考链”。这是事后诊断，不是可提交的预测或 oracle 评价：

| 物种／方法 | free 候选命中 | owner 后命中 | 最终选中命中 |
|---|---|---|---|
| 拟南芥 Q | 3661 | 3014 | 2813 |
| 水稻 Q | 1554 | 1411 | 1267 |
| pooled Q | 5215 | 4425 | 4080 |
| 拟南芥 C0 | 3145 | 3127 | 3127 |
| 水稻 C0 | 1178 | 1171 | 1171 |
| pooled C0 | 4323 | 4298 | 4298 |

Q 的 free unique 为5215／7728（67.4819%），owner 为4425／7728（57.2593%）：2513条参考未在 free 阶段命中，free→owner损失790条，owner→最终损失345条。早期文字汇总误将3661＋1554写成5226，本记录已按终态 receipt 的逐物种实际计数更正；保存的阶段计数和科学gate未发生改变。

C0对应的损失为25条与0条。Q在候选生成阶段覆盖多于C0，但在owner与选择阶段失去较多正确链，最终命中反而低于C0。它把后续问题定位到“生成覆盖、窗口归属、完整链竞争和背景控制”这些相互关联的阶段；未证明移除owner、降低阈值或改变竞争规则会净改善完整输出，更未授权用DEV扫参救本轮门槛。

## 判断、记录与未执行边界

完整 TRAIN 扩大曝光并没有使原Q目标同时获得两物种改进和背景控制。它排除了“这轮完整训练还没跑完”的解释，但不足以断言单纯训练更久、更多数据或某个结构改动一定有效。原模型、末权重、阈值、owner、解码规则及原失败结果均保留；本轮没有加epoch、挑checkpoint、扫阈值、重训救结果或开启新DEV实验。

固定监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3` 的自然定时检查于曼谷时间00:39:34.159确认原五根链完整，00:39:59.300向主对话交接一次。主对话已消费该终态并完成上述数值核验；同一heartbeat的实际TOML为 `PAUSED`，无新对话／fork／cron。CPU receipt由主对话直接处理，不扩展旧监控对象。平台研究goal当前为 `blocked`，本轮未修改其状态，也不声称后台研究迭代仍在运行；长期论文目标未达成。

没有 Git commit／push、公开发布、作业取消或原始数据／checkpoint／冻结结果覆盖。本轮工程运行完成、结果可信、冻结科学门槛失败、独立泛化未验证，这四项判断分别记录。

## 权威来源

Baobab `/home/users/j/jwang/coding-rna`：

- `outputs/M28-FULL-TRAIN-{Q-FIT,C0-FIT}/{summary.json,training.jsonl,step_076536.pt}`；
- `outputs/M28-FULL-TRAIN-{Q-INFER,C0-INFER}/{summary.json,Q.gff3或C0.gff3,实际raw记录}`；
- 原冻结结果 `outputs/M28-FULL-TRAIN-EVAL/summary.json`及两份GFF；
- 新终态 receipt `outputs/M28-FULL-TRAIN-RESULT-CHECK/summary.json`，实际核验入口 `scripts/experiments/M28-METHOD-RESTART/verify_completed_full.py`、`sbatch/M28-FULL-TRAIN-RESULT-CHECK.sbatch`；
- 共同参考 `outputs/M28-COMMON-RULER-R3/result/{species}.references.jsonl`；
- 原生完整结果 `outputs/M28-NATIVE-BASELINE-EVAL/summary.json`，参考身份表 `outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC/references.tsv`；
- 核验日志 `outputs/M28-PILOT-R5/logs/M28FTRESULTCHECK_13361286.{out,err}`。

CPU receipt在读取参考标签前，先重建两臂17380个实际raw窗口并核对候选计数、装配计数、唯一链、推理／评价GFF的ID／坐标／相位；随后核对全部参考身份、重算实际GFF的exact／CDS-base／background和阶段命中，并复算原冻结gate。所有值与保存结果一致，新增比较使用完整原生运行；该检查是确定性自核验，不冒充独立评阅。
