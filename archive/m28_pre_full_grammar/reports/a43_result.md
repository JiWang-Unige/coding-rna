# A43 自然竞争监督：完整 TRAIN 结果

2026-10-07。**科学 FAIL：原五项继续条件全部未通过，两物种完整正确链均为 0。A43 按原停止规则关闭。** 不追加 seed、epoch 或阈值搜索，不重开 A36／A42。作业完成和七项数值证明成立，不等于模型成功。

## 实际结果

下表每格为较强冻结基线 A31_A3 → A43，均为原 1491 个 TRAIN 窗口内的窗口—链配对，不是去重的生物学 gene，也不是完整 raw assembly 的结果。known_nonmatch 与 unknown 分开；后者不作为已知负例。

| 物种 | 已选 | 完整正确 | 已知非匹配 | unknown | ≥5-exon 正确 | 完整参考／长链分母 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | 2193 → 2959 | 1531 → 0 | 347 → 2486 | 315 → 473 | 548 → 0 | 2172／854 |
| 水稻 | 861 → 1975 | 478 → 0 | 251 → 1743 | 132 → 232 | 152 → 0 | 744／232 |
| pooled | 3054 → 4934 | 2009 → 0 | 598 → 4229 | 447 → 705 | 700 → 0 | 2916／1086 |

两个物种“正确严格增加”、两个物种“已知非匹配不增加”、pooled“≥5-exon 正确严格增加”均为 false；实际 `TRAIN_preparation_supported=false`。readout 13424410 与 CPU 13424413 的计数、完整分母及五条件逐项相同，CPU 按实际保存输入重新计数，并非只抄 readout 结论。

本结果否定本次有限离线竞争池拟合在原继续门槛上的收益，不证明所有结构训练必然失败。A43 没有原局部监督，但没有局部单独训练和相同容量的匹配对照，不能将本次退化因果归于这一差异。

## 冻结实验与闭合证据

沿用[准备与一次拟合合同](a43_natural_competitor_preparation.md)：原六条 TRAIN 染色体、1491 个唯一窗口、1536 draws × 3＝4608 次更新、原顺序与 importance、原 R4 gold／unknown 规则。fresh seed0 的 366191 参数表示与十五类发射，固定末 `step_004608.pt`，无研究 checkpoint 续训；预训练 encoder 冻结。有限监督支持 19189 配对，保留原全部正链；不是全语法分区训练或在线挖掘。

自然预测按原完整语法搜索，单位角色势、FP64 累加、null／threshold 0、max128 及原 DNA／ORF／motif／phase／span 规则不变。预测入口不读取 gold、旧候选池或参考标签；全部搜索及预测文件闭合后才连接原 R4 标签。

全部 36 个实际 root（fit、potentials、32 decode 分片、readout、CPU）及其 batch／extern，共 108 条 accounting 行均为 COMPLETED／0:0；STATUS／JOBID、实际配置、来源、完整文件数匹配。以下日期均为 2026-10-07，显式按 UTC 查询，不能混用集群本地 +02 时间。

| 阶段 | 作业 | UTC 起止 | 时长 |
| --- | --- | --- | --- |
| fit | 13424143 | 04:04:47–04:44:29 | 39:42 |
| potentials | 13424395 | 04:44:30–04:54:25 | 9:55 |
| decode32 | 13424400_0…31 | 04:54:27–最末 06:38:52 | 各 1:28:35–1:44:25 |
| readout | 13424410 | 06:38:53–06:44:45 | 5:52 |
| CPU | 13424413 | 06:44:46–06:49:04 | 4:18 |

CPU [核验实现](../scripts/verify_natural_competitor_train.py)的七项实际证明全部 true：原 draw／importance／有限集合 NLL／4608 步；三份有限模型与 AdamW 曝光；末模型发射身份与完整 32 分片；新预测 DNA 语法和直接整窗 FP64 log-probability 评分；同一 A43 势下完整语法目标不低于原自然池重评分目标；原 R4 unknown／较强基线／完整分母；逐物种计数与原五条件重算。

最大直接 log-probability 评分误差为 `2.842170943040401e-14`；完整语法目标减同模型原自然池重评分目标的最小值为 `0.0`。旧池仅用于核验目标下界，不是预测候选输入或新的科学基线。未执行独立神经前向重放或独立完整语法 DP 重放，七证明不扩大到这些范围。

权威产物位于 canonical `/home/users/j/jwang/coding-rna/outputs/`：

- `M28-A43-NATURAL-COMPETITOR-FIT`：实际 4608 行训练记录与三份 checkpoint。
- `M28-A43-NATURAL-COMPETITOR-TRAIN-POTENTIALS`：1491 份发射／DNA 及对应收据。
- `M28-A43-NATURAL-COMPETITOR-TRAIN-DECODE/shard_00` 至 `shard_31`：原 ordinal modulo32、合计 1491 行自然预测。
- `M28-A43-NATURAL-COMPETITOR-TRAIN-READOUT`：完整 `natural_scores.jsonl`、`windows.jsonl`、`summary.json`。
- `M28-A43-NATURAL-COMPETITOR-RESULT-CHECK/summary.json`：实际七证明、逐物种计数与五条件。

固定监控曾在 shell／字段检查失败时过早报告“结构已核验”；该声明不作为证据。主对话已按实际 STATUS／JOBID 文件和配置字段完成一次正确的轻量结构核验，至 07:02:05 UTC 已确认全部 36 根闭合。该核验没有读取 tensor 或预测行内容，不替代上述 Slurm CPU 数值核验。旧 cohort 已消费，原同一心跳暂停，不重报旧作业。

## 下一机制的准备边界

[local_structured_emission.py](../src/m28/local_structured_emission.py)只组合原局部 fine／grouped-region loss 与有限结构 loss，保留各自 mask、分母和 unknown 保护；原 importance 由未来调用者只乘一次。两项系数必须显式传入，无默认研究权重，没有修改 A43 训练或推理路线。

唯一人工检查作业 13425410 于 06:51:05–06:51:14 UTC 运行，root／batch／extern 均 COMPLETED／0:0；[七项检查](../tests/test_m28_local_structured_emission.py)全部 PASS，测试耗时 1.644 秒。实际 private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，STATUS／JOBID 和唯一日志一致，batch MaxRSS 474072K。只用人工 logits、单-exon 竞争几何及模拟 R4 元数据，没有真实 catalog／缓存／checkpoint／模型拟合／optimizer。

这一 PASS 只支持新损失接线。下一科学实验需先固定原 R4 局部监督来源与 mask、两个分母、显式权重和同容量／同曝光的局部单独训练对照，再判断结构竞争是否带来收益；本轮未选择权重或启动新 fit。raw／DEV、独立物种、Setaria、A33 及公开发布仍未放行。完整 raw-FASTA 模型、独立泛化、原生对标优势及论文目标均未达到。
