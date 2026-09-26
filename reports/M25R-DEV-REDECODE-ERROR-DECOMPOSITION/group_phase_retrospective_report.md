# M25R/R4 grouping 与 phase 回顾性损失分解

日期：2026-09-08。状态：本次回顾性分析完成；下述新实验仅为待批准方案，未执行。M25R development NO-GO 不变。

## 结论

严重漏检不是单一 grouping 或单一 phase 问题。epoch 1 的候选生成／分组缺失更大；epoch 3 中，已经选择了完全正确 CDS 链、却被 phase 拒绝的数量明显增加，与候选生成缺失成为两项相近的大瓶颈。错误边界导致的 phase 拒绝确实存在，但不能解释这批正确链的失败。现有材料不足以估算取消 phase 后的最终可恢复链数，更不能证明更换 phase head 一定有效。

唯一推荐下一实验：固定 epoch-1 checkpoint、原 tuple 和既有候选的小规模配对 phase-gate 消融；补回必要 logits，测量解除 phase 拒绝后还有多少链通过原 ORF 和边界阈值。它区分“phase 是可兑现的最终瓶颈”与“解除 phase 后仍被后续检查拒绝”，不重新搜索 decoder，也不训练新模型。

## 执行与证据范围

- 唯一作业：`sbatch --parsable sbatch/M25R-R4-GROUP-PHASE.sbatch`，job `12499139`。
- `private-teodoro-gpu` / gpu035，2 CPU、16 GB、30 分钟上限，无 GPU allocation，`--no-requeue`。
- Slurm COMPLETED，ExitCode 0:0，elapsed 00:03:06，MaxRSS 2507940K。实际 0.0517 节点小时／0.1033 CPU 核小时；即按申请上限计也仅 0.5 节点小时。未重投。
- 28 tests passed in 22.47s（既有15项、勘误5项、本次8项）。每个 epoch/物种/方向的原 stage 逐条复现，6450 条 reference 对账通过；phase terminal lineage 总数对账通过；原 R4 与勘误输入文件元数据未改变。
- 仅使用既有 R4 traces、勘误产物及已批准 A/rice development 数据。未加载 checkpoint、未运行模型推理／训练、未读取 Setaria、未提交或推送 Git、未调用第三方评阅；未创建后台监控。
- 这属于 Agent 自查与确定性对账，不是独立评阅。输出 STATUS 的 `COMPLETED_RETROSPECTIVE_REVIEW_REQUIRED` 是脚本完成状态；本报告完成其结果解释，不将其当作旧框架审批要求。

源目录：`outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4/` 与 `outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4-CASE-CORRECTION/`。

新增实现：

- `scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/group_phase_retrospective.py`
- `tests/test_m25r_group_phase.py`
- `sbatch/M25R-R4-GROUP-PHASE.sbatch`

新增结果：`outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE/summary.json`；各 `epoch_1/`、`epoch_2/`、`epoch_3/` 下的 `reference_decomposition.jsonl`、`phase_lineage_evidence.jsonl`、`summary.json`。以下路径均相对于 Baobab 项目根 `/home/users/j/jwang/coding-rna`（实际 BeeGFS 根 `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna`）。

## 1. reference 分母：可核对的互斥损失分解

每列都是同一6450条参考链；分类按表中逻辑优先级，非因果干预。canonical 为勘误后的6065条，非canonical385条。匹配先保留原 R4 emitted-exact-first / greedy-covering 规则，再额外检查任意 lineage 是否已选中完整 truth chain；保留每条 reference 的 original_stage 和原/新 lineage IDs。

| 分类 | epoch 1 | epoch 2 | epoch 3 |
|---|---:|---:|---:|
| 已输出 exact chain | 1389 | 907 | 499 |
| 已选中正确链，phase 拒绝 | 1149 | 1788 | 2673 |
| 已选中正确链，ORF 拒绝 | 1 | 1 | 0 |
| 已选中正确链，边界阈值拒绝 | 3 | 0 | 0 |
| truth 在完整候选集合中，但选错边界，随后 phase 拒绝 | 169 | 356 | 195 |
| truth 在完整候选集合中，但选错边界，其他/混合终态 | 63 | 38 | 20 |
| 整块规则排除可匹配 truth 的连续 run 子序列 | 108 | 196 | 271 |
| 没有 CDS state 支持 | 105 | 218 | 167 |
| 无覆盖整链的块，CDS 支持分散于多个块 | 1634 | 1216 | 681 |
| 无覆盖整链的块，仅部分块跨度支持 | 254 | 352 | 394 |
| 覆盖块 run 过多，也无兼容连续子序列 | 390 | 394 | 427 |
| 覆盖块 run 不足 | 110 | 107 | 151 |
| run 数相同，但现有边界候选不能组成 truth | 1075 | 877 | 972 |
| 合计 | 6450 | 6450 | 6450 |

将最后7行合并为“现有候选生成／分组不能提供完整 truth”，得到3676 / 3360 / 3063，占所有 missing 的72.6% / 60.6% / 51.5%；正确链 phase 拒绝占 missing 的22.7% / 32.3% / 44.9%。这两个集合互斥；其余是选择错误及少量后置检查。前者是当前 state/candidate 系统的观察性归类，**不是单改 grouping 就能恢复的数量**。

108 / 196 / 271 条子序列可达仅表示：现有 runs 的一个连续真子序列，按原 ±6 motif 枚举，能够包含全部 truth 边界；没有更改 radius、没有生成新输出、没有替代 boundary choice，也未检查其未保存的 phase/score。不能将其作为可兑现收益。

新正确链 phase 数比原 stage 的1145 / 1770 / 2595多4 / 18 / 78。这不是改写原结果，而是找到了原一对一 covering 匹配未归到对应 reference 的 exact chosen lineage（包括原块外 ±6 snapping）；原 stage 每条及总数仍完全复现。

### 重要分层

- 正确链 phase 拒绝，Arabidopsis：971 / 1396 / 2067；rice：178 / 392 / 606，两物种均存在。
- multi-CDS 正确链 phase 拒绝：1134 / 1425 / 2090；single-CDS：15 / 363 / 583。后期失败已不局限于多 CDS 衔接。
- multi-CDS exact：346 / 231 / 112（分母4951）；single-CDS exact：1043 / 676 / 387（分母1499）。
- 328条 span>6144 参考链仍无 exact；其中跨块 CDS 支持为270 / 240 / 149，正确链 phase 拒绝为8 / 17 / 43。这说明长链有多种损失，不能单独归因于窗口长度。
- epoch 3 的2673条正确链 phase 拒绝中1条 noncanonical；其余2672条canonical。前三个 epoch 的正确链 phase 数不能无条件当作“合法可输出链”。
- 各 epoch 使用各自原冻结 tuple；以上趋势不是固定 decoder 参数下的训练因果实验，也不是新的 checkpoint 选择。

## 2. lineage 分母：phase head 错，还是错误边界导致拒绝？

下表分母为所有 phase-rejected lineages，不是6450 reference，不能与上表直接相加。对每条 lineage 检查保存的 phase prefix 及第一次失败，重新用所选 CDS 长度计算 decoder 期望；参考关系优先用完整 exact chain，否则用原 R4 assignment。非exact assignment 不是生物学归属证明。

| 第一次 phase 失败的证据 | epoch 1 | epoch 2 | epoch 3 |
|---|---:|---:|---:|
| 整链 exact，正确 CDS start 上预测与参考 phase 不同 | 1149 | 1788 | 2673 |
| 错误链，但该局部 start 正确、期望与参考一致，预测不同 | 523 | 715 | 751 |
| 错误链改变累积阅读框；head 与局部参考一致 | 85 | 72 | 28 |
| 错误链改变累积阅读框；head 也不符合局部参考 | 37 | 43 | 13 |
| 失败位置不在分配参考的 CDS start，真实目标不可识别 | 134 | 201 | 184 |
| 无可分配参考，真实目标不可识别 | 1350 | 1790 | 2117 |
| phase-rejected lineage 总数 | 3278 | 4609 | 5766 |

可以确认：正确整链的1149 / 1788 / 2673条失败是实际的 phase-class disagreement，不能由该链上游边界错误解释。6450参考链的注释 CDS-start phase 与累计长度约定均一致（disagreements=0）。

另外85 / 72 / 28条在所分配参考下明确是“head 预测符合参考、错误链的累积期望不符”；phase 拒绝不是 head 对这个局部 start 预测错。37 / 43 / 13条同时有错误链期望与 head-reference disagreement，是混合证据，不能分摊独立因果贡献。523 / 715 / 751条只证明局部 phase 不符，不能说修 phase 就能使整链正确。

剩余1484 / 1991 / 2301条的位置/参考不可识别；不能强行归为 phase 错误或 boundary 错误。保存内容仅到第一次拒绝，之后的 phase logits、置信度与最终 score-filter 结果未知。

## 3. phase 失败的形态及结构性限制

正确链在第一段 CDS 就失败：36/1149（3.1%）、939/1788（52.5%）、2066/2673（77.3%）。epoch 1 首次失败主要在第二段（752条）；epoch 3 已转向第一段。期望/预测编码为原类别1/2/3，即注释phase+1：epoch 1 的 `2或3→1` 有1055条；epoch 3 的 `1→2或3` 有2374条。这描述预测偏移，不解释为何发生。

代码证据：`src/foundation_probe/train_generanno_structural_heads.py` 中 phase 为 `Conv1d(hidden+5,4,kernel_size=1)`（约675–689行）；forward 将每个6-mer hidden repeat_interleave(6)，再连接逐碱基one-hot（约750–757行）。因此，同一token内相同碱基输入给phase head的向量相同，输出logits必须相同；但训练标签为 `1+((position_in_CDS+gff_phase)%3)`（164行）。同token同碱基可以要求不同类别，这是当前表征到phase输出的可证结构性别名冲突。

R4 推理按不重叠6144bp窗口切分（`redecode_error_decomposition.py::predict_sequence`），窗口长度可被6整除；本分析按方向化序列的token坐标计算冲突。

- 33252个参考 CDS starts 中10654个（32.0%）在同一CDS、同一6-mer内存在“相同碱基、不同phase标签”的另一个位置。
- 正确链首次phase失败位点的上述冲突：435/1149（37.9%）、482/1788（27.0%）、581/2673（21.7%）。
- 已恢复链全部CDS starts中的冲突：462/2052（22.5%）、220/1302（16.9%）、57/639（8.9%）。这不是匹配对照：CDS序号、物种、链长和选择过程不同，不能作效应量。

**解释边界：**这是参考位点上的结构性冲突统计，不是训练集经 structural_mask 后的实际冲突率，也不是phase准确率上限。冲突要求至少某个位置预测不合标签，但不要求CDS-start这个位置必错；21.7%–37.9%不是“由别名冲突造成的失败比例”。它也不能解释其余全部失败或epoch变化。没有训练干预，不能宣称已证实根因是head架构、LoRA、loss权重或过拟合。

## 4. 唯一最小后续实验（待用户批准，不执行）

### 问题与设计

在一个固定、明确为探索性的epoch-1 lineage小面板上，补回缺失分数，配对比较：A原phase gate；B只在这些固定lineages上旁路phase拒绝，其余 chosen CDS、ORF检查、boundary阈值及score计算完全不变。

H1：正确候选的phase拒绝后有较高可兑现恢复率，phase值得优先干预。H2：多数仍败于ORF/冻结边界阈值，当前phase attrition高估了最终收益。错误候选对照用于观察放开phase会同时放出多少错误链。该实验不检验“6-mer别名是因果根因”，也不采用truth标签替换模型输出。

固定epoch1是顺序上最早的诊断点，不用本实验再选checkpoint。预注册面板目标：64条exact-chosen phase failures、32条nonexact phase failures、32条原emitted controls；按物种/方向比例分层，并在第一组覆盖有/无token冲突，两组均公开实际分层数。每层按稳定lineage ID排序取样，无分数排序、不挑恢复成功者。分层清单及不足处理在推理前冻结；不足则停止报告，不凭结果补样。

推理前从原6144bp非重叠窗口索引确定所有所需窗口；必须包含每条lineage全部必要phase/边界score位置及原score公式所需输入（如需整段，则取整段）。复用原相邻context、padding、tokenization、dtype和原checkpoint；不得截取一个新context替代原窗口。窗口去重总数上限64；若覆盖完整面板超限，停止并报告实际需求，不自动改采样或扩大预算。

同一次forward输出供A/B使用，保存所需logits；不重新生成候选，不重算分组，不修改原R4目录。首先复现面板的已存phase类别/检查前缀、chosen boundaries和emitted controls；发现不一致立即停止，不能用漂移的重放支持消融结论。B仅为隔离输出的分析性反事实，不替代正式decoder，不更新任何冻结成绩。

### 成本

零训练；一个固定checkpoint，最多64个6144bp方向化原窗口（393216bp），一次GPU推理+轻量CPU判定。建议批准硬上限：1张满足原环境显存需求的GPU，最多2小时，8 CPU、32GB主存；提交前仍需核实模型实际显存和加载峰值，不能满足则不上作业。

这些是预算上限，不是已经测得的耗时/显存承诺。历史完整R4作业三epoch含1536训练窗口及整条验证染色体，共84:11:06，不能外推此小面板的精确运行时间。包含加载的2小时若不够则停止，不延时重投。新GPU授权独立于本次CPU授权；本次没有执行该预检、面板生成或实验。

### 读出与预设决策

主要读出：64条exact-chosen phase failures中，旁路后有多少通过全部剩余原检查，附逐条ORF/threshold失败理由；错误链32条同样报告旁路后新输出数；32条emitted controls必须保持不变。分别报物种、首段/后续段、token冲突分层。选择富集面板不能估计整体precision/F1、FPR或gene-count ratio，不能与6450分母拼接。

探索性成功标准（建议预先冻结）：A重放一致且controls不变；至少32/64正确链通过其余检查，同时新增错误链不超过8/32。达到只表示值得另行设计phase修复实验，不是M25R科学PASS或部署门槛。该阈值是新的前瞻性研究决策建议，不是已有冻结gate，需用户批准。

停止条件：输入/重放一致性失败、面板覆盖超过64窗口、资源需求超限、作业失败或2小时超时，立即停止不重投；若完整运行未达上述门槛，停止phase优先路线，交付后置失败分解，不追加训练/调参。即使成功，也只交付结果，不自动开展head改造、扩大样本、grouping实验或Setaria评估。

## 最终判断与剩余不可识别项

本次授权目标已完成：逐reference和逐lineage证据已落盘，原冻结结果保留，损失分解与方案已交付。主要解释为候选生成/分组与正确链phase拒绝并存；后期phase拒绝更突出。完整truth-assisted总数仍为 `NOT_RECOVERABLE_WITH_SAVED_ARTIFACTS`。缺失logits、后续phase检查和边界score通过率均未猜测；本报告不提出独立测试、Setaria泛化或论文确认性主张。新的实验等待用户批准，本次到此停止。
