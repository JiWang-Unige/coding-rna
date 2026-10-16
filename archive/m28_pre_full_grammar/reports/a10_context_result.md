# A10 整窗表示：完整 DEV 结果

2026-10-03。整窗表示相对matched A3改善两物种exact-chain F1，micro增加4.83986个百分点；但拟南芥background span超出固定上限，原科学gate **FAIL**。这是一轮有效的表示机制开发结果，不是原生性能领先、独立泛化或论文PASS。

## 范围与结果可信度

与A3保持同649457参数query模型、seed0从头初始化、原1536draw／1491唯一TRAIN窗／顺序／权重、三遍4608更新和原objective／optimizer。只用末step004608，不选中间权重；同GEN revision的编码从四个6144bp块变为一个24576bp／4096-token窗口。注意力上下文与位置处理同时变化，不是只改变跨块注意力的单因素消融。

原始FASTA直接生成两DEV contig共8690方向窗，无参考／标签／DEV feature cache参与推理；tau0、max_exons128、nearest owner、全局去重、正gain非重叠选择和phase不变。全部7728主参考（5437＋2291）、unsupported／长链／重叠保留；历史6450父参考仅辅助。没有test/Setaria、DEV checkpoint选择、额外seed或阈值扫描。

fit13359631／infer13359725／eval13359727／diagnostic13359766的root、batch、extern全部COMPLETED/0:0，分别31m20s／1h15m22s／29s／14s。评价在读取参考前完成全部实际raw记录、assembly、GFF identity／坐标／phase往返。末CPU receipt13360847为COMPLETED/0:0／13s，进一步从实际GFF与7728行成员表重算：5743唯一预测、重复0、phase全部一致；参考身份／坐标／unsupported及既有A3/native字段未变，所有selected成员符合真实GFF，逐物种及pooled P/R/F1与冻结gate一致。它不做新推理／拟合，不能当独立科研评阅。

## 完整主评价

下表P/R/F1为百分比。背景量是固定reference-background内的预测span bp，不等于生物学假阳性率。

| 范围 | TP／预测／主参考 | Precision | Recall | F1 | background span bp |
|---|---|---:|---:|---:|---:|
| 拟南芥 | 2637／3871／5437 | 68.12193 | 48.50101 | 56.66094 | 84489 |
| 水稻 | 1213／1872／2291 | 64.79701 | 52.94631 | 58.27528 | 181236 |
| micro | 3850／5743／7728 | 67.03813 | 49.81884 | 57.15982 | 265725 |

macro F1为57.46811%。冻结条件要求两物种F1均严格高于原固定A（53.39755%／56.37675%），且两物种background span均不高于原R5_C0（54809／198252 bp）。两项F1和水稻背景通过，拟南芥84489＞54809失败；不能用pooled改善、水稻通过或较弱A3替代这项条件。

| 方法／同7728主尺 | 拟南芥F1 | 水稻F1 | micro F1 |
|---|---:|---:|---:|
| matched A3 | 51.75673 | 53.52360 | 52.31997 |
| 原固定A | 53.39755 | 56.37675 | 54.34270 |
| A10_CONTEXT | 56.66094 | 58.27528 | 57.15982 |
| 原生ANNEVO | 76.26512 | 74.43182 | 75.70250 |
| 原生Tiberius | 74.17714 | 66.58217 | 71.80468 |
| 原生Helixer | 67.45638 | 53.52427 | 62.85117 |

A10水稻F1高于Helixer，但不能称总体或两物种更好；micro仍低于全部三种native。这里是已完成native运行，不用历史cached benchmark替代它们。原生模型及GEN上游训练／validation／未知曝光边界见[native_baseline_result.md](../../../benchmarks/native_baselines/reports/native_baseline_result.md)与[independent_transfer_candidates.md](independent_transfer_candidates.md)；反复使用的这两物种DEV不是独立泛化测试。

## 命中变动与当前瓶颈

A3→A10不是“只多112条”：selected共同3015、新增835、丢失723、两者均未命中3155。拟南芥新增604／丢失486／net+118；水稻新增231／丢失237／net−6。水稻F1提高伴随预测2264→1872，但TP1219→1213；不能说两物种recall均改善。

全部参考的free精确可达从4443→4829（+386），owner从3898→4114（+216），selected3738→3850（+112）。A10未生成2899、owner丢失715、owner非正204、正gain未选60、最终selected3850，合计7728；A3对应3285／545／101／59／3738。free覆盖提升仍伴随更多owner和后续损失，不能据净值认定唯一瓶颈已解决。

ANNEVO精确命中但A10未命中的2033条中，1211未生成、637丢于owner、142非正、43正gain未选。A10与ANNEVO／Tiberius／Helixer共同命中分别3611／3510／3397，native-only2033／1931／1715，A10-only239／340／453。A10命中而三种native都没有精确命中的105条仍是既有注释主参考，不是105个新发现基因，也不是独立泛化证据。

固定描述性分层：1／2–4／5–9／10+参考CDS片段数各1886／2868／1878／1096，selected由A3的1323／1443／730／242变为A10的1272／1532／778／268。单片段反而下降，不能只报告长链收益。span＞24576bp的3条和已记录unsupported的106条全部留在分母；前者两方法均未生成，后者selected4→2。分层只有参考敏感度／成员变动，没有对应子组假阳性分母，不提供子组precision/F1，不证明跨块注意力的独立因果作用。

## 运行成本、决策与未完成工作

本次raw推理wall_loop_seconds为4503.0843，GPU峰值allocated2525017088／reserved2703228928 bytes；包括native编码、query生成、assembly及输出，不是单encoder前向时间。不能与A3的cached推理直接比速度，也未形成与native工具的公平速度排名。

完整整窗表示是可保留的后续研究线索，但本轮按原gate关闭为FAIL：不加epoch、改阈值／背景定义、删难例或用分层重新放行。新A13完整链竞争目标已在首次训练和读取这些A10指标前固定，依据的是A12原TRAIN评分失效；它暂用匹配A6的原特征／生成器／候选池，不把A10更好看的表示偷偷混入。A13仅检验TRAIN准备条件，不能替代完整raw-FASTA DEV；FULL原Q/C0匹配完整训练仍在运行。

论文目标尚未达成：需要完成FULL原生比较、真正符合上游曝光和同源隔离的独立测试、完整推理与成本证据及可识别创新／消融。当前负结果和DEV限制保留，无Git提交或对外发布。

主要证据：outputs/M28-A10-CONTEXT-{FIT,INFER,EVAL,DIAGNOSTIC}，INFER/A10_CONTEXT.gff3，DIAGNOSTIC/references.tsv，outputs/M28-A10-CONTEXT-RESULT-CHECK/summary.json；均在Baobab canonical项目。事前协议：[a10_matched_context_protocol.md](a10_matched_context_protocol.md)。
