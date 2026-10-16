# A26 相位证据固定门槛失败

2026-10-06。原1491 TRAIN窗／16321自然候选的全量拟合、末读出和修复后的CPU数值核验均完成；工程结果可信，原五项科学继续条件只有三项通过。两物种正确链增加，但已知非匹配链也增加，因此关闭本版本，不进入DEV，不扫描阈值或补训救参。

## 完成与一次有界修复

合并fit/readout13392225的root/batch/extern均COMPLETED/0:0，19分14秒；原4608更新、fresh seed0、最后004608、固定生成器及零阈值不变。原CPU13392226 FAILED/1:0、20秒，检查器错误地要求按windows写出的证据行与按records写出的自然读出同序；原失败产物保留。

唯一CPU-only修复13392531的root/batch/extern均COMPLETED/0:0，23秒，batch MaxRSS853600K。改为验证实际materialization的windows顺序，并继续逐ID核对原候选链和28维向量，未排序后忽略身份，也未改容差。实际全部4608 draw／importance、末head有限值、证据geometry／投影分解、候选标签／来源、选择／profiles和原五gate已重算一致。没有模型或backbone前向、候选重生成、optimizer更新或新checkpoint。

规范证据为outputs/M28-A26-PHASE-EVIDENCE-FIT、...-TRAIN-READOUT与...-RESULT-CHECK-REPAIR-R1。具体接线、资源、失败stderr及提交命令见[执行记录](a26_phase_evidence_execution.md)，原事前条件见[协议](a26_phase_evidence_protocol.md)。固定监控已于2026-10-05 17:50:45 UTC真实交接修复终态，主对话已消费；这次来自manual startup，不冒称自然定时触发。

## 实际固定读出

所有计数是候选／窗口配对，不是独立gene数或F1。A6是事前主比较，A21仅辅助比较；候选池和标签相同，未知链始终未作为负例训练。

| 物种或范围 | 指标 | A6 | A21 | A26 |
|---|---|---:|---:|---:|
| 拟南芥 | selected correct | 809 | 1501 | 1518 |
| 拟南芥 | selected known_nonmatch | 98 | 144 | 148 |
| 拟南芥 | selected unknown | 73 | 169 | 175 |
| 水稻 | selected correct | 262 | 495 | 502 |
| 水稻 | selected known_nonmatch | 60 | 104 | 106 |
| 水稻 | selected unknown | 23 | 51 | 53 |
| 两物种合并 | 至少5exon selected correct | 54 | 665 | 675 |

两物种selected correct严格增加和pooled长链严格增加，共三项PASS；两物种known_nonmatch不增加均FAIL。相对A6，A26正确链合并增加949，但known_nonmatch增加96。相对A21，正确链只增加24、known_nonmatch增加6，长链增加10；辅助增益不能替代失败的主门槛。known_nonmatch是相对固定primary参考的非匹配，不等同生物学假基因，unknown也不计为错误。

## 解释与下一研究决定

全量缓存接线检查通过：FP32 hidden18761121792字节，新增证据1827952字节；实际materialization985.651秒，GPU峰值18893957632字节。新投影非零梯度4599步，末权重norm0.15060359137180174。因此本轮不是证据分支未接入或未更新导致的负结果。

本轮支持的有限判断是：这套固定线性28维相位兼容性证据没有消除自然集合选择的非匹配扩张。它不证明相位信息无用、原hidden没有相位，或所有结构化模型无效；共有参数随后共同学习，也不是只有一个输出项变化的严格因果消融。相位／局部状态／路径解码已有先例，不以组件加入本身作创新主张。

下一机制选择优先接续正在运行的A23相对边界、A25起点监督和A27全量遗漏正确前缀诊断，分别判断生成改善与局部排序／路径竞争。A26不进入DEV，不改变旧五gate，不转为救阈值或继续拟合相同版本。完整raw-FASTA性能仍落后原生对标；独立适用范围、相称消融、暴露审计和论文主张均未完成，研究目标保持原范围。
