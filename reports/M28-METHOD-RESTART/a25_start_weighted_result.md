# A25 首起点提示：原 TRAIN 三项门槛失败，关闭

2026-10-06。A25完整拟合、1491 TRAIN窗自然读出及CPU核验均完成；原三项生成覆盖门槛全部失败。关闭A25，不加seed、更新次数或阈值补救，不进入DEV。

## 冻结结果

A25从fresh A17局部上下文结构出发，新增128参数首起点提示及系数1的条件首CDS-base辅助CE，共649585参数；初始共有参数和原四loss与A17一致。它同时改变表征与辅助目标，不是二者的独立因果消融。正式比较基准仍是原A3，而不是较弱的A17。原4608更新、末004608和greedy＋beam4自然池不变。

| 物种／模型 | 自然候选 | 正确候选 | known_nonmatch | 保护／未决 | ≥5外显子正确候选 |
|---|---:|---:|---:|---:|---:|
| 拟南芥 A3 | 8909 | 1751 | 6305 | 853 | 604 |
| 拟南芥 A25 | 11768 | 1698 | 8937 | 1133 | 547 |
| 水稻 A3 | 7412 | 629 | 6261 | 522 | 184 |
| 水稻 A25 | 9451 | 617 | 8099 | 735 | 174 |
| 合并 A3 | 16321 | 2380 | 12566 | 1375 | 788 |
| 合并 A25 | 21219 | 2315 | 17036 | 1868 | 721 |

两物种正确候选分别减少53／12；合并正确候选减少65，长链减少67，known_nonmatch增加4470。原三gate均false，TRAIN_generator_supported=false。不能把候选池膨胀或相对失败A17的小改善当PASS，也不能据此计算最终caller F1。

## 执行与监控提示

fit13392060／readout13392061／CPU13392062的root、batch、extern均COMPLETED0:0，分别59:46／44:16／14:19。STATUS及JOBID匹配。原draw／weight／监督、新增首起点支持、自然R4标签／池／gate、末权重与实际梯度核验通过；未用DEV、sealed、Setaria、新GLM前向或CPU模型前向／新搜索。产物为outputs/M28-A25-START-WEIGHTED-{FIT,TRAIN-READOUT,RESULT-CHECK}。

固定监控回传时将new_head_or_backbone_forward误读成new_head_or_backbone_inference。实际原生提示和summary使用前者，值false；后者不是协议字段。主对话直接核对原字段后排除监控误报，没有修改结果、增加兼容字段或重跑。科学FAIL独立于该监控读取错误。原协议见[a25_start_weighted_protocol.md](a25_start_weighted_protocol.md)。
