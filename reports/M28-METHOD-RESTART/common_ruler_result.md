# M28共同DEV评价：先固定口径，再比较新模型

2026-09-30。Job 13291984 COMPLETED，50s，2CPU/8GiB、0GPU；
14项政策/计分测试通过。未训练新模型，未重跑旧工具。

## 同范围、同参考政策

在既有两条DEV染色体复用M26的CDS-assessable定义，见[预先合同](common_ruler_contract.md)。
自身CDS完整、序列可评估，按每gene最长合格CDS选primary；忽略父级gene/transcript
partial但不忽略CDS本身partial。所有方法统一政策，预测原样保留。

| 物种 | DEV seqid | 历史parent-primary | 新CDS-assessable primary |
|---|---|---:|---:|
| 拟南芥 | NC_003074.8 | 4,151 | 5,437 |
| 水稻 | NC_089041.1 | 2,299 | 2,291 |
| 合计 | 同两条染色体 | 6,450 | 7,728 |

新primary无重复链，phase链一致性检查通过。旧5组方法×两物种exact-chain
全部精确复现；拟南芥新5,437分母复现M26已有结果。
新政策与R2“所有CDS中选最长”在拟南芥2个gene的primary结构不同，水稻0。
正式拟合前需处理训练primary政策；不能把两个政策视为完全相同。

非标准splice 32/38、split终端codon 9/5、同向primary重叠12/9、
长于24,576bp 1/2（拟南芥/水稻）仍留在评价分母，非互斥计数。
不按模型训练eligibility删基因。

## 主结果：unique exact CDS-chain

| 方法 | 拟南芥F1 | 水稻F1 | pooled precision | pooled recall | pooled F1 |
|---|---:|---:|---:|---:|---:|
| M25R-A | 0.3402 | 0.2827 | 0.7555 | 0.2051 | 0.3226 |
| M25R-B | 0.5579 | 0.4820 | 0.6018 | 0.4819 | 0.5352 |
| ANNEVO cache | 0.7628 | 0.7445 | 0.7860 | 0.7305 | 0.7572 |
| Helixer cache | 0.6519 | 0.5711 | 0.6201 | 0.6338 | 0.6269 |
| Tiberius cache | 0.7418 | 0.6658 | 0.7326 | 0.7041 | 0.7180 |

M25R-B新政策TP3,724/预测6,188/参考7,728。它仍低于三种缓存基线，
不能把参考政策变化导致的新TP称为模型改善。pooled为micro聚合，不是物种宏平均。
独立CDS-base union和其他assessable isoform匹配另存在[原始JSON](common_ruler.json)；
后者不自动算primary TP。普通准确率不与constrained指标混写。

## 背景诊断，不是生物学真负例

固定背景为两链向所有已注释gene-feature跨度之外，合计27,956,441bp。
这个定义与“主评价reference之外”不同，因此下列FPR**不能与历史5.61%直接比较**。

| 方法 | 预测span覆盖该背景 | 全链位于背景的数量 | 全背景链/Mb |
|---|---:|---:|---:|
| M25R-A | 0.2845% | 115 | 4.1135 |
| M25R-B | 0.9948% | 179 | 6.4028 |
| ANNEVO cache | 1.1532% | 146 | 5.2224 |
| Helixer cache | 1.9252% | 258 | 9.2286 |
| Tiberius cache | 0.7205% | 165 | 5.9020 |

这些是参考相对的输出位置/密度，未由RNA或蛋白独立证据裁决。
不是precision，也不是新增真基因或生物学假基因的证明。

## 后续用途与限制

这套冻结参考清单作为新C0/B1的共同DEV标尺。它解决同范围计分，
不解决未知backbone暴露、基线历史权重与当前版本差异，或独立泛化。
全部输入是现存sequence-only工具输出，不能声称本轮已重跑各工具最新版本。
参考JSONL与原始输出留在远程；test和Setaria未读取。

实现：[common_ruler.py](../../scripts/experiments/M28-METHOD-RESTART/common_ruler.py)；
Slurm：[M28-COMMON-RULER.sbatch](../../sbatch/M28-COMMON-RULER.sbatch)。
