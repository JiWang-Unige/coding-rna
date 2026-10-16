# A29 搜索截断候选的监督资格结果

2026-10-06。不能自动把淘汰正确路径的cutoff都当负例。原1708个局部通过但遗漏query中，134个cutoff与另一条eligible primary相容，186个受保护或未决，合计320／1708＝18.74%必须保留或跳过负监督；另外1388个仅为参考相对known_nonmatch。该结果支持下一独立方法显式保护竞争路径，但不证明这种训练能提高完整基因注释性能。

## 实际完成与身份

job13393045 root/batch/extern均COMPLETED0:0，gpu034，3min44sec，batch MaxRSS782872K。CPU-only／0GPU／2CPU／8GiB／15min上限，实际入口198.162735秒。四项unittest全部PASS／0.003秒。STATUS COMPLETED、JOBID13393045及原A28 source13392879匹配；cutoff_annotations实际2021行、references545行。原4247 profiles、352 TRAIN窗、545配对身份，原R4 primary几何／保护对象，以及实际保存cutoff标签／参考汇总的重算核验全部通过。

原全量A28来源没有重跑，原GEN缓存仅在CPU读取定向DNA；query assignment未用于资格，模型／GLM前向、候选搜索及optimizer均0。没有checkpoint、DEV或封存Setaria。标签优先primary相容，其余才考虑保护链、unknown chain span和noncallable；保护原因可重叠，不相加作互斥原因。

## 遗漏query的资格

| 首淘汰位置 | other-primary相容1 | 保护／未决-1 | 参考known_nonmatch 0 | query合计 |
| --- | ---: | ---: | ---: | ---: |
| 全局入口池 | 53 | 87 | 613 | 753 |
| 活跃延伸池 | 44 | 53 | 371 | 468 |
| 完成链池 | 19 | 26 | 234 | 279 |
| 保留beam分数界 | 18 | 20 | 170 | 208 |
| 合计 | 134 | 186 | 1388 | 1708 |

所有134个primary相容cutoff都不是当前目标参考的相容节点，而是另一条原eligible primary，不能因不属于该query的teacher分配就称为假gene。拟南芥1322query为125／149／1048；水稻386query为9／37／340。至少5外显子1182query为96／137／949，不能用long-chain query资格比例作独立长基因性能。

221遗漏参考／窗配对中，203个至少一个已观察cutoff为known_nonmatch，116个全部已观察cutoff为known_nonmatch；剩余18个没有这样的已观察cutoff。只读每次首淘汰的最后保留候选，没有标注所有领先竞争路径。116个“全部已观察cutoff为0”不等于其全部领先路径都是负例；18个没有已观察负cutoff也不等于不存在其他可用负例。多个query及重复窗口配对不是独立gene样本。

自然命中对照另有313个未返回query的首淘汰cutoff，资格为27／47／239；每个原324对照参考仍至少一个query返回。对照的返回query不制造cutoff，不并入遗漏分母，不能把本2021条cutoff统一当遗漏query或假gene。

## 下一独立方法的约束

原累计搜索分数的梯度回放组件已通过五项合成数值／梯度测试，详见search_score_replay_preparation.md。正式竞争loss若采用同一gene-set规则，只能对已知不匹配cutoff施加负排序梯度；其他primary与保护／未决必须明确跳过，不能靠matched-query slot标签暗中转负。若需要所有领先路径而不只是单cutoff，需真正采集和重新标注其节点，不能从本报告补推。

这不是A9已失败前缀BCE的救参，也尚未实现在线当前策略采集、rollout或正式结构训练。下一独立模型须固定采集规则、raw-grammar不支持监督边界、margin和系数，再完成实际数据的数值／梯度预检；仍结合A23／A25完整冻结自然读出决定模型路线。不扩大beam、改旧门槛、加seed或epoch来重释旧FAIL。

完整raw-FASTA原生竞争、独立适用范围、论文级因果消融及创新价值仍未证实。工程闭合与监督资格不是科学PASS，known_nonmatch不证明生物学上不存在未注释gene。
