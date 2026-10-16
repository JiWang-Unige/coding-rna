# M28 R3：真实特征、自由候选与联合监督已接通

2026-09-30。已完成固定原始GLM特征提取、两物种B1自由候选/联合梯度联调和同口径基线重评分。**尚无optimizer更新、新模型准确率或优于基线的证据。** 下一步是有限GPU头吞吐、全染色体输出合并和一次冻结的C0/B1首训；不是继续M25R阈值搜索。

## 实际作业

| Job | 内容 | 状态 | wall time | MaxRSS |
|---|---|---|---:|---:|
| 13291587 | 原始GLM特征初次调用 | FAILED 1:0 | 102s | 7,924,784 KiB |
| 13291660 | 单次接口修复、同12窗 | COMPLETED 0:0 | 50s | 3,353,360 KiB |
| 13291835 | 拟南芥B1真实特征CPU联调 | COMPLETED 0:0 | 47s | 1,182,468 KiB |
| 13291836 | 水稻B1真实特征CPU联调 | COMPLETED 0:0 | 50s | 1,093,064 KiB |
| 13291984 | 两物种共同评价器/缓存重评分 | COMPLETED 0:0 | 50s | 1,049,696 KiB |

前两项gpu021各1GPU，后3项gpu035为0GPU。GLM两次合计0.042222 allocated GPUh，未超原15min累计上限。原失败输出保留，没有原地覆盖。

## 原始GLM的具体失败与修复

固定GENERanno1.2B CDS annotator revision
`b0483c23b6b63787b61a6d3a204a9b517d6ba345`，原backbone，无旧LoRA。
首作业第一次前向因 `GenerannoModel.forward() got an unexpected keyword argument 'use_cache'` 失败，未产生成品特征。
读取该固定revision实际forward接口后，确认它是无KV-cache参数的encoder；移除不支持的kwarg，未改权重、样本、attention实现或数据范围。合同原写的`use_cache=False`是初始调用计划，不是成功接口事实。修复目录为
`outputs/M28-GLM-FEATURE-SMOKE-R3/repair1`。

12个训练窗按物种/染色体/链向首个共享随机draw选定，不按标签或模型输出挑选；每窗24,576bp、四个6,144bp块。无special token的四块token拼接与整窗token化完全相同；每窗BF16特征形状4096×2048。
48次block前向总5.6915s，模型加载3.9332s，峰值GPU allocated 2,374,624,768 bytes，12窗缓存201,651,080 bytes。backbone可训练参数和梯度均0；不把局部计时当整染色体吞吐保证。

## B1自由候选与训练分支

[candidates.py](../../src/m28/candidates.py)只接受DNA与模型输出，无参考参数或注释导入。
GT/GC–AG候选携带跨exon未完成密码子，检查跨边界internal stop。
每1kb/类型8端点、每起点/carry/终点类型8合法exon边、每donor8连接、
每事件/carry2条beam、最终8链/kb、最多128exon，均为预先固定工程限额，未扫描。
加性提案beam再交非加性链头重评分，不称全局最优。

每窗自由候选先落盘，再调用监督代码。12窗实际自由候选共1,316条：
拟南芥548、水稻768。它们分别匹配0/7与1/9个可支持训练正例实例；
这是未训练头、孤立训练窗诊断，**不是DEV性能、候选理论上界或架构否决**。

[training.py](../../src/m28/training.py)按窗归一化local fine/grouped似然及端点、
连接、完整链/null的masked类均衡BCE，之后施加保存的抽样权重。
参考链/正连接只在训练分支注入，来源分free/reference_injected；未知和其他isoform
保护不变。局部、端点、连接、链loss分别对对应头及shared层具有非零梯度；
未知候选和endpoint mask梯度为0。六项新测试两作业均通过。

拟南芥实际梯度窗链标签92负/2未知/1正，连接1463负/56未知/10正。
水稻首个梯度窗链标签73负/1正，但基因是单exon，连接1479个均负；
因此该窗不能证明水稻正连接监督接通。后续固定12窗GPU头计时覆盖所有缓存，
需明确检查两物种的正、负连接梯度，而不是把首窗通过外推到所有情形。

## 主评价已固定

见[共同评价结果](common_ruler_result.md)。两物种参考合计7,728条，全部5组
历史预测在原6,450口径精确复现后才重评分。未重跑基线、未改阈值。
训练不支持的合格参考仍留在分母。旧M25R-B仍显著低于现有基线。

当前结果足以继续有限新模型研发，不能证明新颖性、泛化或Nature Communications
发表能力。Setaria和test未使用；GENERanno预训练暴露仍未解除。
全染色体core归属→去重→冲突选择尚未接入正式caller，不能用孤立窗结果替代。

补充：固定12窗的[GPU头吞吐与正连接梯度](head_timing_result.md)已完成，Job13292011为35s、0.009722GPUh，仍零optimizer更新。水稻三个既有窗口含4/21/3个正连接，实际正/负/未知梯度方向均通过；前述首窗局限已由该补充覆盖。R3新增GPU累计0.051944GPUh。

## 可核对产物

- [GLM特征摘要](r3_features.json)
- [拟南芥联调](r3_arabidopsis_integration.json)、[水稻联调](r3_rice_integration.json)
- [工程合同](r3_integration_contract.md)、[有限头计时合同](head_timing_contract.md)
- 远程完整产物：`outputs/M28-GLM-FEATURE-SMOKE-R3`、
  `outputs/M28-B1-INTEGRATION-R3`、`outputs/M28-COMMON-RULER-R3`。
  DNA、特征和候选明细保留远程，不上传公开仓库。
