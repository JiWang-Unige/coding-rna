# 链头精度诊断D1：差异随cuDNN TF32设置改变，暂不采用优化

2026-09-30。Job **13293846 COMPLETED 0:0**，private/gpu035，RTX3090×1、
2CPU/8GiB，申请5分钟，实际35s = **0.009722 allocated GPUh**，
MaxRSS1,403,016KiB。按[运行前合同](chain_precision_contract.md)执行，
一次TRAIN固定输入，8个链头前向＋一个CPU float64链头参照，
**零optimizer更新、无新权重、无DEV/test/Setaria**。原GPU R1失败保持不变。

## 新的观察

同一份shared特征h、同一checkpoint、同一份192条自由候选，原逐链/packed各在
cuDNN TF32开/关下运行两次。matmul TF32始终关闭；作业开始和每格flags均保存。
GPU环境实测PyTorch2.5.1、CUDA12.1、cuDNN90100，未开autocast。
以下原容差固定为atol2e-6/rtol2e-4，未事后改变。

| 链头比较 | nonadditive最大绝对差 | 超限元素 | 最终gain超限元素 |
|---|---:|---:|---:|
| 同为TF32开：packed vs原逐链 | 0.0009055138 | 17/192 | 6/192 |
| 同为TF32关：packed vs原逐链 | 0.0000009537 | 0/192 | 0/192 |
| 原逐链：关 vs开 | 0.0000004768 | 0/192 | 0/192 |
| packed：关 vs开 | 0.0009052753 | 17/192 | 6/192 |

新对照复现了R1相同的17元素/最大绝对差模式，并明确差异在nonadditive分项；
这是新作业的观察，不补造R1缺失的具体key日志。additive、null、additive-only
在四格逐值相等。同格两次重复也逐值相等，但两次重复不足以证明普遍确定性。

相对同h、同原始权重转换的CPU float64链头参照：
原逐链TF32开、原逐链关、packed关的nonadditive最大误差分别
5.32454e-7、5.32454e-7、6.77440e-7，原容差均无超限；
packed开为9.05372e-4，17个超限。参照不是全模型FP64或生物学真值。

四格均选择候选index43，正gain index为0/5/14/43，没有符号或选择翻转。
最近零gain仍约0.8065，因此这不是接近操作点的苛刻检验，更不能保证其他窗口不翻转。

[完整原始summary](chain_precision_D1.json)、
[每格全部分项向量](chain_precision_D1_cells.jsonl)、
[CPU float64参照](chain_precision_D1_cpu_fp64.json)、
[原始stderr](chain_precision_D1_error.log)保留证据。
stderr的GRU权重非连续警告来自每格deepcopy，涉及临时压紧/内存开销；
本作业不据此或7.35s内部循环时间作任何吞吐结论。

## Pro第21轮及解释边界

Pro报告完整读取公开提交7722ff5bd7ec08d62438058cbf41da32be201adb中的
GPU R1报告/JSON/脚本与packed_chain_gru.py，未访问原始GPU产物、未重跑。
它认为一次≤5分钟有限诊断值得，但不支持单窗采用；其答复在D1提交后完成，
并不是提交授权。讨论仍在[同一会话](https://chatgpt.com/g/g-p-6a29d586630481918525796032225f68-ji-wangke-ti/c/6ab777ce-bcd4-83eb-abff-89a8f9361948)。

Pro建议进一步固定逐exon嵌入z、记录GRU末态。本D1固定的是上游h，
z由同权重、同h、同matmul设置重新计算，未独立保存z/GRU末态。
代码中链头除GRU外是pooling、Linear/GELU和读出，未重新执行shared GRU/卷积；
结果足以支持该窗口的cuDNN flag相关链评分路径差异，
但不据此确认唯一底层内核原因，不声称追溯证明R1实际启用了某个TF32内核。
这一限制保留，不另加无终点的定位运行。

## 当前决定

保留“仅为链GRU明确禁用TF32”的实现候选，不直接采用packed或改原R4。
本次没有梯度、AdamW、完整joint loss或吞吐准入结果。未来若继续，只做明确精度
合同下的一次完整步等价/速度验证，保留原逐链参照、原容差与原shared设置；
不能把全模型都改精度后相互接近当作保持原行为。

R4已因规定训练/推理时限未完成正式比较，见[r4终态](r4_terminal.json)。
本诊断没有修复或补齐R4，也不是论文方法收益。下一科学判断仍须完整C0/B1
DEV链P/R、实际自由候选覆盖、背景负担和成本；未获得这些证据前不宣称新模型有效。
