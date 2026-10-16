# A11：完成池组件均值未恢复新的正确长链

2026-10-03。零更新TRAIN探针工程完成、实际预测／原27参考核验通过；事前机制gate FAIL，关闭本读出的DEV扩展。RAW_SUM、COMPONENT_SUM、COMPONENT_MEAN都命中15/27，没有新增≥5段参考。不能把改变候选长度分布当成提高生物学准确率。[事前协议](a11_component_mean_protocol.md)

## 固定执行与数值证据

CPU checks13359870五项PASS后，真实probe13359876仅提交一次：private／1 RTX3090／2CPU／16GiB／20min上限，root/batch/extern均COMPLETED／0:0／1m41s。STATUS、JOBID匹配；8预测及27参考完整。原A3末step004608、8个原A5 TRAIN压力窗、Q16／beam4／完成4／cap128不变；更新0、新头参数0、无新GEN forward／checkpoint／DEV／test。预测循环83.0921s、总84.2850s、peak GPU allocated148045312 bytes是探针记录，不是完整annotation吞吐。

主对话CPU核验13359929 COMPLETED／0:0／14s，证据 `outputs/M28-A6-A7-A11-RESULT-CHECK/summary.json`。直接从全部8真实预测及原A5共同greedy池重建四臂集合，逐一核对27原参考身份、链坐标、实际命中、完成评分公式及DNA合法性；weighted两臂实际计数相同，RAW_SUM每query链集合与原A5完全一致，共同raw分数最大差0。逐物种／exon类计数另从全部27行重算一致。没有按已知参考挑新窗口、移除难例或改变候选。

| 范围 | 参考/window pairs | ORIGINAL | RAW_SUM | COMPONENT_SUM | COMPONENT_MEAN |
| --- | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | 17 | 7 | 7 | 7 | 7 |
| 水稻 | 10 | 8 | 8 | 8 | 8 |
| 1段 | 9 | 8 | 8 | 8 | 8 |
| 2–4段 | 9 | 4 | 4 | 4 | 4 |
| ≥5段 | 9 | 3 | 3 | 3 | 3 |
| 全部 | 27 | 15 | 15 | 15 | 15 |

MEAN对COMPONENT_SUM新增0、丢失0、新增≥5段0。逐物种不低于ORIGINAL、原27身份保留及raw重放条件通过，但严格pooled增加与新增长链两个条件失败；不改条件或用非下降冒充PASS。

## 均值确实改变了完成池，但没有正确长链收益

| 实际读出 | 完成query-chain records | 1段records | ≥5段records | 最大段数 | greedy联合后每窗unique总计 |
| --- | ---: | ---: | ---: | ---: | ---: |
| RAW_SUM | 499 | 322 | 69 | 23 | 104 |
| COMPONENT_SUM | 499 | 323 | 69 | 23 | 105 |
| COMPONENT_MEAN | 499 | 155 | 196 | 62 | 137 |

499是实际终止候选records，不是512满额假定，也不是499个独立基因；同链可由多个query生成。MEAN确实更常保留长链，不能称实现没起作用，但它仍没有命中原12条漏检参考。完成均值与weighted sum共享真实live搜索，不能补回该搜索已剪掉的前缀，也不能从本压力窗结果否定任意其他生成／长度方法。

ORIGINAL／RAW联合池按R4参考相对label1/0/-1为15/88/1；COMPONENT_SUM为15/89/1，MEAN为15/121/1。MEAN增加的是当前主参考不一致候选，不是已证明生物学假基因；-1未知仍保留。这里没有最终GFF、precision／F1或背景测量，不能把15/137当成annotation precision，也不是新的独立泛化证据。

继续边界：不扫EOS系数、均值指数、beam宽度、阈值，不追加训练或进入DEV。原冻结A3、A9 FAIL、FULL与A10均保留；本结果支持结束“只改完成均值”路线，科学目标和其他独立机制分支继续。重要新模型方向需基于真实评分／选择或训练覆盖／表示证据，而不是把相同压力窗反复调到通过。
