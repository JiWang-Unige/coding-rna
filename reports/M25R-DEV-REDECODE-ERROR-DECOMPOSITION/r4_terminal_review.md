# M25R Stage 1 R4：终态审阅与大小写勘误方案

日期：2026-09-08。状态：执行已完成；Pro 复审已完成；CPU 正式勘误尚未实施。M25R development NO-GO 保持不变。

## 执行与证据边界

- Baobab job `12270926`：COMPLETED，ExitCode `0:0`，elapsed `3-12:11:06`；2026-09-04 12:25:58 至 2026-09-08 00:37:04（Europe/Zurich），gpu035。
- 运行代码：`01fbbf15544f18ef0fca19688adc1aa0d3284248`。源产物：`outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4/`。
- 三个 immutable checkpoint、原 grid rows 601/2401/4626 均完成。每 epoch 的 diagnostic、reference_attrition、candidate_lineages、replayed_predictions、structural_validity 及最终 stage1_diagnostic 齐全。
- 原 1,536 sampled training windows；验证范围 Arabidopsis NC_003074.8 与 rice NC_089041.1，共 6,450 reference chains。resolved_inputs 记录 Setaria reads=false、weights_updated=false、threshold_or_decoder_search=false；本次审阅未读取任何 Setaria 文件。
- 四项冻结指标重放 absolute error 每 epoch 均为 0。每 epoch reference ledger 为 6,450 行，exact-chain 数与 aggregate 一致；lineage terminal accounting 完整。
- 这里只进行现有结果取证、报告与方案交付；没有修改模型、诊断代码或 R4 原始产物，没有提交新作业。

## 保留的模型结果

| 指标 | epoch 1 | epoch 2 | epoch 3 |
|---|---:|---:|---:|
| exact CDS interval F1 | 0.120414079 | 0.078302211 | 0.038261788 |
| exact CDS chain F1 | 0.324988301 | 0.237745740 | 0.140227624 |
| pooled intergenic FPR | 0.012468286 | 0.005077354 | 0.003415162 |
| gene-count ratio | 0.325271318 | 0.182945736 | 0.103410853 |
| predicted chains | 2098 | 1180 | 667 |
| exact chains | 1389 | 907 | 499 |
| missing chains | 5061 | 5543 | 5951 |
| false-positive chains | 709 | 273 | 168 |

这些是各 epoch 原冻结 tuple 的诊断，不是新的 checkpoint/tuple 选择。原 count gate 0.80–1.20 仍失败。epoch-1 Arabidopsis FPR 约 0.025668 是物种异质性证据，不应事后替换原 pooled FPR gate。没有 Setaria transfer 或独立测试结论。

## 已证实的诊断大小写错误

`src/screen_anchor/data.py::read_fasta()` 虽在 docstring 声称 uppercase，实际保留输入大小写。R4 `load_species()` 直接传递该结果。生产 decoder 与 trace decoder 在入口 `.upper()`，但 `canonical_truth()` 和独立 `structural_validity_components()` 直接与大写 ATG/stop/GT–AG 比较。独立 evaluator 自带的 uppercase FASTA reader 没有被 R4 使用。

Codex 在 Baobab 对三份原 GFF3 做了全量只读复核：保持 parser、坐标、codon features 与组件函数不变，仅将内存序列转大写。两条验证染色体分别含 4,969,681 与 10,055,400 个小写碱基。

| 检查 | epoch 1 | epoch 2 | epoch 3 |
|---|---:|---:|---:|
| 原序列：准确复现原有效数 | 1338/2098 | 745/1180 | 429/667 |
| 内存 uppercase：有效数 | 2098/2098 | 1180/1180 | 667/667 |
| uppercase 后所有组件失败数 | 0 | 0 | 0 |

因此撤回旧 validity fractions、start/stop/splice failure counts 以及 `SCIENTIFIC_NO_GO_INVALID_STRUCTURES` 的科学解释。原文件保留作为历史记录，不覆盖。100% 合法性仅表示 emitted models 满足被检查的结构规则，不等于准确率、足够 recall 或 caller feasibility。

`canonical_reference_count=2666`、canonical/noncanonical strata、motif reachability 1693/1758/1738 和 truth-assisted 897/577/320 不再用于科学解释。前几项可 CPU 重算；完整 truth-assisted 需要未持久化的 truth-position boundary/phase scores，不能从 candidate coordinates 猜测恢复。

## 保留的诊断及归因限制

- reference 主 stage assignment 不依赖 canonical 标签或 truth-assisted flag；其 trace 已做 uppercase。epoch 1/2/3 的 non_intergenic_block 为 1919/1673/1210，ordered_CDS_runs 为 586/618/800，phase_check 为 1145/1770/2595。没有单项达到 missing chains 的一半；保守结论为 mixed attrition，不是单一 backbone 失败。
- prefilter 2281/1235/700，末端 boundary filter 拒绝 183/55/33。因此仅放松最后一层 threshold 不能解释或解决数量缺口；不授权调参。
- epoch-1 single-CDS exact 1043/1499，multi-CDS 346/4951；epoch-3 为 387/1499 与 112/4951。328 条 span>6144 的参考链三个 epoch 均无 exact recovery。该关联不能单独证明 context length 是因果瓶颈。
- matched-gene strand accuracy 为 0.999379/1/1（分母 1610/1012/548）；exact-matched CDS phase accuracy 均 1（分母 2181/1367/652）。这是成功匹配子集的条件性指标，不是全基因准确率。
- train/validation CDS-region F1：e1 0.9300/0.9120，e2 0.9439/0.9313，e3 0.9526/0.9244。区域信号较好没有转化为完整链恢复。
- train/validation CDS-base phase accuracy：e1 0.3378/0.3368，e2 0.3416/0.3475，e3 0.3416/0.3438；truth-CDS-start phase accuracy：0.6279/0.6271、0.5348/0.5137、0.3189/0.3465。candidate-start expected phase accuracy 为 0.7358/0.5279/0.3146。训练与验证均存在 phase-fit 问题，不能简单解释成只在验证集过拟合，也不能据此隔离 targets、loss、head、LoRA 或 representation 的责任。
- validation boundary AUCPR（start/stop/donor/acceptor）：e1 0.2398/0.1836/0.3594/0.3591；e2 0.2676/0.1645/0.3656/0.3855；e3 0.2543/0.2498/0.5361/0.4449。冻结 threshold 下 exact event recall 约 0.94–0.99，但 precision 很低（e1 start 0.00184、donor 0.00574；e3 start 0.04152、donor 0.08571），不能把高 event recall 当作边界定位已解决。
- raw-head region FPR 与最终 decoded FPR 是不同口径，不能互换。

## Pro 复审与独立核实边界

审阅对话：https://chatgpt.com/g/g-p-6a29d586630481918525796032225f68-ji-wangke-ti/c/6a8c18a3-4370-83eb-b9cf-8924409e2f12

2026-09-08 请求的回复已完成（页面显示思考 6m19s）。Pro 报告 ACCESS_VERIFIED，读取精确提交 01fbbf；明确未读取 Baobab outputs，真实计数引用 Codex 的只读取证。它报告函数级 harness 6 passed/0.03s 与 py_compile exit0，不是完整仓库或真实产物测试，不能作为集成 PASS。

Pro 确认大小写根因、保留 coordinate metrics/主 stage counts、撤回 invalidity、保持 M25R NO-GO，并建议唯一下一步是独立 CPU 勘误。Codex 已直接核实主 assignment 源码未以 canonical/truth-assisted 决定 stage；同意这一影响边界。Pro 的 GO 是技术建议，不代替用户实施/作业授权。

## 唯一建议下一步：CPU 勘误（尚未实施）

1. 仅在诊断 `load_species()` 入口将读入序列 `.upper()`；不改公共训练 reader、模型输入路径、decoder、motif 字典、radius、threshold、weights 或 reference policy。
2. 增加一个窄 `posthoc_r4_case_correction.py`，只读现有 R4 GFF3、lineage/attrition/diagnostic 与 A/rice development FASTA/reference/config。禁止加载 checkpoint 或模型 forward。
3. 独立输出 `outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4-CASE-CORRECTION/`；保存 corrected validity ledger、reference canonical/motif correction 和 summary。原 R4 STATUS/JSON/GFF3/ledgers 不变。
4. 使用 block_span/runs 重建 I/C/G discrete path，先验证 transition totals 3727/3993/4017 完全复现，再重算 uppercase motif reachability。核对每条 exact reference 被 corrected motif 覆盖，以及 exact<=motif<=transition。不能把 Pro 提出的可重建性当作已经验证。
5. 重算 canonical-compatible reference 集合及分层；明确它包含 start/stop、GT–AG、完整 ORF 条件，而不只是 intron motif。
6. corrected truth-assisted total 写 `NOT_RECOVERABLE_WITH_SAVED_ARTIFACTS`，不填猜测值，不为该项重跑 84 小时 GPU。transition/motif 仅解释为冻结 state path/grammar 下的条件性 recall ceiling，不是 representation 或 F1 上界。
7. 聚焦验证：mixed-case diagnostic loading、合法 toy 在归一化后的 validity、同一 logits 大小写 decoder invariance、state reconstruction/每条 exact-chain reachability。运行现有 focused suite 和语法检查；真实 correction 中 counts/原坐标指标/主 stages 必须保持不变。

停止条件：原预测/metrics/stages 改变，重建 transition 不一致，exact chain 不被 motif 覆盖，count/orientation 无法对账，或需要越界访问/加载 checkpoint。仅 truth-assisted 不可恢复不算失败，明确 N/A 即可。成功后写 `POSTHOC_CORRECTION_COMPLETED_REVIEW_REQUIRED` 并停止。

本次交付止于终态审阅与上述具体方案。CPU 勘误实施与运行仍待明确批准；baseline、SegmentNT、新训练、额外 seed、GPU re-decode 和所有 Setaria 操作均未获本次授权。
