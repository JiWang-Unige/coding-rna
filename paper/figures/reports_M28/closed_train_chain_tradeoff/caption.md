# 已闭合 TRAIN 整链选择计数

A36、A37 均减少了相对 primary 参考的已知非匹配，但同时减少正确链和长链配对，因此都只通过原五条件中的两项。这张图用于方法开发中的负结果与权衡说明，不作为完整基因注释准确率或独立泛化结果；A38 未纳入。

## 图注

**Selected chain counts in closed TRAIN experiments.** A31_A3, A36 and A37 are compared on the same 16,321 natural candidates from 1,491 original TRAIN windows, using the original gain-zero selection rule. Panels A and B show selected reference-relative matches (1), known nonmatches (0) and protected unknowns (−1) for Arabidopsis thaliana and Oryza sativa. Panel C shows pooled selected matching chains with at least five exons. The original eligible reference-window denominators are 2,172 and 744, with 1,086 pooled at-least-five-exon reference pairs; they are not replaced by the number of candidates retained by a model. A36 and A37 reduce known nonmatches but fail both species' strict matching-pair increases and the pooled long-chain increase. Counts refer to window–chain pairs rather than unique genes. Unknown is not an error class, and known nonmatch is not a biological false-positive determination. Training inputs and capacity differ: A36 includes positive injection and 366,191 trainable parameters, while A37 is natural-only with 612,324. This is a descriptive comparison, not a causal ablation or genome-level F1 evaluation.

| 方法 | pooled 正确配对 | 已知非匹配 | 未知 | ≥5-exon 正确配对 |
| --- | ---: | ---: | ---: | ---: |
| A31_A3 | 2009 | 598 | 447 | 700 |
| A36 | 1342 | 213 | 230 | 611 |
| A37 | 1827 | 239 | 209 | 639 |

PNG 与 SVG 均使用零起点；两物种 composition 图使用一致尺度，未知类以灰色点纹单列。已检查实际 PNG 的标签、图例和条件线，没有裁切或标签碰撞。SVG 保留矢量与文字，适合后续论文排版；当前仍是负结果图件，不是论文主结果。

## 论文结果段草稿

In an exploratory TRAIN screen, fewer reference-relative nonmatches did not imply improved chain selection. At the fixed zero-gain operating point, the local structured emission model (A36) reduced pooled selected known nonmatches from 598 to 213, but also reduced matching window–chain pairs from 2,009 to 1,342. The jointly trained context and ordered-chain model (A37) selected 239 known nonmatches and 1,827 matching pairs. Relative to the same-pool baseline, matching pairs decreased in both Arabidopsis thaliana (1,531 to 1,138 and 1,396, respectively) and Oryza sativa (478 to 204 and 431). Selected matching chains with at least five exons also decreased, from 700 to 611 and 639. Both models therefore met only the two nonmatch constraints and failed all three matching-pair improvement criteria.

Paired analysis localized most losses of previously selected matching candidates to rejection at the chain-versus-empty operating point. Of 714 matching pairs lost by A36, 665 (93.1%) had nonpositive new gain; the remaining 49 had positive gain but were not selected. For A37, 232 of 279 lost matching pairs (83.2%) had nonpositive gain, while 47 remained positive but unselected. Gains of 47 and 97 matching pairs, respectively, did not offset these losses. These rejected candidates were already present in the natural pool, so their loss was not caused by absence from that pool. The partition identifies a scoring-stage location for further study; it does not establish why their scores became nonpositive.

All counts refer to candidate pairs in 1,491 original TRAIN windows, not unique genes or chromosome-level annotation F1. The common natural readout pool contained 16,321 candidates, and the original reference-window denominators remained 2,916 overall and 1,086 for the at-least-five-exon subset. Unknown candidates were retained as a separate protected category, and known nonmatch was defined relative to the primary reference rather than as a biological false-positive determination. Training supervision and capacity differed between these models, so the descriptive comparisons do not isolate a causal benefit or failure of structured prediction, joint context learning or the null component. They motivate a separate test of absolute chain-versus-empty constraints, not a claim of generalization or superior genome annotation.

该段用于后续论文的探索性方法开发或补充结果，不代替正文所需的完整 raw-FASTA 与独立泛化证据。计数和转移依据分别见 [A36 闭合结果](../../a36_structured_emission_result.md)、[A37 闭合结果](../../a37_context_chain_result.md)；A38 结果仍不纳入此段。

## 来源和复现

数据来自已完成的 A31 readout13398521、A36 CPU13408451 和 A37 CPU13410742。对应原始 summary 快照在 sources/；provenance.json 保存逐物种三臂计数、原始分母和重算五条件。作图前已核对三份完整1491窗结果、共同基线、同池三类数量、选择数加总和原 gate，未重跑原模型或神经前向。

绘图脚本为 scripts/experiments/M28-METHOD-RESTART/plot_closed_train_chain_tradeoff.py。在 generanno 的 Matplotlib3.10.6 环境中，以 --source-dir 指向本目录 sources/，以 --output-dir 指向新的输出目录；已有图件不会被覆盖。在 Baobab 上通过 Slurm 执行。原执行脚本为 sbatch/M28-CLOSED-TRAIN-CHAIN-FIGURE.sbatch。

实际作业13413469的 root／batch／extern 均 COMPLETED0:0，6秒；private0GPU／1CPU／2GiB／3分钟，gpu034，batch MaxRSS147080K，STATUS／JOBID相符，Requeue0／Restarts0。PNG537285字节、SVG103352字节、provenance4679字节，已取回并完成视觉核验。该作业只读三份已闭合汇总，没有 checkpoint、特征tensor、optimizer、DEV、Setaria、A33 或 A38 产物访问；不修改原冻结结果。
