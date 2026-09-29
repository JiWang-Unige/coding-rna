# M27写作包：受限摘要、三图图注与主张矩阵

日期2026-09-29。证据冻结于6350e2bc07d043ba06d41a7b75000502da76b9cb。本文件是写作交付，不是完整论文、已绘制图件、投稿或录用声明。除既有产物的阅读与整理，0新实验、0新前向、0训练、0Setaria访问。

## 定位与标题

**建议定位：定量复现／受限实现研究。** 不以“发现早停绕行”作为原创概念，不以未成功的单边检验建立神经/约束的完整因果分解，不声称选择性诊断。现有资料可以形成可审阅的复现稿底稿；独立贡献及具体期刊适配仍需要作者判断，Nature Communications级主张尚未成立。

工作标题：

**Native gene-annotation responses to paired synonymous and premature-stop substitutions: a registered chromosome-scale pilot**

其中registered指公开冻结协议和代码先于相关实验，并非声称在外部临床/研究注册平台注册。图注中的panel是拟定排版；本轮没有生成图像或补充数据。

## 受限英文摘要

Gene-annotation pipelines may alter predicted structures in response to premature-stop substitutions, but native pipeline responses must be distinguished from biological outcomes and component-level explanations. We conducted a prospectively specified computational pilot using reference-selected sites on human chromosome 22. Twenty-one loci met fixed sequence and annotation criteria; conditional analyses used the loci uniquely recovered in the unmodified genome by each caller: 12 for ANNEVO and 10 for Tiberius, with nine shared loci. Each synonymous or premature-stop single-nucleotide substitution was evaluated independently in the complete chromosome context using the caller's native workflow. All 44 planned mutant runs completed. Every synonymous condition retained the target coding-sequence chain. Each caller produced two strictly assigned premature-stop-specific bypass responses and two unassignable responses. Retaining unknown outcomes in the fixed denominators gave paired-difference identification bounds of [0.1667, 0.3333] and [0.2000, 0.4000], respectively; these are not confidence intervals or population-rate estimates. The callers had no shared confirmed bypass locus within the common set. Four preregistered native reruns reproduced the selected target structures and whole-chromosome chain multisets. In one ANNEVO case, same-DNA cache replay was exact, but releasing a single registered coding-continuation transition did not alter the bypass output or selected-path score components. These results document reproducible, version- and panel-specific structural responses while limiting a particular local explanation. They do not establish biological splice changes, neural-only causation, a general effect of integrity constraints, or a selective risk diagnostic.

## Figure 1. Paired native responses in a fixed conditional panel

**A, Reference selection and conditional execution.** Twenty-one loci were selected from GRCh38 chromosome 22 and GENCODE v49 before mutant prediction. Qualification required a single distinct, evaluable protein-coding CDS chain per gene under the frozen policy; it did not establish uniqueness of all RNA isoforms. The registered internal codon and replacement rules generated independent synonymous and TAA premature-stop alleles at the same genomic position. Each mutant genome differed from the common unmodified chromosome at exactly one registered base; mutations were not accumulated across loci. ANNEVO uniquely recovered 12 target chains in the unmodified input and Tiberius recovered 10; nine loci were shared. WT-inexact loci remain visible as CONDITIONAL_RESPONSE_NOT_RUN rather than being counted as mutant negatives.

**B, All-locus response matrix.** Display all 21 registered rows, in the original genomic-coordinate order, with each caller's WT eligibility and synonymous/PTC category. All synonymous runs in the conditional sets retained the original chain. Among PTC runs, ANNEVO produced 2 BYPASS, 1 PTC_TERMINATION, 1 OUTPUT_SPLIT, 6 OUTPUT_ABSENT and 2 UNASSIGNABLE responses (n=12); Tiberius produced 2, 1, 5, 0 and 2, respectively (n=10). OUTPUT_ABSENT denotes absence of related output after a complete computational run, not biological gene loss. UNASSIGNABLE remains unknown and is not merged with either category.

**C, Fixed-denominator paired differences.** Plot the PTC-minus-synonymous strict-bypass identification intervals: ANNEVO [1/6,1/3] and Tiberius [1/5,2/5]; in the common nine-locus set, both intervals are [1/9,1/3]. These are bounds induced by ambiguous assignment, not sampling confidence intervals or statistical significance claims. Distinct 12- and 10-locus conditional populations cannot rank caller quality. The common set has no locus where both callers have a confirmed bypass: the ANNEVO positive is POTEH and the Tiberius positive is CACNG2. Synonymous and PTC substitutions are different nucleotide replacements, so this design does not isolate a pure biological stop-codon effect.

Data: reports/M27-ALLELIC-INTEGRITY/paired_input_manifest.json, wt_result.json, paired_result.json. The complete 21-row ledger and frozen classifier, not selected illustrations, define the denominators.

## Figure 2. Strict bypass structures and computational reproducibility

**A, Four assigned computational structures.** Compare each WT CDS chain with its independently predicted PTC chain. All diagrams use zero-based half-open coordinates, explicitly show strand and the edited codon, and preserve the complete terminal CDS anchors of the uniquely assigned chain. ANNEVO/POTEH skips the original 174-nt CDS [15695644,15695818), creating a spanning intron [15695485,15698661). ANNEVO/CRYBB2 introduces the 27-nt intron [25227912,25227939). Tiberius/CRYBA4 introduces the 15-nt intron [26625574,26625589), and Tiberius/CACNG2 introduces the 33-nt intron [36566378,36566411) on the reverse strand. Each assigned chain excludes all three edited codon bases through altered splice connectivity and meets the frozen phase/ORF and foreign-coding safeguards. The diagrams are predictions, not observed RNA molecules or validated splice products; the 15-nt GC–AG sequence is an allowed computational motif, not evidence of biological splicing.

**B, Preselected native confirmation.** According to the selection rule fixed before mutant inference, confirm the first positive locus in genomic order for each caller: ANNEVO/POTEH and Tiberius/CRYBA4. Repeat both synonymous and PTC conditions, giving four native runs. Each same-allele repeat exactly reproduces the complete target chain, strand, phase and competing-candidate identity, and the whole-chromosome CDS-chain multiset: 530 chains for each ANNEVO run and 575 for each Tiberius run. Equality refers to an old-versus-new run of the same allele, not equality between synonymous and PTC structures. Repeated computations add no independent locus or biological replicate to Figure 1.

**C, Native workflow context.** For both Tiberius confirmation alleles, the whole-chain filtering sequence reproduces 706→585→575, with unchanged target assignment at each step. Each run makes four observed neural calls, including initial and boundary-reprediction requests; four was observed, not forced by the replay policy. These records locate the selected bypass upstream of the two recorded final filters but do not separate neural predictions, HMM processing and adaptive request generation. Native score recording leaves returned arrays unchanged and is not itself a score intervention.

Data: paired_result.json and confirmation_result.json; detailed selected structures and comparisons are retained in these compact reports. Full native GFF/GTF, request traces and score arrays remain on Baobab.

## Figure 3. A preregistered single-edge release does not reverse POTEH bypass

**A, Exact replay and the local intervention.** Both arms use the same complete 01_PTC chromosome and its native ANNEVO score cache from the confirmation run. The replay arm must reproduce all 530 whole-chromosome CDS chains and the target's original 174-nt exon-skipping assignment before the intervention executes. In the intervention arm, at plus-strand zero-based position15695751 (the third base of the registered TAA), a copy of the A-conditioned transition matrix replaces only CDS1_TA→CDS2 from negative infinity with the ordinary continuation score0. Only this genomic position references the modified matrix; no DNA, emissions, other transition or normalization is changed. The actual target decoding region is [15689950,15720000), with local edited-base index5801.

**B, Native control flow and selected output.** Both arms use the original float32 core and a single min_intron_length=1 pass; neither triggers the native short-intron rerun at20. Candidate regions and recorded target emission arrays match. The released edge is not selected. Both final outputs retain the same bypass chain, and all530 chromosome-wide chains are unchanged. No non-target CDS-chain differences were detected chromosome-wide; outside the target span, no hidden-state differences were detected within the recorded decoding region. The preregistered endpoint—restoration of the original coordinate/phase chain while retaining the PTC—is not met; this is a scientific negative, not a failed job.

**C, Scores of the two actually selected paths.** Sum the native floor-transformed, class-mapped log-emissions over the entire fixed target span [15690077,15719777), including introns: E=−1253.5999851226807, transition sum=−10 and total=−1263.5999851226807 in both arms. Over the complete decoding region, excluding native unscored t=0, the corresponding values are E=−1276.0396003723145, transition sum=−10 and total=−1286.0396003723145. Both emission differences are0, below the prespecified absolute tolerance0.001nats. The two selected paths also have equal totals when rescored on the common released graph. **This does not compare the original reference-chain score against the bypass score:** the original chain was not recovered, and its score was not established by this experiment. Equal scores do not establish unique optimal paths.

The result limits only this single-position, single-edge intervention. It neither proves that all ORF constraints are irrelevant nor attributes bypass to the neural model alone. No alternate edges, combination rescue, second-caller intervention or additional locus was run after the negative endpoint. The complete job used0GPUs and no new neural forward calls.

Data: single_edge_result.json and single_edge_contract.md; code edge_decode.py and edge_run.py. The prospective code/protocol commit is dfed2c440b6af13d3730afeab5acac6384519758. This selected-case intervention is exploratory attribution with a frozen execution rule, not a confirmatory biological experiment.

## 逐主张证据／边界矩阵

| 拟写主张 | 地位与证据 | 允许表述 | 必须保留的限制 |
|---|---|---|---|
| 固定条件面板存在PTC特异严格绕行 | 主要经验结果；44/44原生运行，paired_result.json全21行 | 两工具各有2个明确配对阳性，未知保留后识别下界仍正 | 条件于各自WT精确；不等于发生率/统计显著/同位点普遍复现 |
| 同义对照均保留WT链 | 主要经验结果；12/12与10/10 | 在这些固定替换/目标中成立 | 不代表所有syn都无影响，不是纯生物学stop效应 |
| 预定案例可计算重现 | 有效性控制；confirmation_result.json | 四次同allele重跑的目标与全chr一致 | 不增样本，不证明生物学正确或完整因果归因 |
| 单边释放未逆转POTEH | 局部探索性归因；exact replay gate及single_edge_result.json | 在固定PTC分数下，该干预不足以改变输出 | 非所有边/所有ORF约束检验，不证明原链得分低或NN单独致因 |
| 短内含子和exon skip是真实转录产物 | 不支持；只有计算输出 | 只能称predicted结构 | 无独立RNA、长读长或蛋白裁决 |
| 已得到可区分不可靠修复/合理纠错的风险信号 | 不支持；candidate=false，无双向验证 | 只能列未实现目标 | 不以分数提高/结构合法/基准匹配代替诊断价值 |
| 新通用方法、clean zero-shot或NC级贡献 | 不支持 | 不进入标题/摘要结论 | 训练暴露未知；新颖性受先行工作限制；无独立泛化/决策收益 |

M26的A→B/R3是独立实现资产，参考政策问题属于测量勘误，R8属于有限阴性；不作为本稿的机制证据。完整跨阶段映射见 [证据映射](manuscript_evidence_map_20260929.md)。最接近文献与核查时间见 [新颖性边界](../reports/M27-ALLELIC-INTEGRITY/novelty_scope.md)；本写作包未开展新的文献检索，不声称穷尽相关研究。

## 一次性后续判断

本写作包完成了Pro第14轮建议的下一交付范围，不消耗剩余GPU预算。第15轮Pro已完成一次性写作审阅，只要求收紧隐藏状态比较的地域范围；已修正，不再循环审阅同一骨架，也不把“还不够NC”自动转化为再补实验。

若接受复现型定位，可继续完整稿、真实图件与复现产物索引；若要求必须具有尚未成立的机制/诊断/生物学推进，则本稿按技术/复现报告归档，当前实验线保持关闭。新科学问题需独立立项，不能从本次关闭自动授权。作者名单、投稿目标、声明和正式提交均未决定；本文件不作代作者承诺。

## 第15轮审阅与交付记录

Pro报告完整读取77743384033212e4f0c1b396f8d2b23c2be3016d的本文件及证据映射（2分43秒），未访问Baobab原始缓存、未独立重跑。其唯一实质修订建议是区分全chr CDS链比较与仅目标解码region内的隐藏状态比较；Figure3B已按此收紧。其余数值、条件分母、未知界限和阴性解释未发现改变科学解释的错误。此记录不是外部同行评议，也不代表作者接受投稿定位。

若作者选择复现型完整稿，限定下一制作范围为：正文/方法；由冻结JSON生成三图及源数据表；逐图绑定提交、字段、脚本/命令和公开/远程保存/可重新生成状态的复现索引。当前尚未执行这三个完整稿制作交付，不虚报完成。正式署名、机构声明、投稿选择和提交均未进行。
