# M26 R4：partial 特征层级与 CDS / span-only 分解

冻结于 2026-09-26，运行前。探索性测量诊断，不是独立面板复现。

## 问题与固定输入

NCBI GFF3 的 partial 是 feature 层级属性。gene/mRNA 的不完整不能直接等同 CDS 不完整。本次核对原始 gene、mRNA、CDS 属性，判断 R2 的 partial coding gene 区域是否包含可评价的完整 CDS；再分解实际 CDS 覆盖与仅跨度连接覆盖。

沿用 M26 两条 DEV 染色体、原始参考、五套保存预测及 R2 六分区。原 6,450 complete-primary、阈值、模型、预测不变。Setaria、其他封存数据不访问。无新数据下载、推理、训练或参数扫描。

## 预先规定的分析

1. 按 gene / mRNA / CDS 行分别计数 partial=true、range-only 与未标记；CDS 同一 Parent 的任意行标记均保留。报告属性组合，不能把只查看父级当成 CDS 完整性。
2. 对每个编码转录本报告 gene、transcript、CDS 的独立 partial 标记。CDS sequence compatibility 定义：无 overlapping CDS 区间、第一段 phase=0、长度被3整除、仅 ACGT、ATG 起始、末端 TAA/TAG/TGA、无内部同框 stop。它是技术兼容性，不是表达或功能证明。exception / transl_except 单列，不静默修正。
3. 对 gene 被标 partial 的集合，计算至少有一个 CDS-unflagged 且 sequence-compatible 转录本的基因数。每基因按最长 CDS、再 transcript ID 选择一个候选，只做另名的参考政策敏感性，绝不替换历史主参考。另计任意该类 isoform 的精确匹配。对所有被旧政策排除的编码基因也报告这种候选数，保留排除原因。
4. 对五工具按原六区输出 predicted CDS union、predicted span-only union；二者互斥且和等于 R2 span union。另在 partial coding gene 区分参考 CDS 区域与其余区域的负担。全按碱基并集，不累加重叠转录本。

## 判断与停止

若 partial gene 的负担主要来自 CDS 本身完整的转录本，先把中心问题收窄为 feature-level 可评价性与参考选择偏差，不能继续以“不可评价的不完整 CDS”解释所有 partial 区域。若主要是 span-only，则优先结构连接机制；若主要是真正 CDS-incomplete，才保留原 partial-uncertainty 假说。具体整数和比例全部报告，不用单一事后阈值择优。

原参考数量、R2 区域分子或 CDS / span 恒等式无法对账时停止科学解释；只允许一次不改变科学定义的工程修正，保留失败记录。不依据结果追加阈值扫描或新模型。

预算：一个 CPU 作业，private-teodoro-gpu，2 CPU、8 GiB、15分钟、0 GPU；预计数分钟。输出 outputs/M26-PARTIAL-SEMANTICS-R4/，代码和合成测试公开，原始大型输入不上传。

来源：[NCBI GFF3 官方说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)，partial / start_range / end_range 与 Start and stop codons 小节，2026-09-26读取。
