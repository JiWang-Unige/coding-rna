# M26 R2：参考背景差异归属（R1 后提出，执行前固定）

R1 显示了巨大的 FPR 政策敏感性和工具排序反转。下一次有限 CPU 运行仅解释它，不修改模型或重判旧门槛。输入仍为相同两条开发染色体、五套已保存预测与当前参考。

按以下优先级把每个碱基分到唯一类别：

1. 原 complete-primary exon span。
2. 有 complete-primary 的编码 gene 的其余 gene span。
3. gene 层面标注 partial 的 protein-coding gene。
4. 其余无 complete-primary 的 protein-coding gene（例如只有 partial transcript）。
5. 非编码或其他已注释 gene。
6. 全部 gene 注释之外。

统计每类基因组碱基数、各工具预测 CDS-span 并集覆盖碱基数。每工具、每物种必须精确复现 R1 两种 FPR，才解释区域归属。补充“全部 protein-coding gene 背景”口径，非编码注释仍留作潜在负区域；不默认把所有非编码 gene 都排除后宣称更准确。

如果旧 FP 负担主要落在 partial 编码区域，后续优先研究可评价区域与负区域定义，而非直接发起 repeats/short-ORF 原因实验；否则按实测类别调整假说。无论结果怎样，覆盖已注释区域不证明预测 CDS 链正确，gene-level partial 标签也不证明具体 ORF 真伪。

一次 Slurm CPU：private-teodoro-gpu，0 GPU、2 CPU、8 GiB、10 分钟上限，预计数分钟。输出 `outputs/M26-REFERENCE-BACKGROUND-R2`。一个针对重叠区域优先级、并集不重复计数及分母的合成测试；不重跑 R1、不访问 Setaria、不下载新数据，不扩大到因果或生物学确认。
