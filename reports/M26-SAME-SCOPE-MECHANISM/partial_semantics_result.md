# M26 R4：高 partial 负担主要不是不完整 CDS

2026-09-26。Slurm 13226466，private-teodoro-gpu，gpu034，2 CPU / 8 GiB / 0 GPU；COMPLETED 0:0，60秒，MaxRSS 1,225,616 KiB。3项合成测试通过，stderr 0字节。代码 partial_semantics.py；原始紧凑结果 partial_semantics.json，逐转录本账位于远程 outputs/M26-PARTIAL-SEMANTICS-R4/transcript_semantics.tsv。

## 会改变解释的发现

| 固定开发染色体 | protein-coding genes | gene partial | 至少一条 CDS 无 partial 且 sequence-compatible | 原主参考 |
|---|---:|---:|---:|---:|
| 拟南芥 NC_003074.8 | 5,460 | 1,309 | 1,290 / 1,309 (98.55%) | 4,151 |
| 水稻 NC_089041.1 | 2,307 | 8 | 0 / 8 | 2,299 |

拟南芥 gene-partial 集合中有 2,141 条 CDS-unflagged、sequence-compatible 转录本，其中1,524条 mRNA也标partial、617条mRNA未标partial。另有24条CDS标partial、2条CDS未标但sequence不兼容。水稻8条均为gene/mRNA/CDS都标partial且sequence不兼容。

两份参考所有被标记的gene/mRNA行均有 partial=true 与range属性，不是range-only误解析。CDS检查涵盖同一Parent的全部行，不把父级标记自动传播给CDS。sequence compatibility是完整ORF和基本坐标兼容性，**不证明表达、功能或全长转录本完整**；exception / transl_except在转录本账中保留。

因此，R2“旧FP-span的83.19%在partial coding gene区域”仍是正确区域计数，但不能被解释为83.19%落在不可评价的不完整CDS。当前评价器将gene/mRNA层级的partial排除用于CDS链任务，再将排除区域当负背景；它对两个来源的作用显著不同。原拟南芥头部注明TAIR and Araport，水稻为NCBI RefSeq GCF_034140825.1-RS_2024_06；不能把差异完全归因为物种或模型。

## 五工具的同一参考选择后果

以下均为拟南芥，新增候选基因分母1,290；不加入历史主表：

| 工具 | 精确匹配新增primary CDS链 | 精确匹配任意兼容isoform的排除基因 |
|---|---:|---:|
| A | 196 | 201 |
| B | 373 | 401 |
| ANNEVO | 492 | 581 |
| Helixer | 449 | 523 |
| Tiberius | 485 | 594 |

新增primary按最长CDS、再转录本ID选取，不能将任意isoform计数混入primary指标。水稻新增集合为空。以上是参考已有CDS结构匹配，不是新发现基因，也不是独立生物学证据。B在401个被旧政策排除的基因上精确匹配兼容isoform，这些匹配不应统称已证实假基因；但不能据此把全部2,837条nonexact解释为参考问题。

## CDS 与跨度不是同一种负担

B拟南芥在partial gene互斥分区的1,409,588bp：

- 真正预测CDS：948,924bp；
- 仅预测span连接：460,664bp；
- 其中934,569 / 948,924 = 98.49%的预测CDS，与上述兼容参考CDS的无链向碱基并集重合；
- 1,385,585 / 1,409,588 = 98.30%的预测span，与上述兼容参考CDS跨度并集重合。

无链向覆盖交集不能替代strand-aware精确链判定。五工具的CDS与span-only并集之和在全部10个工具×物种单元精确复现R2，原参考计数不变。

在所有已注释gene以外，B拟南芥107,311bp预测跨度中85,111bp (79.31%)仅为span-only；水稻170,789bp中53,525bp (31.34%)为span-only。这是连接跨度的诊断，不直接证明融合、内含子错误或生物学假阳性。

## 决策

1. 撤下“partial区域=不完整CDS无法评分”的笼统解释；优先按任务层级定义CDS可评价性。
2. 不撤销原固定R1 NO-GO，不改旧6,450分母，不将R3条件容量上界外推到新增参考集合。
3. 该结果首先是一项本项目测量修正，不自动具备论文新颖性，也不构成独立物种复现。
4. 下一项受控检验：相同FASTA/预测下，对照原RefSeq与公开Araport注释，区分相同CDS链的元数据政策差异与真正结构更新。两者有共同来源，只能称来源/表示对照，不称独立真值。独立面板阶段后移，不因找不到缓存无限检索或自动开新模型。

依据：[NCBI GFF3 feature-level partial说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)。此次数据质量Skill使核对从父级标签转向gene→mRNA→CDS对象层级；Baobab Skill保证分析在Slurm中执行。
