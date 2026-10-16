# M26 R5：固定序列与预测的 RefSeq—Araport 来源/表示对照

2026-09-26，下载与分析前固定。由R4实测提出；不替换历史协议。不属于C独立复现：原RefSeq注明TAIR and Araport，两者共享注释来源，且Araport文件有后续更新。

## 输入与范围

- 仅拟南芥NC_003074.8，长度23,459,830；Araport对应Chr3。仅在TAIR10官方FASTA Chr3与现有FASTA逐碱基大写字符串完全相等后作坐标比较，不用同长度或assembly名字替代。
- 保持M26五套预测不变；无新推理、训练、阈值变化、Setaria或其他sealed数据访问。
- 原RefSeq：data/m1_screen/arabidopsis_thaliana/reference.gff3，GCF_000001735.4，原4,151 complete-primary。
- Araport11：TAIR作者公开存档 https://zenodo.org/records/15889110 ，文件Araport11_GFF3_genes_transposons.20240701.gff.gz，16,470,613字节，CC-BY-4.0。存档声明发布2024-06-30，文件日期2024-07-01，记录创建2025-07-14；不称当前最新版本。
- 官方TAIR10 FASTA：https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
- 输入下载到data/m26_reference_pair_20260926/，保留URL/HTTP头；总新增存储限1GiB。不将原始GFF/FASTA上传Git。

## 解析与政策

Araport基因类型来自locus_type，显式映射protein_coding；其他gene/非编码/TE标签不猜为蛋白编码。将Chr3映射NC_003074.8，其余染色体不评分。保留原属性，转换只新增统一gene_biotype字段。精确链比较含strand及全部CDS坐标，末端stop不凭参考匹配结果加减。

每来源各给两种预先固定政策：
1. parent-filter：既有gene/mRNA partial排除、最长CDS再ID选择；RefSeq必须复现4,151和原M26该物种链指标。
2. CDS-assessable：CDS自身所有行无partial/range且通过R4 sequence compatibility；不因gene/mRNA partial直接丢弃；每基因最长CDS再ID选择。ORF不兼容的旧参考移除数量与新增数量分列，不把两种效应混淆。

报告各来源/政策的参考数量、精确链TP/P/R/F1，以及与该政策exon-span背景相对的预测CDS-span分子分母。该背景敏感性不是生物学FPR真值；另外报告固定所有protein-coding gene背景的CDS和span覆盖分子/分母。不挑更有利的政策，不把两来源并集当真值。

## 配对机制检验

在两来源都存在的相同locus、相同strand+CDS坐标链中，报告CDS-assessable及parent-filter eligibility的对应。基因键为RefSeq locus_tag与Araport gene ID，不能只用易变的mRNA ID。一个链可有多个同坐标isoform记录，eligibility按存在至少一个合格记录定义。

按各工具计数因父级筛选而丢失的共享、可评价精确链；共同结构与来源独有结构分别统计。此受控子集只报告对象计数/召回，不能因限制参考而硬造总体precision或F1。

## 停止与资源

FASTA不一致、seqid/坐标/Parent关系不成立、RefSeq原指标不能复现时停止解释。来源具有共同构建依赖已经明确，不宣称独立性。如果效应主要在共享链的父级过滤中且语义对齐后消失，则收束为测量修正，关闭“普遍生物学新机制”叙事；若残余来源争议存在，没有独立证据不得判谁正确。

一个CPU作业：private-teodoro-gpu，2 CPU、8GiB、15分钟、0 GPU。输出outputs/M26-REFERENCE-PAIR-R5。预估几分钟；允许一次不改变科学定义的工程修复，保留失败。后续新推理只在对象、版本、masking和预算另行明确后考虑，不是本协议内容。
