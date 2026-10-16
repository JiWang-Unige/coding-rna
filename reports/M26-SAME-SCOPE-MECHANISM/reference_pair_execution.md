# R5 执行记录：输入不合格，尚无来源比较结果

## 已核实输入

2026-09-26下载到data/m26_reference_pair_20260926/，约51MiB压缩文件。URL见reference_pair_protocol.md；HTTP头保留。TAIR作者Zenodo记录15889110的Araport11.20240701.gff.gz为16,470,613字节，许可CC-BY-4.0；原始文件不改、不上传Git。

TAIR官方FASTA Chr3与本项目NC_003074.8逐碱基大写字符串完全相等，长度23,459,830。两个作业都通过此断言后才进入GFF规范化。不是仅核对长度/名称，也未计算无决策作用的新哈希。

## 两次运行

- 13226530：private-teodoro-gpu，2CPU/8GiB/0GPU，FAILED 1:0，27秒，MaxRSS636,164KiB；2项合成测试通过(14.23秒)。在GFF文本读取阶段遇到非UTF8的0x91。定位为原文件第20,350行curator_summary的引号字节，不涉及坐标、ID或Parent。
- 按预设的一次工程修复机会，对此明确输入用Latin-1无损字节映射读取，规范副本写UTF8；不用忽略/删除整条gene的方式跳过。非ASCII描述不用于评分，原始gzip保存。没有修改科学政策。
- 13226538：唯一重试，独立输出outputs/M26-REFERENCE-PAIR-R5-REPAIR1，FAILED 1:0，11秒，MaxRSS636,124KiB；2项测试通过(3.71秒)。在规范化gene类型时发现locus_type字段缺失，停止于评分前。

原失败日志与不完整规范化文件均保留。首次目录中额外的post_repair_snapshot文件是修复后复制的代码/sbatch，不冒充首次提交快照；初次差异为gzip文本读取未指定Latin-1。

## 为什么不能直接补一个默认值

对Chr3原始gene行的定向核对：6,544行中113行缺规范locus_type属性；101行在其他属性值中含protein_coding类型文字、3行含其他类型文字、9行无这种类型文字。该核对不做任何预测评分，不把文字替代规范字段后继续运行。

实际例子：AT3G01560的类型文字出现在Dbxref值里的“locus_type equals protein_coding”，并不是规范locus_type属性。AT3G06010和AT3G60830的gene行仅剩描述等字段，没有该类型。另有small_nuclear_rna例子，不能将所有缺失类型的gene一律补为protein_coding。描述中还出现未转义的分号，不宜通过放宽解析静默改变科学对象。

本阶段未生成reference_pair.json，因此没有来源×政策F1、工具排序或共享链政策差异的实测结果。成功的FASTA核验与测试不等于来源比较完成。

按冻结的执行边界，此R5不再追加重试。后续只能先处理输入资格：选择可追溯、格式一致的注释发布，或使用同发布的权威gene-type映射；明确记录其与原文件的依赖关系，重新固定输入与规则后再评价。不能根据预测成绩选择来源，也不能把来自同一Araport构建的镜像叫独立注释。

本记录不改变R4结果、R1 NO-GO、Setaria封存及无新训练的边界。
