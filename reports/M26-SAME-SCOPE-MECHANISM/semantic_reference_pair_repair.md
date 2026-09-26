# R7一次工程修复记录

2026-09-26。13227041 FAILED 1:0，38秒，MaxRSS637,148KiB，6项测试通过（18.41秒）。FASTA逐碱基一致已通过；879个R6父节点问题全部连接为transposable_element_gene（879个mRNA对应877个gene-like节点）。没有预测评分结果，失败停在输入资格阶段。

唯一残余项AT3G20365：原GFF feature=pseudogene、locus_type=pseudogene，子特征为pseudogenic_tRNA（Chr3:7,104,309–7,104,392，正链）；同发布功能表AT3G20365.1明确为pre_trna。两者对是否为假基因存在生物学注释分歧，但均明确非protein_coding。原检查器将所有类型差异都计为protein_coding_type_uncertainty，过度泛化了当前endpoint的不确定性。

本合同已规定TE/pseudo单列且不进入coding正例；修复只将这一明确类型对作为非coding告警保留，不宣称同义，不删除原行，不增加coding对象，不改变背景或任何预测分母。所有涉及protein_coding或未知标签的冲突继续阻断。新增对应真实类型对及coding冲突的测试。

按R7合同允许的一次工程修复，独立输出outputs/M26-REFERENCE-PAIR-R7-REPAIR1。最初实际提交代码/sbatch已经在原作业目录保存；原资格JSON中的false和1个不确定项保留，不能回写为true。修复后会重新生成coding_view_qualified及非coding类型争议计数。除此之外来源、坐标、CDS政策、评分和统计规则均不变；不追加第二次修复。
