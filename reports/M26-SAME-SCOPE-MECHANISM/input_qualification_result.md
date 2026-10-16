# R6输入资格结果：缺失类型已补齐，整体资格未通过

2026-09-26。最初13226882在过窄ID正则处FAILED 1:0（3秒，MaxRSS28,084KiB，3测试通过），源表显式uORF不应被当作不合法标准gene model。第五轮Pro商议后另定TYPEVIEW合同，13226976 COMPLETED 0:0（18秒，MaxRSS137,732KiB，4测试通过）。实际运行的是2CPU/2GiB/0GPU，private-teodoro-gpu；尚无预测评分。

## 确定结果

- 原Chr3 GFF中6,544行gene；同发布功能表包含11,789个标准模型、7,624个locus，另有19条明确uORF子特征。19条宿主均由原文Derives_from确认；类型不传播到宿主。
- 113个缺失规范locus_type的gene全部取得唯一类型：110 protein_coding、3 small_nuclear_rna。没有利用描述中的“locus_type equals”文字或模型预测决定类型。
- 6,381个gene的已有GFF类型与功能表相容；另50处差异全部为GFF的mirna与功能表的miRNA_primary_transcript。它们都是非蛋白编码对象层级的表示差异，保留不默改。
- 全部解析出的gene中protein_coding计数5,460。功能表无其他未知格式对象。
- 检查器另报879条mRNA的gene父节点缺失。定向核对AT3G60930：父节点实际以transposable_element_gene显式存在，且有两个mRNA子节点。这是当前检查器只识别feature=gene的限制；单例不能证明879条全部如此，仍需全量父节点连接后判断是否影响编码比较。

原始qualification.json因此qualified=false，未生成合格GFF。这是保守检查器对其支持范围的判断，不应称整个TAIR源文件损坏。没有坐标错误、CDS phase错误或CDS缺mRNA的问题进入本次issues列表；后续仍需把父节点类型完整连接，不能只凭这一列表宣布全部正确。

## 下一决定

本阶段已经回答113个缺失类型能否同发布补齐：能。剩余问题不需要逐基因查询，也不需要第三个注释发布。下一合同应明确所有gene-like父节点（gene、pseudogene、transposable_element_gene）及对象类型分层，验证全部879个报告项；针对protein-coding比较单列类型和结构资格。

若任何coding链的父子关系/类型仍不明，停止受影响比较并报不确定范围；若仅非coding表示差异，不能把它伪装成coding不确定性，也不能为了强行通过把TE或miRNA改成coding。保留全体原始feature行和坐标，单独输出逐节点连接账。之后才是两来源×两政策评分及真正结构差异分解。

本阶段没有读取预测来选择数据、没有工具成绩，未改R1 NO-GO与历史6,450，未访问Setaria。紧凑JSON和测试公开；原始GFF、功能表、全部逐gene账留在远程。
