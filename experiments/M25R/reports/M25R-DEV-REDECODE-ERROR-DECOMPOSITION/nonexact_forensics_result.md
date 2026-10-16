# Phase-bypass 新增9条 nonexact：参考结构诊断

2026-09-08；CPU Slurm 12520779，private-teodoro-gpu，无GPU，1 CPU/4 GiB/5分钟上限；COMPLETED 0:0，耗时6秒。只读既有配对输出及拟南芥/水稻参考GFF；无新推理、训练、阈值调整或Setaria访问。

## 结论

9条不是9个已证实的生物学假阳性：4条精确匹配非主转录本CDS链，2条精确匹配因gene-level partial标记被排除的注释CDS链，2条相对主参考截短，1条与partial参考差一个6 bp剪接边界。这里的精确仅指同链完整CDS坐标链相等，不宣称所有注释属性或生物学真实性已验证。

冻结的 complete-primary 评价定义不变：9条仍是该定义下的 nonexact。本诊断不是事后重定义成功门槛，也不把partial参考或替代转录本添加回冻结主指标。

## 逐条证据

| 物种/链/lineage | 解释 | 参考证据 |
|---|---|---|
| 拟南芥 +342 | 已注释非主转录本CDS精确匹配 | NM_202484.2；主参考NM_001202873.1，AT3G03790 |
| 拟南芥 +4633 | 已注释非主转录本CDS精确匹配 | NM_114434.2；主参考NM_001339214.1，AT3G45650；替代链多一个78 bp内含子 |
| 拟南芥 -6740 | 已注释非主转录本CDS精确匹配 | NM_001337453.1；主参考NM_111168.2，AT3G02980 |
| 水稻 -3380 | 已注释非主转录本CDS精确匹配 | XM_015791979.3、XM_015791980.3；主参考XM_015791978.3，LOC4343229 |
| 拟南芥 +5800 | 已注释CDS精确匹配，但gene标记partial | NM_115416.2，AT3G55590；gene与mRNA均partial=true |
| 拟南芥 -232 | 已注释CDS精确匹配，但gene标记partial | NM_001340090.1、NM_001340091.1、NM_116002.4，AT3G61380；即使NM_001340090.1自身未标partial，gene仍partial=true |
| 拟南芥 +2975 | 相对主参考提前终止，4 CDS对7 CDS | NM_113485.5，AT3G25830；起点相同，终点提前1270 genomic bp；共享3内含子，缺少后3内含子；新增CDS碱基100、漏掉715 |
| 水稻 +2217 | 相对主参考晚起始，9 CDS对13 CDS | NM_001422076.1，LOC4342850；起点后移1672 genomic bp，终点相同；共享8内含子，缺少前4内含子；漏掉768 CDS bp，无新增CDS碱基 |
| 拟南芥 -772 | partial参考上的6 bp剪接边界偏移 | NM_115511.1，AT3G56530/NAC064；基因和mRNA均partial=true；预测中间CDS为[20949383,20949658)，参考为[20949389,20949658)，其余两CDS一致 |

最后三条没有与任何同链protein-coding注释转录本CDS链完全相等的匹配；这不足以证明生物学错误。+2975与+2217是明确的主参考链截短；-772位于已注释partial基因，并非无注释基因组区域。

## 验证与产物

调用既有eval_structure_diagnostic.py的parse_annotation和primary_transcripts，重建6450条complete-primary参考。将反向互补坐标转换为0-based half-open基因组坐标；所有面板B输出均与保存的B GFF链匹配，64条已知exact候选全部匹配原参考ID，避免负链坐标或参考选择错配。

脚本：scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/nonexact_forensics.py。

产物：outputs/M25R-E1-PHASE-NONEXACT-FORENSICS/{cases.json,summary.json,STATUS}；日志logs/M25R9FORENSIC_12520779.{out,err}。

机器输出分类是相对于complete-primary的第一层描述：annotated_nonprimary_exact=6、intron_chain_difference=2、no_same_strand_primary_CDS_overlap=1。这里nonprimary意为不在选定complete-primary集合中，不能将6条全部解释为替代isoform。本文进一步直接检查原GFF的partial标记，细分为4+2+2+1；no_primary_overlap亦不能解释为无注释。原机器产物保留不覆盖。

## 对研究决策的影响

phase一致性硬拒绝既阻断了64条面板主参考精确链，也阻断了至少6条与其他已注释CDS链一致的候选；因此不能以本面板“9条nonexact”否决phase-bypass路线。另一方面，仍存在链截短，且所有9条通过高置信度边界分数与ORF检查，说明这些检查不能单独保证主参考结构精确。

下一项真正决定是否继续的实验应是固定epoch-1、原row601阈值、完整拟南芥/水稻开发集的phase-only A/B配对对照，不是新增head训练。A必须先重现R4冻结结果；B只绕过phase拒绝，其余候选生成、ORF、阈值、去重及评价保持一致。主要决策继续使用原complete-primary指标和原开发门槛；全注释链匹配/partial分类只能作解释性补充，不用于选择阈值。失败或A不一致即停止。

本次仅完成CPU诊断，未启动完整开发集GPU对照；其窗口量、可复用分数覆盖、资源预算及停止条件需先由现有脚本/产物核定，再提交具体授权方案。局部门控结果不能证明完整流水线净收益、head因果缺陷、Setaria泛化或整个路线终态。当前开发NO-GO、Setaria embargo不变。
