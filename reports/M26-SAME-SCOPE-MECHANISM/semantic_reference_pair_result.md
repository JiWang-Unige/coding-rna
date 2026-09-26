# R7结果：CDS语义对齐后无残余结构争议，结束DEV来源对照线

2026-09-26。Slurm13227054，private-teodoro-gpu / gpu034，2CPU/8GiB/0GPU，COMPLETED 0:0，51秒，MaxRSS1,017,736KiB；7项测试通过（12.70秒），stderr0字节。首次13227041的失败及唯一工程修复见semantic_reference_pair_repair.md。此前R5/R6状态原样保留。

## 决定性结果

两来源同一FASTA Chr3/NC_003074.8逐碱基一致，均5,460个protein-coding locus，locus全部对应。CDS-assessable政策下，**7,970条逐locus唯一CDS isoform链完全一致，真正来源独有的可评价CDS链为0，残余结构造成的五工具exact判定翻转均为0。**

5,460个位点分解：

| 情况 | locus数 | 含义 |
|---|---:|---|
| 相同可评价primary CDS | 5,428 | 结构与primary一致 |
| 仅primary选择不同 | 9 | 双方被选链均存在于两来源可评价isoform集合；不是新增/缺失结构 |
| 两来源均无可评价CDS | 23 | 保留不可评价，不强配结构差异 |
| 真正primary结构争议 / 仅一方可评价 | 0 | 本次没有可用于后续“来源争议”验证的对象 |

原RefSeq4,151参考与五工具该染色体chain指标、span FPR全部精确重放。RefSeq新政策增加1,290个gene、移除4个无兼容CDS的旧gene，得到5,437个primary；Araport由5,460移除23个后同为5,437。各来源内部primary另有变化的数量与跨来源9个位点不是同一统计。

## 评分：仅拟南芥这条染色体，不是旧pooled6,450总表

| 工具 | RefSeq parent | Araport parent | RefSeq CDS-assessable | Araport CDS-assessable |
|---|---:|---:|---:|---:|
| M25R_A | 0.348375 | 0.339028 | 0.340170 | 0.340170 |
| M25R_B | 0.554558 | 0.557840 | 0.557925 | 0.559155 |
| ANNEVO | 0.762515 | 0.762502 | 0.762845 | 0.764199 |
| Helixer | 0.645664 | 0.651481 | 0.651944 | 0.652872 |
| Tiberius | 0.740011 | 0.741478 | 0.741771 | 0.743115 |

以上为exact CDS chain F1。语义对齐后工具排序不变：ANNEVO > Tiberius > Helixer > B > A。两来源的小幅分数差异全部由9个位点的primary选择解释；在这些位点B匹配primary为RefSeq0 / Araport6，ANNEVO1 / 8，Helixer1 / 6，Tiberius1 / 8，A均0。任意兼容isoform的精确匹配不存在来源差异。预测数量原样保留，未因TE/pseudo或非coding区域而删预测。

B的policy-span FPR：RefSeq parent12.7439% -> RefSeq CDS-assessable2.7512%；Araport CDS-assessable2.7294%。后两者的少量差异还受所选转录本exon-span背景影响，不是CDS结构真伪证明。全部coding gene背景分母两来源一致，10,842,973bp。该解释性指标不替换历史冻结FPR，也不撤销NO-GO。

8,000条同locus同CDS的共享链中，7,970条两边均可评价；其中1,920条因父级过滤仅在Araport parent政策保留，涉及1,290个locus。B精确匹配其中401个locus。这复现的是既有结构被筛选排除，不是401个新基因。

## 输入与验证边界

- 113个原gene缺失类型以同发布表补齐；扩展gene-like后另7个非coding父节点也由同表补齐，总120。coding参考仍5,460，不因接入877个TE gene-like及203个pseudogene而扩大。
- 879个旧“missing gene”mRNA全部连接到TE gene-like父节点；普通gene的mRNA为9,474。180,753条Chr3原feature行全部保留。
- 50个miRNA表述差异及1个pseudogene/pre-tRNA分歧保留为非coding告警；coding身份不确定为0。新coding_view_qualified=true不回写旧qualified=false。
- 两来源coding CDS链的phase按转录方向均一致。sequence compatibility仍不证明表达、功能或完整转录本。
- 所有大输入/派生GFF/逐位点账留在远程outputs/M26-REFERENCE-PAIR-R7-REPAIR1/；公开紧凑JSON、代码、测试和协议。
- 来源共享TAIR/Araport构建依赖，本结果不是独立复现，也不能外推“所有参考注释都没有结构争议”。基线是旧缓存，预训练暴露与历史版本局限仍在。

## 决策

关闭本DEV来源敏感性作为新论文主张的路线，不再找第三来源，不为不存在的残余结构争议启动100位点验证。数据质量Skill推动了任务层级对齐，但所得主要是本地测量修正，不包装成普遍生物学机制。

保留原R1 NO-GO与R3在原6,450范围下的候选支持分解；本次不把R3上界外推到新参考。Nature Communications目标未完成。下一步应回到有独立可证伪机制问题的研究设计；本阶段结果与拟议方向交由下一轮Pro商议，不以调阈值、增加seed或更大backbone替代新主张。
