# 原 TRAIN 完整链首尾与内部几何诊断

2026-10-13。用已完成的 E、M、等权能量组合和硬共识预测，判断主参考链失配发生在转录首尾、内部 CDS 结构，还是没有同方向 span 覆盖。类别在新结果前固定；它们将帮助选择下一项 TRAIN 改进，不构成因果证明、独立验证或新的 DEV 方法。

## 固定输入与分母

只使用原两物种六条 TRAIN：AT NC_003070.9、NC_003071.7、NC_003076.8；rice NC_089043.1、NC_089044.1、NC_089046.1。参考为13449803保存的22,688条 CDS-assessable 主链 profile，不重读参考 GFF、DNA、势数组或权重。E fit13499817／评价13499823，M fit13510324／评价13510351，组合评价13551150／六份已完成 GFF，硬共识评价13550067。共识仍由两原输出的全部 CDS 几何交集重建，不根据参考筛选。

保留全部参考，包括 unsupported strata、参考重叠和长于 native window 的链。按既有 ruler 用 seqid、strand 和有序 half-open CDS 区间精确匹配；phase 不参与该 ruler，源预测已有独立 DNA／phase 验证。每染色体重放四组完整 exact-chain PRF，物种及 pooled TP／参考数也须与已完成评价一致。这里不重算 base、背景或其他 isoform，亦不重复源 DNA 验证。

## 参考相对类别

每个 arm 对每条主参考只给一个类别，优先级如下。首尾指转录方向中的最外 CDS 边界：正链 start 为最小 CDS 起点、stop 为最大 CDS 终点，负链顺序反转；不是转录本 UTR 边界。候选必须与参考同 seqid、同 strand 且 gene span 有正长度重叠，首尾接触不算重叠。

1. `exact_primary`：全部有序 CDS 区间相同。
2. `same_outer_boundaries_different_internal_CDS`：同一预测拥有两个正确外边界，但内部 CDS 区间不同。
3. `both_boundaries_in_different_predictions`：首尾分别存在于不同预测，没有一个预测同时拥有两个正确边界。
4. `transcriptional_start_only`：只有正确首端。
5. `transcriptional_stop_only`：只有正确末端。
6. `same_strand_span_overlap_neither_boundary`：有同方向 span 重叠，但两个外边界都不匹配。
7. `no_same_strand_span_overlap`：没有同方向 span 重叠。

源输出先核对同 strand span 非重叠，才使用排序起终点的 bisect 索引；参考不作非重叠限制。分类不能判断所有替代 isoform 或真实生物学错误：span 重叠可能仅落在 intron，首端失配也可能同时伴随内部结构失配。数量和比例按 retained references 解释，不使用这些类别生成 subset precision／F1 或标签驱动推理选择器。

## 输出与下一步

报告六染色体、两物种、pooled 七类数量，single／multi CDS 和九条 long-span 的同口径分层，以及 E→组合、M→组合类别转移与 exact gained／lost。结果用于确定下一项 TRAIN 机制实验应优先改进完整链首端、末端、内部结构还是漏覆盖；不会据此调整冻结 DEV、封存数据访问、比较者阈值或当前运行作业。

一个 private-teodoro-gpu／generanno／no-requeue／无 GPU 的2CPU、4GiB、10min作业，先执行四组轻量合成测试，再读取完成源和 profile／预测。测试针对负链边界、内部区间与 phase 的口径、分裂首尾、half-open 接触／异向预测／重叠参考和索引前提。独立输出为 `outputs/M28-FULL-GRAMMAR-E-M-GEOMETRY-TRAIN`，新 JSONL／summary 声明64MiB上限；不写 GFF、不训练、不做 forward，不改旧产物或模型。

完整 STATUS／JOBID／started／summary／reference_geometry.jsonl 和 root／batch／extern COMPLETED0:0 才构成完成。确定性接线错误允许同范围一次修复并用独立新输出保留首次失败；科学负结果保留，不取消作业或自动延长。短作业由主对话直接消费；只有需要外部等待时才复用原固定 watcher 与同一 heartbeat，不创建新的监控对话或自动化。

## 实际作业

五个新叶文件完成 AST／SBATCH 语法核验、精确同步回读后，提交 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-E-M-GEOMETRY-TRAIN.sbatch`，实际返回13555928。提交前 own expanded46，小于并发80；private 节点可用，无覆盖目标或10min维护冲突。使用 completed 源产物，不把已离开 controller 的旧源 ID 写成依赖。提交不是诊断完成，结论等候程序产物与全部 Slurm 终态。

## 完成结果

13555928于08:45:53–08:47:09UTC运行76秒，root／batch／extern全部COMPLETED0:0，MaxRSS1008524K。四组测试通过，程序27.106205秒；STATUS／JOBID／started／66759-byte summary和10625649-byte reference JSONL齐全，总输出小于64MiB。四组既有每染色体 exact-chain PRF重放一致，物种／pooled TP一致，保留全部22,688主参考、5,395条single CDS、17,293条multi CDS及九条long-span。主另行复算六chr→species→pooled七类数量、分层分区、两张转移矩阵的行列边际及 gained−lost，全部一致；没有重跑推理或改旧结果。主已直接消费，未加入 watcher、不补发监控交接。

| 参考相对类别 | E | M | 等权能量组合 | 硬共识 |
| --- | ---: | ---: | ---: | ---: |
| 完整主链精确命中 | 16420 | 16106 | 16527 | 15636 |
| 同一预测首尾正确，内部 CDS 不同 | 1740 | 1858 | 1632 | 1161 |
| 正确首尾分属不同预测 | 114 | 128 | 128 | 63 |
| 只有转录首端正确 | 956 | 943 | 887 | 603 |
| 只有转录末端正确 | 1501 | 1543 | 1482 | 1025 |
| 同向 span 重叠但首尾均不匹配 | 199 | 177 | 177 | 75 |
| 无同向 span 重叠 | 1758 | 1933 | 1855 | 4125 |

组合的6,161条未命中主参考中，无同向 span 重叠1,855条（30.108749%），首尾正确而内部 CDS 不同1,632条（26.489206%），只有末端正确、首端失配1,482条（24.054537%），只有首端正确、末端失配887条（14.397013%）。剩余128条是分裂首尾，177条是有同向 span 重叠但两端不匹配。百分比分母是未命中主参考，不是全部参考或预测。不能由这些描述性类别推断全部缺失源自 ATG、splice、语法限制或模型 score。

| 转移 | 原 exact TP | 保留 TP | 新增 TP | 丢失 TP | 净变化 |
| --- | ---: | ---: | ---: | ---: | ---: |
| E→组合 | 16420 | 16154 | 373 | 266 | +107 |
| M→组合 | 16106 | 15978 | 549 | 128 | +421 |

E→组合新增373条中181条原本是首尾正确／内部不同，96条只有末端正确，68条只有首端正确，24条原无同向覆盖；其余4条来自分裂首尾／双端均不匹配。M→组合新增549条中287条原本是首尾正确／内部不同。两者内部类别分别占新增的48.53%和52.28%，支持“内部结构是重要改进方向”，但不是独立机制证明。组合未保留所有 parent TP：E丢失266条、M丢失128条，需同时报告。

| 物种 | 主参考 | E TP | M TP | 组合 TP | 共识 TP | 组合净 TP 相对 E |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| AT | 17739 | 13111 | 12873 | 13080 | 12611 | −31 |
| Rice | 4949 | 3309 | 3233 | 3447 | 3025 | +138 |

组合的完整链 F1 两物种都略高于 E，但 AT TP反而少31条；此前的AT F1改善不能写成召回改善。硬共识虽然完整链 precision较高，其无同向 span 覆盖的主参考增至4,125条，相比两parent／组合均有明显覆盖代价。该数量只描述已输出链的覆盖，不宣称这些基因在模型内部没有证据或语法不能支持。

下一项机制设计不采用“只加首端权重”这一单因解释。先把本几何分解与正在运行的 R1 missed-reference score／grammar 诊断和4-cell学习率对照共同使用，区分共享语法不支持、能量竞争和实际首尾／内部失配；在原 TRAIN 范围选择有针对性的 objective／表示改进。当前结果不改变冻结 DEV family，不放行新物种、封存数据或独立验证，也不证明 Nature Communications 及以上论文目标已达到。
