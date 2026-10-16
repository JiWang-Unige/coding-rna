# 最终4608步严格raw-FASTA原TRAIN端到端结果

2026-10-09，Asia/Shanghai。最终权重的生产raw-FASTA接入及原六条TRAIN染色体端到端推理／共同尺评价已完成，数值可核验；完整CDS链表现仍有明显缺口，尤其水稻。原TRAIN执行阶段没有DEV／封存数据执行、训练参数变更、增训或对外发布；后续已授权DEV执行另见文末，原论文及原生基线竞争goal未完成。

## 实际执行与工程结论

唯一权重为 `outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt`，三遍／4608步／seed0／epsilon0.1的固定最后checkpoint。严格raw调用只接收原TRAIN完整FASTA，原生冻结GENERanno编码及C0 head后双strand whole-contig exact exon-budget skyline解码；不读reference、R4标签、特征缓存或观测模型，不运行optimizer。窗口、nearest-center所有15角色ownership、cap128、unit role、null0及FP64均保持合同不变。详见[full_grammar_raw_train_contract.md](full_grammar_raw_train_contract.md)。

数组13443323_0–5与CPU诊断13443329的root/batch/extern均为COMPLETED/0:0；六份固定范围输入、双strand25512窗口、完成summary/GFF及诊断产物齐全。输入共156779632bp，21032条unique完整预测链，重复链0。GFF数量、实际TRAIN DNA/codon/phase/cap、坐标和同strand gene-set非重叠检查已通过。诊断先确认全部推理结束及固定checkpoint／caller身份，再打开原六条TRAIN reference；预测没有使用reference。

接入合成检查13443304与评价组合合成检查13443443成功，分别覆盖生产skyline绑定／跨窗完整链／实际codon carry与GFF phase，以及strand敏感exact key／其他isoform另计／背景gene跨度／计数micro汇总。它们不是研究准确率数据，不冒充独立科学评阅。主对话另外用已有summary的TP/预测/参考计数重算P/R/F1及背景覆盖，并核对逐染色体→逐物种→pooled计数相等；无评价公式或聚合不一致，不重复重跑未变验证。

共同尺保持M26 CDS-assessable：每gene在完整、可翻译CDS候选中取最长、平分取transcript ID；不依M28训练eligible、motif、exon cap、长度或同strand overlap缩小参考。完整链key含seqid、strand及所有半开CDS边界，phase另核验；其他assessable isoform命中不计primary TP。不评价UTR、全部isoform或非编码gene。

## 固定主指标

以下P/R/F1均为百分比；pooled由计数合并，而非平均物种F1。

| 范围 | 完整链 TP／预测／参考 | 完整链 P | 完整链 R | 完整链 F1 | CDS-base F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 拟南芥三条TRAIN | 8603／16341／17739 | 52.6467 | 48.4977 | 50.4871 | 93.0848 |
| 水稻三条TRAIN | 1363／4691／4949 | 29.0556 | 27.5409 | 28.2780 | 80.8511 |
| pooled六条TRAIN | 9966／21032／22688 | 47.3849 | 43.9263 | 45.5901 | 90.3601 |

species macro完整链F1=39.3825%，worst-species=水稻28.2780%。不分strand CDS-base union：拟南芥TP/预测/参考=19188173/19588682/21638607、P/R=97.9554%/88.6756%；水稻4775658/5063484/6749983、94.3157%/70.7507%；pooled23963831/24652166/28388590、97.2078%/84.4136%。较高base覆盖不能替代完整链边界、拼接和gene-set重建；二者差距提示完整链重建是当前需要解释的问题，不单凭此归因于某个模块、监督通道或训练量。

| TRAIN seqid | 完整链 TP／预测／参考 | 完整链 F1 % | CDS-base F1 % |
| --- | ---: | ---: | ---: |
| NC_003070.9 | 3454／6645／7135 | 50.1306 | 92.9198 |
| NC_003071.7 | 1995／3835／4293 | 49.0896 | 92.6350 |
| NC_003076.8 | 3154／5861／6311 | 51.8239 | 93.5577 |
| NC_089043.1 | 463／1575／1638 | 28.8204 | 83.0429 |
| NC_089044.1 | 481／1545／1601 | 30.5785 | 81.5135 |
| NC_089046.1 | 419／1571／1710 | 25.5410 | 78.1510 |

背景是所有annotated gene feature跨度任一strand之外的bp；这里的span_FPR只是背景bp被预测链跨度覆盖的比例，不是gene-level分类FPR，也不认证新预测为生物学假阳性。拟南芥背景27874501bp，预测span/CDS覆盖111198/29948bp，span比例0.398924%，纯背景链15（0.538126/Mb）；水稻背景53345031bp，390606/105156bp，0.732226%，纯背景链147（2.755646/Mb）。pooled背景81219532bp，501804/135104bp，0.617837%，纯背景链162（1.994594/Mb）。命中非primary assessable isoform而未命中primary的预测分别718／30、合计748条，未加进上述9966个primary TP。

## 原TRAIN输出分层诊断

CPU-only作业13449803在唯一新root `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE/`成功，root/batch/extern均COMPLETED/0:0；UTC 2026-10-09 02:39:26–02:41:06、Slurm wall100秒、内部profile28.834秒、最大RSS1548228KiB。只重读已完成GFF、原六条TRAIN reference／已物化FASTA及原draw的TRAIN R4元数据，不读checkpoint／feature、不跑native forward、不改预测或训练。没有重跑主评价；新分层exact计数逐染色体／物种／pooled与13443329完全一致。脚本内可手算fixture验证了不同contig／strand、单多CDS段、长span和positive metadata集合划分。

| 范围 | 单CDS段 TP／预测／参考 | 单CDS段 F1 % | 多CDS段 TP／预测／参考 | 多CDS段 F1 % |
| --- | ---: | ---: | ---: | ---: |
| 拟南芥 | 2766／3430／4114 | 73.3298 | 5837／12911／13625 | 43.9931 |
| 水稻 | 673／937／1281 | 60.6853 | 690／3754／3668 | 18.5934 |
| pooled | 3439／4367／5395 | 70.4569 | 6527／16665／17293 | 38.4416 |

水稻多CDS段P/R为18.3804%/18.8113%，单CDS段为71.8250%/52.5371%；拟南芥多CDS段为45.2095%/42.8404%。差距在结构复杂组更显著，但该描述性分层并不能证明某个head、边界、剪接或目标函数是原因，也不授权改语法、重训或筛除难例。

原1536 draws（每物种768、唯一核心窗745/746）中的R4 `positive_ids`按**全局CDS几何key**与固定primary参考相交，拟南芥2093、水稻720条，均属于当前primary参考。这里只检查metadata中的完整链出现，不证明链内所有fine-role观测均完整无冲突；训练实际仍可能有alternative/overlap遮蔽、coarse／partial观测及邻域context／冻结encoder曝光。

| 范围 | positive_ids曾出现：TP／参考 | 完整链 recall % | 未出现：TP／参考 | 完整链 recall % |
| --- | ---: | ---: | ---: | ---: |
| 拟南芥 | 1081／2093 | 51.6484 | 7522／15646 | 48.0762 |
| 水稻 | 215／720 | 29.8611 | 1148／4229 | 27.1459 |
| pooled | 1296／2813 | 46.0718 | 8670／19875 | 43.6226 |

曾出现组中的多CDS段召回率为拟南芥756/1599=47.2795%、水稻107/508=21.0630%。因此缺口不只在metadata未出现的链：不能仅归因于原draw未覆盖完整TRAIN。不过这不是随机曝光对照，未出现组不是“未训练”、held-out或独立泛化，出现/未出现差异也不是训练效果的因果估计；这里不为reference曝光组制造没有合理预测分母的precision/F1。

跨度大于24576bp的固定primary参考为拟南芥2、水稻7、合计9条，本次TP=0且所选预测中也无超过24576bp的完整链。这不说明解码器必然截断：whole-contig语法允许跨窗，此前人工30000bp链接入检查只证明工程支持。本次真实9条全部未完整命中说明当前权重尚无这些长链的成功证据，不能据人工可达性宣称learned长基因准确率；样本仅9条，不泛化成所有长基因的概率结论。

分层 `summary.json`及逐reference的 `reference_chain_profile.jsonl`（8536829bytes）留存，所有参考仍在原主分母，没有按分层重算或替换主指标。此诊断定位了可报告的结构缺口，但没有选择新参数、模型、阈值或新训练方案。

## 原TRAIN几何误差追加诊断

2026-10-09另提交CPU-only作业13449985，执行 `profile_full_grammar_original_train_errors.py`；private-teodoro-gpu、无GPU、2CPU／4GiB／5分钟，新root `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-ERROR-PROFILE/`。它在人工fixture阶段因新代码使用两项CDS、项目实际key要求三项start/end/phase而FAILED/1:0，UTC03:31:39–03:32:02、wall23秒；尚未读取研究产物，不是模型或科学指标失败。失败目录和日志保留，最小修复使用项目实际三项格式并将新增CDS-segment事件明确投影为start/end几何、不混入phase，重跑只用独立root `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-ERROR-PROFILE-R1/`，不覆盖失败或冻结结果。

修复重跑13450002已成功：UTC03:34:47–03:35:08、root/batch/extern COMPLETED/0:0，实际wall21秒、内部8.546秒、MaxRSS1033992KiB，实际AllocTRES只含2CPU／4GiB、无GPU。R1 root有正确JOBID、COMPLETED状态、27246bytes summary及22688行逐参考误差清单。该诊断只复用已完成profile的固定primary参考链清单与六份冻结预测GFF，不再读DNA、reference.gff3、checkpoint或feature，不执行模型或optimizer。内部人工fixture覆盖strand方向端点、内部结构、intron中span重叠但无CDS重叠、两个预测对一个参考的覆盖，以及局部junction命中不等于整链命中；所有新分层的exact-chain TP／预测／参考与原13443329逐染色体、物种、pooled一致，主指标没有变化。

参考链按同contig／strand的正CDS-overlap选最大重叠bp预测；平分依次取更多exact intron、较小CDS端点距离、几何key／ID。所有primary参考仍在原分母，exact hit及计数与原完成诊断逐层核对；该关联不是one-to-one gene匹配，不把它作为新主评价。类别区分无同strand CDS覆盖、同intron链但start／stop／两端不同、两端相同而内部结构不同，以及端点与内部结构都不同。独立报告strand敏感unique CDS start／stop、CDS段及intron-pair事件读数；它们不是gene-level准确率或phase指标。多预测／多参考overlap只提示可能的分裂／融合形态，不证明其生物学含义或模块原因。

| 相对固定primary参考的几何类别 | 拟南芥 | 水稻 | pooled |
| --- | ---: | ---: | ---: |
| 完整链exact | 8603 | 1363 | 9966 |
| 无同strand CDS重叠 | 2597 | 1141 | 3738 |
| 同intron链仅CDS start不同 | 399 | 138 | 537 |
| 同intron链仅CDS stop不同 | 0 | 0 | 0 |
| 同intron链两端均不同 | 8 | 13 | 21 |
| 两端相同但内部CDS结构不同 | 2854 | 763 | 3617 |
| 端点与内部CDS结构均不同 | 3278 | 1531 | 4809 |

12722条未完整命中的primary参考中，8426条（66.2317%）涉及内部结构差异，558条（4.3861%）仅端点不同，3738条（29.3822%）没有同strand CDS重叠；三类计数精确覆盖全部未命中参考。内部结构差异占各物种未命中的比例为拟南芥67.1191%、水稻63.9710%。在多CDS段未命中的10766条中，其占比为8339/10766=77.4568%。因此现有输出缺口不是主要由“同一完整内部链只选错start／stop”这一形态解释；它仍不能证明是head、监督、局部context或global decoding中的哪一项导致。原主评价中的748个非primary assessable isoform命中继续单列；本几何差异也不等于全部为已确认的生物学错误。

| strand敏感unique几何事件F1 % | 拟南芥 | 水稻 | pooled |
| --- | ---: | ---: | ---: |
| CDS start位置 | 79.4107 | 65.2282 | 76.2833 |
| CDS stop位置 | 78.5196 | 58.4647 | 74.0970 |
| 完整CDS段两边界 | 77.7122 | 58.2060 | 73.6413 |
| intron两边界pair | 80.5294 | 61.4471 | 76.6010 |

这些事件全部带seqid和strand，CDS段／intron-pair必须两边界exact；不设置bp容差、也不以片段命中替代整链。共享位置或段在事件集合中去重，因此它们的参考数与gene参考数不同。局部intron-pair F1明显高于完整链F1，但水稻局部结构也仍有明显缺口；不能据此断言编码器已经正确、只剩跨exon组装。CDS start／stop是CDS端点，不是UTR端点、TSS或转录终止位点。

参考逐链bp加权的最大单预测覆盖率为拟南芥86.1019%、水稻66.7997%、pooled81.5200%；同strand全部重叠预测union覆盖为88.5224%、70.7480%、84.3031%。这里保留了重叠参考各自的CDS分母，不与原unstranded union-base F1混用。1355条参考与多个所选预测有正CDS重叠，54条所选预测与多个primary参考重叠；这不是已确认分裂／融合事件计数。无同strand覆盖而存在opposite-strand CDS覆盖的参考仅109／2、合计111条，未将反strand覆盖视为正确gene。该诊断是原TRAIN输出的描述性证据，不是DEV性能、独立泛化或可用于调参的确认性对照。

## 实际运行成本与监控

UTC 2026-10-08 18:13:15六条推理同时开始，各Slurm wall为38:14、24:29、33:54、35:02、32:09、36:18，最迟18:51:29结束。随后CPU评价18:51:31–18:52:29完成，Slurm wall58秒、评价内部43.583秒。整条推理→评价wall为39分14秒，不能与六条单GPU时长之和混淆，也不能把监控回传等待算为计算时间。

每条使用1×3090、2CPU/16GiB；六条Slurm单GPU占用合计12006秒（约3.335 GPU小时，非GPU活跃利用率）。raw主循环合计11799.096秒，其中原生encoding8650.297秒、joint decoding3143.877秒；各task外层入口合计11915.626秒，未包含Slurm包装开销。最大batch RSS7969924KiB（约7.60GiB）；各条PyTorch GPU allocated/reserved峰值2411252736/2554331136bytes（约2.25/2.38GiB），不等于全进程VRAM。六份raw结果的artifact_bytes_before_summary合计29459209bytes，不含物化FASTA或全部诊断文件。未与不同染色体／硬件／计时边界的历史基线作速度比较。

原固定监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3/local`于2026-10-09T02:21:39Z完成原TRAIN的一次实际终态回传，主体仍为gpt-5.6-luna/low；同一 `coding-rna-r5-light-monitor` heartbeat在该次交接后PAUSED并保留，后续DEV已复用启用，无新监控对话或重复自动化。原TRAIN首次定时运行及完整自动回传均已实际观察，不把配置成功等同于回传成功；当前DEV的验证另见其执行记录。监控只核对调度／完成结构，本报告由主对话核对实际指标。

## 证据位置与后续边界

原始产物位于Baobab项目 `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW/<seqid>/`、`outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC/summary.json`及 `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE/`。项目 `/home/users/j/jwang/coding-rna`解析到 `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna`；summary中的resolved checkpoint路径是同一文件，不是更换权重。代码／SBATCH保留于 `scripts/experiments/M28-METHOD-RESTART/` 和 `sbatch/`，未覆盖冻结预测或首训结果。

首训仅原1536 draws／1491唯一核心窗、每遍同一顺序，不是完整TRAIN所有窗口的训练；本次在其原六条染色体整条诊断，不等于独立测试。核心窗之外还可能经邻域context／coarse观测／冻结encoder有曝光，不能把未列入positive_ids的区域称为“未训练”。composite有限context观测目标也不因whole-contig推理而变成chromosome-joint MLE。冻结backbone的上游暴露未知，raw输入不证明clean zero-shot。

既有summary字段 `retained_unsupported_reference_strata`中的 `longer_than_native_window`仅是描述性跨度标记，不是当前whole-contig语法的自动unsupported条件；其他真实motif／终端结构／overlap限制仍保留在共同分母。现有主指标没有预设新PASS阈值，不将数值较差转换成未事先定义的科学gate，也不据此增训或选择其他权重。

下一阶段的[冻结DEV／基线方案](full_grammar_dev_baseline_plan.md)在后续人类明确授权后已启动：两条推理13449922并行，共同评价13449926等待整个数组afterok。仍固定最后权重、原生固定基线、同输入和共同尺；不重跑基线、不增训、不放行sealed/held-out/Setaria。本报告的TRAIN主结果不因此成为DEV结果，也没有在TRAIN与历史DEV指标之间作排名或宣称优于ANNEVO/Tiberius/Helixer、泛化成立、论文成功。原goal继续，科学比较与独立面板证据尚未完成。
