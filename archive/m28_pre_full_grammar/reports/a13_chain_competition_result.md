# A13 完整链竞争：TRAIN 准备门槛 FAIL，未执行 DEV

2026-10-03。固定拟合、末权重 TRAIN 读出及实际候选核验已完成。A13明显增加了正确链的正gain数量，尤其是5+ exon链；但两物种known非匹配的正gain数量也增加，违反事前继续条件。因此本分支停止，不提交DEV、不调整阈值或训练配置，不称科学／泛化／论文PASS。

## 实验与实际执行

[冻结协议](a13_chain_competition_protocol.md)规定：同A6的ChainHead H128／165635可训练参数、seed0 fresh、原1491唯一TRAIN窗／1536draw顺序权重×3／4608更新，原A3 step004608生成器及GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345冻结。原四块6144bp特征、FP32／AdamW参数、相同自然候选及正链注入均不变。仅完整链目标改为null＋已知严格span冲突竞争和负组null锚；它同时改变局部竞争与正负组质量，不是单组件因果消融，也不是损失函数首创。

| 阶段 | Slurm ID／终态 | 实际耗时及范围 |
|---|---|---|
| 固定fit | 13360834／COMPLETED 0:0 | 22分07秒；private／1RTX3090／2CPU／16GiB，45分钟上限 |
| 末head TRAIN读出 | 13360836／COMPLETED 0:0 | 6分34秒；同资源，20分钟上限，afterok:13360834 |
| 条件式raw入口首次检查 | 13360949／FAILED 1:0 | 13秒；0GPU／2CPU／8GiB，测试夹具ROOT不一致，不是拟合失败 |
| 同一检查的最小夹具修复＋实际候选核验 | 13360981／COMPLETED 0:0 | 13秒；同CPU资源、5分钟上限，独立RAW-CHECKS-R2输出 |

上述root／batch／extern已核对；首次CPU检查的extern为COMPLETED 0:0，root／batch为FAILED 1:0。每段STATUS／JOBID对应，原失败日志保留。实际fit在曼谷时间22:44:26–23:06:33，读出23:06:34–23:13:08；时间来自`TZ=Asia/Bangkok sacct`，不是执行主机时区。

fit实际training.jsonl为4608条连续step，epoch顺序各1536；三遍每个位置的window_id／物种／importance／自然和注入候选及1/0/-1计数完全相同。三份checkpoint step001536／003072／004608各2139423字节，评价只用末004608。fit循环1312.2615秒、峰值allocated129449472字节；只是缓存特征下头训练成本，不是端到端原生编码吞吐或方法速度排名。

## 实际配对 TRAIN 结果

核对A12末A6与A13末head的1491条实际windows记录：窗口及顺序、候选链几何、来源、1/0/-1标签和exon class相同；16857候选／窗口对＝16321自然＋536注入，标签2916／12566／1375。逐窗重算正确链的实际known span冲突最大gain，并核对有限分量和`gain=additive+nonadditive-null`。下表计数来自实际记录，不仅复述summary；分母是候选／窗口配对，不是独立基因数。

| 指标，A6→A13 | 拟南芥 | 水稻 | pooled |
|---|---|---|---|
| 全部正确链正gain | 952→1531／2172 | 285→508／744 | 1237→2039／2916 |
| 自然正确链正gain | 858→1350／1751 | 270→481／629 | 1128→1831／2380 |
| known非匹配正gain | 223→689／6305 | 88→329／6261 | 311→1018／12566 |
| 5+ exon正确链正gain | 134→624／854 | 21→173／232 | 155→797／1086 |
| 正确链不高于最大known冲突 | 325→249／1374 | 138→115／600 | 463→364／1974 |

八项原继续条件重算后六PASS、两FAIL：两物种的全部／自然正确链正gain都严格增加，pooled 5+正确链正gain增加，pooled不高于known冲突的数量下降；但两物种的known非匹配正gain均增加。所有分母保留，原`TRAIN_preparation_supported=false`和各布尔条件与实际记录完全一致。

未知／保护标签-1的正gain也由113→243／1375（拟南芥87→184／853，水稻26→59／522）；它不进loss的已知正负监督或本次gate，不能当真负或新基因。known非匹配也只是本协议标签0，不等于确认的生物学假阳性。

## 判断与未执行边界

事实是：本目标提高了正确链越过null的数量，也改善了这里测量的局部冲突比较；同时更多known非匹配越过null。因此它未同时满足“正确链增加且known非匹配不增加”的预声明准备条件。不能把5+链改善单独包装成准备PASS，不能将此TRAIN现象直接换成最终背景span、基因P/R/F1、独立泛化或论文竞争力。两个目标组件同时变化，不能据此断言哪个组件造成上升，或断言改某个组权重必定可解。

本轮不加epoch、不挑早期checkpoint、不扫温度／margin／组质量／阈值，不重训A13救结果。A13完整raw-FASTA推理、GFF、完整7728参考评价和ANNEVO／Tiberius／Helixer同尺对比均**未执行**；两个DEV输出目录不存在。它们不是失败的DEV实验，也不存在A13的DEV F1或科学gate结果。

等待TRAIN读出期间准备了条件式raw调用器和评价源码，复用A6的生成器／四块编码／greedy＋beam4池／tau0／owner／非重叠装配／phase导出。CPU检查覆盖failed TRAIN gate拒绝、raw规范化与双方向坐标／minus phase、零gain选择null、完整网格与特征边界、配对自然池几何／来源／顺序。首次有四项PASS和一项测试夹具失败：临时FIT仍被run_chain_competition.ROOT解释为原项目外路径。仅补齐临时ROOT替换后五项PASS／0.197秒，实际候选核验PASS；原共享GFF parser的非阻塞ResourceWarning保留，未借机改无关解析器。没有原生GEN前向或DEV访问，也没有因软件检查PASS启动被科学gate拒绝的推理。未来DEV协议仅是源码中的预留路径，本轮未冻结或启用正式DEV协议。

FULL原两项完整TRAIN拟合13358083／13358084及其原raw推理／共同评价链保持原参数和停止条件。其结果仍未完成。方法开发和论文目标继续，不把A13的工程完成当成论文交付；新的研究分支须按新结构问题另定协议，不把对本轮冻结目标的救参包装成独立证据。

## 来源与接续

权威数据位于Baobab `/home/users/j/jwang/coding-rna`：

- `outputs/M28-A13-CHAIN-COMPETITION-FIT/{summary.json,training.jsonl,step_004608.pt}`；
- `outputs/M28-A13-CHAIN-COMPETITION-TRAIN-READOUT/{summary.json,windows.jsonl}`；
- 配对来源`outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC/windows.jsonl`；
- 实际记录核验入口`scripts/experiments/M28-METHOD-RESTART/verify_chain_competition_train.py`；
- `outputs/M28-PILOT-R5/logs/M28A13RAWCHECKR2_13360981.out`及对应err保存五项检查和重算结果；原`M28A13RAWCHECK_13360949.err`保存夹具失败。

同一固定监控对话及heartbeat完成原七根cohort的自然终态交接（曼谷时间23:19:12.300），实际PAUSED已核对；A13已消费后只保留FULL原五根。新范围完整保存prompt／ID／固定target／heartbeat类型／15分钟／failed_runs_only／created_at均精确读回；23:25:33.047的同对话人工启动已完成，不能冒充新范围定时运行或完整FULL自动回传。沿用gpt-5.6-luna／low，没有新对话／cron／模型改动。正常状态不应回传；本次人工启动仍误发一次健康消息，已明确记录而不算终态交接。固定对话调度和本机在线限制见[OpenAI官方任务说明](https://learn.chatgpt.com/docs/automations?surface=app)。未提交Git、未对外发布。
