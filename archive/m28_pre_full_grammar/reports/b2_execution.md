# M28 B2 执行记录

本轮执行[b2_first_fit_contract.md](b2_first_fit_contract.md)中已固定的单seed评分位置机制对照。目标是完整拟合、全DEV配对解码与最终评价，不以提交成功、训练loss下降或合法GFF结构判定方法成功。

## 提交前实现与检查

2026-10-01提交前记录：新增fit_b2.py、infer_b2.py、evaluate_b2.py三阶段驱动，以及prefix_inference.py的单次union读出与R5重放核对。复用既有标签、候选Budget、owner、selector、GFF导出和统一评价函数。仅新增头训练；原B1与GLM冻结资产保持不变。

Python AST与三份sbatch的bash语法均通过。两项新增配对测试安排在拟合作业prologue内、任何optimizer更新之前：共同链在不同候选顺序中得到相同gain；控制接受JSON坐标表示但拒绝真实离散变动、参考参与、非有限或超过1e-5的proposal误差。测试时间计入该3小时GPU分配，不重复未改动的6项准备测试。

训练原1,536 draws的顺序及importance weight不变，固定加性控制池/首次截断观察按1491个唯一window ID在拟合内复用。每步读取固定特征、重新计算冻结表示，构造独立的新头计算图；不缓存新头logit或hidden跨optimizer更新。注入完整正链只出现在损失池，实际free及首次剪枝观察单独保存。

全DEV两种解码共用一个最终新增checkpoint。late先重放原B1实际固定池，再与early的候选组成单次最终读出union；每窗保存两池、proposal与新gain，记录控制误差和共同链gain相等。推理不读参考，评价在全部8690窗完成后才接入冻结7728主参考及6450历史副表。

B2继承B1已完成的4608次基座更新，另增加4608次新头更新；新增更新数不是整套模型的训练总数。early与late的曝光和权重匹配，但B2与历史C0/B1没有因复用旧资产而自动变成训练曝光匹配的完整方案因果对照。基座历史训练与特征成本另按R5/R4来源公开保留。

## 调度选择与资源

提交前Slurm记录显示gpu034为DRAIN，gpu035为MIXED，8张RTX3090中已分配7张；CPU已分配36/126，RAM已分配212992/503000MiB。选择private-teodoro-gpu，因为存在合格单卡及CPU/RAM空间，无GPU相关维护预约限制本轮墙时。共享分区不在本次选择中替代GPU型号或数值设置。

全局BeeGFS可用约164TiB，新产物上限仍为2GiB。下列三目录在提交前均不存在，既有159GiB特征与原权重仅只读引用。首次节点查询因scontrol每次只接受一个node参数而返回参数错误；分别查询两节点后获得上述事实，无作业因此启动或重提。

| 阶段 | sbatch | 新目录 | 分配上限 |
|---|---|---|---|
| 拟合 | M28-B2-FIRST-FIT.sbatch | M28-B2-FIRST-FIT | 1 RTX3090 2CPU 16GiB 3小时 |
| 配对推理 | M28-B2-PAIRED-INFER.sbatch | M28-B2-PAIRED-INFER | 1 RTX3090 2CPU 16GiB 8小时 |
| 评价 | M28-B2-EVALUATE.sbatch | M28-B2-EVALUATION | 0GPU 2CPU 32GiB 2小时 |

推理afterok拟合，评价afterok推理；序列资源上限为11GPUh与4独立CPUh，不改变旧实验上限。失败、控制不匹配、超时或产物超限停止并保留证据，不自动重跑、续训、延时、改变阈值/数据或取消别人的作业。Slurm日志沿用outputs/M28-PILOT-R5/logs下独立的新jobname/jobID文件，结果目录与旧R5完全分离。

## 提交与终态记录

此节在实际提交后补入命令与Slurm ID，并在整条流程终态后记录分配成本、产物及原条件是否满足。提交前没有B2正式拟合checkpoint或DEV性能结果。以上调度状态是提交前快照，不作以后重提的事实依据。


## 2026-10-01 实际提交与启动取证

实现于公开提交63730f07d4c19a0850e9e1813a7c89747bd34403固定，随后按以下命令各提交一次，无重提：

```bash
sbatch --parsable sbatch/M28-B2-FIRST-FIT.sbatch
sbatch --parsable --dependency=afterok:13312819 sbatch/M28-B2-PAIRED-INFER.sbatch
sbatch --parsable --dependency=afterok:13312821 sbatch/M28-B2-EVALUATE.sbatch
```

| 阶段 | Slurm ID | 启动取证状态 |
|---|---:|---|
| 拟合 | 13312819 | RUNNING，gpu035，1 RTX3090、2CPU、16GiB |
| 配对推理 | 13312821 | PENDING Dependency |
| 最终评价 | 13312826 | PENDING Dependency |

拟合启动取证时已分配610秒、STATUS为RUNNING；两项配对测试均通过，用时1.692秒。日志已进入第576次optimizer更新、第一遍，565个唯一固定控制窗口已生成。本快照不是最终更新数、总分配成本或模型收敛证据；未完成全部4608步和8690窗时仍不排名。

原独立监控已更新为这三个作业与新完成产物，保留gpt-5.6-luna、low、每15分钟、本机运行；正常运行/依赖等待不向主对话报告，整链完成或发生失败/超时/依赖无法满足/产物缺失时一次交接，成功后暂停该轮监控。仅记录状态和必要错误，不计算性能，不改变参数、预算、数据范围或goal状态。截至本轮配置时，更新启用已确认；新cohort的独立调度检查与完整回传尚未观察，不能用旧R5/Pro的成功回传替代这两项验证。

本机定时监控依赖机器开机、应用运行和连接可用，不保证休眠或断线时即时回传；Slurm作业本身不依赖本机SSH会话持续。最终结果须回到主对话核对工程完成、指标可信、机制条件及后续科学方向，不能由监控作出论文结论。

## 2026-10-01 拟合阶段完成

13312819及batch均COMPLETED/0:0，完整4608次optimizer更新、seed0；仅新增184369参数训练，基座可训练参数为0。原始训练summary逐字段副本见[b2_fit_result.json](b2_fit_result.json)，源文件为outputs/M28-B2-FIRST-FIT/summary.json；它不是配对推理或评价的终态结果。

training.jsonl实际4608行且step连续1至4608，每物种各2304次曝光；三遍各1536 draws、1491唯一窗口，每步基座梯度记录均为0。固定控制观察文件1491行。每遍均有1978次首截断观察曝光、3423个已知排序对及1094个被保护的未知竞争者曝光；三遍和有放回样本不称独立基因或独立重复。

| TRAIN每遍均值 | 第一遍 | 第二遍 | 第三遍 |
|---|---:|---:|---:|
| 完整链BCE分量 | 0.204134 | 0.159944 | 0.139309 |
| 前缀排序分量 | 0.045688 | 0.033282 | 0.024249 |
| importance加权total | 0.224266 | 0.173086 | 0.145381 |

final/rank是在每窗归一后、尚未乘importance的分量；total先把每窗两分量之和乘其原importance，再取全遍均值，不能按表中两列相加重建total。损失下降仅是当前TRAIN拟合变化，不证明收敛、自由链恢复、两物种DEV机制条件或SOTA收益。

Slurm实际分配2803秒，即0.778611 allocated GPUh，1 RTX3090、2CPU、16GiB，低于3GPUh拟合硬限。Slurm原始起止时间为2026-10-01 10:55:49至11:42:32，服务器CEST/+0200；对应08:55:49至09:42:32 UTC。程序内部循环2654.459570秒另列，不能替代分配成本；GPU阶段CPU/I/O已包含于GPUh，不在独立CPU评价成本中重复相加。batch MaxRSS为2193360KiB，GPU peak allocated为116296704 bytes，reserved为159383552 bytes。

三个epoch checkpoint各2226130 bytes；唯一正式权重仍为step_004608.pt，不以DEV选中间checkpoint。summary记载写入自身前合计产物95078370 bytes，这不是整条流程最终存储。拟合没有DEV前向/选择或test/Setaria使用，源summary及TRAIN-only驱动相符。

两项新增配对测试实际通过，1.692秒，时间已计入该分配。独立监控09:46:10 UTC记录拟合完成、配对推理13312821运行及评价13312826等待依赖；本轮主对话只核对已完成拟合，没有读取或排名未完成DEV输出，正常等待仍由该监控负责。整链评价、控制重放及最终机制条件尚未在主对话完成核对，不能用本节代替最终结果。缓存基线的版本及迁移解释范围补在[来源说明](../../../benchmarks/native_baselines/reports/baseline_cache_provenance.md)，旧预测和共同评价器不变。

## 2026-10-01 全范围完成后的逐参考链诊断准备

新增 `scripts/experiments/M28-METHOD-RESTART/diagnose_b2_recovery.py`，只在配对推理和评价均为COMPLETED、同一checkpoint/来源匹配后接入既有主参考ledger及保存的两臂自由池。使用原owner、selector和坐标转换重新汇总；实际free/owner/selected命中及候选数必须与正式评价完全一致，否则不输出解释。

每条主参考保留物种、strand、CDS坐标及unsupported标签，分为not_free、lost_at_owner、owner_nonpositive、owner_positive_not_selected、selected五种描述状态。各阶段同时报告both、early_only、late_only、neither和late→early状态转换，不只报告改善案例。正增益未被选择不自动归因于某个唯一竞争者；这不是新的训练、gate、独立测试或生物学因果证据，不修改原继续条件。主参考分母仍为7728，不增加参考筛选。

两份诊断源码已与准备时内容逐字符核对，前一准备轮的4项无真实数据合成测试已通过；尚未运行真实DEV诊断。无GPU入口为 `sbatch/M28-B2-REFERENCE-RECOVERY.sbatch`，固定2CPU、32GiB、10分钟，输出新目录 `outputs/M28-B2-REFERENCE-RECOVERY`，不得覆盖或重提。先由主对话核对整链结果及Slurm终态成本，只有评价实际CPUh加本诊断最大1/3 CPUh仍不超过原4独立CPUh硬限时才提交；不足则报告，不自动扩预算。诊断文件计入原2GiB新增产物总限。

当前只准备执行入口，没有提交额外作业；三阶段正式代码、权重、冻结协议和既有作业不变。完整结果交接仍由固定监控对话负责。

## 2026-10-01 全DEV正式结果与冻结条件判定

三个作业13312819/13312821/13312826及batch均COMPLETED/0:0，三目录STATUS均COMPLETED。完整4608步拟合与8690方向窗配对推理已结束；每臂四scope文件行数1909、1909、2436、2436，总8690，不以部分窗口排名。主参考仍为7728条（5437/2291），历史6450副表保留。评价已逐窗验证ID/自然网格、重建原参考ledger、回读GFF3身份与数量，并要求selected TP与正式exact指标一致；主对话另按整数分子分母核对F1/macro及gate布尔值，均一致。

全8690窗late离散控制与原R5精确重放，旧proposal最大绝对差0.0（固定容差1e-5）；9092条共同候选窗口记录的最终gain精确相同且使用一次union readout。推理源码只读固定特征中的DNA，不接入参考或TRAIN标签；没有test/Setaria使用。远程Git对象库不含公开freeze提交63730f0（remote HEAD仍01fbbf1且本轮相关源码为未跟踪资产），不能用本机git diff证明快照。已通过GitHub固定提交读取fit/infer/evaluate、四份prefix模块、assembly、common_ruler与evaluate_pilot十份实际文件，与当前远程内容逐字符相等；未fetch、切分支、改源码或提交。该核对不能代替原作业的全过程见证，但修正了Git对象不存在造成的取证限制。

| 同一checkpoint/物种 | free exact | owner exact | final TP / predictions / references | P | R | F1 | 背景span覆盖 |
|---|---:|---:|---|---:|---:|---:|---:|
| late / 拟南芥 | 1889 | 1257 | 1033 / 2261 / 5437 | 45.69% | 19.00% | 26.84% | 6.51% |
| early / 拟南芥 | 1932 | 1649 | 1430 / 8165 / 5437 | 17.51% | 26.30% | 21.03% | 69.81% |
| late / 水稻 | 1315 | 1079 | 874 / 1597 / 2291 | 54.73% | 38.15% | 44.96% | 2.55% |
| early / 水稻 | 759 | 665 | 580 / 3447 / 2291 | 16.83% | 25.32% | 20.22% | 64.49% |

early的拟南芥free仅净增43、owner净增392、final TP净增397，但新增输出5904条使精度和F1下降；水稻free/owner/final TP分别净减556/414/294。两物种背景span均明显增加。预设四项中拟南芥仅free/owner增加通过、水稻四项全部未通过，整轮机制继续条件FAIL。不能以合并TP增加103、拟南芥召回改善、GFF结构合法或训练loss下降替代原条件；不扫阈值/K/seed、不追加训练挽救该版本。这是此冻结off-policy头/预算/DEV的未支持结果，不是学习搜索一般不可能有效的证明。

同头late与early的macro exact F1为35.90%和20.62%，micro-pooled为32.92%和20.79%（late1907/3858、early2010/11612，参考7728）。最终输出unstranded CDS-base F1：late拟南芥49.94%/水稻67.63%，early59.67%/45.86%；pooled56.28%/53.92%，并非原局部head质量。背景全空跨度链/Mb：late20.41/10.52，early118.41/30.70；背景CDS覆盖分别51706/211778bp与715377/2911246bp，背景分母8183405/19773036bp。背景来自所有注释gene span的补集，不是经生物学裁决的真负例；nonexact不能统称生物学假gene。

同口径缓存exact F1（拟南芥/水稻/macro/micro）：C0 44.33/36.73/40.53/42.03%；B1 25.67/40.67/33.17/30.84%；ANNEVO 76.28/74.45/75.37/75.72%；Helixer 65.19/57.11/61.15/62.69%；Tiberius 74.18/66.58/70.38/71.80%。late较旧B1有所改善，但继承B1再新增4608更新，非曝光匹配因果对照；两种B2均未给出超越上述缓存基线的完整方法证据。外部方法和C0/B1没有重新原生执行，不称最新版本SOTA排名或独立zero-shot。所有late3858条与early11612条输出均通过当前CDS结构诊断，不用于事后删输出或论证正确gene。

实际root分配：拟合2803秒、配对推理16546秒（各1 RTX3090），合计19349秒=5.374722 allocated GPUh，低于11GPUh限；评价282秒×2CPU=0.156667独立CPUh，低于4CPUh限。GPU阶段CPU/I/O已计入GPUh，不另重复加到独立CPUh。推理程序16475.301209秒；共享base、union、I/O和控制的stage计时不能当两次独立端到端速度对照。推理batch MaxRSS1701656KiB，GPU peak allocated114814464、reserved146800640 bytes；评价MaxRSS2580008KiB。三目录文件实际95080645+1194195619+30418416=1319694680 bytes（1.229061GiB），低于2GiB，日志另位于R5/logs，不拿du目录元数据代替文件预算。无作业失败、重提、延时或阈值/协议变更。

现在满足既有逐参考链诊断的启动前提：正式评价成本0.156667 CPUh加其最大1/3 CPUh=0.49 CPUh仍低于4CPUh，当前磁盘余量827788968 bytes。只提交原诊断入口一次，保持7728分母、原owner/selector，并要求全部stage计数与正式结果一致；这用于解释阶段变化，不增加gate或新模型推理。A仅保留工程准备，真实A首训/政策/预算仍待另行冻结批准。

## 2026-10-01 逐参考链诊断完成

实际提交 `sbatch --parsable sbatch/M28-B2-REFERENCE-RECOVERY.sbatch` 一次，job13318144及batch均COMPLETED/0:0，STATUS=COMPLETED，无重提。沿用原2CPU/32GiB/0GPU/10分钟入口；source AST及bash语法通过。references.tsv为7728条主参考加一表头，物种/参考ID粒度唯一；全scope的free/owner/selected命中与候选数精确重放正式评价。另以标准库回读TSV核对summary每个交集及状态转换，并核对保留原FAIL条件。仍保留106条带unsupported标签的参考（两物种各53），没有新参考过滤、模型推理或test/Setaria访问。

| primary参考状态（整条链） | late | early |
|---|---:|---:|
| 不在实际free池 | 4524 | 5037 |
| free存在但owner丢失 | 868 | 377 |
| owner增益非正且未选 | 134 | 4 |
| owner增益正但未选 | 295 | 300 |
| 最终exact selected | 1907 | 2010 |
| 总分母 | 7728 | 7728 |

free交集both1528、early-only1163、late-only1676、neither3361；owner交集both1000、early-only1314、late-only1336、neither4078；最终both865、early-only1145、late-only1042、neither4676。因此最终净增103不是只增加103条而无损失：新增1145条同时丢掉1042条。新增来自late原not_free839、lost_at_owner302、owner_positive_not_selected4；丢失去向为early not_free909、lost_at_owner56、owner_positive_not_selected77。909/1042（87.24%）旧TP丢失已经发生在free可达性阶段，不能靠最终gain阈值找回。

拟南芥最终both471、early-only959、late-only562；旧TP丢失中485已not_free、31 lost_at_owner、46 owner-positive未选。水稻最终both394、early-only186、late-only480；旧TP丢失中424已not_free、25 lost_at_owner、31 owner-positive未选。这些是保存池的描述状态，不把正增益未选唯一归因于某个竞争者，也不把off-policy分布偏移直接断言为根因。

实际分配204秒×2CPU=0.113333独立CPUh，MaxRSS2174316KiB。加正式评价后累计0.27独立CPUh，GPU仍5.374722h，均低于原4CPUh/11GPUh硬限。诊断目录文件1885393 bytes；三正式目录加诊断1321580073 bytes，八份本cohort日志63902 bytes，合计1321643975 bytes（约1.231GiB），低于原2GiB新增磁盘限。目录中的详细结果与逐链TSV是解释资产，不是独立验证集或论文机制成功证据。

本轮决策：关闭当前B2前移版本的继续资格，保留late/early全部结果及诊断，不增加训练、改变阈值/K或在DEV上补救。下一方法方向保留原A gene-query/坐标指针备选，先固定可审查的真实标签、容量/损失、操作点、曝光比较及资源协议；A首训仍须既有明确批准门槛，不把合成工程PASS当注释性能。

监控核对发现14:29及14:44两个固定心跳虽追加到原监控对话，但空squeue被误当旧RUNNING/Dependency且未查sacct；这些观察已明确更正。修正后的14:59:01UTC定时心跳实际追加到同一固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3，先查squeue再查sacct，三root均COMPLETED/0:0，完成字段/产物核验通过，并实际向原主对话发送一次接续消息；15:01:26UTC轮次完成。因此配置、定时单次检查、完整自动回传均已验证。其收尾却调用delete而非pause，删除心跳及自动化状态记录。只读本机应用实现发现完成心跳的内置指令要求delete，与本项目暂停保留约定冲突；没有改应用或尝试用低优先级prompt绕过。主对话随后用原名称经正式自动化工具恢复原ID coding-rna-r5-light-monitor，明确PAUSED，原target_thread_id、完整修正prompt和15分钟间隔逐字段匹配，未新建或替换监控对话，模型仍由原gpt-5.6-luna/low对话继承。依据实际定时轮次重建本cohort交接去重记录，旧状态全文没有假称完整恢复。当前没有活动轮询，也没有第二次交接；未来自动结束后能否直接保留暂停仍受产品指令限制，不能把本次人工恢复视作已保证的自动pause。
