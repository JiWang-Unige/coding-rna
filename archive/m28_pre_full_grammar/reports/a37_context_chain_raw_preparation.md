# A37 raw FASTA 注释入口准备

2026-10-06。A37联合上下文训练及完整自然读出仍在运行／依赖等待，本步只准备raw-FASTA接口及人工工程检查，不执行真实DEV。原[训练协议](a37_context_chain_protocol.md)、[自然读出与核验](a37_context_chain_readout_protocol.md)及五项继续条件不变；TRAIN FAIL则关闭模型，不用此入口或人工检查放行，不挑checkpoint、调阈值或增加曝光。

## 自然提议与完整链评分

infer_context_chain.py沿用原双方向窗口、四块GEN编码、尾窗valid-length、nearest-owner、几何去重最大分数、span非重叠装配及GFF3导出。实际复用的annotate_query.encode_window会按原4×6144bp分块，核对每块1024 six-mer及与整体token对齐，原GEN hidden2048维经与TRAIN缓存相同的BF16 materialization后转FP32；没有另改上下文长度或数值口径。每个方向窗口只运行一次GEN编码，原A3末004608与A37末004608共享该输出及DNA，但分别执行自己的上下文。A3只用原query greedy＋beam4／completed4／cap128产生自然链；A37以自己的SharedContext和既有有序exon／junction头计算FP32 gain。旧A3发射或query存在分不参与A37新评分。

保存全部自然候选与来源，不先删非正链；保存exon additive、junction additive、两者和、链nonadditive及窗口null，最终gain＝additive＋nonadditive−null，操作点0不变。组件保存值可按既有CPU分量口径核对，不把A36的15类局部分数或FP64单位累积套用到A37。自然生成及新评分不读取参考、标签或feature cache，没有训练或阈值扫描。

两个组件仍带额外参数、原各4608训练曝光和第二个上下文／整链前向，不是容量或总计算匹配消融。共享GEN只去除重复编码，不解决A3未生成的正确几何；本入口仍仅支持原canonical primary CDS，不支持UTR、全部异构体、已证实的新基因或独立泛化主张。

## 原生执行边界

读取原生模型、GEN、FASTA或创建输出之前，必须核对实际fit13410061、自然13410739、CPU13410742均COMPLETED0:0及对应STATUS／JOBID／summary，固定612324参数／末004608／原自然16321候选、四项实际CPU证明、源作业身份、两物种计数和原2916／1086分母。重算原五条件，不能仅凭TRAIN_preparation_supported字段放行；CPU保存分量核验不是独立神经前向复现。

之后仍须另行固定a37_context_chain_dev_protocol.md及实际数据访问授权。raw接口检查阶段未创建此DEV合同或提交原生推理／评价，不放行独立物种、Setaria、A33待答复例外或公开发布。合同存在仅是科学执行前提，不是人类许可的替代。源TRAIN PASS只允许推进完整原生评价准备，不支持SOTA或投稿主张。

显式FASTA必须至少有一个非空contig，contig ID不能重复，否则不同序列会共用错误的GFF坐标；这些实际支持输入错误会终止而不写完成GFF／summary。中途失败保留诊断输出，不覆盖旧结果或自动重提。

## 针对新接口的检查

五项人工检查分别覆盖实际源身份／较强基线的严格五条件及长链分母；失败门槛或缺少独立合同在原生模型／FASTA／输出前停止；双方向规范化、同一GEN与DNA对象进入两个组件、新nonadditive整链gain进入选择、负链半开CDS坐标／phase及GFF往返；非正链选择null、空自然池和错误精度不改变原候选合同；空FASTA或重复contig不写虚假完成注释。

使用原QueryPointer人工fixture、手工整链分量和模拟源收据，不读取科研checkpoint、真实GEN、FASTA、参考、缓存或DEV，也不创建optimizer。非加性fixture只证明接线，不证明训练后效果；人工特征不证明原生GEN数值。检查不恢复或重跑A36，也不扩展既有集合DP／梯度测试。

SBATCH沿用实际A37 generanno环境（不是旧cluster_config默认环境）、private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟／no-requeue、原排除节点和独立%x_%j日志。唯一短检查输出M28-A37-CONTEXT-CHAIN-RAW-CHECKS；只需短CPU运行和小型临时人工记录，不建新缓存。短作业由主对话直接消费，不加入或替换当前13410061→13410739→13410742的固定轻量监控。

## 首次检查与人工 fixture 修复

四个任务文件已明确同步到canonical，两个Python AST及SBATCH语法通过，原生DEV合同与两个DEV输出根均实际不存在。提交前完成五项实际调度／预约检查，没有相关维护、并发或输出冲突；仅提交sbatch --parsable sbatch/M28-A37-CONTEXT-CHAIN-RAW-CHECKS.sbatch一次，真实ID13411128。

首次root／batch实际FAILED／1:0／0:25，extern COMPLETED／0:0；private0GPU／2CPU／8G、gpu034、batch MaxRSS756056K及STATUS／JOBID相符。五项中四项通过，负链人工用例未选到目标链（内部0.617秒）。失败根与日志保留，不当作通过或训练模型性能。

实际completed_union源码返回list外层、tuple exon，人工manual_scores却把目标链与嵌套list比较，所有目标链因此被赋成负分。这是人工fixture的几何比较错误，不是A37真实JunctionChainHead或坐标装配错误。唯一修改是该人工条件用tuple(map(tuple,c))比较原几何；真实caller、自然池、末模型、训练、阈值、五条件和数据范围不变。

至多一次修复复核采用独立M28-A37-CONTEXT-CHAIN-RAW-CHECKS-REPAIR-R1目录及独立M28A37RAWR1日志，相同generanno／0GPU／2CPU／8G／5分钟，不覆盖首次失败。真实源TRAIN整链仍由原固定监控负责，R1只检查人工接口，不执行原生DEV或optimizer。

## R1 实际复核

仅提交sbatch --parsable sbatch/M28-A37-CONTEXT-CHAIN-RAW-CHECKS-REPAIR-R1.sbatch一次，返回13411198；即刻squeue RUNNING／gpu034。真实root／batch／extern全部COMPLETED／0:0／0:20，AllocTRES为private0GPU／2CPU／8G，batch MaxRSS734872K，STATUS／JOBID相符。五项检查全部PASS，内部0.667秒；未删除或覆盖首次FAILED13411128。

人工目标链的additive−null为负，但nonadditive分量使最终gain28进入原装配和最终transcript；两方向共享一次人工GEN输出和同一DNA对象，未调用A3旧发射评分。负链CDS[(0,5,phase2),(10,14,phase0)]、保存分量的FP32组合及实际人工GFF往返一致。原零／负分保留候选却选null，空自然池和错误精度按原合同处理，门槛／合同失败与空／重复contig不产生虚假完成注释。

原scripts/eval_structure_diagnostic.py:51再次出现unclosed file ResourceWarning，保留stderr；本次GFF往返通过，不把此警告解释为科学FAIL或扩大到旧解析器重构。修复仅是新人工fixture的几何比较，真实caller及训练／读出源码未因人工失败修改。

本次真实调度读取fit13410061仍RUNNING／49:58／gpu034，下游13410739／13410742正常Dependency；同一固定监控15:30的定时turn01a111d6-36f0-7892-acfc-f2ce57c656e1正常结束／error-null／28.074秒，未向主交接健康状态。两个短人工检查均已由主消费并从长期监控中排除，原三根监控／模型／15分钟／冻结五条件不变。真实训练末模型、自然选择门槛、原生完整注释、独立泛化和论文目标仍未完成。

## 完整评价入口准备

evaluate_context_chain.py重建原4个方向scope／8690窗口，核对实际两末模型、原4块特征、612324＋649457参数、源TRAIN作业及零新训练／无参考推理。逐窗与实际A31_A3保存记录对照，要求候选几何、顺序、greedy／beam4来源相同，只允许分数改变；完整保存的FP32 exon＋junction＝additive、additive＋nonadditive−null＝gain以原CPU规则精确核对，不把保存分量核验称为独立神经前向复现。空候选窗口仍必须保留，不能缩小评价范围。

全部新自然记录、owner／去重／span选择、负链CDS坐标／phase及最终GFF往返关闭后才调用reference_context。原7728条primary CDS-chain参考、原背景bp分母及生成器不支持的参考层不删，free／owner／selected阶段的覆盖损失另列；selected阶段TP须与最终exact-chain TP一致。对照A31_A3、A31_A23、A3、A6、A、R5_C0及ANNEVO／Tiberius／Helixer沿用原完整尺，报告逐物种和pooled P／R／F1、物种macro F1与最差物种F1，不用候选窗口计数当基因级F1。A31_A23改变候选生成器，不作为同池机制消融；额外容量、监督曝光和计算不支持因果创新、独立泛化、公平原生运行时间或论文PASS。

五项针对性人工检查覆盖完整网格与模型／源身份；同池来源及精确FP32分量；无参考负链装配／GFF；所有对照的完整参考与背景分母、macro／worst；空候选窗口和部分文件拒绝。只使用人工metadata、手工分量与临时GFF，不读真实checkpoint、GEN、FASTA、reference、cache或DEV，不创建optimizer。

短SBATCH M28-A37-CONTEXT-CHAIN-EVAL-CHECKS沿用generanno／private0GPU／2CPU／8G／5分钟／no-requeue和原排除节点，独立输出与日志；由主直接消费，不加入或替换13410061→13410739→13410742的同一长期监控。完整评价人工检查时仍未创建独立DEV合同或提交原生DEV。新入口和人工通过只属于评价工程准备，不放行TRAIN FAIL。

## 完整评价工程检查实际收据

新脚本／五项人工检查／SBATCH／本记录四个任务文件明确同步到canonical，远端两个Python AST及SBATCH语法通过。原生DEV合同与DEV推理／评价输出根实际仍不存在；新短检查根无冲突。提交前完成五项实际调度与预约检查，没有相关维护、并发或输出冲突，不改其他项目作业。

仅提交sbatch --parsable sbatch/M28-A37-CONTEXT-CHAIN-EVAL-CHECKS.sbatch一次，实际ID13411914。即刻squeue／scontrol为RUNNING／gpu034／private0GPU／2CPU／8G／5分钟／Requeue0／Restarts0；最终root、batch、extern均COMPLETED／0:0／0:20，batch MaxRSS731396K，STATUS／JOBID一致。五项检查全部PASS，内部0.140秒，不读取科研模型／GEN／FASTA／参考／cache或DEV，不新增optimizer曝光。

负链保存组件经gain重建到最终GFF一致，错源／早checkpoint／改变候选几何或来源／错误FP32组合会拒绝，原7728参考及背景分母、十项对照的macro／worst计算和空池窗口完整性均由人工用例验证。原scripts/eval_structure_diagnostic.py:51的unclosed file ResourceWarning保留；GFF往返实际通过，未借此扩大旧解析器修法。

提交前实际fit13410061仍RUNNING／1:11:50／gpu034，自然13410739与CPU13410742仍正常Dependency。短评价检查已由主消费并从长期监控排除；原同一监控、三根对象和五项科学继续条件不变。真实末模型／自然门槛、完整原生注释、独立泛化和论文目标仍待完成，评价入口通过不构成任何科学PASS。

## 后续完整 DEV 执行准备

2026-10-07北京时间提前固定[a37_context_chain_dev_protocol.md](a37_context_chain_dev_protocol.md)及原生GPU／评价CPU两个SBATCH，只有原TRAIN整链闭合并实际五项PASS才提交；不是提前执行DEV或改停止条件。合同保留原两条染色体／8690窗、两末004608、原FP32 gain0、完整7728参考及背景分母。原A31同网格实际Slurm3:13:47及既有零更新时序用于拟定原生1RTX3090／2CPU／16G／7:30上限，CPU0GPU／2CPU／8G／5分钟；这些不是原生速度实测或科学门槛。

此时原fit13410061已实际COMPLETED0:0／1:15:53，匹配STATUS／JOBID／summary、4608训练行及三个各7389798字节checkpoint；数值与分支曝光仍待原CPU核验。自然13410739实际RUNNING／5:56／gpu034，CPU13410742正常Dependency。只读取终态和固定结构，不读取部分自然性能，当前仍无原生DEV作业或输出。只复用原同一三根监控，A36关闭、A33未答复、独立数据及公开发布限制不变。
