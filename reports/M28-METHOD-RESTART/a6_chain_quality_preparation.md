# M28 A6 完整链质量评分准备

2026-10-03。A5固定A3末权重的八窗TRAIN诊断新增7个正确参考窗配对（8→15/27），但候选29→94，且beam不保留全部greedy候选。当前假设：保留同一冻结生成器的greedy+beam4联合池，用整链质量与window-null判别替代将query existence分机械复制到多条路径。这是独立A6准备，不修改A3 FAIL或运行中的全TRAIN Q/C0、原生基线。

## 本次事前固定范围

只做零更新梯度preflight，不开始拟合、选择权重或计算DEV性能。生成器唯一权重为已完成A3 step004608（GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345）。原八条TRAIN压力哨兵及A5已完成自由预测快照只读复用，不重复解码未变化候选。greedy+beam联合去重按完整有向CDS几何排序，保留每链greedy/beam4来源；仅terminated greedy和真实beam completed进入free。beam4/completed4/max_exons128及原DNA语法不变，不使用路径分或参考信息筛选联合池。

先生成并保存全部八窗联合free池，关闭文件后才读取原R4 TRAIN窗口/catalog，旧缓存标签不作为监督。继承eligible positive、alternative/protected、unknown/callability的1/0/-1政策；0是当前主参考相对非匹配，不是生物学真负。仅TRAIN loss加入缺失正参考链，保留free/reference_injected来源；不把注入链算作生成命中，未知-1不进入BCE或负例。

评分器直接复用src/m28/core.py::ChainHead，H128、FP32、seed0：exon特征及additive、ordered exon GRU/nonadditive、window-mean null，完整链logit=additive+nonadditive-null。link模块在本完整链目标中不使用并冻结，不新增局部/link/segmentation损失。不把beam、已有链GRU或本接线当作新颖性证据。

生成器eval且全部requires_grad=false，普通torch.no_grad输出detached hidden供新头反向；不使用不能供head backward保存的inference tensor。固定每窗known-candidate-mean BCE；不按结果切换class balancing、长度归一或阈值。8窗各一次反向，optimizer不存在、更新0、checkpoint写入0，最终直接比较评分头参数与初值，确认未变。

## 具体检查与下一决定

三项合成检查只覆盖实际可达风险：联合池遗漏原greedy或重复计数、TRAIN注入覆盖free来源、保护-1误入负梯度/全未知失去可微零损失。实际八窗检查输入几何/DNA与固定快照一致、free与注入计数、score/loss有限、exon/GRU/additive/nonadditive/null各分支梯度有限且非零、冻结生成器/unused link无梯度、评分头参数未变；对应错误会停止后续拟合并修正具体接线，不改变科学规则。

入口prepare_chain_quality.py及M28-A6-CHAIN-QUALITY-PREFLIGHT.sbatch；generanno/private-teodoro-gpu、1 RTX3090、2CPU/16GiB、3分钟，新目录outputs/M28-A6-CHAIN-QUALITY-PREFLIGHT。八窗已有预测不重复解码，以已完成同环境诊断56秒和原生烟测25秒为小规模运行依据，加启动余量；这不是正式拟合时长估计。只保存free_pool.jsonl、records.jsonl与summary，不保存权重。沿用Slurm动态路由，不更改现有作业、监控、预算或数据范围。

成功仅说明可以进入固定TRAIN拟合方案准备，不能证明训练校准、raw-FASTA推理或注释性能。正式拟合前仍需明确额外曝光、训练顺序、最终checkpoint、raw-FASTA出口及完整7728参考/8690窗口、逐物种F1和背景gate；不能把追加训练后的提升称为匹配曝光的单因素效应，不扫描DEV或开放test/Setaria。重要终态结果追加本文，不创建公开托管、Git提交或论文新主张。

## 实测终态与准备结论

`sbatch --parsable sbatch/M28-A6-CHAIN-QUALITY-PREFLIGHT.sbatch`返回13358869，Slurm root/batch/extern均COMPLETED/0:0，36秒，STATUS COMPLETED、JOBID一致；private-teodoro-gpu/gpu034、1 RTX3090/2CPU/16GiB/3分钟、generanno。三项合成检查通过（0.308秒）。五个新文件的远程精确读回、三份Python AST及sbatch bash语法检查通过，不更改原A3/全TRAIN/原生作业。

八窗全部完成，联合free候选逐窗数为20、16、5、14、13、16、9、11，加总104。参考注入逐窗7、3、0、0、0、1、0、1，加总12，训练池116。自由池标签精确主参考匹配15、参考相对已知非匹配88、保护/未知1；注入链全部为正标签，故训练池27/88/1。每链来源与去重几何在raw记录中保留，注入的12链没有计入free的15条命中。保护-1没有转成负例，0仍不作生物学真负主张。

评分头可训练165635参数、未使用link冻结33153参数；生成器0可训练参数。八窗exon、ordered GRU、additive、nonadditive、null各分支梯度均有限且非零，loss与每链初始logit均有限；生成器与link无梯度，hidden为普通detached tensor而非inference tensor。实际optimizer不存在、更新0、无checkpoint保存，评分头所有参数直接与初值比较完全一致。未训练分数和损失不作模型选择或准确率结果。

实际评分/反向循环8.8341秒，GPU峰值allocated116330496 bytes；首窗6.2994秒含初始化，后七窗0.1766–0.7315秒，不是稳态速度benchmark，也不含正式逐窗beam生成开销。产物summary.json10705 bytes、records.jsonl31410 bytes（8行）、free_pool.jsonl10420 bytes（8行），位于outputs/M28-A6-CHAIN-QUALITY-PREFLIGHT。主对话逐链核对free/注入来源、联合池与训练池计数、1/0/-1守恒、source映射及有限梯度/冻结边界，与summary一致。

结论为TRAIN评分接线/梯度准备通过，不是新fit或科学PASS。下一实施是冻结事前A6拟合与raw-FASTA推理/评价方案，使用上述同生成器联合池和单一seed，保持原标签、完整DEV分母与背景约束；正式训练和DEV仍未执行。A5候选新增提供机制依据，不保证整链质量能在独立注释中获益，且≥5 exon候选遗漏尚未解决。论文目标仍未达成。
