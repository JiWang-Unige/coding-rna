# A31 完整 DEV 正确链分数与选择损失诊断

2026-10-06。原评价把owner后的正确链损失合并为“冲突或非正gain”，不能区分评分没有越过null和正分数链未被选择。本诊断复用既有A6/A7五状态分类及原装配规则，只读取完整闭合的A31输出；不改变模型、候选、owner、阈值、参考或原门槛。

## 来源和执行条件

A31 A3 infer13399434/eval13399454、A23 infer13399441/eval13399455及旧A6/A7分数分类13360116，五个root均须实际COMPLETED/0:0。原生两条链未闭合前不打开A31输出或参考表，不预提交真实诊断。每个来源的STATUS、JOBID与原协议身份须对应。A31 A3已核验同A6自然池；A23不是同池。

只读每个variant已保存的四scope、8690窗完整自然候选及其原分数，以及旧分类保存的全部7728条primary参考。沿用拟南芥5437、水稻2291的完整分母及unsupported strata；不读取FASTA、特征缓存、模型权重、catalog或R4标签，不重新生成标签或预测。

## 分类和核对

按原nearest-owner、重复几何最大分数及同strand span非重叠正gain选择重建。逐参考记录五个互斥状态：not_free、lost_owner、nonpositive_gain、positive_unselected、selected。保留原owner分数与unsupported strata，逐物种、exon层级和合并汇总，并报告相对于旧A6的参考状态迁移。

候选/owner/选择计数须精确复现原inference；五状态合并须复现原evaluation stage计数，selected须等于原完整链TP。A31 A3逐参考free及owner成员关系须与已闭合A6一致；不将此约束施加到A23。任何不一致先视为分析实现或来源错配，不给出科学解释。

positive_unselected只表示原选择器没有选中正分数正确链，不能直接归因为真实基因互相冲突或某一训练缺陷。A31 A3对A6固定池但改变评分机制和head训练曝光，A23还改变生成器；状态迁移不是纯评分因果效应。诊断不增加F1或通过门槛，不做阈值扫描、checkpoint选择、补训练或发布；重复使用DEV不能证明独立泛化。

## 验证和资源

四项纯合成检查覆盖：zero/negative与正分数竞争排除的区别（错误竞争者压过互不重叠正确链，且预测数与参考TP不同）、non-owner不补回与重复最大分数、负方向坐标和owner平局、来源RUNNING时在真实输出读取前停止。检查和真实诊断均为private0GPU/2CPU/8G/5分钟，沿用generanno，日志唯一；实际诊断只解析已保存记录及7728行参考，预计新输出小于10MiB。

当前仅准备并运行合成实现检查。首轮13400815已COMPLETED/0:0、四项PASS，原输出outputs/M28-A31-DEV-SCORE-ATTRITION-CHECKS保留。随后修正首个用例的参考集合，避免把错误竞争者也列作正确参考；只复验改后的合成用例，使用独立outputs/M28-A31-DEV-SCORE-ATTRITION-CHECKS-R1。由主对话消费、不加入原七root监控。真实诊断输出outputs/M28-A31-DEV-SCORE-ATTRITION；五来源闭合后才提交sbatch/M28-A31-DEV-SCORE-ATTRITION.sbatch，提交前核对实时调度与唯一输出。等待仍复用同一固定监控对话和心跳，不改A30/A31运行代码或已有实验产物。

## 实施回执

2026-10-06 07:44 UTC。独立诊断实现、测试、两份SBATCH及本协议已部署canonical remote；本地／远程两个Python AST与shell语法通过，提交前五项实时调度检查和唯一输出核对完成，沿用实际generanno。未改变运行中的A31或已闭合A30代码、参数和产物。

13400815的root／batch／extern均COMPLETED/0:0、37秒，四项PASS/0.033秒，MaxRSS720164K。修正合成参考集合后的13400879同样全部COMPLETED/0:0、16秒，四项PASS/0.003秒，MaxRSS721972K，实际07:40:48至07:41:04 UTC；独立R1输出STATUS／JOBID相符。两者均实际0GPU/2CPU/8G，主对话已消费，不加入研究监控。真实DEV诊断尚未提交或读取；A31两条原生链仍未闭合，当前仅证明限定实现行为，不是科学结果。

## 真实 DEV 诊断执行结果

2026-10-06 09:19 UTC。五个实际source root已COMPLETED0:0，原A31两条链及字面来源身份由主对话直接核对，随后完成最新五项调度检查及唯一输出核对，按原命令提交13402231。实际09:07:06至09:08:20 UTC，root／batch／extern COMPLETED0:0／1:14，MaxRSS765784K，private0GPU／2CPU／8GiB／5分钟；四项合成检查PASS／0.003秒，输出STATUS／JOBID匹配。

references.tsv为7728条数据加表头；保存自然分数、owner、DP选择、全部原stage计数及selected=原TP均精确复现。A31 A3逐参考free／owner与A6成员关系相同，A23不强加同池。主对话独立核对五状态总数、CDS片段层级和状态迁移边际一致，完整结果及限制见[a31_local_path_caller_dev_result.md](a31_local_path_caller_dev_result.md)。没有新模型／GEN前向、R4标注、预测、训练、删参考、threshold scan或Setaria，tau仍0；不把此诊断当额外F1提升、因果归因或新gate。短作业已由主对话消费，不加入已暂停的四root心跳。
