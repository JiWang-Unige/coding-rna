# A27 全量遗漏链诊断执行

2026-10-06。已按[a27_missing_true_prefix_protocol.md](a27_missing_true_prefix_protocol.md)启动独立A27零更新诊断，用于区分原536遗漏的局部必要排序失败与尚未定位的搜索竞争。不是新自由预测、训练或性能gate，A23/A25/A26原冻结配置保持不变。

## 实际样本和实现

直接读取原A24全2916正链／窗metadata与A6原自由池、来源和几何，确定全部536遗漏（拟南芥421、水稻115），396个原全宽TRAIN窗；原metadata第一条自然beam正链对照324个，共860配对／13760 query。遗漏3711exon加对照1320exon，最多80496完整事件。对照不并入遗漏分母，没有对照的72窗不排除。

入口scripts/experiments/M28-METHOD-RESTART/analyze_missing_true_prefix_train.py直接复用冻结原probe_query及其codon/EOS/width4语义，没有新head或decoder。新增聚合保留语法不支持、Q16和reference/query分母；实际保存后重读全部query_profiles与references，核对身份、几何、path组成、首次排除和汇总。数据只来自原六条TRAIN染色体，原GEN revision/A3末004608不变。

新脚本和测试AST、SBATCH语法通过，四个新增文件逐一与远程字节比较相同。只上传相关源码／协议，不同步仓库或提交Git。

## 已启动的作业

实际命令：sbatch --parsable sbatch/M28-A27-MISSING-TRUE-PREFIX-DIAGNOSTIC.sbatch。Slurm13392411，private-teodoro-gpu、1RTX3090/2CPU/16GiB/40min，generanno。17:33–17:34 UTC取证RUNNING gpu034，root目前未终态；既有三训练链正常运行／有效Dependency等待。

提交前五项实际调度检查显示gpu034/035各有一块可用3090，private无新增优先候选、维护不截断本run，目标新目录不存在，BeeGFS余136TiB。旧A5实际3456 exon事件／432 probes／57.356秒，本次按事件与396 cache窗估约25–30分钟，40min包含操作余量；这不是速度benchmark。

实际日志M28A27PREFIX_13392411.err记录原四项真实前缀测试PASS（1.100秒），新增三项聚合检查PASS（0.002秒）。尚不能从测试通过或RUNNING推断全数据探针已经完成。

首次实际数据进度日志已取证：8窗／240探针／40.455秒，16窗／528探针／79.117秒，24窗／47配对／752探针／99.915秒。这证明当前数据loader和实际冻结头探针已推进，不只完成合成测试；尚非396窗全量结果，不据中途记录判断机制比例或调整范围。此进度与40分钟上限没有当前冲突，未延时或重提。

输出目录outputs/M28-A27-MISSING-TRUE-PREFIX-DIAGNOSTIC，要求summary.json/STATUS/JOBID、references.jsonl860行、query_profiles.jsonl13760行齐全；完成时还需root/batch/extern ExitCode0:0，全部来源与实际记录核验一致。

## 监控与下一决定

同一监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor已扩展为M28-A23-A25-A26-A27-TRAIN-20261006：原A23三job、A25三job、A26两job全部保留，新增A27单job，总九job。沿用gpt-5.6-luna/low、15分钟、终态／可行动故障一次交接；不新建、fork、cron或改主模型。

17:33:58.974Z原生保存；完整prompt/id/name/kind/target/cadence/notification/created_at精确读回一致，原created_at保留。同对话manual startup已验证：取证17:34:42.304Z，turn01a10d21-a066-7412-b8db-531ec37c2057 completed/error-null/19.137秒，四项root运行、下游有效Dependency，未回传终态；实际turn_context17:34:30.037Z仍gpt-5.6-luna/low。新cohort定时和最终自动交接尚未观察，不借用旧cohort定时证据。

完整结果未出，不能指认具体错误边界、EOS或长度偏好。若多数遗漏全部query局部排除，下一独立机制优先相应局部边界；若存在较多局部通过query，进一步真实搜索竞争诊断仍必要，不能据本探针直接改变beam或归因。最终选择还须结合A23/A25/A26真实冻结结果，不救参旧FAIL配置。

完整raw-FASTA竞争性模型、独立测试／适用范围、容量匹配机制消融和论文创新仍未完成。本执行只是减少下一迭代的不确定性，不能称Nature Communications目标达到。

## 实际完成与接续

18:00–18:02 UTC主对话取证：13392411 root/batch/extern均COMPLETED0:0，27分05秒，batch MaxRSS1508844K。summary/STATUS齐全：396窗、860配对、13760 profiles、80496事件，实际保存后重读核验PASS；324自然beam对照均满足必要条件。536遗漏中315全部query局部排除、221至少一个query局部通过，5+类157/298局部通过。结果和解释边界已另记a27_missing_true_prefix_result.md。

该221配对仍需真实搜索竞争证据，主对话开始独立零更新回放实现；不改变原A23/A25训练、原A26 FAIL或本诊断的冻结范围。A27已由主对话消费，后续同监控移除该完成对象，不重复报告旧事件。
