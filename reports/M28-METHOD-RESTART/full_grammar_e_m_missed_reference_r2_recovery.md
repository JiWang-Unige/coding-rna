# 原 TRAIN 漏检诊断的超时恢复

2026-10-13。13554511在原1h时限内只完整完成四条染色体，实际TIMEOUT，不能解释为全诊断成功或科学负结果。按用户已授权的常规实施／结果恢复范围，保留失败与全部原产物，独立核验四个已完成scope，只并行重新计算余下两个scope，再合并原22,688主参考。没有延长／取消旧作业、改变数学定义、权重、阈值或数据范围。

## 失败证据

R1诊断09:54:26UTC开始，Slurm于10:54:54达到时限终止：root TIMEOUT／ExitCode0:0，batch CANCELLED／0:15、10:54:55结束，extern COMPLETED／0:0、10:54:56结束。MaxRSS10430816K，低于16GiB申请，不是OOM；stderr明确DUE TO TIME LIMIT。旧STATUS仍RUNNING，源自时限信号未完成收尾，以Slurm为准，不覆盖旧标记伪造成功。

stdout确认NC_003070.9／NC_003071.7／NC_003076.8／NC_089043.1完整双模型scope完成，各有两模型reference_energy.jsonl和summary。第五条NC_089044.1只写出E summary、M未完整；第六条未完整，两者均重新计算，不借单模型半scope完成推定整体完成。已有十二posterior与评价13554510均完整成功，后验排序结果仍有效，R1漏检全结果并不存在。原watcher已真实一次回传该超时，主收到并消费；恢复沿用同一监控对话／心跳，不创建重复监控。

## 不变计算与复用核验

从原诊断入口提取单染色体函数，角色、literal-DNA支持、原GFF／selected-set gain重放、重叠冲突、replacement margin、1e−5近零报告和分类保持原代码；原六scope主入口及最终聚合行为保留。提取前后按确切文本变换核对，除缩进、scope结果返回和原all_rows bookkeeping移动外无计算改动；AST及SBATCH语法核对后才执行。

一个2CPU／4GiB／5min复用核验，读取原master started来源及四个完整scope的八份JSONL／summary，不读DNA／势／checkpoint或做forward。每scope要求全参考几何key（包括species／seqid／strand／全部CDS）唯一且与原13449803 profile直接相等、reference数／exact TP与原E/M评价相等，counts和role terms从逐行数据重算完全相等，两strand selected-set gain重放成立；两模型support及reference-role rows相同，异常improving replacement为0。截断或缺行、partial summary、重复／错配或异常替换即拒绝复用，不用摘要签名代替真实行比较。

资格通过后两个2CPU／16GiB／30min任务各完整计算NC_089044.1和NC_089046.1的两模型／双方向，既有一条rice约8min，30min包括导入、数据读入和余量。最终2CPU／4GiB／10min聚合依赖这两个实际job，重验四个旧scope和两个新scope的所有参考行／role terms／paired口径，逐chr→species→pooled重建22,688主参考与E16420／M16106主TP；异常替换必须为0。没有参考驱动重解码或新预测，没有新训练、DEV／sealed／新panel读取。

## 独立输出与完成条件

新输出`outputs/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE-R2`。完整旧四scope留在R1目录只读，聚合以精确来源路径记录，不复制或重写；新两个scope独立JSONL／summary与STATUS／JOBID／started，聚合独立evaluation／summary。不复用失败job的全局成功标记，不写旧root总summary。新增产物声明128MiB，原256MiB诊断源及失败记录保留。

恢复成功需新两个component和aggregate共3roots／9root-batch-extern全部COMPLETED0:0、markers／非空summary和四份新逐链文件齐全、旧八份完成块通过实际复用核验、全六chr／两species／pooled参考与TP完整重建。原十三个后验／评价成功阶段及新的三成功阶段构成有效完成流程，同时明确保留原13554511 TIMEOUT；不能宣称原十四stage全部0:0。若恢复产生数学counterexample或科学负结果，保留并停止相应解释，不改规则求漂亮结果。冻DEV八／十／六家族及论文goal边界不变。

## 实际执行记录

复用资格13565487实际11:23:00–11:24:57UTC／117s，root、batch、extern全部COMPLETED0:0，batch MaxRSS658068K；程序28.111883s。原四scope共19,377主参考的八完整块全部逐行／角色汇总重建一致，单模型不完整第五scope明确未复用。此资格主对话已核验消费，不是待监控事件，watcher未发送其交接。

通过资格后，task4 NC_089044.1实际13565663、task5 NC_089046.1实际13565660并行提交，均afterok:13565487、2CPU16GiB30min，无GPU。聚合实际13565687，afterok:13565663:13565660，参数顺序为资格13565487／task4 13565663／task5 13565660，2CPU4GiB10min。实际ID来自提交返回值，不按相邻数字推测。11:27UTC两个component RUNNING，聚合等待依赖；提交和RUNNING不等于诊断完成。

## 完整结果与判断

task4实际11:26:57–11:34:54UTC／477s，task5 11:26:56–11:33:55／419s，聚合13565687 11:34:56–11:35:50／54s；三个root及batch／extern共9条全部COMPLETED0:0。两个component的MaxRSS分别2139820K／2232332K，聚合1118492K，程序34.650288s；新STATUS／JOBID／started／四份逐链文件及151002-byte聚合summary齐全。八旧块实际资格与最终全六scope重建均通过，完整保留22,688主参考，替换改进counterexample为0。主另行独立核对分类和single／multi分区、六chr→species→pooled及两模型转移矩阵行列边际，一致。有效流程16个成功正式root／48条记录，另保留旧13554511 TIMEOUT；不是原十四stage全部成功。

| 模型／范围 | 主参考 | Exact TP | 正能量但被冲突链压过 | 低于背景 | 语法不支持 |
| --- | ---: | ---: | ---: | ---: | ---: |
| E／AT | 17739 | 13111 | 2935 | 1585 | 108 |
| M／AT | 17739 | 12873 | 3048 | 1710 | 108 |
| E／Rice | 4949 | 3309 | 961 | 578 | 101 |
| M／Rice | 4949 | 3233 | 982 | 633 | 101 |
| E／Pooled | 22688 | 16420 | 3896 | 2163 | 209 |
| M／Pooled | 22688 | 16106 | 4030 | 2343 | 209 |

E漏6268条，62.156988%正能量被已选冲突链压过、34.508615%低于背景、3.334397%语法不支持；M漏6582条，对应61.227590%／35.597083%／3.175327%。两branch同支持22,479条参考，只有209条不支持，不能将绝大多数漏链归因于语法不可表示。正能量竞争类别的E／M中3565／3713为multi-CDS；这是参考相对分区，不是subset precision/F1或独立长基因结论。全链角色分解不是替换margin的因果解释：E已命中链也有16360条start/stop整体负项，M有14789条，因此“负首尾项导致所有漏链”不成立。

两branch共同命中15636、E-only784、M-only470、共同漏5798；M−E TP为−314，与既有paired相同。E→M丢784中682转为正能量冲突压过、102转为负能量；M-only470由E冲突类423和负能量类47组成。没有改变预测去验证某一角色权重的因果作用。结果支持将下一项固定TRAIN探索放在竞争能量及branch分歧上，并继续等待已运行的LR四格来区分rank／LR；不据此改DEV或放行封存数据。

原watcher11:44:27UTC实际新轮已核验R2整体成功，但先交接POOLPOST，因此R2 handoff仍NOT_SENT_BY_WATCHER；主已独立消费，不补发或伪造YES。数据验证判断为“附TRAIN探索／非因果限制可用于下一项机制决策”，不是独立性能、泛化或投稿验收。
