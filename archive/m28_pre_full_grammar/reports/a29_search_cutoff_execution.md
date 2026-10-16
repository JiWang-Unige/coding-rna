# A29 实际搜索截断候选标签执行

2026-10-06。先固定a29_search_cutoff_protocol.md，再提交CPU诊断；尚未读取全量标签分布。它为下一版搜索竞争损失提供监督资格，不训练模型，也不恢复已失败A9。

## 实施与实际提交

新增analyze_search_cutoff_labels_train.py、四项专用测试及M28-A29-SEARCH-CUTOFF-LABELS.sbatch。标签采用gene-set兼容性：其他eligible primary仍为正，保护／unknown／noncallable为未决，剩余仅为参考相对known_nonmatch。正确目标和首淘汰cutoff分开，尤其保留分数界的DONOR目标与完整STOP竞争链区别。仅原352 TRAIN窗、545参考配对和4247 A28 profiles，不调用模型或GLM、不产生搜索或checkpoint，不读DEV／封存Setaria。

两Python文件AST与SBATCH语法本地通过，四个相关文件上传成功；远程再次AST PASS 2、bash语法通过。四项unittest尚未运行，会在CPU作业入口执行，不能把静态检查写成测试通过。提交前实际五项调度检查完成，没有相关维护截断或新输出目录碰撞；A28源仍RUNNING。

实际命令：

```sh
sbatch --parsable --dependency=afterok:13392879 sbatch/M28-A29-SEARCH-CUTOFF-LABELS.sbatch
```

返回job13393045，generanno，private CPU-only／0GPU／2CPU／8GiB／15min，唯一输出outputs/M28-A29-SEARCH-CUTOFF-LABELS，日志outputs/M28-PILOT-R5/logs/M28A29CUTOFF_13393045.out／err。入口还核对实际源JOBID13392879及原全量A28的全部身份／搜索核验条件；Slurm依赖本身不代替科学来源核验。

首次取证A29为PENDING（Dependency），源A28 job13392879 RUNNING gpu034／13min20sec。A23／A25自然读出13391609／13392061 RUNNING gpu034／42min24sec和20min47sec，下游CPU13391610／13392062有效Dependency。此时没有A29测试或标签结果，也没有A23／A25完整自然读出结论。

## 固定监控

同一个监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3、同一个heartbeat coding-rna-r5-light-monitor更新为M28-A23-A25-A28-A29-TRAIN-20261006，仅原A23三job、A25三job和A28→A29两job。A23／A25结构白名单逐字保留；A28与A29必须整条链成功且各自结构齐全才交接，A28单独完成不触发；任一新可行动故障仍交接。标签比例与机制解释只由主对话读取。

原生保存2026-10-05T18:25:10.465Z／updated_at1791224710465，完整prompt及id/name/kind/target/cadence/notification/created_at精确读回PASS，created_at1790981930511保留。未创建／fork对话、cron、改变模型或平台goal。固定gpt-5.6-luna／low与15分钟不变。按[官方定时任务说明](https://learn.chatgpt.com/docs/automations)继续绑定既有对话；配置成功、manual startup、当前cohort定时运行和终态自动回传分别记录，不借用旧cohort证明当前自动运行。

当前cohort manual startup已实际观察：取证2026-10-05T18:25:55.391Z，turn01a10d50-965f-7771-8467-e28b3cb3a6d8 completed/error-null／19.218秒。A23自然读出13391609已COMPLETED0:0，CPU13391610正在gpu034执行；A25自然读出及A28仍运行，A25与A29的CPU有效Dependency。没有完整链或新故障，不发送健康回执。实际turn_context18:25:47.645Z仍gpt-5.6-luna／low、相同cwd，cursor29057c64-b701-4669-a8f2-a99bf1e621bb:52已消费。当前cohort定时运行及终态自动回传尚未观察。

完整A28／A29现已由主对话直接读取。A29 job13393045 root/batch/extern均COMPLETED0:0／3min44sec／gpu034，batch MaxRSS782872K，四项测试PASS／0.003秒；入口198.162735秒。cutoff_annotations2021行／references545行、STATUS COMPLETED／JOBID13393045及全部原profile／R4／保存记录核验闭合。完整资格结果见a29_search_cutoff_result.md；工程闭合不是性能PASS。

原当前cohort在2026-10-05T18:27:21.633Z实际定时取证过，turn01a10d51-d204-7f80-9ed9-6e25fdf5ae4c completed/error-null／28.466秒，cursor53已消费，健康状态未给主对话发消息。A28／A29终态由主对话直接消费，不虚称已验证自动终态回传。旧journal保留其manual及scheduled记录。

同一原生heartbeat现仅复用A23／A25六root，cohort M28-A23-A25-TRAIN-POST-A29-20261006，两个结构白名单逐字不变。保存2026-10-05T18:34:15.704Z／updated_at1791225255704，完整prompt和id/name/kind/target/cadence/notification/created_at精确读回PASS，原created_at保留。A28／A29及短测试13393178全部移出，避免再次交接已消费事件。

新cohort manual startup已验证：取证2026-10-05T18:35:27.178Z，turn01a10d59-4e6a-7a33-8335-36374d46a0db completed/error-null／18.443秒；actual turn_context18:35:19.059Z仍gpt-5.6-luna／low、相同cwd，cursor54已消费。A23 CPU及A25自然读出仍运行，下游A25 CPU有效Dependency，没有健康回执。新cohort定时和终态回传尚未观察，不借旧定时记录。

A23仍需CPU核验及主对话读取冻结门槛，不能从自然读出完成推断PASS。后续需结合已完成A28／A29和A23／A25完整自然读出，决定下一独立模型与监督保护规则，不据提交或组件测试成功宣称科学PASS或论文创新成立。
