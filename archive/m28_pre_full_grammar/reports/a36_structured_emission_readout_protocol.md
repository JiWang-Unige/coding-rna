# A36末模型自然池读出与实际产物核验

2026-10-06，在读取新模型分数或效果前固定。沿用[A36首次拟合协议](a36_structured_emission_protocol.md)的全部科学配置和五条件，不增加训练、不恢复A33、不开放DEV、独立物种或Setaria，不公开发布。拟合13406144仍在进行；工程完成和TRAIN继续条件分别判断。

## 固定比较

只接受实际13406144完成、4608更新对应的末step_004608.pt；001536／003072仅核对实际保存及optimizer曝光，不读其效果、不选checkpoint。原1536 draw对应1491唯一TRAIN窗，预测只含原A3已完成greedy＋beam4并集的16321条自然候选，顺序、几何和decoder_sources不变。训练中的536参考注入不进入预测。原自然池中不存在的链不会因重评分而出现，本轮不能直接声称恢复遗漏或自由生成。

基线是已实际闭合的A31_A3，原读出13398521与CPU核验13398524。新模型仅从原GEN缓存和DNA前向，所有权重冻结、无梯度、无optimizer、无新GEN前向、无新搜索。FP32发射，原15类角色／相位／start-stop优先级，单位FP64 role-logit减intergenic-logit累积；null0、原span非重叠DP，阈值0，不扫描、归一化或加权。未知候选正常竞争，不用标签过滤预测。

先把所有1491窗自然分数、选择和valid-base FP32发射矩阵保存并关闭分数文件，再接入原A31已经独立验证的自然候选1／0／-1标签。因为几何身份完全不变，不进行新R4标注，不读新catalog或参考序列来评分。原2916 eligible reference/window配对分母不变。标签0仅是primary相对已知非匹配，-1另列，不将其解释成生物学假阳性或新基因。

五个继续条件保持首次拟合协议：Arabidopsis和Oryza自然已选正确配对各严格增加；两物种自然已选known_nonmatch各不增加；pooled自然已选至少5-exon正确配对严格增加。比较均对A31_A3，不退回较弱A6、greedy或不同池。任一FAIL即关闭本模型，不加epoch、不换checkpoint、不改loss／权重／阈值；全部PASS也只支持另行准备raw-FASTA DEV方案，不等于SOTA、机制因果、独立泛化或投稿准备完成。

## 本轮实施与验证

先通过CPU短作业运行四项针对新实施的检查：保存FP32发射的完整窗口log-probability差与原score一致并拒绝错误分数；未知冲突链保留竞争且不记入known_nonmatch；拒绝重排／新增参考几何和缺失配对；五条件严格使用较强基线及分别物种。再核对原1491窗的完整身份、feature-index覆盖、原标签／分母／基线计数，仅读旧metadata，不打开新fit、feature tensor或运行模型。此准备能发现字段错配和评价口径错误，不能证明新模型收益；既有DP穷举和已闭合梯度检查不重复铺设。

GPU读出保存全部1491个发射矩阵（每窗最多24576×15×4字节，合计上限约2.05GiB，连同记录及格式余量计划低于2.3GiB）和natural_scores.jsonl／windows.jsonl／summary.json／STATUS／JOBID。没有新特征缓存或新checkpoint。提交前实际BeeGFS剩余156TiB，足够当前声明产物。

CPU结果核验只读实际产物，不进行新神经网络或backbone前向：逐条核对原1536 draw×3的顺序、importance权重、候选counts与4608真实有限更新记录；三个checkpoint的结构、有限FP32 tensor和AdamW参数／step曝光；全部1491窗保存发射的身份、形状和完整性；独立通过完整窗口FP64 log-softmax差重算16321分数（rtol1e-12／atol1e-8），不走A33 scalar-GRU replay、不放宽容差；再用原DP和标签重算已选counts、原较强基线、分母和全部五条件，与GPU记录逐项相等。此处score数学实现独立，角色映射和选择DP沿用原实现，不冒充独立算法或独立同行评阅。

## 执行资源与停止

三个入口均generanno、Slurm、no-requeue，沿用原排除节点和outputs/M28-PILOT-R5/logs唯一%x_%j日志。短实施检查private／0GPU／2CPU／8GiB／5分钟；自然读出private／1RTX3090／2CPU／16GiB／40分钟；CPU实际结果核验private／0GPU／2CPU／8GiB／30分钟。GPU读出不运行原A31的两个模型和逐窗greedy查询；原A31两模型＋查询完整读出实际1:00:41，原A34单模型396窗读出6:35，故单模型1491窗按约24:48加运行与保存余量设40分钟；这是操作估计，不是匹配速度benchmark。原A31完整CPU核验实际19:28，本轮无需重新构建R4标签／查询，保留30分钟上限。实际运行超过上限按TIMEOUT记录，不自动延长。

短检查实际COMPLETED0:0、四测试PASS及完整metadata核对后才能提交自然读出；读出afterok绑定拟合13406144与短检查真实job ID，CPU核验afterok绑定自然读出真实ID并显式传入。源码也要求实际匹配STATUS／JOBID／summary，不读取未闭合训练权重或部分结果。

任何工程失败保留原输出、错误与Slurm状态；本批次不自动重提、取消、延期或覆盖，主对话先定位具体原因并依既有授权和失败政策决定后续。科学FAIL不是工程失败。正常外部等待仅复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3和heartbeat coding-rna-r5-light-monitor，完整拟合→自然读出→CPU核验闭合或新可行动故障才一次交接；不以fit完成替代结果分析，不读取运行中指标，不重复主轮询。

## 实际执行记录

新源码、四项检查、三SBATCH和本协议共七个任务文件已明确传输至canonical；只传本任务文件，没有全项目同步或Git提交。三个Python AST与三个SBATCH语法均通过，旧数据／输出及无关修改保留。

短检查13408367经sbatch --parsable sbatch/M28-A36-STRUCTURED-EMISSION-READOUT-CHECKS.sbatch提交一次，初始PENDING调度等待，未因此重提或取消。实际root／batch／extern均COMPLETED0:0／0:57，batch MaxRSS743872K，四项检查PASS／0.311秒。完整metadata记录确认1536 draw／1491窗／16321自然候选／1491 feature-index，原身份、标签、较强基线和2916参考配对分母一致；新fit checkpoint／feature tensor未开、模型或backbone前向false、optimizer0、无新R4／DEV／Setaria／A33。此记录不包含新模型性能，不能称科学PASS。

自然读出13408449与CPU结果核验13408451已分别提交一次，命令如下，原拟合13406144没有重提、改参数或读取中间checkpoint效果：

```bash
sbatch --parsable --dependency=afterok:13406144:13408367 --export=ALL,A36_CHECKS_JOB=13408367 sbatch/M28-A36-STRUCTURED-EMISSION-TRAIN-READOUT.sbatch
sbatch --parsable --dependency=afterok:13408449 --export=ALL,A36_READOUT_JOB=13408449 sbatch/M28-A36-STRUCTURED-EMISSION-RESULT-CHECK.sbatch
```

实际提交前再次完成五项活调度／预约检查，private可用RTX3090容量、无相关维护窗口、实际存储156TiB及两独立输出无冲突；不触动其他作业。主即时squeue显示原fit13406144仍RUNNING／gpu034／37:03，两个新root为PENDING／Dependency。scontrol确认GPU读出ReqTRES1RTX3090／2CPU／16G／40分钟、CPU核验0GPU／2CPU／8G／30分钟，均no-requeue／Restarts0／generanno；CPU依赖实际13408449，GPU原SubmitLine含两个依赖，Slurm自动去掉已满足的13408367，仅保留未完成fit。这里是请求资源及排队证据，不是后两项实际allocation或成功结果。

同一固定监控心跳已原生更新为当前整条拟合→自然读出→CPU核验，保存的prompt／kind heartbeat／ACTIVE／固定target_thread_id／15分钟／failed_runs_only均准确读回，没有新线程或改模型。发出一次同线程manual startup核对请求，model／thinking省略以保留原设置；该启动检查、当前全链scheduled运行与实际终态交接分别取证，配置成功不代替执行成功。自然池效应、实际1491发射产物、CPU核对和五条件结论尚未产生。

manual startup现已实际完成：同监控turn01a1114a-0e82-76e2-94d6-9147e4a9b66c正常结束／error-null／56.644秒，真实12:57:26.485Z调度证据为fit RUNNING0:0／gpu034，natural及CPU PENDING0:0／有效Dependency；无产物或性能读取、无健康状态交接。实际turn_context12:57:08.517Z确认原gpt-5.6-luna／low和工作目录不变；不是scheduled验证。原fit-only在12:38的scheduled检查是真实历史证据，但不继承为新全链验证。当时全链定时执行、终态回传及科学结果均未发生。

当前全链scheduled执行现已取证：同固定对话turn01a11153-be4d-7253-bb7a-99624544e381正常结束／error-null／53.226秒，实际13:08:01.147Z状态为fit RUNNING0:0／gpu034，natural及CPU PENDING0:0／有效Dependency，未读取产物、未向主发消息。实际turn_context13:07:43.394Z仍是gpt-5.6-luna／low及原工作目录，定时运行确实追加到原监控对话；与先前manual启动分开验证。心跳继续ACTIVE，终态交接和科学结果尚未发生。
