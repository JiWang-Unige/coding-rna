# 联合学习上下文与完整剪接链评分

2026-10-06，首次实际梯度与拟合前固定。A37独立检验fresh上下文与有序exon／junction整链头联合学习，能否改善原零操作点下完整链对空集及冲突结构的选择。A36已科学FAIL关闭；本轮不续训它、不改其局部发射loss、分数、阈值或停止条件。目标仍为已授权的raw-FASTA基因注释模型迭代，不放行独立物种、Setaria或公开发布。

## 实际依据与结构区别

[A36实际转移](a36_structured_emission_result.md)中，原已选正确配对共丢失714条，其中665条新单位gain≤0、49条为正分未选，另新增47条。这个阶段定位支持研究链对空集的学习机制，但不证明局部发射是唯一原因或非加性头一定有效。

既有[A21](a21_junction_chain_result.md)已经有有序exon／junction GRU及内含子内容，不能说这些表示缺失或把同一头重新命名为创新。A21仅学习头，输入为冻结A3上下文；A36学习fresh上下文及15类发射，但完整链分数是固定单位局部累积。A37让fresh SharedContext的GEN投影、双向GRU、DNA卷积及融合同时接受既有JunctionChainHead的整链梯度，不包含query、link或segmentation头。

模型612324参数＝上下文364256＋整链头248068，feature_dim2048／H128。fresh seed0的一次随机流先构造SharedContext再构造头，不继承A3、A21或A36训练权重；头初始tensor因此不冒称与单独seed0的A21头相同。gain为exon additive＋junction additive＋有序链nonadditive−窗口null，FP32。空集相对gain为0，操作点固定0；不是A36的单位角色分数，也不做尺度匹配／阈值扫描。复用现有有序头的方向化半开池化、累计相位、相触连接及valid-length语义，不改原头实现。

这是组合模型假说，不是容量或总训练曝光匹配的因果消融，也不是新结构目标。可微结构层、神经表示联合学习已有[直接先例](../../../paper/manuscript/M28_complete_chain_competition_prior_work_20261003.md)；只有后续真实完整链收益、机制控制和独立适用范围才能支持论文贡献。

## 固定数据目标与曝光

原1536 draw、1491唯一TRAIN窗、六条原TRAIN染色体、顺序、importance权重和GEN revision不变。GEN冻结，只读原四块6144bp特征和DNA，不做新GEN前向。候选仍是原A3末004608生成器的16321自然greedy＋beam4链，几何与来源不变，没有新搜索。本轮监督只使用自然域，标签1／0／-1为2380／12566／1375；不使用536 TRAIN参考注入。这与A21的自然域相同，但与A36的含注入训练域不同，不能把两者差异单独归因于评分结构。

复用原chain_set_loss，K只含原0／1，G为全部自然标签1兼容集合，loss＝(logZ−gold gain和)／|K|，分区包括所有span兼容known集合与空集。-1不入分区及gold，直接gain梯度为0，共享参数更新仍可能改变未知分数；预测保留所有自然未知竞争。没有local CE、类别权重、辅助损失、margin、温度或新标签。去注入后原72唯一窗失去所有注入gold，空／全未知和零gold draw仍按原目标与步数保留，不删除或补正链。

通过零更新核验后的一次候选拟合固定fresh seed0、1536×3＝4608 AdamW更新，lr3e-4／wd.01／betas(.9,.999)／eps1e-8／clip1。保存001536／003072／004608，科学读出只允许末004608，不选最佳epoch、不追加seed或曝光。具体fit上限依据本次实际load＋forward＋backward测量另行记录，不能用小窗口时序作方法速度排名。

## 零更新核验与失败处理

四项针对新接线的人工检查覆盖上下文→既有头的实际组合、候选顺序不改变单链分数、valid-padding不进入评分、整链loss到上下文／DNA及头分支的梯度、未知加入后known loss／联合梯度不变、空／全未知为零。不重写或重复穷举既有集合DP。

实际preflight按原metadata规则选两物种各四窗，不看fresh分数挑样。核对原自然身份、原标签与较强A31_A3及2916 eligible reference/window分母一致；自然至少5-exon池不是1086参考配对分母。四个上下文及七个头分支须有限非零梯度；没有known多exon候选的窗口，其junction两分支须无梯度或严格零，不因这一真实适用条件替换窗口。两物种各须至少一个known连接覆盖，未知gain直接梯度严格零，原GEN缓存／DNA无梯度，所有新参数前后逐tensor一致，无optimizer或checkpoint。

入口preflight_context_chain.py，唯一输出outputs/M28-A37-CONTEXT-CHAIN-PREFLIGHT；Slurm generanno／private／1RTX3090／2CPU／16GiB／10分钟／no-requeue，原排除节点及唯一日志。只写小型核验记录，不建隐藏缓存或保存新logits。短暂工程失败保留原产物；仅在证实optimizer0时允许至多一次同科学、同资源的明确实现修复至独立R1目录，不能盲目重提、删状态或将科学FAIL当工程失败。

## 固定继续条件

末模型在全部1491窗、16321自然候选上零更新读出，原span非重叠DP及gain0操作点，不按标签过滤或加入参考。较强同池基线仍为已核验A31_A3，不退回A6，也不重判旧A14／A15／A21／A36 FAIL。五条件：两物种已选正确配对各严格增加、两物种已选known_nonmatch各不增加、pooled至少5-exon已选正确配对严格增加。

原基线正确1531／478、known_nonmatch347／251，pooled至少5-exon正确700；eligible reference/window2172／744和至少5-exon854／232不变。unknown另列，不视为生物学假阳性或新基因，known_nonmatch也仅为primary参考相对标签。候选／窗口配对不是独立基因或完整染色体F1。

任一FAIL即关闭A37，不调loss／权重／阈值、不改checkpoint／曝光；PASS只支持另定真实raw-FASTA方案，不能代表SOTA、独立泛化或投稿成功。原自然池缺失几何不能由本重评分模型恢复。需要长期外部等待时仅复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同heartbeat coding-rna-r5-light-monitor；A33例外、DEV／独立面板／Setaria／公开发布边界不变。

## 实际核验与拟合资源

零更新preflight13409938实际root／batch／extern均COMPLETED0:0／0:57，private／1RTX3090／2CPU／16G／gpu034，batch MaxRSS1858748K，匹配STATUS／JOBID／summary。四项检查PASS／1.548秒；原1491窗自然身份／标签／较强基线及2916／1086分母逐项相符，原1536 draw每遍自然支持16806、标签2473／12935／1398、zero-known draw3，唯一窗zero-known2，去注入失去gold的窗72。八个原metadata窗口均未替换；七窗known junction适用，另一水稻窗仅单exon，连接两分支确实为零，其余九分支有限非零，两物种均有known连接梯度覆盖。未知直接gain梯度零、参数前后相同、optimizer0；612324参数和结构与协议一致，没有监督checkpoint、生成器／GEN前向、新搜索／R4／DEV／Setaria／A33。

实际八窗内部load／forward／backward共约19.58秒；第一窗8.171秒含冷启动，其余七窗0.713至3.374秒、均值约1.630秒。内部整段19.614秒，GPU峰值165758464字节。按热窗口均值外推4608步约2:05，再留运行与I/O余量，本次fit上限3小时；这不是随机总体计时，也不是方法速度比较或保证完成时间。训练入口fit_context_chain.py先要求真实13409938闭合及匹配零更新收据，使用同fresh构造和原全部draw。

fit请求generanno／private／1RTX3090／2CPU／16GiB／3:00／no-requeue，三checkpoint连optimizer与日志预计低于64MiB，不建新缓存／logits。记录每个参数实际AdamW step：没有候选／连接的窗口可能不调用相应头分支，不能一概把其参数曝光假定为4608；后续CPU据原候选调用结构重算。原零loss draw仍执行既定optimizer过程，不因unknown保护跳过更新。

首次fit失败保留全部状态；仅明确证实未发生任何optimizer更新的工程失败允许至多一次同科学、同资源最小修复到独立R1目录。已有更新、TIMEOUT或不明更新状态均不自动重提／续训／延长；科学FAIL不属于工程修复。fit终态后仍须全部自然候选读出和实际结果核验，不能以训练完成替代五条件、完整注释或论文成功。

实际fit13410061已提交一次，命令sbatch --parsable sbatch/M28-A37-CONTEXT-CHAIN-FIT.sbatch；canonical新源码AST／SBATCH语法通过，原零更新核验闭合、唯一输出无冲突后才提交。再次完成五项活调度／预约检查，原五个用户作业低于80，没有相关维护，实际存储156TiB足够声明产物。首次squeue／scontrol记录RUNNING／1:44／gpu034，实际allocation为上述1RTX3090／2CPU／16G／3小时，no-requeue／Restarts0，canonical batch command与SubmitLine匹配。此时真实4608更新、三checkpoint、自然读出和科学门槛仍未产生，不把运行状态当作结果。

同一固定心跳已原生更新到当前fit13410061，实际保存的完整prompt、heartbeat类型、ACTIVE、原固定对话、15分钟、failed_runs_only逐项准确读回，created_at保留，没有新监控对话／自动化或模型改动。启动核对、实际定时追加与终态自动交接须各自取证，不继承旧A36阶段的成功声明。主负责下游实施，监控仅读最低调度和终态结构，不读运行中性能。

manual startup已实际完成：同固定对话turn01a111ae-b4b3-79a1-8679-f724be685ab0正常结束／error-null／26.220秒，真实UTC秒级14:47:23调度输出为13410061 RUNNING／gpu034、sacct RUNNING0:0，无性能读取或健康状态主交接。实际turn_context14:47:04.680Z确认仍为gpt-5.6-luna／low及原目录；这是新cohort的手动检查，不是定时检查。定时运行与训练终态／自然读出／科学结论仍待取证。
