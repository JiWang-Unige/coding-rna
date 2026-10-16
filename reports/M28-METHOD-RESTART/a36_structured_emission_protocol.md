# 整链监督的发射表示模型

2026-10-06，首次实际梯度检查与拟合前固定。本轮先实现、核对零更新梯度，再准备一次固定训练；不修改任何旧模型、失败门槛或停止条件。目标仍是用户已授权的raw-FASTA基因注释模型迭代，不将本文件存在视为额外授权，不据此放行独立物种、sealed或公开发布。

## 为什么更换学习目标

[A34](a34_missing_path_emission_result.md)全536个原遗漏配对即使有正确几何，仍有190个非正、179个正分未选；[A35](a35_saved_role_contributions_result.md)定位非正链的coding／intron负总贡献，同时有43个正分竞争失败没有任何负贡献组。故只调整原query搜索不充分，也不能凭这些统计降低intron权重救原分数。

新机制让表示和15类发射直接接受完整链集合监督：复用现有SharedContext的GEN投影、双向GRU、局部DNA卷积、融合与15类线性头，fresh seed0，feature_dim2048／H128／366191可训练参数。GEN本身仍冻结，仅用原四块6144bp缓存；原A3查询生成器不更新，也不把其已训练hidden当成新模型输入。当前新模型不含query、ChainHead或link头。

A14／A15曾在冻结A3表示上训练165635参数的候选ChainHead；本轮训练的是从原缓存和DNA学习的表示与发射，不是给原发射加六个权重。候选集合似然、结构学习和相位发射本身已有先例，见[既有方法边界](../../docs/M28_complete_chain_competition_prior_work_20261003.md)，不据实现差异宣称论文创新。

## 固定科学选择

完整原1536 TRAIN draw／1491唯一窗、六条原TRAIN染色体、原顺序及importance权重不变。监督使用已闭合原A12全部16857候选／窗配对：16321自然、536训练参考注入；标签1／0／-1分别2916／12566／1375。注入仅用于训练目标，不进入自然预测或继续门槛；标签0仍是primary相对已知非匹配，不是生物学真负，-1保护。

原15类角色、相位、start／stop覆盖优先级与A31单位role logit−intergenic logit和不变，FP32网络、FP64累计、null0。新可微函数须逐项复现原固定函数，并证明梯度等于完整窗口log-probability差的梯度。没有长度归一化、额外转移先验、温度、margin、阈值或成分权重。

复用原chain_set_loss：K仅标签0／1，G为全部标签1，loss＝(logZ−gold gain和)／|K|；logZ覆盖K的全部span兼容集合及空集，相触可共存。未知不进入分区和gold，自然选择仍保留未知。空／全未知为可微零，不删除draw；gold不兼容即停止，不删正链。本模型不加原局部CE、query或其他辅助损失；目的是检验整链监督能否同时改善完整链对null和实际竞争结构的相对分数，不把局部数学分解误称因果机制。

一次fit计划为fresh seed0、原1536 draw×3＝4608 AdamW更新，lr3e-4／wd.01／betas(.9,.999)／eps1e-8／clip1，保存001536／003072／004608，主读出只用末004608。不从旧监督checkpoint续训，不挑最佳epoch，不增seed或曝光。生成器既有4608训练与新表示4608训练属于额外模型曝光，参数量、目标和监督亦不同，不能宣称总计算／曝光匹配的因果消融。

本轮新表示只在原冻结自然池评分，不能生成原池中缺失的链；536注入目标可改善表示监督，但不存在候选的正确链不可能仅靠重评分被恢复。成功也须另做原A3生成器＋新表示的真实raw-FASTA协议与自然搜索机制，不能把训练参考注入变成推理依赖。

## 零更新准备

入口preflight_structured_emission.py，输出outputs/M28-A36-STRUCTURED-EMISSION-PREFLIGHT。先复用原A14全1491窗gold兼容／原A12身份核对，再按原metadata规则选两物种各四个梯度窗口，不看新gain／loss挑样例。实际检查project、context、dna、fuse、segmentation五分支有限非零梯度，输入缓存／DNA无梯度，所有state前后逐tensor相等；不构造optimizer，不读监督checkpoint，不写checkpoint。

四个针对新变化的检查覆盖可微score／梯度与完整log-probability差一致、相容多基因目标对应独立known BCE而非单基因softmax、未知冲突链不改变已知loss／梯度、空／全未知保留零梯度。已有集合DP穷举验证不重复铺设。全部模型和重数据运行通过Slurm；private1RTX3090／2CPU／16GiB／10分钟，generanno，原排除节点名单、唯一日志与输出。只有小型记录，不写新缓存或logits。

本preflight不能证明拟合、分数改善、自然注释或泛化。通过后按实际load＋forward＋backward时间制定一次fit的运行上限与读出资源，保持全部科学配置；短暂失败保留输出，至多一次同范围、同资源、同科学选择的明确工程修复使用独立R1目录，不隐瞒原FAILED。不能把失败科学结果当工程重提。

## 拟合后继续条件

末模型在同1491窗全部16321自然候选上零更新评分，只对自然池使用原span非重叠DP、null0、无标签过滤／阈值扫描。基线为实际已闭合并独立核验的A31_A3同自然池、同原R4标签，不退回较弱A6门槛。五条件：两物种自然已选正确配对分别严格增加、两物种自然已选known_nonmatch分别不增加、pooled自然已选至少5 exon正确配对严格增加。原未知选择另列，不归假阳性或新基因。

FAIL关闭本模型，保留原结果，不加epoch、改loss权重／score单位／阈值或另挑checkpoint。PASS仅支持准备下一完整DEV方案，不等于独立适用范围、SOTA或论文成功。DEV已用于原研发，独立面板上游曝光、Setaria封存及公开发布边界不变；A33原FAILED与一次工程例外待人类答复，本轮不读其输出、不覆盖验证器、不重提失败分片。

需要长期外部等待时只复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3与同heartbeat coding-rna-r5-light-monitor，不创建新线程、cron或goal。工程准备、一次fit、TRAIN继续条件、完整raw-FASTA比较、独立泛化和论文主张分别报告。

## 实际执行记录

2026-10-06提交零更新preflight13404309，请求private／1RTX3090／2CPU／16GiB／10分钟，命令sbatch --parsable sbatch/M28-A36-STRUCTURED-EMISSION-PREFLIGHT.sbatch。提交前核对实际集群配置、五项调度／预约信息、源文件和唯一输出无冲突；原8个用户作业未触动。环境沿用本实验系列实际已成功运行的generanno，而非旧框架通用默认coding-rna。仅五个新增任务文件明确传输，canonical三Python AST与SBATCH语法通过。首次squeue证据为RUNNING／0:52／gpu034；这不是检查通过、拟合或科学结果。随后两次ProxyJump连接在banner阶段超时，暂时不能据此核实终态或实际allocation，不把连接故障当作作业失败，不重提或修改SSH配置。

连接中断期间本地新增fit_structured_emission.py，唯一目标outputs/M28-A36-STRUCTURED-EMISSION-FIT。入口只接受13404309实际COMPLETED／0:0和对应完整零更新记录，固定原4608曝光、三个checkpoint与末权重；没有自动提交、重提、挑checkpoint或DEV分支。当时该新入口仅完成本地语法核对，尚未上传、未运行，也没有提交fit；恢复后仍须先取证preflight，再核对运行资源和相称实施检查。当时唯一在集群提交的新对象仍是13404309。

随后原SSH路径恢复，主对话实际sacct核实13404309 root／batch／extern均COMPLETED／0:0、0:58，batch MaxRSS1766132K；AllocTRES确为1RTX3090／2CPU／16G、gpu034／private／10分钟。scontrol已返回Invalid job id，未据此判失败；以sacct和匹配STATUS／JOBID／summary闭合。四项实际检查PASS／1.284秒，八窗五分支梯度全部有限非零、原A31分数误差均0、366191参数前后相同、optimizer0。GPU峰值163786240字节，内部八窗19.853秒；第一窗10.129秒包含冷启动，之后七窗0.645至2.768秒、均值约1.386秒。这些是工程时序，不是速度benchmark或准确率结果。

同一监控仅完成一次manual终态交接并暂停，未声称scheduled验证；主消费证据根始终为/home/users/j/jwang/coding-rna。实际全过程无需重提，SSH配置未改。

按用户既有同范围迭代授权推进一次上述固定fit，SBATCH请求private／1RTX3090／2CPU／16GiB／2:30，generanno；原训练步数、数据与科学配置不变。操作时限参考七个热启动窗口的约1:46模型外推，加约40%运行余量；八窗并非随机计时总体，不能据此承诺实际时长或与原A3拟合1860.993秒作匹配速度比较。新增训练日志／三checkpoint预计低于64MiB，未增加特征缓存或logits。首次fit失败保留全部状态；仅明确证实尚无optimizer更新的工程失败允许至多一次同科学／同资源最小修复、独立R1目录，不能隐瞒原失败或盲目重提。已有更新或TIMEOUT不自动重提、延长、续训或解释未闭合结果；科学FAIL按原五条件关闭。拟合完成仍须新自然池读出和相应实际结果核对，不能把fit终态替代研究成功。

正式fit13406144已提交，命令sbatch --parsable sbatch/M28-A36-STRUCTURED-EMISSION-FIT.sbatch。期间一次单文件SCP因banner超时失败，随后只补传缺失任务文件；canonical fit Python AST及SBATCH语法均通过、实际13404309准备闭合、唯一fit目录无冲突后才提交一次，没有重复Slurm作业。运行终态、4608真实更新、三checkpoint和自然池读出仍待完成；不能称为模型收益或论文结果。长等待仅将同一固定心跳／对话改为13406144，不监控已消费preflight或旧失败依赖。
