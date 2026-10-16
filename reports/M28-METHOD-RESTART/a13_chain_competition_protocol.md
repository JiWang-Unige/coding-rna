# A13 固定完整链／null／冲突竞争目标

2026-10-03，任何拟合前固定。A12原TRAIN末权重诊断已COMPLETED/0:0：2916个eligible正确链／窗口几何对仅1237个gain正，5+ exon正确链仅155/1086为正；自然正确链1128/2380为正，known非匹配311/12566为正。正链在有已知span冲突的1974个配对中有463个不高于最大非匹配分。失败并非只在DEV，支持测试新训练目标；这些数值不是独立基因P/R/F1或泛化证据。

## 唯一变化与固定曝光

原A3 step004608生成器及GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345全部冻结。使用原A6自然free_pool（1491窗／16321候选）与A12记录的相同eligible正链注入、1/0/-1标签、候选几何和顺序；不重新搜索，不按参考生成自然池。读原1536draw顺序／importance_weight_uniform_within_stratum，3遍／4608 head更新。评分头仍为ChainHead H128／165635可训练参数、link冻结；同A6 seed0重新初始化，不继续拟合A6权重。FP32、AdamW lr3e-4、wd.01、betas(.9,.999)、eps1e-8、clip1均不变。沿用原TRAIN特征，无新GEN原生前向、DEV或test/Setaria。

相对A6唯一修改的是完整链目标，包含局部竞争与正／负组质量两项设计，不称其中单一组件的因果消融。自然候选来源、头架构、初始化配方和优化曝光与A6匹配；新的loss/正负组分配并不等于相同优化梯度或相同运行成本。原BCE在重建原draw上正类加权质量15.6234%，但末权重的加权BCE gain导数为正类−119.9028、负类+65.6092；比例不等于梯度占比，不能宣称负例支配当前更新、已收敛或BCE不可能有效。

## 固定目标

记每条链相对window null的gain为g。对每条已知正确链p，N(p)仅包含同窗、实际span严格重叠的已知标签0候选；触碰端点不算重叠，其他正链与-1未知／保护链不入竞争。正项是对p分别计算 `logsumexp(0, g_p, {g_n:n∈N(p)}) - g_p`，再对本窗所有正链均值。0代表真实null竞争，温度1、margin0，不调参；没有冲突时退化为softplus(−g_p)，仍要求绝对正gain。

负项为本窗所有已知标签0的softplus(g_n)均值，包含不与正链冲突的背景候选。若正、负两组均存在，取两项等质量均值；仅一组存在就该组单位质量；全未知／空池是包含null的可微零，仍按原固定draw计数执行AdamW，因此不声称每一步有标签监督。之后乘原draw权重。不是全窗一个正确类的softmax，不竞争多个合法共存基因，不把标签0当生物学真负。绝对负项是为避免“所有链一起升分”充当改善。

完整链判别式竞争及局部候选对比均已有原始先例，见[相关工作边界](../../docs/M28_complete_chain_competition_prior_work_20261003.md)。本轮是机制试验，不声称损失函数首创、精确全局最优或解决生物学假阳性。

## 检查、执行和冻结继续条件

先运行五项针对多正链独立、未知与真实重叠、竞争梯度、null绝对锚及可微零的检查；再直接核验A12实际1491条记录的来源／计数／分量恒等式／冲突最大值与三类主profile，重建原draw计数。从原1491窗按metadata选择八个新梯度哨兵：每物种在同时有已知正／负的窗口中依次选最多正链、最多正链exon、最多实际span冲突、首个正链且有-1保护的窗口，不重复选择，同值按原唯一窗顺序；完全不使用gain、loss或性能选择。一次真实反向验证五个使用中评分分支有限非零梯度、生成器/link无梯度、权重不变；不重跑无关测试。

预提交metadata检查发现旧FULL readiness八个哨兵有七个不在原1491 draw窗内，故本轮不沿用它们；上面的固定新梯度样本始终在原1491窗内，不扩充训练曝光。此修正在任何A13作业执行或新目标结果之前完成，不是按结果更改样本或gate。

preflight通过后正式拟合。独立outputs/M28-A13-CHAIN-COMPETITION-{PREFLIGHT,FIT,TRAIN-READOUT}，入口run_chain_competition.py；private／1 RTX3090／2CPU／16GiB，三段分别5／45／20分钟。原A6三遍head更新实际合计约23分钟，A12全1491窗评分421秒，所申请时窗留操作余量，不是新方法吞吐结果。约数MB新checkpoint和文本，不复制特征。保存head step001536／003072／004608，但评价只用末004608；不挑seed、epoch或阈值，不改变tau0、既有实验或停止条件。

末权重在原1491窗同一候选几何／来源／标签上零更新读出，与A12末A6配对。事前继续条件全部满足才准备单独raw-FASTA DEV协议：两物种各自自然正确链的正gain数和全部正确链的正gain数均严格增加；两物种各自known非匹配正gain数均不增加；pooled 5+ exon正确链正gain数严格增加；pooled正确链不高于最大已知冲突的数量不增加。所有分母及候选池保持相同，禁止删未知或长链。它是in-sample TRAIN准备gate，不是最终注释、背景span、F1或独立泛化。

若gate FAIL，保留负结果，不加epoch、扫margin／温度／权重／阈值或挑另一head救本轮。若PASS，下一阶段仍须另定完整8690 raw方向窗、7728全部主参考、原F1／背景门槛及同尺原生对标的正式DEV执行记录；本轮不提交DEV推理，不把TRAIN通过当论文PASS。FULL与A10原合同、固定等待对话和同一heartbeat不变，没有Git提交或对外发布。

## 预检终态与固定拟合提交（2026-10-03 15:44–15:48 UTC）

预检13360817的root／batch／extern均COMPLETED/0:0，27秒，STATUS COMPLETED／JOBID对应。五项单元检查PASS，0.369秒；八个原1491窗内按事前metadata规则选出的实际哨兵均通过五个评分分支的有限非零梯度检查，真实反向7.4169秒。生成器／link无梯度，整个head权重未变，optimizer0、无新搜索／原生GEN前向／DEV／test。不是训练或科学PASS。

实际loader逐条重建A12的1491条windows记录，核对16321自然候选＋536注入配对＝16857候选／窗口对、标签2916／12566／1375、来源／顺序／链几何、分量恒等式与每条正确链的实际已知span冲突最大分；重算正／非匹配／未知三个主profile并与A12收据完全一致。重放原1536draw得到自然16806／注入555、标签3028／12935／1398、无known标签窗3，与原A6每遍曝光相同。unique与重复draw分母分开，不能当独立基因数。

| 阶段 | Slurm ID／实际状态 | 资源及依赖 |
|---|---|---|
| 零更新preflight | 13360817／COMPLETED 0:0 | private、1RTX3090／2CPU／16GiB、5分钟 |
| 固定fit | 13360834／RUNNING gpu034 | 同资源、45分钟 |
| 末head TRAIN读出 | 13360836／PENDING有效Dependency | 同资源、20分钟；afterok:13360834 |

实际命令为 `sbatch --parsable sbatch/M28-A13-CHAIN-COMPETITION-PREFLIGHT.sbatch`、`sbatch --parsable sbatch/M28-A13-CHAIN-COMPETITION-FIT.sbatch`、`sbatch --parsable --dependency=afterok:13360834 sbatch/M28-A13-CHAIN-COMPETITION-TRAIN-READOUT.sbatch`。提交前按coding-rna-baobab执行五项live调度检查，核实实际存储及新输出不碰撞；private有可用3090，无切断本轮的GPU维护预约。源／单测／三份SBATCH／协议七文件的AST／bash语法及local/remote完整文本一致；本节仅追加已发生执行事实，不变更冻结目标。fit实际STATUS RUNNING／JOBID13360834对应，启动不是4608步完成。

同一监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor复用，gpt-5.6-luna／low／15分钟不变。原十根cohort于15:42 UTC自然触发，A10＋A12终态交接15:44:29.300 UTC被主对话接收，实际PAUSED已读回。随后更新为M28-FULL-A13-COMPETITION-20261003：保留FULL原五根，加A13 fit／readout两根；整份prompt及身份、固定target、ACTIVE、cadence、failed_runs_only、created_at精确读回。15:48:10.857 UTC人工启动检查在同一对话完成，七根运行／有效依赖正常；它不是新cohort的定时运行，尚未观察本cohort自然检查或终态自动回传。正常等待继续由原心跳负责，主对话不重复等待轮询。本机／应用／连接须在线，见[官方任务说明](https://learn.chatgpt.com/docs/automations?surface=app)。

A13目标／gate在首次作业和读取新A10指标前固定；不以A10后验结果改变本轮训练或候选。当前仅有预检PASS及真实拟合启动，TRAIN准备gate尚无结果，未提交A13 DEV。

## 终态关闭（2026-10-03）

上述“尚无结果”为提交时历史状态。fit13360834及末head TRAIN-readout13360836现均COMPLETED/0:0，实际22分07秒／6分34秒；4608条连续step和三遍原draw身份／权重／标签计数已核对。实际1491窗口／16857候选配对保持原几何／来源／1/0/-1分母，八项事前继续条件重算后六PASS、两FAIL。两物种known非匹配正gain由223→689／6305和88→329／6261增加，因此原TRAIN准备gate FAIL，未提交DEV，不改阈值／epoch／温度／margin／组质量救本轮。

完整事实、CPU夹具检查13360949的失败及最小修复13360981、实际记录核验、未执行DEV和固定监控接续见[结果记录](a13_chain_competition_result.md)。条件式raw源码仅属准备，不是原生DEV推理／7728参考或论文PASS；FULL原合同保持不变。
