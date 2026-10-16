# 原 TRAIN 粗粒度观察的信息梯度诊断

2026-10-11。已有部分观察目标不能直接证明粗粒度标签改善完整链预测。原fresh准备13453348中，113,246,208次核心碱基曝光包括111,308,039 fine、524,144 coarse-only和1,414,025 unknown；各占98.288535%、0.462836%和1.248629%。这是重复draw的曝光量，不是独立碱基、基因或样本数。coarse-only虽少，也可能位于重要结构位置，不能只凭比例认定无作用。

## 固定问题与接入

固定使用已完成JOINT24576末9216头（fit13453777），在原seed1 fresh4608计划中选择每物种首个coarse_positions大于0的draw，共两个原TRAIN视图。选择规则在运行前固定，不按梯度、loss、预测结果或DEV选择窗口。复用原fine可行性支持13453408及原联合24576生产编码、nearest-center ownership、有限视图、cap128和FP64分区。参考或标签不输入encoder／head，只有其后观察目标读取标签。

同一份native logits分别计算epsilon0.1的原观察NLL、epsilon2/3的均匀coarse通道NLL，以及仅保留相同fine观察的NLL。第三项仅是计算中的目标副本，原标签、准备计划和冻结模型不改。三者均除以原max(1,fine+coarse)行数，不能因删除coarse信息同时改变归一化尺度。

均匀通道的三个region概率均为1/3，其NLL必须等于fine-only NLL加coarse行数乘log3，并保持相同logit梯度。这是实际实现中绝对观察概率及background校正的核验，不是把均匀通道的常数损失误解为新的学习信号。原epsilon0.1与均匀通道的梯度差描述该固定头、该视图中粗粒度观察的条件信息作用；零或小差值不能外推所有4219个视图、所有参数值或长期优化结果。

## 执行和判定

只对同一固定势函数求logit导数，不对头或backbone求参数梯度，不创建optimizer、不更新权重、不运行decoder、不保存checkpoint／势函数／大特征。保存两视图的细／粗观察计数、三项raw及原尺度归一NLL、均匀通道常数校正误差、gradient最大差及五个互斥角色组的归一梯度L2。所有十五角色（包括background）均保留；这些量不是概率置信度、基因恢复率或生物因果效应。

均匀通道与fine-only常数或梯度关系失败时停止解释，修正具体实现；实际原观察梯度差无论大小都完整记录，不据此改现有训练或冻结DEV。两视图只回答固定点的局部问题，不能宣布粗粒度通道普遍必要或冗余，也不能替代正式性能对照。只有确有额外观察信号及需要检验其长期学习效果时，另固定TRAIN训练对照，不根据已读DEV决定模型或权重。

单24GiB RTX3090、2CPU／16GiB、20分钟，使用现有generanno环境；时限依据原两物种native预检约4分钟，加三套分区和导数的余量。全部科研读取及native计算经Slurm。独立根outputs/M28-FULL-GRAMMAR-COARSE-INFORMATION-PROBE，产物上限256KiB，不覆盖任何原输出。原fresh×joint正式训练13494895及下游作业不改、继续由原唯一心跳监控；本短诊断由主在终态后消费，不创建监控对话或自动化。

本方案在完整DEV结果已读取后固定，只使用原六条TRAIN；不称DEV前预注册，不改原八方法及十方法家族，不访问held-out、sealed、Setaria或新物种。原论文目标仍未完成。

## 实际提交与启动

新增脚本AST和SBATCH语法通过，三个新文件与远程按固定顺序直接逐字节相等，旧训练／推理脚本不改。已有人工数学检查13429599的实际summary为completed、五项通过、零failure／error，包含epsilon0.1和2/3；没有重复运行这些不变检查。原fresh准备收据13453348的两物种fine、coarse和unknown曝光分别重构各56,623,104总行，pooled也完整守恒。

实时Slurm显示private两节点各有空闲24GiB 3090，private与shared均无影响本20分钟运行的预约，账户作业数低于80；新根、脚本和同名作业无碰撞，BeeGFS可用147TiB。单次sbatch返回13495642，UTC2026-10-11 15:30:46提交、15:30:47开始于gpu034；实际AllocTRES为单3090／2CPU／16GiB／20分钟、无dependency、Requeue0及Restarts0，%x_%j日志独立。

启动收据确认原joint末9216权重和本固定配置，选择的是fresh计划additional_step9（拟南芥）及7（水稻），均为各物种首个有coarse-only观察的draw，而不是原预检中不含coarse的首个draw。零optimizer更新、DEV／sealed=false；启动时stderr空。原fresh×joint正式训练13494895未改。

## 终态与局部结果

13495642的root、batch、extern均COMPLETED／0:0；UTC15:32:46结束，调度用时1分59秒，脚本用时99.861秒。STATUS、JOBID、started和summary齐全且非空，总计9,344字节，stderr为0字节。实际收据为source fit13453777／9216，native forward=true、optimizer未创建、零更新、头和encoder不变，原目标不变，未运行decoder、未改预测或保存checkpoint，DEV／sealed=false。此短诊断已消费，不加入心跳。

拟南芥step9的core为NC_003076.8:6721536–6746112负链，fine24,479、coarse1，分母24,480；informative、uniform、fine-only的raw NLL分别3.6163532979612683、4.609605070971156、3.5109927823032194。水稻step7的core为NC_089046.1:23248896–23273472正链，fine24,354、coarse3，分母24,357；对应raw NLL为0.3200948958562435、3.2998502148870443、0.004013348882722312。两者native context均49,152bp，坐标为0-based半开区间。

独立从summary重算原分母、全部归一NLL和五组L2平方和：归一NLL重构误差为0，角色组L2重构最大误差2.17e-18。uniform−fine-only与n_coarse×log3的误差分别−1.73e-13、−7.11e-15，原始logit梯度最大差均0；informative−uniform的归一梯度L2、最大绝对值及五组值也全部0。informative−fine-only还分别与−n_coarse×log0.9相符，误差2.23e-13、4.24e-14。分区使用FP64，求导结果返回native FP32 logits的精度；0是本次实际返回的数值，不是对任意参数或观察集合的数学证明。

结论仅限这两个固定点：粗标签在这两个视图没有提供额外可测logit学习信号，其NLL差表现为常数。这与fine约束已条件确定相关region相容，也可能受当前条件分布集中和返回精度影响；本次没有证明是哪一种原因。因为两个视图合计仅4个coarse位置，不能代表原计划524,144次coarse曝光，不据此启动均匀coarse正式训练、删除观察通道、改论文性能claim或变更DEV家族。

## 后续：每物种最大coarse曝光视图

在上述结果已知后，固定追加每物种原seed1计划中coarse_positions最大的视图；重复draw或计数相同时选最早additional_step。不按loss、梯度、预测恢复率或DEV选取，仍非独立验证或DEV前预注册。轻量计划元数据确定拟南芥step4106为NC_003076.8:9486336–9510912正链（fine17,449、coarse7,127），水稻step4048为NC_089043.1:14401536–14426112负链（fine7,102、coarse17,474）；unknown均0，原分母均24,576。全计划有coarse的draw分别1,301、245个，不能用两个最大值视图替代其总体。

复用同一个脚本，仅增加显式--selection max-coarse路径；默认first选择及其全部计算保持不变，原13495642收据不改、不重跑。新路径仅改变视图选择和独立输出根outputs/M28-FULL-GRAMMAR-MAX-COARSE-INFORMATION-PROBE；头仍为fit13453777的9216，三项NLL、原行数归一、logit导数、uniform恒等关系及禁止更新／decoder／DEV／sealed的边界相同。新增SBATCH仅改job name、输出根及这一参数，继续单3090／2CPU／16GiB／20分钟、产物上限256KiB。

这次检查专门区分“前两点稀少coarse导致未见额外信号”与“更多coarse的固定点也无信号”。若非零，才说明至少这个原TRAIN固定点存在额外logit信号，不说明性能提高、普遍必要或因果机制；若零，仍不证明所有窗口、参数或轨迹均冗余。两种结果都不自动触发正式训练或冻结DEV变更。

新选择路径AST及本地／远程SBATCH语法通过，修改脚本、新SBATCH和本记录在提交前逐字节一致；默认first分支不重跑。新输出根与同名作业无碰撞，实时账户展开作业38个、低于80上限，private两节点各余6个24GiB 3090，无影响本20分钟运行的预约。实际sbatch为sbatch/M28-FULL-GRAMMAR-MAX-COARSE-INFORMATION-PROBE.sbatch，返回13495706，UTC15:41:00提交、15:41:01开始于gpu034；单3090／2CPU／16GiB，20分钟、无依赖、Requeue0／Restarts0。原13494895及其唯一心跳仍不变。

## 最大coarse视图的完整结果

13495706的root、batch、extern均COMPLETED／0:0，UTC15:44:47结束，调度用时3分46秒，脚本189.325秒。实际started和summary均选step4106／4048，source fit13453777／9216；原目标、权重和预测不变，optimizer未创建、零更新，未运行decoder或保存checkpoint，DEV／sealed=false。STATUS、JOBID、started和summary完整且非空，共9,787字节，stderr为0字节；此短诊断不加入心跳。

拟南芥三项raw NLL（informative／uniform／fine-only）为877.8269683386886／7829.822451148193／0.012669810577335738，归一后为0.03571887078201044／0.3185962911437253／5.155359121637263e-7。水稻对应为1900.8611374708503／19197.152100822375／0.0009686358375802229，归一后为0.07734623769005738／0.7811341186858063／3.941389313070568e-8。二者均保留原24,576行分母，不能将uniform的大常数NLL误作更差学习或性能。

独立重算uniform−fine-only−n_coarse×log3分别为−2.73e-12、−1.46e-11，原始logit梯度最大差均0。全部归一NLL重构误差0，五组L2平方和重构最大误差4.17e-17。informative−uniform的归一梯度L2分别0.0010135433970350615、0.005446409732048435；最大绝对值分别4.069007991347462e-5、4.0690104166666664e-5。拟南芥五组差值L2（background／coding／intron／splice／start-stop）为1.27698e-8／2.92114e-4／9.67696e-4／7.41936e-5／3.85049e-9；水稻为3.60515e-3／7.51939e-4／4.00940e-3／1.40675e-4／7.61913e-5。角色组是互斥列的梯度范数，不是因果贡献或错误比例。

informative−fine-only超过“所有coarse均匹配”的−n_coarse×log0.9常数分别126.90990343478336、59.790518230156295，与前两稀少coarse视图的纯常数行为不同。这两个最大coarse视图在当前固定头确有额外可测logit信号，故“原计划coarse普遍无作用”不成立；观察通道的必要性、平均作用、参数梯度和长期完整链性能仍未证明。无需为已有非零信号再增设结构可达性证明或重跑原数学检查。下一项应是另行固定的同起点TRAIN学习对照，而不是改原头、调角色权重或根据DEV选模型。
