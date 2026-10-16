# 原 TRAIN 的粗标签信息追加学习对照

2026-10-12 UTC。Uniform与informative最终9216原TRAIN整链均已完成并由主核验；[完整结果](full_grammar_uniform_coarse_result.md)保留逐染色体、物种、pooled、macro及worst。Informative减uniform的F1为拟南芥−0.908052 pp、水稻+1.536665 pp、pooled−0.342257 pp，不支持一致收益或普遍必要性。既有原两物种六条TRAIN训练／并行范围不变；不访问封存数据、新物种或新面板，不扩展冻结DEV家族。以下保留固定协议和历史执行证据，不再代表仍待完成。

## 问题和唯一干预

只问：在相同已训练起点之后，追加学习中的coarse-only观察信息是否改变完整CDS链恢复？已完成的只读13495642／13495706中，首个coarse-bearing视图无额外logit信号，但每物种最大coarse视图有非零信号。它们不证明性能改善，却使学习对照有实际机制依据。方案在这四个TRAIN固定点和历史DEV结果已知后确定，不称DEV前预注册，也不依据正在运行的fresh×joint最终性能选取。

信息分支复用13494895及其将来的严格final9216原TRAIN输出；不重训。对照分支只把追加fit中的epsilon0.1改为2/3：三个coarse region的观察概率均1/3，保留原标签与绝对观察常数。fine仍是完全相同的硬观察，unknown仍边缘化，完整gene-set竞争、空集和background校正不变。所有draw继续除以原max(1,fine+coarse)，不改成fine-only分母，不添加辅助loss、角色权重、长度或null权重。ε=2/3是数学上的无信息通道，不是搜索或拟合的噪声率。

## 固定学习和读出

两分支都从原fit13434603的step004608头和完整AdamW状态开始，原参数／moment步数均4608；复用原seed1 fresh计划的4608个有序draw、4219个unique窗口、三段各1536步和importance1。使用相同joint24576原生产编码、原BF16→FP32特征物化、nearest-center ownership、原有限视图、所有十五角色、cap128、FP64分区及checkpoint_bp256，冻结encoder且不使用缓存特征／native CDS先验。头、训练模式和随机种子沿用原fresh×joint分支；lr、weight_decay、betas、optimizer eps、clip均不改。

正式追加4608步，唯一科学读出为累计9216／epoch6，不选择6144或7680 checkpoint，不依据TRAIN或DEV丢弃成功分支。原初始化和前4608步已使用informative coarse，因此本对照仅能检验追加学习阶段的边际作用，不是“从零完全没有coarse”的消融，也不能证明通道在完整训练中的必要性。

## 实施前接入和完成条件

先实施独立配置、输出根和严格raw adapter。checkpoint必须明确记录无信息通道及同一sampling／encoder身份；错误通道、协议、累积步数、encoder revision或draw设置不得被当作匹配模型。保持原13494895源文件和输出不变。

native preflight固定使用上述每物种最大coarse视图，从同一004608源各恢复一次独立模型／AdamW；核对uniform loss=fine-only loss+C×log3、logit梯度相同，且两者均用原分母。各只做一次disposable update至4609，验证原optimizer完整恢复、参数有限且更新，更新后训练势函数与严格raw adapter实际逐元素一致。final9216 metadata仅作未保存routing fixture，不冒充正式权重。预检不写checkpoint、不改源权重或moment，不访问DEV／sealed。失败按具体问题停止解释和正式提交，不用另一视图或超参数替换失败。

预检通过后运行一条固定正式fit，再核对4608条实际draw／step／观察计数／原分母、三次checkpoint的有限头与AdamW状态、所有原数据读取边界。六条原TRAIN经严格raw FASTA注释，保留原whole-oriented-contig caller及完整链规则；所有原参考（包括不受支持者）仍留在22688分母。评价在原DNA与两个分支实际caller FASTA逐碱基一致后进行，原参考不进入encoder、head或decoder。

主要读出是informative−uniform的完整链TP、预测链数及count-derived precision／recall／F1，逐染色体、逐物种及pooled micro都保留；P／R／F1另列species macro及worst species。macro为两个物种指标的算术平均，worst为每个分支各取两物种该指标的最小值再作差，不与两个物种配对差的最小值混同。同样保留原CDS-base指标，不能以base收益替代完整链结果。全部成功结果完整保留在本项目记录，不因方向、显著性或物种不一致筛选。单次同起点对照描述这一次追加训练的总作用；不同梯度、裁剪及AdamW轨迹是干预后果，不称独立重复、一般因果机制、普遍必要性或泛化。原始uniform NLL包含entropy常数，不能直接按其大小和informative比较学习好坏。

## 运行、监控和下一决策

实际沿用原fresh×joint资源：预检单24GiB 3090／2CPU／16GiB／20分钟，fit单同类GPU／2CPU／16GiB／18小时；verify为2CPU／4GiB／5分钟，raw六任务最多同时3个、各同类GPU／2CPU／16GiB／2小时，评价为2CPU／4GiB／20分钟，配对为2CPU／4GiB／10分钟。提交前已核对实际private-teodoro-gpu队列／GRES／维护、空闲容量及新路径；原joint24完整fit实耗47490.182秒（13.19小时），支持18小时上限。fit控制器已实际分配gpu034上的1张3090、2CPU、16GiB，其余阶段仍等待依赖，并非已分配。各作业Requeue=0。fit沿用96MiB输出上限，raw各沿用128MiB；不建立新框架或额外缓存。

仅复用原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3和唯一心跳coding-rna-r5-light-monitor，保存了这条新流程及原13494895流程的精确job和完成产物；不新建对话、cron或重复自动化。正常RUNNING／PENDING保持静默，主对话消费终态及科学结果。两个事件分别记录交接：任一整流程首次终态或需行动故障只交接一次并暂停同一心跳；主对话消费后为尚未交接的事件恢复同一心跳。基线完成不等于uniform完成，uniform评价等待基线的有效依赖也不算故障。完整对照结果决定是否有值得继续确认的粗信息性能作用；无论结果如何，不自动放行sealed、改DEV家族、增加seed／步数或宣称论文达标。当前论文goal仍未完成。

## 独立接入实施

新增annotate_full_grammar_uniform_coarse.py及fit_full_grammar_uniform_coarse.py，输出根分别为UNIFORM-COARSE-PREFLIGHT、FIT和VERIFICATION；原fresh×joint五个流程源文件及配置不改。新模型继承原sampling／encoder／原参数容量身份，仅改变追加观察epsilon及独立模型／协议标签，checkpoint严格拒绝错误epsilon、sampling_seed、unique_windows和encoder_blocks_bp。新raw adapter只导入原raw推理实现，复用joint24576完整contig调用，不导入训练／标签读取函数或改变输出规则。

预检实现使用同一原source4608模型的native logits计算uniform和fine-only两项分区，保留原fine+coarse分母，并核对绝对entropy常数及logit导数；随后各做一次独立AdamW更新并检查实际raw适配器同值。两物种仍固定最大coarse且最早draw的选择规则，未改为更容易通过的视图。formal fit和verify复用已验证的原加载、数据接入、抽样与原生产势函数；独立记录全部4608draw、三个checkpoint及同参数AdamW状态。新增四个Python脚本AST及六个SBATCH语法已通过；提交前对应本地／远程源文件逐字节相同。

## Native预检实际结果

13495941的root／batch／extern均COMPLETED、ExitCode0:0；15:56:49至15:58:41UTC，脚本实耗95.805秒，stderr为空。Arabidopsis固定additional step4106，fine17449、coarse7127；Rice固定step4048，fine7102、coarse17474；两者原分母均24576。

| 原TRAIN固定视图 | uniform减fine-only减C×log3的残差 | 两项logit导数最大绝对差 | 更新后训练势／严格raw势逐元素相同 |
| --- | ---: | ---: | --- |
| Arabidopsis最大coarse | −2.7284841053×10⁻¹² | 0 | 是 |
| Rice最大coarse | −1.0913936421×10⁻¹¹ | 3.3306690739×10⁻¹⁶ | 是 |

主对话另按输出值复算entropy恒等式、raw NLL除以原分母及梯度差。实际uniform配置与正在运行的信息基线started配置只有model标签、epsilon及coarse_information三处差异。两次disposable更新均从原4608完整AdamW恢复至4609，未保存checkpoint；源权重／moment未变，encoder冻结，无缓存、DEV或sealed读取。这里验证数学与接入，不是准确率结果；临时9216 metadata fixture不代表存在正式9216权重。

## 实际提交和依赖

下表为16:06:06–16:06:07UTC实际提交；16:07:56UTC控制器回读确认fit RUNNING，其余均PENDING／有效Dependency。

| 阶段 | 实际Slurm job | 已核对的afterok依赖 |
| --- | --- | --- |
| Uniform正式fit | 13496011 | 已完成13495941作为显式参数和启动时严格预检门槛；无Slurm依赖 |
| 最终权重／draw／AdamW核验 | 13496012 | 13496011 |
| 原六条TRAIN raw-FASTA推理 | 13496013_0–5，max3 | 13496012 |
| 原共同尺评价 | 13496014 | 整个13496013数组与信息基线评价13494898均成功 |
| Informative−uniform配对 | 13496015 | 13496014 |

最终权重预期为outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-FIT/step_009216.pt。新infer_full_grammar_uniform_coarse_train.py复用原raw输入和共同评价实现；评价先逐碱基核对两分支所有六条caller FASTA与归一化原DNA，再读取uniform参考并记录paired_input_match.json。新compare_full_grammar_uniform_coarse_train.py核对两边实际job、参考政策／逐染色体几何／全部22688参考及count-derived P／R／F1，保留chain与CDS-base的逐染色体、逐物种、pooled、macro及各分支worst差。RAW实际numeric root尚未分配，不猜arrayID+task。

同一心跳16:07UTC更新成功，实际保存的heartbeat类型、原目标对话、ACTIVE、15分钟、failed_runs_only、created_at及完整7942字符双事件prompt已回读核实。原对话人工turn01a12bbb-a9f1-7392-ba6a-2cc525912b21于16:12:02UTC检查B+C；其后真实心跳输入16:12:45.617UTC包含原automation_id、调度时间及新双事件范围，scheduled turn01a12bbc-f276-79c0-831b-963dc1f40506于16:12:58UTC检查B+C。主对话核对两次实际turn_context均为gpt-5.6-luna／low；两次都发现两条fit RUNNING、其余有效Dependency，无科学内容读取、交接或暂停。配置成功、人工检查和新范围定时追加已分别验证，终态回传仍未发生。当前只等待外部计算，由这一原固定心跳静默监控，不形成主对话重复等待循环。
