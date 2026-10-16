# 完整 CDS 事件后验与原 TRAIN 接入预检

2026-10-13（北京时间）。当前模型已有合法 complete gene-set 的all-path训练分区和MAP调用器，但尚未提供整条预测CDS链的结构边际概率。本轮新增只读后验读出及原TRAIN数值／成本预检，不训练、改温度、筛预测或以历史DEV选择方法。其可能用途是排序与不确定性分析；数学成立、原生可执行、置信度校准和准确率提升是四个不同判断。

## 固定概率与结构范围

在给定一条oriented DNA及原十五角色势的条件下，`Z = ΣS exp(E(S))`仍覆盖全部合法complete gene sets和空集，不限于MAP输出或另枚举的候选池。基因`g=[a,b)`的内部CDS／intron角色决定`E(g)`。`α_B(a)`汇总a之前以BACKGROUND结束的完整外部sets；`β_B(b)`汇总从BACKGROUND开始的b之后完整外部sets。因同方向基因span非重叠，

`log P(g∈S | DNA, scores) = α_B(a) + E(g) + β_B(b) − log Z`。

完整给定set的point probability是`exp(Σg E(g)−logZ)`，空集为`exp(−logZ)`。单gene事件会边际化其他不相冲突的基因，不能用其内部score的sigmoid，也不能一般地把逐gene边际相乘当作完整set概率。内部查询点不重置partition的开放状态；BACKGROUND只约束该gene事件之外的两段。两个方向仍独立，不加入cross-strand排斥。

复用原`_transitions`、实际未完成codon carry、GT／GC-AG与min-intron5、起止密码子连续性、角色precedence及每gene cap128。源码新增独立`full_grammar_gene_posterior.py`，原训练、raw skyline和边界／梯度模块不改。两遍FP64消息只保存被查询边界的BACKGROUND标量；仍扫描完整供给DNA和全部结构状态，不winner-prune。保留log probabilities用于下溢时的排序；仅≤1e−8的P>1取消误差截回0，实质异常停止，不隐藏。

条件分布／动态规划并非新概念，[semi-Markov CRF先例](https://proceedings.neurips.cc/paper_files/paper/2004/hash/eb06b9db06012a7a4179b8f3cb5384d3-Abstract.html)及项目既有Conrad比较继续有效。本轮只检验该具体CDS事件在现有DNA与carry语法中的实现，不仅凭forward/backward或“posterior confidence”主张原创，也不把模型结构支持概率解释成生物学正确率。

## 实施验证与固定原生案例

四组新增CPU人工检查：独立DNA/CDS连接后枚举全部兼容sets，比较Z、每gene事件及exact-set／null概率；相邻两基因与共享STOP竞争（零势时单gene边际分别1/2与1/3）；actual carry／短exon precedence／GC和intron N／128与127上限拒绝；row-shift、null、log下溢和输入／梯度不改写。只跑本轮四组，不重复未变化的底层全套测试。CPU作业2CPU／4GiB／无GPU／5min，generanno／no-requeue；资格由主核验后才提交native。

固定E rank4末9216 fit13499817与M dense末9216 fit13510324，native age均4608，原b048 source及所有权重只读。不是本轮仍在训练的两个LR格，不从中间checkpoint挑选。仅复用E已完成strict raw TRAIN FASTA及其原DNA相等收据，不读参考GFF、标签、训练catalog或feature cache。

两个genomic core：AT NC_003070.9 start14389248、rice NC_089046.1 start14462976，均width24576，原stride12288邻居组成49152bp供给视图；与既有native首draw见证相同，不按后验或准确率挑例。两模型×两物种×双strand共八案例，实际原生joint24576／BF16→FP32／原ownership／原skyline MAP。每例计算所有原MAP基因的事件后验与完整set／null后验，确认gain不改；与原compiled all-path Z比较，再将第一条预测链的角色作为事件约束、外部unknown，独立核对numerator。无预测的案例改检验null事件，必须报告是否真正检验了非空gene，不能以null代替非空证据。

单RTX3090／2CPU／16GiB／30min，先按实际调度选择同硬件的合格private或shared入口，输出小于4MiB的JSON／源收据，不保存势数组、模型或新预测GFF。资源时限是预检上限，不预言full-chromosome成本。结果保留八案例原始log概率／genes／runtime／峰值；不读truth计算F1、AUROC或校准，不定confidence阈值。

## 下一步决策与边界

仅当独立枚举及真实native Z／事件概率／不改MAP验证成立，才根据实测代价准备同original TRAIN整染色体的置信度读出／独立排序评价。49152bp视图后验不是整染色体后验，不能将两者混名；可执行不等于排序有效、校准有效、优于Carbon-A confidence TCN或论文已达标。若原生极慢、全为空或概率异常，保留结果并据实际问题选择优化或停止该读出，不改训练／temperatures／exon cap寻求漂亮成绩。

输出roots为`M28-FULL-GRAMMAR-GENE-POSTERIOR-CHECKS/<actualjob>`与`M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-PREFLIGHT`。实际取得handles后记录；正常长等待沿用原watcher01a0f71f-8845-7432-bd13-bf54c1808fb3及原coding-rna-r5-light-monitor，不创建新chat／automation。旧E/M只作source，不重新监控已关闭链；原LR两链与13545908四格汇总不改。明确接线／调度错误至多一次同科学范围修复，数学／数值负结果不偷偷变定义。范围仍original TRAIN，不增加DEV方法、放行sealed／held-out／Setaria或对外发布。Nature Communications及以上的完整竞争性模型目标仍未达到，不缩减为预检PASS。

## 实际资格结果与原生提交

新增三个Python入口AST及两个SBATCH语法通过，六份新源／协议精确同步并回读后才提交。02:14:43UTC实际调度显示合格3090节点全部GRES已分配，private16张与shared同型节点没有更早空闲卡；CPU有可用资源，维护预留仅11月6日CPU。GPU025实际6张而非旧配置8张，按live GRES取值；没有根据nvidia-smi推断空闲。首次只读`scontrol show node`参数写法被拒绝，改为每node独立查询后完成路由检查，未创建作业或改变调度。原两LRfit仍RUNNING、下游有效Dependency，不作重复结果监控或取消其他项目。

CPU13546934于02:16:32–02:17:24UTC（北京时间10:16–10:17）、52秒，root/batch/extern全部COMPLETED0:0；四组新增独立枚举／竞争／carry与cap／不变性检查全部通过、0error／failure，测试实耗5.501515秒。没有研究DNA、参考、checkpoint、native模型或DEV／sealed读取，无研究更新；主直接消费，watcher handoff=NOT_SENT_BY_WATCHER。数学通过不是原生部署或校准结论。

实际原生13547015于02:22:38UTC（北京时间10:22）提交：`sbatch --parsable sbatch/M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN-PREFLIGHT.sbatch 13546934`。提交前own expanded仍30<80、独立输出根不存在；单RTX3090／2CPU／16GiB／30min／private／Requeue0及完整实参已由scontrol核对。无需将已成功CPU或旧E/M IDs加作controller依赖，程序读取真实资格／source收据。首次队列为PENDING/None，尚未启动或产生原生结果；不得把成功提交称作八案例完成。

同一原heartbeat现包含两个独立等待事件：原LR全链批次QV-LR-FOUR-CELL-13545908与GENE-POSTERIOR-NATIVE-13547015，交接均NO。POST只需自身root/batch/extern、started／summary／STATUS／JOBID及八唯一案例结构字段成立，不依赖LR结束；任一完整终态或需处理失败先向原main一次真实交接并暂停同心跳，另一未交接对象保持NO，由主处理后恢复。原monitor/main／Luna-low／15分钟／通知偏好不变；不新建chat／automation，不让watcher解读概率或更新代码／数据／goal。新scope配置和实际定时追加需分别取证，完整自动回传尚未观察。

## 原生完成及实际可支持结论

13547015于02:32:48–02:34:51UTC（北京时间10:32–10:34）、123秒，root/batch/extern全COMPLETED0:0；STATUS／JOBID／688字节started及23505字节summary齐全，脚本实耗93.189711秒。主核验八唯一model×species×strand案例、固定9216／4608 sources及DATA限制；每例前向与后向Z、原compiled all-path Z和直接事件numerator相同，原MAP链／gain不改。独立复算exact-set／null的log mass及exp、逐gene的exp关系，以及exact-set事件包含于每个gene事件的一致性，通过后由主消费。无权重／预测改写、更新、参考／标签／cache／DEV／sealed访问、阈值或校准测量。

| 模型 | 供给视图与方向 | MAP完整链数 | 两遍后验实耗 秒 | 非空gene边际范围 |
| --- | --- | ---: | ---: | --- |
| E_RANK4 | NC_003070.9 + | 0 | 6.065058 | null only |
| E_RANK4 | NC_003070.9 - | 0 | 5.055360 | null only |
| E_RANK4 | NC_089046.1 + | 1 | 5.135851 | 0.440584 to 0.440584 |
| E_RANK4 | NC_089046.1 - | 2 | 5.129123 | 0.431021 to 0.911912 |
| M_DENSE | NC_003070.9 + | 0 | 5.200542 | null only |
| M_DENSE | NC_003070.9 - | 0 | 5.086896 | null only |
| M_DENSE | NC_089046.1 + | 1 | 5.154990 | 0.517892 to 0.517892 |
| M_DENSE | NC_089046.1 - | 2 | 5.133100 | 0.463160 to 0.952541 |

四个AT案例只有null事件，不将其当作真实非空gene验证。四个rice案例包含每模型三条、共六条被查询的MAP基因；第一条事件在每非空案例都由直接角色约束numerator复核。六个边际为0.431021–0.952541，且对应原MAP基因内部gain很大。这说明相对于null的高gain不等于整链事件概率趋近1，其他合法结构的竞争仍影响概率；不说明哪些gene真实正确，也没有accuracy／AUROC／校准结果。八例使用固定49152bp供给视图、19078结构状态，仅存2／4／6个查询边界标量；不称整染色体置信度。稳态后验约5.05–5.20秒，首次6.065秒含编译影响；GPU峰值allocated最多2561631744字节。完整染色体代价与排序效用仍需实际执行，不把这些局部成本线性换算成完成保证或相对其他工具的优势。

原watcher在02:45:54UTC的真实定时turn01a13326-f82b-7b80-9183-665c756a4d51（Luna／low）收到了完整9487字符新scope，只作合并squeue查询；POST未出现在返回行，却被当作普通等待，未按约定查sacct。这是监控终态识别缺口，不是科研作业失败。本轮main已直接完成sacct与产物／实际数值核验，POST CLOSED／main-consumed／NOT_SENT_BY_WATCHER，不伪记自动终态交接。现仅等待原LR批次，后续规则明确：预期事件末端缺行即查该ID的sacct与markers，不得推断为正常pending；不新建监控或重提已完成POST。
