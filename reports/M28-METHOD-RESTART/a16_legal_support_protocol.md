# A16 合法事件概率质量诊断

2026-10-04，运行新探针前固定。A15按原协议继续，不改其训练或门槛。本轮只诊断冻结A3生成器的训练／解码事件域差异，不拟合新模型、生成新自然预测或进入DEV。

## 当前证据与问题

既有A5原8个metadata-selected TRAIN压力窗、27个参考／窗口配对中，自然beam命中15，遗漏12。八个遗漏在全部query至少一处局部必要名次超出前四，另四个有局部可通过query；A8完成池重排和A11组件均值都未恢复新正确链。仅重排后段候选不能解决全部缺口，也不能据此断言搜索宽度或链长度是唯一原因。

实际 `query_pointer.py`／`query_training.py` 的teacher边界CE在几何合法位置归一化，后续entry下界是previous_end＋1；`query_decode.py`／原A5 probe的自然事件域则受ATG、GT或GC–AG、最短5bp内含子、真实codon carry／ORF和STOP约束，entry下界previous_end＋5。自然结束得分还乘prospective donor-continue或STOP-terminate概率。此差异是真实接口语义，不先当作实现bug或模型失败的因果解释。

本轮量化正确前缀下，原几何pointer概率有多少落在解码合法事件中，及合法事件的联合boundary／EOS概率质量。若非法域承担明显路径分数代价，下一结构假设应研究训练与合法事件域对齐，而不是继续给完成链加事后长度均值；若质量接近1或其他排序缺口仍占主导，不能以此诊断为理由添加新归一化训练。仅对同一前缀常数重归一化不会改变合法事件局部排名，不能声称它直接修复八个局部失败。

## 固定输入及公式

复用已闭合 `outputs/M28-A5-TRUE-PREFIX-DIAGNOSTIC`（13359200，COMPLETED／0:0）的27参考、432个query profile和原 `M28-A5-TRAINED-TRAIN-DIAGNOSTIC` 自然预测；冻结A3 step004608、Q16、原TRAIN四块特征及有向DNA。保留全部8窗，包括零参考窗；全部27参考×16query，不挑最优query、改窗口或过滤复杂链。所有参考均已由A5验证原语法可达；当前几何、合法entry／同种end名次、carry和raw分数必须与保存probe一致，raw总分最大差≤1e-3，否则停止解释。

entry的原几何log-prob记为l(a)，合法起点或acceptor集合E。记录 `log P(E)=logsumexp(l[E])` 与 `l(a*)−log P(E)`。end合法事件以(b,donor或STOP)区分，联合log-weight为原几何end log-prob＋对应prospective EOS log-sigmoid；记录合法联合事件log-mass、目标条件log-prob及不含EOS的合法唯一end坐标mass。后者不重复计算同坐标事件，也不把joint mass当作仅DNA mask mass。

这里的合法域是当前一步的解码支持，不保证每个起点／donor最终可完成整条基因；不是全部可完成CDS路径的全局分区函数。

将上述条件值相加只作路径分数的诊断分解：条件路径分＝原路径分−每步合法log-mass之和，不执行新的自然beam，不与其他路径形成新的top4准确率。另记录teacher与decoder entry下界改变造成的真实目标log-prob差；end下界在二者相同，该下界差为0。全部值有限、每个原概率质量≤1（FP32容差1e-5），条件分数分解误差≤1e-5。

每reference先在其16个query和对应事件上计算均值，再按reference等权汇总，分别报告物种、exon类及原自然beam命中／遗漏；单exon没有next_start，保持缺失而非当0。27是reference/window pairs，432是条件query探针，不是独立基因样本或训练总体估计。强制参考前缀是teacher诊断，不是自然模型输出。

## 实施及下一判断

新增独立 `legal_support.py`、三检查和 `analyze_legal_support_train.py`，不改原A3、A5、A15代码／权重／结果。三检查验证精确质量分解与合法域内rank不变、EOS联合mass不同于单纯boundary mask、条件NLL梯度的非法坐标项确实抵消。真实probe重放原参考身份、合法名次／carry／raw分数；不重复旧A5的未改语法测试。

独立输出 `outputs/M28-A16-LEGAL-SUPPORT-DIAGNOSTIC`，query_profiles.jsonl 432行、references.jsonl 27行、summary／STATUS／JOBID。generanno，private／1RTX3090／2CPU／16GiB／10分钟，依据A5同432条件探针实际88秒留余量，非吞吐基准。先三检查后同一作业内跑数据，optimizer0、无checkpoint、新候选搜索或GEN编码器前向；仅冻结任务头前向。无DEV／sealed／Setaria、其他作业修改、Git提交或公开。

随后一次只读CPU结果核验使用独立 `M28-A16-LEGAL-SUPPORT-RESULT-CHECK`，private／0GPU／2CPU／8GiB／5分钟。核对432个保存profile与原query／参考身份、局部名次／真实carry／raw组成分；按保存事件直接重算条件概率分解、27个reference均值和全部分组汇总，另汇总不含EOS的合法end坐标mass。无新头前向、训练或预测，本短诊断及核验由主对话直接处理，不加入A15三根监控对象。

诊断后直接核对实际事件分解与汇总，再决定是否值得为新合法事件训练写独立固定协议；本探针无准确率gate，不自动启动新训练、调A3／A15或把条件分数升级为胜出结果。新模型仍须完整raw-FASTA、全部参考、逐物种F1／背景、原生ANNEVO／Tiberius／Helixer以及独立测试证据。

## 论文贡献边界

生物学语法参与可微训练不是首创：[Tiberius原作者论文](https://academic.oup.com/bioinformatics/article/40/12/btae685/7903281)已将HMM用于训练和推理，包含阅读框、起止密码子及剪接模式约束。本轮是当前pointer实现的概率域诊断，不是新HMM、CRF或全球最优解码；未来论文贡献须来自可验证的新模型效应、完整输出与适用范围，不能仅靠“加入语法”声明创新。
