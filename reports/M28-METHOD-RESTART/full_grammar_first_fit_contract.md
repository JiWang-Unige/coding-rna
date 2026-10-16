# 完整基因语法观测模型首轮拟合合同

这是新观测目标的首轮真实拟合合同。2026-10-08，用户直接回答“我同意你的方案”，批准此前明确提出的固定epsilon0.1、原抽样与importance、4608步、不用DEV及下列30小时资源方案；这次批准来自人类答复，不由准备产物、配置相同或goal状态推定。实际job13434603于北京时间2026-10-08 10:28:35开始、2026-10-09 01:13:03结束，root／batch／extern全部COMPLETED0:0，4608步及三份权重齐全；private／单3090／2CPU16GiB／30小时配置不变。只读数值及权重核验R1 job13442813已COMPLETED0:0，正式checkpoint仍固定最后4608步。详见[首训结果](full_grammar_first_fit_result.md)。原目标仍是严格raw-FASTA完整gene注释及对ANNEVO、Tiberius、Helixer的性能／适用范围比较；本次TRAIN discovery拟合完成不替代最终整染色体、独立泛化或论文证据。既有负例与准备历史见[准备记录](full_grammar_partition_preparation.md)。

## 模型与训练范围

只使用已授权的原Arabidopsis／rice六条TRAIN染色体及原1,536 draws／1,491唯一窗，连续三遍原顺序，共4,608次更新；每物种每遍768次曝光。原窗口、primary、fine／coarse／unknown、原importance及所有旧结果不变，不按不可行性过滤24窗。暂不扩大为已有25,512个natural TRAIN窗，也不使用DEV、独立物种、Setaria、A33或重开A44；discovery之后仍须推进与论文目标匹配的更完整训练／外部评价。

冻结released GENERanno-eukaryote-1.2b-cds-annotator-preview，revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、k6、6144bp native block、BF16 hidden转FP32正常tensor；每draw原生重新编码，不继承R4/C0/B1/A44监督权重、不持久缓存新features或跨更新复用C0 logits。fresh seed0 M28Core C0，2048→128、366,191个FP32参数，训练shared project／双向GRU／DNA conv／fuse及十五角色输出层，GLM无梯度。

拟合view沿已验证的raw坐标规则：24,576bp core，两侧各一12,288bp stride邻窗，按实际染色体长度裁剪／去重；nearest-containing-center逐碱基归属，平局取较低genomic start，负链只反转片段顺序。无注释派生halo、Gold边界初始化或候选池。内部core边缘传递精确消息，只有整个固定view外缘闭合background。

## 目标与优化设置

每个完整且非重叠gene集合G的分数是原十五角色logit相对background的单位增益和，null=0，FP64累积、cap128、min-intron5、严格ATG／单exon完整terminal STOP／actual coding carry／GT或GC–AG。fine精确角色不变；每个fine未知而coarse已知的位置，正确region类别观测概率为1−epsilon，另外两个类别各epsilon/2；unknown行的观测概率为1。

首轮经上述用户批准固定epsilon=0.1，不训练它、不扫范围、不由24/1491窗口比例估计逐位置误差率，不根据DEV重选。该值是明确的首轮机制试验参数，不是经验证的真实注释错误率。此轮只回答神经评分在同一保留原标签的目标下能否学习；避免同时学习观测误差率而混淆两个机制。推理不使用epsilon或标签，只用训练后的神经分数和原严格语法。

原始NLL为logZ(theta)−logZ(theta, observations)。实现用logits+logq计算相对观测配分后补回sum logq(background)，不丢弃其绝对概率因子。每draw目标为原importance × NLL / max(1, fine行数+coarse行数)。分母与原local_loss的已知观察计数规则相同，但联合结构NLL不等于原局部CE；这个替换本身是新机制，不冒充曝光以外全部因素匹配的因果对照。不加local CE、链BCE、长度权重、温度、平衡类别或其他辅助loss。

沿原优化设置：AdamW lr3e−4、weight_decay0.01、betas(0.9,0.999)、eps1e−8；全C0梯度clip_norm1，batch size1，无scheduler或新增seed。每draw重建原R4目标、原raw view和新C0计算图后更新一次；只在每个原unique core首次使用时核对其实际raw DNA与原缓存，缓存feature／旧标签不作为模型输入。任何原计数／view漂移、不可行、非有限NLL或梯度均停止，不跳过draw或变更规则。

实际准备13429756发现原draw序号612及640重复同一个全unknown窗口arabidopsis_thaliana|NC_003070.9|14770176|-：fine=coarse=0，unknown=24576，分母按原规则为1。这两次曝光每遍保留，观察NLL及其梯度为零；仍按上述原顺序执行AdamW步骤，其已有动量及weight decay可能改变参数，不将其称为来自观察的学习信号。三遍4608次optimizer步骤不等于4608个非零观察梯度步骤。

## 结果与停止条件

每遍保存step_001536.pt、step_003072.pt、step_004608.pt，含模型及optimizer状态；只使用最后4,608步权重作为此轮正式checkpoint，不以TRAIN loss或DEV选择中间权重。逐步记录原draw序号／ID／物种、importance、fine/coarse/unknown数、分母、raw NLL／归一NLL／加权目标、clip前梯度范数及耗时。完整训练验收要求原三遍曝光顺序、4,608条连续日志、三checkpoint及summary齐全、全部原core DNA已核实、GLM仍冻结、C0参数实际改变。

训练loss下降只证明当前目标的TRAIN拟合变化，不证明完整链准确率提升、泛化、SOTA或论文PASS。出现数值／数据矛盾、TIMEOUT／OOM、异常更新或产物缺失则保留失败与已完成checkpoint停止；不自动重提、加时、续训、换epsilon／normalizer／importance／seed或修旧A44。明确纯实现错误最多允许同科学设置的独立最小R1一次，不能以修复名义改变本方案的科学选择。

## 完整染色体输出边界

后续raw caller必须只接收FASTA及最终checkpoint，按既有manifest.starts的完整raw网格原生编码，使用固定base-owned角色势，在整个定向染色体上执行一次完整非重叠gene-set最大分数解码。窗口只是神经特征视图，不能逐窗关gene后owner过滤，也不读取参考／catalog／缓存预测或按参考选择候选。输出gene、transcript、CDS及实际phase，完整ORF／motif／cap和坐标转换逐链核实；零分null保留，无新操作阈值。

已准备annotate_full_grammar_coarse_channel.py，复用contig stitch，当前后端改为独立的full_grammar_viterbi完整状态自动机／分块重算回溯。原13429790四项人工接口检查实际通过，但其旧endpoint编译依赖存在未遵守的A43停止条件，不能由四项结果解除；原分支及失败保持关闭。新表示已由13429852八项实际检查通过（四项独立最优分数／cap／回溯数学及四项新后端raw接线），包括两链30000bp基因跨24576bp窗及负链GFF phase；另有13429908两项高密度49kb测试取得与独立Python相同的整数最优分数。新FP64逐碱基累积及background-first／state-table-arrival同分规则已明确，不宣称旧prefix累积和同分链身份逐位一致；训练配分目标及首训科学配置未因此改变。30427671bp人工染色体资源测量13429995已实际TIMEOUT，20分钟内未完成解码或逐链核验，batch峰值约5.33GiB；这项有界测量按停止条件不加时或重提，染色体资源可行性未证实。尚未读取真实训练checkpoint或运行released GLM推理；不得把人工接口通过或固定49kb上下文支持宣称为整染色体完成，资源超时也不等于训练目标的数学或科学FAIL。六个完整primary超出诊断view的warning仍存在；此轮是有限上下文core事件的composite拟合，不是整染色体所有基因集合的联合MLE或六条完整primary恢复的证明。这些问题须在原目标中继续解决，不通过改变成功定义回避。

## 执行入口与资源

fit_full_grammar_coarse_channel.py prepare只核对当前TRAIN元数据、已完成1491支持／24候选／四窗native证据，生成原1536行拟合计划及拟议loss multiplier；不读raw DNA、feature tensor或模型，无optimizer更新。M28-FULL-GRAMMAR-COARSE-CHANNEL-FIT-PREPARATION.sbatch为generanno/private/0GPU2CPU4GiB2分钟/no-requeue，独立准备根及M28COARSEFITPREP_%j日志，明确实现错误允许同条件独立R1一次，数据／数学矛盾、TIMEOUT／OOM或再次失败停止。

M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT.sbatch为拟议真实训练入口：generanno/private/1×3090至少24GiB、2CPU16GiB30小时/no-requeue，独立FIRST-FIT根、M28COARSEFIRSTFIT_%j日志，新增fit产物上限64MiB。30小时不是旧实验预算扩张或已分配成本；依据四窗native稳态约16秒/draw估计4,608步约21小时并留运行余量，不能保证全队列时长。提交前重新核对活Slurm、预约、存储、所有输出路径及明确首训授权；准备记录或CONFIG相同不是用户批准。没有自动afterok训练提交，准备成功后仍不自动启动首训。

只等待长训练时复用原固定监控对话与heartbeat；正常RUNNING不汇报，终态一次交接后暂停。完整模型评价及任何对外发表／分享另按其既有授权范围处理。

2026-10-08新增人工完整长度skyline资源实测13434312已COMPLETED0:0：30427671bp、2476原raw网格、cap128／checkpoint8192不变，decode含cold JIT537.100558秒，逐选中链检查30.284062秒；57968条完整chain、最大122外显子，全部实际ORF／carry／terminal／splice／cap／非重叠／整数分数重算通过。该人工输入证明了新表示的资源可行性，不是完整长度独立穷举最优oracle、真实GLM势的运行保证或准确率结果。原dense13429995 TIMEOUT仍关闭，正式raw源此时仍用dense；13434348已通过的前瞻skyline接线与本结果支持后续替换，不复活停止分支。训练logZ实现和本合同科学设置不变。

## 已完成的准备结果

13429756实际2026-10-07 17:08:26–17:08:42 UTC，COMPLETED 0:0，gpu034/private、0GPU2CPU4GiB2分钟、Restarts0，batch MaxRSS848884K。1536行原计划与原draw逐行ID／顺序／stratum／importance相同，1491个unique的原view、fine／coarse／unknown及硬支持状态逐值核对通过；两物种各768 draws，unique分别745／746，24个原硬负例全部保留。五个非空产物合计738946 bytes，started／summary配置相同、stdout与summary逐值相同、stderr为空。原生四窗稳态外推20.8066小时，仍不是整队列实测。FIRST-FIT根不存在，实际模型、feature和raw FASTA均未读取，optimizer步骤0；这里只交付了准备，不是训练或科学成功。
