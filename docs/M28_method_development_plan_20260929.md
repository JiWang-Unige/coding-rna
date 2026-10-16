# M28：恢复新方法主线

2026-09-29。当前用户目标：提出可执行纯FASTA ab initio gene annotation的新模型，对标ANNEVO、Tiberius、Helixer，以更好的完整基因性能或经验证的更广适用范围，争取Nature Communications及以上期刊。不是把现有失败结果凑成论文；也不保证必然得到创新或录用。

## 本轮已完成

- 两个远程CPU作业并行完成开发数据盘点；另一个原生dataset计数核对完成。见[结果](../archive/m28_pre_full_grammar/reports/census_result.md)。
- M25R排序后prefix采样导致训练窗口全来自拟南芥，且计数与历史摘要一致。新模型必须比较修正采样、同数据/预算的对照。
- 24576bp半步长窗口可几何容纳6439/6450个开发参考链；尚未含侧翼、motif/候选限制，也未证明模型预测能力。
- 第16轮Pro商议完成（8分58秒）：实际读取73760c7三份指定公开文件，未访问Baobab、未独立重跑。建议首攻B1，A保留为真正不同的备选。不是外部同行评议。

## 首攻：B1，重新学习候选与非加性完整链竞争

输入只有DNA和固定权重。局部GLM特征先冻结，以小型跨块上下文模块连接局部窗口，配合碱基分支保留单碱基定位。

- C0：相同多尺度特征、15类分割与标准结构解码，不能依赖旧run图。
- B1：新的端点/连接候选生成头；链编码器读取exon特征、间距与顺序，学习完整链与错连、截断、邻gene竞争以及“无基因”输出的区分。局部监督与链损失共同训练；离散候选的梯度边界明确，不假称完全可微。
- 同底座、输入样本与顺序、单seed、可比训练预算。头从头训练，不继承受旧DEV选择影响的LoRA。增加数据或修复采样的收益不能归因于B1。
- A仅为备选：gene-query存在性＋坐标指针预测边界序列和结束标记。只有B1被候选可达性限制且首个pilot有证据时再另测；不把两个复杂系统同时堆在一起。

候选主张只是待证伪假设：整链非加性兼容性和空输出竞争，能否把新的高召回候选转成可迁移的精度—召回改进。不能将GLM、HMM、结构化损失或“覆盖≠精度”本身称首创。

## 下一执行阶段：可训练manifest与最小可运行C0/B1

两个工作包可并行，都是新模型的实现准备，不再重复M26/M27审计：

1. 数据包：24576bp候选上下文，计算含侧翼/尾窗的真实覆盖、窗口gene数/exon数及motif候选规模；保存完整链身份和曝光次数、自然背景、partial/重叠/isoform掩码。负例只称参考不一致，排除匹配其他合格isoform的候选。标签不足时补训练数据政策，不扩大backbone掩盖。
2. 模型包：复用冻结GLM加载路径，构建共享特征和C0/B1最小实现；合成坐标/链/空输出测试及小样本梯度/自由推理检查。先测吞吐和缓存，再冻结更新数、候选限额与损失。任何训练用参考注入与完全不注入参考的推理候选必须分开记录。

R2更新（2026-09-30）：共享manifest、局部标签、跨块特征及整链/null核心已通过三个CPU作业与16项针对性测试，见[R2结果](../archive/m28_pre_full_grammar/reports/r2_result.md)。尚未接通真实GLM特征、C0原生解码和B1自由候选，也没有新训练准确率；下一步直接完成这些接口及吞吐测量，不再停留在旧模型审计。训练CDS完整性与历史parent过滤不同，正式主比较需先固定CDS-assessable评价并对全部方法同口径重评分，旧6,450范围保留副表，不按训练eligibility删参考。

R3更新（2026-09-30）：[真实特征与B1联合接线](../archive/m28_pre_full_grammar/reports/r3_integration_result.md)和[两物种共同评价](../archive/m28_pre_full_grammar/reports/common_ruler_result.md)已完成，尚无optimizer更新。主评价固定为7,728条CDS-assessable链；R2 longest-all与评价primary在拟南芥2个DEV gene存在差异，正式拟合前明确训练primary政策。接着完成GPU头吞吐、全染色体输出归属/去重/冲突和一次冻结的公平首训，不增加旧模型审计。

Pro建议48GPUh（特征8、两臂各16、评价8）、RAM≤64GiB、缓存≤200GiB，这是**建议上限，未直接转成提交或承诺**。先冻结更小的工程吞吐测量；实际可完成的公平pilot数据量/更新数与资源合同必须在读取该pilot性能结果之前写定。优先使用私有资源，不以并行为由越过队列、成本或旧数据边界。

## 评价与继续条件

首轮产物是primary protein-coding CDS链和GFF3，不冒称全部UTR/isoform注释。DEV：现有两物种允许染色体；所有长链/候选不可达链保留分母。报告每物种及宏平均exact-chain P/R/F1、参考未注释背景上的FP链/Mb、CDS及gene-span背景覆盖、推理候选整链覆盖、运行和存储成本。背景非生物学真负例，必要时另用公开独立转录/蛋白证据裁决。

模型、checkpoint选择、操作点冻结后，另登记未参与本项目开发选择的非Setaria近/远缘物种作一次性迁移评价；没有合格对象时先只作DEV筛选。不把同源重叠不清的划分或未知backbone暴露称干净zero-shot。ANNEVO等使用同一参考/输入/评价器重评分，正式新面板固定当前适用权重；不能拿原论文不同gene定义的分数拼成总表。

B1若在同预算C0之上，两DEV方向一致地改善完整链而不过度增加FP，再进入独立迁移和扩数据/LoRA。若推理候选漏失而只有注入参考的训练表现好，切换到A假设而不是扫候选阈值。如果两臂未拟合，结论是预算/冻结特征下未检验充分，不声称理论失效。阈值是研究决策，不是期刊录用标准。

## 原始资料与新颖性边界

本轮Exa定向检索共24个返回结果额度（8+6+4+6，非24篇全文精读），覆盖现代模型与早期结构学习两个主题；只将下列原始来源用于立项约束，不声称穷尽文献：

- [ANNEVO](https://www.nature.com/articles/s41592-026-03036-7)：已有远程上下文和MoE/GLM，多物种性能不能靠重复其关键词超越。
- [Tiberius](https://pubmed.ncbi.nlm.nih.gov/39558581/)：已有端到端可微HMM，训练9999bp不等于推理只支持9999bp。其当前多类群能力在正式基线冻结时以实际版本核对。
- [Helixer](https://www.nature.com/articles/s41592-025-02939-1)：已有完整gene输出，亦提出稀疏边界序列方向；A也不能仅凭输出形式宣称首创。
- [SegmentNT](https://www.nature.com/articles/s41592-025-02881-2)：已有基础模型多元素分割和上下文扩展，局部元素指标不等于整链恢复。
- [CRAIG](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030054)、[Conrad](https://pmc.ncbi.nlm.nih.gov/articles/PMC1950907/)：结构化/半马尔可夫判别式gene训练不是新概念。
- [GENERanno1.2B preview模型卡](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview)：本次可读取文本仍未提供可核验完整训练暴露清单。不能用监督未见替代预训练未见。

## 历史隔离

M25R/M26冻结结果和候选上界保留；M27固定先导关闭，复现写作包留作资产，不进入本方法的核心贡献链。Setaria继续封存。旧“等待作者接受复现型定位”已由此次用户明确的新模型目标解除，不能再据此阻塞常规研发。作者、投稿声明和正式投稿仍是后期真实作者决策。

## 学习型搜索的相关工作

2026-10-01补充[基因结构搜索的相关工作与贡献边界](../paper/manuscript/M28_learned_search_related_work_20261001.md)：early update、LaSO、BSO、DAGGER和GeneDecoder已有相关先例。B2仅检验同一新增评分器的搜索位置机制，不将学习前缀或结构图本身称首创；不修改正在执行的冻结合同。

## A备选的合成接口准备

2026-10-01。R5已给出固定自由池的能力边界，A仍是原计划中的备选；本节只做无真实数据的工程准备，不选择A替代当前B2，不启动真实训练或新的DEV性能实验。新增 `src/m28/query_pointer.py` 是gene-query存在性、半开坐标边界指针和结束标记的张量接口，不是完整FASTA注释调用器。

每窗输入有效长度内的冻结表示，新query状态不跨窗口或optimizer更新缓存。query数由调用者显式给出，本轮只用2个合成query及hidden4/8，不冻结真实模型容量。正query到基因的匹配由调用者明确提供；不存在的query不自动标负，标签-1保护未知项。存在性BCE、边界CE及结束BCE分别返回，尚未决定真实拟合系数或匹配政策。greedy接口不接收参考，EOS输出名为terminated，不冒称ORF合法或完整CDS；截断或窗口耗尽不补造结束标记。

验证入口 `sbatch/M28-A-POINTER-SYNTHETIC.sbatch` 只运行9项新合成测试：坐标L末端、多exon事件、截断与EOS、padding、冻结输入/新头梯度、未知保护、空参考背景无伪造目标、标签/坐标错配、多query非互斥。使用既有generanno环境，无GPU、1CPU、4GiB，每次墙时不超过5分钟且本准备累计分配最多300CPU秒（0.083333CPUh），新产物不超过1MiB；不安装依赖，不读FASTA/GFF/缓存/checkpoint，不执行optimizer。确定性的源码错误允许在累计限内修一次并重新核验；超时/OOM或无法在限内修复则保留失败并报告，不延时或改合成目标。目录 `outputs/M28-A-POINTER-SYNTHETIC` 必须为新目录，日志按新job ID隔离。

真实query匹配、query数、CDS/剪接及跨exon codon约束、完整无参考调用器、正式数据/更新数和资源协议仍未完成或冻结。当前B2的模型、权重、数据、预算、三阶段作业及继续条件不变；本工程准备成本单列，不能当作A的科学验证或原B2正式预算的放宽。

实际提交命令为 `sbatch --parsable sbatch/M28-A-POINTER-SYNTHETIC.sbatch`，返回job13316149。该job及batch均COMPLETED/0:0，outputs/M28-A-POINTER-SYNTHETIC/STATUS为COMPLETED；9项新测试全部通过，测试循环0.552秒。没有重复旧测试或发生重提。

Slurm root实际分配24秒、AllocCPUS=1，无GPU，即0.006667 allocated CPUh；不能用0.552秒替代分配成本。batch MaxRSS为435060K；本次STATUS、JOBID及两份Slurm日志共1084 bytes，低于事前1MiB运行产物限。最终两份源码及sbatch均逐字符回读匹配，Python AST和bash语法通过。这仅证明合成接口及新头梯度可运行，未执行optimizer、读取真实特征/参考或得到A的注释成绩。

## A备选的DNA语法接线（工程阶段）

2026-10-01。追加 `src/m28/query_decode.py`：输入仅为同方向DNA和有效冻结特征，沿用B1的连续ATG起始、GT/GC–AG剪接、≥5bp intron、实际未完成密码子carry、编码区无未知碱基/内部STOP及连续末端TAA/TAG/TGA。首外显子和末端STOP外显子至少3bp；不支持跨剪接拆开的起止密码子、非标准遗传密码或非规范剪接。以后评价不能据此删除不支持参考或宣称广泛适用。

每个query直接读取全部DNA合法起点、acceptor和上下文合法末端，不经过B1的每bin/endpoint/link截断候选池；仍是每步单一greedy历史，不保证全局最优、召回或正确gene。末端按boundary log-probability与prospective post-end EOS/continue log-probability联合选择：STOP与DONOR受DNA语法限制，唯一合法STOP可在原始EOS为负时被选；这不是旧裸指针的EOS>0规则。达到exon上限、不存在合法acceptor/末端时保留unterminated，不能回退强行挑其他STOP。返回所有query原始logit及链，不冻结存在性阈值、不去重或匹配参考，不冒称已实现完整FASTA调用器。

事前验证合同：`sbatch/M28-A-QUERY-GRAMMAR.sbatch`仅运行10项新的合成DNA语法测试，覆盖GT/GC、相位相同但实际carry不同的junction STOP、最短intron、exon截断、编码/内含子N、明确不支持的拆分起止、softmask/padding长度对齐、正反链半开坐标以及无候选数截断的boundary+EOS联合选择。只导入先前fixture，不重跑未变化的9项裸指针测试。新目录`outputs/M28-A-QUERY-GRAMMAR`、按job分离日志；generanno环境，1CPU/4GiB/0GPU，单次≤4分钟且本追加阶段累计≤240CPU秒。与前阶段实际24CPU秒合计≤264CPU秒，仍受原准备累计300CPU秒及1MiB运行产物限制。确定性源码错误只允许累计限内一次修复重验；不处理超时/OOM为自动重提理由。不得读取FASTA/GFF/特征缓存/checkpoint、更新optimizer、安装依赖或改变B2冻结协议。

实际提交 `sbatch --parsable sbatch/M28-A-QUERY-GRAMMAR.sbatch`，job13316408及batch均COMPLETED/0:0，STATUS=COMPLETED。10项新测试全部通过，测试循环0.480秒；没有修复、重复旧9项测试或重提。该验证证明合成窗口的合法续接/拒绝、联合末端选择和接口坐标行为，不证明真实gene准确率或A优于B2/C0。

Slurm root分配25秒×1CPU=0.006944 CPUh，无GPU，batch MaxRSS=431036K。本阶段STATUS/JOBID及两日志1193 bytes；与前阶段累计49CPU秒=0.013611 CPUh、运行产物2277 bytes，均低于事前限。两份新源码和sbatch提交前逐字符回读匹配，Python AST及bash语法通过；旧裸指针只更新模块说明以指向新语法接口，其计算实现未改。B2仍按原合同运行和由原固定监控对话负责等待。

下一科学决策仍依据B2完成的公平DEV结果；A尚缺真实query匹配/容量与标签政策、正式first-fit协议、存在性操作点、全基因组ownership/去重以及无参考端到端评估，不把本合成PASS升级为科研PASS或论文主张。

## A备选的query匹配原语（工程阶段）

2026-10-01。追加 `src/m28/query_training.py`，仅训练侧使用显式提供的完整正链。每个正链与每个query用原张量接口重新teacher-force形成存在性、平均边界CE及平均EOS BCE的成对成本，再作矩形线性分配；分配离散且不反传，匹配后重新计算可反传损失。全部正链必须一对一分配，query数不足直接报错，不能按分数删基因。标准线性分配不是新颖性主张，也不读取基因组、参考文件、缓存或checkpoint。

未匹配query的标签必须由调用方显式给出0或-1，损失系数三个数也无默认；本原语不把缺失/partial/其他isoform自动转负。若unmatched=0，匹配成本存在性项是positive BCE减negative BCE，再除以总query数；若unmatched=-1，只用positive BCE再除以正链数。边界和EOS成本各除以正链数，使分配最小化当前teacher_components的完整加权窗口目标，而不是仅最小化漏掉null项的成对损失。空正链保留明确的负标签或未知政策，边界/EOS为可微零，无伪造链。

R4 build_targets的positive_chains、protected_chains、unknown_chain_spans、callable_bases与reference_empty_window仍是后续政策的事实来源。B1针对具体候选跨度判负的规则不能无说明地等同于A对所有未匹配query判负；有partial/alternative/不可调用区域时如何保护，以及真实Q/损失系数/eligibility和解码语法对齐，都须在A正式first-fit前固定。本阶段只验证匹配机制，不接入真实TRAIN或自动选定上述政策，不能当作A真实训练已获批准或B2结果已完成。

事前验证合同：`sbatch/M28-A-QUERY-MATCHING.sbatch`只跑8项新合成测试，含小规模穷举完整损失最小值、null基线/归一一致性、参考顺序、query容量、未知项梯度保护/冻结输入、空参考、padding/边界错配及非有限成本停止；不重复既有9+10项测试。既有generanno已确认安装scipy，不新增依赖；仅1CPU/4GiB/0GPU、单次≤2分钟，本匹配阶段累计≤120CPU秒。前两阶段实际49CPU秒，本阶段最大后累计≤169CPU秒，仍遵守原A准备累计≤300CPU秒和新运行产物≤1MiB；原准备仅一次确定性源码修复机会仍未使用，超时/OOM不自动重提。新目录`outputs/M28-A-QUERY-MATCHING`与job隔离日志；不执行optimizer、不安装、不读FASTA/GFF/真实特征、不更改B2。

实际提交 `sbatch --parsable sbatch/M28-A-QUERY-MATCHING.sbatch`，job13316565及batch均COMPLETED/0:0，STATUS=COMPLETED。8项新测试全部通过，测试循环0.456秒；没有失败、修复或重提，也未重复旧9+10项测试。小规模穷举核对了当前完整窗口目标的最低损失，未知query行梯度为零且合成输入不反传。坐标匹配测试与上一阶段DNA语法测试分开，未接入真实R4标签或证明联合模型训练/自由推理性能。

Slurm root实际19秒×1CPU=0.005278 CPUh、0GPU，batch MaxRSS=493480K；本阶段STATUS/JOBID及两日志1049 bytes。A准备三阶段累计68CPU秒=0.018889 CPUh、运行产物3326 bytes，低于事前300CPU秒/1MiB限；本阶段亦低于120CPU秒限。两份新源码和sbatch提交前逐字符回读匹配，Python AST及bash语法通过；旧裸指针只更新模块说明，计算代码未改。没有optimizer更新、依赖安装或真实FASTA/GFF/缓存/checkpoint读取，B2及原固定监控安排不变。

A真实first-fit仍未启动/批准或冻结。下一步需根据完整B2公平DEV结果确定路线，并在正式A协议中明确Q与正/负类比例、partial/alternative保护、匹配/训练系数、共享底座来源/训练曝光、吞吐预算与无参考全基因组输出；Q变化会改变现有存在性BCE的类质量和归一，不能只把它当不影响实验的容量参数。本合成PASS不能代替独立测试、对标收益或论文主张。

## A备选的TRAIN元数据容量与保护盘点（事前范围）

2026-10-01。本节从前三个纯合成接口阶段进入仅真实TRAIN元数据的准备盘点，不是A首训。只读现有R4 `outputs/M28-LABEL-POLICY-R4/{arabidopsis_thaliana,oryza_sativa}/windows_train.jsonl`、对应catalog中按split=train预筛后所需transcript记录，及原R2 `paired_training_draws.jsonl`；TRAIN染色体allowlist完全沿用B2。不得读DEV参考记录、val窗口、test/Setaria、FASTA/GFF、特征缓存或checkpoint，不执行模型前向/反传、optimizer或性能评价。

新增 `scripts/experiments/M28-METHOD-RESTART/query_training_census.py`区分自然TRAIN窗口、原768/物种有放回draw及其唯一窗口。核对窗口/正链ID粒度、catalog连接、primary资格/FLANK和原draw stratum/weight；失败时停止盘点，不静默删目标或修R4标签。报告完整正链数的最大值/分布、exon规模/短intron、未知跨度与完整受保护非正链，以及参考空窗元数据。容量最大值仅是当前TRAIN正目标下界，不直接选Q；保护率仅是元数据下界，未检查全窗callable DNA，不能把剩余窗口直接定为已知负例或宣称泛化容量。

`sbatch/M28-A-QUERY-TRAIN-CENSUS.sbatch`使用既有generanno、private分区1CPU/4GiB/0GPU，单次与本阶段累计≤120CPU秒；前三阶段实际68CPU秒，本次最大后累计≤188CPU秒，仍受A准备总计≤300CPU秒及运行产物≤1MiB限。新目录`outputs/M28-A-QUERY-TRAIN-CENSUS`只保存compact summary/STATUS/JOBID，日志按job隔离；不重复旧27项测试，不改B2或原监控。原确定性源码错误一次修复机会仍未使用，若发生则只能在累计限内处理；超时/OOM/输入身份或科学限制冲突不自动重提或改数据。

实际提交 `sbatch --parsable sbatch/M28-A-QUERY-TRAIN-CENSUS.sbatch`，job13316799及batch均COMPLETED/0:0，STATUS=COMPLETED。窗口ID粒度、TRAIN allowlist、catalog连接、正链资格/FLANK与原draw stratum/weight检查均通过；未改写任何输入，未解析DEV参考或执行模型。详细计数保留于`outputs/M28-A-QUERY-TRAIN-CENSUS/summary.json`及盘点脚本，不额外生成notebook、托管报告或公开分享。

| 物种 | 自然TRAIN窗口 | 每窗正目标最大值 | 原768次draw最大值 | 自然窗元数据保护 | 原draw元数据保护 |
|---|---:|---:|---:|---:|---:|
| 拟南芥 | 12,548 | 11 | 10 | 8,743/12,548（69.68%） | 567/768（73.83%） |
| 水稻 | 12,964 | 8 | 7 | 3,502/12,964（27.01%） | 237/768（30.86%） |

若要容纳当前全部自然TRAIN正目标，Q<11必不足；原抽样最大值10漏掉三个11目标窗口。这只是TRAIN容量下界，不选定真实Q、不给DEV/其他物种容量保证，也不证明模型可以恢复这些链。保护率是unknown_chain_intervals或完整protected非正链的窗口并集；拟南芥尤其不能默认所有未匹配query负。未读DNA，保护率仍是下界；无上述标记不自动等于callable已知负窗，reference-empty也不是生物学真负证据。检查没有暴露身份/连接质量失败，未进行不相关的时间漂移分析或修标签。

1,536条原draw对应1,491个唯一窗口（拟南芥745、水稻746），包含3,028次完整正链曝光但仅2,813个不同primary转录本ID（2,093+720）；有放回曝光不称独立样本。全部自然TRAIN则是40,023次正链窗口曝光、22,407个不同primary ID。正链最多61/51个exon，最短intron10/31bp；这批元数据正目标没有>128 exon或<5bp intron，不是重新验证DNA ORF的结论。

当前标量匹配实现若取容量下界Q=11，一个原draw pass需要33,308次gene-query teacher pair，不能直接借用B2的每窗用时估计A预算。合理的下一工程步是把同一参考链对所有query的teacher成本向量化并核对目标等价，再作单独吞吐准备；不据此启动A真实首训、重采样或选择操作点。

Slurm root分配7秒×1CPU=0.001944 CPUh、0GPU，batch MaxRSS=301504K；本阶段summary/STATUS/JOBID及两日志6493 bytes。A准备四阶段累计75CPU秒=0.020833 CPUh、运行产物9819 bytes，仍低于300CPU秒/1MiB总限，本阶段低于120CPU秒限。源码和sbatch提交前逐字符回读匹配，AST/bash语法通过；无作业失败、修复或重提。B2冻结合同、主评价分母与原固定监控不变，最终路线仍等完整公平DEV，不把这份容量/标签盘点当性能或论文证据。

## A备选的匹配成本批量化（事前工程范围）

2026-10-01。只将同一完整参考链对全部query的no-grad成本计算批量化，减少逐query重复的pointer/GRU调度。保留标量pointer/advance的原计算路径；二维query状态使用共享坐标keys和批量GRU。teacher_chain_components使用原始teacher几何、逐链的边界/EOS事件均值，不混入DNA解码器的最短intron约束。匹配仍逐正链组装矩形成本，unmatched=0的positive-minus-negative/Q、unmatched=-1的positive/N与边界/EOS各除N不变；分配后仍用原teacher重新构建可反传损失。没有缓存跨optimizer状态、混合未知query政策、默认系数或真实Q选择。

事前验证：重跑受影响的8项匹配测试，新增3项标量/批量pointer与GRU、单/多/长链事件均值及padding测试，不重复不受影响的9项裸指针和10项DNA测试。小张量loss/logit绝对容差1e-6、GRU状态2e-6（rtol=0）；完整窗口分量最大绝对差须≤1e-5，失败不放宽容差。单次合成CPU测量固定seed0、L=24576/H128/Q11、61 exon人工链且最后end=L，随机未训练张量，先标量后批量各一次。Q11仅来自TRAIN容量下界的尺寸准备，不冻结容量。只计teacher成本时间，共享window初始化、底座、线性分配、grad readout/反传、optimizer、拼接和I/O均排除；标量仍含旧存在性BCE/标签检查。单次固定顺序且前置测试可能暖机，不作稳健加速比、GPU吞吐、注释性能或论文效率主张。

新脚本time_query_teacher.py与sbatch/M28-A-QUERY-BATCHED.sbatch使用既有generanno、private分区1CPU/4GiB/0GPU，单次及本追加阶段累计≤180CPU秒；此前实际75CPU秒，最坏累计≤255CPU秒，仍受准备总计300CPU秒/1MiB运行产物限。新outputs/M28-A-QUERY-BATCHED目录及job隔离日志。确定性源码错误最多使用原未用过的一次修复机会且必须守上述累计限；超时/OOM/容差失败不自动重提、放宽或改变尺寸。不读取真实DNA/参考/特征/cache/checkpoint，不执行optimizer、依赖安装或新的DEV性能实验，不改B2及原固定监控。真实A首训/标签政策/公平比较协议仍未批准或冻结。

实际提交 `sbatch --parsable sbatch/M28-A-QUERY-BATCHED.sbatch`，job13318020及batch均COMPLETED/0:0，STATUS=COMPLETED。受影响的8项匹配测试与3项新批量测试全部通过（循环0.888/0.191秒）。完整窗口分量最大绝对差：boundary CE=9.5367431640625e-7、EOS BCE=5.960464477539063e-8，低于固定1e-5；既有完整加权窗口目标、未知项梯度保护及容量停止仍通过。未重复旧9+10项未受影响测试，无失败、修复或重提。

单次固定CPU合成teacher成本为标量2.689103秒、批量0.545663秒，标量/批量=4.928140；这只是此人工长链fixture的工程诊断，既不估计真实GPU训练预算，也不证明完整FASTA注释提速或科学性能。详见outputs/M28-A-QUERY-BATCHED/timing.json及固定脚本。

Slurm root分配32秒×1CPU=0.008889 CPUh、0GPU，batch MaxRSS=587872K；本阶段timing/STATUS/JOBID及两日志3337 bytes。A准备累计107CPU秒=0.029722 CPUh、运行产物13156 bytes，仍低于原300CPU秒/1MiB总限和本阶段180CPU秒限。新源码/sbatch提交前逐字符回读匹配，pointer两处改动及事前范围已回读，AST/bash语法通过。没有optimizer、真实数据模型拟合或DEV新实验，B2冻结设置不变。

## B2全范围判定后的路线状态

2026-10-01。完整B2三阶段及逐参考链诊断已完成，详见[正式结果与诊断](../archive/m28_pre_full_grammar/reports/b2_execution.md)。8690窗late精确重放原控制、共同gain一致、7728主分母完整，数值与阶段计数可用；当前early前移版本未通过预设继续条件。拟南芥/水稻exact F1由late26.84/44.96%降为early21.03/20.22%，背景span覆盖同时大幅增加。拟南芥召回净改善不能替代双物种条件，当前B2不继续扩训/改阈值救结论。

逐链结果显示early最终新增1145条、同时丢掉late的1042条，其中909条已不在early自由池；不是仅最终筛选阈值的变化，也不证明某个未经检验的单一根因。A已完成合成指针、DNA语法、匹配、TRAIN容量/保护盘点及成本批量化，仍非完整无参考全FASTA调用器或经真实训练的模型。下一安全工作是形成具体A首训提案：显式Q及标签保护/负类质量、loss和匹配归一、共享来源/曝光可比性、固定操作点、ownership/去重、真实吞吐预算及不支持参考保留。未经冻结和明确批准不执行A首训；不因为B2失败就改变发表目标或开放test/Setaria。

本轮B2正式加诊断实际5.374722GPUh、0.27独立CPUh、含本cohort日志约1.231GiB；历史成本另留。论文/SOTA/泛化目标仍未达到。固定轻量监控修正后14:59UTC在原对话完成新鲜终态/产物检查并自动回传一次，已验证。本机心跳内置完成指令要求delete，收尾删掉自动化而非按约定暂停；主对话已通过正式工具恢复原ID、原对话绑定及15分钟间隔，状态PAUSED并重建本cohort交接记录。未来自动pause保留尚不能保证，详情与证据见b2_execution.md，不新增对话或重复交接。

## A query输出与公共合并接口的工程准备

2026-10-01。本追加仅将现有无参考DNA query输出接到既有assembly/export路径，使用人工DNA及未训练小张量；不选择真实存在性阈值、不改损失/匹配/梯度政策、不读取FASTA/GFF/缓存/checkpoint或执行optimizer。新增query_output.py要求调用者显式给出有限logit操作点τ，不设真实默认。所有terminated query及其logit−τ，包括非正分和同窗重复链，进入标准free记录；未结束链仅留诊断。公共合并仍按最近中点owner、同链最大分、最大正增益同链向非重叠选择，零分选null，分数不能与C0/B1数值横向比较。

新test_m28_query_output.py仅核验6个新跨接口失败：合成多外显子query→正反链全局坐标/GFF3及phase回读，显式操作点/零分null且free记录不偷删，截断不冒称完整gene，无ATG空输出，nonowner-only丢失与同窗query重复分开，以及非有限操作点/存在性logit停止。只复用已有fixture，不重跑旧9+10+8+3项未改测试。新入口sbatch/M28-A-QUERY-OUTPUT.sbatch为generanno环境、1CPU/4GiB/0GPU、单次与本阶段累计≤120CPU秒；此前107CPU秒，最坏累计≤227CPU秒，仍守A准备总300CPU秒/1MiB运行产物限。输出目录必须新建、日志按job隔离；原一次确定性源码修复机会未使用，若需修复也受同一累计限，超时/OOM不自动重提。不变更B2历史结果、test/Setaria或原已暂停监控。

首训提案见[a_first_fit_proposal.md](../archive/m28_pre_full_grammar/reports/a_first_fit_proposal.md)，它是待批准科学选择而非既有冻结合同。当前pointer的detach与未均衡存在性目标不自动等于提案中的联合梯度/类均衡；未经既有批准门槛不执行新的真实A首训。

实际工程提交为sbatch --parsable sbatch/M28-A-QUERY-OUTPUT.sbatch，job13319189及batch均COMPLETED/0:0，STATUS=COMPLETED，6项新测试全部通过（0.433秒）。中断的前一次提交已用当前squeue、当日sacct、目标目录及本机提交进程核对没有创建job，再作本次唯一实际提交；没有已启动作业失败或重提。没有GPU、真实DNA/标签/cache/checkpoint读取或optimizer；未重复旧9+10+8+3项未改变测试。

原始日志报告一次ResourceWarning：既有eval_structure_diagnostic.py::parse_annotation使用open(path)未显式关闭句柄。人工GFF3的正反链坐标及phase回读均通过，警告不是模型失败；已保留原日志，不为该警告改动冻结评价器或重复整套验证。

Slurm root实际6秒×1CPU=0.001667 allocated CPUh，AllocTRES没有GPU，batch MaxRSS472072KiB。本阶段STATUS/JOBID及两日志1065 bytes；A准备累计113CPU秒=0.031389CPUh、14221 bytes，仍低于原300CPU秒/1MiB限，本阶段低于120CPU秒限。提交前新源码/提案与范围逐字符回读一致，Python AST及bash语法通过。完整原始FASTA调用器、新联合梯度/类均衡训练入口及真实首训仍未完成，首训提案保持待批准。

## A联合梯度与类均衡匹配的工程准备

2026-10-01。继续实现待批准A提案的共享上下文/query联合梯度与类均衡目标，但不启动真实首训。query_pointer.py仅将detach移到旧GeneQueryPointer.for_window的外部输入入口，旧teacher/greedy接口仍冻结原features且数值路径不变；底层QueryPointerWindow保留调用者给出的图。新增query_model.py复用SharedContext及相同C0共同参数构造顺序，再接独立query头；外部GLM/DNA输入固定，新shared可接收query梯度，没有stop-gradient flag或新backbone。

新balanced_query_components明确使用已观察正/负类各半质量的存在性BCE；成本中positive和null基线归一一致，边界/EOS沿用原teacher事件/链均值。已有普通匹配函数和B2实现不改。保守unmatched_query_policy按unknown跨度、完整protected非positive链及有效DNA noncallable保护整窗，padding不参与，正目标保留。joint_query_losses把已有local_loss与上述三项按显式四系数及原window weight相加；没有默认真实Q/系数/阈值、optimizer、抽样或真实数据访问。

本阶段只运行8项新的合成测试：C0共同初始化/固定宽度padding、query-only梯度至shared但不至缓存/未知query、N1/2/3两种unmatched的穷举最低完整损失、每个分配的成本加null常数等于完整目标、空窗负/未知可微零、三种保护原因且不删正例、local+query与importance weight、容量/坐标/非有限停止。因detach入口移动，重跑受影响9项裸接口和8项普通匹配，不重复未变的10项DNA、3项批量或6项输出测试。固定绝对容差1e-6，不放宽。只构造一次未训练D2048/Q16/H128模型统计参数，尺寸来自待批准提案，不据此选容量或称训练。

入口sbatch/M28-A-QUERY-JOINT.sbatch仅用generanno、private 1CPU/4GiB/0GPU、单次与本阶段累计≤120CPU秒；此前113CPU秒，最坏累计≤233，仍守原准备300CPU秒及1MiB运行产物限。新输出目录/按job日志，不读真实FASTA/GFF、TRAIN/DEV元数据、cache或checkpoint，不执行optimizer，不安装依赖、不动test/Setaria、B2历史或已暂停监控。原一次确定性源码修复机会若需要，亦计入相同阶段/总限；超时/OOM或无法守限则留证停止。首训资源与科研选择仍待人类批准，纯合成PASS不证明实际训练/无参考全FASTA或性能。

job13320418及batch FAILED/1:0，STATUS=FAILED exit=1；9项旧裸接口和8项普通匹配均通过，新联合测试7/8通过。失败为完整成本加null常数的断言：4.306426048278809与4.306425094604492相差9.5367431640625e-7；事前合同是绝对差≤1e-6，而测试源码误用places=6（round(diff,6)==0），两者并非同一规则。已读取unittest实际实现确认此语义错误，不改模型、目标、fixture或声明的容差。使用原唯一一次确定性源码修复机会，将本新测试中的三个places=6断言明确改为delta=1e-6；失败原目录与日志完整保留，不把首轮记为PASS。

首轮实际24CPU秒、0GPU、batch MaxRSS535312KiB，STATUS/JOBID及两日志3793 bytes；准备累计137CPU秒及18014 bytes。修复后入口仍为同一sbatch/M28-A-QUERY-JOINT.sbatch，只复验8项新联合测试并完成此前尚未执行的未训练参数统计，不重复已通过的17项旧接口；新目录outputs/M28-A-QUERY-JOINT-R2，剩余墙时硬限96秒×1CPU，首轮加本次最坏≤120CPU秒，准备总最坏仍≤233CPU秒/1MiB。一次修复机会现已使用，后续失败留证停止；真实首训/计时未批准，也不读取真实数据或执行optimizer。

修复后实际提交sbatch --parsable sbatch/M28-A-QUERY-JOINT.sbatch，job13320471及batch COMPLETED/0:0，STATUS=COMPLETED，8项联合测试全部通过（0.352秒），没有再次修复或重复已通过17项旧接口。合成证据支持query-only梯度进入新shared而不进入外部缓存/未知query、共同参数与C0初始化一致、类均衡矩形匹配及null常数与完整目标一致、整窗保护不删正例、原importance weight与局部监督仍接通；不支持真实训练/完整FASTA调用或性能结论。未训练D2048/Q16/H128模型实际583153可训练参数，optimizer更新0。

Slurm root本次实际11秒×1CPU=0.003056 allocated CPUh、AllocTRES无GPU，batch MaxRSS533176KiB，本次STATUS/JOBID及两日志1193 bytes。首轮加复验合计35CPU秒=0.009722CPUh、4986 bytes；A全部准备累计148CPU秒=0.041111CPUh、19207 bytes，均在本阶段120CPU秒与总300CPU秒/1MiB内。调度限制需如实区分：脚本请求00:01:36，但Slurm实际保存Timelimit=00:02:00，故不能声称96秒剩余预算由调度墙时强制执行；本次实际11秒未越限。两轮原日志和失败目录保留，一次修复机会已使用，不再自动重提。本轮没有GPU、真实数据/模型拟合或封存数据访问；已将提案中的“联合梯度/类均衡尚未实现”改为合成已核验、真实训练及全FASTA接线仍待完成，首训提案仍待批准。

## A单轮计时训练与DEV推理入口的工程准备

2026-10-01。新增run_query.py的prepare/fit/infer入口，按现有待批准提案实现固定Q16/H128/seed0、四loss系数1、原1536 draws三遍/4608更新与最后checkpoint。只读R4原cache，真实训练按R4窗口与TRAIN预筛catalog重建目标，不读缓存旧标签，不继承C0/B1/B2监督权重；prepare按每物种四类原draw选定计时窗，同分最小ordinal、缺类别停止，optimizer更新0。fit必须有完成的同配置prepare及两个成本gate；推理只读末A checkpoint、独立完整DEV grid/DNA特征，不调用TRAIN/catalog/评价器，保存全部原始query与完成链供既有assembly/export后续处理。代码存在不意味着首训获得批准；本轮不调用任一真实数据CLI，不创建GPU作业或真实optimizer。原始FASTA编码CLI、真实运行与A统一评价仍未交付。

新增7项合成测试核验新的输入接线失败：原draw/TRAIN split/stratum/weight/positive资格，cache几何与旧标签隔离，固定计时类别及同分ordinal/缺失停止，两个预算gate及prepare身份，最后A checkpoint身份，R4重建目标进入新shared→loss/backward，及新shared→无参考query输出/全query序列化。只用内存人工metadata、零/未训练小有效长张量；不读取真实FASTA/GFF、TRAIN/DEV元数据、特征缓存或checkpoint，不执行optimizer，不重跑旧未变测试，不以合成PASS替代真实CUDA/训练/全8690窗。固定科学配置来自待批准提案，本轮不选参或改变资源协议。

sbatch/M28-A-QUERY-RUNNER-TEST.sbatch为generanno、private 1CPU/4GiB/0GPU、单次120秒；此前准备实际148CPU秒，最坏累计268CPU秒，仍守300CPU秒/1MiB运行产物限。只新建outputs/M28-A-QUERY-RUNNER-TEST及按job日志。原一次源码修复机会已使用，任何本次失败均保留证据，不自动修复重提/放宽测试或加预算。新的生产入口将来仍须先有人类首训批准和实时Slurm检查；本轮不改B2历史、test/Setaria或已暂停监控。

实际唯一提交sbatch --parsable sbatch/M28-A-QUERY-RUNNER-TEST.sbatch，job13321036及batch COMPLETED/0:0、STATUS=COMPLETED，7项新合成测试全通过（1.108秒）。源码/sbatch及事前范围逐字符回读匹配，Python AST与bash语法通过；没有失败、修复或重提，不重跑旧未改测试，也未调用任何生产prepare/fit/infer模式。该PASS覆盖新helper接线，并非真实TRAIN/CUDA计时、4608更新、8690窗或FASTA端到端验证。

Slurm root实际45秒×1CPU=0.012500 allocated CPUh、AllocTRES无GPU，batch MaxRSS976864KiB；STATUS/JOBID及两日志918 bytes。A准备累计193CPU秒=0.053611CPUh、20125 bytes，低于300CPU秒/1MiB；本阶段低于120秒。一次源码修复机会此前已用完，本次没有追加机会。已更新现有首训提案的代码状态，保留原Q/损失/标签/数据曝光/阈值、10.166667GPUh/4CPUh/2GiB预算与科学继续/停止条件。计时gate沿用提案的最大完整TRAIN load→backward与decode分量投影，后者不是包含load/forward的全推理时间保证；生产阶段仍须守原硬限。真实首训/计时及后续评价未批准或提交；原始FASTA编码CLI和A统一评价入口仍待补齐，SOTA/泛化/期刊目标未达成。

## A冻结评价及原始FASTA调用代码的工程准备

2026-10-02。新增evaluate_query.py，只在完整同末checkpoint/配置的8690窗、四scope及自然网格验证后读取原DEV参考；沿用common_ruler/evaluate_pilot选择、owner/去重/selector、GFF3和统计函数。主分母5437+2291=7728，旧6450副表单列；不按A训练资格/语法/长度过滤参考。继续条件直接用整数tp/predicted/reference判断两物种F1均严格改善，且同背景分母下span均不增加；没有round后的gate或论文成功判定。保存实际free/owner/selected恢复、全query/未结束计数及CDS结构诊断，GFF回读同时比较ID/坐标/phase。历史SOTA和R5控制仅标缓存，不称原生重跑或因果消融。

新增annotate_query.py从唯一ID的原始FASTA读取长度、大小写统一/非ACGT为N、两链完整网格，按原revision的四块6144bp/six-mer无特殊token编码并BF16物化后送FP32 A头，走同query/DNA语法/owner及GFF3出口。原生模型只在显式main中从本地HF缓存加载，不下载、不读取标签、reference或特征缓存；输入权重限匹配的A末checkpoint。命令要求显式FASTA/checkpoint、新输出目录和输出容量限，记录全部raw query。此代码不是已执行的原生FASTA验证，也不授权新的数据访问或GLM GPU运行；原生数据/资源预算须另行批准，不挪用现有缓存首训10.166667GPUh/4CPUh/2GiB。

仅新增8项合成核验：FASTA多记录/大小写/歧义与重复/空记录停止、四块/token/hidden对齐、人工FASTA经替身encoder及未训练小头→双链record/owner/GFF、完整DEV结构与checkpoint门槛、A记录阈值/Q/provenance拒绝旧B1、双物种整数gate/背景与主分母、模拟参考loader保留unsupported并拒绝变化、GFF phase不能被坐标相同掩盖。只用人工临时FASTA/GFF、内存元数据/小头和替身编码器；真实参考loader在测试中替身化，不读取真实TRAIN/DEV/genome/cache/checkpoint或执行optimizer/GPU，不调用两种生产CLI，不重跑旧未变测试。

sbatch/M28-A-QUERY-CLOSURE-TEST.sbatch用generanno、private 1CPU/4GiB/0GPU，整分钟60秒硬限；此前实际193CPU秒，本次最坏累计253，仍在原300CPU秒/1MiB内。新目录outputs/M28-A-QUERY-CLOSURE-TEST及按job日志；一次源码修复机会已用完，本次任何失败保存后停止，不自动修复重提/扩大预算/放宽断言。不改B2冻结历史、test/Setaria或已暂停监控；首训提案的科学/资源/批准状态不变。

实际唯一提交sbatch --parsable sbatch/M28-A-QUERY-CLOSURE-TEST.sbatch，job13321229及batch均COMPLETED/0:0、STATUS=COMPLETED，8项新合成测试全通过（0.580秒）。Python AST、bash语法和四文件逐字符回读通过。Slurm root实际17秒×1CPU、AllocTRES无GPU，batch MaxRSS828152KiB；STATUS/JOBID与两日志1428 bytes。本阶段使用人工临时FASTA/GFF、替身编码器和参考loader，未调用任一生产CLI、无真实TRAIN/DEV/cache/checkpoint读取或optimizer更新。冻结解析器再次出现已知未关闭文件ResourceWarning，GFF ID/坐标/phase核验通过；不改冻结评价器、不压制警告、不重跑旧测试，没有修复或重提。

A合成准备累计210CPU秒=0.058333CPUh、21553 bytes，守原300CPU秒/1MiB；0GPU、0 optimizer更新。现有首训提案已同步：prepare/fit/infer、原始FASTA调用和冻结评价代码均已具备，真实CUDA/原生encoder、实际训练权重和完整DEV结果仍未验证。A尚无性能或泛化结论，期刊目标未达成。下一实质步骤是原提案的真实TRAIN计时→4608更新一次拟合→8690窗统一评价，需明确首次批准（原10.166667GPUh/4独立CPUh/2GiB）；原生FASTA另行定数据/预算。本次停止于该批准边界，不继续重复合成审计，不开放test/Setaria，也不启用监控或改旧实验。
