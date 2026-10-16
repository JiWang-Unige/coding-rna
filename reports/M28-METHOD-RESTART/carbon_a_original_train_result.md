# Carbon-A：原六 TRAIN 完整结构比较

2026-10-11 UTC的固定比较。原生 Carbon-A 完整流程和两种预先声明的置信度操作点均已完成；它在本面板上的整链成绩高于我们已完成的候选，但不构成独立泛化或解码器因果证据。后续完成状态于2026-10-13更新：B fresh×joint原TRAIN结果见[独立四格记录](full_grammar_fresh_joint24_result.md)；[C uniform-coarse](full_grammar_uniform_coarse_result.md)、[D codon-energy](full_grammar_codon_result.md)、[E rank4 q/v](full_grammar_joint_lora_result.md)及[M dense q/v](full_grammar_dense_qv_result.md)各自固定TRAIN闭环均已完成。下面保留原固定Carbon-A比较及其观察口径，不回写旧家族或新学习率对照结果。

## 完成证据与共同尺

评价作业 13503419 的 root／batch／extern 均为 COMPLETED 0:0，实际 UTC 23:33:10–23:42:57，9分47秒；[评价收据](../../outputs/M28-CARBON-A-TRAIN-EVAL/summary.json)和 STATUS 为完成状态。六条原 TRAIN 染色体共156,779,632 bp，保留22,688条 primary参考（拟南芥17,739、水稻4,949）。坐标、strand和完整有序CDS区间须全部相同才计为primary TP；使用唯一结构key，而非GFF记录数。逐染色体、分物种和pooled的TP／预测／参考及PRF直接重算一致，species macro和worst也一致。unsupported参考保留，未按训练资格筛除，推理未使用参考。

固定revision为f6789ee79e6ef397de333577ffec1cb734496620；标准原生上下文98,304 bp、重叠6,144 bp、BF16／FlashAttention-2、batch8、dispatch4、GPU Viterbi及其余作者默认约束不变。两种操作点为confidence≥0.1或NA，以及未过滤；没有根据观察成绩选阈值。未使用DEV或sealed数据，更新0，冻结DEV家族未改。

五条原raw作业成功。NC_003071.7的原13503402仍为FAILED 1:0：原生概率和4,172条基因／蛋白输出已完成，独立validator误将GCN歧义密码子推断成A，而作者字面标准表写X。只修自己的翻译校验；CPU恢复13504510的root／batch／extern均COMPLETED 0:0，逐字DNA、双链概率、全部基因的坐标／CDS／phase／蛋白检查通过。原预测与失败记录未覆盖，没有重跑GPU；其接受收据位于独立恢复目录。原wrapper计时缺失保持null，不能用CPU恢复时间代替推理时间。执行与来源细节见[准备及恢复记录](carbon_a_pipeline_readiness.md)。

## 完整链结果

表中均为F1百分数。旧工具来源为[原生基线评价](../../outputs/M28-NATIVE-TRAIN-BASELINE-EVAL/summary.json)，初始及joint来源分别为[4608步诊断](../../outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC/summary.json)和[joint诊断](../../outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC/summary.json)；不回写原历史家族。

| 方法 | 拟南芥 | 水稻 | Pooled | Species macro | Worst species |
| --- | ---: | ---: | ---: | ---: | ---: |
| Carbon-A：confidence≥0.1或NA | 78.958 | 75.259 | 78.172 | 77.108 | 75.259 |
| Carbon-A：未过滤 | 76.745 | 73.530 | 76.063 | 75.138 | 73.530 |
| ANNEVO | 76.831 | 71.179 | 75.542 | 74.005 | 71.179 |
| Tiberius | 75.417 | 63.577 | 72.655 | 69.497 | 63.577 |
| Helixer | 68.122 | 50.294 | 63.592 | 59.208 | 50.294 |
| 已完成joint 24,576 bp候选 | 70.614 | 53.984 | 66.963 | 62.299 | 53.984 |
| 初始grammar：4608步 | 50.487 | 28.278 | 45.590 | 39.383 | 28.278 |

过滤Carbon-A匹配16,896条primary／20,540条唯一预测，pooled precision82.259%、recall74.471%。未过滤匹配17,023／22,072，precision77.125%、recall75.031%。两者均无重复结构记录。过滤口径比ANNEVO／Tiberius／Helixer分别高2.630／5.516／14.580个百分点pooled，且两个物种均高于三工具；未过滤口径比ANNEVO pooled高0.522点，但拟南芥低0.086点，不能称每物种均胜。

过滤Carbon-A比已完成joint候选pooled高11.209点，primary命中多1,970、唯一预测少1,352；未过滤也高9.101点。该比较说明当前完成候选仍未达到此强完整流程对照，而不是证明某一模块导致差距。

## 固定过滤的得失

只读算术诊断13504740的root／batch／extern均COMPLETED 0:0，UTC23:44:55启动并结束，使用已完成评价JSON而非新预测或标签。其[结果](../../outputs/M28-CARBON-A-CONFIDENCE-DIAGNOSTIC/summary.json)保留每条染色体、每物种与pooled的固定端点。

过滤移除1,532条唯一预测，其中127条原primary命中丢失，余1,405条是非primary-exact的benchmark残差；并非均为生物学假基因。没有重复记录移除。拟南芥移除1,229／丢失95个命中，水稻移除303／丢失32个命中。Pooled召回下降0.560个百分点，而F1上升2.108点；六条染色体的F1均上升，但不能据此调阈值或声称独立统计显著。

按先缩小预测分母、再扣除命中损失的顺序，F1端点差为+2.695688−0.587582=+2.108106个百分点。这只是代数核算，不是可实现的中间模型或两个独立因果效应。过滤还将unstranded CDS-base F1从95.859%降至94.147%，reference-relative背景span比例从4.545%降至1.282%；不同终点仍有取舍。

## 与 joint 候选的逐参考链配对

只读配对作业13504864的root／batch／extern均COMPLETED 0:0，UTC23:54:58–23:55:39，41秒。它复用已完成的原primary参考记录13461637及既有GFF，不读DNA、参考GFF或checkpoint，不执行forward／更新／预测修改。六条染色体上两固定Carbon-A口径与joint的TP、唯一预测和重复数全部复现；每条染色体、每物种及pooled的参考交集与分层加总均一致，详见[配对结果](../../outputs/M28-CARBON-A-JOINT24-PAIRED/summary.json)和[执行程序](../../scripts/experiments/M28-METHOD-RESTART/compare_carbon_a_joint24_train.py)。以下为参考级召回计数，不是分层prediction precision或F1。

| 固定 Carbon-A 口径 | 双方命中 | Carbon-A独有 | joint独有 | 双方均遗漏 | 净primary差 |
| --- | ---: | ---: | ---: | ---: | ---: |
| confidence≥0.1或NA | 14,336 | 2,560 | 590 | 5,202 | +1,970 |
| 未过滤 | 14,377 | 2,646 | 549 | 5,116 | +2,097 |

过滤口径下，拟南芥Carbon-A独有1,608／joint独有453，净差1,155；水稻952／137，净差815。未过滤对应为1,677／427、净差1,250，以及969／122、净差847。Joint独有命中不是新基因发现，也不能仅凭真阳性交集构造一个可用的预测合并器。

CDS段数预先分为1、2–5和≥6；基因跨度预先按24,576 bp分层，没有依据成绩选择切点。

| CDS段数 | 参考数 | joint命中 | 过滤Carbon-A命中 | 净primary差 | 召回差 pp |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 5,395 | 4,317 | 4,410 | +93 | +1.724 |
| 2–5 | 9,898 | 6,130 | 7,034 | +904 | +9.133 |
| ≥6 | 7,395 | 4,479 | 5,452 | +973 | +13.158 |

多CDS段参考贡献1,877／1,970（95.279%）的净差；未过滤时贡献1,970／2,097（93.944%）。水稻单段链中，过滤Carbon-A比joint少2个命中，而2–5段及≥6段分别多405和412；后者召回从45.337%到73.383%。因此多段完整结构是当前差距集中的可观察对象，不是已证明的模块失效原因。

超过24,576 bp的参考只有9条（拟南芥2、水稻7），joint命中0、两种Carbon-A口径均命中1；其余22,679条参考贡献过滤／未过滤净差1,969／2,096。总体缺口不能主要解释为参考基因跨度超过joint编码窗口。这个分层不检验长上下文对短基因评分的作用，也不是whole-contig grammar的硬长度上限。两种分层各自覆盖同一全体参考，其贡献不能互相相加。

此结果将作为B/C/D/E固定终态解释的参考：关注多段链恢复与旧命中损失，并同时保留全输出precision／recall／F1。它不授权改阈值、挑checkpoint、凭参考命中合并预测或扩充冻结DEV家族；当前不足以归因于上下文、transition网络、训练暴露或某个decoder组件。

## 非primary残差与已注释其他可评价 CDS 链

2026-10-12 UTC只读既有评价13503419／13507526／13453780的计数，不新开计算作业、不读DNA／注释／checkpoint、不重跑推理或改变预测。三源均使用同一原六TRAIN、22,688条primary参考和38,152条可评价isoform记录（记录数，不是去重后新评价分母）。作者默认窗口为98,304 bp；旧13503419 JSON未保存context字段，该条件来自其已核实的原生执行协议，而非以缺字段推定短窗。

远端实际common_ruler.score先形成唯一结构集合：P为预测、R为primary、I为同源可评价isoforms，结构key为染色体、strand和有序CDS区间，phase另外核验。已有字段exact_other_assessable_isoform_not_primary_TP定义为|(P∩I)−R|；所以唯一预测可分成primary精确匹配、非primary其他可评价isoform精确匹配、两类均不精确匹配。所有方法重复结构计数为0；各染色体和物种加总与pooled一致。下表只重算P−primary TP及再减其他isoform匹配的残差，没有改原F1。

| 固定条件 | 唯一预测 | Primary TP | 非primary残差 | 其他可评价isoform精确匹配 | 两类均不精确匹配 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Carbon-A 98,304 bp／过滤 | 20540 | 16896 | 3644 | 1365 | 2279 |
| Carbon-A 98,304 bp／未过滤 | 22072 | 17023 | 5049 | 1374 | 3675 |
| Carbon-A 24,576 bp／过滤 | 20369 | 16698 | 3671 | 1361 | 2310 |
| Carbon-A 24,576 bp／未过滤 | 21871 | 16840 | 5031 | 1377 | 3654 |
| 已完成joint 24,576 bp | 21892 | 14926 | 6966 | 1267 | 5699 |

| 固定条件 | AT其他isoform匹配 | AT两类均不精确匹配 | Rice其他isoform匹配 | Rice两类均不精确匹配 |
| --- | ---: | ---: | ---: | ---: |
| Carbon-A 98,304 bp／过滤 | 1205 | 1658 | 160 | 621 |
| Carbon-A 98,304 bp／未过滤 | 1208 | 2789 | 166 | 886 |
| Carbon-A 24,576 bp／过滤 | 1208 | 1620 | 153 | 690 |
| Carbon-A 24,576 bp／未过滤 | 1220 | 2658 | 157 | 996 |
| 已完成joint 24,576 bp | 1159 | 3610 | 108 | 2089 |

默认过滤口径的3,644条非primary残差中1,365条（37.459%）精确匹配其他可评价isoform；未过滤为1,374／5,049（27.213%）。短窗口对应1,361／3,671（37.074%）和1,377／5,031（27.370%）。已完成joint输出有1,267／6,966（18.188%）。这说明primary未命中不等同于已证明的假基因，且同源isoform匹配在原共同尺中并未丢失；它们仍不加算为primary TP。

默认固定过滤的1,532条移除量分解为127 primary命中、9其他可评价isoform匹配及1,396条两类均不精确匹配；短窗口1,502条分解为142／16／1,344。前两项是精确几何匹配资格，最后一项不是生物学假结构判定。集合差算术依赖预定过滤输出确为同一未过滤输出的子集，不能当成独立模型或重选阈值的依据。

这里是全输出唯一预测计数，不是逐gene／逐primary参考的损失归因、同gene恢复数、表达或UTR证据，也不证明完整isoform召回或独立生物学有效性。不能将这些总数直接解释为joint丢失primary的原因；先前[LOCAL逐参考isoform诊断](full_grammar_objective_contrast_contract.md#主链变化中的同基因其他可评价-cds-isoform)是不同模型、不同分母和配对单位，不能借用其结果解释当前joint。现有primary排名、阈值、checkpoint、冻结DEV及B/C/D/E均不变。

来源为[默认评价](../../outputs/M28-CARBON-A-TRAIN-EVAL/summary.json)、[短窗口评价](../../outputs/M28-CARBON-A-CONTEXT24-TRAIN-EVAL/summary.json)、[已完成joint评价](../../outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC/summary.json)，以及远端scripts/experiments/M28-METHOD-RESTART/common_ruler.py的既有score实现。

## 结论边界与下一步

这是本项目重新执行的完整raw-FASTA原生流程，不是作者预计算缓存。主权重共享原b0483c23b6b63787b61a6d3a204a9b517d6ba345来源；上下文、原生CDS heads、transition／confidence辅助网络、结构评分和过滤同时不同，未过滤仍保留learned transition。主模型已包含本面板物种／assembly的上游注释暴露，辅助网络各自训练与验证暴露尚未补齐，见[暴露记录](pretrained_exposure.md)。不能称独立encoder、干净zero-shot、纯decoder消融、独立泛化或匹配速度优势。

继续消费B/C/D/E原范围的固定终态及预先声明的机制诊断。只有原已冻结DEV计划中的候选按该计划接续；E和Carbon-A不自动加入。此结果不放行封存／新物种数据，不改变冻结阈值、公开发布或论文claim边界。研究goal保持未完成。

2026-10-12 00:01 UTC调度读回：四个已提交raw数组13494897／13496013／13497945／13499821的ArrayTaskThrottle均已由3提高为6（分别执行scontrol update JobId=<该数组ID> ArrayTaskThrottle=6）。依据实时Slurm：private两节点各8张RTX3090，已分配5及4张，合计尚有7个GPU槽；没有相关维护预约。仅调整既有六任务数组的并发上限，没有取消／重提／新增任务或修改提交脚本；每task仍1RTX3090／2CPU／16GiB／2小时，原命令、checkpoint参数、afterok核验依赖和评价依赖均不变。四组fit仍运行，raw仍为有效Dependency等待；这不是已启动六路推理或已缩短实际用时的证据。原唯一心跳的B/C/D/E对象和结构检查不变，继续安静等待。
