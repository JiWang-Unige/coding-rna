# A29 实际搜索截断候选的监督资格

2026-10-06，在读取全量截断候选标签前固定。下一版若用搜索竞争训练，不能自动把淘汰正确路径的候选当负例：它可能与另一条primary正链相容，也可能是保护异构体或未决区域。本次只给已固定A28的真实首淘汰cutoff加上原TRAIN注释资格，决定后续监督如何避开错误负例，不训练、不重提A9、不改变原搜索。

## 数据与定义

唯一来源是全量A28 job13392879。只在其STATUS/JOBID及原搜索逐query对照、原候选池、正确前缀分数／几何、自然命中对照和保存记录核验全部成功后执行。完整352原TRAIN窗、545参考／窗配对（221遗漏＋324对照）、4247参考/query profiles保持身份与顺序。遗漏1708个局部通过query全部首淘汰；自然命中对照的返回query不制造cutoff，其其他实际首淘汰另列，不并入遗漏分母。

每条首淘汰只读取A28实际保留池的最后一个cutoff，而不是所有超过正确路径的候选。global入口cutoff是START或NEXT_START；延伸池cutoff是DONOR；完成池和分数界cutoff是完整STOP。分数界的正确目标仍是活跃DONOR前缀，不能把它改为提前STOP。累计分数只识别原事件，不重算排序或选择阈值。

复用原run_query.load_train_inputs、build_targets和query_prefix_value_training.compatible，从同六条TRAIN染色体的原定向DNA／R4 catalog取得primary、保护链、chain_unknown_intervals和callable掩码。缓存只在CPU读取和核对原revision／窗口／DNA；不使用旧R2/R3缓存标签，不调用模型或GLM，不读DEV、封存Setaria或新数据。

本诊断是gene-set层的资格，不是A9的matched-query slot标签。cutoff若是任一eligible primary的正确前缀／下一入口／完整STOP，记1；即使它不是该query所分配的目标也不自动负标。否则，保护链相容、跨unknown chain span或其观察span含noncallable任一成立，记-1；剩余记参考相对的known_nonmatch 0。exact primary优先、部分DONOR与完整STOP区分、prefix span与保护OR条件均按已使用的原规则；所有原因flags保留，重叠原因不相加作互斥归因。没有改写A9原标签或任何冻结训练协议。

记录实际source profile ordinal、window／original candidate／query、cohort／物种／exon类、首淘汰stage、cutoff节点、对应primary indices及保护原因。正确目标本身必须仍与该eligible primary相容，否则停止。保存后从实际A28 profile及原R4对象重算每条标签和reference汇总，防止把query数、重复窗口配对或other-primary当独立gene错误。

## 执行与决定

入口analyze_search_cutoff_labels_train.py，唯一输出outputs/M28-A29-SEARCH-CUTOFF-LABELS，cutoff_annotations.jsonl、references.jsonl545行、summary/STATUS/JOBID。四项测试覆盖其他primary仍为set-positive、DONOR与提前STOP、保护／unknown／noncallable、分数界完整cutoff与活跃目标不混淆。

generanno，private CPU-only／0GPU／2CPU／8GiB／15min。352原缓存只顺序读取已有4096×2048特征文件以获取定向DNA，不产生特征或checkpoint；预计文本<8MiB。参考旧CPU汇总约数十秒，考虑缓存I/O估数分钟并留15min；这不是吞吐benchmark。提交时若A28仍运行使用afterok依赖；若已完成，入口仍显式校验原实际job和产物，避免过期依赖。失败保留，只允许有明确证据且同范围的最小工程修复。

报告遗漏／对照分开、逐物种／stage／exon类的三类资格，以及221遗漏配对中至少一个／全部已观察cutoff为known_nonmatch的数目。它不证明所有领先候选都已知负，也不建立训练PASS阈值。若出现other-primary或保护cutoff，未来loss必须显式处理／跳过，不能训练为生物学假基因；即使cutoff全为0，也尚未实现on-policy采集、rollout或结构损失，不称LaSO/DAgger/LOLS已完成。

完成后结合全量A28和A23/A25原自然读出选择下一独立模型，不通过改seed、beam、阈值或追加epoch挽救旧FAIL。完整raw-FASTA模型、原生ANNEVO/Tiberius/Helixer竞争、独立适用范围、论文级因果消融和创新仍是目标，本诊断不是其替代验收。复用同一固定watcher/heartbeat，仅交接终态或可行动故障；无Git提交、作业取消、公共发布或旧框架恢复。
