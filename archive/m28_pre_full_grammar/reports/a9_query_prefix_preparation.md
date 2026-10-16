# M28 A9 query 条件化前缀价值：零更新准备

2026-10-03。用户已授权在 raw-FASTA 注释目标内继续迭代；本轮只执行 TRAIN-only 工程与监督预检，不以未训练头的自由命中决定性能，也不修改已在运行的 FULL/A6/A7。

## 与已经失败的 B2 的关系

B2 在冻结 B1 exon-DAG 上用独立 exon-GRU 的绝对兼容值前移评分；完整 DEV 继续条件 FAIL，水稻自由完整链下降、两物种背景负担上升。保留失败，不覆盖 `src/m28/prefix_value.py` 或重提 B2。学习前缀评分不是凭命名即可成立的创新。

A5 的 A3 指针诊断显示 12 个 missing reference/window pairs 中 8 个在所有 16 queries 的某个合法 entry 原指针局部排名均超过 4。A8 只改变完成链读出，没有新增保留命中，停止扩展。A9 的具体不同点是以冻结 A3 的 query 初始状态、prospective recurrent state、已完成 CDS 的逐碱基平均表示及六项几何，给**全部**合法 ATG/AG entry 在 top-4 之前评分；原 pointer 排名不预筛站点。终点同样读取全部原 ORF 合法事件。beam4、每 query 完成4、exon cap128，GT/GC–AG、intron>=5、实际 codon carry 与完整 start/stop 语法均保留。

新增 `query_prefix_value.py` 头仅 67,745 参数、H128、seed0、FP32；原 A3 `step_004608.pt` 全部冻结。事件为 START/NEXT_START/DONOR/STOP；非加性 logit 不加旧累计 path score，不把均值表示称作消除了所有长度偏差。旧 pointer/EOS 分仅作 receipt，不能给非单调新价值施加提前停止 bound。完成值非正也保存；本轮没有最终注释或阈值选择。

## 本轮监督定义与隔离

只用既有 8 个 A5 TRAIN pressure windows，拟南芥/水稻各4、27 primary reference/window pairs；它们不是独立准确率队列，也不替换将来拟合的原 1,536 draws / 1,491 unique windows。

所有新随机头的 DNA-only 自由预测先完成并关闭文件，再接入当前 TRAIN catalog。监督采样先保存/重放既有 label-free greedy/beam4 前缀，加所有初始 ATG；随后按冻结 A3 原系数与 known-query-mean Hungarian 给每个 primary 分配 query。真实 prefix/entry 与 teacher-history 下 raw-top4 entry alternatives、全部合法 end 形成固定 TRAIN-only 样本；合法生成不可达的真实事件显式记为 grammar-unavailable injection，不删正例、不冒充 free。预检要求正事件总数恰为 primary exon 数的两倍。

query-slot compatibility=1 只指向该 query 分配的完整 primary 的精确前缀/下一个 start；DONOR 必须是 proper prefix，STOP 必须是完整链。另一已知 primary 对当前 query 为 slot-negative，不是生物学假基因。未匹配 unknown query 全部 -1；其他 isoform 的精确兼容前缀、unknown span、noncallable span 保持 -1；assigned 正事件优先保留。known-empty TRAIN 窗口没有虚构正例。

损失是每个实际存在 known labels 的事件类型内 mean BCE，再对这些类型等权平均；不是正/负类均衡。不把 logit 解释为自然正确概率。空池/全 unknown 为可微零；统计各类型正负/unknown 数量，若实际监督极稀疏或保护策略使监督为空，先据此决定是否适合正式拟合。新增头反向必须有限、known 样本时非零；冻结 generator 无梯度，预检前后新头权重逐参数相等。没有 optimizer、checkpoint 写入、新 GLM forward、DEV 或 test/Setaria。

## 运行与后续

入口 `sbatch/M28-A9-QUERY-PREFIX-PREFLIGHT.sbatch`，新输出 `outputs/M28-A9-QUERY-PREFIX-PREFLIGHT`，1 RTX3090/2CPU/16GiB/30分钟，环境 generanno。8窗全站点随机头解码加8窗构造监督/反向；7项合成测试先验证第五起点可进入、非单调价值不被错误截断、carry/软掩码/padding、coding mean、未知保护、loss 及梯度边界。新产物估计<200MiB，不覆盖旧结果。调度按提交前实时队列；失败保留证据，明确实现错误允许一次最小修复的新目录重试，不自动延时或变更科学参数。

准备成功只证明可执行与监督口径，**不是**模型学习成功、free recovery、F1、背景控制或论文资格。B2 的 off-policy 风险尚未解决，不能单靠本轮结构不同断言已解决。正式拟合及同头 early/late 对照尚未在本文件固定或启动；须在更新前固定实际采样、资源、停止条件与继续条件，沿用既有授权自主推进。

## 实际准备结果

`sbatch --parsable sbatch/M28-A9-QUERY-PREFIX-PREFLIGHT.sbatch` 返回13359370；root/batch/extern均COMPLETED/0:0，gpu035单RTX3090、2CPU、16GiB，实际55秒（0.015278 allocated GPUh），batch MaxRSS1,795,712KiB。7项测试通过、0.903秒；8窗新头随机解码18.4017秒，全部返回链DNA合法，监督/反向加读取总程序34.7089秒。新增文件合计13,850,883 bytes，GPU peak223,748,608 bytes。

27 primary pairs/216 exon均保留，432/35,557/33,700正/负/unknown事件，299正事件仅TRAIN teacher注入；8窗新头梯度范数3.26至11.88，frozen generator无梯度，新头逐参数不变，optimizer=0/checkpoint=0。逐身份与前缀回读并独立重算loss最大差2.904e-8；不是训练成功或参考恢复结果。后续单seed拟合与TRAIN机制筛选另在[a9_query_prefix_fit_protocol.md](a9_query_prefix_fit_protocol.md)于更新前固定，准备源码不再覆盖旧B2。
