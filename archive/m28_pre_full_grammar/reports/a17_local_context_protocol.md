# A17 固定局部上下文 query 生成器

2026-10-04，首次运行前固定。用户已授权自主迭代；本轮是新的单配置、单 seed 架构实验，不改 A3/A15，不放行 DEV／test／Setaria，不恢复旧框架或公开研究。

## 要检验的假设

A5 表明原 27 个 TRAIN 压力配对中，八个遗漏在全部 query 有局部必要名次失败，另有完成路径竞争。A16 核验后 start／next-start 合法概率约 97.8%／99.6%，不支持仅靠入口 mask 解释这些失败。原 query 从全窗 attention 汇总初始上下文；没有证据已证明 query 塌缩是根因。

新假设是按位置区分初始上下文，可能使不同 query 学到不同局部基因线索，改善自然完整链覆盖。这不是增加 beam、事后改整链分数、扩大训练或已验证的新颖性贡献。局部区域也可能限制上下文、损害长基因或引入更多错误候选，必须按原自然输出检验。

## 唯一模型改动

对有效长度 L、Q=16，第 i 个 query 的上下文 attention 支持为 `[floor(iL/Q), max(floor(iL/Q)+1, floor((i+1)L/Q)))`，i=0..15。实际完整窗口每区约 1,536bp；不从参考坐标选 anchor，不用 padding 宽度。支持短 contig 时相邻 query 可共享一个实际碱基，不出现空 attention。区域宽度由固定 Q/L 决定，不扫描。

仍用原 query、context key/value、位置编码和 initialize 参数；之后沿用原四头 pre-norm 残差 query 通信、全窗 pointer、GRU、EOS 和 DNA／ORF 语法。局部域仅限制初始 pooling，不限制 Hungarian 正链分配、首个起点、后续 splice 或基因跨度。teacher 和所有自然解码共用新状态入口。

新模型参数仍为 649,457，完整构造顺序／seed0 权重与 A3 逐 tensor 相同，无新参数或监督 checkpoint 继承。保留原源文件和构造器，新增独立模型类；原 global pooling 的临时状态被新初始状态替换，不参与下游损失或解码。

## 固定训练与预检

沿用 A3 的六条 TRAIN 染色体、原四块 D2048 cache／revision、24,576bp 定向窗口和 R4 正／部分／alternative／unknown 标签政策。原 1,536 draws 每遍、1,491 唯一窗、物种各768 draws、原顺序和 importance weights 三遍，共4,608更新。H128/Q16、FP32/batch1、四 loss 系数1、known-query mean BCE、原精确 Hungarian 成本、AdamW lr3e-4/wd0.01/betas0.9/0.999/eps1e-8、clip1 均不改。末 step004608 是唯一主权重，前两 epoch 权重只保留记录。

三项新检查检测：完整参数初始化是否真的匹配；局部支持／短 contig／padding 是否正确；teacher/free 初始状态一致、跨区边界仍可预测且 shared/query 梯度可达。未改的旧损失、匹配、语法和导出不重跑整套测试。

随后原 metadata 方法选取的八个 TRAIN 压力 case 只做梯度预检，optimizer0，逐分支有限且至少有非零梯度，权重不变。预检通过后立即一次从头拟合，不再要求重复批准。非有限数、路由或身份错误停止相应阶段；真实工程问题允许一次有证据的最小修复，科学 FAIL 不追加 epoch、换 seed、改阈值或重提本轮。

## 自然 TRAIN 读出和继续条件

末权重在原全部1,491唯一 TRAIN 窗生成自然 greedy＋beam4 的几何去重联合池，width4/completed4/max_exons128及原 raw 路径分数不改。先关闭所有自由池输出，再打开 TRAIN 参考 catalog 做候选标注；无参考注入、自然生成中的标签或新 GEN 编码器前向。复用闭合 A6/A12 中冻结 A3 的实际同窗自然池，重新计算其原标签并逐项比较，再计算新池标签。

继续条件仅是三项生成机制筛选：两物种自然正确候选／窗口配对分别严格增加，pooled ≥5 exon 正确候选／窗口配对严格增加。known-nonmatch、unknown、全部候选数同时完整报告；它们不是已选基因 FP，也不能因自由池增加自动判断最终 F1。没有结果后再改门槛。

三项未全部通过，保留结果并关闭本架构试验；全部通过只支持为新生成器另写匹配的完整 caller／整链评分协议，不自动进入 DEV 或声称模型成功。真实论文判断仍须完整 raw-FASTA、全部参考、逐物种／最差物种 P/R/F1、背景和原生 ANNEVO/Tiberius/Helixer，独立测试与开放范围另守现有约束。

GPU 读出后一次 CPU 核验：重放4,608实际更新的 draw/weight/监督计数与原 A3 的每遍曝光，检查末 checkpoint 参数有限和结构身份；对全部新自然链直接用原 TRAIN cache DNA/targets 重算标签及三项门槛，不做新模型前向、解码或更新。1,491窗／候选配对不是独立基因或泛化样本，TRAIN机制 PASS 不替代科学 PASS。

## 运行产物和资源

独立输出 `M28-A17-LOCAL-CONTEXT-PREFLIGHT`、`FIT`、`TRAIN-READOUT`、`RESULT-CHECK`，每阶段 STATUS／JOBID／summary。拟合保存完整4,608训练记录及三份 checkpoint；读出保存全部自然 free_pool 和窗口链／标签记录。估计新增结果小于0.5GiB，不复制缓存或原始数据。

generanno，private-teodoro-gpu，预检／fit／readout各1张RTX3090、2CPU/16GiB；预检10分钟，fit45分钟（原A3实际31分钟），readout45分钟（A6原1,491窗自然联合池实际2,006秒／33.4分钟，不是完整DEV）。CPU核验0GPU、2CPU/8GiB、10分钟，可能读取约50GiB既有cache以取得真实DNA，不复制这些特征。新模型的实际速度／资源在运行后记录，不把初始预算当已测结果。

提交前按本机 Baobab Skill 核实实时调度／维护／空间和输出冲突。固定监控对话与同一 heartbeat 仅等待，完成后交给主对话验证解释；不能创建第二监控或让主对话长时间轮询。不取消任何现有作业、Git提交／推送或公开资料。
