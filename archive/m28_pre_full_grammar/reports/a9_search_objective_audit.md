# A9 搜索监督目标与相关工作核查

2026-10-03。结论：当前 A9 是 query 条件化的非加性前缀相容性分类器，用于改变近似搜索的排名；不是已实现的 cost-to-go 学习、DAgger、LOLS 或 SaIL。冻结生成器不等于冻结搜索状态分布，列举全部初始 ATG 也不等于覆盖新策略的深层历史。这个区别影响研究解释与下一版设计，不改变正在运行的 A9、FULL、A6、A7，亦不增加本轮通过门槛。

## 新增原论文证据

使用 Exa 作三个定向检索，每个返回额度5，共15个结果；去重后核对以下三篇原论文的方法，不是15篇全文精读，也不是穷尽性新颖性检索。只向检索服务提供公开方法名称和问题，没有上传未公开代码、数据或结果。DAgger 在这里指 Dataset Aggregation，与既有相关工作记录中2001年的基因识别 DAGGER 不是同一算法。

| 原论文及核查位置 | 已有机制 | 与 A9 的实际区别（本项目比较推断） |
|---|---|---|
| Ross、Gordon、Bagnell，AISTATS 2011，[DAgger](https://proceedings.mlr.press/v15/ross11a.html)，第3节、Algorithm 3.1；[原文](https://proceedings.mlr.press/v15/ross11a/ross11a.pdf) | 迭代执行当前策略（可混合专家），为其访问状态取得专家动作，聚合数据后更新策略。 | A9 的自由轨迹来自固定 A3 的 greedy/beam，后续拟合未按更新中的新头重采轨迹。它不复现该算法，不能借用 learner-induced distribution 下的保证。 |
| Chang 等，ICML 2015，[Learning to Search Better than Your Teacher](https://proceedings.mlr.press/v37/changb15.html)，第2节、Algorithm 1；[原文](https://proceedings.mlr.press/v37/changb15.pdf) | 用学习策略 roll-in，对动作偏离进行完整 roll-out，以终态结构损失生成 cost-sensitive 监督；roll-out 可混合参考和学习策略。 | A9 在固定 teacher 历史上列举 entry/end 事件，再给精确前缀相容性标签，没有对每个替代动作完整续接并计算终态损失。它不是 LOLS 的实现。 |
| Bhardwaj、Choudhury、Scherer，CoRL 2017，[Learning Heuristic Search via Imitation](https://proceedings.mlr.press/v78/bhardwaj17a.html)，第2.2、3.2、3.3节、Algorithm 2；[作者原文](https://mohakbhardwaj.github.io/SaIL/pdfs/bhardwaj17a.pdf) | SaIL 在搜索状态上用 clairvoyant oracle 的 cost-to-go 监督 frontier 选择，并聚合当前混合策略的搜索数据；其主要目标是减少到达可行解的搜索展开。 | A9 不学习该搜索代价，beam 截断也不是该 best-first frontier 过程。原方法的搜索效率目标不能直接替代本项目完整 CDS 链准确率或背景负担目标。 |

SaIL 的会议年份由 PMLR 正式条目核实为2017；PDF抽取首次仅返回摘要及聚合元数据，因此方法事实改由作者托管原文核对，不使用抽取器的 citation count、日期或“preprint”标签作发表事实。既有 early-update、LaSO、BSO 和基因结构图证据见[相关工作记录](../../docs/M28_learned_search_related_work_20261001.md)，这里不重写 B2 的历史合同或结果。

## 对照当前实际源码

本轮从 Baobab 读取下列实现，而非仅依据方法命名：

- `run_query_prefix_value.py::fit()`：先关闭全部1491窗固定 A3 自由轨迹，再读 TRAIN；每个 draw 重新计算冻结表示和固定 matching，但没有用更新中的新增头进行搜索状态采集。
- `query_prefix_value_training.py::collect_free/collect_teacher()`：所有初始 ATG、原自由链事件、真实 teacher 前缀及 raw-top4 entry／全部合法 end 替代事件。替代事件的表示不是完整反事实续接轨迹。
- `query_prefix_value_training.py::compatible/prefix_labels/prefix_loss()`：1 表示与分配 primary 的精确前缀／下一 entry 相容；0 为参考与保护规则下的 slot-negative；unknown 为-1。损失是各存在 known 样本的事件类型 mean BCE 再等权平均，不是完整输出损失或 rollout 回报。
- `query_prefix_value_decode.py::decode_prefix_values/complete_query_values()` 与 `query_prefix_pair.py`：非单调前缀分数排名，完成链一次 union STOP 读出，再在各臂按 geometry 对 query 值取 max；没有可采纳 heuristic 或全局最优保证。

由此可作的推断是：新头可能访问固定控制池与 teacher 历史没有出现过的组合，哪怕 recurrent 状态的生成规则及参数完全冻结。全初始 ATG 只覆盖第一步入口，不证明深层状态已覆盖；一个前缀得到正相容标签，也不证明现有 beam4／cap128 下实际可完成或最终会被 owner／selector 保留。这里指出的是当前设计的证据缺口，不是已观察到的 A9 失败原因。

## 对下一步决策的影响

继续使用[原拟合和 TRAIN probe 协议](a9_query_prefix_fit_protocol.md)，最后 step004608、seed0、采样、损失、资源和原最低机制条件均不变。该条件只决定是否值得进入完整 DEV 准备；不把本次论文对照追加成否决条件，也不改为新的训练算法。

若原 TRAIN 条件通过，未来完整 raw-FASTA 对照仍需在前向前另行固定；必须看完整自由池、owner、最终 exact-chain F1 和背景负担。即使有收益，也只能先称这套系统下的候选生成／搜索配置收益：early 去掉原 raw bound 并改变排名，late 另含 greedy，候选曝光和工作量不同，不能归因为唯一评分位置变量或等预算优势。

若原条件失败，关闭当前头的扩展，不以改类权重、追加步骤、换 seed／beam／阈值挽救。下一版若选择策略分布采样或 rollout 结构损失，那是独立的新科学假设，须先定义 TRAIN oracle、unknown／非目标 isoform 处理、采样与资源，再运行；本记录没有选择或启动它。不能把旧失败泛化成学习搜索不可能，也不能凭文献名称断言新版本必然有效。

对论文贡献的当前可检验表述仍是：在 DNA-only 输入、真实跨剪接 codon carry 和明确参考保护政策下，改善完整 primary CDS 链恢复，并把收益保留到染色体最终输出。组件先例存在；任务特定组合是否新颖、能否胜过 ANNEVO/Tiberius/Helixer、是否具有独立迁移收益，都需要实际证据，未由本次文献核查证明。

本轮只有文献与解释记录。没有新 optimizer 更新、GLM forward、DEV/test/Setaria 读取、作业提交／取消、goal 状态变更、Git 提交或公共发布；等待继续由既有固定监控负责。
