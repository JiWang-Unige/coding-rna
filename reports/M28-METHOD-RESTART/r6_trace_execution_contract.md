# M28 R6 搜索追踪执行合同

2026-09-30，新前向之前固定。此合同补充[r6_candidate_trace_plan.md](r6_candidate_trace_plan.md)，只决定是否实施剪枝前学习兼容性的 B2；不作为新的模型性能或论文主张。

## 方法选择及解释修正

Pro 第22轮实际读取提交 c1ecee5a734a75da10ff928f7f6ed87cfe84813c 的 R5 与 D1 两份 Markdown、两份 JSON，以及 candidates.py、core.py、training.py。未访问 Baobab 缓存、未独立重跑。该回复是设计咨询，不是执行授权或独立验证。

C0 0.935735、B1 0.561102 是最终预测链 CDS 的无链向并集 F1，不是原始局部分类头的质量评估。不能由此归因多任务梯度损害局部表示。R5 两物种一致优势失败，以及仅重排当前最终池不足以达到缓存 Helixer 的能力界，仍成立；不否定重新生成候选的新模型。

选择一个下一方法假设：让专门训练的链前缀兼容性在 beam 和末次截断前生效。旧完整链 BCE 训练的 GRU 不能直接视作前缀评分器。首个方法对照须共享并冻结同一局部表示和端点连接头，两臂同起点、相同新增训练量，预算、归属和输出规则一致；不混入局部表示保护收益。不改 K、seed、R5 步数或阈值。

## 固定输入与观察范围

[r6_trace_sample.json](r6_trace_sample.json) 已在读取 checkpoint 前写入：四个物种链向 scope、每类按 reference_id 字典序取一条，共20 anchors、23个去重窗口。没有缺类，不补选或扩展。每个 nonowner-positive anchor 的 owner 窗与原阳性见证窗均包含在既有窗口中；终止检查核对原保存 gain 为正的非 owner 见证存在。

checkpoint 固定 outputs/M28-PILOT-R5/B1/step_004608.pt、4608步；GENERanno 特征 revision b0483c23b6b63787b61a6d3a204a9b517d6ba345。只读取该样本的既有 DEV 特征及冻结 R5 JSONL；不读 TRAIN 特征、test 或 Setaria，不载入 GLM，不训练、重缓存、生成新 GFF 或全 DEV 排名。

## 输出等价与首次淘汰证据

原 generate 与独立旁路副本接受相同 DNA 和神经输出。参考只供旁路观察，不能改变排序、得分、候选或预算。原 src/m28/candidates.py 和模型实现不修改。

每窗两个实现的候选身份与顺序、links、counts、budget、reference_used 完全相同，proposal 最大绝对差为0；原生成器与 R5 保存记录的前述离散字段完全相同，proposal 最大绝对差不超过1e-5。任何不一致立即失败，不放宽容差，不以失败回放作归因。

观察原始与保留端点、ORF/carry、必要 links、参考 carry 下条件化 exon 排名、实际 entry/donor 精确前缀、完整链末次截断。真实前缀首次在 beam 或 final 被淘汰时，记录同事件同 carry 存活竞争者的完整几何身份、目标名次、cutoff 分差及 start/donor/stop/acceptor、exon 区域均值 log-odds、link 分项。分项须重构真实加性得分；否则追踪失败。只保存所选参考相关事件，不保存全部 beam。

条件化 exon 可得性不是实际执行轨迹；前缀丢失后的条件化观察明确标记 actual_execution_claim=false。单 exon 的 links-ready 是空条件，仍须核对完整 start/stop/exon 支持。首次实际失败沿参考路径顺序定位；entry/complete 本应到达却缺失属于观察错误疑点，不包装成 beam 机制。

## 一次性决策与预算

若两物种样例中分别出现真实可续接前缀，在加性 entry/donor beam 或最终截断被淘汰，则准备 B2 的冻结首训合同。为避免把更晚的另一个阻断误称可恢复，机制见证还要求全部参考端点、必要连接和条件化 exon 边保留，且整链 ORF 语法合法。这是存在性依据，不估计总体归因率。

若未获得两物种见证，不扩大样本或追踪预算寻找支持案例，不自动进入该 B2 版本。若只见端点/exon/语法阻断，应改变候选获取或表示监督；若两物种支持不一致，明确记录并重新选方法，而非改门槛。回放、对照或轨迹有误则诊断未完成，不强作方法判断。

单个 Slurm 作业：RTX3090/24GiB、2CPU、16GiB RAM、15分钟，新增 allocated GPUh ≤0.250000。这是新诊断成本，不挪 R5 余额、不改旧预算。环境 generanno，复用既有初始化。输出 outputs/M28-R6-CANDIDATE-TRACE 必须不存在；日志使用 M28R6TRACE 与 job ID 唯一路径，复用既有日志父目录。compact 全产物 ≤32MiB。

不安排第二个 GPU 作业；超时、OOM、非有限、回放不一致、阳性对照或非 owner 见证丢失、产物超限均停止并保留失败。提交前只做会改变执行决定的语法、合成轨迹对照、实际资源和路径检查。完成后以 Slurm 状态、summary、每窗等价和事件证据共同判定，不能仅凭 COMPLETED 或机制计数宣称研究成功。
