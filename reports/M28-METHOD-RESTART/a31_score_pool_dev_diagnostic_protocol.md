# A31 固定评分与候选池 DEV 诊断

2026-10-06，真实诊断执行前固定。完整A31 DEV已经闭合；本轮只解释已保存的A3部署差异，不新增训练、候选搜索或模型选择，也不重判任何原门槛。沿用同一个A3末004608、全部8690窗、四scope、7728条primary参考、null0及原owner／重复取最大／非重叠DP。

三格固定为原A3 greedy＋原query分数、同一greedy候选记录＋A31保存的局部分数、原A31 A3 greedy+beam联合池＋其局部分数。逐窗核对身份、网格，以及旧greedy去重几何恰等于新联合池中标记greedy的几何集合；按坐标映射局部分数，保留旧greedy的重复query记录，由原装配器执行owner和重复取最大。若完整几何不一致，停止匹配对照，不猜补分数、删掉不匹配窗或将其称为工程错误重跑。旧A3来自缓存推理，新A31来自raw-FASTA；逐窗几何匹配不证明编码数值完全相同，因此评分对照也不能宣称纯同一hidden上的因果效应。

source roots为旧A3推理13357568／评价13357569、新A31 A3推理13399434／评价13399454、已保存共同参考诊断13360116。实际sacct须全部COMPLETED／0:0，STATUS／JOBID与相关协议、末checkpoint须一致。先完成全部候选几何及三格参考无关装配，再读取既有7728行references.tsv；只复用已核验参考，不重开DNA、原GFF、feature cache或R4 catalog。两个原格必须精确复现原完整链TP／预测数／参考数、全部stage计数、原A3 query计数及A31装配计数，否则停止解释。

报告每物种和合并完整链P／R／F1、free→owner→selected参考覆盖、至少五CDS片段正确选择。第一项差固定greedy池改变评分，第二项差固定局部分数扩池；只报告沿此顺序的条件差，不称独立因果份额，不估计交互。缺少union＋query格，理由及评分语义边界沿用TRAIN诊断；不复制query存在logit给beam备选、不直接使用非正累计path log概率，也不新增阈值／归一化／校准。新增中间格重组选择只作诊断，不输出生产注释GFF，不用于选择A3／A23或独立面板。

先三项合成检查：不同候选顺序及重复greedy仍正确映射且复用原重复取最大规则；未知／非参考高分候选仍参与DP而不通过标签删池；缺失／额外greedy几何或身份不一致拒绝匹配。随后private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟、generanno，唯一输出outputs/M28-A31-SCORE-POOL-DEV-DIAGNOSTIC，命令sbatch --parsable sbatch/M28-A31-SCORE-POOL-DEV-DIAGNOSTIC.sbatch。提交前按Baobab skill做实时调度及输出冲突检查。主对话直接消费短CPU作业，不启用已消费A31的心跳或新建监控对话。

预期产物：summary.json、8690行greedy_local_scores.jsonl、四scope的diagnostic_selections.jsonl、JOBID、STATUS及Slurm日志。失败只允许一次不改变科学条件的最小工程修复，原失败产物保留；几何不匹配或原格不复现属于证据边界，不能靠修参数救结果。无GEN／模型前向、optimizer、新参数、sealed／Setaria访问、阈值扫描、公开发布或Git提交。

## 实际执行闭合

2026-10-06 09:36 UTC。真实作业13402574 root／batch／extern COMPLETED／0:0，实际09:31:19至09:32:09 UTC、50秒、MaxRSS898120K，三合成检查PASS／0.047秒。全部8690窗greedy几何匹配，原两个格完整stage和整数指标精确复现，保存中间分数8690行及四scope选择；原科学协议及门槛未改。[结果](a31_score_pool_dev_diagnostic_result.md)解释有序条件差，未据此选择模型、调阈值或提交新训练。该短作业已由主对话消费，未加入旧A31监控。
