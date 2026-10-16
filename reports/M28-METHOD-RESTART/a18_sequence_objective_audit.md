# A18 训练事件均值与完整链分数诊断

2026-10-05，诊断计算前固定。当前A15／A17恢复链保持原协议运行。本诊断只用已完成A5／A16的原八个TRAIN压力窗口、27条参考／窗口配对及432个保存的条件探针，判断事件均值与完整事件和是否实际改变query分配；不追加训练、模型前向、自然解码或DEV，不改变已运行模型。

## 代码事实与待检验问题

`teacher_chain_components`把一条E-exon链的2E个边界CE平均为B，把E个EOS BCE平均为T；A3按已知query平均存在性、按正链平均B和T，Hungarian成本对应其完整窗目标。原beam累加2E个几何边界log-prob和E个EOS log-prob。A11仅改完成池排名并失败，不等于训练过完整链事件和；A16已量化合法域与teacher／decoder下界差，不能据形式差异断言失败原因。

令N为窗内正链数，z为query存在性logit，D为原存在性配对增量：unmatched0时D=−z/Q，unmatched−1时D=softplus(−z)/N。反事实成本固定为：MEAN=D+(B+T)/N；EVENT_RATE=D+(2B+T)/N；SEQUENCE_SUM=D+E(2B+T)/N。EVENT_RATE只分离边界／EOS事件比例，SEQUENCE_SUM再引入链长权重，不扫系数、指数或阈值。这些成本不是模型质量分或新的注释操作点。

三个成本中的D保持不变，因此EVENT_RATE也改变边界项相对存在性项的权重，SEQUENCE_SUM也改变整链结构项相对存在性项的权重；不能把分配变化单独归因于长链竞争。此处只判断固定成本定义的后果，不做因果分解。

用A16保存的teacher entry log-prob恢复B，结束边界和EOS复用原事件；验证每个原raw path分数满足`raw=−E(2B+T)+decoder/teacher entry下界差之和`。A5保存的query logits补齐存在性项，真实R4 TRAIN catalog和原保存DNA重建unknown政策／正链身份，不假定所有窗口都可监督未匹配query为负。保留所有27身份、16queries、原自然命中／局部必要条件和原draw曝光计数；部分压力窗未被原1536 draw抽中，不能称全部为模型训练样本。

## 执行与边界

新增独立CPU脚本`audit_sequence_objective.py`与SBATCH，generanno／private／0GPU／2CPU／8GiB／10分钟，新输出`outputs/M28-A18-SEQUENCE-OBJECTIVE-AUDIT`。三项解析检查检测事件计数倍数、不同链长影响一对一分配及unknown存在性增量；随后重算八个成本矩阵、三组分配和全部原raw分数恒等式。输出windows／references JSONL与summary／STATUS／JOBID，约数百KiB，不读取大特征文件或checkpoint权重，不覆盖任何旧产物。实际job和结果在取得证据后补充。

这是保存的单query FP32探针用float64重组的反事实分配，不是原batch16 CUDA Hungarian逐bit复现，近似相等成本的分配需谨慎解释；不称与原模型保存的实际训练assignment相同。强制参考历史上的局部width4通过是必要条件，不是新的自然命中、F1、背景改善或泛化。

若SEQUENCE_SUM相对EVENT_RATE改变分配，说明链长权重在本固定探针中确实影响query资源分配；仅EVENT_RATE改变则组件事件比例已有影响；都不改变时，不据形式差异自动训练新目标。无论哪种情况，都不能证明相应训练会改善模型。A15／A17 gate及全部历史失败保持不变，新目标若开展须单独固定训练和自然读出协议，不借本诊断自动延长当前训练、放行DEV或公开发布。

## 首次执行与局部恢复

首次CPU作业13388763的root／batch为FAILED、ExitCode=1:0、Elapsed=00:01:17；三个解析检查通过，随后因错误要求每个窗口必须含正链而停止，未写出诊断summary。原八窗包含六个正链窗和两个`zero_positive_unprotected`窗，27条参考身份无缺失；失败不是模型的科学负结果。

R1只修复零正链窗口的空成本矩阵／空匹配，并用保存的greedy query logits支持这两个无参考profile窗口。新增一次零正链回归检查，正链成本公式、三臂、参考分母、原数据隔离及其他校验保持不变。全部八窗仍写入结果，背景窗的assigned-cost=0只表示没有正链匹配项，不是完整背景损失或背景性能。保留首次失败目录和日志，恢复输出使用`outputs/M28-A18-SEQUENCE-OBJECTIVE-AUDIT-RECOVERY-R1`；不重跑任何模型前向。

## R1实际结果

恢复CPU作业13388894通过`sbatch --parsable sbatch/M28-A18-SEQUENCE-OBJECTIVE-AUDIT.sbatch`提交，private／0GPU／2CPU／8GiB／10分钟。Slurm root、batch、extern均COMPLETED、ExitCode=0:0、Elapsed=00:00:25；目录JOBID=13388894、STATUS=COMPLETED，stderr为空。三个解析检查及零正链回归通过；windows.jsonl为8行、references.jsonl为27行，432个原条件探针身份和raw分数分解全部通过，最大分解绝对误差7.11e-15。这是保存值的算术恒等式，不是native batch16前向等价性验证。

EVENT_RATE相对MEAN没有改变任何参考的query分配；SEQUENCE_SUM相对两者改变4／27，均来自拟南芥同一窗`NC_003070.9|4620288|-`。拟南芥4／17、水稻0／10；one-exon 0／9、2–4-exon 3／9、≥5-exon 1／9。所有27条参考配对所属的六个正链压力窗均未被原A3的1536 draw抽中，不可称为训练曝光样本，也不是独立测试或泛化证据。

|参考ID|exons|MEAN／EVENT_RATE query|SEQUENCE_SUM query|三臂局部必要条件|
|---|---:|---:|---:|---|
|rna-NM_101222.4|2|11|13|均通过|
|rna-NM_101224.2|2|13|7|均通过|
|rna-NM_101225.2|2|6|11|均未通过|
|rna-NM_101226.5|18|7|6|均通过|

对保存的11×16成本矩阵，用query子集动态规划重算最小一对一总成本，三臂最优值与保存的SciPy assignment总成本相符（绝对差<1e-12）。对四条变化参考逐一禁止当前最优query后重新优化，MEAN的最小成本增量4.27e-4，EVENT_RATE为8.54e-4，SEQUENCE_SUM为1.74e-3；因此不是保存float64矩阵中的精确并列取舍。SEQUENCE_SUM采用旧MEAN分配的成本7.236933154263542，高于其最优成本7.196604120302874，差0.040329033960668。不同臂的总成本尺度不能直接比较，此检查也不量化native batched前向误差。

动态规划复核只读R1的windows.jsonl；按参考原顺序迭代，状态为已占用query集合S，初值DP₀(∅)=0，转移为`DPᵢ₊₁(S∪{q})=min(DPᵢ(S)+C[i,q])`，q不属于S；禁用检查只去掉指定参考／query转移，最终取所有大小为11的集合最小值。四条变化分别是原矩阵零基行5、7、8、9，没有改变成本或重新搜索候选链。

三臂分配后的局部width4必要条件通过数都为15／27：拟南芥9／17、水稻6／10；one-exon 4／9、2–4-exon 4／9、≥5-exon 7／9。四条变化参考也没有任何一条从局部未通过变为通过或反向变化。原A5自然预测完全复用，本次未增加自然命中、候选链、F1或背景性能读出；optimizer_steps=0、模型／backbone forward=false、checkpoint_written=false、DEV=false、test_setaria=false。

结论：完整链事件和确实能改变这个固定探针中的资源分配，但当前仅一个拟南芥窗发生变化，且没有局部必要条件增益；不能把形式上的目标差异定为既有模型失败原因，也不能据此宣称长链恢复或新模型性能提升。本诊断不选出新的训练赢家，不以这4次分配变化自动开展序列和首训。A15／A17 R1继续按原冻结fit／自然TRAIN读出／gate推进，待实际结果决定后续；本次不改其参数、阈值或数据范围。
