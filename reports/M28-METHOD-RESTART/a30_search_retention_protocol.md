# A30 当前搜索首淘汰竞争损失：首次固定 TRAIN 协议

2026-10-06。独立检验搜索竞争目标，不叠加A23相对边界或A25首起点结构。原四监督保留；额外目标只比较当前自然搜索中首次淘汰的正确路径与有安全负例资格的cutoff。当前尚未拟合，不预称性能收益或论文创新。

## 唯一方法变化

使用fresh seed0原A3 CollaborativeGeneQueryModel，649457参数，无新增head／参数、无监督checkpoint继承。每步先按原四loss及Hungarian成本完成teacher分配；对全部已分配正query，在当前权重下观察原beam。参考链只判定首次淘汰，不进入或保留任何搜索候选。

对raw grammar支持且未自然返回的目标，回放首次淘汰的目标与cutoff原累计log分数。仅完整R4注释判为known_nonmatch 0的cutoff计算relu(1＋cutoff_score−target_score)；另一primary、保护或未决cutoff为零，无负梯度。已返回及raw不支持的参考额外项也为零。后一类保留全部原teacher监督，记录raw支持差异，不删除训练目标。按全部已分配正query平均，包含所有零项；无正query为连接计算图的零。

首次固定margin=1个累计log分数单位，loss系数=1，不搜索系数或margin。总损失为原importance-weighted四loss，加相同importance weight乘额外平均项。用model.query_window取得未detach的共享特征；保留原teacher分配规则。局部start／next_start top4及四种全局首淘汰均纳入同一既定规则，不使用A28旧权重路径离线训练。

geometric log-softmax、GRU／EOS、ORF carry、raw最短内含子5bp、beam4、每query completed4、exon cap128、score-bound停止及自然greedy＋beam联合池均不变。未加长度归一化、score校准、候选注入、解码阈值或完整caller。

## 数据和曝光

沿用原A3的TRAIN loader、两物种六条TRAIN染色体、原1536 draw／1491唯一窗口、GEN固定缓存版本、R4 eligible primary及保护规则。单seed0、FP32／batch1、AdamW lr3e-4、weight_decay .01、betas(.9,.999)、eps1e-8、clip1，原draw顺序重复3次，共4608更新。末step004608唯一主读出，每256步只作工程保存，无早停、重采样、选最好checkpoint或失败后追加曝光。无新GLM前向、DEV、sealed或Setaria访问。

先运行原metadata挑选的八项零更新真实TRAIN预检；不按loss或margin激活情况换案例。原四loss及分配与fresh A3逐case相同、当前beam与原实现的自然完成几何／顺序／累计分数一致、新额外项数值及梯度有限，且固定案例中两物种均实际接通共享和query梯度，才进入正式拟合。所有参数不变、optimizer0。若不通过，区分实施问题与科学支持不足；不伪造非零信号或缩减保护规则。

## 本候选首次固定五项 gate

末模型在全部1491 TRAIN窗完成reference-free原greedy＋beam4自然候选池后，再按原R4标注，CPU复算逐记录及计数，基准为已闭合原A3自然池。

1. 拟南芥正确自然候选数严格大于A3的1751。
2. 水稻正确自然候选数严格大于A3的629。
3. 合并≥5外显子正确自然候选数严格大于A3的788。
4. 拟南芥known_nonmatch自然候选数不超过A3的6305。
5. 水稻known_nonmatch自然候选数不超过A3的6261。

这是A30首次拟合前的新协议，不修改A23／A25原三gate或其已完成判定。错误候选数量约束针对A23暴露的实际候选膨胀，仍不是最终选出gene的FP或精度。五项全部通过才保留；任一失败关闭A30，不用新seed、更新次数、系数、阈值或另一基准救回。同一方法的实施故障仅作边界明确的工程修复，不改变科学配置或覆盖冻结结果。

PASS只支持另行固定完整caller／链质量评价，不自动开放DEV。TRAIN探索不能证明独立泛化、因果机制或优于ANNEVO／Tiberius／Helixer；这仍须后续相应冻结实验。

## 执行和等待

八项预检用唯一outputs/M28-A30-SEARCH-RETENTION-PREFLIGHT，private／1RTX3090 24GiB／2CPU／16GiB／10分钟，沿用实际generanno环境。正式fit资源根据该预检实际计时确定，不能为凑资源上限减少原曝光。长等待只复用原固定监控对话和heartbeat，监控只看调度、完成标记和本候选字面结构字段，不做科学解释或重提作业。
