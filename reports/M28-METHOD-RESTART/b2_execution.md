# M28 B2 执行记录

本轮执行[b2_first_fit_contract.md](b2_first_fit_contract.md)中已固定的单seed评分位置机制对照。目标是完整拟合、全DEV配对解码与最终评价，不以提交成功、训练loss下降或合法GFF结构判定方法成功。

## 提交前实现与检查

2026-10-01提交前记录：新增fit_b2.py、infer_b2.py、evaluate_b2.py三阶段驱动，以及prefix_inference.py的单次union读出与R5重放核对。复用既有标签、候选Budget、owner、selector、GFF导出和统一评价函数。仅新增头训练；原B1与GLM冻结资产保持不变。

Python AST与三份sbatch的bash语法均通过。两项新增配对测试安排在拟合作业prologue内、任何optimizer更新之前：共同链在不同候选顺序中得到相同gain；控制接受JSON坐标表示但拒绝真实离散变动、参考参与、非有限或超过1e-5的proposal误差。测试时间计入该3小时GPU分配，不重复未改动的6项准备测试。

训练原1,536 draws的顺序及importance weight不变，固定加性控制池/首次截断观察按1491个唯一window ID在拟合内复用。每步读取固定特征、重新计算冻结表示，构造独立的新头计算图；不缓存新头logit或hidden跨optimizer更新。注入完整正链只出现在损失池，实际free及首次剪枝观察单独保存。

全DEV两种解码共用一个最终新增checkpoint。late先重放原B1实际固定池，再与early的候选组成单次最终读出union；每窗保存两池、proposal与新gain，记录控制误差和共同链gain相等。推理不读参考，评价在全部8690窗完成后才接入冻结7728主参考及6450历史副表。

B2继承B1已完成的4608次基座更新，另增加4608次新头更新；新增更新数不是整套模型的训练总数。early与late的曝光和权重匹配，但B2与历史C0/B1没有因复用旧资产而自动变成训练曝光匹配的完整方案因果对照。基座历史训练与特征成本另按R5/R4来源公开保留。

## 调度选择与资源

提交前Slurm记录显示gpu034为DRAIN，gpu035为MIXED，8张RTX3090中已分配7张；CPU已分配36/126，RAM已分配212992/503000MiB。选择private-teodoro-gpu，因为存在合格单卡及CPU/RAM空间，无GPU相关维护预约限制本轮墙时。共享分区不在本次选择中替代GPU型号或数值设置。

全局BeeGFS可用约164TiB，新产物上限仍为2GiB。下列三目录在提交前均不存在，既有159GiB特征与原权重仅只读引用。首次节点查询因scontrol每次只接受一个node参数而返回参数错误；分别查询两节点后获得上述事实，无作业因此启动或重提。

| 阶段 | sbatch | 新目录 | 分配上限 |
|---|---|---|---|
| 拟合 | M28-B2-FIRST-FIT.sbatch | M28-B2-FIRST-FIT | 1 RTX3090 2CPU 16GiB 3小时 |
| 配对推理 | M28-B2-PAIRED-INFER.sbatch | M28-B2-PAIRED-INFER | 1 RTX3090 2CPU 16GiB 8小时 |
| 评价 | M28-B2-EVALUATE.sbatch | M28-B2-EVALUATION | 0GPU 2CPU 32GiB 2小时 |

推理afterok拟合，评价afterok推理；序列资源上限为11GPUh与4独立CPUh，不改变旧实验上限。失败、控制不匹配、超时或产物超限停止并保留证据，不自动重跑、续训、延时、改变阈值/数据或取消别人的作业。Slurm日志沿用outputs/M28-PILOT-R5/logs下独立的新jobname/jobID文件，结果目录与旧R5完全分离。

## 提交与终态记录

此节在实际提交后补入命令与Slurm ID，并在整条流程终态后记录分配成本、产物及原条件是否满足。提交前没有B2正式拟合checkpoint或DEV性能结果。以上调度状态是提交前快照，不作以后重提的事实依据。


## 2026-10-01 实际提交与启动取证

实现于公开提交63730f07d4c19a0850e9e1813a7c89747bd34403固定，随后按以下命令各提交一次，无重提：

```bash
sbatch --parsable sbatch/M28-B2-FIRST-FIT.sbatch
sbatch --parsable --dependency=afterok:13312819 sbatch/M28-B2-PAIRED-INFER.sbatch
sbatch --parsable --dependency=afterok:13312821 sbatch/M28-B2-EVALUATE.sbatch
```

| 阶段 | Slurm ID | 启动取证状态 |
|---|---:|---|
| 拟合 | 13312819 | RUNNING，gpu035，1 RTX3090、2CPU、16GiB |
| 配对推理 | 13312821 | PENDING Dependency |
| 最终评价 | 13312826 | PENDING Dependency |

拟合启动取证时已分配610秒、STATUS为RUNNING；两项配对测试均通过，用时1.692秒。日志已进入第576次optimizer更新、第一遍，565个唯一固定控制窗口已生成。本快照不是最终更新数、总分配成本或模型收敛证据；未完成全部4608步和8690窗时仍不排名。

原独立监控已更新为这三个作业与新完成产物，保留gpt-5.6-luna、low、每15分钟、本机运行；正常运行/依赖等待不向主对话报告，整链完成或发生失败/超时/依赖无法满足/产物缺失时一次交接，成功后暂停该轮监控。仅记录状态和必要错误，不计算性能，不改变参数、预算、数据范围或goal状态。截至本轮配置时，更新启用已确认；新cohort的独立调度检查与完整回传尚未观察，不能用旧R5/Pro的成功回传替代这两项验证。

本机定时监控依赖机器开机、应用运行和连接可用，不保证休眠或断线时即时回传；Slurm作业本身不依赖本机SSH会话持续。最终结果须回到主对话核对工程完成、指标可信、机制条件及后续科学方向，不能由监控作出论文结论。
