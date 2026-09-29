# M28 R4运行诊断与隔离的chain-GRU实现验证

2026-09-30（Bangkok）。R4原作业仍运行；本文件不是正式DEV结果，不改R4合同。

## 真实训练负载

从两臂training.jsonl各取前1,536个已完成记录，按字段作算术均值：

| 每窗秒数 | C0 | B1 |
|---|---:|---:|
| load | 0.637860 | 0.628378 |
| shared forward | 0.069947 | 0.067028 |
| proposal | 0 | 0.287864 |
| supervision/backward/optimizer | 0.128487 | 0.784111 |
| 上述分项之和 | 0.836294 | 1.767382 |

源文件为远程`outputs/M28-PILOT-R4/{C0,B1}/training.jsonl`；
前1,536行属于完成的第一遍，不随后续追加改变。
load包括torch.load、CPU转float、DNA one-hot、H2D及同步，不能称纯磁盘计时。
supervision字段还包括标签构造、全部loss、反向及optimizer，不能称纯GRU计时。
B1自由候选JSON序列化/flush在上述分项之外；分项和不是完整墙时或Slurm费用。

B1自由候选mean190.984、median/p95/max192；训练候选mean192.394、
median193、p95 197、max201；参考注入mean1.410。
训练候选数与监督阶段耗时Pearson约0.212，仅是描述性关联，不用于因果归因。
接近候选数量上限不等于完整参考链覆盖高。R3十二窗热缓存计时不能保证此负载。

## Pro第20轮及采用的判断

Pro确认读取公开固定提交`fa6937cd11fc24854bb2c2268ae67f3af9e58dcc`的
r4_pipeline_execution.md、core.py、training.py、fit_pilot.py和r4_pilot_contract.md；
未访问Baobab日志或独立重跑。上述第一遍汇总由本Agent读取日志得出，
不能记成Pro独立核验。讨论保留在[原会话](https://chatgpt.com/g/g-p-6a29d586630481918525796032225f68-ji-wangke-ti/c/6ab777ce-bcd4-83eb-abff-89a8f9361948)。

代码确实每条链单独调用GRU；值得检验批处理是否降低开销，但尚未测得该子模块占比。
只合并chain GRU调用，不同时做exon去重/累计和pooling，不改候选、损失或读出。
[PyTorch GRU文档](https://docs.pytorch.org/docs/2.14/generated/torch.nn.GRU.html)
支持packed变长输入；本机实际版本的行为仍须通过数值测试验证。
工程优化本身不作为论文贡献；R4是否完成及B1是否有效仍由真实终态/完整DEV判断。

## 独立CPU验证范围（运行前固定）

新增`src/m28/packed_chain_gru.py`是隔离原型，不由任何R4 runner导入。
原`core.py`、训练/推断入口、运行参数和权重均不修改，不热替换或续训R4。
单个CPU测试作业：private-teodoro-gpu，0GPU、2CPU/2GiB、最多10分钟，
即最多1/3 allocated CPUh；另列工程用量，不改变R4 GPU/拟合/推断时限。
输出`outputs/M28-PACKED-GRU-CPU-R1`；只用合成特征和临时克隆权重，
不读TRAIN/DEV/test/Setaria数据或现有checkpoint，不生成可评价的训练模型。

测试空候选、单exon、未按长度排序的混合链、共享exon、全未知、
全参考负例、128exon长链、有效区外padding以及201候选/hidden128负载。
FP64固定atol1e-11/rtol1e-9；FP32固定atol2e-6/rtol2e-4，不在失败后放宽。
比较各输出、类均衡链损失、输入特征/各参数梯度、clip前norm、
一步AdamW参数增量及已非零的一阶/二阶矩状态；记录绝对/相对误差。
非链loss/端点连接未被修改；此CPU测试不冒称完整真实joint_losses/GPU验证。
额外构造近零gain，记录正负翻转而非用数值接近掩盖离散操作点敏感性。
一步一致不证明长期训练轨迹逐位一致。

若CPU数值不符，拒绝替换并定位具体差异；若通过，仅说明该批处理语义有初步证据。
真实GPU完整步吞吐仍未测，不能声称提速或据此恢复/扩展R4。
Pro提出的一次至多15分钟真实TRAIN GPU对照是后续提案，不是本次已执行工作，
也不授权新完整拟合。任何下一阶段先记录其范围、成本和停止条件，
保留本轮未完成状态，不将新ID作为原冻结预算的重置。

## 当前验证状态

CPU作业 **13293620 COMPLETED 0:0**，gpu035，2CPU/2GiB、0GPU，
wall17s（34 allocated CPU-seconds = 0.009444 CPUh），MaxRSS597,488KiB。
15 tests passed in13.55s；JUnit核实15项、0失败、0错误、0跳过。
[未改写的原始stdout](packed_gru_cpu_R1.log)保留每例容差、误差和近零gain。

13个含反向/AdamW的案例覆盖6种小形状×FP64/FP32，以及201候选/hidden128的FP32；
另2个是近零gain的前向诊断。最大FP32绝对差：score2.98023e-8、
feature gradient1.39698e-9、parameter gradient1.49012e-8、
参数增量1.49012e-8、optimizer state1.86265e-9；链loss差0。
15个构造案例均未出现正gain符号翻转，但这不保证真实候选中永不翻转。
相对误差分母下限为对应atol，原始日志同时保留该量；不能只看相对误差或只看loss。

判断：隔离实现通过这组CPU数值验证，可以作为未来GPU等价/完整步计时的候选；
目前没有真实TRAIN GPU测量、完整joint loss比较或速度优势。原R4代码未修改，
新原型不在其导入链中，不因本测试通过就热替换或扩预算。
