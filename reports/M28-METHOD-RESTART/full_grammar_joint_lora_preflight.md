# 原 TRAIN 的联合窗口低秩适配器预检

2026-10-12。五组人工检查13499633和原生预检13499640均已完成并经主对话核验。两条固定真实TRAIN view支持结构损失到适配器的梯度、24GiB GPU上的显存可行性和更新后的训练/raw一致性，不构成完整链准确率收益。每物种只用B预检已固定的第一个fresh TRAIN draw，从原最终4608头及AdamW状态独立恢复，各做一次可丢弃更新；不保存checkpoint，不启动正式fit，不解析、保留或使用DEV或封存DNA及其注释。并行B、C、D的代码、配置、作业和冻结评价范围不变。

## 固定容量变化

已缓存的GENERanno revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345`有26个编码层，hidden2048、32个query头、4个KV头。最后4层（零基22至25）各在q_proj和v_proj上添加rank4、alpha4、scale1的低秩残差，共102400个FP32参数。原BF16权重、embedding、RoPE、其他投影、MLP和归一化均冻结；整个编码器保持eval，无dropout或gradient checkpointing。使用原生SDPA的Flash Attention后端，不改双向注意力或位置编码。

每个新投影为原BF16输出，加上FP32计算的`B(Ax)`经BF16舍入后的残差。A按CPU generator seed0在正负输入维度平方根倒数之间均匀初始化，B为零，原模型嵌套于这一初始状态。第一步A梯度应为零、B梯度应非零；这不是梯度失联。后续是否能改善完整链准确率是另一个问题，预检不能回答。

## 可微路径和固定观察

原joint24576编码函数有inference_mode装饰器，因此只训练C0头。新增路径保留完全相同的24576bp窗口、4096个六聚体token、转录方向、非ACGT转N、短窗右侧N补齐、BF16隐藏表示转FP32及最近窗口中心ownership，只去掉冻结特征的梯度屏障。旧函数不修改；raw推理继续调用旧joint caller，但传入已装入同一适配器参数的编码器。

两条真实view仍来自同一fresh plan，不能按新loss换样本。fine／coarse／unknown、epsilon0.1、全部合法gene-set分区、cap128、原fine+coarse分母、FP64分区和importance1保持B的定义。头部AdamW从4608至4609，新的16个适配器张量从0至1；两类参数共同norm1裁剪，沿用lr0.0003、weight_decay0.01、betas0.9／0.999及eps10⁻⁸。更新后状态仅保留在内存中。

## 决定是否可用的验证

五组人工检查覆盖零B时的逐元素嵌套、首步A／B和上游梯度、非零残差计算、只训练固定投影、适配器状态完整加载，以及joint窗口的文字token化、padding和可微BF16→FP32路径。不读取研究DNA、标签、模型或cache；经Slurm CPU作业执行。

原生GPU预检先保存两个固定view的冻结父势函数，再接入适配器。零参数的native势、真实NLL和role VJP必须与父模型逐元素相同；真实归一化结构损失必须到达每个B投影，全部头和适配器梯度有限，原始编码器参数无梯度。确认两类optimizer年龄和实际非零B更新后，核对独立头与适配器内存恢复后的训练/raw势及实际context MAP一致，并运行14bp caller的两个方向。短输入若无CDS，只证明执行、token与方向完整，不证明非空输出坐标或准确率。

GPU资源固定为单24GiB RTX3090、2CPU／16GiB、20分钟、no-requeue。记录每view实际native窗口数、时间、峰值allocated和reserved显存、投影梯度、NLL及更新年龄；不保存特征、DNA、势矩阵或权重副本。失败时保留输出，先定位具体原因；不根据NLL或显存事后替换样本、rank、层数、观察模型或训练参数。

通过后才准备同起点的正式对照及原TRAIN端点；本预检本身不向现有六方法DEV家族添加候选，不证明LoRA新颖性、泛化或论文达标。

## 人工验证结果与原生作业

人工验证13499633的root／batch／extern均COMPLETED、ExitCode0:0，2026-10-11 18:34:45至18:34:51UTC。五组全部通过，错误与失败均为零，脚本实耗0.435749秒，batch MaxRSS473576KiB。未读取研究DNA、参考、checkpoint或原生模型，未使用DEV、封存数据或cache，optimizer更新数为零。

原生预检13499640经afterok13499633提交，科学命令为`python -u scripts/experiments/M28-METHOD-RESTART/preflight_full_grammar_joint_lora.py --checks-job 13499633`。提交前Slurm展开的当前用户作业为44项，低于80项限制；gpu034有1个未分配的RTX3090，gpu035已满，近期reservation不涉及本作业节点。沿用private-teodoro-gpu和固定24GiB资源，不根据GPU闲置改变科学配置。

## 原生预检实际结果

13499640的root／batch／extern均COMPLETED、ExitCode0:0，2026-10-11 18:35:27至18:37:17UTC，实际gpu034／单RTX3090／2CPU／16GiB／20分钟／Requeue0。脚本实耗93.370838秒，batch MaxRSS4545872KiB，stderr为空。两条view各有49152bp context、3个真实joint24576窗口，原已知观察分母均为24576。

| 固定 TRAIN view | 零适配器 raw NLL | 归一化梯度裁剪前 L2 | 更新后 B 的 L2 | GPU reserved 峰值 GiB |
| --- | ---: | ---: | ---: | ---: |
| Arabidopsis NC_003070.9 core14389248 负链 | 0.00149551027688 | 0.0000116859574 | 0.0000287567 | 6.806640625 |
| Rice NC_089046.1 core14462976 负链 | 30.7035429458 | 0.440460235 | 0.0506390440 | 6.810546875 |

两次各从同一原4608头和完整moment、同一seed0适配器初始化独立开始；实际头年龄均4608→4609，全部16个适配器张量年龄均0→1。所有8个B投影的真实归一化梯度非零且有限，A梯度严格为零，符合零B初始化；全部原编码器参数冻结且无梯度。总峰值allocated为6828795904字节（6.359811783GiB）、reserved为7312769024字节（6.810546875GiB），说明这两个固定view可在24GiB GPU上完成本路径，不等于所有训练样本或未来运行均已验证。

零适配器的真实native十五角色势、NLL和role VJP与冻结父模型逐元素相同。更新后独立恢复的头和16个适配器张量，经可微训练路径与原joint raw caller得到完全相同的native势及实际context MAP。两次14bp caller均完成双方向且预测CDS为零，因此不将其作为非空坐标恢复或准确率证据。输出只有started.json、summary.json、STATUS和JOBID，均非空；未保存checkpoint、特征、势矩阵或DNA副本，未使用DEV或封存数据。

该结果使同起点的正式TRAIN对照值得准备，但尚未启动适配器正式fit，也未修改现有六方法DEV协议。每view耗时包含父分区及VJP比较、参数恢复和多次raw／MAP运行，不能据此直接外推正式训练速度。研究goal保持未完成。
