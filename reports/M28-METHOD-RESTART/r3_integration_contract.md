# M28 R3：真实GLM特征与B1自由候选接线

2026-09-30。在用户明确的新模型研发范围内推进有限工程运行；不是正式首训或性能筛选。

## 先固定真实特征 smoke

- 固定本机缓存现有官方GENERanno1.2B CDS annotator revision b0483c23b6b63787b61a6d3a204a9b517d6ba345；仅原backbone，丢弃分类头，无旧LoRA，eval+requires_grad=False，BF16、SDPA、use_cache=False，离线读取，不下载或安装。
- 从R2共享随机抽样顺序中，取每species/chromosome/strand首个窗口，恰好12个；按层选择，不按标签/预测挑例子。只读原train seqids，不用val、test、Setaria。
- 每24,576bp窗口分四个6,144bp局部块，逐块前向。核对无特殊token的四块token拼接与整窗token化一致，1024token/块、4096token/窗；保存对齐后的BF16特征供C0/B1共享。
- 1张24GB RTX3090、2CPU、32GiB RAM、15min（上限0.25GPUh），缓存上限1GiB。私有GPU现满；共享3090可作短作业，具体调度由实时Slurm决定，不以nvidia-smi猜可用。
- 记录模型加载/逐块前向时间、GPU峰值、缓存大小和冻结状态。不执行optimizer step，不测DEV指标，不以吞吐推断新方法有效。
- 失败时保留原输出，最多一次具体工程修复后在独立目录复验；不得改数据、模型revision或变成参数扫描来救性能。若预算/资源不成立则报告并等待资源或缩小非科学实现开销，不扩大正式训练授权。

## 并行CPU实现

执行补充（在候选/梯度结果产生前写定）：两物种各一个CPU作业，private-teodoro-gpu，各2CPU/16GiB/20min；每物种使用上述6个固定窗口。输出为新的outputs/M28-B1-INTEGRATION-R3/<species>。不执行optimizer step。

工程候选预算：每1,024bp/类型最多8端点；每起点与carry状态、每结束类型最多8个ORF合法exon边；每donor最多8个intron链接；每事件/未完成密码子保留2条beam；每kb最多8条送入链头的完整候选，最多128个exon。GT/GC-AG，最小intron总长5bp。本次不扫描任何限额，未训练头低召回不作为增加K的理由。搜索按加性提案分数剪枝，再以非加性链头评分；不是非加性全局最优。

联合loss为局部fine/grouped似然＋端点/连接/链的已知正负类均衡BCE，系数暂均1用于接线。每窗各项先归一化，再施加R2采样权重；正式pilot在吞吐后另行一次冻结。自由候选先落盘，之后才构造参考监督；只有两物种真实小批各头和共享层均有对应梯度、mask零贡献、正负candidate/null方向正确，才判定接线完成。

B1候选API只接受DNA、局部输出与连接分数；不接受参考。先生成并保存剪枝后自由候选，再在训练分支独立标注/注入真实链。携带跨exon未完成密码子，记录端点/连接/beam/完整链剪枝数量；beam仅近似搜索，不保证非加性整链分数全局最优。

端点masked BCE、连接正/参考负/未知监督、局部细分类/分组似然、完整链/null损失分别接线，每窗归一化后乘抽样权重。真实小批确认各头与共享层梯度，unknown/mask不贡献损失，冻结特征不接收梯度。空候选零loss不能记为拒绝能力。

候选数量、loss系数及操作点目前用于工程 smoke，不根据召回扫描；实际吞吐后一次冻结正式公平pilot合同。旧M27关闭，新主评价采用独立CDS-assessable政策且保留模型不支持的合格参考；训练eligibility不删评价分母。
