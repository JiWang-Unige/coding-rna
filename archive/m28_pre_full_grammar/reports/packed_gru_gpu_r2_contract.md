# Packed链头R2：精度作用域修正后的完整一步验证

2026-09-30，执行前固定。R1已失败并停止，D1另取得flag相关的实际证据。
在用户已授权的迭代研发范围内，本次只检验D1支持的一项明确实现修正，
不是原R1无改变重试、放宽容差或续跑R4。原代码/产物/失败结论全部保留。
该修正暂不接入任何正式fit/infer入口，科学目标仍是完整基因模型而非工程优化论文。

## 实现、参照和计算边界

- 参照完全沿用原ChainHead.forward和原shared forward，cudnn.allow_tf32=True、
  matmul.allow_tf32=False；没有把参照一起换成新精度求一致。
- 新候选仅在packed链头调用期间将cudnn.allow_tf32设为False，并在finally恢复。
  链头之外不改shared GRU/卷积、候选生成、pooling、Linear/GELU、读出或矩阵设置。
  前向结束及一步反向/optimizer后均检查原flag已恢复。
  不猜测反向是否等价，直接比较其梯度、clip和AdamW结果。
- 读原step001536及AdamW快照；原paired draws前12窗（Ara6/rice6）、
  原feature/标签/importance weight/loss/候选budget不变。每步重置初态，更新丢弃，
  不生成可评模型，不作连续续训，不读DEV/test/Setaria。
- 1×RTX3090/24GiB、2CPU/8GiB、8分钟硬限（≤0.133333 allocated GPUh），
  private即时可用则优先。预估数分钟完成，输出M28-PACKED-GRU-GPU-R2，
  文本预计<10MiB，不复制特征缓存或权重。
- 此优化分支GPU累计上限仍15分钟：R1实际49s＋D1实际35s＋R2至多480s，
  最坏564s，不把更换ID当作费用重置。CPU验证17s另列。R4的4.763333GPUh另记，
  已关闭未完成，不挪其阶段余额或扩其上限。

## 原容差下的完整步检验

两实现各warmup一步，然后12窗按AB/BA交替，每次从相同checkpoint/AdamW初态。
自由候选先落盘再建参考监督；候选、注入池、标签严格相同。
额外要求shared features逐值相同，避免作用域泄漏。

沿用FP32 atol2e-6/rtol2e-4：逐score key、四项loss/total、shared梯度、
全部参数梯度、clip前norm、一步参数增量、AdamW矩状态和step。
GPU起始版本/flags先落盘；每步trace即时flush；完整记录当前窗口各项误差与
具体失败key后再判定，失败不进入速度阶段。数值门槛或实际自由候选的gain符号/
非重叠选择任一不符均拒绝，不放宽容差、不继续扫flag，也不自动再次修补重试。
近零稳定性只能按实际出现的gain解释；12窗不能保证完整训练轨迹一致。

## 只有通过后才测吞吐

同一GPU两轮×12窗×两实现，按窗和轮反转AB/BA，共48个未插桩完整步；
加warmup与等价共74个被丢弃的一步更新。保持加载/转换/DNA/H2D/shared、
原提案、自由候选JSON写入flush、监督/反向/clip/optimizer、trace写入flush。
复位权重/AdamW与启动耗时在完整步计时外，另报Slurm总wall。
这是固定小批热缓存配对比较，不冒称冷读或完整连续训练成本。

只有数值/离散输出均通过，且两轮逐窗baseline/candidate比值的中位数都>1、
候选完整步中位数低于原版，才保留它用于下一次完整试验的预算估计。
不把小样本计时当普遍性能或论文贡献。数值失败/无整步收益就结束本优化分支；
OOM、非有限、输入错配或超时也停止，无自动重试。新完整模型实验需另冻结
可完成的运行方案，历史成本和未完成状态不抹除。
