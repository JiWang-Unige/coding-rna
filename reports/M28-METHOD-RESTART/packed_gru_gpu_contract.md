# M28 packed-chain GRU：独立GPU工程验证合同

2026-09-30（Bangkok），在执行前固定。已完成的CPU15项测试不能代替GPU等价或速度。
本阶段是用户已授权新模型研发中的有限实现验证，不是新的科学模型筛选，
不以Pro建议作为授权。R4原作业保持原实现/时限/依赖，不热替换、重启或续训。

## 范围与资源

- 1张RTX3090/24GiB、2CPU、8GiB RAM，硬限15分钟（最多0.25 allocated GPUh）。
  private可用且完成时间相近则优先；只按live Slurm路由，不改科学内容。
- 输出`outputs/M28-PACKED-GRU-GPU-R1`，工程成本单列并计入研究累计消耗。
  R4各阶段上限及总8GPUh不重置；本作业不能补齐R4未完成的更新或DEV窗口。
- 只读已保存的B1 `step_001536.pt` 和对应AdamW状态。每一次比较均恢复同一状态，
  完成一步后立即丢弃；不保存任何模型/optimizer，不形成新的训练轨迹。
- 原paired_training_draws.jsonl前12项，已核实Ara6/rice6；顺序不因损失或标签改变。
  只用对应TRAIN特征及R4标签，原importance weight、loss、候选预算和优化器不变。
  不读DEV/test/Setaria，不选epoch、seed、操作点或候选上限。
- 只在本验证进程的模型实例上切换forward函数，不改core.py或R4 runner。

## 数值验证与停止

原逐链GRU对照`src/m28/core.py::ChainHead.forward`；
候选仅替换为已CPU验证的`src/m28/packed_chain_gru.py::packed_chain_forward`。
exon均值/几何/project、加性/非加性/null读出均保留，只有GRU调用批处理。

两个实现先各warmup一步，均从相同checkpoint恢复。随后12窗逐一AB/BA交替：
自由候选先生成落盘，再建参考监督，完整joint_losses四项保持原样。
比较自由候选/连接/剪枝记录、训练候选与标签严格相同；
逐候选各分项score、四项loss及weighted total、共享特征梯度、全部参数梯度、
clip前norm、一步参数增量、AdamW一阶/二阶矩与step。
FP32固定atol2e-6/rtol2e-4，不事后放宽；记录绝对/有atol分母下限的相对误差。
自由候选gain正负翻转和同窗nonoverlap选择变化另记；数值接近不保证离散输出相同。

数值/输入不符、非有限、OOM或超时停止，保留失败产物，不自动重试/扩预算。
任一数值比较失败时不进入性能阶段；选择差异即使处于容差内也须报告，
不能称这些窗口的输出完全等价。一步比较不证明长期训练逐位一致。

## 计时与判断

- 等价阶段另测带CUDA同步的chain forward，作为定位用的插桩读数；
  其整步时间含诊断传输，不能作为速度结论。
- 等价通过后，两轮×12窗×两实现=48个完整训练步，同GPU、每窗/每轮交替AB/BA，
  每步仍重新恢复权重/AdamW。加warmup和等价阶段，共74个被丢弃的一步更新。
- 正式速度计时不插桩链头；包括缓存load/转换/DNA/H2D、shared forward、
  原候选生成、自由候选JSON写入flush、建标签/四项loss/反向/clip/optimizer、
  训练trace写入；不含模型/optimizer复位或初始载入，另报总循环和Slurmwall。
- 此时文件系统缓存已热，只能称配对热缓存实现比较，不能保证连续训练冷读吞吐。
  同时报两实现完整步median/p90、逐窗配对比值、VRAM与actual候选负载。
- 无完整步收益，就不据此安排新完整试验；有收益也仅作为新预算估计依据。
  正式训练必须另外冻结方案并如实累计旧试验/新试验成本；
  不借此让R4半成品进入正式性能表。

相关脚本：`scripts/experiments/M28-METHOD-RESTART/packed_gru_gpu.py`、
`sbatch/M28-PACKED-GRU-GPU.sbatch`。本合同不声称新方法性能、候选召回或期刊创新。
