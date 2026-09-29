# M28 packed-GRU GPU R1：数值验证未通过，未进入速度阶段

2026-09-30（Bangkok）。沿用[运行前合同](packed_gru_gpu_contract.md)，没有放宽容差、
重跑或热替换R4。CPU15项测试通过没有被冒充GPU数值等价。

## 实际执行与失败

Job **13293653 FAILED 1:0**，private/gpu035，1×RTX3090、2CPU/8GiB，
申请8分钟（低于合同15分钟上限），实际49s，即 **0.013611 allocated GPUh**；
MaxRSS1,957,380KiB。命令：

```bash
sbatch --parsable --time=00:08:00 sbatch/M28-PACKED-GRU-GPU.sbatch
```

使用只读B1 step001536及AdamW快照，固定TRAIN首窗
`arabidopsis_thaliana|NC_003070.9|29786112|-`。
两个warmup和两个数值对照各完成一步，四次均恢复相同初态并丢弃更新；
无新权重保存、无连续续训、未读DEV/test/Setaria。尚未到水稻窗口。

第一窗比较ChainHead返回评分分项时触发原定FP32容差
`atol=2e-6, rtol=2e-4`：192元素中17个不符，最大绝对差
0.0009055137634277344，最大相对差0.002578454790636897。
原stderr没有注明失败的具体score key，不能把该差值明确写成最终gain误差；
该记录缺口保留，不靠事后推断填成已观察事实。

已通过的完整窗口数0；48步性能阶段执行0；summary.json未产生。
因此没有完整步/链头速度结论，也未完成所有梯度和optimizer状态比较。
这是实现数值准入失败，不是B1模型性能失败或结构假设被否证。

## 已保存输出能说明什么

第一窗两个实现的完整自由候选字典相同，均192条候选/训练链。
local、endpoint、link loss逐值一致；chain loss为
0.056456465274095535（原版）和0.05645260214805603（packed）。
weighted total为0.39778271317481995和0.39777854084968567。
warmup与正式对照内各实现重复得到相同这些数值。

这些证据把本例差异范围收窄到链评分相关的执行路径，而非抽样或候选改变；
尚未证明具体CUDA/cuDNN原因。loss差小不抵消逐候选分项越界。
[机器可读结果](packed_gru_gpu_R1.json)、
[原始stderr](packed_gru_gpu_R1_error.log)、
[未改写的四步trace](packed_gru_gpu_R1_steps.jsonl)保留原证据。

## 精度线索与当前决定

失败后在同一generanno环境做只读版本查询（登录节点，无GPU计算）：
PyTorch2.5.1、CUDA build12.1，默认cudnn.allow_tf32=True、
cuda.matmul.allow_tf32=False。这不是失败作业内事前保存的flag快照。
[PyTorch v2.5.1官方CUDA说明](https://github.com/pytorch/pytorch/blob/v2.5.1/docs/source/notes/cuda.rst)
明确指出TF32影响使用相关计算的GRU/LSTM；FP32张量类型不等同禁用TF32。
这使精度/内核路径差异成为可检验假说，但目前不能宣称它就是原因。

当前决定：不接入packed实现，不放宽容差，不重跑本R1，不更改R4的精度设置。
若继续定位，应将“同实现可重复性、具体失败分项、固定精度下两实现差异”
写成独立且有限的诊断问题，明确它不是把本次失败改判为通过。
原R4仍按原预算等终态；即使将来数值验证通过，速度和正式方法收益仍各需实证。
