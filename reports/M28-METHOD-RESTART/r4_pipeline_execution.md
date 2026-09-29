# M28 R4：完整DEV推理与汇总链已接入

2026-09-30（Bangkok）。本轮完成的是运行接口与计分验证，**不是新模型性能结果**。
沿用[首训冻结合同](r4_pilot_contract.md)，没有改数据、loss、候选预算、更新数、
checkpoint或操作点。此前“DEV runner待接”的准备项已落实。

## 本轮代码

- [infer_pilot.py](../../scripts/experiments/M28-METHOD-RESTART/infer_pilot.py)：
  只读取最后step004608 checkpoint和同一冻结特征/DNA，覆盖全部8,690个DEV
  方向窗口；缺窗、重复窗、越出两条允许染色体或checkpoint不符则停止。
  不读取参考标签、不注入参考候选。C0保存原生链及gain；B1保存实际剪枝后
  自由候选、完整链score及同候选additive-only score。
- [export.py](../../src/m28/export.py)：将全局半开CDS链导出GFF3，
  GFF起点转1-based，phase按转录方向累计编码长度计算。输出不冒称UTR/全isoform。
- [evaluate_pilot.py](../../scripts/experiments/M28-METHOD-RESTART/evaluate_pilot.py)：
  只有两臂完整拟合与完整推理均完成才运行正式比较。候选文件生成完毕后，才读取
  冻结参考；采用已验证共同assembler，分别计算自由候选→owner→冲突选择的
  真正参考链损失，以及新7,728分母和旧6,450副表、物种宏平均和micro汇总。
  另输出背景span/CDS覆盖、全背景链/Mb、CDS结构诊断与GFF3往返核对。
  参考不一致不等同生物学假基因；B1−C0不单独归因于非加性GRU。

C0固定HMM依赖当前HEAD确为
`37bdd9aa62ddf24fa55941fb827061f7ed49ce53`，`src/HMM.py`无工作树改动。

## 实际验证

**Job13292518 COMPLETED 0:0**，private-teodoro-gpu，2CPU/2GiB、0GPU，
wall56s，MaxRSS768,568KiB；**10 tests passed in52.23s**。
其中新测试覆盖正负链phase/坐标GFF往返、空GFF、DEV缺窗/重复/越界拒绝、
实际链键的逐阶段损失、宏平均与micro区别，并复用共同标尺和assembler测试。
这是针对实现失败的验证，不代表预测正确率或独立科学复现。

## 已调度的完整执行链

下表为UTC 2026-09-29 18:58–19:01附近的实时调度快照；不是终态结果。

| 阶段 | Job | 成功依赖 | 资源上限 | 快照状态 |
|---|---:|---|---|---|
| 共享特征 | 13292334 | 无 | 1 GPU/2CPU/16GiB，2h | RUNNING |
| C0拟合 | 13292356 | 13292334 | 1 GPU/2CPU/16GiB，2h | PENDING Dependency |
| B1拟合 | 13292357 | 13292334 | 1 GPU/2CPU/16GiB，2h | PENDING Dependency |
| C0完整DEV推理 | 13292540 | 13292356 | 1 GPU/2CPU/16GiB，1h | PENDING Dependency |
| B1完整DEV推理 | 13292541 | 13292357 | 1 GPU/2CPU/16GiB，1h | PENDING Dependency |
| 两臂组装/GFF3/评价 | 13292543 | 13292540、13292541 | 0GPU/2CPU/32GiB，2h | PENDING Dependency |

全部使用private-teodoro-gpu。新增两项推理的1h+1h正好位于原合同
“合计GPU推理≤2h”内，整轮GPU上限仍2+2+2+1+1=8h。CPU解码发生在GPU
推理作业内时，等待/解码wall照计GPU分配用量。最终CPU汇总最多4 allocated CPUh，
本轮测试实际112 CPU-seconds，均单列，不用CPU运行冒充GPU吞吐。

实际新增提交命令：

```bash
sbatch sbatch/M28-PIPELINE-TEST.sbatch
sbatch --job-name=M28INFER-C0 --dependency=afterok:13292356 sbatch/M28-PILOT-INFER.sbatch C0
sbatch --job-name=M28INFER-B1 --dependency=afterok:13292357 sbatch/M28-PILOT-INFER.sbatch B1
sbatch --dependency=afterok:13292540:13292541 sbatch/M28-PILOT-EVALUATE.sbatch
```

共享缓存快照已完成4,288/10,181窗，新增71,851,756,292 bytes；有进展、
未报告失败。该数字仅为时点进度，不是最终缓存规模或GPU成本。
原始执行目录为`outputs/M28-PILOT-R4`，其中execution.json记录准确job与依赖；
测试证据在`outputs/M28-PIPELINE-TEST-R4`。缓存、预测明细、GFF3仍留远程。

## 共享缓存终态与训练启动（UTC 2026-09-29 19:27–19:31）

共享缓存job13292334已 **COMPLETED 0:0**，gpu035；wall3,010s，
即 **0.836111 allocated GPUh**，MaxRSS8,057,680KiB。
[完整摘要](r4_features.json)与索引共同确认10,181窗：
1,491个唯一TRAIN窗＋8,690个完整DEV方向窗，其中12个TRAIN缓存复用R3。
新增缓存170,874,523,021 bytes（159.1393GiB），低于180GiB缓存上限。
纯backbone前向2,445.9568s，提取循环2,932.4193s；这些不替代Slurm分配时间。
GPU峰值allocated2,374,624,768B，reserved2,489,319,424B。
本阶段optimizer更新0、backbone可训练参数0、test/Setaria使用false。
上述用量是R4新作业用量；复用的R3运行成本保留在原R3记录，不冒称免费重建。

C0 job13292356已在gpu035执行，日志已记录至少192/4,608次更新；
这只是运行证据，不是DEV性能。B1的缓存依赖已解除，但19:31快照仍排队。
原成功依赖链、checkpoint、样本、loss、操作点和总8GPUh上限均未改变。

依据实时Slurm分配，gpu035仅剩7,384MiB未分配内存，原16GiB请求不利于并行。
先让未启动B1两作业可选private-teodoro-gpu及shared-gpu，并保留RTX3090类型；
后根据C0实际MaxRSS2,254,880KiB、此前B1含反向测量1,782,104KiB，
把B1训练/推理主机内存请求降为8GiB，使其可利用gpu017当时的剩余10,968MiB。
这是调度资源请求调整，不是降低GPU显存、增加预算或改变科学参数。
原job ID不变，没有取消、复制或重跑任务。实际命令：

```bash
scontrol update JobId=13292357 Partition=private-teodoro-gpu,shared-gpu ExcNodeList=gpu023,gpu024,gpu034,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
scontrol update JobId=13292541 Partition=private-teodoro-gpu,shared-gpu ExcNodeList=gpu023,gpu024,gpu034,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
scontrol update JobId=13292357 MinMemoryNode=8192
scontrol update JobId=13292541 MinMemoryNode=8192
```

本机Slurm拒绝了首次以“8G”为MinMemoryNode值的请求，未改变作业；
改用8192（MiB）后，scontrol确认ReqTRES为mem=8G。
排队是否消除以实际RUNNING状态为准，不用调度意图声称已并行。

随后UTC19:33附近的squeue确认两臂**确实并行RUNNING**：
C0 job13292356在private/gpu035，B1 job13292357在shared/gpu017。
B1初始MaxRSS2,452,700KiB，低于8GiB请求；最终训练/DEV结果仍未完成。

## 下一判断

依赖链是普通Slurm afterok，不恢复退役自动研究框架，也不新增审批。
继续核对任务终态和完整产物，再回答B1是否在两物种相对C0改善、改善来自覆盖
还是输出选择、背景负担如何。任何超时/未完成先作为执行状态处理，
不拿半成品排名，不因希望获得Nature Communications结论而改冻结规则。
尚未取得新模型准确率、独立泛化或投稿级主张；test和Setaria继续封存。
