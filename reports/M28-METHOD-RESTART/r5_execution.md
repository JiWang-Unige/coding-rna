# M28 R5执行记录：原实现的完整比较正在运行

2026-09-30。按[执行前合同](r5_execution_contract.md)启动，
不是模型性能结果；也不是被关闭的GRU优化重试或R4中间权重续训。

| 阶段 | Job | 提交后实际状态 | 资源硬限 |
|---|---:|---|---|
| B1从头拟合 | 13294132 | RUNNING，private/gpu035 | 1×3090、2CPU/8GiB，3h |
| C0完整DEV推理 | 13294133 | RUNNING，private/gpu035 | 1×3090、2CPU/8GiB，2h |
| B1完整DEV推理 | 13294138 | PENDING，afterok B1拟合 | 1×3090、2CPU/8GiB，4h |
| 组装/GFF3/共同评分 | 13294139 | PENDING，afterok两臂推理 | 0GPU、2CPU/32GiB，2h |

上表为提交后约42s快照，不是终态。B1拟合与C0推理已确认同时运行。
新增GPU上限9h，CPU独立分配上限4CPUh；旧R4的4.763333GPUh和GRU分支
134GPU秒保留，不因为新ID省略。非GPU工程与更早R3成本继续在各自记录中。

C0复用原已完成4,608更新的step004608，不重复拟合；B1从同seed0/相同draw开始
重新完成4,608更新，原逐链实现不变。两臂都在新目录完整推理8,690窗，
不拼接旧R4部分预测。features与C0目录为明确相对symlink指向已完成资产，
已核实实际解析到R4相应目录，不复制数据、不写入这些目标。

这次代码变化仅为fit接受合同路径元数据，以及评价报告从实际pilot目录和
最终checkpoint路径注明来源。模型、损失、候选、解码、合并和评分函数均未改；
AST/shell语法通过，变更逐行核对。native ANNEVO依赖仍为
37bdd9aa62ddf24fa55941fb827061f7ed49ce53，src/HMM.py无工作树改动。
存储检查显示共享文件系统170TB可用；本轮新增产物限2GiB。

实际命令：

```bash
sbatch --parsable sbatch/M28-PILOT-R5-FIT-B1.sbatch
sbatch --parsable --job-name=M28R5INFER-C0 --time=02:00:00 sbatch/M28-PILOT-R5-INFER.sbatch C0
sbatch --parsable --job-name=M28R5INFER-B1 --dependency=afterok:13294132 sbatch/M28-PILOT-R5-INFER.sbatch B1
sbatch --parsable --dependency=afterok:13294133:13294138 sbatch/M28-PILOT-R5-EVALUATE.sbatch
```

完整调度记录在远程outputs/M28-PILOT-R5/execution.json。
任何失败按合同保留并停止，不延时或自动重提；只有最后checkpoint和全DEV都完整，
才生成正式比较。原R4两个DependencyNeverSatisfied作业未修改或取消。
目前尚无新的DEV准确率、独立泛化或论文主张；test/Setaria继续封存。

## 阶段完成：C0全DEV推理（2026-09-30）

Job **13294133** 的作业及batch均为 `COMPLETED, ExitCode=0:0`，
`infer_C0/STATUS` 为COMPLETED；正常生成summary及四个scope文件。
Slurm原始起止时间为2026-09-29T23:58:41至2026-09-30T01:04:59，
实际分配 **3,978s = 1.105 GPUh**（1×RTX3090、2CPU、8GiB），
低于本阶段2h硬限。batch MaxRSS为1,307,476KiB，stderr为空。
这是已完成C0推理的新增成本，不包含仍在运行的B1或更早研究成本。

| DEV范围 | 方向 | 窗口数 | 窗口内链记录数 |
|---|---|---:|---:|
| Arabidopsis NC_003074.8 | + | 1,909 | 5,807 |
| Arabidopsis NC_003074.8 | - | 1,909 | 6,005 |
| Rice NC_089041.1 | + | 2,436 | 2,905 |
| Rice NC_089041.1 | - | 2,436 | 2,495 |
| 合计 | | **8,690** | **17,212** |

完成后直接解析全部四个JSONL（合计3,090,312字节）：与summary窗口列表逐项核对
window ID及顺序、无重复ID、window_index连续、链与分数等长、分数有限、
每窗reference_used=false，均通过；各scope链记录数与summary一致。
运行代码在推理前核对固定DEV网格，逐窗核对缓存revision及几何，结束前核对8,690窗。
summary确认使用原已完成C0 `outputs/M28-PILOT-R4/C0/step_004608.pt`，
4,608更新，冻结revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345`；
reference_used/test_setaria_used均为false。未拼接R4部分预测。

**17,212是重叠窗口内的链记录，不是去重基因数、精确链数或生物学真阳性。**
这里只确认推理产物完整及计算成本，不是模型性能PASS或独立复现。
B1拟合继续原作业13294132；13294138与13294139仍等待原afterok依赖。
不提前生成C0单臂排名，不改变任何阈值、候选预算、训练步数或资源上限；
完整方案的科学判断继续等待B1最终checkpoint、全DEV推理及共同评分。
