# M27 WT 原生／缓存重放结果

2026-09-26。结论：**两工具全 chr22 缓存重放等价性通过；得到后续条件干预分母，没有突变效应结果。**

## 执行与事前固定

协议与代码在首次推理前公开于提交 17c07711d11e394edacbedad54fc062d3ab1a0f9。执行 [wt_replay_contract.md](wt_replay_contract.md)、[M27-WT-REPLAY.sbatch](../../sbatch/M27-WT-REPLAY.sbatch)，输入/模型/精度/窗口/原生自适应控制器均未更换，未重试。

Slurm 13227533：gpu035，COMPLETED 0:0，15分49秒，1×RTX3090、8CPU、32GiB RAM。分配GPU时间为 949/3600 = **0.2636 GPUh**，不是按GPU利用率折算。MaxRSS 14,726,532 KiB（约14.04GiB）。5项针对性测试全部通过（2.23秒），详见 wt_tests.xml。0训练、0mutant推理、0Setaria访问。

完整运行目录 outputs/M27-WT-REPLAY-R1 保留在 Baobab。公开 [wt_result.json](wt_result.json) 是未改动的紧凑输出；原始分数、输入缓存、两次整染色体预测、完整过滤链轨迹和日志不上传。

## 重放结果

| 指标 | ANNEVO | Tiberius |
|---|---:|---:|
| native / replay CDS transcript数 | 530 / 530 | 575 / 575 |
| chr22、strand、完整CDS链、phase多重集合完全相同 | 是 | 是 |
| native-only / replay-only链记录 | 0 / 0 | 0 / 0 |
| 21位点中的唯一WT精确链 | 12 | 10 |
| WT非精确 | 9 | 11 |
| 重复exact导致的身份歧义 | 0 | 0 |

共同 WT 精确位点为 **9**，不是任意新命中数门槛，也不是跨工具准确率排行榜。

ANNEVO：原生FP32前向生成双方向 float16 HDF5，各为50,818,468×15；文件3,049,112,560 bytes。两个原生解码进程正常退出，未出现发布代码可能吞掉的worker异常；第二次仅使用同DNA和同HDF5，没有第二次神经前向。

Tiberius：native共4次神经请求，initial正/负链各1次，repredict正/负链各1次；重放逐次元数据和实际单热输入完全一致、完整消费请求，神经前向次数为 **0**。初始shape各128×400050×5；重推理分别6×400050×5与8×400050×5。缓存8,641,081,024 bytes（约8.05GiB），未触及20GiB上限。原生两个DNA后过滤的顺序和全部前后链也完全相同：

- 最短coding长度过滤：706 → 585；
- internal-stop过滤：585 → 575。

此验证不等于神经推理跨硬件确定性测试，也不证明 mutant改变原生请求后仍可使用WT分数；本阶段没有运行此类交叉条件。

## 全部21位点出账

“是”表示唯一的完整坐标链、链向和phase与登记参考一致，不是按基因名或邻近距离分配。完整IDs见JSON。

| 登记名称 | ANNEVO WT exact | Tiberius WT exact |
|---|---|---|
| POTEH | 是 | 是 |
| XKR3 | 是 | 是 |
| ENSG00000283809 | 否 | 否 |
| GSC2 | 是 | 是 |
| LRRC74B | 是 | 是 |
| GSTT4 | 是 | 否 |
| CRYBB2 | 是 | 否 |
| CRYBB1 | 是 | 是 |
| CRYBA4 | 否 | 是 |
| NEFH | 是 | 是 |
| ENSG00000285404 | 否 | 否 |
| RFPL2 | 否 | 否 |
| SLC5A4 | 是 | 是 |
| APOL5 | 否 | 否 |
| CACNG2 | 是 | 是 |
| TPTEP2-CSNK1E | 否 | 否 |
| APOBEC3A | 是 | 否 |
| ENSG00000284554 | 否 | 否 |
| ENTHD1 | 是 | 是 |
| SHISAL1 | 否 | 否 |
| TTLL8 | 否 | 否 |

没有替换readthrough、未命名或WT不正确位点。各工具自身 W_m 为12和10，共同9只用于同位点跨工具比较；不能把12/21或10/21外推成全基因组或物种泛化准确率。

## 实测成本与下一步决策

外层完整进程 walltime：ANNEVO prediction 65.96秒、native decode 63.02秒、replay decode 49.57秒；Tiberius native含容器/导入523.44秒，replay含容器/导入200.45秒。Tiberius内层run_tiberius计时分别438.81/117.19秒，不能用内层数代替全部作业成本。

WT输出目录最终11,698,501,778 bytes；准备数据239,212,157 bytes、专用权重目录62,030,362 bytes，合计约11.18GiB，低于50GiB上限。既有SIF和ANNEVO权重不是本阶段新下载。

按本次完整native进程成本线性估算（不是承诺上界）：

- 两工具都运行全部21位点的syn/PTC：新增42×(128.98+523.44)秒，约7.61GPUh；加WT约7.88GPUh，几乎没有8GPUh总预算余量，不应直接默认能完成。
- 若严格限于已事前定义的各工具 W_m 条件estimand：24次ANNEVO+20次Tiberius，约3.77GPUh；加WT约4.03GPUh。未进入W_m的原登记行仍保留并标明“WT不精确，条件响应未测”，不能当阴性。

下一步将向Pro第十轮提出这个明确取舍，同时冻结唯一目标对应、bypass/终止/消失/重启/不可判定的分类与配对差分母。该建议**尚未成为已执行的mutant协议**。不重新筛参考位点，不选择有利WT阈值，不开始分量交叉、训练、新模型或Setaria测试。

## 科学边界

通过的是同输入的工程／测量有效性前提，不是“PTC诱发绕行”证据，也不产生新的Nature Communications中心主张。ACE/SGRF先行工作、synthetic与真实生物学的区别、模型训练暴露限制继续有效。后续若有明确PTC特异响应，还需可复现分量干预及独立真实等位证据才能讨论外推和实用价值；不能凭本次重放或零歧义就宣称发表目标完成。
