# M28 固定 TRAIN 见证的原生编码上下文诊断

2026-10-09。检验原生长上下文是否改变固定4608步head的结构信息。使用已完成13429053的原有八个TRAIN见证，不按本轮性能选窗口，不读取DEV或封存数据。Carbon-A来源及上下文消融见[pretrained_exposure.md](../../../../paper/supplement/pretrained_exposure.md)；作者结果是诊断动机，不是本项目差距的因果证据。

## 固定因素和四种编码条件

固定GENERanno原revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、BF16表示、原C0最后4608步FP32权重、原49152 bp语法视图、三个24576 bp神经读取窗口及原nearest-center所有权。只改变编码器的共同注意力输入，不改head、RoPE、标签、语法、epsilon=0.1或阈值。原有三组目标函数训练不变。

| 条件 | 编码器输入 | head和语法读取区间 |
| --- | --- | --- |
| BLOCK6144 | 每个原24576窗口独立编码四个6144块 | 原三个head窗口，拼接到同一49152视图 |
| JOINT24576 | 每个原24576窗口一次共同编码 | 同上 |
| JOINT49152 | 整个原49152视图一次共同编码，裁切特征 | 同上 |
| JOINT98304 | 原视图左右各增加24576真实TRAIN碱基，一次共同编码，裁切特征 | 同上，不扩大语法或参考分母 |

八个原见证均为内部全长窗口。扩展仅按原坐标和固定长度，不依赖参考边界；遇到不满足真实上下文的源变化即停止，不伪造补齐。JOINT24576对BLOCK6144隔离共同编码；更长条件同时增加编码器可见侧翼，不能将其差异全部解释为注意力长度本身。

## 验证与结果解释

先在CPU上检验正负链非对称特征裁切及裁切后所有权拼接；实际GPU前向再逐窗口比较长输入与局部输入的六碱基token身份。BLOCK6144在每个见证core上必须精确重放原BF16 cache，防止编码资产、数值后端或坐标变化被误判为上下文收益。长输入仅允许原SDPA的Flash kernel；不可用时失败，不转用可能耗尽显存的平方注意力矩阵。

GPU保存八份四条件的FP32势函数，不保存新encoder权重或特征缓存。CPU以同一原DNA、fine/coarse/unknown行计算原结构NLL、局部观测NLL及无参考约束的MAP。BLOCK6144须重放已完成原有限视图诊断的NLL（绝对误差不超过1e-6）及精确MAP角色混淆；不匹配即停止解释。逐物种报告相同已知行分母的NLL、15角色参考行和命中，另列原core重叠且完全位于固定视图内的primary链命中。链计数单位是reference-window pair，可能重复且不含所有背景，不能称为unique-gene F1或full-grid精度。

结果只决定是否值得准备TRAIN长上下文训练／接入；若未变或变差，不基于DEV筛选上下文或改当前三臂。旧head是在短块表示上训练，直接换表示的负结果不能证明长上下文经过适配训练仍无效；正结果也不能证明泛化或论文创新。

## 执行

脚本diagnose_full_grammar_native_context.py有checks、native、readout三个模式；同一sbatch脚本按模式选择独立输出目录。CPU检查5分钟，GPU诊断单RTX3090 24 GB、2 CPU、16 GB内存、20分钟，CPU读取30分钟；均使用已有generanno环境与本地固定权重，不下载依赖或模型。GPU约96 MiB势函数输出，不覆盖原数据／结果。三个阶段按afterok连接，提交后保存实际作业号。

本轮实际作业为13452547 → 13452550 → 13452596，三个root及batch/extern均COMPLETED/0:0。两个人工坐标检查均通过，GPU在gpu035完成全部八见证、32条件，CPU在gpu034完成对应32次结构读取。private分区与shared候选的即时可用资源接近，提交时private gpu035有四张未分配3090，选择private；CPU教学预留不涉及这两个GPU节点。

实际提交命令如下，输出目录分别为M28-FULL-GRAMMAR-NATIVE-CONTEXT-CHECKS、M28-FULL-GRAMMAR-NATIVE-CONTEXT-DIAGNOSTIC、M28-FULL-GRAMMAR-NATIVE-CONTEXT-READOUT。

```sh
sbatch --parsable --job-name=M28NATIVECTXCHK --time=00:05:00 sbatch/M28-FULL-GRAMMAR-NATIVE-CONTEXT.sbatch checks
sbatch --parsable --job-name=M28NATIVECTXGPU --time=00:20:00 --gres=gpu:nvidia_geforce_rtx_3090:1 --dependency=afterok:13452547 sbatch/M28-FULL-GRAMMAR-NATIVE-CONTEXT.sbatch native
sbatch --parsable --job-name=M28NATIVECTXREAD --dependency=afterok:13452550 sbatch/M28-FULL-GRAMMAR-NATIVE-CONTEXT.sbatch readout
```

## 已完成的诊断结果

八个core的旧BF16特征全部精确重放；短块NLL与原最终权重有限视图记录的最大绝对差为0.0，MAP fine/coarse混淆也完全一致。全部长输入与局部输入的六碱基token裁切一致。GPU模式用时66.45秒（不含模块导入），峰值分配显存2.93 GiB，输出约91 MiB；CPU读取用时103.62秒。原生98304 bp在24 GB RTX3090上实际跑通，未修改模型实现、RoPE或任何权重，也没有发生optimizer update。

以下NLL按共同已知core行加总后归一化，越低越好；拟南芥分母91868行，水稻98220行。链命中分别以相同的24和14个reference-window pairs为分母，不是F1。

| 编码条件 | 拟南芥 NLL／已知行 | 拟南芥精确链命中 | 水稻 NLL／已知行 | 水稻精确链命中 |
| --- | ---: | ---: | ---: | ---: |
| BLOCK6144 | 0.007770 | 5／24 | 0.008679 | 5／14 |
| JOINT24576 | 0.006536 | 9／24 | 0.005836 | 5／14 |
| JOINT49152 | 0.006316 | 11／24 | 0.006255 | 6／14 |
| JOINT98304 | 0.006585 | 11／24 | 0.005928 | 7／14 |

在这八个固定见证上，共同编码已带来结构变化，不只是特征距离变化。JOINT49152的拟南芥NLL比短块下降18.72%；JOINT98304的水稻NLL下降31.70%。水稻intron1角色命中从200／1539增至761／1539，仍有显著剩余错误。长度收益不单调：拟南芥最低NLL是49152，水稻最低是24576；不能把最大上下文默认宣布为优胜。

每种长条件仍有2／8窗口NLL升高，分别是拟南芥NC_003070.9的25190400负链窗口与水稻NC_089046.1的3416064负链窗口。水稻NC_089044.1的21467136负链窗口在JOINT24576／49152下NLL降低，但精确链命中由4变为3；JOINT98304为4。物种汇总改善不代表每个窗口改善，NLL改善也不能替代精确链指标。

这些窗口最初为边界／支持诊断挑选，非随机独立面板，四个水稻窗口也均为负链。当前结果支持“6144短块确实损失了这些见证可利用的结构信息”以及继续检验原生上下文，但不能归因全部水稻差距、推广到全TRAIN／DEV或形成SOTA、泛化或创新主张。长条件附带真实侧翼，新增侧翼与共同注意力的作用仍须区分；下游权重保持原短块训练状态，尚未做表示适配。

下一步扩大到原1491个TRAIN有限视图，保留原已完成短块基线、不筛掉失败或unsupported参考，并在完整结果后决定是否实施匹配更新次数的上下文适配训练。现有三组目标函数训练与其冻结DEV方案不变，不按本八见证结果挑DEV checkpoint或修改评价规则。

扩大范围前已读取原training_plan的几何元数据：1491个唯一视图中1487个为49152 bp，另外4个染色体端点视图分别为36864（2个）、40206及45924 bp。后续实现必须按真实染色体范围裁切编码上下文，并保留这些原视图和原参考分母，不能直接套用八个内部见证的固定49152断言或删掉端点。该检查只读取原TRAIN计划元数据，没有新增DNA或评价读取。

权威结果位于远程outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-DIAGNOSTIC/summary.json及八份window_00.pt至window_07.pt，以及outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-READOUT/summary.json、windows.jsonl。CPU已实际加载全部八份势函数并完成32条件读取；没有以提交状态或文件存在替代数值验证。
