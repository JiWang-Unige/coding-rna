# 等权双模型 raw FASTA 调用器接入

2026-10-13。0.5／0.5角色能量组合已接入 raw FASTA 调用器，原生验证13556298的六个方向读出与独立平均能量解码／GFF完全一致，双支实测GPU峰值4.536227GiB。入口直接读取序列和两支末步权重，不读取旧GFF或缓存势；这是工程接入与数值口径验证，不是新训练、调权或独立性能验证。

## 模型与推理规则

输入为 raw FASTA 和两个明确 checkpoint 参数。E为末 head9216／native adapter4608，M为末 head9216／dense qv4608，各自沿用原生 checkpoint 验证和 revision 固定；接线验证使用已完成 fit13499817／13510324。两模型使用独立 backbone，不把一支的 adapter／qv 状态覆盖到另一支。实际单支编码峰值显存2.524／2.562GB，原生双支调用峰值另行测量，不能直接把估计当成实测。

沿用原 FASTA 规范化、24576窗口／12288 stride、nearest-containing-center ownership及负链规则。每条输入 contig 两个方向各重新编码 E/M，在相同 oriented DNA和ownership下逐角色先转FP64再取算术平均。原 whole-contig skyline、unit role、null0、cap128、codon carry、GT／GC-AG和min-intron5不变；不混合旧decoded chains或gene marginals，不设置 score gate。角色数组只存在内存，不保存新全长数组。输出 complete primary CDS GFF3和推理summary，不能宣称UTR、全部isoform、异向不重叠或校准正确率。

## 实际原生接线验证

固定复用此前后验13547015的两个原TRAIN49152bp视图：AT NC_003070.9，半开原坐标14376960–14426112；rice NC_089046.1，14450688–14499840。两视图来自两支既有原TRAIN输入，先直接比较完整原DNA，再裁切，不使用参考边界／catalog／prediction选择位置。另用已有14bp短序列`ATGCGTNAGTTTAA`，FASTA写成小写和R歧义字符形式，核对既有规范化规则和multiple-record入口，不作为生物学结果。

实际调用器一次读取三条FASTA、两个正式checkpoint，执行六个 oriented-contig 读出并写 GFF。随后从原两branch独立重编码相同三个输入，用独立FP64表达式`(E.double()+M.double())*0.5`及既有解码器重放全部链、scope总gain、格式化GFF gain、坐标、ordered phase、ownership和同向span非重叠。这回答接线是否保持原平均能量模型及输出口径；旧四组完整set／FP64／gauge／carry／cap资格13551105复用，不重复数学资格作业。

不重跑全部六染色体 native forward。原六染色体的缓存组合结果已完成，但本次原生接线只覆盖两个固定49152bp视图与短序列，不能称作新调用器已完成全原TRAIN端到端重编码或外部验证。CPU内存峰值和双支GPU峰值以实际作业为准；推理summary包括模型加载和每支fresh native编码时间，不拿 cached redecode时间冒充完整成本。

## 执行范围与命令

一个generanno／private-teodoro-gpu／RTX3090×1／2CPU／16GiB／10min／no-requeue原生验证作业，独立输出`outputs/M28-FULL-GRAMMAR-E-M-RAW-NATIVE`。模型来源是completed产物身份核验，不依赖已离开controller的旧job ID。新输入、GFF、metadata声明256MiB总输出上限；原DNA、GFF、checkpoint和冻结结果不改写。确定性接线失败至多同范围一次修复，独立新输出保留失败；不取消现有作业或自动延长。

调用器入口：

```bash
python scripts/experiments/M28-METHOD-RESTART/annotate_full_grammar_e_m_energy_pool.py \
  --fasta /absolute/path/to/authorized_input.fasta \
  --e-checkpoint outputs/M28-FULL-GRAMMAR-JOINT-LORA-FIT/step_009216.pt \
  --m-checkpoint outputs/M28-FULL-GRAMMAR-DENSE-QV-FIT/step_009216.pt \
  --output-dir /absolute/path/to/new_output
```

这只是接口用法，新的数据仍需在当前授权范围内，经Slurm执行；不授权DEV、sealed／held-out或新物种。对本验证，全部root／batch／extern COMPLETED0:0、STATUS／JOBID／started／summary、非空input及prediction GFF／summary，且六个原生重放全部通过才构成完成。若外部等待超过当前短验证，由原固定监控对话和同一heartbeat接续，不创建新对话／自动化。方法记录按现有协议格式写清工程完成、数值一致与科学性能的区别；没有发布或对外传输。

## 实际接线

两份Python AST和SBATCH语法通过，四个新叶文件精确同步回读一致后，提交`sbatch --parsable sbatch/M28-FULL-GRAMMAR-E-M-RAW-NATIVE.sbatch`，实际job13556298。提交前own expanded45，小于80，目标路径不存在，无10min维护冲突。基于既有单支峰值及可用RTX3090选择private；实际1GPU／2CPU／16GiB／10min／no-requeue，13秒排队后在gpu034运行，未改资源或训练参数。

## 完成证据与原生匹配

13556298于09:04:47–09:06:43UTC运行116秒，root／batch／extern全部COMPLETED0:0，MaxRSS4382624K。STATUS／JOBID／started／5979-byte summary、98354-byte输入FASTA、8100-byte预测GFF和19453-byte推理summary齐全。程序64.291635秒，实际调用器模型加载14.293428秒、包含fresh双支编码和解码的raw loop24.556087秒，共38.849515秒；余下程序时间包括源输入核验和独立重编码／重放，不据此宣称全基因组速度优势。

| 原生输入 | 方向 | 输出完整链数 | Scope gain | 独立平均能量与导出重放 |
| --- | --- | ---: | ---: | --- |
| AT固定49152bp视图 | + | 0 | 0 | 一致 |
| AT固定49152bp视图 | − | 0 | 0 | 一致 |
| Rice固定49152bp视图 | + | 1 | 216.64205093076453 | 一致 |
| Rice固定49152bp视图 | − | 2 | 1182.8977205420379 | 一致 |
| 14bp短序列 | + | 0 | 0 | 一致 |
| 14bp短序列 | − | 0 | 0 | 一致 |

全部六例的whole-oriented-contig gain、输出数、窗口和ownership一致，三条非空输出的全部CDS区间、ordered phase及8-significant-digit导出gain一致，导出后literal DNA／phase和same-strand span非重叠均通过。这里3是完整链输出数量，不是参考TP或新发现基因数。双支实测allocated4870736896bytes（4.536227GiB）、reserved5186256896bytes（4.830078GiB），两backbone全部梯度禁用，模型年龄末head9216／native4608保持不变。

主独立复算六例身份／scope gain／输出总数／实际资源数值、completed源／相同固定identity及所有无mutation／无参考／无DEV／无sealed字段，与两份实际summary一致。全部产物小于256MiB；主已直接消费，未加入watcher、不补发监控交接。旧缓存全TRAIN评价不被改写，新接口尚未重新编码全部六染色体，亦没有独立性能、校准、泛化或新颖性结论。接下来将其作为已核验可执行入口，继续结合R1后验／失配诊断和LR对照选择有科学依据的下一项TRAIN改进；冻结DEV family与封存数据边界不变。
