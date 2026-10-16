# M28原生基线结果与当前模型差距

2026-10-03。六次固定版本原生 raw-FASTA 注释及共同评价已完成。ANNEVO 在这两条 DEV contig 的完整 CDS 链 F1 最高；A3 尚未达到该水平。该表是开发集上的当前版本比较，不是独立测试、干净 zero-shot 证明或论文成功。

## 完成证据与数据范围

原数组13358414_0–5（ANNEVO拟南芥/水稻、Tiberius拟南芥/水稻、Helixer拟南芥/水稻）及评价13358432的Slurm root/batch/extern均COMPLETED/0:0。六元素2026-10-03 10:40:03 UTC启动，最终元素10:46:03结束；评价10:46:04–10:46:24 UTC执行。主对话核对六组STATUS/summary、方法/物种、真实GPU/原生命令完成标记与GFF结构，并由TP/预测数/参考数独立复算precision/recall/F1、背景比例及pooled/macro汇总。

输入仅为NC_003074.8（拟南芥23459830 bp）与NC_089041.1（水稻29936421 bp），原DNA大写化但坐标不变。六次推理均reference_used=false、labels_or_feature_cache_used=false、optimizer_steps=0、test_setaria_used=false；实际每元素只见一张RTX3090。参考只在共同评价阶段打开。两物种固定primary分母5437/2291，共7728条，未删unsupported、长链、重叠或难例。

固定运行资产及隔离Tiberius修复见[native_baseline_preparation.md](native_baseline_preparation.md)；预训练暴露见[pretrained_exposure.md](pretrained_exposure.md)。Tiberius实际2.0.8源码871d497在不可变的实际2.0.7依赖SIF内隔离执行，不能称为发布方原生2.0.8镜像。Helixer a_0080训练列表含两物种；ANNEVO实际植物权重作者声明排除两者但未独立审计同源/验证暴露。GENERanno当前b048来源已有RefSeq CDS监督及两个目标assembly的训练快照归属证据；仍未知的是具体优化窗口、验证划分与同源隔离（来源表述2026-10-13更新）。整表不能统一称为未见物种泛化，也不据快照归属认定具体DEV标签泄漏。

## 同尺完整链与背景

完整链TP要求同strand、相同有序CDS坐标链命中固定primary参考；phase单独诊断，不属于坐标exact-key。precision=TP/预测唯一链，recall=TP/参考primary链，F1=2TP/(预测唯一链+参考primary链)。其他eligible isoform命中不计primary TP；非精确预测也不等于生物学假基因。

| 物种 | 方法 | TP / 预测 / 参考 | Precision % | Recall % | 完整链F1 % | 背景预测span bp | 背景span覆盖 % |
|---|---|---|---:|---:|---:|---:|---:|
| 拟南芥 | ANNEVO | 3941 / 4898 / 5437 | 80.4614 | 72.4848 | 76.2651 | 34221 | 0.4182 |
| 拟南芥 | Tiberius | 3865 / 4984 / 5437 | 77.5482 | 71.0870 | 74.1771 | 29321 | 0.3583 |
| 拟南芥 | Helixer | 3673 / 5453 / 5437 | 67.3574 | 67.5556 | 67.4564 | 102135 | 1.2481 |
| 水稻 | ANNEVO | 1703 / 2285 / 2291 | 74.5295 | 74.3344 | 74.4318 | 289043 | 1.4618 |
| 水稻 | Tiberius | 1576 / 2443 / 2291 | 64.5108 | 68.7909 | 66.5822 | 172110 | 0.8704 |
| 水稻 | Helixer | 1439 / 3086 / 2291 | 46.6299 | 62.8110 | 53.5243 | 874226 | 4.4213 |

背景为所有注释gene feature跨度（任一strand）之外的bp；拟南芥8183405 bp、水稻19773036 bp，合并27956441 bp。summary中的span_FPR实际是背景bp被预测链span覆盖的比例，不是逐基因分类FPR。按该口径Tiberius背景span覆盖最低，不能仅以ANNEVO链F1较高推定所有维度都优越。

| 方法 | 合并TP / 预测 / 参考 | Micro F1 % | 两物种Macro F1 % | CDS-base F1 % | 合并背景span bp | 背景span覆盖 % |
|---|---|---:|---:|---:|---:|---:|
| ANNEVO | 5644 / 7183 / 7728 | 75.7025 | 75.3485 | 93.7935 | 323264 | 1.1563 |
| Tiberius | 5441 / 7427 / 7728 | 71.8047 | 70.3797 | 94.0981 | 201431 | 0.7205 |
| Helixer | 5112 / 8539 / 7728 | 62.8512 | 60.4903 | 92.0312 | 976361 | 3.4924 |

Micro先合并计数再算F1；Macro为两个物种F1的简单均值。CDS-base按不分strand的CDS union计算，参考union为拟南芥6601350 bp、水稻3063207 bp；它不衡量完整基因结构。Tiberius合并CDS-base F1略高于ANNEVO，但ANNEVO完整链F1更高。

六组重复CDS链记录均为0，各组所有CDS-bearing预测记录均phase-consistent。非互斥unsupported strata仍保留：拟南芥split_terminal_codon9 / same_strand_primary_overlap12 / nonstandard_splice32 / longer_than_window1；水稻5 / 9 / 38 / 2。不以这些诊断重定义参考分母。

## 与A3及候选覆盖上限的关系

| 模型/方法 | 拟南芥完整链F1 % | 水稻完整链F1 % | 合并Micro F1 % |
|---|---:|---:|---:|
| A3（已有缓存特征DEV实验） | 51.7567 | 53.5236 | 52.3200 |
| ANNEVO（本轮原生） | 76.2651 | 74.4318 | 75.7025 |
| Tiberius（本轮原生） | 74.1771 | 66.5822 | 71.8047 |
| Helixer（本轮原生） | 67.4564 | 53.5243 | 62.8512 |

ANNEVO相对A3的完整链F1差距为拟南芥24.5084个百分点、水稻20.9082个百分点、合并23.3825个百分点。A3与Helixer在水稻只差0.0006696个百分点，近乎精确持平，不能表述为有意义的性能优势。A3原科学gate FAIL保持；本轮原生表不覆盖其冻结记录，也不把缓存特征耗时与原生整条FASTA耗时直接比较。

A3原自由池精确primary命中为拟南芥3022、水稻1421，合并4443；owner后2590/1308，合并3898；最终选中2519/1219，合并3738。仅用于诊断的乐观上限为2H/(H+R)：假设完美precision、全部精确TP都保留，忽略非重叠与unsupported冲突，故不是可实现模型结果或新gate。

| 范围 | 自由池H | Owner池H | 自由池乐观F1 % | Owner池乐观F1 % | Precision=1时严格超过当前ANNEVO所需最少TP |
|---|---:|---:|---:|---:|---:|
| 拟南芥 | 3022 | 2590 | 71.4505 | 64.5322 | 3352 |
| 水稻 | 1421 | 1308 | 76.5625 | 72.6869 | 1359 |
| 合并 | 4443 | 3898 | 73.0096 | 67.0566 | 4707 |

必要最少TP计算为floor(TP_ANNEVO×R/(Pred_ANNEVO+R−TP_ANNEVO))+1，并非保证可达的目标。A3自由池合并4443条，即使理想选择也只有73.0096%，低于当前ANNEVO75.7025%；拟南芥自由池上限也低于ANNEVO。A3原owner池两物种上限都低于ANNEVO。因此只重排序原greedy池不足以支撑全面超越ANNEVO的主张；必须获得更多可保留的正确候选。合并自由池至少还缺264个TP才能跨过这个理想下界。

已授权的FULL Q/C0在全TRAIN重训生成器；A6冻结A3生成器，以greedy+beam4联合自由池训练非加性质量头，均仍待完整结果。A6实际DEV自由池/owner池覆盖尚未得到，不能用8窗TRAIN预检中的beam增益冒充DEV性能。当前诊断支持检查这两条现有路线，不修改冻结阈值、评价分母或sealed限制，不新增实验。

进一步的[逐参考缺口诊断](native_reference_diagnostic.md)已完成：ANNEVO正确而A3未选出的2143条中，1573条未生成、486条在owner丢失，两类合占96.0803%；缺口集中在多CDS片段链。111条A3独有命中全部是既有primary参考，不是新基因发现。该诊断只分析已有DEV产物，不修改运行中的FULL/A6或冻结规则。

## 运行收据与剩余工作

| 方法 | 拟南芥入口耗时秒 | 水稻入口耗时秒 |
|---|---:|---:|
| ANNEVO | 89.168 | 91.893 |
| Tiberius | 218.260 | 229.527 |
| Helixer | 264.520 | 359.634 |

每元素1RTX3090/8CPU/32GiB/2小时上限；共同评价0GPU/2CPU/8GiB/20分钟上限，实际20秒。入口计时包含GPU probe、进程启动、I/O及输出核对，不含排队；Slurm六元素elapsed分别90/92/219/230/265/360秒。不同模型、CPU/memory/cache管线与并发条件尚未控制，这些数字仅是操作收据，不是公平速度benchmark。

2026-10-03 11:00 UTC固定监控已自动交接本原生分支；主对话已消费并核对指标。该次状态中FULL fits13358083/13358084与A6 fit13358977仍RUNNING，其各自infer/eval为有效Dependency。复用同一监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor，只接续FULL五作业和A6三作业；已完成原生分支退出监控范围。FULL/A6成功后主对话按原协议验证、解释；未完成时不宣布winner、独立泛化或论文目标达到。

## 可追溯实现与原始产物

远程规范仓库：/home/users/j/jwang/coding-rna。本记录及preparation/exposure记录同步到本机项目；不复制原始大产物，不提交Git或对外发布。

- scripts/experiments/M28-METHOD-RESTART/evaluate_native_baselines.py：原生完成边界、方法/物种范围及共同评价汇总。
- scripts/experiments/M28-METHOD-RESTART/common_ruler.py：固定primary选择、CDS链、CDS-base union及背景定义。
- outputs/M28-NATIVE-BASELINES/{ANNEVO,Tiberius,Helixer}/{arabidopsis_thaliana,oryza_sativa}/{summary.json,prediction.gff3,command.log,STATUS,JOBID}：六次真实原生运行。
- outputs/M28-NATIVE-BASELINE-EVAL/{summary.json,STATUS,JOBID}：本表原始计数与汇总。
- outputs/M28-A3-EVALUATION/summary.json：A3既有分数及自由池/owner/selected命中诊断。
- outputs/M28-PILOT-R5/logs/M28NATIVE_13358414_<0-5>.out/err与M28NATEVAL_13358432.out/err：Slurm日志。
