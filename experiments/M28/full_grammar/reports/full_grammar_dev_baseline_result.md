# 最终4608步模型的固定DEV与三基线比较

2026-10-09。固定两条DEV的raw-FASTA推理及四方法共同评价已完成，输入、参考分母和数值可核验。最终4608步M28完整CDS链F1为拟南芥50.8442%、水稻30.0960%，两物种均低于ANNEVO、Tiberius和Helixer；pooled为44.6085%，不支持当前候选具有准确率竞争优势。背景覆盖在部分比较中较低，但不能抵消完整链P/R缺口或据此宣称更广适用范围。本轮没有改参、增训、选中间checkpoint、访问sealed/held-out/Setaria或公开模型／论文，原研究goal未完成。

## 执行与评价条件

按[冻结DEV方案](full_grammar_dev_baseline_plan.md)及人类后续明确授权，唯一候选为 `step_004608.pt`，来自原1536 draws／1491唯一窗／seed0／三遍的有限TRAIN discovery fit。冻结GENERanno revision、原C0、所有15角色nearest-center ownership、完整raw网格、双strand whole-contig exact skyline、cap128、unit role、null0、FP64及同分规则不变；推理没有reference、标签／feature cache、观测模型或optimizer。不是根据DEV数值选择训练步数或运行阈值。

拟南芥NC_003074.8（23459830bp）和水稻NC_089041.1（29936421bp）完整预测先于reference选择。13449922_0–1及afterok共同评价13449926的root/batch/extern均COMPLETED/0:0，STATUS和规定产物齐全；数组显示ID与实际数字JobIDRaw分别为13449922_0→13449923、13449922_1→13449922。评价核对最终模型身份、全部raw网格及两链joint输出、六次历史原生运行的成功终态／固定版本和权重／原生命令及操作点／GFF完成标记，再对candidate、baseline preparation及original source三路FASTA记录逐base直接比较，全部相等后才解析参考。

复用2026-10-03六份固定原生输出，不重新执行基线或宣称其为今日最新版本：ANNEVO 2.3.3／Magnoliopsida，Tiberius实际2.0.8隔离包／基础SIF实际2.0.7／angiosperms，Helixer 0.3.7／land_plant_v0.3_a_0080.h5；源码和完整资产收据见[原生基线记录](../../../../benchmarks/native_baselines/reports/native_baseline_result.md)。同原始FASTA、大写与原生固定预处理，不加mask、RNA、蛋白或参考边界；M28仅在模型输入时将非ACGT转N，不改变长度或参考尺。

共同参考固定为M26 CDS-assessable最长eligible CDS、平分取transcript ID：拟南芥5437、水稻2291、合计7728 unique primary链。所有方法从完整原生预测选每gene最长CDS、ID平分，`complete_only=False`；不按reference匹配、motif、cap、长链、DNA合法性或overlap删除基线预测。exact key包含seqid、strand及全部半开CDS边界，phase另核验。M28 GFF数量、实际DNA／ORF／phase／cap和同strand非重叠符合生产约束；四方法所有parsed prediction records均phase-chain consistent、重复链数均0。其他assessable isoform命中不加入primary TP，不评价UTR、全部isoform或非编码gene。

## 完整链与CDS base结果

以下P/R/F1均为百分比；CDS-base为不分strand的union指标，不替代完整链准确率。

| 范围 | 方法 | TP／预测／参考 | 完整链P | 完整链R | 完整链F1 | CDS-base F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 拟南芥 | M28 final4608 | 2650／4987／5437 | 53.1382 | 48.7401 | 50.8442 | 92.7856 |
| 拟南芥 | ANNEVO | 3941／4898／5437 | 80.4614 | 72.4848 | 76.2651 | 94.8089 |
| 拟南芥 | Tiberius | 3865／4984／5437 | 77.5482 | 71.0870 | 74.1771 | 95.6153 |
| 拟南芥 | Helixer | 3673／5453／5437 | 67.3574 | 67.5556 | 67.4564 | 95.2979 |
| 水稻 | M28 final4608 | 674／2188／2291 | 30.8044 | 29.4195 | 30.0960 | 83.8545 |
| 水稻 | ANNEVO | 1703／2285／2291 | 74.5295 | 74.3344 | 74.4318 | 91.6185 |
| 水稻 | Tiberius | 1576／2443／2291 | 64.5108 | 68.7909 | 66.5822 | 90.8470 |
| 水稻 | Helixer | 1439／3086／2291 | 46.6299 | 62.8110 | 53.5243 | 85.5152 |

| 方法 | pooled TP／预测／参考 | pooled完整链F1 % | species macro % | worst species % | pooled CDS-base F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| M28 final4608 | 3324／7175／7728 | 44.6085 | 40.4701 | 30.0960 | 90.0898 |
| ANNEVO | 5644／7183／7728 | 75.7025 | 75.3485 | 74.4318 | 93.7935 |
| Tiberius | 5441／7427／7728 | 71.8047 | 70.3797 | 66.5822 | 94.0981 |
| Helixer | 5112／8539／7728 | 62.8512 | 60.4903 | 53.5243 | 92.0312 |

M28相对ANNEVO／Tiberius／Helixer的完整链F1差值：拟南芥−25.4209／−23.3329／−16.6122个百分点，水稻−44.3358／−36.4862／−23.4283个百分点，pooled−31.0940／−27.1962／−18.2427个百分点。M28 pooled链precision=46.3275%、recall=43.0124%；它的预测链数7175与ANNEVO7183接近，并非只通过少输出而获得更高完整链precision。较高CDS-base F1不能改写整链比较的负结果。

非primary assessable isoform命中M28／ANNEVO／Tiberius／Helixer为242／495／517／424，继续单列，不补进primary TP。主对话从实际TP／预测／参考计数独立重算逐物种和pooled P/R/F1、结构分层、背景比例、micro汇总、macro／worst及差值，均与实际summary一致；没有更换评价口径或重跑未变评价。

## 背景负担与结构分层

背景固定为全部annotated gene feature跨度任一strand之外，拟南芥8183405bp、水稻19773036bp、合计27956441bp。span覆盖比例不是gene-level分类FPR，也不证明背景预测全为生物学假阳性。

| 方法 | 拟南芥背景span覆盖 % | 水稻背景span覆盖 % | pooled背景span覆盖 % | pooled纯背景链 | pooled纯背景链／Mb |
| --- | ---: | ---: | ---: | ---: | ---: |
| M28 final4608 | 0.4154 | 0.7330 | 0.6400 | 64 | 2.289276 |
| ANNEVO | 0.4182 | 1.4618 | 1.1563 | 146 | 5.222410 |
| Tiberius | 0.3583 | 0.8704 | 0.7205 | 165 | 5.902039 |
| Helixer | 1.2481 | 4.4213 | 3.4924 | 563 | 20.138472 |

M28在水稻及pooled上的背景覆盖较三基线低，但拟南芥span覆盖高于Tiberius；并不存在所有条件下统一胜出。M28 pooled背景span／CDS覆盖为178921／43936bp，纯背景链64；它在真实primary链P/R上的明显损失与背景负担是不同读数，不能选择其中一个作为综合性能优势。

| 方法 | 拟南芥单CDS段F1 % | 拟南芥多CDS段F1 % | 水稻单CDS段F1 % | 水稻多CDS段F1 % |
| --- | ---: | ---: | ---: | ---: |
| M28 final4608 | 73.2209 | 44.0385 | 64.7236 | 20.2067 |
| ANNEVO | 78.9212 | 75.4574 | 77.7966 | 73.2627 |
| Tiberius | 75.7674 | 73.6816 | 71.3608 | 64.8415 |
| Helixer | 77.4725 | 64.3970 | 65.6587 | 49.2979 |

M28的多CDS段缺口在DEV也明显，尤其水稻；这些固定几何分层不是筛选后的新主分母。参考中nonstandard splice为拟南芥32／水稻38，split terminal codon为9／5，同strand primary overlap为12／9，全部保留。超24576bp span参考仅1／2条，M28与Tiberius均无exact命中，ANNEVO与Helixer各命中水稻1条；M28没有此跨度的所选预测。whole-contig语法支持跨窗不等于当前权重已经获得长基因性能，这3条也不能支撑广泛适用范围结论。

## 局部几何补充诊断

CPU-only结果处理复用[原TRAIN诊断](full_grammar_raw_train_result.md)中已定义的strand敏感unique CDS start／stop位置、CDS段两边界及intron-pair两边界指标，作用于全部四方法冻结输出。没有bp容差、不混入phase、不改primary链分母，不把片段事件得分称为gene准确率，也不由该描述性分析选择模型或参数。`profile_full_grammar_dev_local_geometry.py`实际job13450082在新root `outputs/M28-FULL-GRAMMAR-FIXED-DEV-LOCAL-GEOMETRY/`成功；UTC03:57:24–03:57:53、root/batch/extern COMPLETED/0:0、wall29秒、内部12.800秒、MaxRSS1262252KiB，实际0GPU／2CPU4GiB。每组parse的完整链计数与原13449926相同，参考5437／2291不变；主对话从事件计数重算PRF和micro相加也一致。只读取已获授权的DEV输入、reference与保存GFF，未读checkpoint／feature，不重跑原生模型或主评价。

| 范围 | 方法 | CDS start F1 % | CDS stop F1 % | CDS段 F1 % | intron-pair F1 % |
| --- | --- | ---: | ---: | ---: | ---: |
| 拟南芥 | M28 final4608 | 78.7297 | 78.5721 | 77.6309 | 80.8187 |
| 拟南芥 | ANNEVO | 85.4074 | 87.2629 | 91.1308 | 93.4937 |
| 拟南芥 | Tiberius | 84.3570 | 87.0993 | 90.6066 | 93.2731 |
| 拟南芥 | Helixer | 83.4604 | 83.8247 | 86.9628 | 89.4197 |
| 水稻 | M28 final4608 | 68.2001 | 59.9375 | 60.1438 | 62.6517 |
| 水稻 | ANNEVO | 82.4481 | 85.6393 | 89.5724 | 91.9993 |
| 水稻 | Tiberius | 77.3294 | 80.5409 | 87.2528 | 90.2965 |
| 水稻 | Helixer | 70.6845 | 71.6146 | 76.9550 | 79.0596 |

pooled的M28／ANNEVO／Tiberius／Helixer intron-pair F1为75.4385%／93.0567%／92.4089%／86.1935%；CDS段F1为72.4372%／90.6707%／89.6181%／83.8080%。所有事件带contig和strand，共享位置／段去重；例如CDS start参考事件为5436／2290，而不是gene参考5437／2291，gene主分母没有因此缩小。CDS start／stop不是TSS／转录终止位点或UTR边界。

M28在两个物种的四类局部事件上也均低于三基线。因此当前差距不只体现为“已正确局部片段尚未组装成整链”；局部结构准确率本身也未达到这些工具的水平。TRAIN观察中内部结构差异占多数与此相容，但两项描述性分析不能识别encoder、head、监督覆盖、训练规模或全局解码哪个是原因，不能把“局部势已足够，只需换组装”作为已证实前提。

## 运行与监控证据

两条推理均UTC03:12:18开始，拟南芥03:43:01结束（Slurm wall30:43），水稻03:49:13结束（36:55）；共同评价03:49:14–03:49:56（wall42秒、内部26.657秒）。原推理→共同评价端到端wall37分38秒，不包含后续监控交接等待。推理每条1×3090／2CPU16GiB／2小时上限，CPU评价0GPU／2CPU8GiB；没有超时、重提或资源扩大。

M28内部raw主循环为1785.656／2158.582秒，原生编码合计2829.012秒、joint解码合计1113.503秒；Slurm单GPU占用合计4058秒，不等于GPU活跃利用率。最大batch RSS7349780KiB；PyTorch GPU allocated／reserved峰值约2.25／2.38GiB，不代表全进程VRAM。raw结果artifact bytes before summary合计9906014，不含物化FASTA等全部产物。本轮不将不同计时边界的历史基线收据用于速度或资源优越claim。

原监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3/local`完成本cohort结构白名单／非空产物核验并向原主一次交接，原心跳 `coding-rna-r5-light-monitor`实际PAUSED／heartbeat并保留，无新监控对话或重复自动化。本cohort原对话中的定时轮次已实际观察；终态交接轮次为UTC03:51:18–03:53:08。监控曾把上一轮03:48:26标签错误沿用到终态摘要，早于最后作业结束；主已在操作journal按实际turn metadata修正时间窗口，计算耗时只采用Slurm Start/End，科学指标不受该记录错误影响。监控没有读取或解释指标，本报告的科学比较由主完成。

## 结论资格与后续边界

这两条DEV已被多轮开发使用，不是独立测试；固定本次checkpoint不能抹去历史开发依赖。冻结GENERanno上游species／annotation／同源暴露unknown；ANNEVO植物权重作者声明排除拟南芥／水稻但未经独立审计，Helixer权重包含两物种，Tiberius包含水稻，见[预训练暴露记录](../../../../paper/supplement/pretrained_exposure.md)。严格raw FASTA输入与clean zero-shot是不同条件。当前比较可以回答这个固定开发面板下的真实性能差异，不能回答无偏跨物种泛化或因果机制优势。

工程接入及固定DEV比较完成，不意味着研究或论文成功。当前候选尚未显示相对三原生工具的准确率竞争力，也没有更广适用范围的实证；这是有限TRAIN discovery fit的负面开发结果，不是对所有GLM／结构化基因模型的否定。继续原模型与发表goal需要解释局部结构和整链缺口、准备与假设相对应的后续机制对照，并取得充分训练和独立面板证据。现阶段仍不改训练参数、不增训、不重开封存数据或失败旧路线、不修改冻结指标和阈值、不公开未发表研究。

原始结果位于Baobab `outputs/M28-FULL-GRAMMAR-FIXED-DEV-RAW/<seqid>/` 及 `outputs/M28-FULL-GRAMMAR-FIXED-DEV-COMPARISON/summary.json`；项目 `/home/users/j/jwang/coding-rna`解析到 `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna`，resolved checkpoint路径指向同一最终权重。脚本／SBATCH及所有冻结和失败结果保留，不自动Git提交或发布。
