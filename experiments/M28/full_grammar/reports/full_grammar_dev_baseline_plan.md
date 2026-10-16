# 最终4608步模型的冻结DEV／基线评价准备

2026-10-09，Asia/Shanghai。状态：**固定DEV推理、三基线共同复评及只读局部几何诊断均已完成，原监控一次交接后暂停**。实际结果见[DEV比较结果](full_grammar_dev_baseline_result.md)：最终候选两物种完整链F1均低于三基线，工程完成不等于科学竞争优势。本轮人类原话为“授权，这种不再需要我批准。甚至我允许你并行执行多项实验来更快速完成goal”；据此执行此前准备的固定DEV推理、三条既有原生基线共同复评及常规实施、验证与结果处理。此记录转述当前对话授权，不以文件存在推定权限。最终4608步严格raw-FASTA原TRAIN六条推理及共同尺诊断已完成，见[full_grammar_raw_train_result.md](full_grammar_raw_train_result.md)。本轮不改训练参数、不增训、不放行sealed/held-out/Setaria，不依TRAIN或DEV数值选择中间权重、修改评价尺或公开未发表研究。

## 冻结对象与比较问题

问题是“最后4608步coarse-channel模型在既有同物种DEV染色体上的严格raw-FASTA完整CDS链表现，及相对三条固定原生工具的准确率／背景负担”，不是新独立测试、clean zero-shot、因果消融或最终论文成功。当前原TRAIN诊断也不能替代该比较；独立物种及更完整训练不在本阶段。

唯一候选为 `step_004608.pt`，来自固定seed0/原1536 draws/1491唯一窗/三遍/epsilon0.1首训。全部权重、冻结GENERanno revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345`、C0结构、原窗口网格、nearest-center base context、双strand whole-contig exact skyline、cap128、unit role、null0、FP64及已有tie policy原样冻结。不会对DEV选择step1536/3072、调epsilon、温度、tau、role权重、cap、grid或decoder。观测模型和训练标签不进入推理。无需重新训练三条基线。

## 输入与评价分母

后续仅既有两条DEV contig：拟南芥NC_003074.8/23459830bp、水稻NC_089041.1/29936421bp。所有方法接受相同完整原FASTA记录，不加softmask、RNA、蛋白、参考边界或其他证据，不按基因裁剪。M28原caller的大写/非ACGT转N策略原样声明，不改变长度或坐标；基线保持其既定原生预处理并记录差别。

共同primary reference使用既有M26 CDS-assessable longest-eligible-CDS/tie-transcript-ID尺，已冻结分母为拟南芥5437、水稻2291、合计7728 unique primary CDS chains。后续先完成全预测和GFF，然后只按授权范围重放这个尺。身份/坐标/分母或GFF phase不一致就报告并停在结果核验，不以删unsupported、partial标签改写、motif/cap/长链/overlap筛选使其匹配。完整CDS可评估性与gene/UTR partial分别处理，不把M28训练eligible当测试参考筛选条件。

预测统一从完整原生输出选每个predicted gene最长CDS、平分取transcript ID，`complete_only=False`；不按与reference的匹配、DNA合法性或模型可预测性删除基线预测。完整链key含seqid、strand及所有半开CDS边界；exact phase单列诊断，不偷换为只匹配一个边界、片段或CDS-base结果。其他assessable isoform命中另计，不加入primary TP。UTR／所有isoform／非编码gene能力本轮不评价，也不宣称已有覆盖。

## 固定基线与复用边界

复用2026-10-03六次原生运行13358414_0–5及共同评价13358432所记录的固定资产，不称其为未经重新查询的“今日最新版本”。具体证据在[native_baseline_preparation.md](native_baseline_preparation.md)、[native_baseline_result.md](native_baseline_result.md)；准备阶段未重新读取DEV产物，获授权后的共同评价将按下述边界核对并复用。

| 工具 | 固定源码／权重 | 原生运行操作点 |
| --- | --- | --- |
| ANNEVO | 42c920f70f52f96c7c06c665a639ec47d9304fff／v2.3.3，Magnoliopsida | batch32，8CPU，普通植物模式 |
| Tiberius | 871d4979568ae09e56de7eca5d7489046ee833ea／实际2.0.8隔离包，基础SIF实际2.0.7，angiosperms | batch4，seq_len200070，softmasking=false，clamsa=false |
| Helixer | d17bb496b8842542590644b2437798711c0125d5／0.3.7，land_plant_v0.3_a_0080.h5 | batch16，subsequence64152，既定原生HelixerPost默认与deterministic |

后续复用预测前须核对六组实际成功终态/STATUS/运行收据、输入seqid/长度/原FASTA直接相等及固定版本、GFF结构/坐标、无外部证据。这里使用真实单元格、FASTA直接比较和已有资产身份，不新建无用途hash。若证据不匹配，先报告具体缺口；不得静默换为历史缓存表、更新版本或为更好成绩重跑。科学比较优先复用合法固定输出，不因新模型出现自动重跑所有工具。

## 固定读数与解释

复用 `common_ruler.select_reference/score`；逐物种列完整链TP/预测数/参考数及P/R/F1、不分strand CDS-base union P/R/F1、非primary assessable isoform命中、重复链数、背景span/CDS bp覆盖、纯背景链/Mb。背景固定为所有annotated gene feature跨度任一strand之外，覆盖率不是gene-level分类FPR，也不等于已确认的新假阳性基因。

合并表由计数micro汇总，同时给两物种macro与worst-species，报告每物种相对每基线的差值（百分点）及背景负担，不以单个pooled胜出写成全面优越。单／多CDS段、超过24576bp原生窗口的跨度、same-strand primary overlap及非标准motif等参考分层保留，描述性分层不能被用作重选主分母。长链和多CDS段不是全染色体语法的自动unsupported条件；超cap128、模型语法实际不允许的motif/终端结构及同strand重叠与这些描述性长度分层区分记录。当前不新增科研PASS阈值，也不恢复已失败旧路线的gate；所有结果如实记录，负结果不自动触发参数调整。

耗时分开记录原生encoding、joint decoding、输入物化、模型加载、Slurm wall和GPU峰值。新M28计时与旧工具计时只有起止点、输入范围和硬件一致时才可作直接效率比较；否则标为历史运行成本，不能用缓存feature时长宣称raw速度优势。训练成本单列，不混入纯推理表。OOM/超时/未完成是工程／资源边界，不能当0准确率或静默排除。

## 已准备的执行接线

新增 `infer_full_grammar_fixed_dev.py`及 `evaluate_full_grammar_dev_baselines.py`，不修改已完成TRAIN的生产caller、权重或评价结果。DEV runner只支持上面两条固定contig，保留原FASTA的大写记录供直接比较；仍由原caller在模型输入前将非ACGT转N。后续评价将候选、原生基线和原始source FASTA记录逐base直接相等比较，参考尺仍读取原始序列，与既有原生比较的共同尺一致，不通过模型归一化改参考可评估性。

后续CPU评价先核对两份候选及六份固定原生运行：完成字段、最后checkpoint、模型身份、双strand完整grid、whole-contig评分与tie policy；固定方法／物种／版本／权重收据；原生实际命令与batch、seq_len及softmask操作点；原生完成标记和GFF范围／坐标；Slurm实际数字JobID对应root/batch/extern COMPLETED/0:0；然后直接比较三路输入。全部通过才解析reference。任何缺口会停止评价而非删预测、换权重、重跑基线或读取其他数据。

共同评价同时输出四方法的逐物种计数、primary exact-chain／CDS-base／背景指标、phase诊断、其他isoform及重复链数，pooled按计数合并并列macro／worst-species。单／多CDS段及span>24576bp仅作结构分层，不修改参考。实际运行收据随表保留，性能差值为每物种M28相对三基线的完整链F1百分点，不自动生成科学PASS或速度优势claim。

两个SBATCH为 `M28-FULL-GRAMMAR-FIXED-DEV-RAW.sbatch`（2元素，最多并行2，每元素1×3090／2CPU／16GiB／2小时）与 `M28-FULL-GRAMMAR-FIXED-DEV-COMPARISON.sbatch`（0GPU／2CPU／8GiB／20分钟，整个候选数组afterok）。授权后已按这两个固定操作点提交，见下方实际记录。提交前核对队列、维护、BeeGFS可用空间及正式输出目录不存在；private-gpu空闲3090足够并行两元素。资源选择基于原TRAIN实际24–38分钟／contig及既有共同尺资源，不是完成时间保证。两种作业沿用原成功generanno环境，不安装依赖，不调用基线原生推理。

本轮正式产物为 `outputs/M28-FULL-GRAMMAR-FIXED-DEV-RAW/<seqid>/`（JOBID、STATUS、started.json、原input.fasta、summary及prediction/{summary.json,full_grammar.gff3}）和 `outputs/M28-FULL-GRAMMAR-FIXED-DEV-COMPARISON/`（JOBID、STATUS、summary.json）。这些目录在准备检查时未创建，由本轮获授权的作业首次建立，不覆盖历史TRAIN或原生基线结果。若实际身份或输入不一致，先保留事实、明确问题，不静默放宽合同；常规实现错误可在已授权冻结范围内作最小修复并保留失败产物，但不得改变科学参数、重选候选、自动扩展资源或放行其他数据。

`check_full_grammar_dev_preparation.py`只使用人工短序列、人工receipt和人工调度行，检查新接线会影响科学解释的错配：最后权重／标签边界／coverage、固定基线物种与实际版本、失败或缺失Slurm步骤、未完成时不得进入reference选择、strand-sensitive结构分层，以及driver保留原IUPAC记录再调用不变caller。它以mock替代研究数据读取和模型运行；不重复旧caller的训练或完整语法检查。准备阶段仅提交此人工检查，未打开历史DEV receipt／GFF／DNA或label；其终态成功只能证明人工接线检查，不是实际DEV就绪或模型表现。

人工检查13449873已完成：UTC 2026-10-09 02:56:43–02:57:12，root/batch/extern COMPLETED/0:0、实际wall29秒、MaxRSS861492KiB；五类人工接线检查通过。唯一检查root `outputs/M28-FULL-GRAMMAR-DEV-PREPARATION-CHECKS/`有JOBID、STATUS和summary，明确research DNA/labels、历史DEV receipts/predictions、checkpoint和native forward均未读取／执行，optimizer steps=0。三份新增Python AST及三个SBATCH语法通过；Mac只做语法与文本操作，人工运行经Slurm。这些是授权执行前的准备证据，未经改动不重复检查。

## 本轮实际执行与等待交接

人类明确授权后已提交下列固定作业。实际Slurm显示两条推理于2026-10-09 UTC03:12:18（北京时间11:12:18）开始运行；不将后续记录更新时间当作计算起点。

- DEV推理数组 **13449922**，任务0为拟南芥NC_003074.8，任务1为水稻NC_089041.1；首次调度快照两任务均RUNNING于gpu034。
- 共同评价 **13449926**，显式 `afterok:13449922`，首次快照PENDING/Dependency；只在两条候选完整成功后核对六份历史基线输出、输入与实际Slurm终态，然后读取固定DEV参考并评分。
- 等待仅复用监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3/local` 和心跳 `coding-rna-r5-light-monitor`，15分钟、gpt-5.6-luna/low；正常RUNNING、PENDING及依赖等待不回报主对话，整条流程终态或需要处理的失败一次交接后暂停本轮心跳。

实际完成：两条推理分别UTC03:43:01／03:49:13结束，共同评价03:49:14–03:49:56完成，全部root/batch/extern COMPLETED/0:0及正式产物齐全。输入、权重、原生版本／操作点、参考分母及数值已由实际CPU评价和主对话核验。后续只读局部几何job13450082同样成功，使用原TRAIN已定义的事件指标，不改变本冻结主评价或候选。完整科学比较见结果文件；研究goal尚未完成。

UTC 03:18后，原心跳更新成功；本机实际保存字段当时为ACTIVE、heartbeat、原目标对话、15分钟与原通知策略。原监控对话的启动轮次 `01a11eab-4028-7260-bed5-1945bb914c34` 实际turn_context为gpt-5.6-luna/low，未改变主对话设置。首次独立检查于UTC 03:20:30完成并写回同一journal；本轮定时检查也已实际观察：轮次 `01a11eb9-0d5c-7fb2-b4ac-bc8a9c78ae51` 于UTC03:33:26–03:33:55运行，heartbeat wrapper对应同一automation与当前DEV cohort，实际model/effort为gpt-5.6-luna/low，并更新同一journal；正常运行时未向主对话发消息。完整终态核验与一次回传现已观察：原监控轮次 `01a11ec9-6939-74d0-a47d-a9c24920bca0` 于UTC03:51:18–03:53:08关闭本cohort，随后同一心跳PAUSED并保留。监控沿用旧03:48:26时间标签的记录错误已按实际turn metadata修正；运行成本仅使用Slurm Start/End，不受此标签影响。

## 资格与下一步

这两条DEV已经被多轮研究使用，当前checkpoint在此次比较前固定不能抹去历史开发依赖。GENERanno为发布CDS注释checkpoint的冻结backbone，其上游物种/annotation/同源暴露仍unknown；ANNEVO植物权重作者声明排除拟南芥/水稻但未独立审计；Helixer当前权重明确包含两物种；Tiberius包含水稻。详见[pretrained_exposure.md](pretrained_exposure.md)。raw FASTA输入与各方法未见物种是不同条件，不能统一称clean zero-shot。

原TRAIN共同尺诊断已完成：拟南芥完整链F1=50.4871%，水稻=28.2780%，pooled=45.5901%；CDS-base F1分别93.0848%、80.8511%、90.3601%。这是有限原draw首训后的训练染色体描述性读数，不是DEV、独立泛化或基线比较；不能以loss下降或base指标替代完整gene准确率。该结果不触发任何自动改参／增训，也不改变本固定候选。

本次固定DEV推理、三基线共同复评、必要结果核验和局部几何描述性分析已完成，不读封存数据、不公开论文或模型。沿用本固定输入、最后权重、原生基线和共同尺，预测先行／标签后置，并仅复用原监控对话与同一心跳；本cohort已一次回传后暂停，不再检查旧作业。独立面板与论文claim所需证据另按真实授权和科学隔离要求推进。原论文goal继续active，本阶段完成不意味着整个goal完成，也不由负结果自动增训、改参或重新定义成功。
