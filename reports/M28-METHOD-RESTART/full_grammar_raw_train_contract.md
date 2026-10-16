# 最终4608步 raw-FASTA 接入与原TRAIN诊断

2026-10-09，Asia/Shanghai。用户本轮明确恢复原goal，并授权先完成严格raw-FASTA接入、原TRAIN端到端诊断，再准备冻结DEV与基线评价方案。原论文及ANNEVO/Tiberius/Helixer比较目标不缩小；平台goal已经active，不新建、改写或以工程完成标记整个goal完成。现阶段不改训练参数、不增训、不放行sealed/held-out/Setaria数据，也不执行DEV。

## 固定推理

仅用首训最终 `outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt`，严格核对步数、三遍、观测模型身份和冻结encoder revision；不选择中间权重。C0结构、原24576/12288网格、四块6144bp原生GENERanno编码、BF16物化后FP32 head、十五角色原logit减background、unit role、null=0、cap128均不变。不使用标签、epsilon观测模型、特征缓存或参考候选。

生产入口改为已有且验证过的exact exon-budget skyline：同stage/实际carry、较低exon预算且分数不差的max-only精确dominance，FP64逐base积累，8192bp检查点重算trace。它不恢复已停止的dense、A43或A44链。相同最优分数的链可按skyline已有tie policy选择，不声称与dense的平分路径逐链一致。

每个base的全部15角色来自最近且包含该base的原窗口中心；平分取较小genomic start，两strand相同，负链只映射行坐标、不换角色列。每条oriented全染色体一次联合complete nonoverlapping gene-set解码；不是窗内闭合再owner过滤。输出primary complete CDS链及GFF3 gene/mRNA/CDS，不主张UTR、全部isoform或跨strand不重叠。

## 唯一输入范围

| species | 原TRAIN seqid | bp |
| --- | --- | ---: |
| arabidopsis_thaliana | NC_003070.9 | 30427671 |
| arabidopsis_thaliana | NC_003071.7 | 19698289 |
| arabidopsis_thaliana | NC_003076.8 | 26975502 |
| oryza_sativa | NC_089043.1 | 27474823 |
| oryza_sativa | NC_089044.1 | 23931887 |
| oryza_sativa | NC_089046.1 | 28271460 |

Slurm任务按固定seqid/长度在保留DNA之前过滤species FASTA，仅物化当前一条TRAIN序列；丢弃其他记录，不解码、不保留非TRAIN DNA。与既有raw caller一致，大写化且非ACGT转N，不改变坐标/长度。推理不打开reference GFF、split_val或catalog。全染色体输出覆盖原TRAIN范围，而不仅是被原1536 draws曝光的1491个窗；这不增加训练曝光。

## 验证、运行和产物

入口检查 `13443304` root/batch/extern均COMPLETED/0:0；严格加载最后权重，实际生产绑定下两项既有测试/六次skyline调用通过：跨原生窗的30000bp完整基因双链与实际codon carry/含N内含子/负链GFF phase。该检查无原生GLM、研究DNA或标签读取，不是准确率结果。旧prospective检查与旧dense结果不覆盖。

诊断接入检查 `13443443` root/batch/extern均COMPLETED/0:0，STATUS/summary完成。以新诊断的实际imports、GFF writer/parser、共同尺及aggregate做两组合成检查：一组验证坐标/phase往返、strand敏感的完整链、其他isoform另计、背景按gene跨度而非CDS；另一组验证跨染色体micro先合并计数而非平均F1。人工fixture的pooled完整链TP/预测/参考为2/5/3，F1=0.5（不是两条fixture F1的均值），CDS-base计数为24/36/27；这些仅是可手算测试值，不是任何模型的TRAIN结果。检查不读研究DNA/标签、checkpoint或DEV，不跑原生GLM；生产推理和评价代码均无需修改。脚本 `check_full_grammar_train_diagnostic.py`、SBATCH与唯一产物 `outputs/M28-FULL-GRAMMAR-TRAIN-DIAGNOSTIC-CHECKS/summary.json`保留供复现，不重复运行未变检查。

正式数组六元素，各private-teodoro-gpu、1×3090至少24GiB、2CPU/16GiB、12小时、no-requeue，最多同时6条。成功首训相同generanno环境/本地固定encoder缓存，不安装依赖。此前30Mb人工skyline537秒、约5.3GiB RSS仅用于资源准备，不能保证真实神经分数frontier或耗时。原生全染色体编码预计小时级；12小时是每条固定上限，不自动延长、改语法或分块近似。每条输入约20–31MB、预测输出上限128MiB；整个新cohort预留1GiB以上存储，不落盘全部神经分数。

唯一推理root `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW/<seqid>/`：JOBID、STATUS、started.json、input.fasta、summary.json、prediction/{summary.json,full_grammar.gff3}。summary覆盖双strand全grid、最优分数、skyline统计、编码/解码时间及GPU占用。唯一CPU诊断root `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC/`：JOBID、STATUS、summary.json。日志为 `outputs/M28-PILOT-R5/logs/M28FULLTRAINRAW_<array>_<task>.out/.err` 与 `M28FULLTRAINDIAG_<job>.out/.err`。CPU诊断只在整个数组afterok后启动，2CPU/4GiB/20分钟、0GPU。

只有所有六元素root/batch/extern成功、所有规定产物与完成字段齐全，再加CPU诊断成功，才算端到端工程完成。运行中无final summary完全正常。具体实施错误可以在同一科学范围内做一次最小修复并新root重试，保留失败；OOM/TIMEOUT/真实科学问题先报告，不滚动加时或自动改协议。

## TRAIN端到端读数

全部预测先完成，再由独立CPU诊断打开且只解析六条TRAIN reference，复用既有M26 CDS-assessable尺：逐gene在完整、可翻译CDS候选中选最长，长度平分取transcript ID；基因/UTR partial不等于CDS partial。不得按M28训练eligible、motif、cap、长度或same-strand overlap删除参考。实际phase不一致会明确报错，不默默缩分母。

先核对raw身份、最终checkpoint、双strand全grid、GFF导出数量/坐标/actual DNA与phase/cap/同strand gene-set不重叠，再计算unique完整CDS链P/R/F1（key含seqid/strand和所有CDS边界，phase另诊断）、不分strand CDS-base union P/R/F1、背景跨度与CDS覆盖、纯背景链/Mb、非primary但命中其他assessable isoform、重复链数。背景是所有annotated gene feature跨度任一strand之外的bp；span_FPR是背景bp覆盖比例，不是gene-level分类FPR。报告逐染色体、逐物种、pooled计数、species macro与worst-species，不以合并值替代物种失败。unsupported参考分层保留。

这是训练染色体上的描述性端到端读数，不能证明独立泛化、原生基线优势或论文成功；即使未经draw曝光的同TRAIN染色体区域，也不是独立held-out测试。不设置新准确率gate或根据本次结果改变训练参数。

首训是有限原context下的composite观测目标，不能因为本轮做整染色体joint推理而称为chromosome-joint MLE；旧context外参考警告保留，不由本次推理覆盖。后续方案见[full_grammar_dev_baseline_plan.md](full_grammar_dev_baseline_plan.md)。终态及实际结果见下面追加记录；不依新结果更改本合同的冻结科学条件。

## 等待和接续

只复用原监控对话 `01a0f71f-8845-7432-bd13-bf54c1808fb3/local`，主对话 `01a0332e-e0a9-7aa2-90f9-159405701833/local`。旧心跳已被删除且当前无重复coding-rna自动化；本轮用户明确允许只重建一个绑定原对话的heartbeat。保持15分钟、quiet RUNNING/PENDING、单次终态交接、不动主对话模型、不新建/fork监控。监控只做轻量调度/完成结构检查，主对话负责指标和科学解释。交接后暂停并保留该固定心跳用于下一cohort，不删除。配置、单次执行、定时运行和完整回传分别如实核验。

## 本轮执行记录

生产接入检查13443304成功后，`sbatch --parsable sbatch/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW.sbatch`返回数组13443323；`sbatch --parsable --dependency=afterok:13443323 sbatch/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC.sbatch`返回13443329。六元素已在gpu034/gpu035同时RUNNING，CPU诊断PENDING/Dependency。输入总计156779632bp，原全grid双strand总计25512个窗口；无训练更新。首条已观察到真实native encoding进度，尚无整条预测、GFF或科学指标，不能提前记端到端完成。

提交前fresh五项Slurm检查及存储/输出不碰撞检查通过：private两节点配置各8卡、已用5/4，共7张3090可分配；无影响该GPU窗口的maintenance reservation，BeeGFS约152TiB可用，当前用户作业数远低于80，数组6远低于100。已有TE作业未改动。选择private而非等待共享特定GPU队列，不更改科学命令；正式资源仍为每元素1GPU/2CPU16GiB/12小时。登录节点和Mac未运行Torch/NumPy/native forward。

旧ID `coding-rna-r5-light-monitor`缺失且无重复后，工具创建成功复用了这个可用ID并绑定同一原监控对话，而非创建新对话。已保存heartbeat/ACTIVE/15分钟/failed_runs_only/目标对话；旧对话最新实际turn曾漂移为gpt-6.1-sol/low，已按用户指定恢复仅该监控对话为gpt-5.6-luna/low，新turn `01a11cb9-7480-7d92-9a2a-f30e6b0a638d`确认luna/low，主模型未改变。同一对话的初始化检查于2026-10-08T18:15:01Z完成：六推理RUNNING、诊断PENDING，未发正常进度接续；它早于journal物化时读到memory缺失，现已补齐current-cohort记录。

首次真实定时回合 `01a11cc7-bf12-7e52-a96e-b5d7636baec8`已在同一原监控对话于2026-10-08T18:30:15Z触发、18:30:29Z完成；实际user message为带同一automation ID的heartbeat，turn context仍为gpt-5.6-luna/low。它执行最小squeue与独立sacct，六推理RUNNING、诊断PENDING/Dependency，无终态或需用户处理情况，未发主对话接续。首次回合未写journal，由主补齐实证并在同一心跳明确每轮更新本cohort记录及正常状态仅DONT_NOTIFY格式；不改变对象、间隔、模型或科研范围。已区分配置成功、人工初始化成功、首次定时运行成功；本cohort完整自动终态回传及后续记录修正的实际执行仍待观察。主对话不另外轮询等待Slurm，不暂停或缩小原goal。

## 2026-10-09终态接续（以上RUNNING为当时记录）

原监控于2026-10-09T02:21:39Z实际完整回传：数组13443323_0–5及CPU诊断13443329全部root/batch/extern COMPLETED/0:0、规定产物完整；同一heartbeat已PAUSED并保留，journal记Handoff YES。本cohort定时检查和完整自动终态回传均已证实，不新建对话／重复自动化。主对话独立核对sacct、实际summary与指标公式／计数汇总；端到端工程完成。

六推理实际在UTC 2026-10-08 18:13:15–18:51:29完成，CPU诊断18:51:31–18:52:29，整流程wall39分14秒。拟南芥完整链F1=50.4871%、水稻28.2780%、pooled45.5901%，不分strand CDS-base F1分别93.0848%／80.8511%／90.3601%。工程完成与完整gene科研成功不同；全部结果与限制见[full_grammar_raw_train_result.md](full_grammar_raw_train_result.md)。冻结DEV／基线方案已准备但未执行；训练参数、最后权重选择及封存数据边界不变，原goal未完成。

原TRAIN输出的补充分层13449803（2CPU/4GiB/5分钟、0GPU）root/batch/extern COMPLETED/0:0，实际wall100秒。单／多CDS段、原1536 draw中R4 positive_ids完整链几何出现／未出现及span>24576bp分层完成，所有exact计数与原共同尺完全相等，不读DEV或模型，不改预测。多CDS段F1拟南芥43.9931%、水稻18.5934%，曾出现metadata正链组的完整召回率51.6484%/29.8611%；9条长span参考均未命中。它们定位描述性结构缺口，不是因果消融、新科学gate、独立测试或改参依据。原 `longer_than_native_window`字段只表示长度分层，whole-contig语法并不自动排除长链。唯一新root `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE/`保留summary及reference_chain_profile.jsonl，未覆盖任何冻结结果。
