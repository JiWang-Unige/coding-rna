# 原 TRAIN 原生基线比较

2026-10-09。本轮补齐 ANNEVO、Tiberius、Helixer 在原六条 TRAIN 染色体上的固定版本原生结果，用于判断当前模型在训练范围内与成熟调用器的差距。比较不用于证明独立泛化，不改变原八方法 DEV 家族或另列十方法扩展，不据此在 DEV 选择参数。

## 输入和运行资产

直接只读复用 `outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW/<seqid>/input.fasta`，与原4608步 TRAIN 推理使用相同文件。原输入已转为大写，非 A/C/G/T 字符统一为 N；此次不重新提取基因组，不引入软掩码、参考、蛋白或 RNA 证据。完成前不打开真实参考 GFF。

| 物种 | 原 TRAIN seqid | 长度 bp | 固定 primary 参考数 |
| --- | --- | ---: | ---: |
| 拟南芥 | NC_003070.9 | 30427671 | 7135 |
| 拟南芥 | NC_003071.7 | 19698289 | 4293 |
| 拟南芥 | NC_003076.8 | 26975502 | 6311 |
| 水稻 | NC_089043.1 | 27474823 | 1638 |
| 水稻 | NC_089044.1 | 23931887 | 1601 |
| 水稻 | NC_089046.1 | 28271460 | 1710 |

复用原 DEV 对比已验证的隔离资产及操作点：ANNEVO 源码42c920f、Magnoliopsida权重、batch32；Tiberius 隔离2.0.8源码871d497及实际2.0.7依赖 SIF、angiosperms权重、batch4、seq_len200070、unmasked且无ClaMSA；Helixer 0.3.7源码d17bb49、land_plant_v0.3_a_0080权重、batch16、64152 subsequence、deterministic、原生后处理默认。详细身份与暴露限制沿用[native_baseline_preparation.md](native_baseline_preparation.md)和[pretrained_exposure.md](pretrained_exposure.md)，不把各工具的上游暴露统一称为未见。

`run_native_baseline.commands` 仅增加明确 FASTA 参数。未传参数时原 DEV 命令保持不变；TRAIN runner 通过只读单文件容器挂载指定共享 FASTA，其他参数逐项保持原样。原生 GFF、日志、GPU 探测和权重收据分别保存在新目录，不覆盖旧 DEV 或原 TRAIN 结果。

## 评价与决策

18组预测全部完成后，才用原 M26 CDS-assessable、最长 eligible CDS、transcript ID tie-break 共同尺打开参考。固定分母为拟南芥17739、水稻4949、合计22688。保留 unsupported、重叠、长链与其他困难参考，不为原生模型或本模型单独删例。逐染色体参考数、eligible isoform 数、unsupported 分层、CDS union 参考碱基与背景碱基必须和原4608步 TRAIN 评价一致。

输出完整链 precision／recall／F1、CDS-base union、背景 span、重复链、其他 eligible isoform 命中、phase一致性诊断，以及逐物种、pooled、两物种 macro 和 worst-species 汇总。phase不是坐标 exact-key，不因原生输出不符合本模型的128-exon上限而删预测。相对原4608步控制的差值是 TRAIN 描述性差距，不是目标函数因果证明；入口计时不作公平速度 benchmark。

下一判断依据是：成熟调用器在相同 TRAIN 输入上是否仍明显领先，以及该差距在两物种、CDS 碱基和完整结构层面是否一致。原有六条训练分支及四格汇总继续按自己的依赖执行，此比较不等待它们，也不改变候选、训练、DEV启动条件或封存数据边界。

## 执行和终态

预检使用 private-teodoro-gpu，0GPU／2CPU／8GiB／10分钟，检查合成范围、共享输入挂载和原参数不变、真实 parser 汇总、完整性闸与分母不删，以及原六份 FASTA 实际长度／编码和既有资产是否齐全。

预测数组为18格、最多6格同时执行；每格1 RTX3090／8CPU／32GiB／30分钟。原两条 DEV 的同工具运行约90–360秒，最大 TRAIN 染色体30.43Mbp与原 DEV 最大29.94Mbp接近，30分钟留有余量。提交时依据实际可用 RTX3090 选择最早完成的分区；不因经济成本延迟。预计原生临时文件及输出为数十GiB，使用独立目录，保留错误产物，不取消其他作业。

评价使用 private-teodoro-gpu，0GPU／2CPU／8GiB／20分钟。只在18元素全部 afterok 后运行；任何失败或依赖不可满足都交回主对话处理，不自动重提。完成证据要求每元素 Slurm 成功终态、原生完成标记、真实 GFF 和 summary／STATUS，及评价终态、22688分母和完整汇总，不把 GPU probe 或预检通过当作准确率结果。

输出目录分别为 `outputs/M28-NATIVE-TRAIN-BASELINE-CHECKS`、`outputs/M28-NATIVE-TRAIN-BASELINES` 和 `outputs/M28-NATIVE-TRAIN-BASELINE-EVAL`。长等待仅扩展原固定监控对话和心跳；不建立第二个监控。作业 ID 与实际终态在提交后补齐。

首个预检13461185于15:03:08–15:03:33 UTC运行，root/batch FAILED1:0，extern COMPLETED0:0。12项合成测试全部通过；随后检查脚本错误地从 TRAIN runner 读取原 runner 的资产常量而触发 AttributeError，未进入真实 FASTA 核验，未运行任何原生模型。修复仅改为原 runner 的正确常量路径；首轮产物和日志保留。修复后的预检使用独立 `outputs/M28-NATIVE-TRAIN-BASELINE-CHECKS-R2`，不覆盖首轮，也不改权重、命令参数、参考尺或科学范围。

修复后预检13461217的root/batch/extern均COMPLETED0:0，15:08:19–15:08:45 UTC，共26秒；12项合成测试、六份实际输入的scope／长度／ACGTN编码和固定资产均通过，无GPU前向、真实参考、DEV或封存访问。中间一次SSH在banner阶段超时，重新取证确认无R2输出目录及同名作业后才提交，没有重复任务。

正式命令为 `sbatch --parsable --dependency=afterok:13461217 sbatch/M28-NATIVE-TRAIN-BASELINES.sbatch`，返回数组13461306；`sbatch --parsable --dependency=afterok:13461306 sbatch/M28-NATIVE-TRAIN-BASELINE-EVAL.sbatch` 返回13461307。提交后首六格实际RUNNING于shared-gpu/gpu021，其余12格为JobArrayTaskLimit，评价为有效Dependency；这只是运行收据，尚无新的准确率结论。提交前用户展开队列58项，加19项后77，未超过80；共享gpu021有7张空闲24GiB RTX3090和充足CPU／内存，私有gpu034已满、gpu035仅2张空闲；无相关GPU维护预约，BeeGFS可用约164.5TB，选择共享分区以缩短完成时间。

原固定心跳于15:12:49.862 UTC扩为六条未完流程，原监控对话、gpt-5.6-luna／low、15分钟间隔、通知偏好及此前五条依赖均保留，保存配置已精确回读。该更新只增加13461306／13461307和其结构白名单；合成检查及已消费失败不交给监控重复处理。后续实际单次检查、定时检查和终态回传须分别核实，不用配置保存成功代替运行证据。原论文goal仍未完成。

新范围独立检查已实际追加到原监控对话，turn01a1213b-af31-77b0-83cc-fdcaecbe1b76，15:15:22–15:15:59 UTC；实际turn_context为gpt-5.6-luna／low。这是原主触发的manual检查，不是定时运行。15:15:50的调度核实显示新数组0–5成功完成、6–11运行、12–17受并发上限等待、评价有效Dependency；其余五条流程仍正常，未读产物／指标、未交接或暂停。

更新后实际定时检查也已追加到原对话：turn01a1213d-1fb4-74f1-9038-fd50fff0a17d，15:16:56–15:17:22 UTC，实际turn_context为gpt-5.6-luna／low；heartbeat输入和真实调度命令均包含新13461306／13461307。15:17:13核实六条流程正常、未读产物／指标、未交接或暂停。配置保存、独立manual检查和新范围scheduled检查均已分别验证；整链终态及最终科学解释仍未完成。

本轮整链现已闭合：18格13461306及CPU评价13461307的全部实际root／batch／extern均COMPLETED0:0，评价UTC15:31:19–15:32:40。原监控在实际scheduled turn01a1214b-52ef-7ea3-991b-efb7ac6653b2回传结构终态并于15:36:06UTC暂停同一心跳，其他五cohort仍未消费。原主已经验证并解释共同尺结果，见[native_train_baseline_result.md](native_train_baseline_result.md)；追加短CPU逐参考诊断13461637亦直接完成消费，见[native_train_reference_diagnostic_result.md](native_train_reference_diagnostic_result.md)。原生准确率比较已完成，但原论文goal未完成。
