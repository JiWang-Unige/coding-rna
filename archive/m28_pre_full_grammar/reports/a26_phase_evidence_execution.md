# A26 相位状态证据执行记录

2026-10-06（本次取证UTC 2026-10-05 17:23）。实验按[a26_phase_evidence_protocol.md](a26_phase_evidence_protocol.md)固定：fresh seed0，总248096参数，新增28；原A3 generator/15类局部分支冻结，原1491 TRAIN窗/16321自由候选、1536 draw三遍4608更新、末004608、原集合NLL及unknown保护、阈值0/五gate均不变。不是A21救参，不改变A23/A25，也没有DEV/sealed/Setaria、新候选搜索或公开发布。

## 已完成的工程证据

预检13392186 root/batch/extern均COMPLETED/0:0，38秒，gpu035。五项实际测试PASS，1.060秒；原metadata八TRAIN窗检查18.493秒、optimizer0。八窗全部自然候选证据有限，fresh共有A21参数逐tensor相同，零初值项下完整scores和原集合loss精确一致。八窗新增投影梯度均有限且非零，覆盖两物种；局部generator无梯度，head参数未变，unknown独立gain直接梯度0。水稻一个实际case仅单exon已知候选，因此两个junction分支应无梯度；该真实支持边界按协议处理，不硬要求全部分支每窗都激活。

证据目录outputs/M28-A26-PHASE-EVIDENCE-PREFLIGHT，实际summary.json、STATUS/JOBID及logs/M28A26PREFLIGHT_13392186.err。预检不是全1491窗缓存验证，也不是性能结果。

## 固定运行已启动

合并GPU fit/末TRAIN readout：13392225，private-teodoro-gpu、1RTX3090/2CPU/16GiB/60min，17:22取证RUNNING gpu035。CPU结果核验13392226，afterok:13392225、0GPU/2CPU/8GiB/10min，PENDING有效Dependency。提交前检查五项当前调度/预计启动/维护状态、实际storage与目标目录不存在；未取消历史作业或扩大资源。

运行入口scripts/experiments/M28-METHOD-RESTART/run_phase_evidence_chain_set.py，复用本次八窗预检的new_head、evidence_inputs、score_loss、zero_update_witnesses，而不沿用旧头收据。runner AST及两SBATCH语法通过，相关新增源码逐文件与远程字节比较相同；未提交Git。

独立阶段目录：
- outputs/M28-A26-PHASE-EVIDENCE-FIT
- outputs/M28-A26-PHASE-EVIDENCE-TRAIN-READOUT
- outputs/M28-A26-PHASE-EVIDENCE-RESULT-CHECK

原FP32 hidden保留18,761,121,792字节，新增固定28维证据只需1,827,952字节；不保留1491窗全部15类logits。GPU job先执行显存预算检查、materialization及八窗新头live/cache score/loss/梯度等价检查，均在optimizer0阶段；之后才优化。提交前已完成实现和基于原实际pipeline峰值的显存估算，全量新缓存实际峰值与新live/cache收据由本GPU job零更新阶段记录，尚未完成，不能把八窗预检称为全缓存验证。

每256工程checkpoint，主读出只有末004608。训练额外记录新投影实际梯度支持、末权重norm；CPU核验原4608顺序/importance/计数、末head有限/结构、八窗等价收据、全部候选几何与证据行顺序，并用保存28维向量和末投影权重直接重算phase term及nonadditive分解；不跑新GLM、生成器或head推理。CPU重算原selection/profiles和五gate，辅助配对已闭合A21，不以A21门槛取代A6。

## 等待与未完成项

同一watcher 01a0f71f-8845-7432-bd13-bf54c1808fb3与heartbeat coding-rna-r5-light-monitor已扩展为M28-A23-A25-A26-TRAIN-20261006，保留A23三阶段13391608/09/10与A25三阶段13392060/61/62，新增A26两job。原gpt-5.6-luna/low、15分钟、quiet terminal-only交接不变，17:22:47.369Z原生保存及完整prompt/id/kind/target/cadence/notification/created_at读回一致。没有新建、fork或cron。

扩展cohort单次同对话startup已完成：取证2026-10-05T17:23:23.065Z，turn01a10d17-414b-7940-a9e4-d2234c21946a completed/error-null/24.391秒。随后新cohort定时SAME CHAT取证亦已实际完成于17:24:17.441Z，turn01a10d18-22ec-79c0-8a94-50536846ffc0 completed/error-null/32.180秒；A23/A25 fit与A26合并root均运行，下游有效Dependency，未回传终态。该定时turn实际turn_context17:24:08.091Z仍为gpt-5.6-luna/low。配置成功、单次检查、定时运行和最终自动交接分别取证，只有最终终态交接尚未观察；不是沿用旧cohort的定时证据。

当前尚无A26 optimizer完成、实际末TRAIN选择或gate结论。FAIL按冻结规则闭合，不扫参数救参；PASS也只支持下一独立固定raw-FASTA caller/DEV协议。完整raw-FASTA同输入原生ANNEVO/Tiberius/Helixer比较、独立适用范围、机制消融与论文级创新证据仍未完成。工程验证或TRAIN机制收益均不能称为泛化或Nature Communications级模型成功。

## 训练终态和结果核验顺序修复

2026-10-06，本次取证UTC 2026-10-05 17:45–17:48。合并13392225的root/batch/extern均COMPLETED/0:0，19分14秒，batch MaxRSS6410828K。原FIT与TRAIN-READOUT的STATUS/JOBID和summary齐全；末004608已写出，原4608更新完成。全1491窗缓存与八窗live/cache接线检查已经执行，不再是仅有预检：实际GPU峰值18893957632字节，新投影非零梯度4599步、末权重norm0.15060359137180174。它们是训练与接线证据，不是性能成功。

原CPU核验13392226 root/batch FAILED/1:0、20秒；extern COMPLETED/0:0，STATUS FAILED exit=1。stderr指向run_phase_evidence_chain_set.py第233行：保存证据的window_id列表与list(originals)顺序不相等。训练materialization明确按windows.items()写出phase_evidence.jsonl，natural readout按records.items()写出windows.jsonl；两个字典来自不同原始顺序，训练和读出均以同一window ID取证据。原检查器将两种遍历顺序错误地当成同一种，不据此推断训练错配，亦不跳过身份检查。

唯一有界修复仅把证据顺序断言改为与实际materialization的list(windows)比较，随后仍逐ID检查原候选完整链顺序、28维向量、末权重投影、全部组成分、4608 draw/importance、选择与五gate。没有放宽容差、排序后掩盖重复、修改模型／checkpoint、重新生成候选、训练或读取DEV。新增CPU入口M28-A26-PHASE-EVIDENCE-RESULT-CHECK-REPAIR-R1.sbatch，输出同名新目录，保持generanno/private/0GPU/2CPU/8GiB/10min。原失败目录与日志不改；原训练／读出仍只读。runner AST及新SBATCH语法通过，正式修复核验待实际完成，不能用源码判断替代全量核验。

原读出报告暂示两物种正确链增加且known_nonmatch也增加，事前五gate只有三项通过；在修复核验完成前只把它视作待复核的原始读出，不称已验证科学结果。固定heartbeat已临时PAUSED并精确回读，保留同ID、target、prompt、间隔、通知及created_at；这是主对话直接消费旧CPU故障，不是自动终态回传或平台goal暂停。核验提交后或终态消费后只用同一watcher接续仍在运行的A23/A25/A27，不新建监控。

修复核验唯一提交为 `sbatch --parsable sbatch/M28-A26-PHASE-EVIDENCE-RESULT-CHECK-REPAIR-R1.sbatch 13392225`，取得13392531；17:49 UTC立即取证RUNNING gpu034，elapsed15秒。两个改动文件与远程逐字节比较一致，实际输出先不存在，提交前五项调度／维护检查通过。原GPU作业和末权重不重提，原失败CPU记录保留。保存证据与读出的实际首个ID分别为拟南芥86016正链窗和29786112负链窗，进一步确认两种遍历顺序确有差异，而不是从报错猜测。

复用同一heartbeat加入新CPU13392531及REPAIR-R1白名单，旧13392226退出当前对象；cohort为M28-A23-A25-A26R1-A27-TRAIN-20261006。原生完整prompt、ID、名称、heartbeat、target、15分钟、通知与created_at已精确读回一致；新cohort manual startup已向原固定监控对话发送，尚待实际执行收据。这个配置不替代CPU完成证据，也不算自然定时运行验证。

修复终态13392531的root/batch/extern均COMPLETED/0:0、23秒、MaxRSS853600K，STATUS/JOBID一致，全部原数值核验及新增materialization顺序断言通过。实际五gate为三PASS／两FAIL，TRAIN_preparation_supported=false；[完整结果](a26_phase_evidence_result.md)据此关闭A26版本，不进入DEV。没有再次修复、训练或改参数。

本修复cohort的manual startup现场取证完成于17:50:15.846Z，原对话turn01a10d2f-ddeb-7790-8bf7-d61cb8d66086最终completed/error-null、74.838秒，实际turn_context仍gpt-5.6-luna/low。监控于17:50:45Z真实发出A26结构终态，主对话收到并消费；随后同一heartbeat暂停并读回PAUSED。配置、manual执行与终态交接均已验证；本修复cohort未观察自然定时触发，不能用17:39的旧cohort记录替代。A23读出、A25拟合及A27诊断仍运行，下一监控只接续这些原对象。

终态消费后，同一heartbeat已缩为M28-A23-A25-A27-TRAIN-REMAINING-20261006，全部A26退出，原七个剩余root／依赖／白名单不变。17:52:30.018Z保存，完整prompt及固定身份／间隔／通知／created_at精确回读；manual startup实际于17:53:14.903Z取证，turn01a10d32-9dc7-7be2-b152-c688eb1cc637 completed/error-null、23.268秒，仍gpt-5.6-luna/low。A23读出／A25拟合／A27诊断运行，其余有效Dependency；没有健康回传或新监控。新remaining cohort定时触发和终态交接尚未观察，按[官方同对话定时任务机制](https://learn.chatgpt.com/docs/automations)复用原对话，不新建独立cron。科学结果只见上述结果记录，不由监控解释。
