# M28 全 TRAIN 等曝光训练对照协议

## 实验问题与边界

本轮检验在充分覆盖已有 TRAIN 后，学习式完整基因 query 生成相对局部分类加固定 C0 解码是否仍有优势；同时记录扩大覆盖与旧 pilot 的差别。它是两个新权重的开发集实验，不是修改、延长或救回正在运行的 A3。

特征沿用已通过 12 个逐元素 BF16 回放的官方冻结 GLM，标签仍为 R4，六条 TRAIN 染色体与所有未知/部分/alternative 保护政策不变。自然窗口每遍 25,512 个、三遍总计每臂 76,536 次更新。没有新物种、封存/test/Setaria、DEV 模型选择、经济成本 gate 或对外发布。两臂是等数据与等曝光对照，不是等参数、等目标函数或等运行时长。

## 两个固定模型

Q 臂复用 A3 的 CollaborativeGeneQueryModel：D2048、H128、Q16、649,457 个参数，四头 pre-norm 残差 query 通信，dropout0、无 FFN。损失为原 A3 的 local 加 known-query-mean 存在性 BCE、边界 CE、EOS BCE，四系数均为1；匹配目标与未知保护不改。不是新的结构创新，能否成为论文模型取决于后续性能与机制证据。

C0 臂复用 M28Core(C0)：同 SharedContext 与15类 local_loss，366,191 个参数。两臂从相同 seed0 重新初始化，共享上下文及 segmentation 的初始权重逐元素相同，不继承旧监督 checkpoint；训练差异包含 query 结构及附加目标。

## 固定输入与优化

每遍按已生成的 training_order.jsonl 完整读取每个自然窗口一次，采用独立 Python random.Random(0) 的固定打乱顺序，三遍重复同一顺序；不以各臂随机抽样制造不同曝光。

权重为 N_total×原 stratum draw 数/(1536×该 stratum 自然窗数)，固定使用 importance_weight_original_stratum_objective。保持原抽样对应的染色体/链向期望风险、两物种各占一半权重；不声称裁剪梯度/AdamW 更新与旧抽样一致。

两臂同 batch1、FP32、AdamW lr3e-4/wd0.01/betas0.9和0.999/eps1e-8、clip1，三遍，无学习率/seed/epoch/阈值搜索。标签重新来自 R4 的 TRAIN-only metadata/catalog，旧缓存里的历史标签字段不能用于监督。检查所有 feature geometry、窗口、stratum、正例资格及 join；不解析 mixed catalog 中 DEV 行。

## 执行顺序

CPU readiness 在特征仍生成时独立完成：七项针对范围、完整顺序、权重、实际损失路由、共同初始化和旧 DEV 缓存拒绝的检查；生成全 TRAIN 的固定最密正例、最多正例外显子、零正例无保护、正例有保护哨兵（每物种四类）。

全部 16 特征分片及 CPU finalize 成功后，GPU smoke 在同一全 TRAIN 特征上运行两臂。每臂固定训练顺序前64窗，加八个元数据选择的极端/保护哨兵，共72次真实前向和反向；无优化器更新、不保存权重。检查有限损失、有限梯度、GPU容量，记录 prefix 平均耗时与保守运行规划所需证据；不能将吞吐烟测当作性能结果或科学 PASS。

只有真实 GPU smoke 成功后才提交两臂从头固定拟合。操作时窗按真实烟测及旧拟合记录确定；旧 A3 4,608步为31:16，旧 C0 4,608步约3,278秒，线性外推本轮分别约8.66小时、15.13小时，仅为规划参考，不是性能保证。资源预计每臂1 RTX3090、2CPU、16GiB，长任务使用 private，实际提交前重新核实调度/维护/输出冲突。2026-10-03 原生烟测成功后已提交两组新拟合，具体资源与依赖见执行记录。

## 评价与决定

只使用最后 step_076536 权重，禁止根据 DEV 选择 epoch checkpoint。Q 保持原 DNA 语法、tau0、max_exons128、nearest-owner/去重/正 gain 同链向非重叠选择；C0 保持现有固定 DNA 语法和原生 HMM 解码，不能因全 TRAIN 训练换用 oracle 或新阈值。两臂完整无参考 DEV 推理完成后才进行统一参考评价。

主参考仍 7,728 条（拟南芥5,437、水稻2,291），8690个方向窗口；unsupported、长链、同链向重叠参考不删。旧工具结果只作历史比较，不声称完成最新原生对标或独立泛化验证。

预先固定候选推进条件：Q 两物种 exact-chain F1 均严格高于本轮等曝光 C0，两物种背景 span 均不高于本轮 C0；同时完整报告逐物种/micro/macro、CDS-base union、free→owner→selected损失与 query 冗余。旧 A/A2/A3、旧 C0 的比较另列并标明曝光差异；不能把单物种/仅召回改善或未达背景条件包装成成功。不通过则保留结果，不扫阈值、延长曝光或换 checkpoint 救本轮。通过后仍需最新原生 comparator、冻结独立测试及与期刊目标相称的其他证据，不能立即宣称发表条件达到。

## 执行记录

训练入口 run_full_train.py，原始 FASTA 推理入口 infer_full_train.py，共同评价入口 evaluate_full_train.py。新目录包括 M28-FULL-TRAIN-READINESS、GPU-SMOKE、Q-FIT/C0-FIT、Q-INFER/C0-INFER 和 EVAL。保留所有旧缓存、checkpoint、冻结输出和不相关 Git 改动；不提交 Git。

2026-10-03 07:48 UTC 已核实的执行事实：

- CPU readiness13357989：COMPLETED/0:0、29秒；七项检查通过（1.279秒）。全窗25,512及唯一正链22,407完成 TRAIN-only join；两物种每遍有效权重分别12756.000000002032、12756.00000000133。八条固定哨兵已写入 witnesses.jsonl，最密窗11/8条正链，最多正例外显子窗73/51个外显子；不是模型性能结果。
- GPU smoke13358001：PENDING Dependency，实际提交命令为 `sbatch --parsable --dependency=afterok:13357945:13357989 sbatch/M28-FULL-TRAIN-GPU-SMOKE.sbatch`，private-teodoro-gpu/1 RTX3090/2CPU/16GiB/10分钟。Slurm已满足ready依赖，只保留未满足的finalize依赖；尚未进行原生烟测、没有新拟合结果。
- 特征数组13357939_0–3已经COMPLETED/0:0；其余分片继续运行/按并发上限排队，汇总13357945等待全数组成功。只查首个已完成分片的summary作运行取证：1595个窗=1491个原生新增+104个旧缓存复用，409.30秒、backbone可训练参数0；不能据一个分片推定全缓存完成。
- 已完成两份新增Python AST、两份新增sbatch语法与五份新文件远程精确回读。原A3训练/推理实现、参数、权重和评价口径未修改。

等待复用同一监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同一heartbeat coding-rna-r5-light-monitor，gpt-5.6-luna/low、每15分钟。当前cohort M28-A3-FULL-TRAIN-SMOKE-13358001，包含原A3链和全TRAIN特征→汇总→smoke链。任一整链到达完整终点、已使主对话有下一授权动作时，只交接一次并暂停，明确另一条仍未完成的作业；新拟合不由监控启动。配置成功、人工检查、当前cohort定时及完整回传分别记录。

监控已修复实际数组识别问题：sacct使用JobID列匹配13357939_0–15，保留JobIDRaw仅作内部编号（例如_0对应13357940）；不凭父数组或一个分片成功判定全量完成。

## 原始 FASTA 推理和评价入口

infer_full_train.py 固定加载对应臂最后 step_076536 权重，校验完整拟合、配置和冻结编码器版本。正式 DEV 路径只物化已有两条允许的染色体，对全部自然窗和两条链向执行官方四块 BF16 编码，再送入新 FP32 头；不读取参考、标签或 DEV 特征缓存。Q 不调用 ANNEVO HMM，C0 明确保留原生 HMM 研究依赖及其原许可边界；公开发布另行处理。

所有 free 候选及分数写入逐窗 JSONL，组装仍使用原 nearest-owner、全局去重、正 gain 同链向非重叠选择，再写 primary CDS GFF3。evaluate_full_train.py 在两臂终态推理齐全后才调用冻结参考评价，检查8690窗、四完整scope、原GFF与重新组装的坐标和phase，记录逐物种/pooled指标、free→owner→selected真实参考损失和候选冗余。旧 A/A2/A3/R5 C0 单列为不同曝光的历史结果，不用于本轮等曝光 gate；尚未完成的旧结果标记为不可用，不能当成失败或成功。

历史 common_ruler.py 的 reference_nonexact_primary 字段实际上是非精确匹配的预测数量 len(pk-rk)，不是漏掉的参考数量。新报告将它正确命名为 prediction_nonexact_primary；指标数值与旧 A3 实现不改动。合成 GFF 检查还触发旧 parse_annotation 的未关闭文件 ResourceWarning；该警告未使检查失败，没有据此扩大修复范围或改变运行中的旧实验。

## 2026年10月3日正式拟合启动

08:08 UTC 已核实：

- 特征数组13357939全部16元素COMPLETED/0:0；汇总13357945 COMPLETED/0:0、16秒。根索引完整25,512窗，24,021新增/1,491只读复用，12个stratum；编码器与标签边界不变，DEV_forward_used/test_setaria_used均false。
- GPU smoke13358001 COMPLETED/0:0、68秒。Q、C0各72次真实前向/反向，144条记录损失与梯度有限，优化器更新0；Q峰值GPU分配312,109,568字节，C0为165,360,640字节。这是训练就绪证据，不是准确率。
- prefix平均耗时C0 0.48325秒/Q 0.25901秒，76536步外推约10.27/5.51小时，不含AdamW且不是保证；旧完成拟合线性外推15.13/8.66小时。据两类证据申请C0 18小时、Q 10小时的private连续时窗，不改变更新数或科学停止条件。
- 初始推理检查13358066 COMPLETED/0:0、11秒，八项检查通过；增加不同曝光历史参考口径检查后，以新只读结果目录执行13358075，COMPLETED/0:0、12秒，九项检查通过（0.437秒）。保留初次结果；未反复审阅未变化的补丁。
- 首次工程检查提交因不存在的private-teodoro-cpu分区被Slurm拒绝，没有创建作业；改为项目规定的private-teodoro-gpu、零GPU后成功。没有科研重试、原始结果覆盖或原A3参数改动。

| 阶段 | 作业 | 资源时窗 | 成功依赖 |
| --- | --- | --- | --- |
| Q固定拟合 | 13358083 | 1 RTX3090 2CPU 16GiB 10小时 | 13358001与13358075 |
| C0固定拟合 | 13358084 | 1 RTX3090 2CPU 16GiB 18小时 | 13358001与13358075 |
| Q原始FASTA DEV推理 | 13358085 | 1 RTX3090 2CPU 16GiB 3小时 | 13358083与13358075 |
| C0原始FASTA DEV推理 | 13358086 | 1 RTX3090 2CPU 16GiB 3小时 | 13358084与13358075 |
| 共同评价 | 13358090 | 0GPU 2CPU 8GiB 20分钟 | 13358085与13358086 |

实际提交命令均为 `sbatch --parsable --dependency=afterok:<表中依赖ID以冒号分隔> sbatch/M28-FULL-TRAIN-<阶段>.sbatch`；对应文件为Q-FIT、C0-FIT、Q-INFER、C0-INFER、EVAL。所有作业使用generanno，日志在outputs/M28-PILOT-R5/logs，以作业名和ID区分。输出目录均此前不存在；存储可用175,257,609,568,256字节，无新增经济成本gate。

即时Slurm核实两fit已在gpu034并行RUNNING，各日志已出现实际step64/76536，训练循环约19.2秒；推理和评价等待成功依赖。这只证明训练已经开始，不是完成拟合或模型优势。原A3推理13357568仍运行、评价13357569等待，不取消、不改参数。

## 固定监控接续

同一监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3、同一heartbeat coding-rna-r5-light-monitor已更新为cohort M28-A3-FULL-TRAIN-FIT-13358083-13358084；保留gpt-5.6-luna/low、每15分钟、ACTIVE，不新建或fork对话。旧特征/烟测终态已由主对话直接验证并消费，不伪称监控曾发送交接，也不继续监控旧分片。

任一完整A3链，或完整两fit→两raw infer→共同eval链到达终点，或出现需处理的真实失败，监控只向原主对话交接一次并暂停，明确另一条仍运行的链；主对话负责科学验证和后续授权工作。配置、人工检查、当前cohort定时检查及自动终态回传分别记录；旧smoke cohort08:00同对话定时运行已验证，不能充作新fit cohort的定时或回传证明。

## 原生推理接口检查

2026-10-03，本轮额外检查真实原始DNA编码与新头自由解码的接线，避免仅凭缓存训练烟测或合成接口检查推定原生推理已经可用。入口smoke_full_train_native_inference.py直接复用infer_full_train.py的annotate_contig、window_record及GFF出口；不修改正式拟合或末checkpoint主入口，不伪造训练步数/权重身份。

固定使用已有replay_sentinels.jsonl中每物种第一条TRAIN哨兵，共两个24576bp原始片段，各走正反链；Q/C0合计八个方向窗口。只读两条已允许TRAIN染色体的FASTA和既有窗口几何，编码器仍为原revision/四块BF16，头为固定seed0未训练初始化、eval且requires_grad false。无参考/标签/特征缓存或checkpoint读取，无optimizer更新、权重保存、DEV前向或test/Setaria访问。

检查对应真实失败：原生编码器/新头/解码接口若不兼容、候选数与分数错配或非有限、完整窗口缺记录、Q16/tau0改变、输出坐标/phase回读改变，则本检查失败并据日志修正具体实现，不能把错误留给数小时后的正式推理。空输出对未训练头是允许结果，不把候选数或结构合法性当准确率；全染色体与正式训练末权重的性能仍待原成功依赖链。

调度入口M28-FULL-TRAIN-NATIVE-INFERENCE-SMOKE.sbatch：generanno/private-teodoro-gpu/1 RTX3090/2CPU/16GiB/3分钟，新目录outputs/M28-FULL-TRAIN-NATIVE-INFERENCE-SMOKE，日志M28FTNATIVE_<job>.out/err。现有拟合/推理/评价作业的资源、依赖、权重和科学规则不改变；本次没有新基线运行、权重下载、Git提交或公开发布。

### 实测终态与范围

原生接口检查13358115已COMPLETED/0:0，Slurm分配25秒、1 RTX3090/2CPU/16GiB，STATUS COMPLETED。固定两个TRAIN原始片段的正/反链，在两臂合计八个方向窗口全部走通；C0/Q分别6/4条最终CDS链，两份GFF3坐标与phase回读完全一致。这些是未训练头的工程输出计数，不是准确率。

本检查native循环11.4432秒，C0/Q头均零可训练参数、optimizer更新0，backbone可训练参数0；峰值GPU分配分别2411252736/2457264640 bytes。reference_or_labels_used、feature_cache_used、checkpoints_read_or_written、DEV_forward_used及test_setaria_used均为false；没有保存或修改正式权重。证据为outputs/M28-FULL-TRAIN-NATIVE-INFERENCE-SMOKE/summary.json。完整DEV以及真实step076536的最终预测质量仍待原五作业依赖链，不能由八窗烟测推出。

A3原小曝光整链13357566/13357568/13357569已完成并由主对话验证：事前条件失败，结果保留在[a3_result.md](a3_result.md)。不追加A3轮次或改变本full-TRAIN协议。固定同一监控对话及heartbeat已缩小为cohort M28-FULL-TRAIN-13358083-13358084，仅等待Qfit13358083/C0fit13358084/Qinfer13358085/C0infer13358086/eval13358090；模型gpt-5.6-luna/low和15分钟间隔不变。此前08:15定时检查是旧合并cohort证据，不能代替新cohort定时或自动回传验证。
