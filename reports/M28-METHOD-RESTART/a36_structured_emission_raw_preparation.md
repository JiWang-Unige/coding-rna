# A36 raw FASTA 注释工程准备

2026-10-06。A36整链TRAIN结果现已闭合，原五条件仅2/5通过，科学FAIL并关闭模型，详见[a36_structured_emission_result.md](a36_structured_emission_result.md)。以下raw-FASTA caller、四项人工序列检查、完整评价接口及两个Slurm入口是此前完成的工程准备，不放行失败模型，也不是实际原生注释结果。训练13406144、自然读出13408449和CPU核验13408451的原源码及结果均保留。真实DEV、独立物种、Setaria和公开发布均未启动。

## 候选生成和发射评分

infer_structured_emission.py沿用既有原生四块GEN编码、双方向、真实尾窗长度及nearest-owner／几何去重最大分数／span非重叠装配与GFF3导出。每个方向窗口只编码一次，冻结A3生成器和A36表示共享该GEN输出及DNA，但分别运行各自SharedContext。自然候选仍由原A3末004608的greedy＋beam4／completed4／cap128从输入DNA生成；新A36末004608的15类发射负责单位role-vs-intergenic FP64 gain，旧A3发射不参与新分数，阈值仍0。没有参考注入、推理标签、feature cache、长度归一化、温度或权重调节。

两个表示组件有额外参数、训练曝光和前向开销，不能宣称与A31容量或总计算匹配的消融。GEN共享只消除重复编码，不省去第二个表示前向；本入口也不解决原A3未生成的正确几何，仍只支持canonical primary CDS，不支持UTR、全部异构体或独立泛化主张。

显式FASTA为空时不能生成仅有header的虚假complete记录，重复contig ID也不能把不同序列并入同一GFF坐标系。新入口直接拒绝这两个输入错误，不改旧caller或增加通用防护框架；中途失败保留已生成诊断文件，不写成功summary。

## 原生执行条件

入口在模型、原生GEN、FASTA或输出创建前要求本轮三个实际源作业COMPLETED0:0、对应STATUS／JOBID／summary；重算两物种对A31_A3的原五条件，保留2916参考配对分母，且两份自然读出／CPU记录完全一致。实际CPU曝光、checkpoint、自然身份、保存分数和标签／选择核验均须成立。源证据只能证明计算准备，不自行授予新数据访问许可。

之后仍要求另行固定a36_structured_emission_dev_protocol.md。本轮不创建该真实DEV合同、不提交原生注释或评价；FAIL按原停止条件关闭，不以本工程接口或合成测试放行，不挑早期checkpoint、改阈值或重训。现有A33待答复例外不部署、不恢复。

## 四项针对性检查

人工记录检查两物种较强基线的严格门槛和源CPU身份，防止用伪造PASS字段或错作业放行；失败门槛与缺少单独合同须在读取原生模型／FASTA或创建输出前停止。人工双方向序列检查非ACGT转N、负链半开CDS坐标／phase及GFF往返，且同一GEN／DNA对象确实进入两个表示、新发射gain28进入装配，旧A3发射不能被调用。最后检查零分选null、负分保留自然池、无完成候选为空，候选几何和来源不因A36分数改变。

使用原QueryPointer人工fixture及手工15类输出，不读取训练模型、真实GEN、参考、缓存或DEV，不执行optimizer；人工特征不证明原生GEN的数值或新模型准确率。已有装配／DP／梯度测试不重复扩大。四项检查拟经Slurm运行，generanno／private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，no-requeue、唯一输出M28-A36-STRUCTURED-EMISSION-RAW-CHECKS及M28-PILOT-R5/logs/%x_%j。

## 执行记录

四个新增任务文件已明确传输到canonical，两个Python AST及SBATCH语法通过；本地不具备全部原模块，实际检查只在远程generanno运行，没有为了导入测试而同步整项目或安装依赖。实际五项Slurm调度／预约检查无硬约束，输出无冲突，独立DEV合同尚不存在。

sbatch --parsable sbatch/M28-A36-STRUCTURED-EMISSION-RAW-CHECKS.sbatch已提交一次，返回13408563，请求private／0GPU／2CPU／8GiB／5分钟／no-requeue。实际root／batch／extern全部COMPLETED0:0／0:47，batch MaxRSS742248K；AllocTRES确认CPU2／8G／无GPU，gpu034，STATUS／JOBID相符。四项测试全部PASS／内部0.858秒：两方向规范化、新发射gain28进入最终选择、负链CDS[(0,5,phase2),(10,14,phase0)]和GFF往返一致；两个表示共享同一编码和DNA对象，旧A3发射未被调用，非正和空自然池保持原规则。没有真实模型／GEN／FASTA参考或DEV读取，不构造optimizer；这些是工程结果，不是模型准确率或自然完整注释结果。

原scripts/eval_structure_diagnostic.py:51再次出现既有unclosed file ResourceWarning，四项测试和GFF往返通过；保留原stderr，不将警告隐藏、视为科学FAIL或扩大到旧解析器修复。本短作业已主消费，不加入长期心跳。原fit仍RUNNING，两个后续root为正常Dependency，未读取中间权重或指标，也没有取消任何旧作业。

同一固定监控已于13:08:01.147Z实际完成当前三根定时检查并保持安静；同线程turn01a11153-be4d-7253-bb7a-99624544e381正常结束、实际上下文仍gpt-5.6-luna／low。这是本批次scheduled执行证据，不是终态交接或科学PASS。同一心跳继续只覆盖13406144→13408449→13408451，不加入本人工检查、不扩大数据或改变goal。

## 完整评价工程准备

evaluate_structured_emission.py将复用既有完整common ruler，先要求实际TRAIN三根闭合及原五条件PASS、另行固定DEV合同、对应原生推理终态／JOBID，再核对两物种四方向8690窗、同一原生四块GEN上下文、A3生成器和A36发射各自末004608／原4608曝光及两个组件参数量。它从保存自然记录重建候选→owner→selected，并检查实际GFF坐标／phase／数量；全部闭合后才打开参考，不能只读取成功标记或GFF header即评分。

逐窗对照已闭合A31_A3的输入身份、候选几何、次序和greedy／beam来源，只有分数允许改变。A31_A23改变生成器，作为另一性能对照而非同池机制消融。全部7728完整primary参考及背景分母保持原口径，包括模型不支持或不能生成的层；报告逐物种及合并exact-chain P／R／F1、macro和worst-species F1、CDS碱基指标、结构诊断及free→owner→selected覆盖损失。比较包括ANNEVO／Tiberius／Helixer、A31_A3／A31_A23和既有A3／A6／A／R5_C0，全部要求同一完整ruler。

新增四项人工检查针对新接口：完整网格与两个组件末checkpoint／曝光／共享编码／实际TRAIN source身份；A31_A3同池几何／来源必须一致而A36分数可不同，不能把A23或A6文件当作此基线；负链已保存选择与GFF相位重建不触发参考入口；完整参考／背景约束确实接入全部比较和macro／worst物种统计。使用人工记录和原网格metadata，不读真实权重、GEN、FASTA、参考、DEV结果或cache，不进行optimizer。真实科学评价尚未执行，不新增DEV合同或提前放行。本步只新增评价脚本、四项检查和0GPU／2CPU／8GiB／5分钟／generanno的Slurm短检查入口。

### 完整评价检查的实际执行

三个新增源文件及本记录已逐个明确传输到canonical，两个Python AST及SBATCH语法通过；提交前完成五项实时Slurm调度／预约检查，无输出冲突或实际硬约束。sbatch --parsable sbatch/M28-A36-STRUCTURED-EMISSION-EVAL-CHECKS.sbatch仅提交一次，返回13408892；即刻squeue为RUNNING／gpu034。实际root／batch／extern全部COMPLETED0:0／0:15，batch MaxRSS748992K；AllocTRES为CPU2／8G／无GPU，STATUS／JOBID相符。四项测试全部PASS，内部0.084秒；检查在独立人工记录上确实拒绝错checkpoint／曝光／网格／候选来源／完整分母，允许新发射改变分数，并复原负链非空CDS选择、phase和实际GFF。

原scripts/eval_structure_diagnostic.py:51的unclosed file ResourceWarning再次出现，保留stderr；未把此既有警告解释成科学失败，也未扩大修复旧解析器。当前结果只证明新评价接口的这些工程路径，不证明完整8690窗实际原生推理、模型准确率或独立泛化。本短作业已由主对话消费，不加入长期心跳。

固定监控13:38:00.911Z的取证为fit13406144 COMPLETED0:0、natural13408449 RUNNING／gpu034、CPU13408451等待正常Dependency；原五项TRAIN条件尚未获得整链CPU核验结果。监控对话与同一心跳继续只覆盖这三个源作业，不增加轮次、改阈值、替换checkpoint、打开新DEV或提交真实评价。TRAIN PASS之后才按原约定另行固定原生DEV合同；若FAIL则关闭本模型，不用以上工程PASS替代科学门槛。

## 原生执行入口准备

两个新增SBATCH分别调用已通过检查的infer_structured_emission.py和evaluate_structured_emission.py。保持generanno、no-requeue、原排除节点及唯一%x_%j日志；目标为两个尚不存在的M28-A36-STRUCTURED-EMISSION-DEV-INFER／DEV-EVAL目录，不预建目录或覆盖旧结果。Python成功返回后才写实际JOBID与COMPLETED；失败保留Slurm／stderr及已创建的诊断输出，不自动重提、延长、续训或写成功标记。CPU入口要求显式A36_DEV_INFER_JOB，同一真实ID须用于afterok依赖和--source-job，不能猜测ID或误接TRAIN／旧DEV作业。

准备资源为GPU入口private／1RTX3090／2CPU／16GiB／6:30，CPU入口private／0GPU／2CPU／8GiB／5分钟。实际sacct显示同8690窗A31_A3原生13399434完成3:13:47、batch MaxRSS3101384K；其CPU13399454完成0:43、MaxRSS1414820K。A36额外一个SharedContext前向，按既有A34单模型396窗6:35线性估算额外8690窗约2:24:28，总约5:38:15，再留约15%操作余量至6:30；这不是新A36速度实测或公平benchmark。CPU复用已实际完成的同一评价实现，增加基线身份和对照记录核对；5分钟保留启动及新记录处理余量，不构成更改科学门槛。

本步只准备并传输两个执行入口，不提交真实DEV作业、不创建DEV合同、不读取新数据。实际整链TRAIN条件通过后，主对话才固定完整原生DEV合同、再次核对活调度／维护／存储／输出冲突并提交一次GPU入口；其后CPU绑定返回的真实ID。任一TRAIN条件FAIL仍关闭模型；原固定监控仅覆盖13406144→13408449→13408451，新入口没有job ID，不能加入当前监控或标记为已执行。
