# A31 raw FASTA 注释与完整评价工程准备

2026-10-06。A31完整caller和评价接口已实现，各四项Baobab CPU合成检查全部通过。真实1491窗TRAIN读出和CPU科学核验仍独立进行，不开放DEV、sealed或Setaria。这一步减少TRAIN门槛通过后的实现等待，不改变[a31_local_path_caller_protocol.md](a31_local_path_caller_protocol.md)的评分配方或五项条件。工程测试不能替代原生编码、完整caller评价或竞争性结果。

## 完整注释接口

`infer_local_path_caller.py`按既有raw-FASTA管线处理两方向、四块原生GEN编码和真实尾窗长度。非ACGT转N，冻结A3或A23末004608，不载入新的评分头。原greedy＋beam4／completed4／cap128自然完成池按DNA重新生成，使用同一固定局部角色发射差之和计分，无阈值、长度先验或温度调整。所有自然候选和来源，包括非正分数，先存入JSONL；之后按原nearest-owner、重复几何最大分数和span非重叠选择装配并输出GFF3。

原1491窗准备重用已保存池、没有新搜索；本raw接口实际执行时必须从输入DNA重新生成自然池，二者不可混写为同一工程操作。这里仅做合成输入测试，没有加载真实训练权重或原生GLM，没有读真实参考或缓存。合成特征不验证四块原生编码的实际数值。

## 科学执行条件

入口先要求原TRAIN读出与CPU核验完成且源作业对应，再从两物种selected计数和eligible参考重算固定五gate。所选variant须单独通过，另一variant失败不自动否定它，也不被借来替它通过。原分数、自然池、冻结权重和零更新均须与实际核验记录对应。

以上条件在原生权重、编码器、FASTA读取或输出创建前检查。通过后仍须单独固定完整raw-FASTA／DEV协议；协议文件只是科学合同条件，不是人类授权证据。本轮不创建该DEV协议或提交原生注释作业。失败variant按原停止条件关闭，不救阈值或另换checkpoint。

## 四项工程检查

检查所选variant的门槛重算及两物种参考分母；失败门槛或缺少完整DEV合同不能读原生模型／FASTA或创建输出；两方向规范化、负链半开坐标及GFF3 CDS相位往返，固定局部分数确实进入最终装配；零分选择空集、负分不提前删池，以及无完成候选的空池。使用原QueryPointer合成fixture和明确构造的15类局部输出，不是新模型性能。

## 实际检查结果

本机两个新Python文件AST和SBATCH语法检查通过，窄范围上传完成。提交前完成五项实时Slurm检查；private有可用CPU资源，没有相关维护截断、并发限制或目录碰撞。命令sbatch --parsable sbatch/M28-A31-LOCAL-PATH-CALLER-RAW-CHECKS.sbatch返回13398776，private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，generanno，唯一outputs/M28-A31-LOCAL-PATH-CALLER-RAW-CHECKS。

root／batch／extern全部COMPLETED0:0，实际1:26，batch MaxRSS718344K，STATUS／JOBID相符。四项检查全部PASS，内部2.395秒。合成多外显子负链的固定局部分数28实际进入装配；半开坐标CDS[(0,5,phase2),(10,14,phase0)]及GFF3往返一致。零分最终选择为空集，负分仍保留自然候选，无ATG完成池为空。未读取真实GLM、参考、DEV、sealed或Setaria，optimizer0。

复用GFF解析器scripts/eval_structure_diagnostic.py:51再次产生既有unclosed file ResourceWarning；本次往返比对通过。此文件句柄清理问题已在A21工程结果记录，本轮未扩大到解析器重构或为相同警告重跑。

当前主研究监控仍只跟踪A30三根和A31两根；本短工程作业已由主对话消费，不加入心跳或新建监控，不扩大科学数据范围，也未改运行中的TRAIN源文件。

## 完整评价接口

`evaluate_local_path_caller.py`只接受两物种双方向共四个固定scope、8690窗、原生四块GEN上下文及对应末004608的完整DEV推理记录。它先从全部自然候选重建原装配并核对计数、GFF坐标和相位，再读取参考计算结果。评价复用原两物种7728条完整primary参考链及全背景分母，报告exact-chain P／R／F1、CDS碱基指标、结构诊断和候选到装配的覆盖损失；不会将参考分母缩成模型可生成的子集。

A31-A3与A6使用同一冻结A3自然生成器，评价逐窗核对几何、次序和greedy／beam4来源；分数允许不同。A31-A23的生成器不同，候选池不要求相同，只核对输入网格和编码上下文。两者不能混写为同池消融。比较接入既有A3、A6、A、R5-C0和ANNEVO／Tiberius／Helixer，逐项要求同一完整参考及背景分母。原pilot gate仅作辅助报告，不作为A31新的科学继续条件。

此接口尚未执行真实DEV评价。A3对照同时改变完成池和分数，A6有额外评分头训练暴露，A23改变生成器及训练后的共享权重；重复使用的DEV不能证明独立泛化，也不能支持匹配容量的因果结论或论文成功。

## 完整评价工程检查结果

本机两个新Python文件AST及SBATCH语法通过，三个文件窄范围上传完成后，重新完成五项实时Slurm检查并确认唯一输出未存在。命令sbatch --parsable sbatch/M28-A31-LOCAL-PATH-CALLER-EVAL-CHECKS.sbatch返回13398850，private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，generanno，唯一outputs/M28-A31-LOCAL-PATH-CALLER-EVAL-CHECKS。

root／batch／extern全部COMPLETED0:0，实际1:48，batch MaxRSS726992K，STATUS／JOBID相符。四项检查全部PASS，内部7.669秒：拒绝不完整网格、错误上下文或checkpoint；A3同池几何必须一致而A23异池允许；从保存候选重建选择和GFF且不调用参考入口；完整7728参考及背景分母不能被替换。全部使用合成记录，没有加载真实模型、读取真实DEV参考或执行原生编码。既有GFF解析器的unclosed file ResourceWarning再次出现，但往返核对通过，未为相同警告扩大修复。

此短作业已由主对话消费，仍不加入研究心跳。固定监控对话的同批次04:20及04:35定时运行已核验，仍使用gpt-5.6-luna／low，每15分钟检查A30和A31原五根。所选A31 variant的TRAIN五项条件经实际CPU核验通过后，才另行固定并执行完整raw-FASTA／DEV协议；本轮未创建该协议或提交原生推理、评价作业。
