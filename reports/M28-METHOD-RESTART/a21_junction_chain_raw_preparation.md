# A21 raw FASTA 注释与完整评价准备结果

2026-10-05。A21注释接口及完整评价代码各通过四项Baobab CPU合成检查，但随后真实TRAIN门槛3 PASS／2 FAIL，候选关闭，见[实际结果](a21_junction_chain_result.md)。没有执行真实raw FASTA或DEV评价；这份工程准备不放行科学FAIL。拟合、候选、标签、4608步曝光和五项门槛均未改变，见[冻结协议](a21_junction_chain_protocol.md)和[执行记录](a21_junction_chain_execution.md)。

## 接口与执行条件

`scripts/experiments/M28-METHOD-RESTART/infer_junction_chain.py` 复用原四块原生编码、冻结 A3 生成器、greedy＋beam4 完整链几何并集、零 gain 阈值、nearest-owner 及非重叠集合组装，只将评分替换为实际 `JunctionChainHead`。原始序列先转大写，非 ACGT 转 N；两个方向分别编码，尾窗保留真实长度。JSONL 保存全部自然候选、来源和分数，不在记录前丢弃非正候选；另保存外显子 additive、连接 additive、整链 nonadditive 和窗口 null，便于后续解释。

入口在 checkpoint、生成器、原生编码器、FASTA 读取及输出目录创建前，要求 fit、TRAIN-readout、CPU result-check 全部完成且关联同一 fit/readout 作业。它从两物种实际选择数重新计算原五条件，并与 CPU 核验对应，拒绝只声明 PASS、缺物种、非末004608权重或改变曝光。通过后仍须先准备独立 `a21_junction_chain_dev_protocol.md` 才能原生执行；协议文件本身不是人类批准证据。本轮没有生成该 DEV 协议或提交原生注释作业，保留原 PASS 仅支持准备下一协议的边界。

## 四项实际检查及一次修复

`tests/test_m28_junction_chain_raw.py` 使用合成 DNA、原有脚本化 QueryHead fixture 和真正的 H4 `JunctionChainHead`，不是虚构模型成绩。检查覆盖：从选择数重算门槛；门槛失败或 DEV 合同缺失时不读取任何原生输入或生成输出；双链序列规范化、实际连接分量进入 gain、负链半开坐标及 CDS 相位的 GFF3 往返；零 gain 选择空集而仍保留完整自然池，以及负 gain 不被提前删去。

首次命令 `sbatch --parsable sbatch/M28-A21-JUNCTION-CHAIN-RAW-CHECKS.sbatch` 返回13390546，root／batch FAILED1:0／16秒，extern COMPLETED0:0。前三项通过；第四项把 JSON 的列表坐标与 Python 元组直接比较，出现 `[[[0,9]]] != [[(0,9)]]`。坐标值未改变，这是测试表示比较错误，不是模型或科学门槛失败。原目录、SBATCH 和日志保留。

只将期望候选的外显子元组转为 JSON 对应列表，未修改调用器、候选生成、分数或训练。唯一修复命令 `sbatch --parsable sbatch/M28-A21-JUNCTION-CHAIN-RAW-CHECKS-REPAIR-R1.sbatch` 返回13390593；root／batch／extern全部 COMPLETED0:0／14秒，STATUS／JOBID匹配，四项检查全部PASS／0.157秒。两次均为 private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟上限，generanno 环境。代码部署以直接内容比较核对，没有增加哈希或第二套监控。

复用的 GFF 解析器 `scripts/eval_structure_diagnostic.py:51` 对 `open(path)` 发出未显式关闭文件的 ResourceWarning；本次坐标／相位往返比对通过。该既有文件句柄清理问题已记录，本轮未扩大到解析器重构。

## 完整评价接线结果

`scripts/experiments/M28-METHOD-RESTART/evaluate_junction_chain.py` 要求完整四方向scope和8690窗、固定末权重及原生四块上下文；显式FASTA子集不能作为完整DEV结果。它逐条直接比较A6与A21的自然候选几何、顺序和decoder来源，核对实际gain与外显子／连接／整链分量及null，按原规则重建选择，再核对GFF3结构与总数。完成这些步骤后才读取参考注释。

评价接入原7728条完整参考链和背景分母，不改为模型可达子集；输出逐物种与pooled结果、macro F1、背景负担、结构诊断及候选到最终选择的损失。比较源保持A6、原A、R5 C0与已完成的ANNEVO／Tiberius／Helixer原生结果，保留未支持参考层及原生训练曝光限制。原两物种F1严格超过A且背景span不超过C0的pilot门槛仍单独计算；它不等于原生竞争力、泛化或论文成功。A21相对A6同时改变目标和头容量，不把配对候选池解释为容量匹配因果消融。

命令 `sbatch --parsable sbatch/M28-A21-JUNCTION-CHAIN-EVAL-CHECKS.sbatch` 返回13390788，root／batch／extern全部COMPLETED0:0／14秒，STATUS／JOBID匹配，四项检查PASS／0.060秒。private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟上限，generanno环境。检查使用合成记录及现有常量，不读取真实DEV参考或运行新模型：拒绝部分网格及错误编码上下文，拒绝候选来源／坐标变化及不一致分量，完成合成GFF3重建与数量核对，拒绝参考或背景分母改变。复用解析器仍有上述ResourceWarning，往返比对通过；没有新工程重试。

## 未完成的科研工作

尚未装载 A21 实际训练末权重做原生 GPU raw FASTA 推理，未读 DEV 标签、sealed 或 Setaria，没有新基因 F1、泛化或论文结论。合成接线检查不能替代原生编码上下文、完整窗口覆盖和独立评价。后续只有实际 TRAIN 流程经 CPU 核验且原五条件全部通过，才准备单独完整 DEV 协议；FAIL 则按原停止条件关闭，不用新阈值、epoch 或早期权重救参。

固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor已成功交接A21终态并实际暂停；主对话消费后只恢复A17对象，不再检查关闭的A21。这些短CPU合成检查由主对话处理，没有加入心跳或建立重复等待循环。
