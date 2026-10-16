# 原 TRAIN 双头 CDS 证据支路的匹配续训

2026-10-09。用户已解除原TRAIN内改模型及增训的阶段限制。本轮检验两个冻结原生CDS分类头的逐碱基证据，经28系数角色校准后能否改善完整CDS链恢复；不把预训练分类器本身或一般后处理包装为本项目创新。

## 固定比较和输入

以原13434603末4608步C0为共同起点。新模型保持原366191参数结构，新增零初始化Linear(2,14,bias=False)，将两列六位置CDS减NON_CODING的FP32分数映射到14个非背景角色增量，背景role0不变。冻结原b0483c23b6b63787b61a6d3a204a9b517d6ba345编码器和两个BF16原生头，仍按6144 bp编码，原DNA分支、24576读取窗、nearest-center所有权、完整基因语法和epsilon=0.1均不改。

匹配对照复用已经运行的OBJECTIVE-CONTRAST/STRUCTURED_CONTINUE末9216步，不另训重复控制。新支路与该控制都使用原1536 draws／1491 unique、原顺序和importance、fine优先／coarse／unknown分母及完整基因集合观察NLL；都从同一旧末权重增加3遍、4608次更新，正式末权重为累计9216步。新28系数也是4608次更新，但其AdamW step从0起，不虚构已有4608次暴露。

恢复全部旧C0参数的原AdamW moments及step4608，再追加新系数参数组，AdamW lr 3e-4、weight decay 0.01、betas 0.9/0.999、eps 1e-8，所有366219参数共同clip norm 1。新增参数改变共同裁剪和梯度流，因此这里只匹配样本、目标及更新次数，不声称有效优化幅度完全等同。无辅助损失、head选择、上下文长度选择、阈值调整或DEV选模。

## 接入和实施验证

先完成[八个原TRAIN见证的零支路重放和独立梯度检查](full_grammar_native_cds_prior_preflight.md)。再固定各物种原抽样序列第一个draw，各独立恢复同一旧权重和优化器，执行一个一次性更新：旧参数step4608→4609、新系数step0→1、两个冻结头及编码器无梯度；真实更新的神经权重仅在内存序列化后严格加载，通过原三输入raw neural接口重放相同预测。正式raw入口拒绝4609步一次性权重。这两次更新不保存科学checkpoint，不带入正式续训。

严格raw caller只导入推理模块和固定模型定义，不导入训练计划、参考catalog、标签、观察损失或旧feature cache。使用原raw FASTA读取、BF16冻结编码器、两个原生头、末权重28系数、原owner拼接及exact exon-budget skyline双链整染色体解码。输出仍为完整primary CDS链，不扩张为UTR／所有isoform。

验证成功后才提交新正式续训，1 RTX3090 24GB、2CPU、16GiB、18小时。资源估计沿用原结构续训单遍17323秒／1536 draws，28系数矩阵仅附加小幅计算；以一次性更新实测核对，不把估计称为完成保证。输出独立为outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-FIT；每遍保留工作checkpoint，正式仅末009216，不覆盖旧结果，不自动恢复或取消作业，新增fit文件总上限96MiB。

## 后续判断

末权重经完整draw顺序／数值算术／两组AdamW各自步数核验后，运行原六TRAIN染色体的严格raw FASTA推理，再用原共同评价尺读参考，对照旧4608锚点与STRUCTURED_CONTINUE9216。至少报告AT、Rice和pooled exact-chain precision／recall／F1，CDS base、背景假阳性、重复预测及不支持参考层。训练NLL或见证梯度不是完整链性能。

本支路尚无新DEV执行方案，不进入既有三目标DEV冻结实验，不按已知DEV表现调支路。未放行sealed／held-out／Setaria，无公开发布、Git提交或推送。即使TRAIN改善，也不能直接证明独立泛化、基线竞争力或论文目标完成。

长等待仅复用原监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3及心跳coding-rna-r5-light-monitor，不新建对话或自动化；正常调度状态静默，完成或需处理的终态一次交接。

2026-10-09 10:23 UTC：八见证初始化预检13453019已按实际终态和产物核验通过；随后提交两次独立更新预检13453049，shared-gpu／gpu021、1 RTX3090、2CPU、16GiB、5分钟。本阶段尚未提交正式支路续训或其评价。

## 更新预检实际结果

13453049的Slurm root／batch／extern均COMPLETED 0:0，UTC10:23:11–10:24:03、52秒。两个固定draw分别为AT NC_003070.9|29786112|-和Rice NC_089046.1|26886144|-，原始结构NLL分别22.866312936和0.001662265；加原importance／known行数后为0.001003254和9.01836e-8。裁剪前梯度范数0.474438及1.71780e-5，未据此更换样本或改变学习率。

两个独立模型的全部旧AdamW参数均4608→4609，新支路均0→1。更新后的神经权重在内存序列化、strict加载后，通过真实原生双头三输入raw neural适配器逐元素重放直接四输入模型；正式caller拒绝非末4609权重。两个一次性更新没有保存checkpoint、未带入正式训练，旧源文件不变、冻结编码器及双头无梯度、DEV／sealed未使用。稳态第二draw实测11.68秒，与原结构训练速率相近，保留18小时正式上限。

后续流程采用afterok：正式fit→CPU全部更新／checkpoint核验→六任务raw FASTA数组（max3并行）→整个数组成功后的CPU共同评价。核验2CPU／4GiB／5分钟；raw每任务1RTX3090／2CPU／16GiB／12小时；评价2CPU／4GiB／20分钟。后三类输出分别为M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-VERIFICATION、NATIVE-CDS-PRIOR-TRAIN-RAW和NATIVE-CDS-PRIOR-TRAIN-DIAGNOSTIC；均置于outputs，带M28-FULL-GRAMMAR前缀，不覆盖原预测或默认评价路径。此时只有预检已完成，正式链的任务ID和终态须另记。

## 正式链已启动，尚无性能结果

路径／元数据绑定检查13453172的root／batch／extern均COMPLETED 0:0，UTC10:30:22–10:30:35、13秒，两项合成检查通过；验证新fit配置被末9216 caller接受、旧／非末身份被拒绝，以及原评价函数收到显式的新raw路径／末权重／模型身份。未读取DNA、模型权重或研究指标。

| 阶段 | 实际Slurm任务 | 当前状态 |
| --- | --- | --- |
| 正式支路续训 | 13453177 | UTC10:31:16起RUNNING，gpu035、1 RTX3090 |
| CPU全部更新核验 | 13453180 | PENDING，有效afterok:13453177 |
| 六TRAIN整染色体raw | 13453182_0–5 | PENDING，有效afterok:13453180；max3并行 |
| CPU共同评价 | 13453184 | PENDING，整数组afterok:13453182 |

资源核对时gpu034全部8张GPU已分配、gpu035分配6／8；shared-gpu最长12小时不足本轮18小时合同，故采用可立即分配的private-teodoro-gpu。同型号训练启动后取证了实际资源和TimeLimit；未取消、抢占或改动其他任务。BeeGFS剩余151TB，不因本轮小权重和六份raw输入触发存储阻碍。未将调度开始或先前工程验证称为训练完成。

2026-10-09T10:33:42.090Z已更新原同一ACTIVE heartbeat，绑定原01a0f71f-8845-7432-bd13-bf54c1808fb3，15分钟、gpt-5.6-luna／low设置保留；原三目标、全1491上下文和本支路三条独立cohort均保留，任一需处理的终态交接后暂停，由原主消费并恢复同一心跳继续其余对象。配置精确读回成功；本次更新后的独立首检、真实定时追加和最终自动交接须分别记录，不能用更新前10:33检查代替。按[官方定时任务说明](https://learn.chatgpt.com/docs/automations?surface=app)，本机定时检查仍需要电脑开机、应用运行和连接可用，不保证休眠／断线即时执行。

原监控对话已完成本次三cohort手动首检，实际turn01a1203b-7aaf-7a92-a0cc-7339540aad03，UTC10:35:31–10:35:53，实际gpt-5.6-luna／low。10:35:37调度检查确认三个cohort均未终态：旧结构／混合续训仍RUNNING；局部目标raw继续运行；全1491上下文的六GPU已完成、六CPUreadout运行；本支路fit运行、后三阶段有效Dependency。未读指标、未交接或暂停。此证据确认同一对话首检成功，不等于本次更新后的真实定时追加或最终自动回传；后两者待后续观察。主仍需完成训练末态、完整更新核验、raw输出和共同评价的结果解释，原论文goal未完成。

真实三cohort定时追加现已观察：原对话turn01a12047-531b-7950-b274-dda06b82dc9c，heartbeat wrapper actual current_time_iso=2026-10-09T10:48:27.673Z、automation_id仍为原coding-rna-r5-light-monitor，输入含三cohort全部实际job；turn_context为gpt-5.6-luna／low，10:48:46完成DONT_NOTIFY。10:48:31检查均为正常运行或依赖等待，未读科学指标、未交接或暂停。这与手动首检分别取证；最终终态交接仍待完成。本轮新准备的[支路冻结DEV方案](full_grammar_native_cds_prior_dev_plan.md)及人工检查没有修改正在运行的训练、三目标DEV或监控对象。
