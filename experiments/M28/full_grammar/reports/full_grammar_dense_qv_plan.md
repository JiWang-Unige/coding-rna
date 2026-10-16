# 原 TRAIN 全秩 q/v encoder 探索

2026-10-13（北京时间）核验。M最终头9216／native q/v4608、六raw和共同评价／配对整链完成：pooled完整链F1为72.199933%，比匹配B高1.897309 pp，六染色体均改善；低于单fit低秩E的72.958322%，也仍低于Tiberius pooled和ANNEVO。完整结果按计数核验后记录，不作纯rank因果、独立泛化或论文成功判断。E/M都保持TRAIN-only，冻结B/C/D DEV不变。以下保留固定设计与执行历史。

## 固定科学设置

仍只原Arabidopsis和rice的六条TRAIN染色体及原fixed seed1 fresh4608 draws／4219 unique。起点为原13434603末4608 C0及AdamW状态，不接续已完成B9216，不继承验证中的可丢弃更新。原标签、coarse ε0.1、unknown边缘化、已知行分母、importance1、24576 bp joint编码、nearest-center ownership、FP64完整结构NLL及skyline／cap128不变，不加native CDS先验、缓存或新数据。

仅原b048 GENERanno最后四层22–25的q_proj／v_proj八个无bias权重矩阵改为独立FP32主副本，前向按原输入dtype投为BF16；其余encoder参数冻结、eval模式、原RoPE和SDPA flash不变。实际矩阵几何q为2048×2048、v为256×2048，新增可更新权重18,874,368，连同原366,191参数头为19,240,559。共享原发布资产和原起点文件只读，不覆盖；整体encoder已经适配，不能称整体冻结。

头AdamW恢复原age4608、lr3e−4；原生q/v新组age0、lr1e−5；weight decay0.01、betas0.9／0.999、eps1e−8及共同clip1。全秩矩阵、低秩残差、学习率与数值舍入方式不同，本探索不是单独的rank／容量因果消融。发现阶段仅一个固定fit／seed；性能不足是结果，不以DEV或中间checkpoint调参。

## 当前实施与接入验证

首先运行四项CPU人工检查：原BF16线性输出一致／FP32主副本梯度；修改后的投影与AdamW dtype／age；状态重载及错误身份／shape／dtype／非有限拒绝；八个固定native q/v矩阵的几何、参数量及其他参数冻结。它们不读研究DNA、注释、checkpoint或发布模型。

CPU通过后，GPU只选择B/E相同的每物种第一fresh TRAIN视图，各自从原4608头及原q/v独立做一次可丢弃4609／native1更新。检验原生初始势函数、NLL和角色VJP与冻结parent逐张量相同，八矩阵梯度非零且有限、其他参数无梯度、原moments正确接续；更新状态重新装载后训练／严格raw势与context MAP一致，并检验14 bp短输入双链覆盖。报告实际显存及耗时；不保存checkpoint，不跑全染色体准确率，不读DEV／sealed／Setaria。两个可丢弃更新不是正式训练结果。

独立文件为src/m28/joint_encoder_dense_qv.py、tests/test_m28_joint_encoder_dense_qv.py及preflight_full_grammar_dense_qv.py；复用已通过的joint tokenizer／编码函数与原抽样／目标实现，不改既有LoRA入口、资产或其他实验。CPU输出为outputs/M28-FULL-GRAMMAR-DENSE-QV-CHECKS/<job>，GPU为outputs/M28-FULL-GRAMMAR-DENSE-QV-PREFLIGHT。CPU 2CPU／4GiB／5分钟，无GPU；GPU单RTX3090／2CPU／16GiB／20分钟，沿用generanno原环境、不安装依赖，均no-requeue。提交前读回真实Slurm、限制、维护与输出路径。

## 已完成的接入验证

CPU13509581于05:26:12–05:26:28UTC完成，root／batch／extern均COMPLETED0:0，四项检查全部通过。GPU13509644于05:30:02–05:32:20UTC完成，同样三项终态全部COMPLETED0:0，脚本107.776秒；实际沿用generanno和单3090，无依赖安装。两个原固定TRAIN视图为Arabidopsis NC_003070.9起点14389248负链和rice NC_089046.1起点14462976负链，各读取49152 bp／三个native窗口，各自从原头4608／native0独立恢复至可丢弃头4609／native1。

两者初始native势、绝对NLL及角色VJP与冻结parent逐张量相同；八个q/v矩阵均有有限非零梯度和非零FP32更新，其他encoder参数冻结且无梯度。更新状态恢复后实际训练与raw势逐元素一致、context MAP链和gain一致，14 bp人工输入的双链覆盖通过。最大allocated6,346,762,240 bytes，reserved6,765,412,352 bytes；两视图端到端验证分别45.587和33.033秒，其中包括额外parent／VJP／重载／raw检查，不能直接外推正式单步速度。未保存checkpoint或覆盖原权重／moment，未访问DEV／sealed。这是接入和可行性证据，不是完整链准确率结果。

## 正式读出与停止边界

CPU和原生接入已通过且由主核验，独立fit／核验／六raw／共同评价／配对入口已实现，接下来只运行新增接口的六项人工CPU检查，不重复已通过的四项底层检查。它们验证最终checkpoint身份与八个完整FP32矩阵、实际native装载及前向公式、头／qv双组不同学习率及age／moment登记、raw收据和配对分母／macro／worst计算；不访问研究数据或真实权重。正式设计仍追加4608、只最终头9216／native4608读出，6144／7680仅作日志完整性记录，不选checkpoint。共同评价对齐原全22,688 primary参考及两物种／pooled／macro／worst、CDS-base和背景端点，所有unsupported保留；配对对照为已完成B，不依据正在运行E或新DEV选参数。

正式fit保持单3090／2CPU／16GiB、20小时上限，与既有encoder反传分支资源相同；native预检已证明约6.35GB allocated可行，但其额外复算使耗时不能作为单步估计。首个实际训练进度用于一次速率可行性判断，不改变科学参数或自动延长。严格核验2CPU／4GiB／5分钟；raw六任务最多并发6，各单3090／2CPU／16GiB／2小时；评价2CPU／4GiB／20分钟；配对2CPU／4GiB／10分钟，均no-requeue和原generanno。三次checkpoint各保留FP32 q/v及两个AdamW moment，独立fit产物上限由96MiB改为768MiB以容纳实际全秩状态；它不是改变训练样本或更新数。新输出根和旧C/D/E完全隔离，不写原发布资产或原4608状态。

正式训练已经提交并运行，仍需最终权重与整链结果；不能以资格验证或任务提交代替最终研究结果。只更新原唯一监控对话与同一心跳，记录M fit、verify、六raw、eval、paired和整链完成产物；首个终态一次交接后由主消费。不因单个探索点胜过TRAIN对照改变冻结DEV或作独立泛化／发表claim，不新建监控或自动化、不取消／重提既有C/D/E。原始数据、旧权重、预测和冻结结果保留。

## 最终接口检查和实际提交

2026-10-12 09:00UTC首个长程速率判断：fit13510324仍RUNNING，最新完整日志头step5568／native960；832至960的128步耗时1306.649秒，约10.208198秒／step。按这段速率和尚余3648步推算总fit约13.105707小时，20小时相对余量约6.894293小时，不需要调度或科学参数变更。这只是运行风险估计，不是保证ETA或准确率；E已完成结果不会改变M固定设置。

13510272于06:05:00–06:05:34UTC完成，root／batch／extern均COMPLETED0:0，六项新增接口检查全部通过，检查本体0.561秒；未读取研究DNA／参考／真实checkpoint／发布模型或执行模型更新。底层四项未重复运行。

06:06:56UTC提交前private gpu034／035各8张3090、各已分配2张，原own expanded job32，低于80限制；唯一维护预留是11月6日CPU节点，scratch可用146TiB，五个新输出根不存在。20小时连续fit使用private，既有D/E不动。

| 阶段 | 实际job | 有效依赖 |
| --- | --- | --- |
| 固定追加4608 fit | 13510324 | 已完成13509644／13510272为显式输入和启动资格核验 |
| 最终权重／draw／AdamW核验 | 13510325 | afterok13510324 |
| 原六条TRAIN raw | 13510326_0–5，max6 | afterok13510325 |
| 共同评价 | 13510351 | afterok整个13510326数组；程序核验已完成B13494895及其结果 |
| 对B配对评价 | 13510352 | afterok13510351；已完成B评价13494898仅作为只读输入 |

最初共同评价提交还声明了afterok13494898，Slurm拒绝Job dependency problem、未创建评价job；只读scontrol确认旧B ID已从控制器移除。fit／verify／raw已经成功创建，均未重提；移除过期的Slurm依赖后仅补交评价和paired，程序中B完成状态、输入一致性及科学配对核验保留。一次scontrol逗号列表查询不被支持而返回Invalid job id，随后逐个ID查询确认五个实际job及依赖；这不是作业失败。

06:09UTC实际控制器核对fit RUNNING gpu034、1张3090／2CPU／16GiB、Requeue0，四个下游PENDING／有效Dependency、raw throttle6。started.json保留原source13434603、epsilon0.1、seed1、24576、头lr3e−4／qv lr1e−5与固定9216／4608目标；实际前两更新已写日志，不能用其初始化／编译耗时直接判定全程速率。原唯一心跳已扩展到D/E/M，完整25554字符提示、原target／heartbeat／ACTIVE／15分钟／failed_runs_only／created_at均回读一致；本次新范围的定时追加及M终态回传尚未发生，不宣称已验证。

## 完成后的只读 TRAIN 预测交集诊断

M结果已知后，另做一项明确的探索性计数诊断：比较固定低秩E评价13499823／fit13499817和全秩M评价13510351／fit13510324的原六TRAIN预测，对同一保存primary profile13461637列出共同命中、M新增、M丢失、共同遗漏，并沿用原N的CDS段数1／2–5／≥6和span≤24576／>24576分层。分层规则在本次执行前固定，但诊断选择是post-outcome，不称独立确认。它回答M较E少314个净TP背后的命中交集，不选择模型、合并预测或改变训练设置。

只读取两个完成评价JSON、两份已对同一原DNA直接相等的输入收据、12个既有GFF及已保存references.jsonl，不读FASTA／真实reference GFF／权重，不native forward、不更新或修改预测、不DEV／sealed。复用既有四格／分层计数及GFF重放函数，要求实际hit／unique prediction／duplicates和全部chr／species／pooled重现完成结果，原22688及unsupported全部保留。CPU2／4GiB／10分钟／无GPU／no-requeue，独立输出outputs/M28-FULL-GRAMMAR-DENSE-QV-LORA-RECOVERY；实际job与结果另记于full_grammar_dense_qv_result.md。AST与SBATCH语法已通过，不以静态检查代替真实结果。


## 最终整链与只读交集实际结果

2026-10-12UTC fit13510324于19:21:26完成，verify13510325于19:24:20完成；六实际raw roots13536720／13536753／13536793／13536911／13536912／13510326分别按原seqid顺序完成，最后20:08:36；eval13510351于20:10:36、paired13510352于20:10:58完成，十roots及batch／extern共30条COMPLETED0:0。主核对末9216／native4608、固定draw／参数／moment、六raw双链／身份和共同DNA／参考／计数算术。[完整结果](full_grammar_dense_qv_result.md)保留M对B正差、对E与强基线的负差及背景权衡，不选择6144／7680或改变冻结DEV。

只读post-outcome交集作业13538281于20:25:02–20:25:29UTC完成，27秒，三Slurm步骤COMPLETED0:0；复用原N计数／分层和GFF重放函数，四格边界与实际12个既有GFF／22688参考重放闭合。M／E共同15636，M新增470／丢失784／共同遗漏5798，净少314；段数分区净少42／113／159，9条长spanE命中2／M1。只有探索性reference-recall含义，不是子组F1、rank因果、ensemble或独立泛化；无新训练／推理／DNA／真实reference GFF／checkpoint／DEV／sealed访问。此短任务由主直接消费，不产生新后台监控或伪造watcher通知。

原watcher20:14UTC真实一次M终态交接、handoffYES并保留完整33966字符提示暂停同一心跳；主现已消费M与短CPU13538281，全部等待对象关闭，原监控对话／自动化保留暂停供以后复用，不再轮询旧作业。当前只有研究结果与论文主张未达成，不能把本批成功关闭写成论文goal完成。

