# 逐候选空集模型的原生 DEV 注释条件方案

2026-10-07北京时间，在原生DEV访问前固定，当前仅准备且未执行。A38 fit13412630仍运行，自然13412701和CPU13412702等待依赖；尚未产生可信的完整末模型读出或科学判据。只有原TRAIN链闭合并满足全部五条件才允许执行本方案，不能仅凭afterok或本文件存在放行。A37和A36仍科学FAIL关闭，不改变源训练协议、数据隔离或原评价分母。

## 执行前提和模型身份

原三个root须实际COMPLETED0:0、STATUS／JOBID及各summary一致；主核对实际4608更新、三个checkpoint、原1491窗／16321自然候选和四个CPU证明，再按[a38_logistic_context_protocol.md](a38_logistic_context_protocol.md)重算较强A31_A3五条件：正确配对各严格高于1531／478，known_nonmatch各不高于347／251，pooled≥5-exon正确严格高于700。原2916／1086参考窗口分母不变，A37配对不替代门槛。任一FAIL关闭A38，原生DEV不执行，不选另一权重、改阈值或追加曝光。

infer_logistic_context.py绑定fit13412630／natural13412701／CPU13412702，readout检查13412682、preflight13412601及末step004608；在模型、GEN、FASTA或新目录访问前再次要求完整源身份、数值证明和五条件成立。A38的mean-known BCE配置不可被A37集合似然或旧source ID替代。生成器仍为A3末outputs/M28-A3-FIRST-FIT/step_004608.pt、649457参数；评分器为outputs/M28-A38-LOGISTIC-CONTEXT-FIT/step_004608.pt、612324参数，两者原4608更新、eval及无梯度／新增optimizer。

工程检查只用人工future-PASS收据和小型伪输入：A38身份／目标／源ID绑定，较弱A37或缓存true不能绕过原五条件、失败前不加载科研数据或建目录，新评价身份与保存score kind接线。复用已验证的A37raw／完整评价逻辑，只变更新模型／源身份和输出名；不重复未变的坐标、搜索、GFF、DP或完整分母人工测试。人工PASS不是实际TRAIN通过。

## 原始输入和完整注释

拟执行仅原--dev：拟南芥NC_003074.8共23459830bp、水稻NC_089041.1共29936421bp，来自原data/m1_screen物种genome.fa；不物化其他染色体、独立物种、Setaria／sealed或A33例外。沿用24576bp宽、12288bp stride、固定尾窗和两个方向，4scope／8690窗。GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、原四块6144bp、BF16 materialization转FP32、非ACGT转N、方向化及valid-length保持不变。

每方向窗只编码一次GEN，原A3提议和A38评分复用同一特征与DNA，各有自己的上下文。自然提议仍greedy＋beam4／completed4／cap128及原ORF／内含子／相位规则；完整gain为exon additive＋junction additive＋nonadditive−null，FP32、操作点0。旧query／segmentation分数不进入新gain。没有reference、标签、feature cache、参考注入、候选补正或label筛选；全部自然链及来源、非正gain、四分量和null须保存，原owner／几何去重／span DP不变。

唯一输出outputs/M28-A38-LOGISTIC-CONTEXT-DEV-INFER含4个方向JSONL、完整summary及A38.gff3。新天然池须逐窗匹配A31_A3同网格已保存的几何、顺序和greedy／beam4来源；旧93195是待核对的实际对照数，不预填新结果。模型前后权重相同、所有原窗口及最终GFF完整才允许COMPLETED。接口保留显式--fasta用法，但本方案不授权自动换数据，只拟执行原--dev。

## 完整评价和结论边界

evaluate_logistic_context.py首先核对新推理完整网格、源模型／TRAIN身份及原生作业终态，配对实际A31_A3记录，核对保存FP32分量恒等式并重建装配／GFF。全部预测、空池窗口、源记录和最终CDS坐标／相位闭合后，才打开原参考。不能以缺失窗口或生成器不支持的结构缩小评价。

沿用全部7728条primary CDS-chain参考及原物种背景bp、alternative／unsupported层；报告逐物种和pooled exact-chain precision／recall／F1、macro及worst-species F1、CDS碱基／背景错误、长链／结构诊断及free→owner→selected损失。十项比较为A38、A31_A3、A31_A23、A3、A6、A、R5_C0、ANNEVO、Tiberius、Helixer。A31_A3同提议，但容量和训练曝光不同；A31_A23变了生成器，不作同池因果消融。原native基线来源限制保留，不和cached方法作公平原生运行耗时排名。

这仍是单seed、重复两物种DEV，不是独立泛化、发表成功或BCE创新证据；TRAIN candidate/window不当独立基因F1。若A38通过TRAIN后完整DEV仍劣于强基线，保留真实差距而不救阈值／曝光。本方案不放行新物种、Setaria、claim或公开分享。

## 资源和停止条件

准备检查为generanno／private0GPU／2CPU／8G／5分钟／no-requeue，唯一M28-A38-LOGISTIC-CONTEXT-RAW-CHECKS。拟原生GPU为private1RTX3090 24GiB／2CPU／16G／5:30／no-requeue；依据原A31同网格3:13:47，以及A38preflight去冷启动的七窗load／forward／backward共4.861秒、均值0.694秒外推额外约1:40:36，合计约4:54，再留约36分钟操作余量。这个包含反向／缓存I/O、小样本和不同候选负荷的加法不是原生基准、方法速度主张或完成保证；实际TRAIN自然读出闭合后若计时证据改变，可在保持科学命令不变时调整操作请求并记录依据。

CPU完整评价拟private0GPU／2CPU／8G／5分钟／no-requeue，沿用原A31实际43／69秒及原组件核对余量。预计新保存文本／GFF<64MiB，不缓存GEN。实际执行前检查活调度／预约、必要存储和唯一输出；当前不提交原生GPU／评价CPU，不打开实际DEV模型、DNA或参考。

只有主核验五条件全部PASS且当前工程检查闭合后，才提交sbatch/M28-A38-LOGISTIC-CONTEXT-DEV-INFER.sbatch；不能预先把TRAIN CPU的afterok当科学门槛。随后CPU以真实新GPU ID绑定afterok和A38_DEV_INFER_JOB，调用DEV-EVAL.sbatch。SBATCH不预建推理入口要求不存在的目录，源JOBID和完成标记仅在实际成功后写入。工程失败保留状态／日志，不自动重提、延长、续接、取消或覆盖，不将科学负结果当工程错误。

长等待复用原固定监控对话及heartbeat；本次短准备由主消费，不改当前13412630→13412701→13412702监控对象，不创建新监控或新等待循环。论文和独立验证总目标仍未完成。

## 实际接线核验

新增infer_logistic_context.py、evaluate_logistic_context.py、新目标身份测试、本条件方案及3个SBATCH明确同步canonical，7个文件逐字一致，Python AST和bash语法通过；没有修改正在运行的A38 fit／自然读出／CPU脚本，未修改A37的关闭结果或旧源码。

仅提交sbatch --parsable sbatch/M28-A38-LOGISTIC-CONTEXT-RAW-CHECKS.sbatch一次，真实13412872。初次squeue为RUNNING／4秒／gpu034，scontrol实际private0GPU／2CPU／8G／5分钟，Requeue0／Restarts0及canonical Command／SubmitLine相符。最终root／batch／extern均COMPLETED0:0／16秒，batch MaxRSS724392K，STATUS／JOBID匹配；4项针对新身份与条件入口的检查PASS／0.123秒，无旧坐标／DP／GFF测试重复。

检查使用构造的future-PASS收据、旧人工fixture和小型伪输入，只核对mean-known BCE、A38源ID／末权重路径、实际五条件重新计算、失败前无模型／GEN／FASTA／目录访问、新score kind及完整评价身份接线。没有打开科研checkpoint、feature cache、实际GEN、原始DNA或参考，没有科研模型optimizer、原生DEV或性能产物；人工有效收据不证明真实fit或科学PASS。

两原生输出根DEV-INFER／DEV-EVAL实际仍不存在，未提交原生GPU／评价CPU。当前训练与结果链仍由同一心跳监控，13412872已由主消费并排除。五条件和三个真实源根闭合前不提交原生作业；实际全部PASS后拟GPU命令为sbatch --parsable --dependency=afterok:13412702:13412872 sbatch/M28-A38-LOGISTIC-CONTEXT-DEV-INFER.sbatch，再按新GPU真实ID提交对应afterok／A38_DEV_INFER_JOB的CPU。该依赖是执行顺序，不替代科学判断。
