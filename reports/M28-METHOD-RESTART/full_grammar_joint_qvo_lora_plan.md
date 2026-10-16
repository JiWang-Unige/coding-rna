# 原 TRAIN 注意力输出投影适配对照

## 当前接线状态

2026-10-14 04:59 UTC：15个新叶文件与Baobab逐字相等，7份Python AST和7份SBATCH语法通过；四组实际CPU接线检查尚未运行，native资格及正式训练均未提交，也没有本批输出。实时展开本账户运行中作业102个、等待作业1038个，不满足下述提交前并发<80条件，因此不新增提交、不取消其他作业、不放宽条件。接下来仍是先四组CPU检查、再固定两视图native资格、主消费通过后才正式fit；静态检查不能替代这两级接线结果。

2026-10-14（北京时间）。已完成native学习率四格未产生更好的设置，且两种参数化的学习率敏感性方向相反。下一项固定TRAIN模型对照保留原rank4 q/v高lr设置，新增最后四层注意力输出投影o的同rank适配，检验额外输出映射能力是否改善完整链预测；不是继续扫描lr，也不把未测新结构当成性能改善或算法首创。

## 固定干预和来源

原q/v位于0-based encoder层22–25，各rank4／alpha4／scale1，原8投影／16FP32张量／102400参数。新增四个bias-free BF162048×2048 o_proj的rank4 FP32 A/B：A4×2048，B2048×4，每个16384参数，新增65536，总12投影／24张量／167936参数。原q/v的seed0初始化流逐张量保持；o使用独立CPU seed0均匀±1/√2048流，所有B初始0。仍为原BF16 projection加FP32残差舍入至BF16，不merge或改变released权重；所有其余encoder权重冻结。头366191，总训练参数534127。新增容量、原生optimizer新age和共同clip的影响属于整个干预，不能隔离成纯o或参数数量因果效应。

只从原fit13434603末head4608／native0的C0及保存AdamW恢复，不接续已完成E9216或任何可丢弃更新。原fresh4608有序draw／4219 unique、seed1、三pass1536、joint24576／nearest-center ownership、原fine/coarse/unknown ε0.1／已知行分母／importance1全部不变。head lr与新native组lr均固定3e−4；AdamW decay0.01、betas0.9/0.999、eps1e−8、共同norm1不变，逐步记录真实两组lr。唯一正式读出为累计head9216／native4608；6144／7680只为完整记录，不选checkpoint。没有额外loss、role权重、codon／native-CDS prior或阈值。

仅原两物种六TRAIN；结构语法、十五角色、null0、cap128、codon carry、GT/GC-AG、min-intron和FP64 whole-contig MAP均不变。原E/M／四格／mean-pool权重、势和预测只读；不修改8／10／6冻结DEV家族，不读取sealed/held-out/Setaria/新物种/新panel，不公开发布或Git提交推送。

## 接线检查必须回答的问题

先做四组2CPU4GiB5min／无GPU新增检查：实际十二模块／24张量／167936参数与canonical顺序，旧q/v初始A/B逐张量不变且新增o初始输出严格等于原base；四个o的非零B梯度、原base冻结与命名restore；新末步checkpoint／实际loader／raw接口正确恢复24并拒绝旧16或错o形状／age；新optimizer24参数的names／moments形状／head/native age／lr。只人工fixture，不读研究DNA／label／reference／weights，不重复未变的旧LoRA公式、joint tokenize、grammar或paired ruler全套测试。

通过后，一个1×24GiB RTX3090／2CPU16GiB20min可丢弃native预检。固定复用原B/E/M每物种首个fresh draw，两次各独立从head4608/native0恢复；原始released模型与零q/v/o的实际十五角色势、NLL和role VJP逐张量相同。全部12 B梯度有限非零、A梯度初始0，所有released base冻结；真实head/native lr正确，分别到4609／1且全部o B非零更新。命名restore后训练／raw势和whole-context MAP一致，原14bp短输入两方向接入。无checkpoint写出；不继承这两次更新，也不把短输入零预测或接线资格当成全染色体准确率。

## 完整正式端点和资源

native资格由主消费通过后，单fit20h／1RTX3090／2CPU16GiB，保留原fresh顺序和所有4608步；不自动延长、补训或改变科学设置。独立CPU最终记录核验2CPU4GiB5min，读三checkpoint的head/native moments／age／全部24命名张量和完整log；源weights/moments只读。六TRAIN raw数组每项1RTX3090／2CPU16GiB2h／max6，全部双方向native fresh编码，输入仅原FASTA＋正式step009216；恢复24张量、梯度禁用，参考／label／cache／观察模型均不进入推理。

全部raw完成后，原共同评价2CPU4GiB20min先直接比较新格、B和规范原DNA，再访问参考，保留完整22688主参考及unsupported／overlap／long-span、phase、base/background和其他assessable isoform。对B配对2CPU4GiB10min，六chr／两species／pooled／macro／worst完整报告；与已完成原E q/v及ANNEVO/Tiberius/Helixer/Carbon-A的相同参考口径由主使用已有JSON比较，不重新跑基线或根据成绩删样本。负结果和较差物种不省略，既有paired计算复用，新增接线不改变ruler。

均private-teodoro-gpu／generanno／no-requeue，源SBATCH的实际环境覆盖历史cluster_config默认环境差异；GPU短资格若live shared明显更早且实际已知可用24GiB兼容GPU，则只改调度路由。首训与raw路由以提交前实际GRES、并发<80、storage、maintenance和输出不存在为准。fit上限96MiB、各raw128MiB，新增总预算小于2GiB，不保存新全长势。原q/v实测约13–14.5h，以20h覆盖小幅新增反传；时间只为资源计划，不保证或声称速度优势。

独立输出前缀M28-FULL-GRAMMAR-JOINT-QVO-LORA-{TRAIN-CHECKS,PREFLIGHT,FIT,VERIFICATION,TRAIN-RAW,TRAIN-DIAGNOSTIC,PAIRED}。完成条件为实际root／batch／extern终态0:0、全部STATUS/JOBID/required产物及真实source身份匹配；数组由同一sacct行的JobID与JobIDRaw映射逐task核验，不推算numeric child ID。接线/调度明确错误最多一次同科学范围修复、新输出保留失败；OOM/TIMEOUT/非有限科学更新或负准确率不改参数自动重训。

正常长等待只复用原watcher01a0f71f-8845-7432-bd13-bf54c1808fb3和heartbeat coding-rna-r5-light-monitor，gpt-5.6-luna/low／15min／原通知策略。向研究主对话回传必须裸threadId="01a0332e-e0a9-7aa2-90f9-159405701833"、hostId="local"，省略model/thinking，无main或urn前缀。只有未消费完整事件或需处理失败才真实回传一次，主消费后闭合，不能重查旧LR或已closed事件；不新建对话/自动化/cron，不改goal状态或冻结规则。

工程接入与TRAIN新模型成绩分开。即使本对照改善，也不授权独立测试数据或竞争性／普适性／论文发表claim；原Nature Communications及以上完整论文目标保持，仍需更强性能或独立适用性证据、机制与最终稿件。
