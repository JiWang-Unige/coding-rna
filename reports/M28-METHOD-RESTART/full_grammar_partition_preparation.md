# 完整 CDS 语法配分函数原型

2026-10-07。本原型准备完整语法的可微归一化，为研究现有有限候选训练与完整 DNA 语法推理的支持域差异提供同语法接口；它不是新的拟合、监督损失或性能结果。当前 A44 仍使用原固定有限已知兼容集合目标，本原型不接入其训练、预测或科学门槛。

## 归一化空间与实现

`src/m28/full_grammar_partition.py` 对原完整 CDS 语法允许的所有非重叠 gene 集合求 FP64 log sum exp，包含一次空集合，允许相邻 gene。分数仍是15角色相对 intergenic 的单位增益和，无角色权重、长度归一化、温度或新参数。输入只有定向 DNA、有效长度势和显式 exon cap；没有候选池、参考、标签、beam 或先取 Viterbi winner。

语法沿用 ATG、最终 STOP、GT／GC–AG 与至少5 bp intron，保留 exon 数和实际未完成密码子、coding 不允许 N／内部 STOP、intron 允许 N。Donor 成熟后减去相同 intron-role prefix，仅在相同 exon 数与实际 carry 的未来等价状态间用 logaddexp 累加历史质量；到 acceptor 时加回 prefix。完整 gene 的结束质量逐位置进入此前非重叠集合的 prefix 一次，防止重复计入空集合或先前完成集合。Torch FP64 tensor 保留对输入角色势的梯度，没有逐值 Python detach。

这与原 Viterbi 只留最大分数的递推不同；诊断计数表示被处理的状态／边提案，不表示被枚举的 gene 或历史数量。原型没有性能优化，计算图与缓存的真实窗口资源可用性仍需独立证据，不能从小 DNA 正确外推到染色体。

## 人工验证

七项检查使用既有独立拼接 DNA／密码子穷举，不调用新 DP 或 ORFIndex 来枚举真值；再穷举全部 span 兼容集合，直接按完整链角色计算 logZ 和全15角色梯度。八个短 DNA、cap 1／2／4 覆盖固定整数、非整数 FP64 和非整数 FP32 势；logZ 比较至小数点后10位，FP64 梯度 rtol 1e−9／atol 1e−10，FP32 梯度 rtol／atol 1e−6，区分累加后转回 FP32 的舍入。

具体验证还包括：两条 donor 历史汇入同一 carry 状态时全零 logZ 应为 log3 而不是只保留一个赢家；六个相邻 gene 应有64个兼容集合；无合法完整 gene 时应有可微零和零梯度；独立写明的 spliced 角色数组应产生 posterior 加权梯度；单碱基中间 exon、实际 carry 导致的内部 STOP 与 exon cap 不得合并；共同逐行 logit 平移保持 logZ／梯度，正负大分数保持有限。原 Viterbi 目标仅用于 logZ 的下界核验，不作为配分函数真值。

唯一短入口 `sbatch/M28-FULL-GRAMMAR-PARTITION-CHECKS.sbatch` 使用 generanno／private-teodoro-gpu／0 GPU／2 CPU／8 GiB／2分钟／no-requeue，预计小于1分钟；独立根 `outputs/M28-FULL-GRAMMAR-PARTITION-CHECKS` 和唯一 `M28FULLZCHECK_%j` 日志，输出小于1 MiB。提交前检查实际队列、并发、维护窗口、空根与语法。只一次；明确人工实现或 fixture 错误允许同输入／资源／口径、独立 R1 根一次最小修复；第二次失败、数学矛盾或 TIMEOUT／OOM 停止，不改训练、放宽阈值或延长时限。

## 监督与科学边界

原型没有观察事件配分函数、R4 部分标签条件或 gold likelihood。完整语法 logZ 直接减去不完整注释的 gold 分数，可能把未标注 gene 当负例；未知区域如何边缘化、窗口边界如何建模、分母和损失尺度均不能由本人工检查替研究协议选择。本原型也没有校准生物概率、独立泛化、新颖性或较 ANNEVO／Tiberius／Helixer 优势证据。

实际首次训练仍需其明确科学合同及有效批准，不因人工数学检查通过而自动提交，不替代 A44 原两组五条件或救活原 A36／A42／A43 FAIL。没有打开 TRAIN／DEV／独立物种、研究权重、feature、原始注释或 catalog；原固定70根监控保持不变，本短检查由主处理。

## 首次检查与唯一 fixture 修复

10:25:26 UTC 五项活调度／预约、精确空根／日志和 bash 语法检查通过；10:26:25 UTC 账户展开71项，加本次为72，小于80。两 Python 和 SBATCH 内嵌 Python AST 通过，新四文件与 canonical 逐字一致。实际命令为 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-PARTITION-CHECKS.sbatch`。

13426765实际于10:26:27–10:26:34 UTC运行7秒，private／generanno／gpu035／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，root／batch FAILED1:0、extern COMPLETED0:0，batch MaxRSS469684K。7项运行0.935秒，6PASS／1FAIL；失败为独立手写 spliced fixture 错把完整空间当成 null 加一条链。相同 DNA `ATGCGTTAGTTTAA` 还包含合法单 exon `(0,9)`，不能从归一化空间删除。该用例在错误手写断言之前，实际 logZ 和全梯度对独立穷举的比较已通过，不是原型分区或梯度比较失败。

原失败根、STATUS／JOBID、唯一日志以及当时完整模块／test／SBATCH／协议快照保留。唯一R1仅为手写期望加入原 `(0,9)` 单 exon 的角色、分数和 posterior 梯度，将其空间写成 null 加两条互斥合法链；DNA、所有势、原型算法、七项检查、容差和资源不变，不通过改变语法或遗漏竞争路径消除失败。使用独立 `M28-FULL-GRAMMAR-PARTITION-CHECKS-R1.sbatch`、`outputs/M28-FULL-GRAMMAR-PARTITION-CHECKS-R1` 和 `M28FULLZCHECKR1_%j` 日志；本次消耗唯一允许的 fixture 修复，第二次失败停止。

## 唯一 R1 实际完成

10:29:45 UTC 五项活调度／预约、R1 精确空根与日志、bash 语法通过；原型模块与失败根源码直接 cmp 相同，修正 test／内嵌 Python AST 和四文件 canonical 逐字一致。10:30:01 UTC 账户展开71项，加本次为72，小于80，实际执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-PARTITION-CHECKS-R1.sbatch`。

13426777实际10:30:01–10:30:08 UTC闭合7秒，root／batch／extern均COMPLETED0:0，实际 private／generanno／gpu035／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS455852K；Command／SubmitLine／WorkDir及唯一日志相符。7项人工检查PASS0.292秒，STATUS／JOBID／summary身份、tests7／failures0／errors0、原失败13426765与仅fixture修复字段相符。完成后主核对全部零optimizer、无研究输入／模型／候选搜索／标签／DEV／独立数据和未接入训练／推理标志，原失败根和四份快照保持，未重跑已闭合检查。

当前结论仅为既定短人工 DNA 上的完整语法 logZ 与全角色梯度得到独立穷举支持，包括历史质量累加、兼容集合计数、实际 carry、短 exon 与 cap。没有新增研究训练，尚无实际窗口计算图资源、部分注释观察事件、边界状态、注释性能或创新性结果。原 A44 的完整 TRAIN 解码及其配对读出／数值核验／机制诊断仍由原固定70根监控接续，不能用本次数学准备放行其 raw 或修改 FAIL。

## 固定完整窗口资源探针

数学检查13426777闭合后，只运行一次24,576 bp、cap128的完整前向与输入势反向，不用短窗替代训练长度。DNA为原人工seed0生成49,152 bp的首窗口；势为原seed41首窗口的CPU FP32整数−3至3，各15角色。原prototype算法、单位增益和null0不改；不读研究DNA、标签、模型、GEN或feature，不创建参数／optimizer或新的监督损失。

入口 `probe_full_grammar_partition_window_synthetic.py`／`M28-FULL-GRAMMAR-WINDOW-PROBE.sbatch` 使用 generanno／private／0GPU2CPU16GiB2分钟／no-requeue；2 CPU与16 GiB沿用原配对fit的主机量级，不是GPU训练可行性测试。完整耗时未知，2分钟是该单窗口探针的固定观察上界，不由短DNA秒数推算。唯一根为 `outputs/M28-FULL-GRAMMAR-WINDOW-PROBE`、唯一日志 `M28FULLZWIN_%j`，输出小于1MiB。started记录完整输入身份；只有实际前向结束才写forward记录，反向完成且24,576×15梯度有限、共同row-shift梯度残差≤1e−5时才写completed summary，区分阶段完成与整个探针完成。

若OOM／TIMEOUT或未完成完整前向／反向，保留终态与已有阶段文件，不缩短DNA、不降cap、改势、加预算或重提，不把失败变成小窗PASS。该结果只决定现有逐标量CPU autograd表示能否直接接入，不能否决完整语法训练的科学假设或推断不同计算后端。若成功，也只证明这一个人工窗口的执行和梯度有限，不代表完整真实TRAIN、观察事件正确或性能收益。仅命令提交拒绝且无job时可同命令纠正；已创建job不重试。本独立资源探针不消耗／重开已结束数学fixture的R1，不改变A44的70根监控。

## 完整窗口探针实际超时

10:45:35 UTC 五项活调度／预约、空根与唯一日志、bash语法通过，新Python AST与三文件canonical逐字一致。原型模块与原失败数学根的源码直接cmp相同。10:47:26 UTC账户展开73项，加本次为74，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-WINDOW-PROBE.sbatch`。

13426855于10:47:27 UTC开始，10:49:42 UTC终态TIMEOUT／Reason TimeLimit，root记录ExitCode0:0但不是成功；batch为CANCELLED，extern为COMPLETED0:0。实际private／generanno／gpu035／0GPU2CPU16GiB2分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir与唯一日志相符。Slurm实际elapsed135秒包括终止阶段，batch为137秒，不能改称2分15秒科学预算。STATUS为FAILED exit137，JOBID相符；仅started.json578 bytes和一行“完整前向开始”stdout，forward.json与summary.json均不存在。

实际输入完整24,576 bp／cap128／原seed0、41／15角色FP32，Torch2.5.1；motif为start417／donor3098／acceptor1494／stop1121。只有完整输入建立，不存在完成的logZ、反向或梯度结果。MaxRSS在sacct中缺失，不能将TIMEOUT解释为OOM或给出内存峰值。stderr的时限取消消息后还出现pthread_create与threadpool_create资源错误，证据保留；不把后置消息自行归因为模型错误或改成作业重试理由。

该逐标量CPU autograd原型未在这次固定完整长度与资源观察条件下展示可用性，因此不直接接入训练，不重提、不缩短窗口或降cap。七项小DNA数学证据仍成立，但不代替完整窗口执行；完整语法训练的科学假设、GPU或不同计算表示仍未检验，A44及原历史FAIL不变。

代码审阅发现一个可处理的扫描：每个新的entry／实际carry cache key在legal_ends遍历全部排序事件，再逐一排除首STOP或歧义之后不可能合法的端点。既有 `exon_event_range` 只给必要ORF端点范围，范围内仍必须执行完整ORF检验，并保留原事件顺序及short carry。这是下一节隔离索引版的依据，不是实际加速证据。

## 隔离端点索引版准备

`src/m28/full_grammar_partition_indexed.py` 只在建立合法外显子边缓存时使用必要ORF端点索引。首STOP的终止端点仍保留，携带密码子导致STOP前的短中间exon仍保留；范围内逐一运行原ORF检验。空集合、所有完整非重叠集合、原评分算术顺序、logaddexp、state key、cap及六项原递推计数不变，没有先取winner或删去低分路径。原模块、A44及训练／推理调用者不改。局部AST归一化比较已确认递推差异仅为事件范围和三个整数诊断计数；这不能替代实际运行。

`tests/test_m28_full_grammar_partition_indexed.py` 在原七项人工检查上逐用例要求新旧logZ与全输入梯度逐值完全相同、六项原计数相同，并继承独立拼接DNA穷举真值；另加390 bp人工稠密序列、cap128的新旧logZ／全梯度／计数比较，要求实际检查的端点减少。新增 `indexed_events_inspected` 是运行计数，`original_all_event_loop_iterations` 是按同一cache key数与事件数推算的原循环次数，不是原版实测耗时，不作为加速比。390 bp用例不证明24,576 bp资源可用性。

八项短检查入口 `M28-FULL-GRAMMAR-INDEXED-CHECKS.sbatch` 使用同generanno／private／0GPU2CPU8GiB2分钟／no-requeue；独立根 `outputs/M28-FULL-GRAMMAR-INDEXED-CHECKS` 与 `M28FULLZIDXCHECK_%j` 日志，输出小于1MiB。runner要求实际tests8／零失败，避免空suite通过。只运行一次；若为明确实现或人工fixture错误，只允许保持输入、容差、资源和完整空间的一次最小R1；数学矛盾、TIMEOUT／OOM或第二次失败停止，不借用此检查改变科学合同。

仅八项实际闭合后才允许一次新索引表示的完整窗口资源探针：仍为原24,576 bp／cap128／seed0、41／0GPU2CPU16GiB2分钟，不缩短输入、降低cap或增加资源，另用独立文件与输出根。它不是13426855原逐标量全事件扫描的重提；原超时记录保持关闭。完整前向／反向阶段标记和梯度条件不变；若仍超时或OOM，保留结果、不重试。此准备不打开研究输入，不定义部分注释观察事件或gold likelihood，不启动新拟合，不改变A44原70根监控。

## 索引版人工检查实际完成

本次五项活调度／预约通过，原模块／ORF端点索引／原测试与canonical逐字一致；新source、test、SBATCH和协议四文件逐字一致、AST及bash语法通过。11:14:04 UTC账户展开14项，加本次为15，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-INDEXED-CHECKS.sbatch`。

13427003实际11:14:05–11:14:14 UTC闭合9秒，root／batch／extern均COMPLETED0:0，private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS501036K。八项PASS3.115秒；STATUS／JOBID／summary身份、原数学13426777、新表示、tests8／failures0／errors0及所有无研究数据、无optimizer、无调用者改变字段核对相符。原型与索引版在实际人工用例上logZ、全梯度和六项原计数完全相同；稠密cap128人工用例确实减少了端点检查，尚无完整窗口耗时或内存收益结果。

据此准备唯一完整窗口入口 `probe_indexed_full_grammar_partition_window_synthetic.py`／`M28-FULL-GRAMMAR-INDEXED-WINDOW-PROBE.sbatch`。执行前要求13427003实际八项闭合；输入、势、cap、FP64算术、完整前向／反向记录和有限梯度／row-shift残差≤1e−5条件沿用原探针，仅调用隔离索引版并写明表示与原超时job身份。独立根 `outputs/M28-FULL-GRAMMAR-INDEXED-WINDOW-PROBE`、`M28FULLZIDXWIN_%j` 日志，资源仍为0GPU2CPU16GiB2分钟／no-requeue；没有改动原探针或研究训练。

## 索引版完整窗口实际超时

11:16:01 UTC五项活调度／预约通过；新探针、SBATCH及报告三文件与canonical逐字一致，Python AST／bash语法、精确空根与日志通过。11:16:37 UTC账户展开10项，加本次为11，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-INDEXED-WINDOW-PROBE.sbatch`。

13427015实际11:16:37–11:18:43 UTC运行，root为TIMEOUT、scontrol Reason为TimeLimit，batch为CANCELLED、extern为COMPLETED。sacct root ExitCode0:0、scontrol root0:15、batch0:15均如实记录；root的0:0不表示成功。实际private／generanno／gpu034／0GPU2CPU16GiB2分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir与唯一日志相符。root elapsed126秒、batch127秒含时限终止阶段，固定观察上界仍为2分钟。

11:19:36 UTC核验时仅JOBID、started.json677 bytes和原STATUS；started完整记录24,576 bp／cap128／seed0、41／CPU FP32势／FP64递推／Torch2.5.1与原相同motif计数。forward.json和summary.json均不存在，没有完成的logZ、反向或梯度。stderr只记录时限取消；batch MaxRSS15931056K，约15.2 GiB，但无OOM终态。原STATUS仍为RUNNING，表明时限终止时shell退出trap未刷新文件，不能据该旧标记当作继续运行；调度终态优先，原文件与日志保留、不覆盖为成功或重提。

结论是必要端点索引虽在八项人工用例保持完整配分、全梯度与六项原计数，却仍未在固定完整窗口／资源条件下完成前向。它不能直接接入训练，不缩短窗口、降cap、改势、增加资源或重试，也不能据此否定完整语法训练的科学假设。逐标量autograd计算图的存储和逐状态运算是下一种计算表示需要处理的具体问题；下一步可准备保持全部路径空间的批量运算／显式反向表示，重新核对logZ、全梯度及状态语义后才考虑新的完整窗口探针，当前没有实现或宣称其可用性。

原A44完整70根流程继续由同一固定监控处理；本轮没有新拟合、研究输入访问、raw解码或科学性能结果，没有改变历史FAIL或放行DEV／独立数据。论文目标仍未实现，工程准备不替代严格ab initio全gene结果、独立适用范围与机制证据。

## 紧凑显式反向准备

`full_grammar_partition_reverse.py` 用此语法专用的紧凑数值节点替代逐提案Torch计算图。节点只有null、完整外显子增益、带符号的intron prefix与两历史logaddexp；同一count／实际carry状态仍累加全部历史。合法边仍经必要端点索引和完整ORF检验，评分复用已验证的FP64紧凑margin类；原六项逻辑计数不变。达到cap后不创建本来就不能进入结果的非terminal提案节点，但仍保留其原逻辑计数，没有剪枝有效路径或改变支持域。

反向按节点逆序将logaddexp posterior分到两个父节点，将外显子和intron-prefix伴随量累加，再对角色prefix／按phase的coding prefix作反向cumsum与端点角色替换，最后恢复相对intergenic增益对原15列logits的导数。Torch自定义autograd节点保持输入连接，原logZ仍为FP64；上游梯度先在FP64相乘，最后按原FP32／FP64输入类型返回。这一实现只支持CPU的一阶梯度，不宣称GPU、二阶梯度或新监督损失可用。

`reverse_tape_numeric_bytes` 仅统计实际节点数组的数值元素字节，不包含array预留容量、合法边描述、缓存／dict、评分prefix、反向伴随数组或进程基线，不能当作峰值内存。完整窗口可用性仍须完整前向与反向实际闭合。

九项人工检查包含原七项独立拼接DNA穷举／兼容集合真值，额外比较390 bp／cap128的新旧完整logZ、全梯度和六项计数，并核对FP32／FP64下−0.7、2.3非单位上游梯度与真实autograd连接。新旧logZ比较至小数点后10位；梯度沿用原独立穷举的FP64 rtol1e−9／atol1e−10、FP32 rtol／atol1e−6，不要求不同反向累加顺序位级相同，也不根据运行结果放宽容差。独立手写spliced角色posterior梯度、短exon／carry／STOP、cap、null及逐行logit平移等原断言仍对新结果执行。

唯一入口 `M28-FULL-GRAMMAR-REVERSE-CHECKS.sbatch` 使用generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立根 `outputs/M28-FULL-GRAMMAR-REVERSE-CHECKS`、`M28FULLZREVCHK_%j`日志，runner要求实际tests9／零失败，输出小于1MiB。明确实现或人工fixture错误只允许同输入／容差／资源的一次最小R1；数学矛盾、TIMEOUT／OOM或第二次失败停止。只有九项实际闭合后才考虑同24,576 bp／cap128／seed0、41／0GPU2CPU16GiB2分钟的唯一新表示资源探针，不重提两个已关闭超时job。本准备不读研究输入、模型、标签、catalog、DEV或独立物种，不接入任何训练／推理调用者。

## 显式反向人工检查实际完成

11:32:07 UTC五项活调度／预约通过，两source／test及内嵌Python AST、SBATCH语法、新四文件canonical逐字一致、精确空根与日志核对完成。11:33:20 UTC账户展开8项，加本次为9，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-REVERSE-CHECKS.sbatch`。

13427085实际11:33:21–11:33:28 UTC闭合7秒，root／batch／extern均COMPLETED0:0，private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS516172K。九项PASS1.499秒，STATUS／JOBID／summary身份、tests9／failures0／errors0、原数学与索引数学source job、仅CPU一阶以及全部无研究数据／optimizer／调用者修改字段相符。显式反向在这些人工用例上得到原独立穷举与Torch递推的logZ／全梯度支持，包括非单位grad_output；不将短检查秒数推算为完整窗口可用性。

据此准备唯一 `probe_reverse_full_grammar_partition_window_synthetic.py`／`M28-FULL-GRAMMAR-REVERSE-WINDOW-PROBE.sbatch`，执行前要求13427085实际九项闭合，独立根 `outputs/M28-FULL-GRAMMAR-REVERSE-WINDOW-PROBE` 与 `M28FULLZREVWIN_%j` 日志。仍用原24,576 bp／cap128／seed0、41／FP32势／FP64分数及原完整前向／反向标记、有限梯度与row-shift残差≤1e−5条件，0GPU2CPU16GiB2分钟／no-requeue，输出小于1MiB；无研究模型／数据或optimizer。TIMEOUT／OOM／未完成时保留结果、不得缩短输入、降cap、加资源或重试。A44读出13426152的真实TIMEOUT独立按其冻结协议处理，本人工检查不替代A44的科学比较或放行raw。

## 显式反向完整窗口前向完成但反向超时

11:36:59 UTC五项活调度／预约通过；探针／SBATCH／报告三文件canonical逐字一致，Python AST、bash语法、精确空根与日志核对完成。11:38:02 UTC账户展开8项，加本次为9，小于80，提交 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-REVERSE-WINDOW-PROBE.sbatch`，返回13427113；实际11:38:03 UTC进入Slurm，11:38:04 UTC开始。private／generanno／gpu034／0GPU2CPU16GiB2分钟／no-requeue／Restarts0，实际Command／SubmitLine／WorkDir及唯一日志相符。

完整24,576 bp／cap128／原seed0、41前向实际82.341525秒完成，11:39:32.711515 UTC写入forward.json504 bytes，FP64 logZ10992.597775197577。六项原计数为合法exon边276003、状态exon扩展32010632、donor提案28297444、成熟状态提案1011974、terminal提案3506456、最大active状态2667；专用tape为68298643个节点，数值元素1707466075 bytes。前向自报process峰值2161980 KiB，不能将仅节点字节当作进程峰值。

13427113仍在11:40:13 UTC以TIMEOUT结束，Reason TimeLimit；root sacct ExitCode0:0但scontrol0:15，batch CANCELLED0:15、extern COMPLETED0:0。root elapsed129秒、batch131秒含时限终止阶段，固定上界没有延长；batch MaxRSS2489312K，约2.4 GiB，stderr仅时限取消，无OOM终态。11:40:45 UTC核验时已有started.json745 bytes和forward.json，但summary.json不存在，反向未完整闭合，没有完整窗口梯度有限／row-shift检查结果。STATUS仍为旧RUNNING，调度终态优先；所有原文件与失败日志保留，不重提。

当前证据首次证明此人工完整窗口的前向可以在原观察条件内完成，且记录的进程峰值低于原逐标量Torch索引表示的记录；两次均未完成相同完整前向与反向，不给出公平全流程加速比或完整训练内存收益。九项小DNA一阶梯度证据仍成立，但不能代替完整窗口反向。下一计算表示应处理6800万节点逆序循环的执行开销，再做同路径／logZ／全梯度核对；不缩短窗口、降cap或继续重复已超时的纯Python表示。

A44读出超时及其未运行的数值／机制下游另见配对fit协议的实际终态节；这两条证据线均无科研性能或独立泛化PASS，未启动新拟合或实际raw／DEV评价。论文目标保持原完整范围。

## 顺序编译反向准备

generanno实际元数据为Python3.10.18、Numba0.62.1、llvmlite0.45.1、NumPy2.2.6、Torch2.5.1；仅用标准库读取安装元数据，没有在登录节点导入这些计算库或安装依赖。新 `full_grammar_partition_compiled_reverse.py` 直接继承已验证的原前向，自定义backward仅将专用tape的逆序节点累加循环交给nopython；关闭fastmath、parallel和磁盘编译cache，不改方向／节点顺序或使用并行归约。角色prefix／phase展开和上游梯度相乘仍为原表达式，支持范围仍仅CPU FP32／FP64的一阶梯度。

静态AST已确认编译循环除array字段映射为同一数值view名称外与原循环完全相同，子类没有覆盖原forward。原显式反向模块、原Torch／索引版和所有超时根不改；不调用或重开旧A42编译Viterbi，不接入A44或任何训练／推理调用者。静态比较不能替代真实编译和梯度验证。

十项人工检查包含原九项独立拼接DNA穷举／原Torch logZ及全梯度／cap128／非单位grad_output；额外要求与纯Python显式反向在FP32／FP64、−0.7／2.3上游梯度下完整forward和全部八项tape／原计数完全相同、全输入梯度保持原容差，且实际nopython signature非空，禁止隐式Python回退。logZ仍比较至小数点后10位，梯度仍为原FP64 rtol1e−9／atol1e−10、FP32 rtol／atol1e−6，不依结果放宽。

唯一 `M28-FULL-GRAMMAR-COMPILED-REVERSE-CHECKS.sbatch` 用generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立根 `outputs/M28-FULL-GRAMMAR-COMPILED-REVERSE-CHECKS`、`M28FULLZCREVCHK_%j`日志，runner要求实际tests10／零失败，输出小于1MiB。明确实现或fixture错误只允许同输入／容差／资源的一次最小R1；数学矛盾、TIMEOUT／OOM或第二次失败停止。只有十项实际闭合后才准备同24,576 bp／cap128／seed0、41／0GPU2CPU16GiB2分钟的唯一新表示完整窗口探针；包含新进程冷编译成本，不借短检查进程的warm状态，不重提三个已关闭超时表示。本准备不读取研究输入、模型、标签、catalog、DEV或独立物种，不定义观察事件／新loss，不修改A44已停止的执行链和固定监控。

## 顺序编译反向人工检查实际完成

11:55:58 UTC五项活调度／预约通过，新source／test／内嵌Python AST、SBATCH语法、原显式模块与canonical相同、新四文件逐字一致、精确空根与日志核对完成。11:58:23 UTC账户展开8项，加本次为9，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-COMPILED-REVERSE-CHECKS.sbatch`。

13427223实际11:58:25–11:58:45 UTC闭合20秒，root／batch／extern均COMPLETED0:0，private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，batch MaxRSS688912K。十项PASS5.946秒，STATUS／JOBID／summary身份、tests10／failures0／errors0、原三个数学source job、仅CPU一阶及全部无研究数据／optimizer／调用者修改字段相符，实际nopython signature数1。原独立穷举／Torch／纯Python反向的全部规定比较通过，没有放宽容差或重试。

据此准备唯一 `probe_compiled_reverse_full_grammar_partition_window_synthetic.py`／`M28-FULL-GRAMMAR-COMPILED-REVERSE-WINDOW-PROBE.sbatch`，独立根 `outputs/M28-FULL-GRAMMAR-COMPILED-REVERSE-WINDOW-PROBE` 和 `M28FULLZCREVWIN_%j`日志；仍为原24,576 bp／cap128／seed0、41／FP32势／FP64分数／0GPU2CPU16GiB2分钟／no-requeue，输出小于1MiB。要求13427223实际十项闭合；新进程初始nopython signature数必须0，反向计时包含首次JIT，不在计时前预热。

新完整前向还必须与13427113相同完整人工输入的已有forward.json逐值核对FP64 logZ与全部八项原递推／tape计数，差异即停止。原完整反向完成、全输入FP32梯度有限及row-shift残差≤1e−5条件不改；只有这些条件实际完成才写summary。若TIMEOUT／OOM／未完成则保留已有阶段证据、不缩短窗口、降cap、改势、加资源或重试。该探针不重开纯Python反向失败、不提交新训练，也不放行A44未完成的科学门槛。

完整窗口另核对unit-role logZ导数的数学边界：非intergenic列导数为各角色占用概率，应非负；同一位置所有非intergenic角色之和为覆盖概率，不能超过1。用原FP32人工梯度容差1e−6吸收舍入，任一违反就不写成功；这用于检测大图反向的角色／区间累加错误，不是新的训练或预测阈值，也不说明生物概率校准。

## 顺序编译完整窗口前向反向实际完成

12:06:44 UTC五项活调度／预约通过，探针／SBATCH／报告三文件与canonical逐字一致，Python AST、bash语法、精确空输出根与唯一日志核对完成。12:07:45 UTC账户展开8项，加本次为9，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-COMPILED-REVERSE-WINDOW-PROBE.sbatch`，返回13427267。首次状态读取发生SSH握手超时，未据此判定作业失败或重新提交；随后只读取同一job及其产物。

13427267实际12:07:46–12:09:31 UTC闭合105秒，root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU16GiB2分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir及唯一日志相符。batch MaxRSS2869492K，约2.7 GiB。完整输出根有JOBID13427267、STATUS COMPLETED、started.json834 bytes、forward.json504 bytes及summary.json2513 bytes，无缺失；stderr末尾无输出。

输入仍是原24,576 bp／cap128／DNA seed0生成49,152 bp后取首窗口／potential seed41／全部15列FP32整数−3至3／FP64分数／unit role／null0，motif计数与前次相同。完整前向81.784235秒、logZ10992.597775197577，八项计数均与13427113已完成的原前向完全相同：合法exon边276003、状态exon扩展32010632、donor提案28297444、成熟状态提案1011974、terminal提案3506456、最大active状态2667、tape节点68298643、tape数值元素1707466075 bytes。forward自报process峰值2230152 KiB。

完整反向7.302722秒，包含本进程首次nopython编译；signature数实际从0变为1，没有提前预热、磁盘cache、fastmath或并行归约。全24,576×15输入梯度FP32且有限；row-shift梯度最大绝对值1.1920928955078125e−7，小于原1e−5界。最小非intergenic角色梯度−2.9388516360057126e−12，最大同位置非intergenic占用和1.0000001192092896，均在原FP32 1e−6舍入容差内。summary自报process峰值2812872 KiB，与Slurm batch MaxRSS分别记录，不将tape字节或前向峰值当作完整进程峰值。

summary实际确认完整前向／反向、全梯度有限、角色posterior边界及完整前向与原计数一致；optimizer_steps_this_run为0，model_or_backbone_forward、research_inputs_opened、labels_or_catalog_opened、DEV_used、independent_species_used、observed_mask_or_gold_likelihood_implemented、training_or_inference_caller_changed均为false。十项数学检查13427223与此完整窗口资源探针均为实际完成，不是只提交或阶段标记。

结论是这一人工完整窗口的无条件完整语法配分及CPU一阶反向首次在原资源观察条件内闭合。三个旧超时表示和原失败产物保持关闭、未重提；旧表示没有相同全流程完成计时，因此不给公平加速比。单个人工窗口不证明研究输入、多窗口拟合、GPU训练可行性或模型性能；部分注释观察事件、窗口边界与损失单位仍未定义／验证，不将无条件logZ直接当作R4监督似然。

A44原读出TIMEOUT及数值／机制下游未执行的状态不变，仍遵守原停止条件；科学比较为UNEXECUTED而不是FAIL或PASS。没有新训练、实际raw／DEV或独立物种访问，也没有新的论文主张或对外发布。下一步应解决与真实标签语义一致的观察事件／边界和有效监督似然，而非继续重复已闭合的人工资源检查；论文目标尚未实现。

## 部分角色观察似然准备

`full_grammar_observed_partition.py` 定义同一完整非重叠gene集合语法上的观察事件：每个有效DNA位置给定允许角色集合，gene集合的逐位置角色必须全部属于该集合。无条件空间、原FP64评分、实际carry、cap、完整起止与GT／GC–AG条件不变；观察空间只排除与已知角色不相容的集合，不根据分数、候选池或Viterbi winner删路径。观察logZ由同一紧凑tape和已验证的顺序编译一阶CPU反向计算。

原 `labels.py` 的实际策略是fine已知时用单个角色；fine未知但region已知时用原REGION_COLUMNS；两者未知时不监督。新适配函数严格按此优先级给允许角色，未知行保留全部15类，忽略的padding只裁去、不成为背景。它不从positive chains重建负例、不合并同工型、不改变原build_targets。本轮只用人工target验证这个映射，不读取原R4 catalog／标签或研究DNA。

观察事件必须约束完整gene span里的intron和gene之间的intergenic，不只是exon。递推在穿过禁止背景的位置时清除此前已完成集合的背景延伸质量；新完成到当前位置的gene质量仍能进入。intron历史在成熟时检查donor至成熟位置，活跃历史经过禁止该intron角色的位置时失效，避免不同donor历史合并后保留不相容路径。exon允许性使用原角色覆盖优先级的整数prefix计数，保留first三碱基、stop覆盖、单碱基middle和三碱基terminal的实际规则。

未归一化观察负对数概率为 `logZ_all − logZ_observed`，不是完整logZ减不完整gold分数。全部未知时两者应完全相同且loss／梯度为0；部分未知区内所有与观察相容的合法额外gene仍被边缘化。没有选择损失分母、fit系数、窗口importance用法或新训练配置；单纯映射原角色标签不等于已经选定可用的真实R4训练目标。

边界保持当前完整gene语法：若允许角色与任何完整集合都不相容，明确抛出InfeasibleObservation；不自动删观察、添加截断gene、改背景或跳过窗口。原R4包含跨窗／不完整gene的粗粒度监督，因此实际观察支持是否可行及边界条件仍需科学合同解决。此实现不声称已解决截断边界、部分标签偏差、生物概率校准或真实TRAIN可用性。

十项人工检查用独立拼接DNA穷举全部完整chain及非重叠子集，对观察logZ、全部输入角色梯度和NLL梯度作比较；保持原logZ小数点后10位、FP64 rtol1e−9／atol1e−10、FP32 rtol／atol1e−6及−0.7／2.3上游梯度。检查还要求全部未知时原前向／八项计数完全相同与零loss／梯度、六个邻接gene的部分观察保留16个集合、显式完整spliced角色只留一个集合、fine／coarse／unknown／padding原优先级、gene前后及间隙观察、成熟前与成熟后intron历史过滤、短exon／实际carry／cap和不可行边界；实际nopython signature必须非空。

唯一入口 `M28-FULL-GRAMMAR-OBSERVED-CHECKS.sbatch` 使用generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立根 `outputs/M28-FULL-GRAMMAR-OBSERVED-CHECKS` 与 `M28FULLOBSCHK_%j` 日志，runner要求实际tests10／零失败，输出小于1MiB。明确实现或人工fixture错误只允许同输入／容差／资源的一次最小R1；数学矛盾、TIMEOUT／OOM或第二次失败停止。失败不触发新fit或修改旧终态；通过也只证明人工观察事件数学，不放行真实R4、DEV、独立物种或原A44下游。原无条件模块、十项原检查和已完成完整窗口不重复运行／修改；所有训练／推理调用者保持不变。

## 部分角色观察似然人工检查实际完成

五项活调度／预约通过；新四文件与canonical逐字一致，原无条件source、compiled reverse、labels、角色渲染与margin source保持读取时的原内容，Python及内嵌runner AST、bash语法、十项实际test方法数、精确空输出根与唯一日志核对完成。12:30:44 UTC账户展开8项，加本次为9，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-OBSERVED-CHECKS.sbatch`，返回13427476。首次调度读取已为COMPLETING，随后按同一ID读取终态，没有重提。

13427476实际12:30:45–12:30:59 UTC闭合14秒，root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir及唯一日志相符，batch MaxRSS644944K。十项PASS4.569秒；JOBID、STATUS COMPLETED和summary.json1306 bytes齐全，summary身份、tests10／failures0／errors0、原数学13427223／完整窗口13427267、nopython signature数1及人工范围字段相符。

实际通过了独立穷举的观察logZ、FP32／FP64全输入梯度与非单位上游NLL梯度比较，全部未知时原logZ与八项计数完全相同、loss／梯度严格为0，部分未知六gene用例保留16个相容集合，完整spliced观察仅保留一个集合。背景前后／间隙、成熟前与成熟后intron历史筛选、short exon／carry／cap、fine优先于coarse、unknown／padding及不可行边界均得到实际人工检查支持。没有放宽容差、改变fixture或使用R1。

观察似然算子现已实现并得到人工数学支持；summary中observed_likelihood_implemented为true，loss_normalizer_or_fit_coefficients_selected、actual_research_inputs_opened、original_R4_targets_or_catalog_opened、model_or_backbone_forward、research_candidate_search、DEV_used、independent_species_used、training_or_inference_caller_changed均为false，optimizer_steps_this_run为0。此结果不是实际R4可行性、多窗口NLL资源、新拟合、注释提升或创新性证据。

下一步应在保持原TRAIN隔离和标签生成器不变的前提下，准备真实R4观察支持及跨窗／不完整gene的零更新诊断，先取得哪些实际标签与完整gene语法不相容的证据，再定边界和监督单位；不能自动删标签、跳过不可行窗口或据本次PASS放行新首训。A44原读出TIMEOUT／科学门槛UNEXECUTED及固定监控PAUSED均不变；论文目标仍未实现。

## 原R4八窗观察支持零更新诊断

本诊断检查原R4角色观察是否属于当前完整gene语法的可行事件，不比较模型性能。原监督接线13425772和观察数学13427476的Slurm root／batch／extern均已核对为COMPLETED0:0；复用13425772实际summary中的原八个事前固定TRAIN证据窗及顺序，每物种四类，不重新选择容易通过的窗口，也不以八窗充当全部1491窗的比例估计。

新 `probe_full_grammar_R4_observation_support.py` 调用原run_query的load_train_inputs及validate_cached：加载原1536draw／1491唯一窗，混合catalog在JSON解析前按TRAIN筛选，沿原六条染色体allowlist连接窗口／catalog／缓存。只打开八份原定向DNA＋GEN feature缓存，feature字段仅由原validate_cached检查形状，不用于前向；不读旧缓存标签，按原build_targets重新生成fine／coarse／unknown／N／padding目标。八窗实际fine与grouped-region位置数必须与原summary相符，任一漂移停止；原R4标签生成器、窗口、draw及importance不改。

全部有效位置的十五列势固定为CPU FP32零，完整评分仍FP64／unit role／null0／cap128。只计算原观察事件的logZ：若可行，数值是相容完整非重叠gene集合数的对数；若抛出已定义的InfeasibleObservation，记录不可行并继续其余固定证据窗。这是诊断预期结果，不将不可行窗口删出未来训练或改变观察，也不把工程COMPLETED等同于全部观察可行。其他来源、几何、数值或实现异常立即停止；不计算模型NLL、不做反向／optimizer，不写预测或checkpoint。

每窗保留原左右边缘允许角色、跨左／右窗和不完整primary transcript的原几何、原chain未知区及实际fine／coarse／完全未知数。完整gene集合首位置只能是背景0或起始13，末位置只能是背景0或终止14；原观察排除这些两类角色时，单独标记充分的边缘不可行原因。未触发此必要条件不代表边界或整个观察可行；跨窗／不完整metadata与不可行并存也不单独证明完整因果。

唯一入口 `M28-FULL-GRAMMAR-R4-SUPPORT-PROBE.sbatch` 使用generanno／private／0GPU2CPU16GiB15分钟／no-requeue，独立根 `outputs/M28-FULL-GRAMMAR-R4-SUPPORT-PROBE`、`M28FULLR4SUPPORT_%j` 日志。15分钟是八个实际窗口新诊断的固定观察上限，参考已完成单个人工完整前向81.8秒及八次缓存读取留出余量，不是保证或未来拟合投影；只跑观察配分，不重复原无条件完整窗口。产物started、逐窗即刷新的windows.jsonl及仅八窗全部检查后写出的summary，小于1MiB。

只执行一次，TIMEOUT／OOM、真实来源矛盾或实现异常保留结果停止，不自动重提、加时或修改数据；不可行事件本身保留为诊断结果，不触发修复／新fit。此零更新研究诊断沿现有原TRAIN权限，不放行raw、DEV、独立物种、Setaria、A33或对外发布，不重开A44原TIMEOUT读出及下游。边界模型、损失单位、系数和首训合同仍须在真实证据后另行确定；论文目标保持严格raw-FASTA完整gene性能或独立适用范围及机制证据的原范围。

## 原R4八窗支持实际完成及确定边界冲突

12:43:20 UTC五项活调度／预约通过，新三文件与canonical逐字一致，原观察算子／labels／run_query保持原内容，Python AST、bash语法、精确空输出根与唯一日志核对完成。12:44:23 UTC账户展开8项，加本次为9，小于80，执行 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-R4-SUPPORT-PROBE.sbatch`，返回13427583。实际12:44:24–12:45:38 UTC闭合74秒，root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU16GiB15分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir及唯一日志相符，batch MaxRSS903300K。stderr末尾无输出。

JOBID13427583、STATUS COMPLETED、started.json1409 bytes、windows.jsonl9754 bytes及summary.json14850 bytes齐全。实际八行与summary records逐字段相同，身份／顺序与原八窗一致；fine／grouped-region位置数均与13425772相符。summary计时9.677805秒从TRAIN metadata加载结束后开始，含八窗缓存／标签／观察前向，不等于Slurm总74秒。观察前向各窗0.0351–0.5312秒；只有带原角色约束的零势观察事件，不是实际模型NLL、无条件配分或拟合耗时。

原八窗有四窗可行、四窗不可行：拟南芥1／4可行，水稻3／4可行。四个不可行窗均有单个边缘角色直接排除完整集合的充分矛盾，且原metadata有相应跨窗primary；这些transcript的CDS_complete_supported均为true，即染色体上的完整支持不等于本窗口中的完整起止。

| 物种及原证据类别 | 定向窗口身份 | 支持 | 确定边缘角色障碍或零势log集合数 |
| --- | --- | --- | --- |
| 拟南芥 新增已知最多 | NC_003071.7 12582912 + | 不可行 | 左端仅角色3 CDS2 |
| 拟南芥 新增未知最多 | NC_003076.8 10321920 − | 不可行 | 左端仅角色2 CDS1 |
| 拟南芥 正链exon最多 | NC_003070.9 25190400 − | 不可行 | 右端仅角色4 intron0 |
| 拟南芥 span冲突最多 | NC_003070.9 19611648 − | 可行 | log集合数72.96382371476446 |
| 水稻 新增已知最多 | NC_089043.1 17276928 − | 可行 | log集合数1.9459101490553132 |
| 水稻 新增未知最多 | NC_089046.1 3416064 − | 不可行 | 左端仅角色2 CDS1 |
| 水稻 正链exon最多 | NC_089046.1 565248 − | 可行 | log集合数0.6931471805599453 |
| 水稻 span冲突最多 | NC_089044.1 21467136 − | 可行 | log集合数0.6931471805599453 |

对应跨窗primary的完整定向span分别为[−390,1042)、[−1738,2473)、[24429,26414)和[−862,9771)，有效窗口均为[0,24576)。完整窗口gene集合首位置只能背景0或起始13、末位置只能背景0或终止14；因此这些原精细角色观察与当前边界语法不相容，直接令观察事件概率为0，不能作为有限训练NLL。此为支持域／边界矛盾，不是原R4标签文件错误、数据破坏或模型性能FAIL；放宽边缘条件后是否还有其他支持障碍尚未检验。

尤其水稻“新增未知最多”窗有原chain_unknown区[0,9771)，但fine 24576／coarse0／完全未知0，左端角色2仍为已解局部标签。链未知与局部角色未知不是同一对象，不能拿chain_unknown区默默删除这9771个已知角色来使配分可行。四个可行窗也只证明原观察具有非空完整集合支持；零势log集合数不是拟合后的概率、gene准确率或“所有参考均被模型预测”。

summary中原R4重新生成、原TRAIN catalog读取、八份原缓存读取与研究语法观察支持评估均为实际事实；cached_features_used、model_constructed、model_or_backbone_forward、research_weights_read、optimizer_constructed、checkpoints_written、original_R4_policy_changed、observation_rows_deleted_or_boundary_repaired、old_cached_labels_used、loss_normalizer_or_fit_coefficients_selected、DEV_used、independent_species_used、A33_outputs_opened及training_or_inference_caller_changed均为false，optimizer_steps_this_run为0。没有新训练、原A44恢复、额外样本选择或预测写出。

结论是保留原R4精细／粗粒度标签的完整语法监督必须先解决跨窗边界支持，不能把不可行窗移除后将剩余四窗冒充完整TRAIN新拟合。下一工作应准备保持原标签的边界条件及其数学证据，明确外部上下文／开放边界与完整raw输出的关系，再固定训练损失单位及首训合同；不自动修改原观察、cap、原R4文件或历史停止条件。八窗属于事前接线证据，4／8不是全部1491窗的不可行比例；完整TRAIN可行性、模型收益及独立泛化仍未成立，论文目标保持ACTIVE／NOT ACHIEVED。

## 原八窗完整primary集合与边界补齐证据准备

本轮保持13427583原八窗、顺序、全部原training_primary及所有原fine／coarse／unknown观察，问一个更具体的问题：这些固定primary的完整链能否组成非重叠集合，并在核心窗口逐位置满足原观察？扩展范围取核心与原primary完整CDS span的并，仅作为存在性证据坐标；该范围由TRAIN注释推导，不是raw-FASTA推理规则、选定训练padding或开放边界概率模型。

实际两物种genome.fa及其解析后目标均无现成.fai，data内也未发现TRAIN专用FASTA／索引。因此不扫描包含非TRAIN序列的未索引FASTA、不创建索引或放行新raw数据。仅复用原1491 TRAIN缓存的DNA；参考坐标从原窗口ID与原窗口metadata连接，方向沿原rc。原八份缓存全部保留，其他缓存只在其能补齐当前缺失位置时读取；重叠每个碱基必须一致，不是仅比较长度。metadata显示四个边界冲突窗中只有拟南芥NC_003071.7窗存在可补齐的相邻TRAIN缓存，其余三窗外部DNA覆盖不足；实际读取结果仍以作业内校验为准。

新 `probe_full_grammar_R4_completion_witness.py` 对所有原primary检查原cap128、首末exon至少3bp、intron至少5bp、完整跨度及相互span重叠。序列可得的链逐一用独立拼接CDS检查ATG／末STOP／无内部STOP／无歧义，再与原ORFIndex实际carry逐exon延伸核对；GT／GC–AG沿原DNA直接核实。缺外部DNA的边界链只报告历史CDS_complete_supported及几何，不冒充当前DNA验证，窗内完整链仍用原core DNA验证。没有删primary、改cap、改标签或搜索较易通过的参考子集。

沿原candidate_roles渲染完整链，只比较原core中所有允许角色；单链只检查自身span，非重叠完整集合另外检查gene之间的背景。primary跨度重叠时不按遍历顺序覆盖出伪集合，记录具体冲突且完整集合角色比较为未执行。只有完整DNA、所有链合法、无span重叠且core全部观察相容时才记录DNA_verified_completion_witness；其他结果区分原集合重叠、链语法不支持、观察不相容及DNA未取得。固定primary集合失败不证明所有可能gene集合都失败，历史metadata相容不等于当前DNA证据。

唯一入口 `M28-FULL-GRAMMAR-R4-COMPLETION-WITNESS.sbatch` 为generanno／private／0GPU2CPU16GiB5分钟／no-requeue，独立根 `outputs/M28-FULL-GRAMMAR-R4-COMPLETION-WITNESS` 与 `M28FULLR4COMPLETE_%j` 日志。原metadata加载前次约64秒，本轮不计算配分，只直接检查38条固定primary和最多必要邻窗缓存，5分钟为一次执行观察上限；产物started、逐窗windows.jsonl及八窗闭合后summary，小于1MiB。只执行一次，来源／坐标矛盾、实现异常、TIMEOUT／OOM均保留证据停止，不自动修复重提；预期的集合冲突或上下文缺失是诊断结果，不当作程序失败。

无模型、feature前向、配分／NLL、反向、optimizer、新拟合、padding决策、raw／DEV／独立物种／Setaria／A33访问或预测写出；原A44 TIMEOUT停止条件与固定监控PAUSED保持不变。接下来按实际结果确定边界数学所需的支持，而不是将注释推导的补齐范围带入最终caller或将本诊断包装为性能提升。

## 原八窗完整primary集合与补齐诊断实际完成

13:16:07 UTC五项活调度／预约通过；新三文件与canonical逐字相同，报告旧部分先与canonical逐字比较以保留并发改动，Python AST、bash语法、精确空输出根与唯一日志通过。13:18:28 UTC账户展开8项，加本次为9，小于80，唯一提交 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-R4-COMPLETION-WITNESS.sbatch`，返回13427817。实际13:18:28–13:19:54 UTC闭合86秒，root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU16GiB5分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir及唯一日志相符，batch MaxRSS895632K，stderr末尾无输出。

JOBID13427817、STATUS COMPLETED、started.json1198 bytes、windows.jsonl24331 bytes及summary.json35881 bytes齐全。实际八行与summary records逐字段相同，身份／顺序与started及13427583原八窗一致，原fine／coarse／完全未知数均不变。原1536draw／1491唯一TRAIN窗、每窗全部primary身份及数量均保留。summary计时14.737442秒从原TRAIN metadata加载结束后开始，不是Slurm总86秒。

全部八窗的原完整primary集合均无span重叠；按完整链几何渲染后，所有原core允许角色均相容，逐窗集合及单链观察不相容位置数都为0。原38条primary中35条获得当前缓存DNA的完整语法验证，直接拼接CDS与原ORFIndex实际carry verdict一致，首末／内部STOP、无歧义、原GT／GC–AG、intron≥5及cap128均通过；其余3条跨窗链缺外部DNA，未给当前完整语法判定。38条的历史CDS_complete_supported均为true，但不将其替代这3条缺失的实际序列验证。

五窗有当前DNA支持的完整集合存在性证据：原已可行的四窗加上拟南芥NC_003071.7 12582912 +窗。后者原闭窗事件不可行，实际用原相邻负链缓存 `window_00789.pt` 补入左侧390bp，扩展基因组范围[12582522,12607488)，定向长度24966bp，core offset390。与原正链core共有的12288bp逐碱基一致，且全部原core DNA保持不变。原10条primary的完整链均合法、非重叠、全部原21768 fine／173 coarse／2635 unknown行相容；跨窗rna-NM_128482.4完整span为扩展坐标[0,1432)，5exon。没有删标签、primary或改阈值，证明此窗确有保留原角色观察的真实DNA完整链补齐，而不是通过静默mask制造可行。

另外三个原闭窗不可行窗仍未取得完整补齐DNA，不记作FAIL或已解决：

| 原窗口 | 缺失基因组DNA半开区间 | bp | 当前未验证完整链 |
| --- | --- | ---: | --- |
| 拟南芥 NC_003076.8 10321920 − | [10346496,10348234) | 1738 | rna-NM_122723.4 |
| 拟南芥 NC_003070.9 25190400 − | [25188562,25190400) | 1838 | rna-NM_105396.4 |
| 水稻 NC_089046.1 3416064 − | [3440640,3441502) | 862 | rna-XM_015762904.3 |

这三窗的完整链几何与原core观察相容且无span重叠，窗内完整链已由core DNA验证；只有边界链的外部序列未取得。原八缓存加唯一必要邻窗共9份artifact，features未用于前向；未扫描FASTA、创建索引、填造外部DNA或访问其他数据划分。完整集合证据5／8不是全部1491窗的比例，也不是原四个冲突全部解决；新直接证据只覆盖其中一个冲突窗。

summary实际确认optimizer_steps_this_run为0；model_constructed、optimizer_constructed、model_or_backbone_forward、research_weights_read、checkpoints_written、original_R4_policy_changed、observation_rows_deleted、primary_transcripts_dropped、external_DNA_fabricated、unindexed_FASTA_scanned、raw_FASTA_opened、partition_or_NLL_computed、training_or_inference_caller_changed、DEV_used、independent_species_used及A33_outputs_opened均为false。没有新训练、原A44恢复、预测／checkpoint写出或监控变更。

本轮取得一个真实边界补齐存在性证据，并排除这八个固定primary集合的额外角色／span矛盾；未定义边界先验、未计算扩展观察logZ／梯度、未选损失单位或新首训合同，也没有模型性能／泛化／论文就绪证据。注释派生范围只用于诊断，不能作为无参考推理padding规则。下一步准备保留跨窗角色的边界条件及其配分数学，区分原core完整语法与新的边界语义，再考虑实际TRAIN可行性；原A44停止条件不变，论文目标保持ACTIVE／NOT ACHIEVED。

## 完整语法边界消息数学准备

`full_grammar_boundary_messages.py` 在同一给定全局DNA及15角色势上，将完整非重叠gene集合配分改写为逐碱基有限状态消息。内部切点不是新的完整gene边缘：只在给定全局序列起点初始化背景状态、终点取背景状态质量；消息穿过任意内部分块，跨块gene仍须在全局DNA上完整起止。没有均匀开放边界先验、Gold初始化、参考候选注入或截断gene输出；给定全局DNA的范围也没有由该模块选择。

消息保存gene内exon计数、实际未完成密码子字符串、ATG／STOP已读部分、donor GT／GC、intron最短长度成熟阶段和acceptor AG阶段。先消费明确起止密码子的三个角色；普通exon末碱基发donor角色，middle exon首碱基发acceptor角色。长度3的first exon仍保留start，长度3的terminal exon全部为stop，长度1的middle exon由donor覆盖acceptor，沿原candidate_roles优先级。intron需donor2＋至少1内部碱基＋acceptor2，即原5bp下限；actual carry在每次完成密码子时拒绝内部STOP，终止STOP必须同一exon内连续。cap附在BoundaryMessage中随分块保留，不在每块重置。

每个状态累加所有到达历史的log质量，不保留单个赢家、不枚举候选池；评分仍为每个最终角色的logit减同位置intergenic logit，FP32／FP64输入、FP64累积、unit role／null0。原区间配分、观察算子、compiled reverse和所有研究caller保持不变；新的native Torch表示目前仅用于人工数学准备，没有完整实际窗口或染色体资源证据。

若全局DNA分块转移记为K，各块得分消息连续相乘，只有全局两端为背景边界。核心观察只限制其原允许角色，外部DNA与势仍是实际输入，外部角色保持未知；相应未归一化NLL为同一全局空间的logZ_all减logZ_core_observed。梯度应穿过核心以及外部上下文势，不能detach外部消息后冒充同一目标。该接口不获取外部DNA、不决定有限halo／padding、不处理缺失序列，不选损失分母、importance规则或fit系数；部分序列可组成消息，但未完成gene的状态不接受为全局完整输出。

新十项检查针对新状态表示可能改变科学支持或梯度的失败：完整gene集合logZ及FP32／FP64全部角色梯度对独立拼接DNA穷举及原观察算子；所有单切点及逐碱基切块的值／梯度完全相同；literal motif／实际carry／intron成熟状态；核心首尾均在gene内时的外部前后缀支持、NLL与完整上下文梯度；六gene的未知空间64集合及观察后16集合；跨块cap和短exon覆盖优先级；跨块／跨exon内部STOP拒绝；全部未知时零loss／梯度及common-row-shift不变；GC donor／最短intron／背景约束；全局不完整gene拒绝和空块恒等。保留原值小数点后10位、FP64 rtol1e−9／atol1e−10、FP32 rtol／atol1e−6、上游−0.7／2.3；不改原fixture或放宽容差求通过。

唯一入口 `M28-FULL-GRAMMAR-BOUNDARY-CHECKS.sbatch` 使用generanno／private／0GPU2CPU8GiB2分钟／no-requeue，独立根 `outputs/M28-FULL-GRAMMAR-BOUNDARY-CHECKS` 与 `M28FULLBOUNDCHK_%j` 日志。读13427476／13427817已闭合summary仅用于来源接线，不重跑其检查或读研究DNA、features、标签、权重。产物started与仅十项实际PASS后写出的summary，小于1MiB；参考此前人工检查14秒总时长，2分钟为本次固定观察上限，不是完整TRAIN投影。明确实现或新人工fixture错误最多允许同输入／资源／容差的一次最小R1；数学矛盾、TIMEOUT／OOM或第二次失败停止，保留失败证据。

此准备不修改R4观察、真实TRAIN窗口、历史结果、A44原TIMEOUT停止条件、DEV／独立数据权限、已有监控或训练／推理caller，也不提交新首训。若通过，下一步才评估该边界表示的实际窗口资源及既有TRAIN上下文可用性；缺外部DNA的三个证据窗仍不能由人工PASS补证，论文完整gene性能与适用范围目标保持原范围。

## 完整语法边界消息人工验证实际完成

13:44:31 UTC五项活调度／预约通过；新模块、十项检查、SBATCH及报告四文件与canonical逐字相同，报告旧部分先比较以保留并发改动，Python及内嵌runner AST、bash语法、实际test方法数10、空输出根与唯一日志通过。13:47:57 UTC账户展开8项，加本次为9，小于80，唯一提交 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-BOUNDARY-CHECKS.sbatch` 返回13428297。实际13:47:58–13:48:30 UTC闭合32秒，root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir及唯一日志相符，batch MaxRSS643956K。

JOBID13428297、STATUS COMPLETED、started.json174 bytes及summary.json1501 bytes齐全。summary身份、两个已闭合source job13427476／13427817、tests10／failures0／errors0与范围字段相符；实际stderr有十个test名称逐项ok、Ran 10 tests in 5.180s及OK，没有跳过或R1。summary计时14.266265秒从suite发现前开始，包含检查模块加载，不等于unittest的5.180秒或Slurm总32秒。

人工完整集合logZ及全部FP32／FP64输入梯度匹配独立拼接DNA完整chain／非重叠子集穷举及原观察配分。三个短DNA的全部内部单切点及逐碱基切块，logZ与整段值逐值相同、全部梯度rtol0／atol0相同；ATG／STOP切中、实际未完成密码子C／CT、donor／intron成熟／acceptor状态、跨块exon计数及原短exon角色优先级都获得检查支持。跨exon内部STOP、少于5bp的intron、禁止背景观察和全局未完成gene仍被拒绝；没有通过开放全局末端改变原完整gene语法。

在人工 `ATGCGTTAGTTTAA` 的核心区[3,12)上，核心首角色为donor7，末端在终止密码子内部；单独按完整闭窗语法时原观察不可行。沿真实给定的人工前后缀状态消息后，完整全局观察只留一个合法集合，core NLL及全部输入梯度均匹配独立穷举与原全局观察NLL，核心外前缀和后缀梯度实际非零。这里的前后缀是该人工DNA中的实际字符，不是缺失研究DNA的代用品；角色未删除，外部消息未detach，也没有选择一个均匀边界先验。

六个邻接人工gene的无条件集合数64、核心部分观察后16，NLL为log4，未知及核心外gene仍边缘化而非变成负例。全部未知观察时loss与全输入梯度严格为0；原FP32共同row shift下值及梯度不变。完整值／oracle比较继续采用原小数点后10位、FP64 rtol1e−9／atol1e−10、FP32 rtol／atol1e−6及上游−0.7／2.3，没有放宽容差或重提。

summary实际确认synthetic_only、CPU_first_order_verified、完整全局语法logZ／全梯度对照、所有内部切点一致、实际carry／motif／cap／角色优先级、core NLL／完整上下文梯度、全局不完整拒绝及未解gene边缘化为true，source_metadata_summaries_read为true。optimizer_steps_this_run为0；boundary_prior_chosen、context_extent_or_padding_rule_selected、loss_normalizer_or_fit_coefficients_selected、research_DNA_features_labels_or_weights_opened、model_or_backbone_forward、original_R4_policy_changed、observation_rows_deleted、DEV_used、independent_species_used、A33_outputs_opened及training_or_inference_caller_changed均为false。

结论是逐碱基边界消息及核心观察似然在这些人工完整集合上获得数学支持，可以保留跨内部切点的状态而不强制每块完整起止。native Torch在完整24,576bp／cap128或实际染色体上的资源尚未验证，真实TRAIN前后缀及势的来源也未接入；三个缺外部DNA窗口状态不变。下一步应准备可执行的完整长度边界消息表示及资源证据，再确定实际TRAIN上下文合同，不能将人工短窗PASS当作新首训、性能改善、独立适用范围或创新性证据。原A44停止条件与所有历史失败保持不变，论文目标仍为ACTIVE／NOT ACHIEVED。

## 编译边界消息及检查点反向准备

本轮新增 `src/m28/full_grammar_boundary_compiled.py`，原native模块和全部研究caller不变。状态表从原 `_transitions` 对ACGT及非ACGT等价类作结构闭包，固定cap、实际密码子carry、motif阶段及角色优先级；不按分数、winner、参考或候选池裁剪。每块输入、输出为同一固定状态表上的FP64 log mass，缺失状态为负无穷；自定义一阶CPU autograd对势及传入消息同时求导，不在内部切点关闭gene或detach消息。全局仍只从BG0开始、在真实给定DNA末端读BG。

Numba前向使用原状态排序、原edge顺序的稳定logaddexp，fastmath及parallel关闭。仅保存固定间距前向检查点；反向逐块按同一运算重建alpha并传播完整状态／角色质量。checkpoint_bp默认256只决定存储及重算，不能作为context范围、边界先验、padding或损失规则。角色势仍为原logit减背景logit，背景梯度相应减全部角色占用；不选择loss normalizer或fit系数，明确仅一阶CPU。

新十项人工检查针对新表示的实际风险：48个完整集合值／全FP32及FP64梯度对独立拼接DNA集合穷举、native消息和原观察配分三方；每个有质量的输出状态及带正负权重的向量readout梯度；在ATG、donor、intron成熟、acceptor、实际carry及STOP切点传入消息的梯度；全部内部单切点及逐碱基分块；原core[3,12)观察及完整上下文NLL／外部梯度；随机部分角色观察；cap、短exon覆盖、内部STOP、GC donor及最短intron；64／16集合边缘化、未知零loss及row shift；检查点间距1／2／3／5／256的值及全梯度不变并验证实际nopython；全局不完整拒绝、空块恒等、非法输入和零质量状态的未定义readout拒绝。原oracle值小数点后10位、FP64 rtol1e−9／atol1e−10、FP32 rtol／atol1e−6、上游−0.7／2.3不变。

唯一新入口 `M28-FULL-GRAMMAR-BOUNDARY-COMPILED-CHECKS.sbatch`，独立根 `outputs/M28-FULL-GRAMMAR-BOUNDARY-COMPILED-CHECKS`，日志 `M28FULLBCOMPCHK_%j`，generanno／private／0GPU2CPU8GiB2分钟／no-requeue。仅从实际完成13428297的summary核对来源；不读研究DNA、features、标签、权重，不做model forward、优化或新首训。十项实际PASS后才写summary；源文件、AST、bash及唯一根／日志核对后单次提交。明确实现或新人工fixture错误最多允许同输入／资源／容差的一次最小R1；数学矛盾、TIMEOUT／OOM或第二次失败停止，保留原证据。

数学实际通过后，才对同一原完整合成输入（DNA seed0先生成49,152bp后取前24,576bp，势seed41整数−3至3、15角色FP32、cap128、单位势、null0及FP64累积）准备新表示的单次完整资源探针。拟用同0GPU2CPU16GiB2分钟、checkpoint256，新根／日志，首次前向及反向JIT不预热。全长logZ对已实际完成13427267的10992.597775197577采用事前固定rtol0／atol1e−8（完整FP64不同归约），短数学容差不变；全24576×15梯度必须finite、row-sum≤1e−5、nonBG占用≥−1e−6且逐行≤1+1e−6。新FSA表示不虚报为旧八种递归／tape计数相同；单独记录实际状态、检查点／重建块bytes和进程／Slurm RSS。失败不缩短窗口、降cap、延时或重提。此准备不是实际全长资源PASS、真实TRAIN上下文合同、性能或适用范围证据。

## 编译边界消息数学检查实际完成

14:17:57 UTC五项活调度／预约核对，展开账户7项、加本次为8小于80；旧报告先与canonical逐字比较，新增模块／tests／SBATCH／报告逐文件比较均最终exit0，保留所有既有内容。一次模块比较曾因SSH banner exchange超时失败，只重做读文件比较后成功，没有重复提交作业。AST、内嵌runner AST、bash、实际方法数10及唯一输出根／日志通过。14:20:44 UTC单次 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-BOUNDARY-COMPILED-CHECKS.sbatch` 返回13428495，即时queue为PD(None)。Slurm实际14:20:45–14:21:15 UTC共30秒（该次未设置TZ的sacct／scontrol显示CEST16:20:45–16:21:15），root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU8GiB2分钟／no-requeue／Restarts0，实际Command、SubmitLine、WorkDir及唯一日志吻合，batch MaxRSS695804K。

JOBID9 bytes、STATUS COMPLETED10 bytes、started.json192 bytes、summary.json1609 bytes齐全且身份13428495／source13428297一致。stderr十个实际test逐项ok、Ran 10 tests in 9.037s及OK，summary tests10／failures0／errors0、计时9.322670秒；后者含suite加载，不是Slurm30秒。新进程forward／reverse nopython signatures均由0变1，首次编译确实执行；没有R1、改fixture或放宽容差。

全部48个完整集合值／FP32及FP64全势梯度对三方oracle通过；每个有限输出状态、正负权重readout、所有所测motif／carry／splice切点的传入消息梯度对native通过。全部内部单切点及逐碱基分块保留完整历史；检查点间距1／2／3／5／256的值及全梯度逐值相同。实际检查点bytes与逐块重建alpha bytes和分块geometry相符，并不据此声称它们等于进程RSS。core[3,12)单独闭窗仍不可行，完整给定上下文的NLL／全梯度与独立集合及native相符，核心外前后缀梯度非零；cap、短exon角色、GC donor、最短intron、跨exon STOP及全局完整约束保持。64／16集合边缘化、未知零loss／梯度、row shift不变和全部人工部分角色观察通过。

summary确认上述数学／一阶CPU／synthetic范围及来源summary读取；optimizer_steps_this_run为0。boundary_prior_chosen、context_extent_or_padding_rule_selected、loss_normalizer_or_fit_coefficients_selected、research_DNA_features_labels_or_weights_opened、model_or_backbone_forward、original_R4_policy_changed、observation_rows_deleted、DEV_used、independent_species_used、A33_outputs_opened及training_or_inference_caller_changed均为false。结论仅为编译、检查点反向和消息组合的人工数学支持，全长资源、真实context及研究训练效果仍未由此验证。

## 编译边界消息完整合成窗口资源探针准备

新 `probe_compiled_boundary_window_synthetic.py` 和 `M28-FULL-GRAMMAR-BOUNDARY-COMPILED-WINDOW-PROBE.sbatch` 明确绑定实际数学13428495及实际原完整基准13427267，来源身份／固定输入／完成状态不符就停止；全长DNA motif计数须为原417 start／3098 donor／1494 acceptor／1121 STOP。首次前向计时包含结构状态闭包及Numba JIT，首次反向包含其JIT，不做预热；保存started、forward（包括不符时实际值）及仅通过后写出的summary。不保存24,576×15梯度文件、不读研究样本，产物小于1MiB。

唯一新根 `outputs/M28-FULL-GRAMMAR-BOUNDARY-COMPILED-WINDOW-PROBE`、日志 `M28FULLBCOMPW_%j`，generanno／private／0GPU2CPU16GiB2分钟／no-requeue，checkpoint256；其余完整输入、势和事前容差按上一准备段不变。14:22:32 UTC再次五项活调度／预约通过，账户仍7项，新根与日志为空。只在源文件核对后单次提交，TIMEOUT／OOM、数学不符或运行失败不重提、不减长度／cap、不延长。短数学证明及完整资源门槛均不是新首训；完整基准未保存全梯度，因此本探针不宣称完成了两种实现全长梯度逐值对照，也不比较不同表示的旧八种递归／tape计数。

## 编译边界消息完整合成窗口实际完成

新probe／SBATCH／报告的AST、bash及canonical逐字核对通过，14:24:32 UTC展开账户7项、加本次8小于80，输出根与唯一日志为空；单次 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-BOUNDARY-COMPILED-WINDOW-PROBE.sbatch` 返回13428516，即时queue已R／gpu034。显式TZ=UTC查询确认实际14:24:33–14:24:55共22秒，root／batch／extern均COMPLETED0:0；private／generanno／gpu034／0GPU2CPU16GiB2分钟／no-requeue／Restarts0，实际Command、SubmitLine、WorkDir及唯一日志一致，batch MaxRSS696908K。

JOBID9 bytes、STATUS COMPLETED10 bytes、started.json864 bytes、forward.json419 bytes、summary.json3104 bytes全部非空。started字段逐项保留在summary，forward与summary.forward逐字结构相等，stdout三份JSON分别与started／forward／summary完全相等，stderr为空；数学来源13428495、完整原基准13427267、输入geometry及所有范围字段核对通过。前向／反向nopython signatures各从0变1，没有预热、R1、输入缩短、cap下降或时长扩展。

同一完整24,576bp／cap128输入，实际前向4.657314779秒（含状态闭包及首次JIT），反向3.419833745秒（含首次JIT）。新logZ10992.597775197613对原10992.597775197577差3.637978807091713e−11，低于事前rtol0／atol1e−8。全24576×15 FP32梯度finite，共同row shift梯度残差1.1920928955078125e−7低于原1e−5，nonBG角色梯度最小0，逐行nonBG占用最大1.0000001192092896满足原1e−6容差。没有通过改变语法或观察来取得完成。

实际固定表19,078状态、97行检查点、检查点14,804,528 bytes；反向最大重建alpha块39,224,368 bytes。前向进程峰值594,184 KiB，完整进程峰值646,060 KiB；Slurm batch MaxRSS696,908K为另一实际调度指标，不能将前两项数组bytes或进程RSS与之混用。相同固定合成输入的单次冷启动记录为：

| 指标 | 已完成原表示13427267 | 新消息表示13428516 |
| --- | ---: | ---: |
| 前向秒 | 81.784235207 | 4.657314779 |
| 反向秒，含首次Numba编译 | 7.302722387 | 3.419833745 |
| 完整进程峰值RSS，KiB | 2,812,872 | 646,060 |

这三项是同一实际合成input／角色势／cap及0GPU2CPU的资源记录，不是重复benchmark、真实基因组训练吞吐、GPU可行性或性能比较；原完整梯度没有保存，因此不宣称新旧全长梯度逐值一致，数学一致性来自已完成的人工三方／消息／检查点检查，本全长梯度只取得finite及占用不变量证据。新FSA状态数也不等于原递归／tape的八种计数。

summary实际full_window_forward_and_backward_completed、all_input_gradients_finite、complete_fixed_input_logZ_matches_original_with_declared_tolerance、all_input_role_posterior_bounds_verified及冷JIT字段为true；optimizer_steps_this_run为0。boundary_prior_chosen、context_extent_or_padding_rule_selected、loss_normalizer_or_fit_coefficients_selected、model_or_backbone_forward、research_inputs_opened、labels_or_catalog_opened、DEV_used、independent_species_used、A33_outputs_opened、observed_full_window_NLL_executed、original_R4_policy_changed、observation_rows_deleted、training_or_inference_caller_changed、baseline_full_gradient_equality_checked及old_recurrence_or_tape_counts_compared均为false。

本轮已取得可执行的完整长度跨块消息表示及一阶CPU资源支持，不再仅停留于短人工数学。真实TRAIN上下文及势尚未接入，完整长度core-observed NLL也尚未运行；下一步应从已获证、实际可得的固定TRAIN上下文做零更新观察支持／NLL诊断，再按原授权确定研究训练合同。注释派生补齐范围只能作为存在性诊断，三个缺外部DNA的窗口仍未补证，不能把它们默认为可行或用虚构字符补齐。所有原失败、A44原TIMEOUT停止条件、DEV／独立数据边界及冻结科研规则保持不变；完整gene性能、广泛适用性、原创性及论文目标仍未达成，本轮是PROGRESS而不是科研PASS或goal完成。

## 原八窗真实上下文零势似然诊断准备

新 `probe_full_grammar_R4_context_nll.py` 保留13427817的全部八个固定TRAIN见证及顺序，绑定原观察13427583、编译数学13428495和完整资源13428516的实际完成记录。沿原Q.load_train_inputs的1536draw／1491唯一窗、两物种六条TRAIN染色体allowlist与混合catalog解析前split筛选加载；仅允许13427817已打开的九份原缓存，GEN feature只由原validate_cached检查geometry、不用于前向。使用原build_targets重新生成全部core fine／coarse／unknown观察并核对原位置数和primary IDs，不读旧缓存标签、不改窗口／draw／importance或目标生成器。

上下文只重建13427817已经记录的实际缓存片段、诊断坐标、core定向offset、重叠DNA及缺失区间；原注释推导的CDS补齐坐标仅用于本次存在性／似然诊断，不成为训练或无参考推理padding规则。五个实际可得上下文计算同一给定全局DNA上的logZ_all−logZ_core_observed，其余三窗完整保留并记录NOT_EXECUTED_EXTERNAL_DNA_MISSING，不能据此把不可获取的支持判为数学FAIL或从训练集删除。核心外实际DNA／势保持在可微消息中，外部角色全部未知；不选择开放边界先验、不detach、不填字符，也不扫描无索引原始FASTA。

沿用原诊断的全零FP32十五角色势、单位角色、null0、cap128、FP64累积和checkpoint256，不做模型／optimizer／权重读取。结果是未加权、未归一化的合法集合计数似然及在零势处的完整输入NLL梯度，不是已选定的训练损失单位或学习后的模型损失。五窗logZ及NLL必须finite且≥−1e−8；所有context×15梯度必须finite，common-row-shift残差≤1e−5。NLL梯度是两种角色后验占用之差，允许负值；每个nonBG角色差及逐行covered概率差绝对值应≤1+1e−6，不能误用无条件配分的非负概率要求。记录core及外部前后缀的实际梯度L1，不把非零值预先包装为结果。

四个未扩展上下文的观察logZ，直接与13427583已完成的72.96382371476446、1.9459101490553132、0.6931471805599453、0.6931471805599453按事前rtol0／atol1e−8核对，不重跑旧算子。原Arabidopsis首窗保留完整24,576bp core观察，在已有实际DNA的390bp前缀上检查新的全局支持、NLL及外部梯度；其原闭窗不可行记录仍是原空间的事实，新的诊断空间不能改写它。

唯一入口 `M28-FULL-GRAMMAR-R4-CONTEXT-NLL.sbatch`、根 `outputs/M28-FULL-GRAMMAR-R4-CONTEXT-NLL`、日志 `M28FULLR4CTXNLL_%j`，generanno／private／0GPU2CPU16GiB5分钟／no-requeue，沿已实际运行的原completion环境及资源路由。原completion总86秒、单次完整新资源22秒；按五对配分／反向及TRAIN元数据留本次固定5分钟上限，不是原A44超时的延长或新首训。14:39:02 UTC五项活调度／预约核对通过，账户展开7项小于80，根及唯一日志为空；源文件AST、bash及canonical核对后单次提交。明确新runner实现错误最多允许同输入／资源／容差的一次最小R1；实际数据漂移、数学矛盾、TIMEOUT／OOM或第二次失败停止。产物started、八行windows及仅本次诊断完整执行后写出的summary，小于1MiB，不保存完整梯度、checkpoint或预测。工程完成时仍明确八窗中仅五窗NLL可执行；本轮不放行DEV／独立物种／A33、改变R4或caller、取消原作业或重启A44。

## 原八窗真实上下文零势似然诊断提交

新runner AST、SBATCH bash及runner／SBATCH／报告三个文件的canonical逐字比较通过，报告原内容先核对并完整保留。14:43:10 UTC展开账户7项、加本次8小于80，唯一根及日志为空，单次 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-R4-CONTEXT-NLL.sbatch` 返回13428648，即时queue为PD(None)。后续读取确认该同一job在gpu034实际RUNNING、elapsed1:41，started.json2328 bytes已写出，错误日志当时为空。这是实际启动记录，不是五窗NLL／梯度或八窗诊断已经完成的证明；须待终态、逐行结果及summary核对后判断。

## 原八窗真实上下文零势似然诊断实际完成

13428648实际14:43:12–14:45:02 UTC运行110秒，root／batch／extern均COMPLETED0:0；scontrol提交时间为14:43:11 UTC，与本机提交调用返回记录14:43:10 UTC有一秒差。private／generanno／gpu034／0GPU2CPU16GiB5分钟／no-requeue／Restarts0，Command、SubmitLine、WorkDir及唯一日志吻合，batch MaxRSS1089656K。JOBID9 bytes、STATUS COMPLETED10 bytes、started.json2328 bytes、windows.jsonl16204 bytes、summary.json24165 bytes全部非空。

逐行八条记录与summary.records完全相等，保留来源八个window ID及顺序、全部fine／coarse／unknown位置数、实际片段及缺DNA区间；实际九个打开的缓存与原13427817及started允许列表相等。stdout十四份JSON为五个计算开始、八个逐窗记录及一份终summary，逐项与产物相符；stderr为空。原TRAIN元数据加载8.256855秒，模块导入后总51.698077秒，进程峰值992784 KiB；前两项不含模块导入，不能等同Slurm110秒或batch RSS。前向及反向nopython signatures各0→1，没有预热、R1或容差／数据改动。

五个实际可得上下文均完成观察支持、两类配分、未加权未归一化NLL及全部context×15 FP32输入梯度。零势结果如下，NLL数值不能跨观察量不同的窗口排名为模型性能：

| 原序号与core窗口 | context bp | logZ all | logZ core observed | NLL |
| --- | ---: | ---: | ---: | ---: |
| 0 NC_003071.7 + 12582912 | 24966 | 1153.396869138 | 123.598396688 | 1029.798472450 |
| 3 NC_003070.9 − 19611648 | 24576 | 946.989225616 | 72.963823715 | 874.025401901 |
| 4 NC_089043.1 − 17276928 | 24576 | 1198.610704548 | 1.945910149 | 1196.664794399 |
| 6 NC_089046.1 − 565248 | 24576 | 1266.810719671 | 0.693147181 | 1266.117572490 |
| 7 NC_089044.1 − 21467136 | 24576 | 1274.238462772 | 0.693147181 | 1273.545315591 |

原序号0保持core原21768 fine／173 coarse／2635 unknown行，在先前实际证实的[12582522,12607488)上下文、core定向[390,24966)、内部切点390计算。观察logZ finite、NLL1029.7984724501464；全部梯度finite，core梯度L1为38248.75821282812，实际390bp前缀梯度L1为81.71318615517703，后缀为空所以梯度为0。它取得了原core观察在实际DNA前缀上的配分支持及非detach外部梯度，而非仅人工数学或primary存在性。原闭窗不可行仍是原范围的事实；这里的上下文范围仍是注释推导的诊断坐标，不是已经合法化的训练／raw推理padding或边界先验。

四个原未扩展窗口的观察logZ对13427583实际结果均满足事前rtol0／atol1e−8；最大差仅约4.26e−14，未重跑旧实现。五窗NLL梯度common-row-shift最大残差1.7881393432617188e−7≤1e−5，最大nonBG单角色后验差0.9999710321426392，最大逐行covered后验差1.0000001192092896，均满足原1e−6界；NLL梯度是后验差而不是非负概率。首窗NLL前向5.076094秒、反向4.653044秒含首次编译；其余四窗各前向0.918764–0.957220秒、两类配分合成的反向3.592597–3.726098秒。这些是零势、无模型的单次CPU运行记录，不能推出GEN模型训练吞吐或GPU可行性。

原序号1／2／5仍分别缺原记录1738／1838／862bp外部DNA，状态均为NOT_EXECUTED_EXTERNAL_DNA_MISSING，NLL、梯度及新范围观察支持未计算，context_observation_feasible为null。summary实际available_context_NLL_and_gradient_windows为5、missing_context_NLL_NOT_EXECUTED_windows为3、all_original_witnesses_accounted为true、all_original_context_NLLs_executed为false；不能把五窗可计算写成八窗或1491窗训练支持成立，也不能把三窗未执行写成数学失败。

summary确认本次synthetic_only为false、原TRAIN DNA及catalog读取和R4 targets重建为true，optimizer_steps_this_run为0。cached_features_used、old_cached_labels_used、research_potentials_or_weights_read、boundary_prior_chosen、training_context_extent_or_padding_rule_selected、loss_normalizer_importance_rule_or_fit_coefficients_selected、model_constructed、optimizer_constructed、model_or_backbone_forward、checkpoints_written、original_R4_policy_changed、observation_rows_deleted、primary_transcripts_dropped、external_DNA_fabricated、unindexed_FASTA_scanned、raw_FASTA_opened、DEV_used、independent_species_used、A33_outputs_opened及training_or_inference_caller_changed均为false。

本轮取得真实固定TRAIN上下文上的观察似然／全势梯度证据，并在一个原边界冲突实例中验证真实前缀参与反向传播；这推进了完整gene集合目标，不是新模型训练或性能PASS。下一步应核实剩余原TRAIN序列的合规可得路径及同方向势的来源，补齐三个未执行实例，再确定不依赖参考边界的实际训练／推理上下文合同。缺序列不能变成删除窗口的理由，也不能用反向邻窗feature直接翻转替代未经核实的同向GEN表征。原A44停止条件、所有历史结果与DEV／独立数据边界不变，论文目标仍ACTIVE／NOT ACHIEVED。

## 原八窗 TRAIN FASTA 上下文补齐准备

本轮从原始两物种数据目录的 `genome.fa` 重建同一八个固定TRAIN诊断上下文，补齐原缓存缺失的1738、1838及862bp；13428648的五窗已执行／三窗未执行记录保持不变。新runner绑定同一四份实际来源及13428648，原1536draw／1491唯一窗、六条TRAIN染色体、八个ID及顺序、九份缓存、全部core观察和primary IDs均不变。上下文继续使用原注释推导的core与primary CDS并集，只用于数学／数据诊断，不将它选作训练或无参考推理padding规则。

原FASTA沿用未修改的 `manifest.read_allowed`。原配置中的物种路径先按Q.ALLOW筛成Arabidopsis及rice，长度只从这六条TRAIN记录的sequence-region头元数据取得；不调用会读val split的旧 `census.metadata`，不打开split_val／test文件，不生成全基因组FASTA索引。读取器会流式经过同一原文件中的其他记录行，但在记录门处丢弃，不解析、拼接、保留或使用非TRAIN序列。原配置中Setaria等字符串不成为文件访问。本轮为同一已授权TRAIN数据的读取实现准备，不是新数据获取或封存数据放行。

每窗将定向FASTA core的全部24,576bp与原缓存逐碱基比较；五个已有完整缓存上下文还比较其全部字符，包括已有390bp前缀。原缓存片段、缺失区间及DNA证据仍单独保留。只对原先DNA不可获取的三个primary重新运行几何、拼接CDS与实际carry语法检查；其余35条已获证DNA通过上述字符一致性转移证据，不宣称重新逐条审计。实际DNA不符、语法拒绝或数据漂移即停止，不填字符、不改label、不删行。

八窗全部运行相同全零FP32十五角色势、单位角色／null0、cap128、FP64累积、checkpoint256的未加权未归一化core-observed NLL和完整输入梯度。沿用finite、非负计数／NLL及原梯度界；五个已有结果作为新FASTA接入的正控，三个logZ／NLL量分别与13428648按事前rtol0／atol1e−8比较。记录三个新外部片段实际梯度，不要求非零以改写结果。只保存八条小型定向context、逐窗统计和summary，不复制整条染色体、不保存梯度或模型；本轮无权重、模型、optimizer、feature前向或拟合。

新入口 `M28-FULL-GRAMMAR-R4-TRAIN-FASTA-CONTEXT.sbatch`、独立根 `outputs/M28-FULL-GRAMMAR-R4-TRAIN-FASTA-CONTEXT`、唯一日志 `M28FULLR4FACTX_%j`，沿原generanno／private／0GPU2CPU16GiB／5分钟／no-requeue。原五窗作业110秒，固定5分钟覆盖新增原FASTA流及三个配分／反向；提交前核对五项活调度、账户80上限和根／日志碰撞。本runner明确实现错误最多允许同输入／资源／容差的一次最小R1；实际数据或数学矛盾、TIMEOUT／OOM、第二次失败停止。原A44超时及依赖停止保持关闭，不取消、重提或延长旧作业。即使八窗通过，也只是固定TRAIN上下文上的零势数学支持，不能推出1491窗支持、新模型性能或论文PASS。

15:17:19 UTC五项活调度及预约核对完成，无本次private短CPU任务的预约冲突；账户当时七项。runner AST及SBATCH bash通过，新入口及输出根不存在，报告原82,804 bytes与canonical逐字相同后完整保留；仅同步本runner、SBATCH和追加的现有报告，提交时再核对展开账户及唯一日志。

## 原八窗 TRAIN FASTA 上下文实际完成

本次八窗全部取得真实原TRAIN DNA上的有限观察似然和完整角色势梯度；三个先前缺缓存片段的案例不再处于本次未执行状态。单次提交 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-R4-TRAIN-FASTA-CONTEXT.sbatch` 于15:18:58 UTC返回13428966；即刻核对同一job为gpu034上的RUNNING。实际运行15:18:59–15:21:01 UTC共122秒，root／batch／extern均COMPLETED0:0，private／generanno／0GPU2CPU16GiB／5分钟／no-requeue、Restarts0，Command、WorkDir及唯一日志吻合。没有R1、重提、容差或数据改动，batch MaxRSS1,270,380K。

JOBID9 bytes、STATUS COMPLETED10 bytes、started.json2,683 bytes、windows.jsonl24,976 bytes、summary.json42,074 bytes、contexts.jsonl203,706 bytes，总273,458 bytes。八条逐窗结果与summary.records完全相等，八条context的ID、顺序、定向core区间及字符长度与相应结果／原诊断坐标一致；stdout二十份JSON为八个计算开始、三个新增primary DNA检查、八个逐窗结果和一份终summary，逐项吻合，stderr为空。原1536draw／1491唯一窗、六条TRAIN染色体、八窗全部位置数／primary IDs、原九份缓存及原缺DNA区间都保留。

原FASTA实际尺寸Arabidopsis121,164,952 bytes、rice140,582,944 bytes；读取器只解析、拼接和保留原六条TRAIN记录，未调用census.metadata或打开val／test split文件。每个24,576bp core均与原缓存逐碱基一致，五个先前完整上下文的所有字符也一致。五个已有案例的logZ_all、logZ_core_observed及NLL与13428648的实际值完全相等，最大绝对差为0，满足事前1e−8正控界；这是新FASTA接入的结果一致性，不宣称重新证明旧全梯度逐值相等。原13428648的五窗已执行／三窗未执行和原闭窗四个边界不可行结果仍是各自原空间中的事实。

原序号1、2、5分别补齐实际1738bp定向前缀、1838bp定向后缀、862bp定向前缀。新增检查的rna-NM_122723.4／gene-AT5G28400、rna-NM_105396.4／gene-AT1G67280和rna-XM_015762904.3／gene-LOC4351603依次为9、9、16个CDS exon，全部几何失败及直接DNA失败为空，实际carry编码判定和完整语法合法性均为true。原35条已有DNA判定经确切DNA一致性转移，新增三条取得实际DNA支持；没有把历史metadata_complete当作DNA检查，也没有逐条重跑那35条。新增三窗结果为：

| 原序号与core窗口 | context bp | logZ all | logZ core observed | NLL | 外部片段梯度L1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 NC_003076.8 − 10321920 | 26314 | 1359.807540773 | 173.311464656 | 1186.496076117 | 前缀27.046406701 |
| 2 NC_003070.9 − 25190400 | 26414 | 1274.504895490 | 113.973101281 | 1160.531794209 | 后缀94.096705450 |
| 5 NC_089046.1 − 3416064 | 25438 | 1384.540933983 | 64.729884088 | 1319.811049895 | 前缀59.027576361 |

八窗全部context×15 FP32梯度finite；共同row-shift最大残差1.7881393432617188e−7≤1e−5，最大nonBG单角色后验差0.9999945163726807、最大covered后验差1.0000001192092896，均满足原1e−6界。新增三个片段的NLL梯度实际非零，说明外部真实DNA位置上的角色势实际参与本次观察目标的反向传播，不是detach消息；原首窗390bp前缀梯度L1仍为81.71318615517703。NLL是未加权未归一化的零势集合计数负对数概率，不能按这些窗口数值对模型或物种性能排名。

TRAIN metadata加载3.518352秒、原FASTA读取18.007296秒，模块导入后总76.545648秒，进程峰值1,172,784 KiB。首窗前向4.580897秒／反向4.507573秒含冷编译，三新窗分别前向1.068969／1.032280／1.003477秒、反向4.003704／3.908112／3.794984秒；这不是训练吞吐或GPU可行性证明。前向／反向nopython signatures各0→1。

summary的all_original_context_NLLs_executed为true、context_NLL_and_gradient_windows及contexts_written均为8，optimizer_steps_this_run为0。raw_FASTA_opened、unindexed_FASTA_streamed_with_TRAIN_record_gate及nonTRAIN_FASTA_lines_streamed_but_discarded为true；nonTRAIN_DNA_parsed_retained_or_used、split_val_or_test_files_opened、census_metadata_called、FASTA_index_created及DEV／独立物种／A33使用均为false。原R4规则／观察行／primary集合未改，模型、权重、optimizer、feature前向、checkpoint、训练／推理caller均未接入；边界先验、训练context／padding及归一化／importance／fit系数仍未选择。

本轮补齐了固定八窗的真实DNA支持和零势似然反向传播，原三窗缺序列不再阻塞本次诊断；不能外推到全1491窗训练支持或完整gene性能。下一步转向不依赖注释边界的训练／raw推理上下文合同及实际同方向角色势，不继续重复零势正控。原GEN使用6bp token和6144bp block，新增1738及862bp前缀相对原core token格均余4，不能将原core cache、反向邻窗feature简单拼接或翻转当作同方向全context表征；需由同一实际定向FASTA和明确token格重新建立势的来源。原A44停止条件与所有冻结数据边界不变，论文目标仍ACTIVE／NOT ACHIEVED，本轮是进展而非新模型性能或科研PASS。

## 原 TRAIN 原生 GLM 与完整语法接入准备

本轮不重跑零势正控，改为同一八窗的真实原生GLM、fresh seed0 C0十五角色势与完整语法NLL反向传播。新诊断候选view仅依赖原core坐标、原24576bp宽度／12288bp步长和原FASTA染色体长度：core两侧各一原步长邻窗，起点按原始长度裁剪、去重，内部按已存在context_segments的nearest-containing-center／较低genomic start tie逐碱基归属；不得用CDS边界、标签或先前诊断补齐坐标选view。此一次固定view是接入准备，不是已经批准的拟合上下文、padding或全染色体raw性能合同。

新full_grammar_context_view只增加这个坐标函数及与原逐碱基归属相同的可微拼接。负链只反转片段顺序，局部行仍为各自同方向RC，十五类列／phase不翻转。四个人工正负链／不等长重叠检查将新拼接逐值对照原stitch_contig_emissions，并要求拥有行VJP为1、未拥有行为0；它们随同一Slurm接入作业运行，不另起数学作业。原raw stitch有no_grad装饰，不能直接用于训练；新拼接使用切片和torch.cat，保留全部被选行的图。

编码器沿已有cache_full_train.load_native／annotate_query.encode_window：原发布GENERanno-eukaryote-1.2b-cds-annotator-preview、固定revision b0483c23b6b63787b61a6d3a204a9b517d6ba345、离线现有权重、BF16隐藏后FP32、每窗四个6144bp块及原6bp token，不更换模型、权重或tokenizer。每个邻窗由实际同方向FASTA重新编码；原native函数返回inference tensor，在其作用域之外clone成普通无梯度输入供可训练projection保存，不把反向邻窗cache翻转替代GEN。仅以原八core缓存核对实际DNA和geometry，不用cache feature或旧标签；原build_targets重建全部fine／coarse／unknown行。

头仍是原M28Core C0的fresh seed0／366191参数／FP32共享project、双向GRU、DNA conv、fuse和十五类segmentation，不读任何监督checkpoint，不继承A44或恢复其失败链。无optimizer、更新、checkpoint、局部CE、importance或损失归一系数；本轮只运行实际神经势下的未加权未归一化logZ_all−logZ_core_observed，cap128、单位角色／null0、FP64累积、checkpoint256及原NLL梯度界保持不变。非零神经势的observed logZ允许为负，只要求两配分finite及NLL≥−1e−8，不能沿用零势计数的全部非负要求。

八窗必须全部有限且保留原观察位置数，完整context势梯度finite；逐window完整NN输出的VJP应与所拥有context行逐值相同，未拥有及padding行严格零；五个原共享参数组梯度均finite且非零，fresh模型逐值未改，发布GLM冻结无梯度。原1491窗只作同一view的几何及完整primary span是否超出上下文的metadata census，后者不用于扩大view或删标签、不冒充全cohort NLL支持；发现外部primary将保留为下一步全cohort支持核实的警告。

独立根outputs/M28-FULL-GRAMMAR-GLM-CONTEXT-PREFLIGHT、唯一日志M28FULLGLMCTX_%j，沿原原生GPU smoke的generanno／1RTX3090 24GiB／2CPU16GiB／10分钟／no-requeue。八个core最多24个原生窗口／96块，完整context通常49152bp；已完成的八窗零势122秒及原原生10分钟smoke路由只用于操作余量，不保证本次吞吐。只保存started、八行数值／梯度统计和summary，小于1MiB，无features／势／研究checkpoint落盘。提交前完成五项活调度、实际资源／预约／账户80上限、唯一根日志及明确四个文件比对；实现错误且无更新最多同输入／资源／容差另根R1一次，真实数据／数学矛盾、TIMEOUT／OOM或第二次失败停止。所有DEV／独立物种／A33和原A44继续条件不变。

15:34:37 UTC五项活调度及预约核对完成；private的gpu034仍有一块3090未分配、gpu035八块未分配，shared有普通Resources／Priority等待，无本短任务的维护冲突，按原优先及近期可完成性选private。账户展开七项，加本次八项小于80，三个新实现路径、唯一输出根及日志均为空；两Python AST、SBATCH bash通过，原报告91,511 bytes与canonical逐字一致后完整保留。仅同步新模块、runner、SBATCH和原报告追加段，不复制框架或改旧科研结果。

15:36:38 UTC精确四文件canonical逐字相符后单次 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-GLM-CONTEXT-PREFLIGHT.sbatch` 返回13429053；即时核对同一job在gpu034实际RUNNING，elapsed45秒，1RTX3090／2CPU16GiB10分钟／no-requeue／Restarts0，Command／SubmitLine／WorkDir及唯一日志吻合，JOBID与STATUS RUNNING已写入，stdout／stderr当时均空。该状态不是原生模型／梯度已完成的证据，须以终态及八条实际结果判定。

## 原 TRAIN 原生 GLM 与完整语法实际完成

13429053实际15:36:39–15:40:33 UTC运行234秒，root／batch／extern均COMPLETED0:0；实际原private／generanno／gpu034／1RTX3090／2CPU16GiB10分钟／no-requeue／Restarts0，batch MaxRSS3,741,160K。单次执行，没有R1、容差或输入更改。JOBID9 bytes、STATUS10 bytes、started.json2,053 bytes、windows.jsonl17,014 bytes、summary.json35,352 bytes非空，总54,438 bytes；逐窗八行与summary.records完全相等，stdout十七份JSON为八个计算开始、八个逐窗结果和一份终summary，逐项一致，stderr为空。

实际同向原生GLM编码24个窗口／96个6144bp块，BF16隐藏转FP32的features均为[1,4096,2048]、requires_grad false、inference-tensor false。八个49152bp上下文均保留原core定向[12288,36864)、全部fine／coarse／unknown位置数，并与原core缓存DNA逐字符一致；上下文由FASTA长度／原坐标选择，不读CDS来扩范围。四项正负链／不等长人工归属检查的forward和VJP与原归属逐值一致；实际24个NN输出的拥有行VJP与全context梯度逐值相同，未拥有及padding行严格零。

八窗全部配分／NLL及完整角色势梯度finite，NLL均非负；common-row-shift最大残差1.2079385669494513e−7≤1e−5，最大nonBG角色后验差1.0、最大covered后验差1.0000001192092896，满足原1e−6界。五个共享参数组全部finite且非零，以下为八窗实际未归一化NLL梯度L2范围，不能据此后验选loss系数：project 1,061,015.377–1,316,968.210；context GRU179,766.490–239,475.376；DNA conv11,692.865–14,632.646；fuse261,712.782–319,669.006；segmentation242,900.179–292,838.493。fresh模型逐值未改，发布GLM冻结且无梯度，optimizer_steps_this_run为0。

| 原序号 | 原生神经势未归一化NLL | core外前缀角色梯度L1 | core外后缀角色梯度L1 |
| --- | ---: | ---: | ---: |
| 0 | 28911.646491639 | 38.409569415 | 143.586521708 |
| 1 | 28486.503216370 | 18.539315820 | 81.411205977 |
| 2 | 26442.393240467 | 77.647574498 | 0.425614346 |
| 3 | 26942.070861000 | 76.365699265 | 199.867386395 |
| 4 | 28024.428195249 | 56.867343790 | 75.648935499 |
| 5 | 25753.901364036 | 83.536488323 | 63.356120284 |
| 6 | 25773.453862641 | 22.027971773 | 178.619517610 |
| 7 | 29298.241752393 | 262.243891228 | 228.416875132 |

metadata加载7.388726秒、原TRAIN FASTA读取14.994130秒、原生encoder加载9.048679秒，模块导入后总177.430026秒；进程峰值7,693,604 KiB与Slurm采样batch MaxRSS分开记录。GPU峰值allocated2,714,082,816 bytes／reserved2,829,058,048 bytes。首窗配分前向5.754757秒／反向15.128694秒含冷编译，其他窗前向2.134915–2.290956秒／反向12.551722–13.215940秒；这是一次真实接入的CPU语法与GPU神经头组合资源证据，不是已更新训练吞吐保证或模型性能比较。

1491窗几何census实际Arabidopsis745窗／2234原生窗口、水稻746窗／2237原生窗口，合计4471原生窗口。完整primary span超出本候选view的有Arabidopsis2窗、水稻4窗，六个window–primary配对对应五个独特transcript（拟南芥rna-NM_001333295.1在两窗出现）。它们不是八个已检验实例，metadata超界本身不证明core事件不可行，亦不能据此删窗或扩大view。summary的full_cohort_observation_support_proven为false；本轮只取得八个原生神经实例的完整梯度接入证据。

actual_released_GLM_forward、actual_fresh_neural_forward、all_eight_NLL_role_and_shared_parameter_gradients_completed及new_model_weights_unchanged为true。未构造optimizer、读取监督checkpoint、保存研究权重／features／势、更改R4／删行／删primary、选择边界先验／批准拟合context／归一化／importance／fit系数、使用候选池预测或DEV／独立物种／A33、重开A44；其原失败与停止条件保留。本轮推进实际模型目标的可训练接入，不是训练成功、完整gene性能或论文PASS。

## 固定上下文六个超界案例支持核实准备

下一步只检验上述全部六个实际warning的原core角色事件，绑定13429053的完整实际来源、同一候选view及原1536draw／1491窗，不重跑八个神经正控、不另选样本、不从CDS扩大上下文。原scope仍两物种六条TRAIN染色体；原cache只核对这六个实际core DNA／geometry，原R4 builder保留全部观察。metadata的完整primary超界与core事件是否有任一合法补全分开记录；有限支持也不等于原超界完整链可在该view表示。

全零有限势仅用于这六个新实例的支持存在性：实际同向DNA、cap128、单位角色／null0、FP64和checkpoint256及同core允许角色不变，只计算观察配分，不重复无条件配分、反向或模型。任意有限神经势改变质量但不改变同一图的支持；因此zero-potential finite／InfeasibleObservation能回答本问题，但不能替代模型训练／完整链性能。不可行将作为诊断结果保留，不修label、不填DNA、不删行或事后改固定view。没有无条件NLL、参数、checkpoint或权重访问。

入口M28-FULL-GRAMMAR-CONTEXT-EXTENT-WARNINGS.sbatch、独立根outputs/M28-FULL-GRAMMAR-CONTEXT-EXTENT-WARNINGS、唯一M28FULL6SUP_%j日志；沿generanno／private／0GPU2CPU8GiB5分钟／no-requeue。原CPU raw上下文进程峰值约1.2GiB，新任务不加载GLM参数且仅六次观察forward，8GiB／5分钟为固定余量，不延长旧任务。输出started、六行windows和summary小于1MiB。提交前重新核对五项活调度、账户上限、根日志及三个明确文件；仅新runner明确实现错误允许同输入／资源／容差另根R1一次，数据／数值漂移、TIMEOUT／OOM或第二次失败停止。任何结果都不放行全1491窗支持、首训、raw／DEV／独立物种或A44下游，论文目标保持ACTIVE／NOT ACHIEVED。

15:47:15 UTC五项活调度／预约核对完成，无本次短CPU任务的维护冲突；账户展开七项，新增一项不超过80，两个实现路径、输出根及唯一日志均不存在。新runner AST及SBATCH bash通过，原报告96,338 bytes与canonical逐字相同后完整保留，再同步新runner、SBATCH与原报告追加结果及六窗诊断方案，不改变已完成神经job的输入或源码。

## 固定上下文六个超界案例实际完成

15:50:47 UTC单次提交 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-CONTEXT-EXTENT-WARNINGS.sbatch` 返回13429243；实际15:50:47–15:51:39 UTC运行52秒，root／batch／extern均COMPLETED0:0。原private／generanno／gpu034／0GPU2CPU8GiB5分钟／no-requeue／Restarts0，Command、SubmitLine、WorkDir及唯一日志相符；batch MaxRSS1,123,296K。没有R1、重提、扩大view、修改观察或更改容差。

JOBID9 bytes、STATUS COMPLETED10 bytes、started.json2,624 bytes、windows.jsonl4,691 bytes、summary.json11,160 bytes共18,494 bytes，五个预期产物均非空，无额外研究产物。逐窗六行与summary.records逐值相同，ID／顺序／genomic context／超界primary IDs与13429053全部六个warning及started完全相符；stdout十三份JSON为六个开始事件、六个逐窗结果及一份终summary，全部一致，stderr为空。六个原core与对应原缓存DNA逐字符相同，core全部fine／coarse／unknown行保留，六份缓存只用于DNA／geometry核对，没有使用旧feature或旧标签。

六个固定49,152bp实际TRAIN DNA上下文的原core观察配分全部有限，均为FEASIBLE_CORE_EVENT，无InfeasibleObservation；每个定向core仍为[12288,36864)。表中数值是全零势允许补全集合的log计数，不是神经NLL或模型性能。

| 原序号与core窗口 | fine / coarse / unknown位置数 | 零势observed log计数 | 完整超界primary CDS span |
| --- | ---: | ---: | --- |
| 0 Arabidopsis NC_003070.9 − 17719296 | 24401 / 4 / 171 | 1055.356117446 | rna-NM_001333295.1 [17732581,17758194) |
| 1 Arabidopsis NC_003070.9 − 17756160 | 24448 / 71 / 57 | 1196.831917061 | rna-NM_001333295.1 [17732581,17758194) |
| 2 rice NC_089043.1 + 13443072 | 23267 / 1309 / 0 | 1117.461097297 | rna-XM_066303951.1 [13429657,13444381) |
| 3 rice NC_089043.1 − 8232960 | 24576 / 0 / 0 | 1287.898363428 | rna-XM_066304365.1 [8213178,8236894) |
| 4 rice NC_089044.1 + 3452928 | 24576 / 0 / 0 | 1191.679102996 | rna-XM_015759047.3 [3430159,3458698) |
| 5 rice NC_089046.1 − 12976128 | 24576 / 0 / 0 | 1205.240045403 | rna-XM_015763105.3 [12996390,13013965) |

模块导入后总35.927036秒；首窗6.414135秒含冷编译，其他五窗1.304115–2.200012秒。这里只执行六个原core事件的观察配分，没有无条件配分、反向、模型／backbone前向、权重访问或参数更新。原1536draw／1491唯一TRAIN窗和两物种六条TRAIN记录保持不变；读取原FASTA时其他记录行流经I/O后丢弃，未解析、保留或用于非TRAIN DNA。没有DEV／独立物种／A33／A44执行链访问。

结论是六个超界warning没有造成这六个固定view内的core角色事件支持为零；不能把它们当作需要删窗、删primary或按CDS扩范围的证据。相反，每条原完整超界primary仍有CDS在候选view外，complete_outside_primary_chain_representable_in_view明确为false；有限core支持只证明存在某个合法补全，不证明该补全就是完整原primary。八个原生神经接入实例与这六个warning的支持证据均已取得，但不能外推其余原训练窗，full_cohort_observation_support_proven和complete_original_primary_coverage_proven仍为false。

下一步保留候选坐标view及全部原1491窗，转向全cohort原观察支持核实，不重跑八个神经接入或六个已完成warning、不按结果改标签或重选窗口。正式拟合之前仍需明确loss归一化／抽样importance与拟合上下文合同；全染色体raw完整链解码与训练core事件不能混为同一验收条件。本轮两个作业推进可训练接入及实际支持诊断，未执行首训，不是新模型性能或论文PASS。原A44 TIMEOUT／下游未执行事实保持关闭，论文目标仍ACTIVE／NOT ACHIEVED。

## 全部原训练窗口观察支持核实准备

本轮覆盖原1536draw中的全部1491唯一TRAIN窗，不重采样、不删窗。八个CPU分片按原window插入顺序ordinal modulo8划分，每片186或187窗；其中1477个新实例只计算固定原生坐标候选view中的零势观察配分，13429053的八个有限神经观察配分及13429243的六个有限零势配分转移支持存在性证据，不再次运行这些十四个配分。八个神经来源的质量不是零势log计数，逐行明确区别来源，绝不混合配分数值作为性能指标。读取来源实际完成标记、原配置／ID及每窗view／观察位置数；原FASTA、R4目录和builder保持原冻结输入。

每个core均与其原缓存DNA逐碱基核对，全部原fine／coarse／unknown行由原R4 builder重建，原primary集合及超界IDs保留。缓存会整体读取以核对DNA／geometry，但其中feature和旧标签不作为本轮输入。上下文仅由原坐标／24576bp宽度／12288bp步长与原TRAIN长度决定，正负链core偏移分别映射，不按注释边界扩大。沿cap128、单位角色／null0、FP64累积及checkpoint256，只执行观察forward；真实InfeasibleObservation保留为INFEASIBLE_CORE_EVENT并继续其余窗口，不是实现失败、删标签或重新调view的许可。

新入口M28-FULL-GRAMMAR-TRAIN-CORE-SUPPORT.sbatch为generanno／private／0GPU2CPU8GiB／10分钟／array0-7%8／no-requeue。六个已完成新支持forward首窗含冷编译6.414秒、其余1.304–2.200秒；每片不超过187窗，十分钟覆盖约412秒forward、原FASTA／缓存读取和导入余量，不保证实际吞吐。新独立根outputs/M28-FULL-GRAMMAR-TRAIN-CORE-SUPPORT含八个shard_0..7，各自JOBID、STATUS、started、windows及summary。第二个入口M28-FULL-GRAMMAR-TRAIN-CORE-SUPPORT-RECONCILE.sbatch只在整个array afterok后运行，generanno／private／0GPU2CPU4GiB／2分钟／no-requeue；它要求八个实际完成标记、逐行与summary一致、1491窗各出现一次、1477新计算／8神经复用／6零势复用，并输出总summary及两物种支持分母。不读模型、预测、DEV或新域；总小型产物预计小于5MiB，不保存序列、feature、势、梯度或checkpoint。

只有全1491个原core事件都实际取得有限支持，full_cohort_observation_support_proven才为true；完整超界primary覆盖、首训、raw完整gene性能和论文PASS仍是不同结论。出现不可行时保留逐窗证据并判断原core目标的真实冲突，不重新提交以改变科研结果。明确实现错误最多同输入／资源／容差另根最小R1一次，实际数据／数学漂移、TIMEOUT／OOM或第二次失败停止，不修复或延长旧A44链。

16:02:59 UTC五项活调度及预约核对完成，private短CPU任务无维护冲突；账户展开七项，加八个array任务及一个汇总后十六项小于80，array8小于100。两个private节点可由Slurm安排CPU，不占GPU；唯一根、三个实现路径及两类日志均不存在。新runner AST及两SBATCH bash通过，报告原106,557 bytes与canonical逐字相同并保留，仅同步本runner、两个SBATCH及此报告追加段。

16:05:04 UTC单次 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-TRAIN-CORE-SUPPORT.sbatch` 返回13429430，八个任务均16:05:06 UTC在gpu034实际RUNNING，array throttle8／0GPU2CPU8GiB10分钟／no-requeue／Restarts0及Command／SubmitLine／WorkDir／唯一日志吻合。16:05:35 UTC单次 `sbatch --parsable --dependency=afterok:13429430 sbatch/M28-FULL-GRAMMAR-TRAIN-CORE-SUPPORT-RECONCILE.sbatch 13429430` 返回13429444；实际PENDING普通Dependency，整个13429430_*依赖、0GPU2CPU4GiB2分钟及原路径相符。四文件canonical逐字一致后提交，没有改动运行中的实现。

## 原训练观察不可行案例定位准备

运行中已有834条逐窗结果时发现13个INFEASIBLE_CORE_EVENT；这是部分集合的实际负例，不是全1491的最终分母或实现失败。最先八个负例的原完整primary span均在固定view内，不能仅凭先前六个超界警告解释它们。保持全array继续完成所有窗口；不进入把这些零支持事件直接作有限NLL的拟合，也不删除观察或改R4规则。

新增定位只在13429444实际完整汇总之后运行：按原ordinal选每物种／正负链的第一个不可行事件，最多四个，不择优或重复通过案例。相同actual TRAIN DNA、原core允许角色、cap128／unit／null0／FP64／checkpoint256下，使用原_state_table、_forward_messages及_advance定位第一个全部状态质量为零的碱基；用实际DNA合法transition的角色与原允许角色的空交集逐项验证。若prefix一直有支持而最后background无质量，则单独记录外侧闭合失败，不误称局部标签冲突。记录原core坐标、fine／coarse类、实际DNA、前一步状态、可行角色及该位置原primary的既有完整CDS DNA／motif检查；只在原固定范围内检验，不从注释扩上下文。局部反例不外推全部失败机制。

入口M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-TRACE.sbatch、独立根outputs/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-TRACE及唯一M28FULLFAILTRACE_%j日志，generanno／private／0GPU2CPU8GiB5分钟／no-requeue，afterok:13429444。只保存started、最多四行诊断及summary，预计小于1MiB；不读取模型、不反向／更新、不填DNA／改标签、不放行DEV／独立物种／A33／A44。明确实现错误最多同输入／资源／容差另根最小R1一次；数据／数学矛盾、TIMEOUT／OOM或第二次失败停止。前一轮1491支持筛查中的真实不可行不是可修复失败。

16:10:38 UTC重新完成五项活调度／预约核对，private短CPU作业无维护冲突；账户展开十六项，加定位十七项小于80，新runner／SBATCH／根／唯一日志均不存在。定位原始版本AST及SBATCH bash通过；补入原primary DNA检验后再次作与新编辑相称的AST语法核对。原报告109,772 bytes已精确同步并完整保留，仅转移本定位runner、SBATCH及现有报告追加段。

## 全部原训练窗口观察支持实际完成

原1491窗已经全部核实：1467个core事件有有限支持、24个零支持，不能对全集合宣称原硬观察NLL可直接拟合。13429430八个task／batch／extern全部COMPLETED0:0，共24个规定Slurm身份；八个任务16:05:06 UTC同时开始，16:10:57–16:11:15 UTC结束，elapsed351–369秒。原private／0GPU2CPU8GiB10分钟／no-requeue条件未改，无R1／延时。batch MaxRSS734,692–752,760K。13429444／batch／extern均COMPLETED0:0，16:11:16–16:11:44 UTC运行28秒，0GPU2CPU4GiB2分钟，batch MaxRSS846,508K。

八片windows与各summary.records、started身份／原ordinal modulo8逐项相同；合并后的1491行与总summary.records逐值相同。每片stdout恰为其全部逐窗结果、所有新实例的开始事件及一份终summary，汇总stdout恰一份无records的实际summary；九份stderr全部为空。44个预期小型文件均非空，总4,375,113 bytes，1491份原core缓存恰各核对一次，全部原FASTA core DNA与缓存逐字符相同，所有core fine／coarse／unknown分母保留。新观察forward1477次，其中1453有限／24不可行；原八个神经质量和原六个零势质量及ID／来源逐值复用，无再次运行其配分。

| 原TRAIN物种 | 原唯一窗口 | 有限支持 | 不可行 | 正链不可行 | 负链不可行 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 745 | 732 | 13 | 6 | 7 |
| Oryza sativa | 746 | 735 | 11 | 5 | 6 |
| 合计 | 1491 | 1467 | 24 | 11 | 13 |

全部24个反例的原完整primary span均在各自固定view内；此前六个完整primary超界案例仍全部有core支持，超界集合总数仍为六。不能把这24个反例归咎于这六个warning，亦不能以再扩大halo替代定位。full_cohort_support_evaluated=true而full_cohort_observation_support_proven=false，是完整诊断的科学负结果，不是执行失败；不能丢弃24窗后把1467窗当作原完整拟合队列。逐窗失败ID及原标签分母完整保留在总summary与各shard，不计算性能排名或改原冻结规则。

本轮没有模型／权重／优化器／反向／研究checkpoint，也未使用旧cache标签／feature、DEV／独立物种／A33或重开A44。该结果改变下一步：先定位原硬角色观察与真实DNA／完整基因语法的具体冲突，再设计保持原证据的目标，不直接把零支持事件送入有限NLL拟合。完整gene原始FASTA性能及论文目标仍未达到。

## 四个原训练不可行案例定位结果

16:15:21 UTC单次 `sbatch --parsable --dependency=afterok:13429444 sbatch/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-TRACE.sbatch` 返回13429503；实际16:15:22–16:15:50 UTC运行28秒，root／batch／extern全部COMPLETED0:0，private／generanno／gpu034／0GPU2CPU8GiB5分钟／no-requeue／Restarts0及Command／SubmitLine／WorkDir／唯一日志吻合。batch MaxRSS1,202,564K，模块导入后14.404772秒。没有R1或条件更改。

五个预期产物JOBID9 bytes、STATUS10 bytes、started.json2,984 bytes、windows.jsonl6,413 bytes、summary.json16,474 bytes非空，共25,890 bytes。四行结果与summary.records逐值一致，stdout四行诊断及一份终summary相同，stderr为空。ID／原ordinal及view与完整1491源的确定性每物种／每链首个反例相符；只打开对应四份原core缓存，原DNA、core fine／coarse／unknown行及所有观察均保留。

四例均在原core内部出现首次零质量prefix；前一步实际仍有合法状态，但原允许角色与当时所有DNA合法transition角色的交集为空，不是外侧closure才出现的问题。以下坐标均为原染色体0-based单碱基位置；角色编号沿原十五类顺序。

| 原ordinal及物种／链 | 首次冲突genomic base | 原粗标签允许角色 | 实际DNA合法角色 | 该处原primary及完整DNA判定 |
| --- | ---: | --- | --- | --- |
| 44 Arabidopsis + | NC_003070.9:10351420 | intron [4,5,6] | [3,9] | rna-NM_001084155.2／AT1G29630，非标准剪接；拼接coding carry合法 |
| 190 Arabidopsis − | NC_003070.9:10457438 | CDS组 [1,2,3,7,8,9,10,11,12,13,14] | [6] | rna-NM_102727.2／AT1G29870，非标准剪接；拼接coding carry合法 |
| 871 rice + | NC_089043.1:27347930 | intron [4,5,6] | [2,8] | rna-XM_026020294.2／LOC4347925，非标准剪接及2bp内含子；拼接coding carry合法 |
| 997 rice − | NC_089043.1:26938035 | CDS组 [1,2,3,7,8,9,10,11,12,13,14] | [0] | rna-XM_066303997.1／LOC4347870，内部终止密码子；拼接coding carry不合法 |

四例首次冲突处fine均已为IGNORE −100；原coarse仍分别为intron2或CDS1。相应原primary的CDS_complete_supported和eligible_positive原本均为false，DNA／geometry判定也不满足当前完整语法，未在诊断中改变这些字段。具体证据包括每个真实前一步state／carry、局部DNA、原CDS exon坐标、允许角色空交集及原primary拼接DNA检查。该结果不等于这些注释生物学上错误：它说明当前严格ATG／完整ORF／GT或GC–AG及min-intron5的gene-set universe与这些硬粗观察不兼容。

因此，原local_loss对不支持完整CDS的primary仍提供的粗粒度监督，不能无条件等同于完整基因语法中必须精确满足的角色事件；fine已忽略并不保证coarse事件可行。这四个真实反例把问题定位到观测建模与语法支持的适配，而非GLM梯度接入、注释halo扩展或Slurm故障。它们不覆盖其余20个失败机制，diagnostic_covers_all_failed_mechanisms=false；完整负例清单24条仍保留，不丢弃／修label，不把14个已通过实例或1467个有限实例替代原1491窗。

下一轮先核实全部24个反例的机制覆盖，再据此设计能保留原标签和全部窗口、同时与完整gene语法相容的观测目标；这不是已经选定误差率／loss系数、改变R4、放行首训或声称模型性能改进。当前论文目标仍ACTIVE／NOT ACHIEVED，原A44 TIMEOUT及未执行下游保持关闭。本轮所有新增作业已终态并由主消费，没有仅等待外部作业的阶段，固定监控heartbeat继续PAUSED，未创建新监控或重复goal等待循环。

## 全部原训练反例机制核实准备

本轮绑定已完成13429444的全部24个原TRAIN负例及13429503的四个定位，按原ordinal逐一覆盖，不另选反例或删除窗口。四个已有trace仅核对当前原DNA／允许角色及原primary直接DNA证据后复用；其余20个由同一固定cap128自动机定位首次空prefix或末端background闭合失败。全部core仍与原缓存逐碱基相同，原R4 builder的fine／coarse／unknown分母及固定view不变。重复出现的同一transcript只去重机制统计，不减少24个窗口；其DNA／geometry及原完整性字段若不一致即报告漂移。数据保留原两物种六条TRAIN染色体，不使用DEV、新域或监督权重。

入口M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-CENSUS.sbatch、独立根outputs/M28-FULL-GRAMMAR-TRAIN-SUPPORT-FAILURE-CENSUS和唯一M28FULLFAILALL_%j日志；generanno／private／0GPU2CPU8GiB5分钟／no-requeue。原四窗定位28秒，新增20次相同forward及RNA检验，没有反向、模型、optimizer或参数更新；小型started／24行／summary预计小于1MiB。明确实现错误仅同输入／资源／容差另根最小R1一次；真实数据／数学矛盾、TIMEOUT／OOM及第二次失败停止，旧负例和A44停止记录不改。

16:26:53 UTC重新核对五项活调度／预约；账户展开作业数及新根／两个实现路径／日志碰撞核实后，仅同步新runner、SBATCH和本报告追加段。新AST及SBATCH bash通过，原报告118,953 bytes与canonical逐字相同并完整保留。同期准备的观测概率仅是独立候选数学验证，不选择正式拟合误差率、loss系数或更改R4标签。

## 全反例核实的实际中断与候选观测似然准备

全反例核实实际提交13429578，16:31:06–16:31:57 UTC运行51秒，root及batch FAILED 1:0、extern COMPLETED 0:0，batch MaxRSS 1,203,600K，无重排。实际错误是既有DNA_checks的“Actual-carry and direct concatenated-CDS grammar verdicts disagree”；尚未产生完整summary，不把已输出的部分定位当作24例全部完成。先定位该具体不一致，只有明确实现错误才使用约定的一次同条件R1；若为真实数学／数据矛盾则停止这一支，不改旧结果或放宽检查。

并行准备独立的候选粗观察通道：fine仍是原精确角色，unknown观测概率为1，每个已观察coarse行在三个region类别上归一化，正确类别概率1−epsilon、另两类各epsilon/2。所有原标签及完整gene-set推理语法保持不变；这是改变观测概率模型的候选，不是宣称原硬事件已可行。epsilon仅由调用者给定，0.1和2/3是人工数学例子，未选择正式拟合参数。NLL为logZ(theta)−log(sum_G exp(S_theta(G))q(y|G))。现有自动机能量逐行减背景logit，因此用logits+logq算出的相对logZ必须补回sum_i logq_i(background)，否则似然和epsilon梯度都错误。

人工验证仅使用一／两个可穷举完整单exon ORF的16／27bp DNA，独立枚举null及所有非重叠gene子集。五项验证覆盖FP64数值与全部logit／epsilon梯度、core外观测为1、行公共平移不变、unknown单位似然、epsilon=2/3时保留每个coarse行的log3常数、原硬冲突变为候选有限事件以及中心差分。枚举／解析梯度绝对容差1e−9，差分1e−6，不读研究DNA、模型或checkpoint，无optimizer更新。入口M28-FULL-GRAMMAR-COARSE-CHANNEL-CHECKS.sbatch，独立根outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-CHECKS、唯一M28COARSECHCHECK_%j日志，generanno/private/0GPU2CPU4GiB2分钟/no-requeue。明确实现错误仅同输入／资源／容差另根R1一次；数学失败、TIMEOUT／OOM或再次失败停止，不通过更改容差或人工fixture掩盖失败。

16:32:05 UTC五项活调度及预约重新核实，账户7个展开作业、上述新根／实现／日志均无碰撞；两个Python AST及SBATCH bash已通过。只同步三份新增数学验证文件及本报告追加段，逐文件核对内容后单次提交。这仍是候选目标的数学准备，不是研究训练、模型选择或论文性能结果。

三份数学实现及报告逐文件cmp相同后单次提交13429599，实际16:37:09–16:37:20 UTC、11秒，root／batch／extern均COMPLETED 0:0，gpu034/private、0GPU2CPU4GiB2分钟、Restarts0，batch MaxRSS682,212K。五项人工验证全部通过（5 tests、0 failures、0 errors，4.6223秒），stderr是预期unittest逐项ok及OK，而非运行错误。这只验证候选观测似然及一阶梯度的数学实现，不提供真实TRAIN效果或正式拟合依据。

失败13429578已实际输出21行，下一未输出窗为原ordinal1352、rice NC_089046.1:26824704正链；剩余原ordinal1389及1416也未完成。为判断错误是否来自既有helper把拼接ORF与终止密码子必须位于单exon内的语法条件混同，准备仅检查1352窗三个原primary的逐exon实际carry、拼接CDS、原geometry和helper原错误。不重做21次forward，不改DNA／label／grammar。独立入口M28-FULL-GRAMMAR-CDS-CHECK-DISAGREEMENT.sbatch、根outputs/M28-FULL-GRAMMAR-CDS-CHECK-DISAGREEMENT、日志M28CDSDISAGREE_%j、generanno/private/0GPU2CPU4GiB2分钟/no-requeue、小型summary，失败则保留证据停止，不自动重试；五项调度、账户7作业及新路径无碰撞已重新核实，AST／bash通过。

独立诊断单次13429614实际16:39:23–16:39:43 UTC、20秒，root/batch/extern COMPLETED 0:0，Restarts0、MaxRSS1,106,884K。1352窗前两个原primary直接／carry均合法，第三个rna-XM_015764666.3原CDS_complete_supported及eligible_positive均为false：拼接984bp CDS为ATG…TGA且无内部stop，但末CDS exon仅2bp（oriented [34100,34102)），上个exon留下实际carry T，本exon是GA。因此拼接标准ORF通过，而当前严格语法明确不支持跨exon terminal STOP；原geometry也已报告split_terminal_codon。既有helper漏把这个既有条件纳入直接判定，错误地要求两者等价，并非实际自动机与其语法矛盾。一次R1只在直接报告加入split_terminal_stop后再与原ORFIndex比较，保留该primary不合法结论及所有原标签／语法／负例；源helper和FAILED根不覆盖，沿相同24窗、cap128、资源和容差，在独立R1根执行，第二次失败停止。

同时准备原24负例的候选似然／一阶梯度真实TRAIN诊断，源为实际13429444全1491窗及13429599五项数学通过。仅zero CPU FP64 potential、诊断epsilon0.1、原固定view及原R4全部fine/coarse/unknown；无GLM或C0、无权重读取或更新。逐窗核对core原DNA、view／观察行数，计算全／候选观测logZ、背景补偿、NLL、所有role梯度及epsilon梯度；检查有限性、概率非负、每行公共平移梯度和epsilon概率导数界。若fine硬事件仍不可行，明确记录STILL_INFEASIBLE_EXACT_FINE_EVENT，不能删fine或改语法。只有全部24有限后，才能结合原1467硬可行事件讨论候选的数学支持，不等于全1491的数值NLL、神经梯度接入、正式fit或性能提升。

入口M28-FULL-GRAMMAR-COARSE-CHANNEL-TRAIN-PROBE.sbatch，独立根outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-TRAIN-PROBE、唯一M28COARSECHTRAIN_%j日志，generanno/private/0GPU2CPU8GiB5分钟/no-requeue；24次双forward及一阶reverse预计数分钟，小型summary／24行小于1MiB。实现错误允许同条件另根R1一次，数学／数据问题、TIMEOUT／OOM或再次失败停止。两项新入口的五项调度、账户7作业及路径／日志无碰撞已核实；仅同步四份新增实现和本报告，逐文件cmp，两个Python AST及SBATCH bash通过后单次各提交，不恢复A44。

## 全部硬观察反例机制实际结果

上述五文件逐字cmp相同后并行单次提交，Slurm实际身份为候选TRAIN诊断13429625、一次机制R1核实13429626。两项16:42:16 UTC开始，gpu034/private，各0GPU2CPU8GiB5分钟、Restarts0；13429626于16:43:09结束，53秒，root/batch/extern COMPLETED 0:0，batch MaxRSS965,912K。全部24行、started配置／ID及stdout与summary逐值相同，stderr为空；五个预期文件JOBID9、STATUS10、started6,895、windows40,973、summary87,683 bytes非空，共135,570 bytes。FAILED源的前21例定位除耗时和新job来源外逐值相同，四例原trace复用；修正仅使实际跨exon terminal STOP的直接检查与既有语法正确对齐。

24例均在原core首次失去支持，无仅外侧closure失效。20例是原观察排除全部DNA合法转移，4例是此前被硬观察约束的状态在当前DNA上已经没有合法转移（两例非ATG起始、一次N位点、一次跨exon terminal STOP）。首次fine为IGNORE的23例、fine背景0的一例；首次coarse为intron2的14例、CDS1的8例、background0的一例、unknown−100的一例。故“所有首次冲突都是非完整primary的粗标签”不能推广到全部24例；背景冲突原ordinal202、AT负链NC_003070.9:12575890，fine/coarse均0，而实际前一步被迫的coding只能发出[1,7]，该位置没有primary跨越。

首次冲突所涉及原primary共23个去重transcript；直接DNA失败计数非标准splice19、非ATG2、ambiguous CDS1、internal STOP1、split terminal STOP1、无terminal STOP1（类别可重叠）；geometry有min-intron不足3例及split terminal codon1例。保留这些原primary的既有完整性及eligible字段，不能由诊断宣布注释生物学上错误，亦不删除24硬失败。原1467有限／24不可行硬观察群体结论仍不变。

## 候选神经目标接入准备

在候选24窗实际完整终态、所有NLL／一阶梯度有限之后，拟对13429503既定的每物种／每链首个四个硬负例接入原native GLM→fresh seed0 C0→base-owned完整语法。只将loss换为已验证候选观测概率，沿原raw view、原标签及strict inference grammar；epsilon0.1仍是诊断值。检查四窗真实神经NLL、所有role及epsilon梯度、原五组shared／segmentation参数梯度、原owned行VJP，以及encoder冻结和fresh C0 state逐值未改变；无optimizer、checkpoint或正式训练参数选择。入口preflight_full_grammar_coarse_channel_glm.py及M28-FULL-GRAMMAR-COARSE-CHANNEL-GLM-PREFLIGHT.sbatch已本地AST／bash通过，尚未同步／提交。预期private/1×3090至少24GiB、2CPU16GiB5分钟/no-requeue、独立输出根及M28COARSECHGLM_%j日志；此前原八窗native诊断实际177.43秒、encoder加载9.05秒，四窗预计两分钟左右加少量余量，输出小于1MiB。提交前必须确认上游完整条件及新活调度／路径；不以仍运行或部分结果代替这些条件。

## 全部原硬负例的候选似然实际结果

13429625实际16:42:16–16:47:15 UTC、4分59秒，root/batch/extern全部COMPLETED 0:0，MaxRSS1,138,152K、Restarts0，五分钟限时内正常结束，没有取消、延时或重提。五个预期文件JOBID9、STATUS10、started4,581、windows37,014、summary54,435 bytes非空，共96,049 bytes；24行与summary.records、配置／ID、stdout逐值一致，stderr为空。实际原ordinal及ID恰为全1491源中的全部24个硬不可行事件，13个AT／11个rice，无选择或丢弃。

24/24候选NLL、所有context角色梯度及epsilon梯度均有限，role及epsilon梯度均非零，仍不可行0。保留fine493,051、coarse89,691、unknown7,082行，共589,824个原core碱基；所有标签逐值未改、无coarse行删除，strict gene-set inference grammar未改。NLL范围1,115.3624–2,601.5200，行公共平移梯度最大5.621e−16，epsilon梯度范围−620.00001至14,369.41381；这些是未训练zero potential的目标数值，不是模型性能或可比较的准确率。包括首次fine背景冲突的一例也有限：候选概率放松此前的粗观察约束，不需改该背景fine标签。完整核算after-import283.0674秒；今后同类24窗完整反向不能照搬过紧五分钟估计，应依据此次实际耗时及运行余量准备新执行计划，不修改本次记录。

数学推论的范围明确：原1467个硬可行事件各存在满足fine及coarse的完整gene set；候选在任意0<epsilon<1时给每个已观察coarse类别严格正概率，因此这些事件仍有正支持。加上本次其余24个实际有限事件，原1491个core事件在同一固定raw-coordinate view及当前strict grammar中均有候选正支持。这里没有重跑1467个正例的候选数值NLL，不能称为全1491数值／神经执行完成；也不证明六个完整primary超界warning已解决、整染色体协议成立、正式fit可放行、泛化或论文PASS。原硬观察24失败仍保留为原目标的负结果。

候选native接入的上游24/24条件现已满足。16:48 UTC重新核对五项活调度／预约、账户7个展开作业及新实现／输出／日志无碰撞；Slurm确认gpu034为8×3090／24G每卡、已用7，gpu035为同规格8卡全空闲，shared已有GPU等待队列，故沿private/3090而不改科学输入。仅同步两份新增native实现及本报告后逐字cmp并单次提交。明确实现错误最多同输入／资源／容差独立R1一次，数据／数学问题、TIMEOUT／OOM或再次失败停止；不选择epsilon、normalizer、importance或fit coefficients，不恢复A44。

## 候选神经接口实际提交

两份native实现及报告与canonical逐字cmp相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-COARSE-CHANNEL-GLM-PREFLIGHT.sbatch返回13429660。实际16:50:07 UTC开始，Slurm确认M28COARSECHGLM、gpu034/private、1×3090/2CPU16GiB5分钟、Restarts0、独立命令／WorkDir／M28COARSECHGLM_13429660日志身份。首次检查为RUNNING 17秒，不把提交、运行或部分窗口代替四窗及fresh state／encoder冻结的完整终态核实。本轮没有首训／DEV／新物种／A33或A44重启，论文目标保持ACTIVE／NOT ACHIEVED；固定监控仍复用原对象，不创建新对话或重复心跳。

## 候选神经接口实际结果

13429660实际16:50:07–16:51:56 UTC、1分49秒，root/batch/extern全部COMPLETED 0:0、Restarts0，batch MaxRSS3,729,876K。四个预定负例全部取得有限native NLL、所有角色及epsilon梯度、project/context/dna/fuse/segmentation五组参数非零有限梯度，四个owned-window行VJP均逐值相同。released GLM仍冻结且无参数梯度，fresh seed0 C0完整state逐值未变，未构造optimizer、更新参数或写checkpoint。五个预期小型文件JOBID9、STATUS10、started5,556、windows7,270、summary20,377 bytes非空，共33,222 bytes；四行、配置／ID、原四个选择ID、stdout与summary逐值一致，stderr为空。

native NLL范围24,744.1927–29,984.7219，FP32角色公共行平移梯度最大1.262e−7、角色导数绝对值最大1；epsilon导数2,103.9218–5,465.7035。after-import91.3162秒，GPU峰值allocated2,714,082,816／reserved2,829,058,048 bytes。未归一化fresh C0的project梯度L2约1.02e6–1.37e6，有限及接通不等于正式优化步长已经合适；epsilon0.1、normalizer、importance和fit coefficients仍未选择。这里也不能拿native NLL与zero-potential NLL跨窗口或跨口径作性能排名。

本轮实际完成三类证据：五项独立人工似然／梯度检查、原24硬负例全部候选NLL及一阶梯度、四个既定负例的native网络反向接入。它们支持候选观测目标保留全部标签／窗口且无需改变strict inference grammar即可消除固定诊断view的零支持问题；不证明参数已经学到更好的gene注释。下一轮优先落实可执行拟合协议及原队列训练路径，明确正式观测参数、归一化／采样和完整raw协议的剩余边界；没有新改动、失败或具体疑点，不重复本轮已通过的人工／24窗／四窗检查。原24硬观察负结果、FAILED13429578和A44停止记录完整保留，DEV／新域／A33仍未使用。论文目标保持ACTIVE／NOT ACHIEVED；本轮全部新增作业已终态并由主消费，固定heartbeat保持PAUSED，没有只等待外部作业或重复goal等待循环。

## 原队列首轮拟合入口与元数据准备

本轮新增[首轮拟合草案](full_grammar_first_fit_contract.md)及独立prepare／fit入口，不重跑已经完成的数学、24窗候选或四窗native检查。live sacct重新确认13429625、13429626和13429660均COMPLETED 0:0；原报告135,096-byte前缀与canonical逐字相同后续写，历史证据未改。

新草案只使用原六条TRAIN染色体及原1536 draws／1491 unique，连续三遍原顺序共4608次更新，全部24个原硬负例保留。fine、coarse、unknown及固定raw-coordinate view不变，但明确改变观测概率模型：拟固定epsilon0.1，不学习噪声率、不扫参数或根据DEV选择。每draw使用原importance × 完整集合NLL / max(1,fine+coarse core行数)，分母沿实际原local_loss的计数规则；AdamW3e−4／wd0.01／clip1沿原设置，不增加辅助loss。此为可审查科学选择，不是旧硬观察已PASS、机制公平对照或已经批准的正式fit。

fit_full_grammar_coarse_channel.py prepare仅核实已完成1491／24／4窗来源与原TRAIN元数据，生成1536行原draw计划，统计已知观察分母、importance multiplier及原负例曝光；不读raw FASTA、feature tensor、模型或checkpoint，不构造optimizer。新CPU入口M28-FULL-GRAMMAR-COARSE-CHANNEL-FIT-PREPARATION.sbatch为generanno/private/0GPU2CPU4GiB2分钟/no-requeue；新根及M28COARSEFITPREP_%j日志无碰撞后单次提交。明确实现错误最多同科学输入／资源／容差另根R1一次；数据／数学矛盾、TIMEOUT／OOM或第二次失败停止。

另一个FIRST-FIT入口仅准备而未提交：generanno/private/1×3090至少24GiB、2CPU16GiB30小时/no-requeue，三遍checkpoint和新fit文件上限64MiB，无自动afterok、续训或旧A44重试。约21小时只来自四窗native稳态的外推，不是全队列实测。完整染色体raw caller尚未交付，六个primary超出view的warning仍存在；本轮有限上下文core事件的拟合不能宣称整染色体联合MLE或模型性能／泛化成立。

提交准备前五项活调度／预约检查完成，账户展开7项低于80；旧A44两项依赖未满足作业仍保留而不取消。新实现／文档、两个唯一输出根及日志检查无碰撞。AST／bash只验证新增入口语法，实际CPU准备将检出来源字段、原抽样顺序、view和观察计数漂移；若失败则不启动fit。此准备产物不推定用户首训授权。DEV、独立物种、Setaria、A33和对外发布范围均未扩大，goal保持ACTIVE／NOT ACHIEVED。

## 原队列拟合元数据实际结果

五份新增／续写文件与canonical逐字相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIT-PREPARATION.sbatch返回13429756。实际17:08:26–17:08:42 UTC运行16秒，Slurm root/batch/extern均COMPLETED 0:0，Restarts0；gpu034/private、0GPU2CPU4GiB2分钟与独立命令／工作目录／日志身份一致，batch MaxRSS848884K。元数据核算after-import2.340927秒，无raw FASTA／features／模型读取，无optimizer或checkpoint。

1536行training_plan与原paired_training_draws逐行核对ID、原顺序、物种、stratum和importance；1491个unique与原完整支持来源逐值核对view、core坐标、fine/coarse/unknown及硬支持状态。原两物种各768次曝光、unique745／746、24个原负例曝光完整保留。拟南芥importance sum764.0790550，fine曝光18317883／coarse106068／unknown450417；水稻importance sum749.0897519，fine18676744／coarse76120／unknown121504。已知观察分母范围1–24576，拟议importance／分母范围1.120891758e−5–1.333333333；不同importance并未重归一或改抽样。四窗稳态最大16.25516194秒外推4608步20.806607小时，不是全cohort实测。

发现两次零已知观察draw（原序号612及640）重复同一拟南芥负链窗口NC_003070.9|14770176：fine=coarse=0、unknown24576、原importance4/3、分母1。它们不删除，原观察NLL／梯度为零；原顺序AdamW仍执行，动量／weight decay可能作用，不能声称每个optimizer步骤都有观察学习信号。此事实已补充到同一拟合草案，未改变科学配置或另跑已通过诊断。

五个预期非空产物JOBID9、STATUS10、started1961、summary3975、training_plan732991 bytes，共738946 bytes。started／summary配置和job identity相同，1536序号连续，stdout仅一行且与summary逐值相同，stderr为空；FIRST-FIT根仍不存在。准备工程已完成，尚未提交新科学观测模型的首轮真实拟合；不把原24硬事件FAIL、六个primary extent warnings或整染色体raw caller未完成事实改为PASS。全部新增作业已终态，无只等待外部作业阶段，固定监控保持PAUSED；原目标保持ACTIVE／NOT ACHIEVED。

## 完整序列输出入口与人工接口检查

本轮推进不依赖尚待确认新首训的raw接口：annotate_full_grammar_coarse_channel.py只接收FASTA、最终checkpoint和输出位置／大小上限。核对最终epoch3／step4608、固定神经模型／released GLM revision、完整语法cap128与新观测目标身份，严格加载C0 state；不运行训练readiness，不读catalog、标签、候选池、feature cache、旧预测或A44链。已批准的原队列拟合准备13429756仍COMPLETED 0:0，FIRST-FIT根不存在；此代码不是已有训练权重或真实推理结果。

每个contig按原manifest.starts的24576bp／12288bp网格（包括移位的最后完整窗）在两链原生编码。复用stitch_contig_emissions逐碱基nearest-containing-center／lower-genomic-start归属、十五角色列不置换；调用实际compiled emission lattice在整个定向contig上联合最大化完整非重叠gene-set单位角色增益，null0、FP64、实际carry、cap128，不逐窗闭合gene后owner过滤、不加beam或阈值。两链独立解码，不宣称跨链非重叠。复用genomic_chain／as_transcript／write_gff3输出完整CDS gene／mRNA／exon／CDS，核实输入DNA、GFF phase、原motif、min-intron5及边界；仍不含UTR或all isoforms。

新检查只运行四个明确接口例，不重复旧配分／候选24窗／native四窗检查：人工30000bp完整基因跨24576bp窗，在正负链均需由whole-contig solve恢复并导出正确坐标；14bp的实际未完成codon与intronic N需导出负链CDS phase2／0；40000bp原网格最后起点15424需逐行符合独立owner计算且不交换角色列；最终checkpoint身份需拒绝中间步／旧loss／错误cap等实际可达输入。native encoder和head均为stub整数FP32势，没有released模型forward、真实checkpoint、研究FASTA或optimizer；实际Numba complete-gene-set解码及GFF3文本写入／读回会执行。通过只证明新接口接通，不能证明真实权重推理、完整染色体资源、性能或论文成功。

唯一入口M28-FULL-GRAMMAR-RAW-CHECKS.sbatch、独立根outputs/M28-FULL-GRAMMAR-RAW-CHECKS及M28FULLRAWCHK_%j日志，generanno/private/0GPU2CPU4GiB2分钟/no-requeue；人工长例仅低motif密度，2分钟含首次实际JIT余量，不能外推高密度染色体耗时。AST／bash和新根／文件／日志检查通过，五项活Slurm／预约核对完成、账户仍7个展开作业，报告139919-byte原前缀与canonical相同后续写。四项检查全通过才算本次接口准备完成；明确纯接线实现错误最多同输入／科学设置／资源／容差另根最小R1一次，数学矛盾、TIMEOUT／OOM或第二次失败停止。新首训未提交，不读取DEV／独立物种／Setaria／A33，不重开或取消A44，不创建后台监控；所有历史失败保持原样。

## 完整序列输出人工接口实际结果

五份新接口／测试／提交入口／报告与canonical逐字相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-RAW-CHECKS.sbatch返回13429790。实际17:18:22–17:18:49 UTC运行27秒，root/batch/extern全部COMPLETED 0:0、Restarts0；gpu034/private、0GPU2CPU4GiB2分钟，独立命令、工作目录与M28FULLRAWCHK_13429790日志身份一致。batch MaxRSS1038232K，after-import14.944219秒／进程peak RSS918096KiB。四项实际运行、failures0／errors0，stdout四个test均ok且最终summary逐值相同，stderr为空。

人工30000bp完整基因在正负链分别恢复CDS [3000,33000)／[3864,33864)，单位role增益300000，GFF3坐标3001–33000／3865–33864及phase0正确；没有逐窗关闭后过滤丢弃。带实际codon carry及intronic N的14bp负链例恢复两个CDS [0,5) phase2／[10,14) phase0，GFF3分别1–5／11–14。两例的相反链全zero势返回null。40000bp原网格的末窗[15424,40000)及两链全部owner行与逐碱基独立计算相同，角色列未交换。最终checkpoint身份例拒绝旧loss、中间步／epoch、错误cap／epsilon／revision；这里只验证人工state元数据，不替代实际权重文件的严格加载或完整训练产物核实。

四个非空预期产物JOBID9、STATUS10、started597、summary1145 bytes，共1761 bytes；started／summary配置和job identity相同。已完成代码／测试交付，未加载任何研究checkpoint或released GLM、未读研究FASTA、标签／catalog或feature cache、未做optimizer更新。新raw caller只在每条定向contig上联合选择，不增加跨链非重叠、UTR或all-isoform主张。旧compiled模块开头的“never imported by raw callers”说明现已不适用，核对原文件与canonical相同后仅改说明为新raw接口可复用，递归／浮点／tie／grammar代码逐行未改，无须重复算法测试。

当前下一步是实际染色体规模的可执行资源表示与已提出的首训协议确认，不能以人工低motif密度例缩小原完整gene研究目标。整染色体运行、真实checkpoint推理及模型准确率／跨物种收益尚未成立，原24硬观察FAIL、六个完整primary extent warnings及A44停止条件均保留。FIRST-FIT根仍不存在，没有重新请求常规启动／成本许可或用本次PASS推定首训授权；goal保持ACTIVE／NOT ACHIEVED。本轮新增作业已终态并由主消费，固定监控继续PAUSED，不创建重复等待循环。

## 原 endpoint 编译分支停止条件与新回溯表示

核实历史资源时发现上轮raw接线遗漏了一个有效停止条件。a43_natural_competitor_preparation.md第171–191行及实际13424814／13424849均显示：原emission_lattice_compiled endpoint／前驱arena分支首轮编译接口失败，唯一R1七项通过、一项fixture ERROR，随后修正fixture但按停止条件没有再运行；完整768bp随机矩阵验证及49kb回放未执行，分支已关闭。live sacct再次确认两作业FAILED1:0、28／25秒。13429790的四项人工接口通过是真实结果，但不解除该旧停止条件；上轮使用这一分支不能作为当前可放行raw依赖。本轮将当前raw入口移除其依赖，原失败、未执行条件及文件全部保留，不提交A43 R2或旧49kb回放。旧模块仅修正文档说明，不改算法或以新名字重提同一原型。

新full_grammar_viterbi.py是真正不同的逐碱基自动机最大分数／回溯表示，复用已实际13428495验证的完整boundary状态图，非原endpoint边缓存、成熟donor表或历史arena。每步对同一完整gene语法的每个状态取最高到达分数，固定FP64的逐碱基role-minus-background累积；null0、完整ATG／单exon终止STOP、实际coding carry、GT或GC–AG、min-intron5及cap128不变。它是max-semiring解码，不用于替代logZ或训练loss；原训练观察及首训CONFIG均未改。

仅保存每8192bp边界的FP64状态向量，回溯时重算单块uint32前驱edge token；不保留整条contig的所有历史前驱节点，不创建完整FP64 gain／prefix矩阵。表示成本为(ceil(n/C)+1)×S×8的边界分数及最多C×S×4的单块前驱，加原输入与选中角色路径；这些是具体数组bytes而非完整进程RSS保证。8192接近30Mb规模两项内存的平衡点，只控制计算表示、不截断gene。全局仍BG初始化及闭合，内部不强制完整起止。新tie规则为先background continuation、再state-table到达顺序，以保证全zero返回null；不宣称与旧endpoint prefix减法版FP64舍入及其他同分链身份逐位相同。

新增四项人工数学验证：独立穷举完整gene集合的最优分数与选中chain／直接角色分数（整数FP32和非整数FP64、actual carry、GC、intronic N及cap）；1／2／3／5／7／8192 checkpoint span的同答案／分数及实际内存公式；zero-null、行公共平移与相邻完整gene；人工128-exon完整chain及cap127排除。同时只因后端变更回归原四个raw接口例，不重跑原训练或配分检查。两项nopython主kernel必须实际执行；分数oracle绝对1e−9级检查，不通过扩大容差救援。

唯一M28-FULL-GRAMMAR-VITERBI-CHECKS.sbatch为generanno/private/0GPU2CPU8GiB3分钟/no-requeue，根outputs/M28-FULL-GRAMMAR-VITERBI-CHECKS与M28FULLVITCHECK_%j日志。相比上轮27秒接口检查，余量覆盖新JIT、128-cap状态图与单块约600MiB回溯，不把这项人工上限当染色体时长。新增源AST／bash通过，五项fresh调度／预约及新文件／根／日志无碰撞、账户仍7个展开作业。纯编译／接线实现错误最多同科学输入／资源／容差另根R1一次；最优分数／完整路径宇宙矛盾、TIMEOUT／OOM或第二次失败停止，不回落到已关闭endpoint分支。全部输入人工、无真实DNA／checkpoint／released GLM或optimizer；真实大长度资源仍待这项新表示证明之后测量。

## 新最大分数自动机实际检查与高密度资源准备

七份相关新源／接线／测试／SBATCH／报告逐文件cmp相同后，单次提交13429852。实际17:34:36–17:35:02 UTC、26秒，root/batch/extern全部COMPLETED 0:0，batch MaxRSS1548348K。四项新数学及四项后端变更后的raw接线共八项全部通过，failures0／errors0；after-import11.592855秒、进程peak RSS1448020KiB，forward与trace各两个实际nopython签名。四个预期非空文件JOBID9、STATUS10、started547、summary1192 bytes共1758 bytes，started／summary身份及配置相同，stdout最后summary逐值相同，八项均ok，stderr为空。没有A43 R2、研究输入、实际checkpoint、native GLM或optimizer更新。这项结果替代当前raw入口的停止分支依赖，不改写13429790真实通过和依赖遗漏的历史。

下一项仅测新checkpoint自动机的高motif密度资源表示：沿已完成独立Python A42的固定49152bp人工IID ACGT（DNA seed0），三个24576bp窗起点0／12288／24576、原nearest-center归属，分别全zero与各窗独立FP32整数[-3,3]（seed41），两项cap128。旧完成13422815_0／1的实际独立最优增益为0／21526；只读取它们的已完成summary及身份，不重跑旧Python或关闭的A43 endpoint。新表示必须取得完全相同的整数最优分数、直接合法chain／实际carry／角色增益／非重叠及null／逐碱基上界，同时记录实际checkpoint和单块前驱bytes。新同分规则可能选不同chain，不要求旧链身份或batch专属计数相同，不以比较fixture计数来代替数学判定。

独立M28-FULL-GRAMMAR-VITERBI-SCALE.sbatch、根outputs/M28-FULL-GRAMMAR-VITERBI-SCALE/case_0及case_1，唯一M28FULLVITSCALE_%A_%a日志；generanno/private/0GPU、每项2CPU8GiB3分钟、array0–1%2/no-requeue。三分钟为新JIT与高密度状态双遍的有界首测余量，不依据旧endpoint耗时保证新kernel；小型summary／chains／records总预计小于1MiB。17:43 UTC五项活调度及预约已核实，账户7个展开作业、两个新增实现及新根／日志均无碰撞，旧A44依赖作业保持原样。新实现AST／SBATCH bash及逐文件同步后单次提交；纯实现错误最多同人工输入／科学设置／资源／容差另根R1一次，最优分数或合法集合矛盾、TIMEOUT／OOM或第二次失败停止。这是新表示的必要资源／独立答案验证，不是旧关闭分支回放，也不是研究模型性能或完整染色体可行性证据。

## 新表示的高密度独立最优分数与资源实际结果

新runner／SBATCH及两份既有报告、两份复用人工fixture helper与canonical逐文件cmp相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-VITERBI-SCALE.sbatch返回13429908。首次live检查已终态，两个array task实际17:46:01–17:46:15 UTC各14秒，root/batch/extern全部COMPLETED 0:0、Restarts0，gpu034/private、每项0GPU2CPU8GiB3分钟，独立命令／WorkDir／%A_%a日志身份相同。batch MaxRSS分别1099156K／1070616K；进程采样峰值1182480／1168792KiB，与Slurm采样口径分别记录、不互相替代。

两项新解码耗时含cold JIT为6.850750／6.864056秒，输入准备0.147142／0.138677秒。全zero返回null及gain0；整数势选中88个完整非重叠gene，gain21526、逐碱基upper142551。实际motif counts为ATG801、donor6102、acceptor2964、STOP2251，contexts归属为[0,18432)／[18432,30720)／[30720,49152)。两个最优整数增益均与已完成独立Python 13422815_0／1完全相同，选中链的实际拼接ORF、contiguous terminal codons、splice、coding carry、cap、非重叠、直接角色增益及总和全部复核；不要求同分链或旧batch计数身份，未重跑旧Python或A43停止原型。

实际新自动机S19078，forward／trace各一个nopython签名；每项checkpoint7行共1068368 bytes、单块前驱625147904 bytes、选中角色路径49152 bytes、输入FP32 logits2949120 bytes，forward加回溯重算98304 bases。回溯可跨8192边界，块边界不是基因关闭或截断点。具体数组与(ceil(n/C)+1)S8／CS4公式相同，但没有把单块数组成本当作完整RSS，也没有把这一长度的RSS视作30Mb上限。

每项五个预期非空产物及started／summary／records／stdout身份全部核对一致，stderr为空；case0为JOBID11、STATUS10、started1021、records825、summary2643 bytes，共4510；case1为11／10／1071／18934／64944 bytes，共84970。旧原24硬观察失败、六个view超界warning、FAILED A43与A44 TIMEOUT及未执行链均完整保留，FIRST-FIT根仍不存在。当前raw入口不再依赖关闭后端，新表示的八项人工数学／接线和两项高密度49kb独立答案资源验证均实质完成。

尚未完成整染色体长度运行、真实训练checkpoint严格加载／native推理或注释性能比较。下一步直接检验这同一固定max表示的染色体长度资源，不再重复已通过的短例／49kb；新首训科学协议确认仍是独立未决事项，不通过人工解码PASS推定许可。目标保持ACTIVE／NOT ACHIEVED；本轮新增任务均终态且已消费，固定heartbeat继续PAUSED，没有仅等待外部作业或重复监控循环。按既有报告续写并保留历史，不生成新的发布或云分享产物。

## 染色体长度人工资源测量准备

本轮live sacct和summary再次确认13429852及13429908两项均COMPLETED0:0，八项检查及两项独立整数最优增益保持通过；不重跑这些短例。读取已完成13429660的原TRAIN长度元数据，六条染色体最长是拟南芥NC_003070.9的30427671bp。因此新资源输入固定为同长度的人工IID ACGT（DNA seed0），不是实际研究DNA；沿原manifest.starts完整24576bp／12288bp网格含移位末窗，各窗FP32整数[-3,3]（seed41）及原nearest-center归属。只有一个高密度非零势实例，不增加seed、科学数据或重复全zero对照。

完整新max automaton源码、cap128、FP64逐碱基累积、null0、同分策略及checkpoint8192均不改。在整个30427671bp定向contig上做一次forward和分块回溯，不截断gene或按窗关闭；记录跨checkpoint／owned-context边界的选中gene、实际两遍bases、数组bytes、RSS和耗时。独立检查每条选中完整链的DNA ORF／carry、连续terminal codons、splice、cap、直接整数角色增益及全局非重叠／总和。按不相交gene span调用原独立verifier，不为每个gene分配全染色体角色向量；逐碱基上界也按65536bp小块核算，避免诊断自身构造全长FP64副本。当前没有此长度的独立穷举最优答案，不能把合法链／分数复核冒称30Mb独立oracle验证。

唯一M28-FULL-GRAMMAR-VITERBI-CHROMOSOME.sbatch、根outputs/M28-FULL-GRAMMAR-VITERBI-CHROMOSOME、M28FULLVITCHR_%j日志，generanno/private/0GPU2CPU8GiB20分钟/no-requeue。49152bp的6.86秒包含首次JIT，不能作为稳定每bp速度；20分钟是同表示高密度长度扩展的有界首测，不承诺真实时长或失败后加时。表示公式预计3716行checkpoint共567150784 bytes、单块前驱625147904 bytes、FP32输入1825660260 bytes，须由作业实算核对；不把这些数组数值当完整RSS保证。新增输出上限64MiB，只写started／selected_chains／summary及状态，不持久保存人工DNA、角色势或模型。17:53 UTC五项活调度／预约核对、账户7作业及新源／根／日志无碰撞，BeeGFS可用154TiB。

新AST／SBATCH bash和逐文件同步后单次提交；仅明确纯编译／接线错误允许同输入／科学设置／资源／容差另根R1一次，合法集合／分数矛盾、TIMEOUT／OOM或再次失败停止，不降低cap、改变分块或重开A43／A44。新首训仍未启动，原训练协议科学确认不由这项资源测量代替；六个原view extent warnings及原硬观察24失败保持不变。这项完成只支持染色体长度的人工势资源可执行，不等同真实raw-FASTA神经推理、性能／泛化或论文目标达成。

## 染色体资源实际提交与固定等待交接

三份新入口／SBATCH／报告及未改max backend与canonical逐文件cmp相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-VITERBI-CHROMOSOME.sbatch返回13429995。实际18:01:17 UTC提交、18:01:18开始，gpu034/private、0GPU2CPU8GiB20分钟/no-requeue、Restarts0，实际Command／SubmitLine／WorkDir／M28FULLVITCHR_13429995日志身份相同。首次live检查为RUNNING28秒；人工原网格2476窗及全长FP32输入1825660260 bytes准备完成，18:01:37.808496进入完整解码。preparation12.500520秒与人工逐碱基upper88050390只是准备事实，不是完成trace或模型效果；尚无最终summary符合运行中协议。

当前仅外部资源作业等待，按用户授权直接复用固定监控01a0f71f-8845-7432-bd13-bf54c1808fb3/local和同heartbeat coding-rna-r5-light-monitor，不新建/fork/cron、改模型或维护主循环。工具更新后实际保存kind heartbeat、target_thread_id原监控、15分钟、ACTIVE、failed_runs_only及原created_at1790981930511，新2021字符prompt仅当前13429995结构白名单与一次交接；旧A44工程故障历史完整保留，不继承其Handoff sent YES。

实际手动初始化turn01a1178d-c8a0-7132-bf58-775816d764b2已completed/errornull/35639ms/cursor155，18:08:50.438Z的真实turn_context核实原gpt-5.6-luna/low/cwd。监控实际18:09:13Z查询squeue／sacct确认13429995 RUNNING/gpu034/0:0，无输出／日志／chain读取或主交接。配置保存及一次手动检查已成立；本cohort定时追加、完整终态产物及自动回传尚未验证，不能拿旧批次定时观察代替。完成或可行动故障时监控一次发送原主并暂停同heartbeat，原主解释资源结果；短暂SSH故障不等于实验失败，监控不批准首训、改goal、取消／重提／加时。

读取原主最近六turn后核实：首次拟合方案已经明确询问新观察模型首训（固定epsilon0.1、原1536 draws／1491 unique／importance、三遍4608步、原已知观察分母、不用DEV）。随后raw接线、新max表示／49kb资源及本轮染色体资源三个连续goal回合没有人类协议答复，自动goal续行不等于此项新科学批准。此前B2启动／费用授权范围不同，不据准备记录、签名或人工PASS推定许可。相同未决边界已达到三回合阻塞审计阈值；主已完成可独立推进的当前准备，下一次真实神经拟合需要原已提出方案的一次确认。这一阻塞不是因正常RUNNING、不把论文目标缩小为工程测试，也不取消或暂停已授权13429995及其固定监控。

## 染色体长度资源测量实际超时

作业13429995于北京时间2026-10-08 02:01:18开始、02:21:23结束（Slurm UTC 2026-10-07 18:01:18–18:21:23），root为TIMEOUT／ExitCode0:0／20分05秒；batch为CANCELLED／0:15／20分06秒，extern为COMPLETED／0:0／20分06秒。stderr明确DUE TO TIME LIMIT；root的0:0不覆盖TIMEOUT状态，不是正常完成。batch MaxRSS5586704K，约5.3279GiB，TotalCPU20:04.209；没有Slurm OOM证据，不将峰值内存或退出码包装为可行性PASS。

stdout仅started和chromosome_decode_started两行：2476个原网格窗口及1825660260-byte人工FP32势已准备，进入decode_full_grammar_viterbi后未输出chromosome_decode_finished。selected_chains.jsonl和summary.json均不存在，故独立逐链核验尚未开始；实际失败阶段在完整解码调用内，现有日志不能细分前向、块回溯或返回前重建。三个遗留产物为JOBID9、started1139、STATUS8 bytes，共1156；stdout1235、stderr115 bytes完整保留。STATUS仍为RUNNING，是终止后未刷新的旧标记，不能据此把已被Slurm终止的作业继续当作live；不覆盖旧标记或伪造完成summary。

源码给出一个具体的复杂度解释：_step在每个base遍历完整状态表，cap128时S19078，checkpoint只降低历史前驱存储，并不免除前向和重算的状态循环；完整两遍的状态遍历规模为2×30427671×19078约1.161×10^12。这是执行结构的事实，不是本次已完成的实际迭代计数或独立profiling结论。kernel为parallel=False，不能凭增配CPU承诺同比加速。20分钟测量只证明这一固定资源下未完成，不证明解码数学错误、GLM训练失败或任何科学性能负结论；此前八项人工验证与两项49kb独立最优分数结果仍保留，但不能外推为染色体资源成立。

按本次明确TIMEOUT停止条件，这一20分钟资源测量已停止：不加时、重提、降低cap、改分块救援或恢复A43／A44。保留源码、failed根及全部旧证据；当前真实染色体执行仍缺可验证结果，后续需要另行确定合乎原完整gene目标的执行方案，不用逐窗关闭或缩小成功定义回避。新观测模型FIRST-FIT根仍不存在，首训科学确认仍未收到，goal保持BLOCKED／NOT ACHIEVED，原论文目标未改；本次工程交接不授权训练或恢复goal。

本cohort固定监控的定时运行与故障回传已实际验证：同一原监控session在北京时间02:23:19（UTC18:23:19.906）的heartbeat输入明确automation_id coding-rna-r5-light-monitor及当前13429995；实际turn_context18:23:19.958仍为gpt-5.6-luna／low／原cwd。它从Slurm和必要stderr取得TIMEOUT后只发送一次原主交接，主已收到并独立消费；保存heartbeat实际PAUSED、固定target／15分钟／failed_runs_only不变。此前“定时／终态回传未验证”只属初始化时点历史；此次已证实定时追加和工程故障回传，不是完整资源测量成功。没有剩余本cohort live作业或后台重复轮询。

## 不同的精确 max 状态支配表示

北京时间2026-10-08的续行从已消费13429995 TIMEOUT出发，不加时、重提或修改原dense kernel。新的full_grammar_skyline.py仍使用同一complete boundary状态表，但把每碱基完整S遍历改为有序可达frontier与外显子预算skyline；不是beam、阈值剪枝、逐窗闭合、降低cap或重开A43 endpoint。原raw caller保持既有backend，训练logZ／coarse observation／FIRST-FIT入口不改。

精确支配关系是：同一DNA位置，状态具有相同stage及actual coding carry，较少已用外显子的分数不低于较多者时，后者可以删除。对高count的每一条未来路径，低count可逐边产生相同角色、保持相同stage／carry及count差，cap剩余额度更大；完成STOP回到BG后差消失，因此保持全局最优值。FP64每碱基同一增量的舍入加法是单调的，这不是把prefix差换成其他累积口径。不同carry即使phase相同也不合并，较高count但更好分数仍保留。等分时低count支配可能改变旧dense同分链身份；新tie明确background continuation、state-table到达及等分低count支配，需在不同checkpoint span下保持自身答案。此关系只用于max，不能用于sum/logZ、训练梯度或gene集合概率。

实现以两个dense分数数组作稀疏累加器，仅清除上轮活状态；候选索引排序后按stage／count／carry顺序删除被支配者。checkpoint仍存完整S向量，单块uint32前驱数组只写实际到达的状态，回溯只读选中活路径，校验边及边界分数与最终逐碱基总分。结构proof在cap1／2／4／128的真实图上核对相邻count每个ACGTN转移及角色包含关系，传递覆盖全部count；若失败不声称精确剪枝。新增实际活状态访问数、edge proposals、dominated states和最大frontier计数用来判断资源收益，不用旧dense潜在循环规模冒充实测计数。当前仍没有速度或染色体可行性PASS。

六项新检查只针对新表示：实际状态图的role-preserving simulation；等分／较差count删除与较好count及不同actual carry保留；五个短人工序列、caps1／2／4、整数FP32和非整数FP64的独立完整gene集合oracle及合法chain；checkpoint1／2／3／5／7／8192自身同答案与数组bytes；全zero null、行平移、相邻gene及实际稀疏访问；人工128-exon最高分链完整保留且cap127排除。无旧solver重跑、真实DNA、released GLM、训练checkpoint、optimizer或DEV。新入口M28-FULL-GRAMMAR-SKYLINE-CHECKS.sbatch固定generanno/private/0GPU2CPU8GiB3分钟/no-requeue，新根outputs/M28-FULL-GRAMMAR-SKYLINE-CHECKS及M28FULLSKYCHECK_%j日志。纯编译／接线错误最多同输入／资源／科学设置／容差另根R1一次；数学／路径宇宙矛盾、TIMEOUT／OOM或第二次失败STOP。仅实际六项通过后才考虑新表示的49kb独立最优分数与资源测量，不提交新的整染色体作业来掩盖13429995负结果。

原20分钟dense测量停止、原24硬观察FAIL、六个primary extent warnings、A43／A44失败与未执行证据均保留。新首训科学确认仍未收到，自动goal续行只使当前goal ACTIVE／NOT ACHIEVED，不构成该首训批准；本轮在不依赖它的已授权精确表示方向推进。固定监控仍PAUSED，没有创建重复对话、心跳或后台工作。

## 精确 skyline 实际数学结果与独立资源首测

五份新源／测试／runner／SBATCH／报告与canonical逐字一致；新根及日志无碰撞，五项fresh调度／预约检查无当前private限制，账户五项展开作业低于80。实际M28源与连续已完成作业使用generanno，沿用该环境而不采用cluster_config.yaml历史coding-rna默认。单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-SKYLINE-CHECKS.sbatch返回13434294。北京时间09:19:43–09:20:00（UTC01:19:43–01:20:00）实际17秒，root／batch／extern全COMPLETED0:0、Restarts0，gpu034/private、0GPU2CPU8GiB3分钟及命令／工作目录／日志身份相同。

六项新测试全通过、errors0／failures0；两个forward与两个trace nopython签名实际执行，after-import10.78277055秒、进程peak RSS620064KiB，Slurm batch MaxRSS713104K分别记录。真实完整图ACGTN逐角色simulation成立，独立短例最优增益和合法chain成立，128-exon最高分链未被错误删除；这是新精确剪枝的数学／人工实现证据，不是染色体资源或模型效果。四个非空文件JOBID9、started570、summary870、STATUS10 bytes共1459；配置／job identity与stdout最后summary逐值相同，stderr0 bytes，已由主消费。没有运行旧dense、endpoint或任何研究输入／optimizer。

下一项仅针对新skyline表示测固定49152bp高motif密度人工输入：原DNA seed0、24576bp／12288bp三个窗口、nearest-center owner、case0全zero及case1 FP32整数[-3,3] seed41、cap128／checkpoint8192／逐碱基FP64／null0不变。只读已完成独立Python13422815_0／1的最优分数0／21526及固定输入身份；不重跑旧solver。每项需直接核实完整chain DNA／ORF／motif／carry／cap／非重叠／整数角色分数与独立最优值，并报告实际frontier访问、dominated states、最大活状态、cold-JIT总时长和RSS。潜在dense N×S访问数只作执行结构参照，不能称旧算法的已测计数或把比值当实际时间加速。

新M28-FULL-GRAMMAR-SKYLINE-SCALE.sbatch、唯一根outputs/M28-FULL-GRAMMAR-SKYLINE-SCALE/case_0及case_1、M28FULLSKYSCALE_%A_%a日志；generanno/private/0GPU，每项2CPU8GiB3分钟、array0–1%2/no-requeue，两项小型产物预计小于1MiB。新runner AST／SBATCH bash与fresh五项调度及碰撞核对后提交；纯编译／接线错误最多同输入／资源／科学设置／容差另根R1一次，数学／合法集合矛盾、TIMEOUT／OOM或第二次失败STOP。未改变raw caller后端或首训协议；这项49kb资源PASS也不能解除旧dense20分钟TIMEOUT停止条件或证明真实染色体／训练／泛化。

## skyline 高密度实际资源与完整长度下一步

三份新入口／SBATCH／报告与canonical逐文件cmp相同后单次提交13434301。两array实际北京时间09:22:07同时开始，09:22:21／09:22:23结束（UTC01:22:07–01:22:21／23），14／16秒，root／batch／extern全部COMPLETED0:0、Restarts0；gpu034/private，每项0GPU2CPU8GiB3分钟，Command／WorkDir／原%A_%a日志身份成立。case0／1进程peak RSS1192340／1206264KiB，Slurm batch MaxRSS1116392／1086516K，采样口径分列。

最优整数gain0／21526与已完成独立Python精确相同，全zero为null，非零势88个选中完整gene、upper142551。原三个context segments、实际DNA ORF／terminal／splice／carry／cap／非重叠／逐角色直接整数分数及sum全部核实；同分chain身份不作旧版比较。仍S19078／residual groups154，checkpoint7行1068368 bytes、单块前驱625147904 bytes、FP32输入2949120 bytes与两遍98304 bases不变。case0实际forward活状态访问1758027／edge proposals1894062／dominated19151／max active56；case1为10153962／10964392／133011／1040。相对原结构N×S=937721856的访问规模，约减少533／92倍，但排序和其他开销未包含在这个比值里，不能说实际时间加速533／92倍，也不能保证长序列frontier上限1040。

两项decode含cold JIT为7.567559／7.948051秒；前一dense49kb实测6.850750／6.864056秒仍完整保留，新小长度总时长并未更快，不声称已解决整染色体超时。两项数学上的精确剪枝及实际活状态规模减少，给出与旧dense不同的执行机制；是否转化为完整长度可执行须直接测量。每项forward／trace各一nopython签名。五个非空产物case0 JOBID11／started1075／summary2996／STATUS10／records1049 bytes共5141；case1为11／1125／65301／10／19162共85609。started／summary／records与stdout结构逐值相同，stdout2095／2169 bytes、stderr0，已消费。

新profile_full_grammar_skyline_chromosome.py只替换实际后端为已通过graph／oracle的skyline、对应新prerequisite和活状态计数；与原driver的diff确认30427671bp人工IID ACGT seed0、原2476窗FP32整数seed41、owner、cap128／checkpoint8192／单位角色／null0／FP64口径和独立有界逐链核验完全不变。依旧只读原13429660六TRAIN长度元数据，不读研究FASTA／labels。新唯一root outputs/M28-FULL-GRAMMAR-SKYLINE-CHROMOSOME、M28FULLSKYCHR_%j日志、generanno/private/0GPU2CPU8GiB20分钟/no-requeue、64MiB小型产物上限。20分钟是这个不同执行表示的有界完整长度首测，不是给旧13429995加时或重提；不依据cold49kb时长作确定性线性预测。不重复新短例、49kb或旧solver以等待结果。

新driver AST／SBATCH bash、fresh五项调度／预约、唯一文件／root／log及storage通过后单次提交。仅纯接线／编译错误可同人工输入／资源／科学设置／容差另根R1一次；数学或合法集合矛盾、TIMEOUT／OOM或第二次失败STOP，不降低cap、改checkpoint、增加walltime救援或恢复旧dense／A43／A44。无完整染色体独立穷举oracle，需分清最优路径数学证明、选中路径复核与实际资源是否完成；无模型性能／泛化证据。raw caller尚未切换、FIRST-FIT仍未启动，新科学确认不由这一工程首测替代。只有外部等待时复用原固定轻量监控，不创建新对话／心跳或主轮询循环；goal不标完成。

## skyline 完整长度实际提交与固定交接

三份新driver／SBATCH／报告与canonical逐文件cmp相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-SKYLINE-CHROMOSOME.sbatch返回13434312。实际北京时间2026-10-08 09:25:45提交、09:25:46开始（UTC01:25:45／46）；首次squeue／scontrol／sacct为RUNNING10秒，gpu034/private、0GPU2CPU8GiB20分钟/no-requeue、Restarts0，Command／SubmitLine／WorkDir及M28FULLSKYCHR_13434312日志身份相同。BeeGFS可用154TiB，账户原五项作业加一低于80，原输出根及log无碰撞。

stdout started身份／配置符合新skyline，旧dense13429995明确仅是stopped source而非重跑；2476窗人工输入准备13.013713秒、1825660260-byte FP32 logits及upper88050390与原人工输入口径一致。09:26:05.629714（UTC01:26:05.629714Z）进入完整解码，初始化时尚无decode_finished／selected链／summary，不能称完成、最优染色体结果或资源PASS。原raw caller未切换，首训未启动，所有历史负结果及未执行事实不改。

只复用原固定监控01a0f71f-8845-7432-bd13-bf54c1808fb3/local与heartbeat coding-rna-r5-light-monitor。更新前read_thread核对该对话idle及已消费旧13429995 TIMEOUT；journal仅在顶部加入新13434312 cohort，Handoff sent NO、connection failures0，全部旧history保留。automation_update实际ACTIVE，保存字段精确核对kind heartbeat／原target／15分钟／failed_runs_only／created_at1790981930511及新2166字符prompt，未改模型、创建或fork对话、cron或第二个主轮询。当前prompt只查最小Slurm状态，COMPLETED后结构白名单核验，正常RUNNING静默，仅终态／可行动故障向原主一次接续后暂停。

同一监控初始化消息已送达，未传model／thinking override；配置成功与初始化实际运行、定时运行、终态回传分开记录。原监控历史实际gpt-5.6-luna／low已核实，不把旧定时证明替代本cohort。当前仅外部资源等待，由固定监控负责，主不维持循环，不通过工程等待解除新首训科学确认或把原论文目标标完成。

本cohort实际初始化turn01a1191f-6e50-74d1-8ac7-61ca6170218e已追加到原监控；真实turn_context时间UTC2026-10-08T01:27:33.402Z，cwd／model gpt-5.6-luna／effort low与原设置相同。journal记录实际09:29:38（UTC01:29:38Z）一次squeue／sacct查询，job13434312 RUNNING/gpu034/0:0，未读产物、日志或科学内容，Handoff sent仍NO、无主接续。保存配置与单次检查成功已成立，本cohort定时追加和最终自动回传尚未验证；不以历史13429995的定时／故障交接替代。主不再查询运行中job，后续只消费固定监控一次终态交接并按已授权范围继续验证。

## skyline 的原 raw 接口前瞻接线

上一续行的六项结构／数学检查及两项49kb独立答案与新完整长度提交是实际进展；本续行不重跑这些用例，也不因goal继续推定新观测模型首训许可。原固定监控journal仍为当前13434312 Handoff sent NO，主不启动第二个资源等待循环。为提交新的短接口检查而做的必要账户调度检查，实际09:34:04观察13434312仍RUNNING8分18秒；这是提交决策的live状态，不是额外轮询或资源结果。

check_full_grammar_skyline_raw.py仅在测试进程临时把原raw模块的decoder绑定到实际decode_full_grammar_skyline，不改正式raw源、训练配置、checkpoint身份或公开接口。复用原两个会执行decoder的人工case：36864bp上下文中30000bp完整gene跨24576bp原生窗，在正负链恢复坐标及GFF3；14bp actual coding carry／intronic N的负链双exon phase2／0。预期共六次实际新skyline调用，每次cap128及精确支配元数据须吻合。原Head／native encoder为stub，十五角色FP32势不变；执行真实whole-contig decode、actual DNA检查和GFF3写入／读回，不执行旧dense／endpoint kernel。未变化的独立owner和最终checkpoint身份两项不重跑。通过只证明前瞻接线，不能称正式后端已部署、真实checkpoint推理或完整长度资源PASS。

新增唯一M28-FULL-GRAMMAR-SKYLINE-RAW-CHECKS.sbatch及root outputs/M28-FULL-GRAMMAR-SKYLINE-RAW-CHECKS、M28FULLSKYRAW_%j日志；generanno/private/0GPU2CPU4GiB2分钟/no-requeue，沿用原人工raw检查资源而不以低motif长gene耗时预测染色体。新AST／bash、当前remote配置／worktree和相关源相同、fresh五项调度／预约及文件／root／log无碰撞后单次提交；账户六项展开作业加一低于80，小型JSON产物预计小于64KiB。仅纯编译／接线错误可同输入／资源／科学设置／容差另根R1一次；实际完整chain／score／DNA／GFF phase矛盾、TIMEOUT／OOM或第二次失败STOP。没有新增监控对话或更新现有完整长度cohort，不恢复A43／A44／旧dense测量、不读研究FASTA／label／checkpoint／DEV／独立物种／A33、无optimizer更新。原24硬观察失败和六个primary extent warnings仍保留。

## skyline 的原 raw 接口实际检查

三份新增runner／SBATCH／报告与canonical逐文件cmp相同后，单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-SKYLINE-RAW-CHECKS.sbatch返回13434348。实际北京时间2026-10-08 09:35:53–09:36:14（UTC01:35:53–01:36:14）21秒；root／batch／extern全部COMPLETED0:0、Restarts0，gpu034/private、0GPU2CPU4GiB2分钟，Command／SubmitLine／WorkDir及独立M28FULLSKYRAW_13434348日志身份成立。两项tests实际运行、errors0／failures0；after-import7.877587秒、进程peak RSS1450348KiB、Slurm batch MaxRSS1548304K。

六次调用均进入新skyline kernel，forward／trace各一个实际nopython签名。正负链的30000bp人工完整gene均保留原定向[3000,33000)／gain300000，经原caller输出genomic CDS [3000,33000)／[3864,33864) phase0及正确GFF3 1-based坐标；两次相反链zero势均null。14bp carry／intronic N例目标负链恢复定向((0,4),(9,14))／gain140，经原caller输出CDS [0,5) phase2／[10,14) phase0，另一链null。真实DNA／phase检查和GFF3导出断言通过，没有把跨原生窗完整gene拆开再owner过滤。36864bp调用的实际checkpoint6行915744 bytes、single-block predecessor625147904 bytes；长gene最大frontier5、carry例7，各例dominated state均0，因其低motif结构不制造外显子预算竞争。这项通过不能替代之前高密度dominated-state／独立最优值检查，也不能外推30Mb耗时。

四个非空产物JOBID9、started707、summary6591、STATUS10 bytes共7317；started／summary配置与job identity同一、stdout末尾JSON与summary逐值相同、stdout5509／stderr0 bytes。测试临时binding结束后正式源仍未改，summary production_raw_backend_changed=false；本次实际证明新decoder可以通过原raw模块完成两项接口用例，不声称正式入口已切换。未运行native GLM、读取真实训练checkpoint或研究DNA／标签，无optimizer更新，不重新运行已通过的owner／checkpoint身份项。

本回合已取得新的实际接口证据，原科研目标仍未完成；当前剩余完整长度资源由固定监控13434312 cohort接管，不加时、重提、反复主轮询或新增监控。若它按原边界完成，再合并这一接线结果处理实际raw默认依赖；若TIMEOUT／OOM／数学或再次失败，遵守新skyline资源停止条件，不能凭短例通过放行。首训、真实checkpoint推理、准确率／适用范围比较及论文级证据尚未完成，原新观测模型科学确认仍未收到。历史24硬观察FAIL、六个extent warnings、A43／A44／旧dense失败与未执行链不改。

## 首训未决边界与资源等待分离

上一回合实际完成13434348两项新decoder的raw／GFF3接线检查，属于进展，不是仅重复状态。本回合读取canonical首训合同及最近八回合，原明确问题仍是完整语法新观测模型首训：固定epsilon0.1、原1536 draws／1491 unique及原importance、三遍4608更新、原已知观察分母、不用DEV；没有新的直接人类答复，FIRST-FIT根仍不存在。B2启动及不考虑费用的既有批准不自动覆盖新的观测概率／联合NLL机制；自动goal续行也不是该项科学授权。准备元数据、native四窗接入、完整状态图／max数学、49kb独立答案、raw前瞻接线均已完成，不能继续用同一通过用例或新准备文件替代真实拟合。

此前BLOCKED恢复后，新skyline图／数学／49kb／完整长度提交、前瞻raw接线及本回合构成三个连续续行；前两回合完成了不同的已授权准备，但未解除同一首训未决边界。当前无需该答复的主准备已完成，下一次真实训练需要此前已提出科学方案的一次直接确认。满足首训阻塞审计阈值，平台update_goal实际返回主goal状态blocked／NOT ACHIEVED；原因是强制首训授权尚未补齐，不是计算慢、资源费用、普通RUNNING或需要更多审阅。原Nature Communications及以上论文、严格raw-FASTA完整gene注释与ANNEVO／Tiberius／Helixer性能或适用范围比较目标不变。

当前同一heartbeat实际仍ACTIVE／kind heartbeat／原target／15分钟／failed_runs_only，13434312 Handoff sent NO；主没有收到其终态交接，不据运行标记或历史snapshot判成功／失败。只保留原固定监控，既有作业与停止条件不改、不取消或重提、不创建新等待循环。终态结果仍按原授权一次交接和验证；工程交接不会放行新首训、恢复goal或改科学主张。本cohort定时／完整回传尚未由主验证，机器／应用／连接条件仍适用。

## 2026-10-08 人类批准首训与完整长度实测消费

用户在上述明确首训问题后直接回答“我同意你的方案”。这次答复批准新观测模型固定epsilon0.1、原1536 draws／1491 unique／importance／标签、三遍4608步、已知观察分母、无DEV及private单3090／2CPU16GiB／30小时／no-requeue方案。原首训批准缺失已解除，按相同合同直接执行，不再以旧准备文件或平台goal状态要求重复批准；不改变原raw-FASTA／比较基线／论文目标、封存数据或公开发布边界。

主实际查询sacct，13434312 root／batch／extern全部COMPLETED0:0，2026-10-08北京时间09:25:46–09:35:34（UTC01:25:46–01:35:34），9分48秒；batch MaxRSS5583388K，进程peak RSS5518316KiB，约5.3GiB，采样口径分列。仍是30427671bp人工IID DNA seed0／原2476窗FP32整数势seed41／cap128／checkpoint8192／单位role／null0，未读研究DNA、标签、模型或checkpoint。preparation13.013713秒、decode含cold JIT537.100558秒、独立逐选中链复核30.284062秒；gain13259529，pointwise upper88050390，57968条完整gene、最大122外显子／5589bp span，3497条跨checkpoint、2298条跨owner上下文。所有实际ORF／carry／contiguous terminal STOP／splice／cap／非重叠及整数分数由runner重算；没有完整染色体独立穷举最优oracle。

实际forward active visits6459799335、edge proposals6990537444、dominated80110767、maximum active1624；checkpoint3716行567150784 bytes、single-block predecessor625147904 bytes、input1825660260 bytes、两遍60855342 bases。五个非空产物JOBID9／started1261／summary3362／STATUS10／selected_chains16191378 bytes，共16196020，低于64MiB；选中行数57968，started／summary配置与job身份一致、stdout最后JSON与summary逐值相同、stderr0。此结果是新精确执行表示在这一人工完整长度输入上的完成，不给出旧dense的确定倍数加速、不外推真实模型势、不声称准确率或科研PASS；原dense TIMEOUT和A43／A44关闭事实保留。

固定监控的定时回合已实际追加到原对话，但查询用squeue && sacct；作业离开squeue后Invalid job id令sacct被短路，三轮未取得终态，实际没有成功自动回传。主已直接消费终态并暂停同一heartbeat，随后将修正为两个调度查询独立执行，再复用同一对话／同一自动化监控FIRST-FIT；不将配置成功冒称完整定时交接成功，不重复旧cohort手工消息。

## 新观测模型获批 FIRST-FIT 实际启动

合同与本记录同步canonical并逐文件cmp一致，fit源cmp一致／AST通过、SBATCH bash通过；canonical配置虽仍写旧coding-rna默认环境，此M28链沿实际已完成native13429660和准备13429756的generanno，不往base安装。fresh五项Slurm／预约显示private gpu034有3个、gpu035有8个空闲3090，账户五个展开作业加此项低于80，无相关GPU维护预约，BeeGFS可用153TiB；FIRST-FIT根及同名活作业均无冲突，独立%x_%j日志目录可写。没有重跑已通过的native／准备／数学诊断。

单次sbatch --parsable sbatch/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT.sbatch返回13434603。实际北京时间2026-10-08 10:28:35提交及开始（UTC02:28:35）；首次即时查询PENDING0:0／None，随后主启动核验RUNNING1:30、gpu034。root／batch／extern RUNNING0:0；Command／SubmitLine／WorkDir、输出M28COARSEFIRSTFIT_13434603.out/.err、AllocTRES private／1×nvidia_geforce_rtx_3090／2CPU／16GiB、TimeLimit30小时、Requeue0／Restarts0、无dependency均吻合。独立根outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT内JOBID13434603、STATUS RUNNING、started.job_id13434603／preparation_job13429756／configuration与原准备逐值相同；training.jsonl已创建，启动核验时可见行数0，stdout／stderr为空，不推定已完成optimizer更新或数值拟合。日志按64步flush，空日志不等于失败；正式验收仍需4608连续步骤、三个checkpoint及summary。

同一固定监控01a0f71f-8845-7432-bd13-bf54c1808fb3和heartbeat coding-rna-r5-light-monitor转到当前唯一13434603 cohort；旧skyline终态已消费、不补发重复交接。automation实际保存kind heartbeat／ACTIVE／原target／15分钟／failed_runs_only，prompt逐字读回相同、created_at1790981930511保留。原监控实际新turn01a11958-4ec9-7242-a0a2-b47d691b31c8已追加，turn_context UTC02:29:40.238仍原cwd／gpt-5.6-luna／low；主模型未改。prompt修正两状态查询独立执行，终态只按文件／身份／步数／布尔白名单交接，不读loss、checkpoint tensor或解释性能。人工初始化实际UTC02:30:23（北京时间10:30:23）成功，独立执行的squeue及sacct均返回RUNNING／gpu034／0:0，命令exit0；该同对话turn已completed／error null／53787ms、cursor164，journal标当前人工检查YES，无产物／指标读取或进度交接。配置成功及一次人工检查已验证；当前cohort的定时检查和完整终态回传尚未验证，不冒称训练完成。机器／应用／SSH可用条件不变。
