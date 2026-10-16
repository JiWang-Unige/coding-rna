# 原 TRAIN 共同编码上下文适配

2026-10-09。固定4608权重的全部1491视图比较显示，24576共同编码在两物种均降低汇总NLL并提高有限视图完整参考链配对命中，但每物种仍有退化视图；49k／98k收益不单调。本分支检验在相同起点、同新增4608更新下，对完整24576读取窗口的原生特征适配是否改善整染色体CDS链恢复，不将固定权重诊断当作新模型性能。

## 固定对照与数据

唯一改变是每个原24576读取窗口由四个独立6144块改为一次原生共同编码。采用原b048 GENERanno、原bidirectional SDPA flash、六碱基token、BF16隐藏状态后FP32头；不改变RoPE或backbone参数，不利用两个native CDS头，不加入新窗口。24576是原神经读取单位，训练各邻窗口和全染色体推理各窗口的可见输入相同；不直接将49k／98k的中心共享视图条件当作每窗口raw推理，避免混淆其上下文归属。短FASTA窗口在定向后按原规则右补N到24576，只有真实有效碱基进入势函数及DNA解码，不伪造真实染色体侧翼。

从原final4608 C0及全部AdamW状态独立恢复，新增三遍原1536 draws、1491 unique，共4608更新，正式末权重累计9216。原物种／染色体／方向顺序、importance、fine／coarse／unknown、known-row分母、结构coarse-channel epsilon0.1、24576核心及左右各一个12288步长邻窗、nearest-center ownership、cap128、单位角色、null0、FP64分区函数不变。优化器仍AdamW lr3e-4、weight decay0.01、betas0.9／0.999、eps1e-8、clip1，C0仍366191参数、seed0。已有STRUCTURED_CONTINUE是同起点／相同更新次数／同抽样的6144匹配对照，不重复启动新控制。

原两物种六TRAIN限定不变，不读取DEV、sealed、held-out或Setaria。不把本分支自动并入已冻结八方法DEV家族；进入DEV前单独固定比较方案，不能用DEV反馈选择context、checkpoint或损失权重。backbone上游监督来源的既有局限仍有效，原TRAIN改善不证明独立泛化或原创机制。

## 接通验证

人工输入检查全4096个六碱基token、短输入定向后N padding、完整长度及有效行截断；新checkpoint入口接受本分支末9216、拒绝原6144或中间权重。对人工短序列核验新双链raw调用器仍使用完整原DNA和原skyline规则，不读取参考。

原生预检先用固定原4608头在全部八个旧见证重放JOINT24576势函数，要求精确相同；不挑性能好窗口。随后对原第一draw的各物种各独立恢复一次，做可丢弃step4608→4609更新，核对有限NLL、所有C0梯度、冻结backbone、AdamW及参数变化；原生特征在inference mode之外clone。相同更新后参数经新raw adapter严格加载，重算同输入势函数须一致。adapter的末步metadata仅为人工路由fixture，预检不保存科学checkpoint，不进入正式更新。预检单RTX3090／2CPU／16GiB／15分钟，上限以现有原生预检和八见证重放路径为依据，实际通过后才启动fit。

## 正式训练与端到端评价

单独输出M28-FULL-GRAMMAR-JOINT24-FIT，私有单RTX3090／2CPU／16GiB／18小时；沿用现有结构续训的DP和更新上限。共同编码是同一原生窗口的一次4096-token forward，具体速率由预检记录，不以两个样本推算保证全队列完工。每1536次保存checkpoint，仅step_009216正式；日志512KiB缓冲、每64步刷新、fit文件96MiB上限，不保存全范围特征缓存或改写运行中的源代码。

完成后CPU核验全部4608日志的原顺序、importance、观测数、分母、目标算术及三checkpoint的严格加载、权重／AdamW矩有限和逐参数步数，2CPU／4GiB／5分钟。通过后六原TRAIN染色体raw数组最大并行3，每任务单RTX3090／2CPU／16GiB／12小时；按相同stride、owner、整定向染色体联合skyline解码，reference和标签不进入caller。预测输出每染色体128MiB上限，保留固定输入FASTA，无势缓存。

评价CPU整数组afterok，原M26最长CDS primary／ID tie尺和全部可评价参考、unsupported层不变，2CPU／4GiB／20分钟。报告每染色体、两物种、worst species及pooled的链precision／recall／F1、CDS碱基、背景和重复链，与原4608及同更新次数STRUCTURED_CONTINUE分开比较；不是有限视图pair命中，也不是独立测试。不因指标负结果改条件或重训。实现故障保留原产物，只允许一次同科学设置的明确修复。

原固定监控对话及唯一心跳继续复用，已结束准备不反复轮询；新正式fit及整条后续流程只作为一个cohort加入原监控，正常运行静默，完成或需处理时由主分析。

## 输入和调用器检查结果

CPU检查13453753的Slurm root／batch／extern均COMPLETED 0:0，实际UTC12:06:49–12:07:10。四项人工测试全部通过：完整及1／14 bp短输入的4096-token保留、补N与有效行截断，inference-mode之外clone可参与梯度；checkpoint拒绝原6144编码及中间权重；14 bp反向链的跨内含子codon carry及GFF相位；40000 bp末窗不对齐六碱基token时双链nearest-center ownership及角色列保持不变。这里只用了人工输入和stub／mock，未读研究数据或checkpoint，未运行原生编码器，也未做优化器更新。

使用`sbatch sbatch/M28-FULL-GRAMMAR-JOINT24-PREFLIGHT.sbatch 13453753`提交原生预检13453758，单RTX3090／2CPU／16GiB／15分钟、private-teodoro-gpu。正式fit须等待实际八见证重放、原生短FASTA和两次独立可丢弃更新全部通过；不以人工检查代替原生证据。

原生预检13453758的root／batch／extern现已COMPLETED 0:0，UTC12:07:53–12:09:31，98秒、gpu034。全部八个旧见证的JOINT24576势函数精确重放；真实14 bp短FASTA双链caller完成，保留完整有效DNA且无参考输入。原各物种第一draw独立step4608→4609更新通过，NLL和梯度有限、所有C0梯度存在、backbone冻结，恢复后的raw adapter势函数与训练势函数逐元素精确相同。原生更新只作接线验证，均丢弃，未保存科学checkpoint或更改原起点。人工末9216metadata只用于加载路由fixture，不是已完成末9216训练。

## 正式作业依赖

预检实际通过后已提交fit13453777 → CPU核验13453778 → 六染色体raw数组13453779_0–5（最多并行3）→ 整数组afterok评价13453780。四阶段分别使用本方案的18小时GPU／5分钟CPU／每task12小时GPU／20分钟CPU资源，不因初始调度PENDING／None判作失败或重复提交。

实际命令为`sbatch sbatch/M28-FULL-GRAMMAR-JOINT24-FIT.sbatch 13453758`；后续依次为`sbatch --dependency=afterok:13453777 sbatch/M28-FULL-GRAMMAR-JOINT24-VERIFY.sbatch 13453777`、`sbatch --dependency=afterok:13453778 sbatch/M28-FULL-GRAMMAR-JOINT24-TRAIN-RAW.sbatch 13453777`和`sbatch --dependency=afterok:13453779 sbatch/M28-FULL-GRAMMAR-JOINT24-TRAIN-EVAL.sbatch 13453777`。后续指标尚未产生，正式训练／推理完成和科研成功分别判断；与新抽样候选的未来DEV范围另见[单独固定的开发比较](full_grammar_sampling_context_dev_plan.md)。

主实际核对fit13453777从UTC12:11:14在gpu034运行，原生起点和配置正确、stderr为空。CPU核验依赖fit、raw依赖核验，评价的实际Slurm依赖为afterok:13453779_*，涵盖全部数组任务，不是单个父项成功。已经在UTC12:13:05.129更新原唯一心跳，保存ACTIVE／heartbeat／原target_thread_id／15分钟／failed_runs_only／原创建时间及完整四cohort prompt逐项相符，未新建监控或改模型。新批次独立首检与真实定时追加的证据另行记录，不从配置成功推定实际执行。

原监控对话的新批次手动独立首检已于UTC12:14:21实际检查四条cohort，12:14:39完成；turn01a12095-d6e1-7f83-908a-1cd4dab0c825实际模型仍gpt-5.6-luna／low。五个尚未结束fit运行，其余下游为有效Dependency；LOCAL分支局部已消费但三目标整体仍未完。没有读指标或发送终态交接，该次是手动检查，不是定时追加。

新四cohort规则的真实定时追加现已核验：原监控对话turn01a12097-e53d-7ac3-a01f-3374e488fb91，实际UTC12:16:27.975–12:16:48.320，输入包含heartbeat及全部四个新JOINT24作业ID，模型仍gpt-5.6-luna／low。12:16:32实际检查四条cohort，五个fit运行、后续有效Dependency，未读指标或交接。配置保存、手动独立检查、真实定时追加三项分别得到证据；最终终态自动回传仍待流程结束。固定原心跳继续ACTIVE，不新建或替换等待机制。
