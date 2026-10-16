# M27 原生配对结果：严格绕行已观察到，尚未完成机制归因

Takeaway: Both native callers show PTC-specific structural bypass in the fixed M27 panel, without biological attribution.

在固定的整条人类22号染色体上，同一密码子的同义替换与人工早停替换产生了不同的原生注释响应。两工具的全部计划推理均已完成，所有同义条件均保留原有目标链。ANNEVO的12个位点中有2个明确绕行、2个不可归属，Tiberius的10个位点中也是2个明确绕行、2个不可归属。保留未知病例后的配对差值识别下界仍为正，但这不是置信区间或人群发生率。共同9个位点中，两工具的明确绕行不在同一个基因，不能据此声称同位点机制一致。下一步只按事前顺序确认每工具一个案例的原生重现性，不把本结果称为真实剪接、失活或论文新机制。

## 1. 科学问题与范围

结果核对日期2026-09-27。依据 [事前合同](paired_native_contract.md)（首个mutant前公开于c284cd0919884df91e4fa47d85020fbfa5ac5364），本轮测量各工具WT唯一精确集合内的原生配对总响应。固定GRCh38 chr22、GENCODE49、uppercase/no-softmask、单SNV、原版本/权重/精度/自适应流程；0训练，0Setaria访问。原21个参考登记位点全部保留，不因WT或mutant结果补选。

ANNEVO条件分母12，Tiberius10；共同WT精确9只用于同位点比较。每份mutant独立来自同一WT并通过全序列单位置差异与正负链密码子验证。严格BYPASS仍要求唯一单链保留完整双端CDS锚，整个编辑密码子进入一个内含子，剪接连接改变、ORF与phase有效且无约定的外来编码混杂；未放宽规则。

## 2. 执行状态与终点

Slurm job13227849（1×RTX3090/gpu035、8CPU、32GiB、6h上限）为COMPLETED 0:0；全部44份执行receipt均为COMPLETED且归属同一job，68个模型/解码子进程returncode均为0。24次ANNEVO及20次Tiberius完成，0技术性未完成；主错误日志为空。前置17项测试已在CPU准备阶段通过，不把测试当生物学验证。

分配时间3:43:20，即3.722222 GPUh；MaxRSS14,826,788KiB（约14.14GiB）。连同WT的949秒，M27累计3.985833 GPUh，8GPUh总预算尚余4.014167 GPUh。逐例记录之和不代替分配时间。

24份ANNEVO临时HDF5均按事前政策在成功解码/分类后删除，共73,178,701,440 bytes；WT缓存、完整预测和日志仍保留，临时分数可通过相同输入重算。20份Tiberius轨迹全部complete/native，共80次真实神经前向，每例实测4次请求；4是观察值而非事先强制条件，mutant缓存写入0。M27六个声明目录终态合计13,430,714,955 bytes（12.51GiB）；逐例清理前最大记录16,473,234,377 bytes（15.34GiB），不是连续采样得到的全程存储峰值。

## 3. 主要结果

### 固定分母与未知病例

| 工具/集合 | n | 配对00 | 配对01 | 配对10 | 配对11 | 未知配对 | PTC−syn识别界限 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ANNEVO自身WT精确集合 | 12 | 8 | 0 | 2 | 0 | 2 | [0.166667, 0.333333] |
| Tiberius自身WT精确集合 | 10 | 6 | 0 | 2 | 0 | 2 | [0.200000, 0.400000] |
| ANNEVO共同WT集合 | 9 | 6 | 0 | 1 | 0 | 2 | [0.111111, 0.333333] |
| Tiberius共同WT集合 | 9 | 6 | 0 | 1 | 0 | 2 | [0.111111, 0.333333] |

配对编码依次是PTC、syn；10表示PTC明确bypass而syn明确非bypass。界限来自每条件[0,0]/[1,1]/[0,1]的固定分母运算，不是抽样置信区间，不以完整病例子集作分母。不同12/10集合不能直接排名工具优劣。

所有syn均为WT_CHAIN_RETAINED（ANNEVO12/12，Tiberius10/10），但PTC有多种不同响应：

| PTC类别 | ANNEVO n=12 | Tiberius n=10 |
|---|---:|---:|
| BYPASS | 2 | 2 |
| PTC_TERMINATION | 1 | 1 |
| OUTPUT_SPLIT | 1 | 5 |
| OUTPUT_ABSENT | 6 | 0 |
| UNASSIGNABLE | 2 | 2 |

其余预定类别均为0。OUTPUT_ABSENT是该完整运行中无目标相关输出，不是已证明生物学基因消失；UNASSIGNABLE不改记阴性或消失。

### 全部21行登记账

| ordinal | gene | 编辑坐标 chr22，1-based | ANNEVO PTC | Tiberius PTC |
|---:|---|---:|---|---|
| 1 | POTEH | 15695752 | BYPASS | OUTPUT_SPLIT |
| 2 | XKR3 | 16799808 | UNASSIGNABLE | UNASSIGNABLE |
| 3 | ENSG00000283809 | 18906691 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 4 | GSC2 | 19149726 | PTC_TERMINATION | PTC_TERMINATION |
| 5 | LRRC74B | 21053395 | OUTPUT_ABSENT | UNASSIGNABLE |
| 6 | GSTT4 | 24000153 | OUTPUT_ABSENT | WT_NOT_EXACT/NOT_RUN |
| 7 | CRYBB2 | 25227916 | BYPASS | WT_NOT_EXACT/NOT_RUN |
| 8 | CRYBB1 | 26607928 | OUTPUT_ABSENT | OUTPUT_SPLIT |
| 9 | CRYBA4 | 26625580 | WT_NOT_EXACT/NOT_RUN | BYPASS |
| 10 | NEFH | 29485794 | OUTPUT_ABSENT | OUTPUT_SPLIT |
| 11 | ENSG00000285404 | 31783939 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 12 | RFPL2 | 32193137 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 13 | SLC5A4 | 32232963 | OUTPUT_ABSENT | OUTPUT_SPLIT |
| 14 | APOL5 | 35726806 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 15 | CACNG2 | 36566405 | UNASSIGNABLE | BYPASS |
| 16 | TPTEP2-CSNK1E | 38300016 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 17 | APOBEC3A | 38961482 | OUTPUT_ABSENT | WT_NOT_EXACT/NOT_RUN |
| 18 | ENSG00000284554 | 39025618 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 19 | ENTHD1 | 39765284 | OUTPUT_SPLIT | OUTPUT_SPLIT |
| 20 | SHISAL1 | 44296746 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |
| 21 | TTLL8 | 50033351 | WT_NOT_EXACT/NOT_RUN | WT_NOT_EXACT/NOT_RUN |

非该工具WT精确集合者没有进行条件响应推理，不进入12/10分母，不记mutant阴性。公共原始JSON的conditional_response=PLANNED是原计划标记；实际完成由每条件response、execution receipt及终态技术计数确认。

### 明确绕行的实际结构

以下坐标为0-based half-open，长度直接来自已分配链的CDS边界：

| 工具/基因 | 编辑密码子被排除的方式 | 被排除的区间 |
|---|---|---|
| ANNEVO / POTEH | 原174nt内部CDS块被跳过；两侧原锚保留 | 原CDS [15695644,15695818)，新跨越内含子[15695485,15698661) |
| ANNEVO / CRYBB2 | 原内部CDS中出现27nt内含子 | [25227912,25227939) |
| Tiberius / CRYBA4 | 原内部CDS中出现15nt内含子 | [26625574,26625589) |
| Tiberius / CACNG2 | 原内部CDS中出现33nt内含子 | [36566378,36566411) |

这些是计算输出，不是RNA证实的exon skipping或新内含子。原合同没有事后加入更长内含子门槛；15/27/33nt结构照实呈现，不能由ORF合法推导生物学合理。

## 4. 解释与下一判断

两工具自身及共同集合上的识别下界均>0，支持“这个固定条件面板中存在PTC特异的严格结构绕行”。它不等于发生率估计、统计显著性、所有基因/物种普遍性，也不证明神经模型或ORF约束各自的因果贡献。

共同9中ANNEVO明确阳性为POTEH，Tiberius为CACNG2；前者在Tiberius是split，后者在ANNEVO不可归属。不存在两工具在同一共同位点均明确bypass的病例。因此，两工具均有响应不等于同位点机制复现。CRYBB2与CRYBA4分别只属于一个工具的WT精确集合。

Tiberius所有PTC目标的预定类别在两个末端过滤前后均未改变，没有记录到唯一归属目标链被这两步删除；两个bypass在过滤前已存在并保留到最终输出。这只排除“最终过滤制造该单链连接”作为这些输出的解释，不区分神经分数、HMM与自适应重推理。ANNEVO没有对应过滤轨迹，不能猜测其6个OUTPUT_ABSENT的来源。

## 5. 证据与产物

- 结果：[paired_result.json](paired_result.json)，原始冻结分类器终态输出的精确副本；全部21行和逐条件链/phase/坐标/分类均保留。
- QC/实现：[paired_tests.xml](paired_tests.xml)、[paired_input_manifest.json](paired_input_manifest.json)、paired_response.py及paired_run.py；测试在首mutant前完成。没有NGS reads，因此FASTQ、比对覆盖和湿实验QC不适用。
- 执行：[paired_native_execution.md](paired_native_execution.md)；完整命令、日志、原GFF/GTF、44个execution.json和20个trace.json均在Baobab的outputs/M27-PAIRED-NATIVE-R1/<tool>/<case>/。原始输入/模型/全预测不公开导出。
- 来源与边界：[WT结果](wt_result.md)、[原生配对合同](paired_native_contract.md)、[先行研究](novelty_scope.md)。本报告为Agent对实际产物的解读，不是独立重跑或Pro实验证据。

## 6. 新颖性与生物学限制

ACE/SGRF已描述完整性假设与损伤绕开；早期Helixer已有神经扰动分析。Pro第11轮在读取固定bb5f773的文献说明和合同后，进一步把M27定位为已知现象在现代原生流程中的有限复现。其提出的“神经证据与完整性约束冲突能否选择性识别不可靠修复”仍是未验证候选，不能从本结果推出。Pro未查看运行病例、未重跑；SGRF补充S1直读受限，其样本数未获第二次独立核验。

人工syn/PTC是两个不同碱基替换，不是纯粹终止密码子生物学效应。位点数小、参考先选、条件于WT精确、没有真实等位的独立转录/蛋白证据，未知训练暴露不称干净zero-shot。原生结构完整不证明基因功能完整，推断的内含子也不是已验证剪接。

## 7. 下一步

按首mutant前规则，有限重现确认对象已经由坐标顺序确定：ANNEVO的POTEH（ordinal1，ENSG00000198062.16）和Tiberius的CRYBA4（ordinal9，ENSG00000196431.4）。不是挑最漂亮、最长或最符合故事的例子，也不换成共同位点。

向Pro提交本完整阶段后，下一次最多进行这两个工具×目标的原生syn/PTC四次整chr22重算，检验目标响应与结构能否重现；不修改本次结果或分母，不直接展开全部分量组合。只有确认后另定的有限分量检验产生可定义、可否证的风险信号，才讨论合理纠错与真实损伤的双向验证。当前阶段完成不等于论文或Nature Communications目标完成。
