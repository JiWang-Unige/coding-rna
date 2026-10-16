# A31 完整 DEV 局部路径评分结果

2026-10-06。两条原生raw-FASTA注释及评价均已闭合。A31 A23在既有两条重复使用DEV contig上的合并完整CDS链F1为67.8133%，较原A3提高15.4934个百分点，高于本轮固定Helixer的62.8512%，仍低于Tiberius的71.8047%和ANNEVO的75.7025%。这是实测开发进展，不是独立泛化、清洁zero-shot、统计显著性或论文竞争性通过；两个variant的原始结果都保留，不由这次DEV数值选择新测试模型。

## 闭合证据与共同口径

主对话在09:05 UTC取得原始Slurm终态。以下时间统一为UTC；每个root、batch、extern均COMPLETED／0:0，输出STATUS／JOBID匹配。

| 阶段 | Job ID | 实际开始至结束UTC | Slurm用时 | batch MaxRSS |
| --- | --- | --- | --- | ---: |
| A31 A3原生推理 | 13399434 | 05:37:41至08:51:28 | 3:13:47 | 3101384K |
| A31 A3评价 | 13399454 | 08:51:29至08:52:12 | 0:00:43 | 1414820K |
| A31 A23原生推理 | 13399441 | 05:39:11至08:49:27 | 3:10:16 | 3099912K |
| A31 A23评价 | 13399455 | 08:49:28至08:50:37 | 0:01:09 | 1400260K |

沿用[事前固定DEV协议](a31_local_path_caller_dev_protocol.md)：拟南芥NC_003074.8／23459830bp、水稻NC_089041.1／29936421bp，正负方向四scope，各variant均8690窗，真实尾窗及同一四块GEN上下文。两者分别使用原A3、A23末004608，均为4608次原更新；局部15类角色logit相对intergenic0求和、所有角色权重1、null0及FP32／FP64配方不变。本次没有新参数、训练更新、阈值扫描、feature cache或Setaria访问。

CPU评价在打开参考前重建实际自然池、owner、DP选择及非空GFF坐标／phase，保存候选计数和最终预测数完全一致。A31 A3逐窗几何、顺序和greedy／beam来源与原生A6相同；A23的同池字段为JSON null，不要求异生成器保持同池。共同primary分母为拟南芥5437、水稻2291，共7728，unsupported、难例及原背景范围均保留。旧6450条parent-primary指标不是本表口径。

完整链TP要求同strand、相同有序CDS坐标链；phase另行诊断，不属于exact-key。Precision=TP／预测唯一链，Recall=TP／primary参考链，F1=2TP／（预测唯一链+参考链）。Micro先合并整数，Macro才取物种均值。主对话由两份原summary的整数独立重算全部对照P／R／F1、pooled／macro、stage TP及背景比例，均一致。保存的五状态诊断又复现原候选、owner、选择和评价计数。

## 完整链性能

下表同一两条DEV、同一7728条primary，而非原论文数字拼接。原生对照资产与暴露边界见[native_baseline_result.md](../../../benchmarks/native_baselines/reports/native_baseline_result.md)及[pretrained_exposure.md](../../../paper/supplement/pretrained_exposure.md)。

| 方法 | 拟南芥F1 % | 水稻F1 % | 合并Micro F1 % |
| --- | ---: | ---: | ---: |
| A3 | 51.7567 | 53.5236 | 52.3200 |
| A6 | 49.1465 | 45.0061 | 47.9253 |
| A31_A3 | 67.5419 | 63.0786 | 66.1671 |
| A31_A23 | 69.4210 | 64.1419 | 67.8133 |
| Helixer | 67.4564 | 53.5243 | 62.8512 |
| Tiberius | 74.1771 | 66.5822 | 71.8047 |
| ANNEVO | 76.2651 | 74.4318 | 75.7025 |

| Variant／范围 | TP／预测／参考 | Precision % | Recall % | F1 % |
| --- | --- | ---: | ---: | ---: |
| A31 A3／拟南芥 | 3524／4998／5437 | 70.5082 | 64.8152 | 67.5419 |
| A31 A3／水稻 | 1465／2354／2291 | 62.2345 | 63.9459 | 63.0786 |
| A31 A3／合并 | 4989／7352／7728 | 67.8591 | 64.5575 | 66.1671 |
| A31 A23／拟南芥 | 3621／4995／5437 | 72.4925 | 66.5992 | 69.4210 |
| A31 A23／水稻 | 1465／2277／2291 | 64.3390 | 63.9459 | 64.1419 |
| A31 A23／合并 | 5086／7272／7728 | 69.9395 | 65.8126 | 67.8133 |

A23相对A31 A3的合并提高1.6462个百分点：拟南芥TP增加97，水稻TP仍1465，但水稻预测数减少77。A23相对Helixer的合并优势4.9622个百分点，两个物种均数值较高；相对Tiberius仍差3.9914、相对ANNEVO仍差7.8892个百分点。A31 A3在拟南芥仅比Helixer高0.0855个百分点，不将其小差异描述成稳健优势。原pilot对A及C0的辅助门槛均PASS，不能替代上述强基线差距或作为论文PASS。

## 碱基与背景负担

CDS-base按不分strand的CDS union计算，不衡量整基因连接。背景是所有注释gene feature跨度之外的bp，拟南芥8183405、水稻19773036、合并27956441；这里的背景比例不是逐基因分类FPR。

| Variant／范围 | CDS-base F1 % | 背景预测span bp | 背景span覆盖 % |
| --- | ---: | ---: | ---: |
| A31 A3／拟南芥 | 91.9099 | 39231 | 0.4794 |
| A31 A3／水稻 | 89.9445 | 166472 | 0.8419 |
| A31 A3／合并 | 91.2769 | 205703 | 0.7358 |
| A31 A23／拟南芥 | 93.2850 | 48171 | 0.5886 |
| A31 A23／水稻 | 90.0690 | 170553 | 0.8626 |
| A31 A23／合并 | 92.2631 | 218724 | 0.7824 |

A23并非所有维度都优于A31 A3：背景span更多。两个variant全部保存预测均通过现有ORF／splice／phase结构检查，没有重复坐标链；这表示满足调用器的结构约束，不表示所有预测在生物学上正确。

## 五状态损失及下一机制依据

事前[五状态诊断协议](a31_dev_score_attrition_protocol.md)在全部五个source root闭合后执行，作业13402231：private／0GPU／2CPU／8GiB／5分钟，实际09:07:06至09:08:20 UTC，root／batch／extern COMPLETED／0:0，用时1:14、MaxRSS765784K。四项合成检查PASS／0.003秒，STATUS／JOBID相符，references.tsv为7728数据行加表头。诊断没有模型／GEN前向、新预测、R4标注、训练、参考删例或阈值变化。

每条primary参考恰落入一种状态，下表每行合计7728；selected与上表TP完全一致。主对话还独立重算每物种、CDS片段层级与A6状态迁移的边际合计，均复现原计数。

| Caller | 未进入自由池 | Owner丢失 | 非正gain | 正gain未选 | Selected |
| --- | ---: | ---: | ---: | ---: | ---: |
| A6 | 1761 | 175 | 2967 | 157 | 2668 |
| A31_A3 | 1761 | 175 | 261 | 542 | 4989 |
| A31_A23 | 1654 | 198 | 255 | 535 | 5086 |

A31 A3与A6逐参考free／owner成员关系相同，正确链的非正gain损失从2967降至261，selected从2668升至4989，同时正分未选增加至542。至少五个CDS坐标片段的2974条参考中，A6仅选144条，A31 A3选1694条，A23选1792条；这一层的非正gain数分别1768／11／8。这是保存输出的状态差异，不是匹配容量、匹配训练曝光的纯评分因果实验，也不将CDS片段数直接当作基因bp长度。

A23仍未命中的2642条参考中，1654条未进入自由池，占62.6041%；535条正gain未选，占20.2498%；255条非正gain，占9.6518%；198条owner丢失，占7.4943%。正分未选不等于真实基因彼此冲突，现有诊断没有确定其唯一原因。相对A31 A3，A23的free参考命中在拟南芥4148→4269，水稻1819→1805；不能把合并增益表述成两物种生成覆盖均改善。

下一机制优先针对完整多片段正确候选的生成／保留，以及正分候选间的竞争。非正gain已不是当前最大损失，因此不重复旧质量头或靠DEV threshold扫描救成绩。[三格保存输出DEV诊断](a31_score_pool_dev_diagnostic_result.md)13402574随后COMPLETED／0:0：固定A3的greedy＋query、greedy＋local、union＋local合并F1依次52.3200%／56.4037%／66.1671%。评分变化使非匹配primary预测减少1060、TP净减少10；固定局部分数扩池使TP增加1261、至少五CDS片段正确链增加705。全部8690窗几何匹配，但缓存与原生hidden数值未证明相同，仍只作条件解释；新训练机制另行固定设计，不升级为独立确认，也不执行因A30原五gate FAIL而停止的A32。

## 适用范围与可追溯产物

可用于限定DEV内部判断，不能用于清洁zero-shot、独立泛化、显著性、因果机制或论文验收。当前只评两条contig，不是两个完整基因组。GEN上游暴露unknown，Helixer已知训练暴露，ANNEVO作者排除声明与Tiberius物种暴露按固定资产分别保留。A3主对照同时变池与分数，A6多了质量头训练曝光，A23还改变生成器与拟合共享权重；不是纯单组件效应。

两个原生入口运行段分别11580.02／11399.20秒，峰值GPU allocated2455143424／2455144960字节，实际均单RTX3090／2CPU／16GiB；只是本次操作收据。旧缓存A3、原生三对照的CPU／memory／流水线条件不同，未作公平速度benchmark，不据峰值或运行段称效率优越。

原固定watcher在09:03:20.685 UTC完成一次实际结构交接后，原heartbeat已PAUSED；这只停止已消费作业的等待，不暂停研究goal。三对照短人工FASTA预检此前只提出、没有文件或作业提交；正式结果到达后优先消费主结果，未把它冒称运行完成。本结果未提交Git、公开发布或打开独立面板／Setaria。

原始数据位于canonical远程/home/users/j/jwang/coding-rna：outputs/M28-A31-A3-LOCAL-PATH-DEV-INFER／DEV-EVAL、outputs/M28-A31-A23-LOCAL-PATH-DEV-INFER／DEV-EVAL各自summary.json、四scope JSONL、GFF3、STATUS和JOBID；诊断在outputs/M28-A31-DEV-SCORE-ATTRITION/{summary.json,references.tsv,STATUS,JOBID}。实现为evaluate_local_path_caller.py及diagnose_local_path_caller_dev_attrition.py，日志在outputs/M28-PILOT-R5/logs的对应job文件。保存Markdown未进行渲染QA，不作版式验证声明。
