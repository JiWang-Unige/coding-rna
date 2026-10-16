# 同更新数 TRAIN 参考链与竞争链评分诊断

2026-10-11。本方案在读取本轮 native-prior、fresh、joint DEV 指标之前固定。原 TRAIN 中联合24576编码的末9216模型 pooled 链F1为66.963%，同更新数结构续训为55.745%；两者仍低于ANNEVO和Tiberius。已有全范围评分分解仅覆盖原4608模型，且没有所选竞争链的角色项。这里重放两个已完成末9216模型的实际生产势函数，区分剩余遗漏的背景负分、竞争压过和固定语法不支持，并补齐竞争双方的角色分数。它不是新训练、新操作点、重新预测或因果模块验证。

## 固定模型与范围

对照为STRUCTURED_CONTINUE final9216（TRAIN评价13451625，原6144独立编码）；候选为JOINT24576 final9216（TRAIN评价13453780）。两者来自共同4608起点，使用同原1536 draws、1491 unique views和4608新增更新。只读取这两个已完成正式末权重、原六条TRAIN的冻结raw输入与预测、原共同尺22,688条primary参考记录，以及既有完成收据。保留17,739条拟南芥和4,949条水稻参考，包括209条语法不支持参考；不访问DEV、held-out、sealed、Setaria或新物种。

两条生产接入的checkpoint验证、BF16编码后FP32头、完整24576窗口与12288 stride、nearest-center ownership、双链定向、15角色单位权重、null0、cap128和FP64逐base累积均不变。不用参考改编码窗口、选择权重、输出边界或阈值。参考仅在势函数生成后用于固定路径分数计算。重放不运行新的decoder，不写GFF，不保存大特征或整染色体势函数。

## 评分与配对口径

沿用已核验的实际DNA状态图支持检查及完整角色路径。每个支持参考的gain为整个CDS外侧span内各角色logit减背景logit的FP64和，包含内含子。其同strand、span正重叠的全部冻结所选链构成删除集合；参考gain减删除集合gain为合法单链替换margin。相接不算重叠，反strand不冲突。不支持的参考保留并报告，角色项为空。

所选每链gain须在原八位有效数字GFF舍入界加1e-5内复现，整strand所选gain须与原raw收据在relative1e-8／absolute1e-5内一致。1e-5仅隔离报告中的近零／tie，不是新的操作点。若存在能严格提高冻结所选集合的合法参考替换，停止科学解释，先核对评分重放、坐标和decoder一致性。

保存参考、所删竞争集合和二者差值的15角色gain及base-row数。四组互不重叠：coding body为1–3，intron body为4–6，splice coding boundary为7–12，start／stop为13–14。角色分解须重构原gain和替换margin；base-row是链暴露量，不是unique genomic coverage。多个竞争链的角色项相加，不把只与参考局部相交的部分截断成新分数。

以固定参考的物种、seqid、strand和完整CDS坐标配对两个模型，报告全染色体、两物种及pooled的类别迁移和角色margin，另保留双方均命中、仅对照命中、仅候选命中、双方均遗漏四个固定结果组。TP、预测、参考须重放原评价；候选净命中增量应为2,230。分组是事后结果描述，不是独立因果层、subset precision／F1或新基因发现。不同头的gain尺度未校准，跨模型的绝对gain变化不解释为正确概率或可比较的置信度。

## 执行与下一步

直接复用现有生产编码和已核验的评分／状态图函数，仅新增两个固定模型的接入、竞争角色项和配对汇总。先CPU人工检查新增的多竞争链角色加和、margin重构、空删除集合及错误和的拒绝；无研究数据或原生forward。通过后原六条TRAIN各一个GPU任务，顺序重放两个模型，整数组afterok CPU汇总。单task单24 GiB RTX3090、2CPU／16GiB、90分钟上限；CPU检查及汇总为2CPU／4GiB／5分钟，使用已有generanno环境。既有4608评分诊断最长28分钟，候选原raw的encoding约26分钟且本次省去decoder；两个模型顺序重放并留小幅运行余量。整个产物预计小于256 MiB，不保存checkpoint、DNA副本或大缓存。

输出为独立`outputs/M28-FULL-GRAMMAR-MATCHED-ENERGY-CHECKS`及`outputs/M28-FULL-GRAMMAR-MATCHED-ENERGY-DIAGNOSTIC`，不覆盖历史诊断、预测或冻结结果。每task保存两个模型的reference_energy.jsonl、summary、started和STATUS／JOBID；CPU保存全范围summary及简明comparison。正常执行与明确同科学设置的实现修复按既有自主授权推进；科学负结果不触发重训或规则修改。

结果只判断9216模型的剩余缺口是否仍由实际评分解释、候选相对匹配对照恢复／丢失的参考对应何种角色差值。不能由分数算术证明上下文的机制因果、独立泛化或论文竞争成功。冻结八方法DEV和独立十方法扩展均保持不变，不新增候选、不根据DEV反馈修改训练或本方案。长等待只复用原监控对话与唯一心跳，不创建新监控或重复自动化。原论文goal继续保持未完成。

## 实际检查与提交

方案于2026-10-11 UTC13:30:40固定，依据是已消费TRAIN结果；当时本轮DEV作业已结束，但主尚未读取native-prior／fresh／joint新DEV指标。固定不等于DEV运行前预注册，读取指标后没有改变这里的模型、范围、评分或分组。

人工CPU检查13494221于13:35:19–13:35:32执行，root／batch／extern均COMPLETED／0:0，stderr空、STATUS／JOBID／summary匹配。多竞争角色项加和、margin重构、空删除集合、错误和拒绝、相接／多重交叉、四hit组与不支持参考保留通过；research_data_read／checkpoint_read／native_forward／DEV／sealed均false、更新0。原生产geometry返回tuple而原JSON保存list，实际GPU前仅把ownership比较改为相同表示；未改科学设置或已核验算术。

提交GPU数组前再次由sacct和完成收据确认检查成功。尝试对该已完成检查建立afterok被调度器以Job dependency problem拒绝，没有产生GPU作业；确认通过后直接提交实际数组13494264_0–5，max6，每项单RTX3090／2CPU／16GiB／90分钟。汇总13494270仍依赖整个新数组afterok:13494264_*，0GPU／2CPU／4GiB／5分钟；没有放松数组或产物完成条件。

六项均于13:41:51 UTC启动在gpu034，实际numeric roots依task0–5为13494265／13494266／13494267／13494268／13494269／13494264；JOBID和started均保存numeric root。主核对全部started的原seqid／长度／物种、两个固定末权重、source evaluation13451625／13453780和本协议，reference_used_by_encoder_or_head=false、更新0、DEV／sealed=false。启动时六stderr均空，CPU汇总PENDING有效整数组Dependency。这是实际接入启动，不是评分回放或诊断成功；终态和科学结果待整批完成。

等待仅由原固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3和唯一心跳coding-rna-r5-light-monitor负责。13:49:00.613 UTC保存本批唯一事件M28-MATCHED-TRAIN-ENERGY-13494264-13494270；主直接读回ACTIVE／heartbeat／原target／15分钟／原通知策略／created_at未变，prompt与提交文本逐项相等。首次人工请求turn01a12b3a-24a9-7601-b31c-91db1860fb84实际运行13:49:53–13:50:25，turn_context确认gpt-5.6-luna／low；13:50:18核对六GPU RUNNING、汇总PENDING有效Dependency，不读科研内容、不交接、不暂停。上一DEV批13:40:58的真实自动交接已消费，不再轮询。

本批真实定时追加已于14:04:15–14:04:41 UTC在同一原监控对话运行，turn01a12b47-4bde-72a3-9c7e-184c897739ec。主直接读取实际heartbeat输入：automation_id为原唯一心跳，调度时间14:04:15.197Z，输入包含本批唯一事件；实际turn_context确认为gpt-5.6-luna／low。14:04:34检查六GPU仍RUNNING、CPU汇总PENDING有效Dependency，无科研内容读取、交接或暂停。配置、人工首检及本批定时追加分别验证；本批终态回传和完整评分诊断仍待完成，未记为成功。

## 完整终态与评分结果

六个GPU任务及CPU汇总13494270的全部21个root／batch／extern均COMPLETED／0:0。GPU最晚在2026-10-11 UTC14:39:41结束，汇总于14:39:41–14:40:29执行。主已逐项核对六目录STATUS／JOBID／外收据与十二臂收据、真实numeric roots、两个固定source fit／evaluation、末9216、原seqid／长度、零更新及DEV／sealed=false；产物齐全。每条冻结所选链的分数通过八位有效数字GFF舍入界检查，24个strand所选总分与原raw收据差值均为0，最大单链舍入差为0.0000498786。未发现任何严格提高已选集合分数的合法单参考替换；这是一致性诊断，不是另一次全局MAP证明。

主独立重算六染色体TP／预测／参考加总及全部三scope的PRF、四hit组分母、类别迁移的行列边际、96个角色均值差与按supported人数加权的pooled均值。全部与原评价和保存summary一致，未改预测、阈值、训练或评价规则。全范围参考22,688（拟南芥17,739／水稻4,949），209条单链语法不支持参考完整保留；unsupported移到其他类别为0，近零／replacement tie类别本批为0。原始路径及逐链字段位于outputs/M28-FULL-GRAMMAR-MATCHED-ENERGY-DIAGNOSTIC，各臂reference_energy.jsonl由已完成Slurm汇总配对，主未在登录节点重新处理这些大记录。

| Scope | Structured TP | Joint TP | Gained | Lost | Net |
| --- | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis | 10,837 | 12,284 | 2,060 | 613 | +1,447 |
| Rice | 1,859 | 2,642 | 896 | 113 | +783 |
| Pooled | 12,696 | 14,926 | 2,956 | 726 | +2,230 |

## 剩余遗漏与类别迁移

| Scope | Structured 背景负分 | Joint 背景负分 | Structured 正分被竞争压过 | Joint 正分被竞争压过 | 双方单链不支持 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis | 2,629 | 1,667 | 4,165 | 3,680 | 108 |
| Rice | 1,727 | 819 | 1,262 | 1,387 | 101 |
| Pooled | 4,356 | 2,486 | 5,427 | 5,067 | 209 |

Joint的7,762条遗漏由2,486背景负分、5,067正分被全部所选冲突链压过、209单链不支持构成；前两类7,553／7,762=97.307%。这里说的是单链可表达，不是所有参考可同时共存。既有全参考兼容性oracle还包含19条额外联合约束遗漏，不能把本表的97.307%直接解释为全部可恢复。

2,956条新命中中，2,123（71.820%）来自旧模型的正分竞争遗漏、833来自背景负分。另有1,085条从负分变成正分竞争遗漏，但没有成为exact hit。726条旧命中被丢失，704转为竞争遗漏、22转为负分。背景负分类别减少1,870，不等于恢复1,870条参考；竞争遗漏只净少360，水稻反而增加125。这比只比较两个模型的类别总数更准确。

## 配对结果组的角色 margin

以下是支持参考的逐链平均原始score单位，margin始终为完整参考gain减完整所选冲突集合gain；不是对局部重叠片段重计分。四项用未舍入值精确重构total。

| 固定结果组／评分模型 | 支持链数 | Coding body | Intron body | Splice boundary | Start／stop | Total margin |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 新命中／旧Structured | 2,956 | +5.1171 | −29.0444 | +2.9378 | +2.4963 | −18.4933 |
| 丢失命中／Joint | 726 | −3.0811 | −11.0849 | +5.5680 | −0.0960 | −8.6940 |
| 双方遗漏／Joint | 6,827 | −7.2952 | −34.7842 | −2.4944 | −2.8508 | −47.4246 |

双方遗漏的Joint内含子项在拟南芥为−20.1887，在水稻为−67.7967；水稻splice项仍为正+3.4700。内含子项是这些组平均margin的最大负项，但它是跨不同长度链的未归一化分数和，不能据此断言每条漏链的主要原因，或推导降低内含子权重的因果效果。不同训练头的gain尺度未校准，跨模型绝对值变化不是概率、confidence改善或效应大小。新命中在Joint下的margin为0是已选同链与自身相减的恒等式，不是机制验证。

本结果可用于探索性TRAIN评分核算。它支持继续区分表示／学习得到的评分与语法容量，而不支持因果、独立泛化、端点重调或论文成功。单fit、事后hit分组、有限训练视图及注释预训练暴露边界仍在；不增加训练、不重跑caller、不放行封存数据，冻结DEV家族不变。正在执行的fresh×joint第四格属于其独立既定协议，不是本评分诊断产生的新调参。

## 终态交接与剩余监控

原监控对话在2026-10-11T14:49:55Z发送本事件唯一真实定时终态交接，主已收到并完成上述科学消费。该scheduled turn01a12b70-f440-7ac1-bbf9-70b40ea4af07结束后原心跳PAUSED；主于2026-10-11T14:56:34.801Z通过同一automation工具恢复ACTIVE，仅继续M28-FRESH-JOINT24-TRAIN-13494895-13494899。保存值已精确读回：heartbeat、原target、15分钟、通知策略及created_at不变，原B完整结构白名单保留，不再检查本事件。B原范围的人工和定时检查已有实际证据；本次B-only恢复后的追加尚未单独观察，不另发重复首检。没有新监控对话、重复自动化、cron或模型变更。论文goal仍未完成。
