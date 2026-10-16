# 原 TRAIN 双模型漏检评分诊断

2026-10-13。在完整 CDS 后验排序之外，需要识别 E rank4 与 M dense 的原 TRAIN 漏检究竟对应何种评分失败。本轮将复用同批完成后的整染色体角色势，对全部22,688条原主参考链做只读诊断：区分固定语法不支持、参考链弱于背景，以及合法正得分参考被已选结构竞争压过。结果将指导下一步应检验评分学习还是结构覆盖；不能仅凭分组算术归因于rank、学习率、编码器容量或生物学创新。

## 固定模型与参考范围

固定 E fit13499817、M fit13510324，head9216／native4608，原六条 TRAIN 染色体与现有后验批次相同。只复用 `M28-FULL-GRAMMAR-GENE-POSTERIOR-TRAIN` 双strand完整 FP32十五角色数组，不重新编码、加载checkpoint、训练、改变预测或新增数组。CPU在原汇总评价13548940完成后才启动，保持先完成十二个后验及原参考排序、再做诊断的访问顺序；不改任何现有25个作业或参数。

参考来自已完成13449803的原 TRAIN 主链profile，沿用M26 CDS-assessable最长eligible CDS／transcript ID tie政策。逐染色体reference／TP与 E13499823、M13510351以及现有后验评价精确复现；保留 unsupported 和长于24576窗口的参考。长span不是whole-contig语法自动不支持。原profile的1536-draw R4 metadata只描述历史source，不用它定义 E/M 暴露、训练eligibility或独立测试，不据历史DEV结果选择机制。

## 固定评分与可替换事件

沿用实际carry、GT／GC-AG、min-intron5和每gene cap128。对每条参考验证单链在literal DNA生产语法中是否可达；对合法链按原角色precedence计算 `gain = Σbase (role_logit − background_logit)`。使用既有生产支持／逐base FP64前向累积／角色分解helper，不建设另一套解码语法。

在同strand原选set中，删除所有span与参考span重叠的已选gene，再插入这一条合法参考。原选gene互不重叠，触碰端点不算重叠；因此 `replacement margin = reference gain − sum(conflicting selected gains)` 是该具体替换的精确评分变化。这个替换不要求整个参考set联合可达，不把错误参考的人工插入当新预测或新的召回结果。跨strand不加排斥。

既有1e−5只用于数值近零报告，不是科学阈值或新operating point。分类保留 exact hit、outside grammar、negative／near-zero versus background、positive but outscored、near replacement tie，以及异常的positive improving replacement。最后一类若存在，保留逐链证据和失败summary并停止解释：它可能说明原MAP／输入重放有问题，不能包装成模型性能改善。

每gene得分重建必须与原GFF八有效位得分在relative6e−8／absolute1e−5内一致；各strand选set总得分也重建原scope收据。保存参考gain、冲突gain、margin及十五角色分解；编码body、intron body、splice coding boundary、start／stop组互斥。按六染色体、两物种及pooled报告全部分类、单／多CDS段和长span层，以及 E/M共同命中、E-only、M-only、共同漏检和分类迁移。原主链TP变化保持与既有评价一致。跨模型原始gain尺度未校准，组均值只作各自描述，不用于rank或LR因果推断。

## 必要验证与运行

新增四组检查对应具体失败：paired hit／分类迁移口径；unsupported和长span保留以及缺行／重复／role几何不一致拒绝；端点触碰与多冲突角色margin；实际carry／motif／cap支持与角色总分。资格同时读取原 E/M source和profile的完成metadata，不读取参考profile正文、研究DNA、真实参考GFF、checkpoint或native。它只证明接线和口径，不证明真实漏检分布。

资格2CPU／4GiB／5min，无GPU。正式一个2CPU／16GiB／1h的CPU作业顺序处理六染色体和两模型，private-teodoro-gpu／generanno／no-requeue，仅afterok13548940；不新增GPU用量。源角色数组已经包含在既有44GiB批次声明，新增逐链JSONL／summary声明256MiB上限，不保存第二份数组或DNA。正式完成需root／batch／extern全部COMPLETED0:0、STATUS／JOBID／started／总summary和十二model×seqid逐链文件／summary齐全；全部22,688参考的E/M逐项几何和支持一致、原评价重建、没有异常提高原set得分的合法替换。

确定性接线错误至多一次同范围修复并保留失败输出；数学或科学负结果不改定义求漂亮成绩。监控仍复用原watcher和heartbeat，原LR事件不变，新诊断作为POSTTRAIN原25阶段后的末端扩展；部分后验、排序或诊断完成不代表完整扩展事件成功。正常状态静默，只在完整终态或需处理失败向原main一次接续。不开新chat／automation／cron，不改冻结DEV、sealed／held-out、新panel、goal状态或公开发布边界。论文竞争性模型和独立验证目标仍未达成。

## 实际运行记录

新Python入口／测试AST与两份SBATCH语法通过，五份新源／协议精确同步并回读后才提交资格。CPU13549411于04:02:51–04:03:15UTC（北京时间12:02–12:03）、24秒，root／batch／extern全COMPLETED0:0，STATUS／JOBID／summary齐全。四组新检查、原E/M completed source与主profile metadata分母接线通过，0errors／failures、程序4.868949秒；没有研究DNA、真实参考GFF、profile正文、checkpoint或native读取。主已消费，watcher handoff=NOT_SENT_BY_WATCHER，不补发资格交接。资格通过不是实际漏检分布或机制结论。

正式13549479于04:08:00UTC（北京时间12:08）接受，命令 `sbatch --parsable --dependency=afterok:13548940 sbatch/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE.sbatch 13549411 13548940`。04:09:23UTC实际PENDING有效Dependency；controller核对2CPU／16GiB／1h、无GPU、private、Requeue0、完整实参与独立log成立。原13548940仍正常等待十二后验，own expanded55加本任务仍小于80；不取消、修改或重提原25作业。

原POSTTRAIN等待批次现扩展为原25阶段→评价13548940→诊断13549479，实际26roots／78root-batch-extern；13548940单独成功为中间阶段。原LR13545908仍独立，两个handoff均NO。相同heartbeat已保存完整13,434字符新prompt，ID／heartbeat／ACTIVE／原watcher目标／15分钟／通知偏好及created_at精确回读一致；新诊断结构白名单只核验文件和完成身份，不让watcher读取gain、margin、分类或解释科学结果。

实际定时turn01a13379-5f59-7573-aa8a-d85329c1bdc5于04:15:54–04:16:09UTC追加到原watcher，完整13,434字符instructions一致、当前turn context为gpt-5.6-luna／low。仅查LR末端和新26POST IDs：两个末端均有效Dependency、现有评价仍等待、部分GPU编码与六个CPU后验正常运行；没有读取科学数字或向main交接。扩展scope配置与一次真实正常状态检查成立，终态消失处置及完整自动回传仍未验证。尚无真实漏检评分结果，完整竞争性模型和论文目标仍未完成。
