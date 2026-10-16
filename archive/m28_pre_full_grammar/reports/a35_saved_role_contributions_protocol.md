# 原固定角色分数的成分诊断

2026-10-06，实际读取logits前固定。A34已闭合，536个遗漏配对中190个非正、179个正分但未选、167个在参考辅助池获选。本次只分解原得分以定位负贡献来源，给下一独立模型目标提供依据；不改评分、训练、原选择器、阈值或A33停止条件。

## 原始总体和成分

源为实际GPU13403916与独立CPU13403920，全396窗／860配对、原536遗漏和324对照，全部保存FP32 logits。保留原已核验身份、几何、A34状态和gain，不选高分／低分样例，不重新读取feature cache／catalog或做GEN／模型前向。只使用旧六条TRAIN染色体，不访问DEV、Setaria、sealed、独立物种或A33输出。

原角色1至14互斥完整分为coding1–3、intron4–6、donor7–9、acceptor10–12、start13、stop14。使用原candidate_roles相位和短exon边界优先级；每组记录实际role base数及原role logit−intergenic logit贡献和，FP64累计。组和必须以rtol1e-12／atol1e-8复现原A34 gain。没有转移／长度先验、新权重、归一化评分或zeroing成分后的替代预测。

每参考记录最负的已存在组，仅在该组贡献<0时定义；同分按上列固定顺序。不存在组的每base贡献为null，不伪装成有证据的零值。汇总各原A34状态、两物种、非正组中的exon类别；pooled每base值用总贡献／总role base数，不平均不等长度参考的均值。参考／窗和base两种分母分开，重复窗口不当独立基因。

最负组是原数学贡献的描述，不是因果消融；负总和可能来自多组，组规模与每base强度需一起解释，不由统计直接归因类别不平衡、内含子长度或某一个生物机制。所有原状态与实际保存记录须完整复现；不计算性能gate、F1或新的gene调用。

## 执行范围

入口diagnose_saved_role_contributions.py，唯一新目录outputs/M28-A35-SAVED-ROLE-CONTRIBUTIONS，860行references／summary／STATUS／JOBID。三项小检查覆盖六组完整分数分区、缺角色的null每base语义、不同role base数的加权汇总。运行保存后重读实际全部记录并核对原190／179／167／324分母，不只是内存汇总。

private-teodoro-gpu／0GPU／2CPU／8GiB／5分钟，generanno，显式CUDA_VISIBLE_DEVICES为空；日志M28A35ROLEPARTS_JOBID.*。原A34 CPU完整多候选复算实际2:55，本次只拆860条参考，无新模型或自然池打分，5分钟为操作上限，不是速度benchmark。新文本预计低于16MiB，不写新logits／checkpoint，不覆盖任何原产物。故障保留输出，至多一次同范围、同资源时限的最小工程修复使用独立R1目录；科学参数不变。仅提交不算完成，终态和860真实记录／分数组和必须闭合；需要外部等待时只复用原固定监控／同一心跳。

## 与原生结构评分的区别

本次只读核对实际固定Tiberius2.0.8源码871d4979568ae09e56de7eca5d7489046ee833ea：site/tiberius/hmm.py的HMMBlock继承AnnotationHMM；eval_model_class.py的make_default_hmm使用VITERBI，predict_vit把神经预测和核苷酸输入交给该HMM。当前A31是局部角色相对intergenic的贡献和，再对既有完整候选做区间选择，不等价于该隐状态序列／转移模型的完整评分。该代码差异不证明性能差距的原因，不移植Tiberius默认长度去救冻结结果，也不把已有HMM／结构解码思想宣称为本项目创新。

## 闭合记录

CPU13404094实际COMPLETED／0:0，3:11，全部860记录及原分数成分和核对通过、最大误差0，三项检查PASS。主对话直接消费终态并独立汇总保存记录；同一心跳保持PAUSED。详细数字、解释与结论边界见[结果记录](a35_saved_role_contributions_result.md)，原固定协议未改。
