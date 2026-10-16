# 局部角色输入模型的TRAIN负结果

2026-10-07北京时间。A40拟合、自然读出及七项CPU核验全部完成，工程与数值证据闭合；原五项继续条件通过3项，科学FAIL关闭。水稻正确配对增加1不能抵消拟南芥和长链的失败，不进入DEV，不挑checkpoint、加曝光或救阈值。论文总目标仍未达到。

## 较强同池基线与原条件

粒度为原1491唯一TRAIN窗的自然链配对，不是独立gene或完整染色体F1。原16321候选、原标签、gain0及span非重叠DP均保持，预测没有参考注入。

| 原继续条件 | A31_A3 | A40 | 判定 |
| --- | ---: | ---: | --- |
| 拟南芥正确配对严格增加 | 1531 | 1484 | FAIL，减少47 |
| 拟南芥known_nonmatch不增加 | 347 | 120 | PASS |
| 水稻正确配对严格增加 | 478 | 479 | PASS，增加1 |
| 水稻known_nonmatch不增加 | 251 | 106 | PASS |
| 合并至少5-exon正确配对严格增加 | 700 | 636 | FAIL，减少64 |

完整reference/window分母仍为2172和744，合计2916；至少5-exon为854和232，合计1086。原自然池正确配对1751和629、长链正确配对合计788不替代参考分母。拟南芥未知选择170，水稻47，合计217；known_nonmatch是primary参考相对已知非匹配，未知不计为错误或新基因。

合并已选2406，其中正确1963、已知非匹配226、未知217；较A31_A3的2009、598、447分别变化−46、−372、−230，长链700至636。较已闭合A39，正确1974至1963、已知非匹配232至226、未知225至217、长链653至636。这里只描述同池选择组合，不能用非匹配下降称整体性能改善，也不能从不同容量和初始化随机流的单seed比较推断角色输入的单一因果。

## 实际工程及数值证据

fit13418488、natural13418506、CPU13418507的root、batch、extern均COMPLETED0:0，耗时1:08:53、19:12、7:08。fit和natural实际private1RTX3090、2CPU、16GiB；CPU实际private0GPU、2CPU、8GiB。batch MaxRSS依次2108744K、1699516K、866148K，STATUS、JOBID及源作业相符。CPU于2026-10-06 22:34:13 UTC结束。

fit实际4608行，三个checkpoint均10014552字节；natural_scores和windows各1491行，恰1491个非空window_000000至window_001490局部发射文件。七项实际证明全部true：原draw、importance及4608更新；固定等权复合loss算术；三student checkpoint及真实optimizer曝光；三冻结源逐tensor等于原A3并排除于optimizer；保存角色输入及原A31分数；FP32分量算术和自然身份；原标签DP、较强基线、完整分母及五条件。保存原A31分数最大误差0。

实际SharedContext与window null的参数step为4608，其他非连接头为4602，junction两分支为4425，不把所有参数说成4608。主逐项重算两物种与pooled的五个严格整数条件，并核对读出与CPU的完整counts、分母、配置及实际参数曝光一致。CPU没有新模型前向；这是保存记录和决策核验，不冒充独立神经重放或生物学验证。

原证据根为outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-FIT、M28-A40-ROLE-EVIDENCE-CONTEXT-TRAIN-READOUT及M28-A40-ROLE-EVIDENCE-CONTEXT-RESULT-CHECK，原[事前协议](a40_role_evidence_context_protocol.md)和所有产物保留。固定监控原生暂停并读回同ID、heartbeat、原固定对话、15分钟和通知字段；主在下次定时检查前手动消费终态，不声称完成了自动终态回传。

## 下一步的研究范围

只读转移诊断13420717已完成：2026-10-06 22:50:51至22:51:25 UTC，root、batch、extern均COMPLETED0:0，34秒；private0GPU、2CPU、8GiB，batch MaxRSS573760K。两项新adapter检查PASS，0.009秒。原1491窗、16321候选、三臂counts和2916/1086分母复现；主核对两比较、两物种及pooled、四标签集合共24组的六类互斥转移守恒。

相对A31_A3，正确配对保留1841、丢失168、新增122，净−46；丢失中131为新gain非正、37仍正分但原DP未选。至少5-exon正确保留607、丢失93、新增29，净−64；丢失中72非正、21正分未选。相对A39，正确保留1858、丢失116、新增105，净−11；长链保留597、丢失56、新增39，净−17。原unknown继续单列，不计为错误。低分丢失与正分未选同时存在，仅为保存分数的描述性定位，不证明角色输入的单一因果，也不授权改变DP、阈值或loss。诊断未打开模型、特征tensor、TRAIN参考catalog、新标签或DEV，没有前向和optimizer更新。证据根outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-ATTRITION。

A41的零初始化冻结发射修正在读取此结果前已固定并通过零更新预检。该独立模型保留初始A31最终势，只学习固定上下文的1806个仿射边际参数，不是将旧A40分数混合、继承其student或救A40。其完整接线、一次固定拟合及原五条件见[A41协议](a41_emission_residual_protocol.md)。接线13420720实际通过后，A41已一次提交fit13420820及afterok natural13420821、CPU13420824；这是独立模型运行，不改变A40的FAIL。raw-FASTA竞争性能、独立泛化和论文主张继续未完成。
