# GENERanno 原生输出的比较口径

2026-10-03。静态核对作者固定29a93db37ec8b9a44e41502de5dbb4e5cf630830的[CDS调用器](https://github.com/GenerTeam/GENERanno/blob/29a93db37ec8b9a44e41502de5dbb4e5cf630830/src/tasks/downstream/cds_annotation.py)。

源码事实：`process_sequences_on_gpu`返回各head二元注释；`write_outputs_for_input/write_parquet`保存双链标签FASTA和Parquet数组，不输出多外显子Parent分组或CDS phase。`calc_acc`的exact_f1检验连续CDS片段边界及区间内容，不是整条剪接链。

项目判断：`common_ruler.py::score`比较完整有序CDS坐标链，二元mask不能直接提供该输出。原调用器的链指标应标为不适用，不记零分；既有ANNEVO/Tiberius/Helixer主比较不变。原发布head以后可另作CDS-base／片段对照；若接组装器，须明确署名GENERanno＋该组件，不能称作者原生链输出或按参考补连接。

这是调用器范围核查，不证明模型内部没有结构信息、性能差或所有GENERanno版本均无链输出。本轮未运行原head、生成新成绩、改现有方法／分母、访问test/Setaria或提交作业。原权重身份与上游未知暴露仍见[来源记录](pretrained_exposure.md)。
