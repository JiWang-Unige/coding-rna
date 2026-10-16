# M28 R3 GPU头吞吐：完成，但尚未训练

2026-09-30。Job **13292011 COMPLETED 0:0**，gpu021 RTX3090，
2CPU/16GiB，wall35s，MaxRSS1,782,104KiB，allocated0.009722GPUh。
没有optimizer step。12个现存TRAIN缓存两臂各前向/反向一次，无GLM重提取。

C0有366,191个可训练参数，B1有565,495个；GLM不在可训练参数中。
FP32头峰值allocated C0 165,360,128B、B1 184,058,880B；此处不含已离线缓存的GLM权重。
C0首窗forward5.19s/backward4.29s含GPU/GRU首次初始化，后续forward约0.064–0.067s、
supervision+backward约0.117–0.123s（另有一次0.492s）。
B1每窗forward约0.065s、proposal0.102–0.138s、supervision+backward0.205–0.500s
（首窗1.329s）。全部逐窗值见[JSON](r3_head_timing.json)。

C0先读取缓存，I/O 0.097–1.433s；B1后运行，缓存已热，I/O约0.024–0.027s。
**这个I/O差异不能解释成B1本身更快**。小窗计时也不是训练后候选分布或全染色体保证。

实际两物种都覆盖了正/参考负连接：拟南芥两窗10/42个正连接，水稻三窗4/21/3个正连接。
对应正连接logit梯度为负、负连接为正、未知为0，两物种均通过。
因此补齐了此前水稻首个单exon窗未测试正连接的局限；不靠替换样本或新增挑例完成。
自由候选先落盘再构建参考监督，所有backward梯度有限，未修改模型参数。

本R3新增GPU累计：GLM两次152s + 本作业35s =187s = **0.051944 GPUh**。
待执行的仍是全染色体caller和冻结首训，不是扩大工程smoke或宣称方法效果。
