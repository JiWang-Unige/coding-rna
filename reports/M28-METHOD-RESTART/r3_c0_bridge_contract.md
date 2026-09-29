# M28 R3 C0 原生接口 smoke

2026-09-30，工程接口测试，不是新模型训练/性能实验。

R2内部15类采用沿CDS递增的phase0/1/2；ANNEVO固定依赖37bdd9aa62ddf24fa55941fb827061f7ed49ce53的HMM.py实际把NN输入的各三相位组解释为0/2/1。define_columns状态分组的顺序不能代替viterbi_decoding中的实际输入映射。因此在C0入口显式交换四个三相位组的1/2列。R2标签不改，原文档误称“follows ANNEVO emission mapping”的注释改正；此前没有向原生解码器传入这些标签或报告新模型性能。

范围：合成GT/GC-AG×三种剪接相位、纯背景、同窗多个gene、正反坐标和尾窗归属；读取已存在原生HMM依赖，不下载/复制/发布其源码，不调用其神经网络。该依赖保留自己的非商业科研许可证。测试最小intron内部状态参数1（不是声称真实intron总长1bp），不宣称是正式首训解码参数。

资源：private-teodoro-gpu，2CPU/8GiB/10min，0GPU；只写新的outputs/M28-C0-BRIDGE-R3。失败则修具体接口，不启动首训。它不替代B1自由推理、实际GLM对齐或吞吐测量。
