# M27 WT 重放前的实际实现备注

2026-09-26，参考准备提交 67f5c9a 之后的代码核实。第九轮 Pro 商议已完成；采纳决定与 WT 执行边界见 wt_replay_contract.md。本文件不是模型结果或完整干预协议。

## Tiberius 的坐标不能只看 Sequence.start

实际容器 Tiberius 2.0.5 / bricks2marble 0.0.6。读取的文件：

- /usr/local/lib/python3.12/dist-packages/bricks2marble/struct/fasta.py：原始 Sequence 包含 name/start/end；copy() 重建 Sequence，不能假设所有额外 evidence 自动保留。
- /usr/local/lib/python3.12/dist-packages/bricks2marble/tools/annotate.py，_annotate：初次预测后由 _find_mismatches 决定边界重推理；每条序列将相关相邻块的两侧半窗口拼接、flatten，重建一个 name 相同但起点默认 0 的 Sequence，再 resample。原始窗口索引保留在控制器局部变量 repred_index / repred_sequence / repred_strand 中，不是新 Fasta 对象的基因组起点。

后果：不能直接用重推理对象的 start/end 给变异定位。也不能把拼接多个不相邻窗口的序列当作连续染色体。这不改变 native 预测，而是限制我们怎样正确编写干预。

最小 WT 验证不需要解决突变坐标：保持 native 自适应控制器，只在每次 lstm_prediction 边界记录原生输入（含方向、shape、实际单热编码）和未改动分数；缓存重放必须逐次输入一致且不允许缓存缺项时悄悄调用神经模型。然后比较完整 chr22 的规范化输出链以及每次请求。WT 通过仅证明同条件缓存重放，不证明任意 mutant 交叉条件可用同一缓存。

后续若执行突变交叉，必须选定科学 estimand，再决定记录控制器的真实窗口映射或使用明确标注的固定窗口局部 HMM 诊断。不能为了省缓存而静默关闭 repredict，也不能把 native 自适应、后过滤、剪接语法等所有贡献统称纯 ORF grammar。

## 输入与资源的已知事实

已核实的源码 mammalia_nosofttmasking_v2.yaml 指定无 softmask、无 ClaMSA、默认序列长度 400050；这些是发布配置而不是待扫参数。主程序还指定 native 后过滤最短 coding 长度 200 及删除 internal-stop 链。

在已核实的 checkout 默认 model_weights 位置和容器默认 model_weights 位置，没有对应目录；项目顶层也未见 tiberius_nosm_weights_v2。这里只报告检查过的位置，不能声称整个账户没有权重。后续需明确固定模型的下载/装载路径及版本，不应调用隐式下载并误以为命中了历史权重。

准备输入 229 MiB、ANNEVO Mammalia 既有权重 56 MiB；未来单份双链15类float16全 chr22 分数约2.84 GiB，Tiberius若保留float32同形分数约5.68 GiB，实际请求量还受padding/repredict影响。新模型执行前先固定 WT 缓存上限、时间上限、精度/输入策略和失败规则。当前没有运行任何模型，不将这些估计写成实测资源。
