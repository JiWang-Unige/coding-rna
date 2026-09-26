# M27 准备结果：21 个参考选择的配对位点，尚无模型实验结果

日期：2026-09-26。结论：**参考材料与实际容器接口已取得；可进入有预算的 WT 重放设计，但目前不能声称交叉解码已验证或等位完整性假说成立。**

本阶段按[准备合同](preparation_contract.md)执行，不延续已关闭的 M26 R8，不重开 Setaria。完整紧凑清单在 [preparation.json](preparation.json)；原始 FASTA/GTF（229 MiB）及完整运行日志仅保留在 Baobab。

## 得到了什么

固定 UCSC hg38 chr22（50,818,468 bp）及 GENCODE 49 comprehensive CHR（文件头 Ensembl 115，2025-07-08），不是最新注释。原 FASTA 保留，另生成大写输入，供未来 sequence-only / no-softmask 模型使用。

chr22 共 1,747 个全部类型 gene、11,615 个 transcript。5,267 个 transcript_type=protein_coding 中 4,710 个通过完整 CDS 资格检查；546 个部分 CDS/selenocysteine，另 11 个由于 stop 几何、非规范 ORF/歧义序列或起始特征不合格。没有把后者改标为完整结构。

下表分母为具有 protein_coding 转录本的 434 个基因，不是 chr22 全部基因：

| 参考资格 | 基因数 |
|---|---:|
| 多个不同的可评估 protein-coding CDS 链 | 349 |
| 还存在不可评估的 protein-coding transcript | 15 |
| 没有可评估 protein-coding CDS 链 | 7 |
| 单一合格链，但没有满足位置/边界/配对条件的内部 Tyr 密码子 | 42 |
| 合格配对位点 | **21** |

21 个位点全部保留，没有修改规则凑到 24 个，没有运行或读取任何预测。每个配对都是同一个基因组碱基的 synonymous Tyr 与 TAA 替换；每次干预应各自作用于完整 chr22，不把多个位点合并成一个突变染色体。正链 10 个、负链 11 个。基因、转录本、参考 CDS 链、1-based 碱基坐标及正向基因组 REF/ALT 见 JSON。

这是受限的 Tyr-to-stop、单一可评估 protein-coding 链面板，**不是**随机基因组位点样本或所有 RNA isoform 唯一的集合。面板也没有按 MANE/CCDS 或预测成功率筛选；包含 readthrough 名称/未命名基因等需保持可见的参考复杂性。未改变目标链 GT/AG 二核苷酸并不证明真实剪接不受影响。它只能先测试算法响应，不证明生物学错误或致病性。

## 接口可行性与尚未解决的问题

ANNEVO 固定源码 37bdd9aa62ddf24fa55941fb827061f7ed49ce53（v2.3.2）已有 Mammalia 权重；原生 prediction.py 与 decoding.py 分离，HDF5 提供两方向 15 类神经概率。

实际 Tiberius 容器核实为：

- Tiberius 2.0.5；Python 3.12.3；TensorFlow 2.17.0+nv25.2；bricks2marble 0.0.6；hidten 0.0.1。
- lstm_prediction(self, inp_chunks, clamsa_inp=None, batch_size=None)。
- hmm_prediction(self, nuc_seq, lstm_predictions, batch_size=None)。
- 原生 predict_function / repredict_function 均分别调用上述两个路径。

因此接口层面允许独立提供分数与序列，但**端到端分离仍待验证**。Tiberius 会根据边界进行重推理，最后还按实际 DNA 删除含内部 stop 的链。必须控制这些路径和缓存坐标；只替换一次 HMM 输入不能自动获得“神经 vs ORF grammar”的完整归因。重推理窗口随干预变化时，固定 WT 缓存未必覆盖新请求。

本次没有加载任何模型，没有下载新权重，0 推理、0 训练、0 GPU。容器导入出现 cuFFT/cuDNN/cuBLAS 重复注册信息，但导入和接口检查成功；不能据此推断 GPU 推理兼容性已通过。

## 执行及失败记录

| 作业 | 状态 | 时间 | MaxRSS | 结果 |
|---|---|---:|---:|---|
| 13227362 / M27-ALLELIC-PREP-R1 | FAILED 2:0 | 50 秒 | 450,488 KiB | 5 项测试通过、参考准备和 JSON 完成；容器 cwd 软链接不可达，接口步骤未执行 |
| 13227380 / M27-ALLELIC-PREP-R2 | COMPLETED 0:0 | 68 秒 | 1,228,000 KiB | 仅重跑容器接口读取，成功 |

两作业均 private-teodoro-gpu、2 CPU、8 GiB、0 GPU，合计约 0.066 分配 CPU 小时。第一次日志、失败状态、已产出参考文件保持不变。唯一工程修复将实际 BeeGFS 目录绑定为容器 /work；未重复下载/筛选、未改变科学条件。原准备脚本同步修正路径，单独修复脚本保留重试过程。

源码：[prepare_reference.py](../../scripts/experiments/M27-ALLELIC-INTEGRITY/prepare_reference.py)、[runtime_probe.py](../../scripts/experiments/M27-ALLELIC-INTEGRITY/runtime_probe.py)、[测试](../../tests/test_m27_prepare.py)、[准备作业](../../sbatch/M27-ALLELIC-PREP.sbatch)、[接口重试](../../sbatch/M27-ALLELIC-PROBE-REPAIR.sbatch)。

## 下一决策

不采用缺乏依据的“>=20/24 共同 WT 正确”门槛，也不在看见预测后改门槛。先与 Pro 商议 21 位点的可判断范围，预先确定各工具及共同 WT 正确的条件分母、效果分类、无信息与无效应的区分，以及 native/分离 WT 重放一致性要求。

完整单变异方案有 1 WT + 42 mutant 染色体/工具。单份 ANNEVO 双链 float16 神经分数约 3.05 GB（2.84 GiB），不能把 43 份全量分数永久缓存并仍称满足 50 GiB 上限。需要逐变异处理并保存足以重放的差分缓存，或在第一次模型运行前明确更小的初始可执行范围；不能用联合编辑或裁剪上下文偷换问题。

文献新意边界见合同：ACE/SGRF 已有相近旧概念；即便先导出现绕行，也不等于 Nature Communications 级贡献。未来只有可重复分量干预、独立真实等位证据和实用纠正原则成立后，才考虑扩大。发表目标仍未完成。
