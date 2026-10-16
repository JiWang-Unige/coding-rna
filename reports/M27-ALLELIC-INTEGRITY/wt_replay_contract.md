# M27 WT 原生／缓存重放合同

日期：2026-09-26。状态：FROZEN_BEFORE_FIRST_WT。本阶段只有原序列，不进行任何 mutant 推理或训练。

## 决策与 Pro 第九轮

Pro 已读固定提交 67f5c9a 的准备报告、合同与 JSON（ACCESS_VERIFIED，未重跑），支持保留全部 21 位点、不重新筛 MANE/readthrough，并撤下没有依据的共同 WT 命中数门槛。我们接受先做原生响应、后做分量归因的收窄。

首要 estimand 留给后续单变异阶段：在每工具自己的 WT 精确恢复集合 W_m 内，PTC 与 synonymous 配对替换的绕行比例差；共同集合只用于同位点跨工具比较。两种替换的碱基不同，因此不是纯粹的生物学 stop 效应。算法绕行也不自动是生物学错误。ACE/SGRF 的旧概念不能成为新意。

本次具体决策只有：固定工具能否在整条 chr22 上精确缓存重放，哪些已登记位点具有唯一 WT 精确链，以及实测成本是否允许后续逐变异。WT 成绩只决定可评价条件分母，不用于选位点/权重/阈值。

## 不变的输入与模型

- 输入：准备作业生成的 data/m27_grch38_chr22/chr22.unmasked.fa，50,818,468 bp；不裁剪上下文、不合并编辑、不输入外部 RNA/蛋白/掩蔽轨道。
- 参考：outputs/M27-ALLELIC-PREP-R1/preparation.json 的全部 21 个已登记位点；GENCODE 49 CDS 链已显式合并 terminal stop。
- ANNEVO：源码 37bdd9aa62ddf24fa55941fb827061f7ed49ce53（v2.3.2），既有 saved_model/ANNEVO_Mammalia.pt。原生 FP32 神经计算及发布代码写出的 float16 HDF5，不另行量化。发布窗口 102400+两侧12800、无 overlap_pred；batch_size=2、num_workers=2；解码原 min_intron_length=20、min_prot_length=100、8 CPU。
- Tiberius：现有 2.0.5 SIF，TF 2.17.0+nv25.2 / bricks2marble 0.0.6 / hidten 0.0.1。使用发布 mammalia_nosofttmasking_v2 配置对应权重 URL https://bioinf.uni-greifswald.de/bioinf/tiberius/models/tiberius_nosm_weights_v2.tar.gz 。已下载 29,805,917 bytes，源 Last-Modified=2025-10-24 17:16:14 GMT；此为固定发布权重，不是“已证明未见 human”的模型。运行前解压至专用 refs/weights/tiberius-m27-nosm-v2，不隐式切换/更新权重。
- 已读权重配置：inp_size=5、output_size=15、clamsa=false、hmm=false、units=372、numb_conv=3、numb_lstm=2、pool_size=9。使用原生默认 HMM（不是借口宣称端到端训练机理）；seq_len=400050、batch_size=2、no_softmasking，其他 native 参数保持默认。
- 保留 native 自适应分块、边界 repredict、合并规则和 DNA 后过滤，不为了重放固定窗口。不修改两项目第三方源码。

## 执行与一致性

串行执行 ANNEVO 原生两步 prediction.py→decoding.py，再以同 HDF5、同 DNA 重放 decoding.py；不再运行第二遍神经前向。检查 HDF5 双方向均为 chr22×15，保留原始文件。发布解码器会吞掉 worker 异常，因此日志出现 Process failed 就停止，不把部分输出当成功。

随后 Tiberius 一次 native WT，记录每个初始/repredict 的 lstm_prediction 请求、原始单热输入与原始精度分数；只旁路保存，不改返回值。再在新的进程中通过同一 native 控制器重放：逐次元数据与实际神经输入完全相等才可取缓存；缺失、额外或未消费请求均停止，绝不偷偷计算新神经分数。重放 neural_forward_calls 必须为 0。

记录 Tiberius 两个 native 后过滤（最短 coding 长度 200、internal-stop 检查）各自前后全部链；顺序和链均须重现。重推理拼接对象起点不是基因组坐标，本阶段只验证同 DNA 请求 identity，不据此实现 mutant 定位。参见 wt_implementation_notes.md。

全 chr22 比较 transcript CDS 链、strand、phase 的多重集合，忽略生成 ID 和日志顺序；重复链数量也须相等。不是仅比较21个位点，不是宽松 IoU。两发布工具自身 CDS 包含终止密码子，不能再盲目加3 nt。Tiberius过滤前后链与完整请求账也须相等。仅通过此步才报告 replay_qualified=true。

## WT 分母和对应规则

以 chr22、strand、完整终止密码子包含的 CDS 坐标链及连续 phase 精确对应到登记参考。只有唯一一个 exact prediction 的位点进入 W_m；零个记 WT_not_exact，多个同链预测记 ambiguous_duplicate_exact，不任选一个。同名基因不当身份依据，也不按最近邻/最大重叠强配。全部21位点出账，旁边基因或 readthrough 输出不得误记作目标基因的变异响应。

本阶段不做 genome-wide gene-accuracy 排名，也不把少量 WT 命中解释为模型泛化。重叠、其他 RNA isoform、readthrough 等留作后续解释层；若后续 mutant 无法唯一对应，标注不可判定，不能据此补选新位点。

预先固定后续分支：W_m 为空表示该工具在此面板无条件分析信息；共同集合为空不否决各工具分析。“无效应”只能在完成配对干预后判断，零净差有正反抵消时必须列出。可重放的 PTC 特异绕行才进入机制确认；两工具配对差都为正才考虑跨工具扩展，仅一个工具出现就限制主张。这些不是显著性、群体发生率或发表门槛。

## 资源和停止

作业 M27-WT-REPLAY-R1：1×RTX3090 24GiB，8 CPU，32 GiB RAM，2小时 walltime，private-teodoro-gpu；占未来先导总上限8 GPUh中的最多2 GPUh。输出与日志均在独立 outputs/M27-WT-REPLAY-R1。数据准备与缓存全体新增存储上限50 GiB；Tiberius缓存单独限制20 GiB，按原始float精度保存。native/重放串行，所有原始输入和缓存保留。

预期十几分钟到约一小时，真实峰值未知，以本次计时、MaxRSS和磁盘统计为准。超时、资源超限、模型装载失败、worker错误或重放不等价即停止，不自动调 batch/窗口/精度、换模型或补跑。合成 cache/input/phase/重复对应测试必须先通过。此合同不运行一个 mutant，也不释放 Setaria。下一阶段要根据实际成本重新核算43个完整输入/工具，并完成 Pro 商议后才执行。
