# M27 原生确认结果：两案例完整复现，尚未完成分量归因

2026-09-27。状态：COMPLETED / 4_OF_4_TARGET_AND_WHOLE_CONTEXT_REPRODUCED。

## 核心结果与解释边界

按首mutant前规则选定的ANNEVO/POTEH与Tiberius/CRYBA4，各自syn/PTC均完成原生整chr22重算。新旧目标完整CDS链、strand、phase、候选多重性/竞争身份与分类全部一致；四次同allele重复的全染色体链多重集合也完全一致，没有新增或丢失链。Tiberius两次末端过滤前后的完整链及真实请求元数据也一致。

这确认两个预定案例具有计算重现性，并建立可供后续检验的神经缓存；**不是新的独立阳性，不扩大12/10或共同9的统计分母，不证明真实剪接、生物学损伤被掩盖或选择性风险诊断。**没有替换神经分数，没有实施分量干预，本mutant的同DNA缓存重放仍未验证。原44例分类、未知和界限均保持不动。

## 事前协议与执行

[确认合同](confirmation_contract.md)、代码与19项CPU测试在首次确认模型运行前公开于0a8436d5f12b936a430623353a7facab4a8d1ea3，原配对结果固定43e366a。Pro第12轮建议有限四次确认与旁路记录合并，未访问原始产物或独立重跑；没有因为结果选新案例或改分类。

CPU准备13233925 COMPLETED 0:0，12秒，0GPU，19测试通过（4.79秒）。GPU命令：sbatch --parsable sbatch/M27-NATIVE-CONFIRM.sbatch；job13233933及batch/extern均COMPLETED 0:0，gpu034，1RTX3090、8CPU、32GiB、45分钟硬上限，无重试。实际19分39秒（1179秒，0.3275分配GPUh），MaxRSS14,817,172KiB（14.13GiB），主stderr 0 bytes。4个case execution均完成，6个子命令returncode均0。

## 逐例复现

| 工具 / 固定对象 | allele | 目标响应 | 新旧全chr CDS链 | 差异 |
|---|---|---|---|---|
| ANNEVO / POTEH，ordinal1 | syn | WT_CHAIN_RETAINED，同一链 | 530 / 530 | 0新增、0丢失 |
| ANNEVO / POTEH，ordinal1 | PTC | BYPASS，同一174nt内部CDS跳过 | 530 / 530 | 0新增、0丢失 |
| Tiberius / CRYBA4，ordinal9 | syn | WT_CHAIN_RETAINED，同一链 | 575 / 575 | 0新增、0丢失 |
| Tiberius / CRYBA4，ordinal9 | PTC | BYPASS，同一15nt新内含子 | 575 / 575 | 0新增、0丢失 |

比较单位是同一allele的新旧运行，不是syn与PTC互相相等；忽略自动生成的预测ID，但不忽略结构、多重性、相交候选或竞争锚。完整目标键与逐项比较见原始紧凑 [confirmation_result.json](confirmation_result.json)。

Tiberius两allele都完整复现过滤链706→585→575，目标在每步过滤前后同一身份/类别；每例实际4个神经请求（2 initial、2 repredict）与4次真实前向，形状/精度/上下文元数据同旧运行。4是观察值，不是硬编码断言；本次总8次真实前向，不是缓存代替推理。

CRYBA4的[26625574,26625589)内含子序列为GCCTAACCCGCCGAG，GC–AG边界，完整编辑TAA位于其中；这属于工具允许的模式，但15nt结构没有独立转录验证。原生输出照录，不事后加最小长度筛选。POTEH保留双端锚而跳过[15695644,15695818)的174nt原CDS。合法ORF和可重现不等于生物学正确。

## 缓存与资源账

每份ANNEVO原生HDF5为3,049,112,560 bytes，两个strand均50,818,468×15、float16；两份均保留。每份Tiberius缓存为8,641,081,024 bytes（约8.05GiB），含真实请求的输入/输出，float32；各低于20GiB硬上限，写入前受M27总量剩余约束。记录原样返回原生数组，没有改变网络输出或精度。

确认目录outputs/M27-NATIVE-CONFIRM-R1为23,388,900,090 bytes。按原六个目录加确认目录核算，M27新增数据合计36,819,615,045 bytes（34.29GiB），低于50GiB。没有删除原始产物。公开Git只含代码、协议和紧凑结果；HDF5/神经数组、FASTA、完整GFF/GTF及日志留在Baobab。

分配GPU账：WT949秒＋原配对13,400秒＋本次1,179秒＝15,528秒，累计4.313333GPUh，8GPUh内剩余3.686667GPUh。0训练、0Setaria访问。

## 下一决策

两个目标及全局上下文均稳定，没有本轮发现的全局复现障碍；但不能跳过同DNA缓存等价对照。只读 [实现边界](component_feasibility.md) 显示ANNEVO的DNA条件转换、短内含子重解码和分数过滤，以及Tiberius的核苷酸发射、动态拼接重推理和末端过滤必须区分。

下一轮Pro商议只决定一个可执行的有限检验：保持神经证据，定位实际必要的显式完整性步骤，并检查它是否逆转对原结构与绕行结构的神经发射得分偏好。若无法可靠定位或只再确认已知约束行为，收束为受限实现归因，不扩真实变异/物种/干预组合，也不宣告Nature Communications目标完成。
