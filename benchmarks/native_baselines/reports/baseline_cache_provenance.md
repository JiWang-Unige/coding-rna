# M28 缓存基线来源与迁移解释边界

2026-10-01。本记录说明当前M28评分表中的ANNEVO、Helixer和Tiberius预测来自哪些既有资产，以及现有记录能支持什么版本和数据暴露结论。当前比较可以继续使用这些固定预测；它们不是本轮重新运行的最新版本基准，也不能整体称为各方法均物种未见的zero-shot测试。本次只读取配置、运行收据及已允许评价的来源字段，没有新基线推理、权重下载或数据放行。

## 评分链路与实际缓存

当前共同评价器[common_ruler.py](../../../pipeline/m28/scripts/common_ruler.py)读取下列M12B缓存，在固定拟南芥NC_003074.8和水稻NC_089041.1上重评分；主分母为5437和2291条CDS-assessable primary链。原M12B的CDS-span指标不是这里的exact-chain F1，不能把历史metrics.json中的数值直接替换当前共同评价结果。

| 方法 | 拟南芥与水稻缓存文件 | 读取到的配置或命令 | 仍未由收据证明的内容 |
|---|---|---|---|
| ANNEVO | outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/predictions/{species}.gff | Magnoliopsida，ANNEVO_Magnoliopsida.pt，batch32；原驱动为annotation.py。 | env.txt未给出执行时源码commit和权重完整版本身份，不据当前checkout反推旧执行版本或预训练隔离。 |
| Helixer | outputs/M12B-SAMEPANEL-BASELINES-HELIXER/predictions/{species}.gff3 | land_plant_v0.3_m_0100.h5，64152bp subsequence，batch16；容器路径为helixer-docker_latest.sif。 | latest文件名不是不可变发布版本；当前所读env.txt未给出当时容器/HelixerPost的明确版本，不能称最新工具成绩。 |
| Tiberius | outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/predictions/{species}.gtf | 原驱动用--singularity --model_cfg angiosperms，batch4。缓存预测日志记录chunk200070、Softmasking=False、postprocessing=True。 | M12B YAML列出的tiberius_2.0.5.sif不能单独证明实际用过该容器；入口按已安装tiberius包版本选择容器。实际容器及其内部解析的配置/权重版本尚未由这次收据重建。 |

三份env.txt均记录2026-06-17；分别记载ANNEVO job8982049/task0、Helixer job8982048/task2、Tiberius job8982069/task1。这里保留原字段，不把job与task拼成未经确认的Slurm数组身份。本次没有用旧code-review gate或validate_goal的状态作为现在的研究授权，也没有恢复已退役脚本。

## Tiberius配置提示的物种暴露

本次读取的host checkout为8c49fd029909eaa6098160247363e68e5ecde698，model_cfg/angiosperms.yaml无本地改动；同一文件已在[官方固定提交](https://raw.githubusercontent.com/Gaius-Augustus/Tiberius/8c49fd029909eaa6098160247363e68e5ecde698/model_cfg/angiosperms.yaml)核对。它列出tiberius_version 1.1.6、配置日期2026-04-15、Softmasking=False，训练材料包含OsativaKitaakev3.1及Sviridisv4.1。

这是配置作者对该模型的训练材料说明，不是本项目独立审计过的完整训练数据清单。入口传给容器的是原--model_cfg angiosperms简写，因此容器内可解析其自身配置；不能仅凭host文件就把该列表当作旧缓存运行所用权重的精确暴露记录。也不能把配置中的模型兼容版本、host包版本、容器版本当作同一个字段。

可据此作出的研究判断是：现有材料不足以支持“当前水稻比较对所有基线都属物种未见”的说法。若使用这个已记录的angiosperms模型作未来迁移比较，水稻不是该模型作者声明的训练物种未见对象；具体Kitaake材料与本项目参考是否共享序列或同源簇仍未核实，不能直接宣称样本泄漏。训练列表出现Sviridis只是模型元数据读取，不改变本项目Setaria继续封存的边界。

## 对后续比较的影响

B2事前合同中的历史缓存表照常保留，不换分母、不替换预测或根据新成绩挑选较弱基线。缺少完整版本收据不使已经固定预测的同口径重评分失效，但会限制可发表的复现和zero-shot表述。

如完整方法达到进入新独立面板的条件，冻结协议应分别记录每个方法的实际源码/容器/权重身份、容器内生效配置、推理操作点和训练暴露资料，并使用适用于目标类群的版本。模型元数据不完备时明确称暴露未知，不能写作未见。工具无需与M28相同训练预算，但方法内部训练因果消融需要额外曝光匹配，两种公平性不能混同。

是否以更高完整链性能还是更广适用范围作为主贡献，要由独立结果决定。若选择未见物种迁移，应按各方法训练和基础模型暴露分别标注资格；不同基线见过不同物种时，不能强行把整张表统称干净zero-shot。这里只提出解释边界，不选新物种、运行新面板或更改B2。

## 本轮核对范围

本轮读取了M12B配置、三份env.txt、Tiberius既有拟南芥预测日志、旧sbatch的只读Git对象、当前Tiberius入口及官方固定配置。没有执行旧availability写入脚本或任何旧Auto Research流程；没有重新计算大权重摘要。Tiberius实际旧容器版本和ANNEVO/Helixer旧二进制精确身份仍是已标注缺口，不能由文件存在、配置名或今天的环境版本补成已验证事实。
