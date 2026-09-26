# 冻结epoch-1跨block组装可达性：结果与下一范围

2026-09-08。已完成用户批准的单次CPU诊断；未启动新推理或训练。

## 结论

保留原CDS run、允许跨block连接的候选图尚未被可达性上界排除。6450条冻结参考中3585条存在满足canonical ORF的完整路径，最乐观unique-chain F1上界为2×3585/(6450+3585)=0.7144993，高于必要门槛0.55。这个数是假定只输出全部可达真链、零假阳性的上界，不是模型成绩，不同时证明其他门槛可达。

1634条原跨block碎裂参考中685条存在合法路径，947条缺少至少一个边界候选，2条边界各自存在但不能由同一原run组成相应CDS片段。由此判断：值得提出一次无参考的组装实验，但跨block连接不能修复全部碎裂，也不能据此宣布GO或Setaria放行。

## 执行与检查

- Job 12520954，private-teodoro-gpu，2 CPU、16 GiB、30分钟上限，无GPU；COMPLETED 0:0，耗时2分38秒，MaxRSS 2487308 KiB。
- 提交命令：sbatch --parsable sbatch/M25R-E1-BLOCK-ASSEMBLY.sbatch。
- 使用既有generanno环境，未安装依赖。无checkpoint加载、模型forward、训练、阈值搜索、Setaria文件访问或Git提交。
- 6450条逐参考记录完整；1389条既有精确输出全部通过新图可达性检查；原碎裂类别1634条、旧motif-reachable 3616条均核对一致。三个针对有序路径选择的断言在作业内通过。没有失败或重投。

## 候选图的准确范围

事件类型仍为原TRANSITIONS；每个候选CDS的左右端点必须由同一个原始CDS run的两端在±6 bp内支持。开始/内部/结束CDS分别使用start/acceptor与donor/stop事件；片段必须正长度，连接遵循原至少4 bp的GT-AG间隔及递增run顺序。可跨block和跳过run；不设置参考驱动的距离阈值。phase和边界分数过滤作乐观放宽，完整路径另检查原canonical ORF。

参考仅查询路径是否存在，没有生成预测GFF或进行模型选择。因允许重组选取原输出，此图上界不沿用phase-only的固定709个nonexact输出数。

这不是任意组装器的上界：将多个原run合并为一个CDS、在原run内部切分CDS、改变事件类型或候选位置均不在定义内。不能用本报告否决这些未测试机制。

## 全量分解

| 类别 | 全开发集 | 原跨block碎裂组 |
|---|---:|---:|
| 缺少边界候选 | 2834 | 947 |
| 边界存在但同run片段不成立 | 30 | 2 |
| 片段存在但整链顺序/间隔不兼容 | 0 | 0 |
| 完整路径存在，未加ORF过滤 | 3586 | 685 |
| 完整路径存在且canonical ORF合格 | 3585 | 685 |

3586与3585之差为一个非canonical参考，不能合并为同一计数。单边界可达总数仍为3616；新图约束额外排除30条，ORF再排除1条。

合法路径中，拟南芥2506条（水稻1079条）；跨block见证路径分别268条、417条。词典序最小run-ID见证路径最多跨4/12个block，最多跳过2/3个run；可达参考最大跨度7821/18800 bp。这些是已保存见证路径的描述，不是最少block数，也不是全候选图歧义度或待采用的限制参数。全量85条见证路径跳过run，其中66条来自碎裂组。

## 下一步：明确的无参考机制提案（未执行）

建议只提出一个固定epoch-1的候选CDS有向无环图解码器，不新增head训练：节点为本次定义的run支持CDS候选，边为方向一致且可维持完整ORF的donor→acceptor连接，允许跨block。使用原模型区域与边界分数进行路径排序；phase硬拒绝不恢复，以累计CDS长度构建合法phase。选择路径时不得读取参考注释。

真正实施前仍须冻结路径评分公式、跨基因融合抑制规则及多路径输出规则；本次oracle不提供这些选择，也不能从本次真链跨度分布直接设距离阈值。必须先写清这些机制，再申请执行，而不是在验证集上边试边调。

所需信息是候选图节点/边的原模型分数。旧R4未保存全量logits，158窗口面板不足以给完整图打分；应先统计所需唯一原窗口覆盖。最坏情况是17384窗口的一次epoch-1前向，全11通道float16约2.19 GiB，按面板速度粗估纯前向10.35小时，实际时间与内存需按实现确定。此为预算边界，不是已批准的资源请求；没有提交GPU作业。

后续必须同时报告原冻结complete-primary指标、计数/FPR、结构有效性及split/fusion错误；其他注释链匹配只作解释性补充。一次固定机制结果若不满足冻结开发要求则停止，不自动换阈值、增训练或访问Setaria。

当前授权任务已完整交付。整体研究目标仍未达到终态：phase-only已排除，但无参考组装尚未测试；不能把工程成功或oracle可达性当作研究成功。

## 产物

- scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/block_assembly_reachability.py
- sbatch/M25R-E1-BLOCK-ASSEMBLY.sbatch
- outputs/M25R-E1-BLOCK-ASSEMBLY-REACHABILITY/references.jsonl、summary.json、STATUS、JOBID
- logs/M25RBLKCPU_12520954.out、logs/M25RBLKCPU_12520954.err

所有原始数据、checkpoint和冻结结果保留；没有新增后台监控。
