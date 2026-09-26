# M26 R7：限定protein-coding对象的同序列来源配对

2026-09-26，评分前固定。本阶段承接R6已完成的输入分类，不重写R5/R6失败，不仅换目录重投旧合同。沿用R5科学问题，但明确新增同发布gene_model_type、gene-like父节点支持和miRNA类型差异分层；不更换来源、不依预测选择规则。

## 对象资格

输入仍为Zenodo15889110的Araport11.20240701 GFF、同发布20240630功能表与TAIR10 FASTA；对照原RefSeq GCF_000001735.4。仅Chr3/NC_003074.8，23,459,830bp，运行再次逐碱基核验同序列。无第三个参考。

显式gene、pseudogene、transposable_element_gene均可作为父节点，原feature名与坐标/strand/phase保留。将879个旧“missing gene”报告项逐一连接，报告实际父节点feature分布和完整mRNA父节点账；不能仅凭例子认定全部TE。TE与pseudo从不改成普通protein_coding gene。

gene类型：明确GFF与同发布标准model唯一类型一致则保留；缺失规范字段才以明确功能类型补齐。保留19 uORF对象，不给宿主传播类型。50个mirna vs miRNA_primary_transcript差异单列，保留原GFF类型，不宣称两个对象名称等价；它们均明确非coding，所以只对protein-coding资格判定不构成不确定。若出现涉及coding身份的冲突、其他未知类型或结构ID/Parent/坐标错误，停止评分并保留问题清单。

派生GFF保留全部Chr3 feature行及所有结构属性，描述原文留在原始文件。mRNA唯一gene-like Parent、CDS唯一mRNA Parent、父子坐标范围与strand、CDS phase取值必须成立；并按转录顺序核验链内phase一致性（第一段可非0，两来源不必phase相同），问题单列后停止本次评分，不静默删除。未知coding对象不静默丢弃。此资格不等于生物学真值。

## 固定比较

五套预测为M26原A/B、2026-06保存ANNEVO/Helixer/Tiberius；不新推理、不调参、不修改输出。两来源分别：
- parent-filter：原gene/mRNA partial排除、最长CDS再ID primary；原RefSeq4,151和5工具该染色体chain指标及span FPR必须精确重放。
- CDS-assessable：CDS本身所有行无partial/range，且通过R4固定sequence compatibility；每gene最长CDS再ID。ORF兼容不是表达/功能证明；保留不兼容计数和旧集合增减。
- strand-aware exact CDS链包含terminal stop，坐标不依参考成绩加减；phase不混入chain身份。
- 报告参考数、预测数、TP/P/R/F1及原政策span背景分子分母，另报全部coding gene背景CDS/span覆盖。非coding区域不会因接入gene-like父节点就自动被屏蔽，FPR不是生物学假阳性真值。

历史6,450及NO-GO保留。旧模型版本/历史FASTA身份局限不因重新评分消失，不称最新完美公平SOTA。

## 决策性分解

1. 同locus、同strand+完整CDS坐标的共享链：两政策资格对应、受父级过滤影响的链/基因数及工具恢复数。只报对象计数，不对共同子集硬造总体precision。
2. 语义对齐后逐locus primary：
   - 相同primary；
   - primary不同但双方选中的链都存在于两来源可评价isoform集合：仅primary选择差异；
   - primary不同且不能由这种选择解释：primary结构争议；
   - 仅一来源CDS可评价、双方均不可评价，单列，不冒充一对冲突结构。
3. 即使primary相同，也计数任意CDS-assessable isoform集合差异，防止primary表掩盖其他结构分歧。
4. 同时比较raw CDS集合，区分CDS自身资格差异与真正缺失的结构；未对应locus不强配。对共享locus中真正来源独有可评价CDS，输出逐工具primary/任意兼容isoform精确匹配判定翻转数与分母。
5. 每组输出整数分母、每工具匹配各来源primary/任意兼容isoform的locus数及可追溯逐locus账。类型未知对coding比较若为0，其不确定界限为0；若非0则本合同停止评分，不报告确定排名。

## 终点

若对齐后coding链集合/工具比较稳定，关闭本地来源敏感性作为新论文主张的路线，不找第三来源救故事。若仍有实质结构争议，先量化它能否改变比较，随后才考虑最多100位点、模型身份盲化的独立证据验证。来源有共同构建依赖，不称独立复现，不因一方注释/模型同意就判真。

一个Slurm CPU批次，private-teodoro-gpu，2CPU/8GiB/15分钟、0GPU，预计数分钟。outputs/M26-REFERENCE-PAIR-R7/；限定只允许一次不改变对象或政策的明确工程修复，保留独立失败输出，禁止根据成绩改规则。输入资格不合格或原指标不能复现则停止科学解释。不得修改原始输入、旧输出、checkpoint或Setaria。

第六轮Pro针对R6结果和本合同草案商议后确定方向；其意见不是独立重跑或实验授权。执行授权来自用户的逐阶段自主推进请求。
