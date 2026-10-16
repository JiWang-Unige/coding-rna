# M26 R6-Q：同发布注释类型输入资格核验

2026-09-26，评分前固定。R5两次失败保留，不再重试；本阶段是reference-only输入准备，不打开预测或产生工具指标。第四轮Pro商议正在审阅4c5f90a；本准备不预设其对后续论文主路径的意见。

## 唯一问题与决策

同一个TAIR数据发布中明确提供的gene_model_type，能否无歧义补齐Chr3 GFF缺失的locus_type，并与已有规范类型一致？一致且结构可解析才产出合格输入；冲突/未知逐项保留，不默认coding，不根据模型输出处理。

## 固定来源与规则

- https://zenodo.org/records/15889110 ：TAIR_Data_20240630，同记录GFF日期20240701；原始压缩GFF保留。
- 同记录 Araport11_functional_descriptions_20240630.txt.gz（2,513,392字节），表头name、gene_model_type。下载链接：https://zenodo.org/records/15889110/files/Araport11_functional_descriptions_20240630.txt.gz?download=1 。读取前两列，不从描述判断类型。
- 范围仅原DEV Chr3 -> NC_003074.8，23,459,830bp。原R5已逐碱基核验TAIR FASTA与项目序列一致。本资格阶段不改序列。
- gene ID由GFF显式ID取得；功能表模型去除末尾isoform后缀后对应locus。同一locus类型须一致；表与已有GFF规范locus_type冲突则不合格。
- 明确GFF类型原样保留，缺失时才使用同发布功能表唯一类型。pseudo/非编码/TE均不改为coding；不从是否有CDS推断类型。
- 全部Chr3原始feature行及坐标、strand、phase、feature名称保留。派生文件仅保留评分所需结构属性ID/Parent/Name/locus_tag/locus_type/partial/start_range/end_range/exception/transl_except/pseudo/pseudogene，并新增统一gene_biotype。自由描述中非规范分号不作为属性；原文件不动。不是丢弃gene或结构行。
- 检查gene唯一、mRNA单一已声明gene Parent、CDS已声明mRNA Parent、坐标范围/链向与CDS phase。失败阻止派生合格GFF。
- Latin-1用于这份已知非UTF8源文件的字节映射；类型/ID为ASCII。描述不进入指标，不将其转换结果称为文本勘误。

## 输出与边界

outputs/M26-ARAPORT-INPUT-R6/ 下：input_qualification.json、gene_type_resolution.tsv，且仅qualified=true时生成araport_chr3.qualified.gff3。状态表示资格分析完成，不等同输入必然合格。单独保存逐gene映射、类型冲突和缺失，Git只公开紧凑计数及代码/协议，不公开原始数据。

一次CPU作业，private-teodoro-gpu，2CPU/2GiB/5分钟、0GPU。预计秒级，使用已验证generanno环境中的Python/pytest。本阶段不自动重试。若输入仍不合格，记录具体缺项，再由阶段商议选择是否继续来源对照；不无限找映射。后续评分必须另行确定使用本输入的合同，不冒充原R5完成。

Setaria封存；0新推理、0训练、0参数搜索。数据质量检查只服务“能否安全比较这个来源”这个决定。
