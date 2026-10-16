# M26：同范围性能与条件机制诊断

日期：2026-09-26。状态：执行前协议；开发性、探索性研究，不是独立测试或期刊主张。用户已明确授权转向机制与公平系统评测、GitHub 公开发布及 Pro 商议。M25R 固定路线 NO-GO 保持不变，发表目标并未完成。

## 这一次运行要改变什么决定

先回答：R1 的内部改进在同一坐标、参考和评价器下是否仍落后于既有工具？“候选覆盖不足”和“输出非精确”能否被逐条身份核验分开？结果决定后续是研究一个可迁移的机制，还是仅整理实现特有的失败；本次不据结果调整解码器。

## 固定输入与比较范围

- 拟南芥 NC_003074.8（23,459,830 bp；4,151 complete-primary 参考）、水稻 NC_089041.1（29,936,421 bp；2,299），合计 6,450；仅开发染色体。
- 参考来自各自 `data/m1_screen/<species>/reference.gff3`。最长 CDS、转录本 ID 打破平局；排除参考 partial，与 R1 一致。不选择 canonical-only 参考。
- M25R A/B：`outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY/*_predictions.gff3`。
- ANNEVO/Helixer/Tiberius：既有 M12B 完整预测，路径由脚本固定；不新训练、不新推理、不挑表现最好的版本。每物种先限定 seqid 后解析。无覆盖则报缺失，不能记为零分。
- baseline 历史配置是 `configs/M12B-SAMEPANEL-BASELINES.yaml`，历史执行命令是对应 EXTERNAL sbatch；它们引用相同原始 FASTA 路径，但不能仅凭同路径证明历史内容逐字相同。输入版本/预训练暴露的不足明确保留，不宣称参数量、训练预算、上下文或零样本条件相同。
- 预测采用最长 CDS primary，保留预测 partial，不做基于参考的过滤。读取坐标后不自动补 stop codon；分别输出每工具在实际 FASTA 上的末端 stop、ATG、长度模 3 分布，用来发现格式约定差异。若差异使 exact-chain 比较失真，应将该比较标为待对齐，预先说明修正规则后再重评分，不能据哪个版本分高择优。

## 预定输出

1. pooled 与每物种表：exact strand-aware CDS interval/chain F1、TP/P/R、普通 CDS-base F1、R1 原定义 gene-span intergenic FPR、gene-count ratio。复用冻结指标函数；必须逐项复现 B 指标到 1e-12，否则停止解释。
2. 单列 FPR 参考政策敏感性：以全部已注释 gene span（含 partial 和非编码 gene）定义背景；不替代原指标。不能把完整 primary 之外的所有区域称作真正无基因。
3. 每个预测按优先级互斥分类：精确 primary；精确其他完整 isoform；同 intron chain 但端点不同；同链参考跨度重叠；仅反链重叠；无完整 primary CDS-span 重叠。后几类是参考差异，不是已证实假阳性、fusion、pseudogene 或新基因。
4. 由完整参考 ID 与 oracle JSONL 联接，计算 B exact∩reachable、B exact−reachable、reachable−B exact。先核验全集 6,450 与唯一键；若 B exact 超出 oracle，报告定义不一致而不使用“图恢复率”。可达性仍是当前图的 reference-assisted 条件诊断，不能从基线最终 GFF 推断其内部候选图。
5. 输出小型 JSON/TSV 和每链分类账；所有方法、物种和类别均保留，不挑案例。短 ORF、重复区和生物学支持在本次尚未归因。

## 资源、失败与停止

一次 CPU Slurm：private-teodoro-gpu，0 GPU，4 CPU，16 GiB，40 分钟上限，预计数分钟；环境沿用 generanno。独立目录 `outputs/M26-SAME-SCOPE-MECHANISM-R1`，不覆盖 M12B/M25R。先跑三个针对链方向、坐标、FPR 分母和 nonexact 分类的合成测试，不冒充真实结果。

运行、范围或解析失败时先保留失败证据，只允许科学定义不变的明确工程修复；指标不理想不是重跑/调参理由。缺输入、不同 stop-codon 约定或无法确认的历史版本需要如实报告，不能拼表补零。一个条件图的 oracle 结论不能提升为跨模型普遍机制。

本阶段不访问 Setaria，不训练，不扫描阈值/权重，不取消其他作业，不下载 RNA-seq、不使用参考辅助预测。后续独立转录/蛋白证据及跨物种确认面板，需要先完成数据来源、独立性、统一分母及资源的具体方案；Pro 提供建议，不提供授权或实测证据。
