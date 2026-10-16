# M28：独立迁移候选的公开元数据登记，尚未冻结测试面板

2026-10-03 14:33 UTC。当前两物种DEV用于反复开发选择，不能代替独立迁移。本轮只核对4个非Setaria候选的公开组装及注释元数据，为后续冻结模型后的一次性评测准备；不下载FASTA／GFF，不读取其标签内容，不运行推理或选择最终物种。候选是否曾参与本项目历史开发、同源隔离及各比较资产的上游暴露仍须分别核实，当前没有任何一项获得“干净zero-shot”资格。

## 实时 NCBI reference 记录

以下为2026-10-03查询时 `reference_only=true` 返回的具体RefSeq身份，不是永久latest别名。长度是NCBI报告的assembly总长度，可能包含细胞器／unplaced序列；protein-coding gene数是发布元数据，不是未来CDS-assessable primary分母。

| 候选 | taxid／样本 | RefSeq／配对GenBank | assembly名与level | 总长度bp | RefSeq注释、日期／protein-coding genes |
| --- | --- | --- | --- | ---: | --- |
| Brachypodium distachyon | 15368／Bd21 | [GCF_000005505.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005505.3/)／GCA_000005505.4 | Brachypodium_distachyon_v3.0，Chromosome | 271,163,419 | Release103，2018-03-16／25,535 |
| Brassica rapa | 3711／Chiifu-401-42 | [GCF_000309985.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000309985.3/)／GCA_000309985.3 | CAAS_Brap_v3.01，Chromosome | 352,829,610 | Updated Release103.20201202，2020-12-07／41,403 |
| Solanum lycopersicum | 4081／Micro-Tom（BioSample元数据） | [GCF_036512215.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036512215.1/)／GCA_036512215.2 | SLM_r2.1，Chromosome | 832,773,323 | GCF_036512215.1-RS_2024_09，2024-09-20／28,611 |
| Amborella trichopoda | 13333／Santa Cruz75、female | [GCF_054954075.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_054954075.1/)／GCA_054954075.1 | AtrichopodaSantaCruz75HAP2 v2.0，Chromosome | 700,343,102 | GCF_054954075.1-RS_2026_09，2026-09-10／16,734 |

这四项都有RefSeq注释元数据，但不据此保证primary链完整性或标签正确。表中assembly、annotation release与样本必须作为一个组合：不能拿SL4.0／Heinz1706标签配SLM_r2.1／Micro-Tom序列，也不能把旧AMTR1.0标签配新Amborella haplotype。旧版本可另作候选，但不能因结果不理想而替换版本。RefSeq／GenBank pairing不是两套序列完全相同的保证：Brachypodium、番茄及Amborella的返回记录明确列细胞器增补差异。未来核基因任务需要预先登记contig范围，不把表中总长度直接用作核背景分母。

候选用途是覆盖现有两DEV之外的植物谱系，而不是按预期高分挑对象：Brachypodium可用于禾本科迁移候选，Brassica用于十字花科候选，番茄用于另一植物谱系，Amborella用于远缘探索候选。前两项分类见[NCBI Bd21项目](https://www.ncbi.nlm.nih.gov/bioproject/32607)及[NCBI Brassica WGS记录](https://www.ncbi.nlm.nih.gov/Traces/wgs/AENI03)；Amborella的相对系统位置由[原始基因组项目](https://ncbi.nlm.nih.gov/bioproject?LinkName=bioproject_bioproject&from_uid=238126)描述。这里未测量距离、未审计同源重叠，不给候选正式“近／远缘”统计标签。

## 比较资产的资格仍需逐个确定

当前下游头只在本项目允许的拟南芥／水稻TRAIN拟合；这不抹去上游encoder暴露。GENERanno 1.2B preview的物种／assembly／注释监督暴露仍unknown。ANNEVO当前Magnoliopsida权重的作者说明仅给出排除模式物种的陈述和部分例子，不能推断这4项均已排除。Helixer a0080需逐候选核对实际训练列，不能用论文train_val池或不同权重列替代。既有来源边界见[pretrained_exposure.md](../../../paper/supplement/pretrained_exposure.md)。

Tiberius当前固定[angiosperms配置](https://github.com/Gaius-Augustus/Tiberius/blob/871d4979568ae09e56de7eca5d7489046ee833ea/model_cfg/angiosperms.yaml)明确含 `BrapassprapavarrapaVT123v1.1`：Brassica rapa已有作者披露的物种训练暴露，即使本表是不同cultivar／assembly也不能写成species-unseen。该配置未列本表另外三个名字，不等于已审计所有预训练、验证或同源来源。配置target声明为Mesangiospermae；因此Amborella先保留为适用范围待核实的远缘探索对象，不能直接放入公平主表后用外推劣势宣称普遍优越。

动物范围暂不登记正式面板：当前植物下游头未验证动物，ANNEVO Magnoliopsida也不是动物对照资产。更广适用能力需相应固定资产及独立协议，不能由增加物种名称替代实验。

## 真正进入一次性迁移前的必要决定

模型／最后checkpoint／操作点先固定，再核实候选是否未参与本项目开发选择及既有同源隔离要求，并冻结最终assembly+annotation release、核contig范围、primary政策和每个native comparator适用资产。长链、重叠、unsupported及不可达参考按共同评价口径保留，不按新模型可预测性删分母。标签资格与原DEV政策差异必须在看预测性能前说明。

冻结后raw-FASTA caller只读DNA和固定权重，关闭预测后再接入参考；每物种及worst-species报告整链P/R/F1、背景负担和候选可达性，不拼接原论文不同参考定义的分数。本候选登记没有改变当前A9/FULL/A6/A7/A10的gate、数据或停止条件，不自动解封test／Setaria；最终面板资格不足时仍只能作DEV开发比较。

## 来源与执行记录

使用本机NCBI Datasets Skill原脚本，按官方v2 REST路径 `genome/taxon/<taxon>/dataset_report`，只GET元数据。Brachypodium请求于14:29:51.714 UTC成功，其余三项于14:31:05.890–14:31:06.457 UTC成功；各返回1条reference记录、`sources.kind=evidence/supports_claim=true`，表中链接取返回canonical_url。没有保存raw响应或数据包。表内记录链接支持组装／注释事实，不作为模型泛化证据。

Skill示例的 `genome/taxon/assembly_descriptors` 路径返回404，不作为空数据或生物学结论；根据[官方REST说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)改用实际dataset_report路径后成功。运行时原先未直接暴露requests，复用已安装bundled pip vendor依赖完成路由，没有安装包、改base环境或修改Skill脚本。未以检索失败推断缺少基因组或注释。

## 公开暴露与项目历史的定向资格核对（14:56 UTC）

**本轮改变下一步判断的结论：上述四项均不能作为当前Helixer权重的干净未见物种面板。** 这不排除将来预先冻结后进行实际迁移比较，但必须按各方法暴露状态分层，不能称所有方法都在species-unseen条件下比较。

Helixer补充表S29将当前 `land_plant_v0.3_a_0080.h5` 对应到 `seed_009_017`；S23该列明确列短柄草及番茄为训练，S18的其余train_val物种按S23脚注用于验证（每物种随机800 subsequences）。因此白菜和Amborella是作者披露的验证暴露，不是训练，也不是未见。下表是对具体资产列与脚注的解释，不能只从S18的train_val池推断全部训练。[原始补充PDF：S18印刷页66、S23页81、S29页85](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-025-02939-1/MediaObjects/41592_2025_2939_MOESM1_ESM.pdf)

| 当前候选 | Helixer作者数据ID／版本 | 对本a0080权重的暴露 |
| --- | --- | --- |
| Brachypodium distachyon | Bdistachyon／v3.1 | 训练 |
| Brassica rapa | BrapaFPsc／v1.3 | 验证，不是训练 |
| Solanum lycopersicum | Slycopersicum／ITAG3.2 | 训练 |
| Amborella trichopoda | Atrichopoda／v1.0 | 验证，不是训练 |

作者这些版本不同于本候选登记，不据此断言核序列／标签完全相同，也不能据不同cultivar或新assembly改写为species-unseen。Tiberius固定配置已重新读取：白菜训练暴露维持，上述其他三个名字未列，仍不是已审计完整未见证明；Amborella的范围问题维持。ANNEVO四项个别排除状态和GENERanno上游暴露仍unknown，不用Helixer清单替它们补齐。

仅在远程 `docs/`、`reports/`、`configs/` 的md／yaml／yml／json／tsv／csv中定向搜索候选全名、缩写、cultivar、accession及常见词根，并阅读实际命中。除本登记外，发现白菜旧版 `GCF_000309985.2` 在 `configs/m13_distance_generalization_panel.yaml` 的backup_species_if_scaffold_caveat_blocks、`reports/M13_CLOSE_PLANT_FREEZE.json` 及 `docs/10_findings.md:185` 被登记为M13备选；该协议正式测试的是A. lyrata／水稻。这些文本不能证明白菜已实际运行，也不能证明从未运行。其余三个候选在此次受限记录搜索中未发现既往使用，不把“未找到”当作全历史、数据目录或同源独立资格。未读取任何候选FASTA、GFF、封存预测或标签；正式面板尚未选择／放行。

本次公开定向检索1次、请求5条，按Exa Skill计数sources_reviewed=5，不表示5份均已验证；依据是原始PDF和固定Tiberius配置。PDF网页截图首次超时，未据此改结论；复用此前从同一官方URL取得的PDF，以实际渲染页检查S18/S23/S29列关系，没有安装依赖或下载新数据。初次远程搜索误含不存在的 `config/` 返回exit2；纠正为实际三个目录后exit0，失败不作为阴性证据。正文按用户偏好的证据、范围、未证实事项写法记录，本补充不授权新实验，也不改变A11或其他活动分支的科学gate。

## 替代候选的资产资格与组装事实（15:54–15:56 UTC）

旧四项仍保留为曝光已知的迁移候选，不改写其资格。为寻找作者明确未用于Helixer训练／验证的对象，定向读取同一官方补充PDF的S18印刷页66–67、S23页81及S29页85；实际渲染表和原文一致。下列三项在S18为test (used)，并非train_val，当前a0080是论文用于植物测试的资产。因此可记录“作者披露的Helixer测试物种”，不能称全模型共同未见、全历史无同源暴露或从未用于任何方法开发。[Helixer原始补充表及数据划分方法](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-025-02939-1/MediaObjects/41592_2025_2939_MOESM1_ESM.pdf)

| 新候选／taxid／样本 | 固定RefSeq／配对GenBank | 组装／level／总长度bp | 公开注释release／日期／protein-coding genes |
|---|---|---|---|
| Vitis riparia／96939／Riparia Gloire de Montpellier，isolate1030 | [GCF_004353265.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004353265.1/)／GCA_004353265.1 | EGFV_Vit.rip_1.0／Chromosome／500105807 | Release100／2020-06-01／26189 |
| Hibiscus syriacus／106335／Baekdansim，isolateYM2019G1 | [GCF_006381635.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006381635.1/)／GCA_006381635.2 | ASM638163v2／Scaffold／1826769743 | Release100／2021-01-21／64995 |
| Nicotiana attenuata／49451／UT | [GCF_001879085.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001879085.1/)／GCA_001879085.1 | NIATTr2／Chromosome／2365682703 | Release100／2016-11-28／34094 |

三项accession／assembly名与S18作者测试身份相同；本次注释日期取实时NCBI记录，不据此证明与作者当时GFF完全相同。上表gene数是发布统计，不是本项目未来primary参考分母。Vitis与Hibiscus的RefSeq／GenBank配对分别有增补质体差异，不能交换标签与序列。

当前固定[Tiberius angiosperms训练配置](https://github.com/Gaius-Augustus/Tiberius/blob/871d4979568ae09e56de7eca5d7489046ee833ea/model_cfg/angiosperms.yaml)未列这三个物种，但列同属的V. vinifera及N. tabacum；未列名单不足以证明未见、无同源或所有上游来源干净。另发现Coffea arabica和Phoenix dactylifera虽为Helixer的test物种，却分别以Carabicav1.1、Phoenix_dactylifera_CF_009389715出现在当前Tiberius训练清单中；故没有把两者提升为“共同未见”候选。ANNEVO新Magnoliopsida的[固定说明](https://github.com/xjtu-omics/ANNEVO/blob/f1ccd9c36bd6dce3e32b25dafe50c3f5809acd0e/docs/boundary_aware_model.md)只给排除模式物种的总体陈述，并未逐项证明这三者；其个别暴露仍unknown，GENERanno同样unknown。

实际组装难点也须保留：Vitis报告174 scaffolds、19 chromosomes，Hibiscus1788 scaffolds／contig N50 2560905 bp；Nicotiana虽为Chromosome，仍报告37194 scaffolds／contig N50 64201 bp。公开assembly总长度不是核基因评测背景分母，Chromosome标签也不是全部序列都已染色体化的保证。[Vitis记录](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004353265.1/)、[Hibiscus记录](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006381635.1/)、[Nicotiana记录](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001879085.1/)

将来若冻结整个基因组评测，需按组装元数据事前确定全部核contig范围，保留合法unplaced／短scaffold、难例和零gene区，不按新模型可预测性或标签gene密度只选大染色体。当前未据长度、BUSCO或作者预测分数选择最终面板，也不宣称这些新候选覆盖了禾本科或动物；它们尚非全面适用范围证据。

在canonical项目docs／reports／configs的md／yaml／yml／json／tsv／csv中搜索这三个全名／accession／assembly／cultivar，未找到既往使用记录（rg exit1为无匹配）；该受限搜索不覆盖原始数据、旧作业和全部Git历史，不能作为独立性证明。未读取新DNA、GFF、蛋白、预测性能或任何Setaria封存内容；项目历史资格和同源隔离仍未完成。

三次NCBI Datasets原Skill脚本调用均成功，只GET `genome/accession/<上述accession>/dataset_report`，各一条记录，retrieved_at为2026-10-03T15:54:42.492–15:54:42.497 UTC，sources.kind=evidence／supports_claim=true；上表链接来自canonical_url，raw_output_path=null。复用bundled pip vendor requests（原runtime未直接提供requests），没有装包、改base或修改Skill。写法沿用用户科研记录的事实、来源、范围和未决事项；检索写作参考只用于风格，不引入其中旧TE事实。

决策：这三项可继续做明确暴露分层的独立迁移准备，而非重复旧四项全部train_val的错误前提；最终checkpoint／面板／注释政策尚未冻结，没有下载、推理或评测提交。A13两根和FULL原五根仍按各自事前合同运行，科学gate不变。只有满足本项目独立资格及统一原生口径后，结果才可能支持“下游监督物种外适用”；清洁全模型zero-shot仍另需上游manifest，不以未知推定阴性。

## 核序列范围与短 scaffold 支持

2026-10-06。三个替代候选的完整sequence report已按固定RefSeq身份分页读取，仅含序列身份、长度和组装类别等元数据，没有下载DNA、GFF或预测。若将来采用完整核基因组评价，不能把“Chromosome级组装”理解成只保留已命名染色体。

| 固定候选 | 核序列数 | assembled molecule／unplaced scaffold | 核长度bp | 小于24576bp的核序列 | 最短核序列bp |
| --- | ---: | ---: | ---: | ---: | ---: |
| [Vitis riparia GCF_004353265.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004353265.1/) | 174 | 19／155 | 500105807 | 0 | 99850 |
| [Hibiscus syriacus GCF_006381635.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006381635.1/) | 1788 | 0／1788 | 1826769743 | 8 | 7536 |
| [Nicotiana attenuata GCF_001879085.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001879085.1/) | 37194 | 12／37182 | 2365682703 | 29351 | 332 |

这里核序列对应返回记录的Primary Assembly。三个对象的unplaced记录也将assigned_molecule_location_type记为Chromosome，故不能单凭这个字段认定序列已定位到染色体。完整核范围须保留这些scaffold、短序列及零gene区域，不按长度、标签密度或预测表现删输入；目前只是准备范围，不是最终测试面板放行。

Nicotiana的unplaced序列合计1539651764bp，占核长度65.08%；仅处理12条已组装染色体将遗漏大多数核序列。其20163条核序列短于单个6144bp编码块，最短332bp，因此短contig的真实有效长度、padding、两方向尾窗及跨contig隔离是实际输入需求，不能由两个长DEV contig的完成记录代替。[固定Nicotiana序列报告所属组装](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001879085.1/)

Vitis另外返回non-nuclear叶绿体NC_039680.1／160972bp，全部175条记录合计500266779bp；Hibiscus另有NC_026909.1／161019bp，全部1789条合计1826930762bp。核长度与包括这些增补叶绿体的序列总和要分开；不把paired GenBank的差异或组装发布总长度默认为未来核背景分母。Nicotiana本次37194条均为Primary Assembly，不由报告未列细胞器推论该物种没有细胞器。[Vitis](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004353265.1/)、[Hibiscus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006381635.1/)、[Nicotiana](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001879085.1/)

A31现有infer_local_path_caller.py已支持显式--fasta、多contig及双方向，不需要再造调用包装层。原raw合成测试包含9／15bp的输入与替身encoder，但不证明真实GEN编码器在332bp等短contig上已运行验证。独立迁移前应补这一实际接口条件的原生工程预检，保留所有输入长度；不靠筛去短scaffold掩盖未验证用法，也不把工程PASS当作泛化成绩。当前A31原生DEV、评分配方和冻结gate不变。

本次NCBI Datasets原Skill脚本调用采用官方sequence_reports路径、page_size1000、不设chromosome或role过滤；Vitis／Hibiscus／Nicotiana分别完整1／2／38页，实际175／1789／37194条，逐accession无重复且与首分页total_count一致。后续页可不返回total_count；初次汇总错误读取末页该字段而停止，修正为保留首分页总数后完成余下两个对象，未将部分页当完整人口。最终证据retrieved_at分别08:21:36.575422／08:22:39.796773／08:23:31.845464 UTC，sources均evidence且supports_claim=true，checked_sources为空，raw_output_path=null。复用bundled pip vendor requests，没有安装依赖、保存raw响应或改Skill。

## 项目记录与可达 Git 历史资格

本次把原docs／reports／configs定向搜索扩展到scripts／sbatch／src的文本代码及作业脚本；三个名称、常见缩写、RefSeq／GenBank accession、assembly和cultivar只命中本候选登记，没有发现其他当前记录。另对上述目录的全部12个可达Git commit作名称／RefSeq／assembly别名的变更搜索，返回0且无命中。该历史覆盖2026-08-24初始导入至2026-09-04末commit01fbbf15544f18ef0fca19688adc1aa0d3284248；新近M28资料主要尚未提交，已由当前文本搜索另覆盖。

这个结果减少了现有项目记录中的资格疑点，不是全历史无使用或同源隔离证明：外部备份、不可达Git对象、原始数据和全部历史作业未由这次查询覆盖。三个候选仍待最终独立资格与同源条件确认，GENERanno及ANNEVO的个别上游暴露仍按原记录保留unknown；当前没有干净全模型zero-shot的结论，没有打开Setaria、sealed或任何候选预测，也没有按未来性能选择面板。

## 原生短 contig 工程预检范围

2026-10-06，运行前固定。只构造332bp与7536bp的两条人工DNA，后者把一个连续ORF放在6144bp编码块接缝附近；不是候选物种序列或注释，不评价人工ORF的生物学预测准确率。复用原FASTA reader、真实四块GEN编码器、冻结A3／A23末模型及原A31 annotate_contig、自然greedy＋beam池、owner／装配和GFF出口。原两个TRAIN准备gate及实际readout身份在加载前核对，不绕过源条件，不改变当前DEV代码或冻结科学协议。

两物种候选的最短长度用作接口需求，不据此选择最终面板。实际预检检查两个独立contig的两方向共四scope／每variant，真实4096×2048特征及24576×5 DNA张量、有效长度、有限hidden、原自然候选坐标不越出真实长度、装配计数与GFF ID／坐标／phase回读。冻结头state_dict逐tensor前后直接相等，无梯度或optimizer更新。允许人工DNA没有选中gene，不将非空输出当效果门槛；不读真实FASTA／GFF、特征缓存、DEV输出或Setaria。

独立预检script和SBATCH使用generanno、private1RTX3090 24GiB／2CPU／16GiB／5分钟，outputs/M28-TRANSFER-SHORT-RAW-PREFLIGHT及按job分离日志，预计输出低于1MiB，不覆盖旧产物。短任务由主对话消费，不加入A31固定四root监控。确定性工程错误只作一次最小修复并用独立R1输出保留首轮失败；超时／OOM或科学条件不符不自动重提。真实调用不读取面板数据，完成也只证明上述人工长度的工程接通，不是泛化、纯未见暴露、性能收益、公平速度或论文PASS。

## 原生短 contig 工程预检结果

2026-10-06 08:47 UTC，主对话已消费。作业13401642在gpu034／private-teodoro-gpu运行，实际08:41:57–08:42:55 UTC；root／batch／extern均COMPLETED／0:0，Slurm用时58秒，batch MaxRSS3071700K。STATUS为COMPLETED、JOBID为13401642，stderr为0字节。独立输出outputs/M28-TRANSFER-SHORT-RAW-PREFLIGHT／summary.json完整，没有触发修复或R1重提。

A3及A23各保留两个contig×两个方向的四个scope、四份单行JSONL，共八行；每个短contig仅一个[0,真实长度]窗口，没有把padding扩为输入坐标。真实固定revision编码器与冻结末004608模型均已调用，观察到有效长度332／332／7536／7536、输入4096×2048特征和24576×5 DNA、有限24576×128 hidden；每个头的state_dict前后逐tensor相等，无梯度或optimizer更新。7536bp正方向实际自然候选包含[6140,7535)链，来自greedy与beam4，确实跨6144bp块接缝；该候选的存在不是正确注释证据。

两个variant都只在各正方向产生一个自然候选，全部local gain非正，最终选择均为0。负方向空候选仍有保存记录，两个GFF均保留332／7536bp sequence-region。GFF回读和空预测集合一致，但没有非空gene，因此本轮不能据summary中的gff_ids_coordinates_phase_roundtrip=true宣称非空gene的ID、坐标或phase已获原生验证，也不能评价人工ORF命中率。脚本运行段18.75秒、峰值GPU allocated2410934272字节只记为这次小规模接口运行事实，不与其他方法作速度／显存优越性比较。

结论：真实短输入、跨块候选坐标、双方向／多contig记录和冻结头接通的工程条件通过；没有读取真实候选物种DNA／注释、DEV输出、特征缓存或Setaria，没有增加训练、新参数或改变现有推理代码／科学gate。两个人工长度不能代替整个碎片化基因组测试，完整面板、同源隔离及比较资产资格仍未冻结，准确率、独立泛化和论文主张仍未验证。当前A31原生DEV四根继续由同一固定heartbeat监控，13401642不加入长期监控；正式DEV结果关闭后才执行已准备的saved-output五状态诊断，不由此预检选择模型。
