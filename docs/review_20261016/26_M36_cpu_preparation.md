# 26 M36-BOOTSTRAP-ADAPT：CPU 准备（协议 22 §7 四项；2026-10-21/22）

范围：协议 `22_M36_bootstrap_adaptation_protocol.md` §7 列出的四项不依赖基模的准备——(1) miniprot 独立 env；(2) 10 个训练物种蛋白 → Dryas / Eschscholzia **标签集 A 染色体**的比对、过滤与 S2 种子 GFF；(3) 种子 GFF → M34 标签接口与单元测试；(4) 碱基级"无基因覆盖"边缘后验查询与测速。**没有训练、没有用 GPU、没有对真实目标物种做任何模型推理**；没有读取 Dryas / Esch 的参考注释（`reference.gff3`）作任何用途；评估集 B 染色体没有被拼接、建索引、比对或写出任何文件；没有读 DEV、Setaria、sealed。Baobab 未提交，Modal 未用。ont 上只用 CPU 核 24–31、nice 10、同时 ≤ 3 个进程（miniprot 单进程 8 线程）。协议 22 的判据未改，§9 无需追加（本文件只是准备记录）。

代码：`experiments/M36-BOOTSTRAP-ADAPT/{scripts,tests}`（本地与 ont `~/coding-rna-reorg` 同路径，`push.sh` 用 tar/ssh 同步；本地只新增，未 commit）。数据：ont `~/coding-rna/data/m36/{proteins,genomes_A,miniprot_raw,seeds}`；日志：ont `~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT/logs/`。

## 0. 结论先行

| 项 | 结果 |
|---|---|
| §7-1 miniprot | `~/envs/miniprot`（micromamba 独立前缀）miniprot 0.18-r281；既有 env 未动 |
| §7-2 蛋白 → 种子 | 10 物种 **241,717** 条查询蛋白（每基因最长异构体；拟南芥/水稻只用 M34 训练+验证染色体）；miniprot 比对 Dryas A 85 min、Esch A 114 min（拥挤机器）；按协议 S2 规则过滤后 **Dryas 3,674 条种子（训练染色体 3,014，28.3 / Mb）、Esch 2,309 条（训练 2,000，12.9 / Mb）**，中位 identity 0.79，≥6 段占 50%；评估集 B 未被拼接/索引/比对/写出 |
| §7-3 标签接口 | `m36_seed_labels.SeedLabeler`：与 `m34_labels.Labeler` 同签名，S1/S2/S3（种子精细标签 + 其余 IGNORE + 可选伪背景）与 F1（其余全背景）四臂；**6 个单元测试通过**，其中 F1 臂与 `m34_labels` FULL 在短柄草真实窗上逐元素相等 |
| §7-4 无基因边缘后验 | `m36_nogene_marginal.no_gene_marginal`：用冻结 scaled-mass 内核在每个边界查询背景态前/后缀质量，P(碱基 i 无基因) = exp(α_B(i)+β_B(i+1)−log Z)；**与穷举全部基因集路径的暴力计算逐碱基一致（<1e-9），log Z 与冻结后验一致**；200 kb 在拥挤核上 62 s（312 s/Mb），按 M33 生产记录（65–80 s/Mb）一条 20 Mb 染色体 × 一条链 ≈ 25 min，M36 全部 A 染色体两链 8 进程 ≈ 2–3 h，可行；未对真实目标物种运行 |
| 资源 | ont CPU 核 24–31、nice 10、≤ 3 进程；磁盘增量 **1.4 GB**（蛋白 118 MB、A 基因组+索引 0.8 GB、miniprot 输出 270 MB、种子 10 MB）；Baobab 未用、Modal 未用、GPU 未用 |
| 协议 | 判据与各臂定义未改；§9 无需追加 |

## 1. miniprot 环境（§7-1）

- `micromamba create -p ~/envs/miniprot -c conda-forge -c bioconda miniprot` → **miniprot 0.18-r281**，独立前缀 `~/envs/miniprot`，没有动 generanno / plantcad2 / 其他 env。
- 调用统一走 `scripts/m36_miniprot.sh <species> [threads] [query.faa] [tag]`：`taskset -c 24-31 nice -n 10`，先把该物种 **A 染色体**的 `data/m1_screen/<sp>/inputs/<contig>.fasta` 拼成 `genomes_A/<sp>_A.fa`（Dryas：NC_141366.1 / 368 / 370 / 372 + `NC_141374.1__val20M`；Esch：NW_028238198.1 / 420 + `NW_028238472.1__val20M`；脚本里写死，B 染色体名不出现在任何脚本中），`miniprot -d` 建索引，再 `miniprot -I --gff --outc=0.9 --outs=0.97 --outn=5` 输出 GFF3（含 `##PAF` 行）并 gzip。`-I` 把最大内含子设为 3.6·√(参考长度)（Dryas A ≈ 40 kb，Esch A ≈ 52 kb）；`--outn=5/--outs=0.97` 只保留每条蛋白接近最优的 ≤ 5 条比对（同位点去重只需要最高 identity 的那条，取 5 条是为了让不同蛋白在同一位点竞争）。

## 2. 蛋白来源与过滤规则（§7-2）

### 2.1 查询蛋白（`m36_proteins.py`）

- 来源：M34 的 10 个训练物种 的 RefSeq 蛋白（ont `data/m34_train/raw/*/…_protein.faa.gz` 与 `raw/legacy_proteins/`），md5 已在 19 号文档核过。**不含**目标物种自己的蛋白、不含面板蛋白。
- 规则（先于比对固定）：每个蛋白编码基因只取 **CDS 最长的一个异构体**（从同一 RefSeq GFF 的 `CDS.protein_id` → `mRNA` → `gene` 解析；并列取 ID 字典序小者）；去掉 `partial=true` 的 CDS；去掉细胞器 contig（GFF `region` 的 `genome=chloroplast|mitochondrion|plastid`）；拟南芥 / 水稻只取 **M34 训练 + 内部验证染色体**上的基因（Ath NC_003070.9 / 003076.8 / 003071.7，Osa NC_089043.1 / 089046.1 / 089044.1），即不碰 DEV（NC_003074.8 / NC_089041.1）与 split_test 染色体；8 个新物种取全部核基因（含未定位 scaffold）。序列去掉末尾 `*`，非字母字符改 `X`，长度 < 30 aa 丢弃。头部 `>{species_key}|{protein_id}`。

结果（`data/m36/proteins/queries_stats.json`，212 s）：**241,717 条查询蛋白**，`queries_10sp.faa` 118 MB。

| 物种 | GFF 中蛋白 | 基因（取最长异构体） | 写出 | 因染色体限制丢弃 | 细胞器 | partial | 残基数（M） |
|---|---:|---:|---:|---:|---:|---:|---:|
| 拟南芥（仅 3 条 M34 训练/验证染色体） | 48,265 | 17,756 | 17,728 | 16,900 | 118 | 65 | 7.2 |
| 水稻（仅 3 条） | 13,438 | 4,980 | 4,980 | 6,236 | 163 | 13 | 2.3 |
| Populus | 51,840 | 29,451 | 29,451 | 0 | 98 | 68 | 13.7 |
| Brachypodium | 37,892 | 25,367 | 25,367 | 0 | 81 | 87 | 11.4 |
| Solanum | 44,391 | 28,440 | 28,440 | 0 | 126 | 46 | 12.3 |
| Phoenix | 54,077 | 25,768 | 25,768 | 0 | 0 | 62 | 12.0 |
| Phaseolus | 43,271 | 29,488 | 29,488 | 0 | 114 | 82 | 13.0 |
| Citrus | 40,427 | 23,363 | 23,363 | 0 | 120 | 73 | 10.8 |
| Daucus | 49,321 | 32,885 | 32,885 | 0 | 120 | 74 | 15.3 |
| Dioscorea | 33,298 | 24,247 | 24,247 | 0 | 84 | 231 | 10.4 |
| 合计 | | | **241,717** | | | | 108.4 |

（拟南芥 17,756 → 17,728 的差为 28 条 < 30 aa 的蛋白；水稻 4,980 个基因 = M34 的 2 条训练染色体 3,377 + 验证染色体 1,601 个主编码基因（19 号文档 §2.4），`m1_screen/oryza_sativa/reference.gff3` 只含 M28/M31 用到的 5 条染色体。）

试跑：2,000 条短柄草蛋白 → Dryas A（128 Mb，索引 341 MB）**20 s 实时 / 41 s CPU（8 线程，峰值 0.65 GB）**，输出 gz 0.97 MB、1,183 条比对；据此全量 241,717 条约 40 min（Dryas）与 ≈ 1 h（Esch，209 Mb），输出每物种 ≈ 120 MB。

### 2.2 S2 过滤（`m36_seed_filter.py`，协议 22 §2 的 S2 行）

对每条 miniprot 比对（`##PAF` + mRNA/CDS/stop_codon）依次要求：

1. **全长覆盖** ≥ 0.95：PAF 的 (qend − qstart) / qlen；
2. **identity** ≥ 0.70（mRNA 的 `Identity`）；
3. **无移码、无内部终止**（PAF `fs:i:0`、`st:i:0`；新版 mRNA 的 `Frameshift=0;StopCodon=0`）；
4. **有终止密码子**：miniprot 的 `stop_codon` 特征紧接末端外显子（与 CDS 合并，得到含终止子的完整 CDS 链；终止子被内含子打断者按"split_terminal_codon"落入第 5 条拒绝）；
5. **合法链**：冻结 `manifest.chain_checks`（ATG 起始、终止子、长度 mod 3、无内部终止、GT/GC–AG 剪接、仅 ACGT、首末外显子 ≥ 3 bp）全部通过，外显子数 ≤ 128（语法上限）；
6. **同位点取 identity 最高**：同一 contig、同一链上 CDS 跨度重叠的候选，按 (identity, coverage, CDS bp) 降序贪心保留。

反义链重叠的种子彼此保留（标签按链分别生成，与 M34 一致），只在统计中报告数量。输出 `seeds_S2_<sp>.gff3`（gene/mRNA/CDS，CDS 含终止子，属性含 source_protein / source_species / identity / coverage / n_exons）、`seeds_S2_<sp>.json`（接口直接读）与 `seed_stats_S2_<sp>.json`（漏斗、按 contig / 训练-验证段 / 段数分层 / 来源物种 / identity 分位数）。

### 2.3 种子统计

全量运行（ont，2026-10-21 23:26 → 10-22 02:46，核 24–31、nice 10、8 线程；期间机器负载 70+/32 核，CPU 利用率仅 1.3–1.5 核）：

| 目标 | A 染色体总长 | miniprot 实时 / CPU / 峰值内存 | 输出（gz） | 比对数 |
|---|---:|---|---:|---:|
| Dryas A（4 训练 + val20M） | 126.5 Mb | 5,120 s / 7,278 s / 2.7 GB | 136 MB | 168,750 |
| Esch A（2 训练 + val20M） | 174.5 Mb | 6,814 s / 8,600 s / 3.2 GB | 133 MB | 162,787 |

过滤漏斗（`seed_stats_S2_<sp>.json`；过滤本身每物种 15 s）：

| 步骤 | Dryas | Esch |
|---|---:|---:|
| miniprot 比对（≤ 5 条/蛋白，≥ 0.9 覆盖，≥ 0.97 最优分） | 168,750 | 162,787 |
| 覆盖 ≥ 0.95 | 141,081 | 133,463 |
| identity ≥ 0.70 | 43,366 | 38,151 |
| 无移码 / 内部终止 | 42,331 | 37,370 |
| 有终止密码子特征 | 28,998 | 24,697 |
| 冻结 `chain_checks` 全通过（失败：非 ATG 起始 9,240 / 8,358；非 GT/GC-AG 剪接 224 / 163） | 19,590 | 16,213 |
| 同位点去重（丢弃 identity 较低者 15,916 / 13,904） | **3,674** | **2,309** |

种子分布：

| | Dryas 训练染色体（106.5 Mb） | Dryas val20M | Esch 训练染色体（154.5 Mb） | Esch val20M |
|---|---:|---:|---:|---:|
| 种子数 | **3,014**（28.3 / Mb） | 660 | **2,000**（12.9 / Mb） | 309 |
| 按段数 1 / 2–5 / ≥6 | 380 / 1,121 / 1,513 | — | 278 / 685 / 1,037 | — |
| 按 contig | NC_141366.1 1,033；368 761；370 643；372 577 | | NW_028238198.1 1,093；420 907 | |
| identity 分位 5/25/50/75/95% | 0.709 / 0.745 / 0.791 / 0.848 / 0.931（全部 A） | | 0.708 / 0.744 / 0.791 / 0.858 / 0.945 | |
| 来源物种（前三） | Populus 1,399、Citrus 1,310、Phaseolus 485（禾本科 14、拟南芥 68） | | Populus 527、Citrus 512、Phoenix 313（Daucus 263、Solanum 203、Dioscorea 198） | |
| 剪接基序 | GT-AG 22,112、GC-AG 301 | | GT-AG 13,661、GC-AG 137 | |
| 反义链重叠种子 | 1 | | 0 | |

观察（描述性，不读参考注释）：
- 种子密度相对 M30 报告中评估集 B 的参考基因密度（Dryas 9,508 / 102.8 Mb ≈ 92 / Mb；Esch 12,367 / 170.5 Mb ≈ 73 / Mb，14/18 号文档）约为 **31%（Dryas）与 18%（Esch）**——这是 S2 的"部分标签"规模，介于 M30 的 k=10% 与 k=25% 之间（Dryas 偏近 25%，Esch 偏近 10%），且多外显子基因占比高（≥6 段占一半），单外显子基因明显偏少（蛋白比对对短的单外显子基因的 identity 门槛更难达到）。
- Dryas（蔷薇目，fabids）种子几乎全部来自 Populus / Citrus / Phaseolus 等真双子叶蛋白，禾本科蛋白只贡献 14 条；Esch（基部真双子叶）来源更分散，单子叶（Phoenix、Dioscorea）贡献 22%。与协议的"跨物种、不用目标物种自身蛋白"一致。
- 最大的损失发生在 identity ≥ 0.70（丢 69–71%）与"非 ATG 起始"（miniprot 把蛋白第 1 位对到非 ATG 密码子，或起始段未对齐）；去重丢弃的 1.6 万 / 1.4 万条是同一位点被多个物种的同源蛋白重复命中，正是预期的冗余。
- 首次运行的过滤脚本把 miniprot 0.18 的 `stop_codon`（已包含在末端 CDS 内）当成 CDS 之外的特征，导致 0 条种子；修正后重跑（比对结果不变，只重跑 15 s 的过滤）。

接口冒烟（真实 Dryas / Esch A 训练染色体，`target_labeler(sp, 'S2')`）：GFF3 与 JSON 读入逐条相同；两物种全部种子对基因组重跑 `chain_checks` 0 条不合法；各取 40 个种子窗，S2 臂的已知位点与 F1 臂逐碱基相同，F1 在全部 ACGT 碱基上无 IGNORE，S2 窗内 IGNORE 比例平均 0.83（Dryas）/ 0.90（Esch）。

## 3. 种子 GFF → M34 标签接口（§7-3，`m36_seed_labels.py`）

- `SeedLabeler(species, genome_fa, lengths, seeds, background=None)`：与 `m34_labels.Labeler` / `CatalogLabeler` 同一组属性与方法（`seqs`、`lengths`、`grid`、`window`、`local`、`targets(contig, strand, a, b, treatment)` → `(dict(seg, region), win, core)`），因此 `m34_train_common.example(lab, treatment, row)` 和 M34 训练器可直接消费；`seed_cells()` 给出"种子优先"抽样所需的 (contig, strand, 窗) 单元（与 `C.gene_cells` 同构）。种子行由 `seed_row()` 构造成目录格式（`training_primary=True, CDS_complete_supported=True, overlapping_primary=False, eligible=True`），精细标签由冻结 `src.m28.labels.build_targets` 生成，与参考基因完全同一代码路径。
- 四种臂：
  - **S1 / S2 / S3 / S**：窗内种子位点得到完整精细标签（CDS 相位、donor/acceptor、start/stop、内含子相位；`region` 1/2），**其余位置 seg 与 region 都为 IGNORE（边缘化）**；可选伪背景区段（`background={(contig, strand): [(a, b)]}`，来自 §4 的 P(无基因) ≥ 0.99 区段）标为背景 0/0；种子位点优先于伪背景。
  - **F1 / F / FULL**：种子位点同上，**其余全部背景**（0/0），除 N 碱基与窗尾填充外无 IGNORE——即 Carbon-A 式填背景对照。
- 臂的种子集合：`arm_seeds(arm, s1_rows, s2_rows)`：S1/F1 → 自训练种子；S2 → 蛋白种子；S3 → `merge_seed_sets(S2, S1)`（同链跨度重叠时丢 S1，S2 优先，返回丢弃数）。`target_labeler(species, arm, s1, s2, background)` 在 M30 的 A 训练染色体上直接构造（验证伪 contig 需显式指定；B 染色体名会触发 `assert_not_B` 异常）。
- 读入：`load_seed_gff`（任何含 mRNA/CDS 的 GFF3，CDS 含终止子，相位按翻译顺序重算）、`load_seed_json`。伪背景：`pseudo_background(p_oriented, strand, L, threshold=0.99, min_len=1)` 把某条链的每碱基 P(无基因) 转成基因组坐标区段。

### 3.1 单元测试（`tests/test_seed_labels.py`，ont generanno env）

`6 passed in 29 s`（合成部分 5 个 + 数据部分 1 个）：

| 测试 | 内容 | 结果 |
|---|---|---|
| `test_gff_round_trip` | 合成 6 个两链多外显子基因写成 GFF3 / JSON 再读回：id、链、CDS 三元组（含重算相位）逐一相等 | 通过 |
| `test_seed_roles_and_unknown_mask` | 合成 contig 全部网格窗 × 两链：S2 臂在每个种子跨度内的 `seg`/`region` 与"由 GFF 链直接推出的角色"（CDS 相位 1–3、内含子 4–6、donor 7–9、acceptor 10–12、start 13、stop 14；region 1/2）逐碱基相等；跨度外全部 IGNORE；窗尾填充 IGNORE；种子碱基无 IGNORE | 通过 |
| `test_F1_has_no_ignore_and_matches_S_on_seeds` | F1 臂在全部有效碱基上 seg、region ≥ 0（**无 IGNORE**）；与 S1 臂在种子位点逐碱基相同；其余位置为背景 0 | 通过 |
| `test_pseudo_background_overlay` | 给定 (contig, 链) 伪背景区段后，S1 臂在这些区段为 0/0，其外仍 IGNORE；'−' 链区段的方向/基因组坐标换算正确；`pseudo_background` 从概率向量得到同一区段 | 通过 |
| `test_merge_S3` | S3 合并：与 S2 同链重叠的 S1 种子被丢弃（计数 1），其余保留；`arm_seeds` 的 F1 用 S1 种子 | 通过 |
| `test_equals_m34_full_on_brachypodium` | 用 `m34_labels.Labeler`（短柄草 8 Mb 验证段，真实 GFF）目录中的 eligible 基因作为"种子"喂给 `SeedLabeler`：在两链各 4 个只含 eligible 基因的窗上，**F1 臂的 `seg`/`region` 与 M34 FULL 逐元素完全相等**（同一 `core` 序列）；S2 臂在种子位点等于 FULL，其余 IGNORE | 通过 |

## 4. 碱基级"无基因覆盖"边缘后验（§7-4，`m36_nogene_marginal.py`）

- 推导：冻结自动机（`full_grammar_boundary_messages._transitions`）中只有背景自环 B→B 发射角色 0，任何基因都是一次 B→…→B 的闭合游走，所以
  log P(碱基 i 不被任何基因覆盖) = α_B(i) + gain[i,0] + β_B(i+1) − log Z，其中 gain[i,0] ≡ 0，α_B(i) / β_B(j) 就是冻结 `_prefix_queries` / `_suffix_queries` 在查询位置报告的背景态前/后缀质量。**新增的只是"在 0…n 的每个边界都查询"**：用冻结的 scaled-mass 内核 `_scaled_prefix_queries` / `_scaled_suffix_queries`（M33 后验同款，块重标定防下溢）各跑一遍，一次前向 + 一次后向，代价与 M33 的基因后验相同；不改任何冻结文件，不剪枝，不引入阈值。
- 接口：`no_gene_marginal(dna_oriented, logits[n,15], max_exons=128)` → `log_partition, log_p_nogene[n], p_nogene[n], alpha_B, beta_B`；`segments(p, 0.99, min_len)`、`to_genomic(...)`。命令行 `run <run_dir> <seq> <strand> <tau> <out_dir> [--fasta]` 读 `m33_infer` 格式的 `pot/<seq>_<strand>.npy`，写 `nogene_<seq>_<strand>.npy`（float32 P）与 `.json`（log Z、P ≥ 0.99 的伪背景区段的方向/基因组坐标、覆盖 bp、秒）；`bench <L>` 合成测速。**没有对 Dryas / Esch 的任何势运行**（基模未定，M34 之后才会有 pot）。

### 4.1 测试（`tests/test_nogene_marginal.py`）

`3 passed in 101 s`：

| 测试 | 内容 | 结果 |
|---|---|---|
| `test_brute_force_cap3_and_cap128` | 5 条 72–96 bp 合成序列（植入 1 或 2 个含/不含内含子的合法基因），对冻结自动机 **DFS 穷举全部 B→…→B 路径（= 全部基因集及其外显子结构）**，按路径能量求 log Z 与每碱基"角色 0 质量"，与 `no_gene_marginal` 比较：log Z 相对误差 < 1e-9，每碱基 P(无基因) 绝对误差 < 1e-9；上限 3 与生产上限 128 结果一致。路径数 7–24，其中两基因集 6–10 个 | 通过 |
| `test_against_frozen_posterior_and_chains` | 6 kb 合成问题：log Z 与冻结 `complete_gene_set_posteriors` 一致（< 1e-6）；对 skyline 解码的 18 条链，链内每碱基满足 P(无基因) ≤ 1 − P(链)（最高 P(链)=0.998） | 通过 |
| `test_segments_round_trip` | 阈值区段、min_len、'−' 链到基因组坐标的往返、与 `pseudo_background` 一致 | 通过 |

### 4.2 速度（整条 20 Mb 染色体的可行性）

- 实现代价：状态数 19,078（上限 128），每碱基一次 `_advance`，前向 + 后向各一遍；内存：势 [20 M,15] float32 1.2 GB + 双精度 gain 2.4 GB + 三个长度 20 M 的 float64/int64 数组 0.5 GB ≈ **4–5 GB / 进程**（M33 已对 68 Mb 的 Esch B 染色体跑过同类两遍扫描，内存可行）。
- 合成测速（`bench 2e5`，ont 核 24–31 **在负载 72/32 核、同核上 miniprot 并发的条件下**，nice 10）：200 kb 用 62 s = **312 s/Mb → 20 Mb 外推 ≈ 104 min / 链**；同一输入冻结 `complete_gene_set_posteriors`（只查两个端点）用 157 s，两者 log Z 相差 1.2e-11（553.4549）。合成势下平均 P(无基因)=0.975，P ≥ 0.99 的碱基占 93.6%（随机势、背景偏置 +0.5 的人工数字，不代表真实模型）。2 Mb 的测速在同样拥挤的条件下被我中止以让出核给 miniprot（下次在空闲时段补测）。
- 实际参考：M33 用同一对内核（scaled-mass 前向 + 后向）对整条染色体算基因后验的生产记录为拟南芥 NC_003071.7（19.7 Mb）1,527–2,045 s / 链、Esch NW_028238309.1（68 Mb）4,218–4,465 s / 链，即 **≈ 65–80 s/Mb**（空闲核）。本查询不增加扫描次数，只多写两个长度 n 的数组，因此一条 20 Mb 染色体 × 一条链在空闲核上 ≈ 25 min，拥挤时 ≈ 1.5–2 h；M36 需要的 Dryas A（4 条训练 + 20 Mb 验证段 ≈ 126.5 Mb）× 2 链 + Esch A（≈ 174.5 Mb）× 2 链 ≈ 600 Mb·链，按 8 个进程并行（M33 的 `m33_post_queue.py` 同款）约 1.5–3 h，可行。

## 5. 待 M34 阶段一结束后的启动清单

1. 基模：`M34_STAGE1_DONE` 出现后取 T10 胜出 checkpoint（骨干由 M35 决定：`M34_BACKEND=ntv2|pc2`），确认其 `m34_infer.py` 能对 M30 注册的 Dryas / Esch contig 推理（`m30_common.register()` + `m34_train_common.register()`）。
2. 零样本推理 **只在 A 染色体**（train 4+2 条 + 两个 `__val20M` 伪 contig）上生成 `pot/`、`chains/`（M33 格式）；B 染色体的推理推迟到 §4 评估时一次完成。
3. S1 种子：`m33_posterior.py` 同款后验（τ=1，M33 冻结阈值 log P ≥ −1.400；另记 P ≥ 0.9）→ 链转成 seed rows（`seed_row(seqid, strand, CDS(含终止子，相位重算), id)`），写 `seeds_S1_<sp>.json`；伪背景：`m36_nogene_marginal.py run …` 每条 A contig × 两链（Dryas ≈ 5×2，Esch ≈ 3×2 任务，§4.2 速度），`pseudo_background(…, 0.99)` → `background` 字典。
4. S2 种子已就绪（`seeds/seeds_S2_<sp>.json`）；S3 = `arm_seeds('S3', s1, s2)`；F1 用 S1 种子 + `treatment='F1'`。
5. 训练：复用 `m34_train_lora.py` 的 draw → `example` 路径，把 `arm_labelers` 换成 `target_labeler(species, arm, s1, s2, background)`（每个目标物种单独训练），抽样"种子优先 + 均匀"各半（`seed_cells()` + A 训练染色体网格），512 步、batch 8、头 lr 5e-5 / LoRA 1e-4、每 128 步存 checkpoint（协议 §3 固定）；内部验证伪 contig 只作描述。
6. 评估：B 染色体一次性推理 + `m33_evalB.py` 口径（exact-chain F1、分层、1 Mb 块配对 bootstrap 2,000 次），闭合比例按协议 §4，判据 §5；B 的 GFF 只在评分脚本读取。
7. 记录：每臂种子数 / 覆盖参考基因比例（评分阶段事后统计）、GPU 时数；偏离只在协议 22 §9 追加。

## 6. 复现

```
# ont, CPU
cd ~/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT
taskset -c 24-31 nice -n 10 ~/envs/generanno/bin/python scripts/m36_proteins.py
bash scripts/m36_miniprot.sh dryas_octopetala 8 ; bash scripts/m36_miniprot.sh eschscholzia_californica 8
taskset -c 24-31 nice -n 10 ~/envs/generanno/bin/python scripts/m36_seed_filter.py dryas_octopetala
taskset -c 24-31 nice -n 10 ~/envs/generanno/bin/python scripts/m36_seed_filter.py eschscholzia_californica
taskset -c 24-31 nice -n 10 ~/envs/generanno/bin/python -m pytest -q -s tests/
taskset -c 24-31 nice -n 10 ~/envs/generanno/bin/python scripts/m36_nogene_marginal.py bench 2e6
```
