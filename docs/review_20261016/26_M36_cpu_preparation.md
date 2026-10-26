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


## 7. GPU 流水线实现（2026-10-26）

范围：实现并在 ont **仅 CPU**（核 25–31，nice 10）上测试协议 22 的 GPU 阶段全部脚本（§5 启动清单 1–6 项）：A 染色体零样本推理 → S1 种子与伪背景（CPU）→ 7 次适配训练 → 评估集 B 一次性推理 → 打分与 §5 判据。驱动 `m36_queue.sh` 已在 ont 后台启动，**正在等待 `P2B_DONE`，没有使用 GPU**。**没有训练、没有任何 GPU 作业、没有读取任何 B 集输出或分数**（B 的 GFF / 参考只在 `m36_score.py` 里读取；测试里读取的是 M30 的 Carbon-A 预测与 M33 / M30 既有输出，都是已打分的旧结果）。判据、臂定义、超参数、数据划分、基模均未改动（基模按协议 22 §9 为 T10；Z 臂复用 P2b 的 T10 输出）。M34 / M30 / M33 / M35 文件**一个字节都没改**：M36 训练器在进程内 monkey-patch 包装 `m34_train_lora.main`，所以不存在需要证明"未使用时行为不变"的 M34 钩子；测试验证了补丁退出后 M34 的函数与 `CFG` 原样恢复，补丁内其他臂（如 `T10`）仍解析到 M34 定义。

### 7.0 结论先行

| 项 | 结果 |
|---|---|
| 基模 / Z 臂 | `~/coding-rna/outputs/M34-NTV2-SCALEUP/pc2_T10/model_003072.pt`（md5 `280d4efbf7a252e3f0b8adb8babb6531`；队列启动时与 `p2b/model_T10.json` 核对，不同则拒绝）；Z = P2b 的 `p2b/T10`（同一模型、同一路径，`pot/` 保留），**不重跑** |
| 实现 | 12 个新脚本 + 5 个新测试文件（§7.1）+ `m36_common.py` / `tests/conftest.py` 小改；本地与 ont `~/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT` 同路径，已 push、md5 核对 |
| 测试（ont，CPU） | 69 项全部通过（新增 60 + 既有 9），其中包括真实 `m34_train_lora.main` 循环的 CPU 冒烟（4 个臂 + 断点续跑 + 指纹拒绝）、队列干跑、Carbon-A 数值复现、判据各分支、假根端到端打分（§7.6） |
| GPU 时数（现实估计） | **≈ 58 h（43–78 h）**，协议估计 28–30 h；差别几乎全在训练步时（pc2 实测中位 41.4 s/步，协议按 NT-v2 的 18–25 s 估），见 §7.5 |
| 队列 | `m36_queue.sh` 已于 2026-10-26 17:22 在 ont 后台启动（PID 1492157，`setsid nohup`，日志 `~/coding-rna/logs/m36-queue.out`），`queue.jsonl` 最后一行为 `waiting`（原因：无 `P2B_DONE`；P2b 进程 334609 / 3174339 仍在运行），GPU 上没有 M36 进程 |
| 需用户知悉的取舍 | 7 次训练里只有 **1 个 F1 对照**，协议没写在哪个物种上做，默认 Dryas（`M36_F1_SPECIES`）；GO-强 / NO-GO / "CI 不排除持平"的读法；见 §7.2。这些只影响判据的**解读代码**，不影响任何训练或推理 |

### 7.1 文件（`experiments/M36-BOOTSTRAP-ADAPT/{scripts,tests}`）

| 文件 | 作用 |
|---|---|
| `m36_common.py`（改） | 增加 `M36_OUT` / `M36_M34_OUT` 环境覆盖（测试与干跑用假根）；基模 / Z 路径；S1 常量（M33 冻结 τ=1、thr −1.400270636、P≥0.9、CDS≥96 bp）；训练超参（协议 §3）；训练顺序 |
| `m36_s1_seeds.py` | S1 种子：零样本解码链中 M33 后验 log P ≥ thr 的链（另写 P ≥ 0.9 敏感性集，不被任何臂使用）→ 种子行（相位重算；冻结 `chain_checks` 不过的丢弃并计数）；伪背景 = `pseudo_background(P(无基因), 0.99)`，并与 `m36_nogene_marginal` 自带区段交叉核对；写 `seeds_S1_*.json`、`background_*.json`、`seed_stats_S1_*.json`、`READY_<species>` |
| `m36_cpu_queue.py` | CPU 任务队列（M33 `m33_post_queue` 的模式）。`seeds` 模式：每个 A contig × 链的 M33 后验 + 无基因边缘（`taskset 24-31`、nice 10、8 进程、内存预算 40 GB；某 contig 的 pot / chains 一出现就开算，大者优先）；一个物种全部完成后自动构建 S1 种子。`bpost` 模式：B 集各臂的后验（描述性），完成后删除该臂 `pot/`（Z 的不删） |
| `m36_train_inputs.py` | 各臂种子 / 背景读入与标签器：S1 = S1 种子 + 伪背景；S2 = 蛋白种子、无背景；S3 = `arm_seeds('S3')`（S2 ∪ S1，S2 优先）+ S1 伪背景；F1 = S1 种子 + 其余全背景 |
| `m36_schedule.py` | 抽样序列（协议 §3）："种子优先 + 均匀"各半，`default_rng(1)`，几何同 M34 `draw_schedule`（均匀半：A 训练 contig × 链 × 24,576 网格窗；种子半：均匀选种子，再在覆盖其 CDS 跨度的 ≤5 个网格窗里均匀选一个，取种子所在链）；512 步 × 8 = 4,096 次；S1 与 F1 种子集相同 → 同一序列 |
| `m36_train.py` | 导入 `m34_train_lora.main`（不改），进程内补丁 `T.arm_spec / T.arm_labelers / T.forbidden_contigs / BK.load / BK.check_init / CFG['lr_final']`。`BK.load` 外包一层，把**基模的 LoRA 载入**新建的 LoRA（M34 训练器自己是零初始化 LoRA；协议 §3 要求"继续训练头 + LoRA"）；头从基模载入。lr 头 5e-5 / LoRA 1e-4，余弦到 1e-6，5% warm-up（26 步），clip 1，batch 8，512 步；可恢复快照每 64 步，模型文件每 128 步；首次运行把基模 md5、种子文件 md5、调度 md5 写入 `m36_arm.json`，之后任何不同的重启退出码 11 |
| `m36_supervise.sh` | `m34_supervise.sh` 的 M36 版：任何中断后从 `ckpt_latest.pt` 续跑；空闲显存 < 7.5 GB 等待；连续无进展且非 CUDA/OOM 时 FATAL |
| `m36_infer.sh` | `m34_infer_retry.sh` 的 P2b 包装版（`m34_p2b_infer.py`：M30 注册后 Dryas / Esch 与 `__val20M` 可解析，解码子进程里同样）；按 (contig, 链) 续跑；`--keep-head` 保留 `pot/`；显存 < 5 GB 等待 |
| `m36_gate.sh` | P2b 闸门逻辑：`P2B_DONE` 存在，无 `m34_train/infer/p2b/stage/queue/supervise/smoke` 进程，无属于 coding-rna 的 GPU 计算进程（同账号 TE 项目的 GPU 进程不计）；每个 GPU 阶段之前再查一次 GPU 进程 |
| `m36_sentinel.py` | 一次性哨兵：首次创建并记录身份（模型 md5 等）；之后只有相同身份才放行（可续跑），不同 → 11，已完成 → 10 |
| `m36_score.py` | 评分：`main`（未过滤 exact-chain F1，判据）、`post`（后验过滤，描述）、`carbonA`（复现 M30 数值） |
| `m36_queue.sh` | 驱动（§7.2） |
| `m36_cpu_smoke.py` | 在 CPU 上跑**真实的** `m36_train.main → m34_train_lora.main` 循环（伪编码器，§7.6） |
| `tests/test_m36_{s1_seeds,schedule,train,queue,score}.py` | §7.6 |

### 7.2 流水线、顺序与协议未写明之处的取舍

`m36_queue.sh`（`flock` 单实例；事件 JSONL `~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT/queue.jsonl`；终止哨兵 `M36_DONE`）：

1. **闸门**：`P2B_DONE` 存在 + 无其他 coding-rna GPU 作业，连续 3 次、间隔 60 s 通过；再核对基模 md5 = P2b `model_T10.json` 的 md5、`pc2_T10/best_step.txt` = 3072、`p2b/T10/DONE_*.json` 存在。
2. **A 零样本**：基模对 8 个区域一次推理（Dryas `NC_141366.1 / 368.1 / 370.1 / 372.1` + `NC_141374.1__val20M`；Esch `NW_028238198.1 / 420.1` + `NW_028238472.1__val20M`；每链 **301.0 Mb**，简报写 ≈335 Mb，按 `lengths.json` 实算为 301.0），`pot/` 保留 → `zeroshot_A/`。B 染色体名不出现在 A 推理或任何训练输入里。
3. **S1 / 背景（CPU，与 GPU 重叠）**：A 推理开始前，队列就在后台启动 `m36_cpu_queue.py seeds`，每个 (contig, 链) 的 pot 与链一落盘就算后验（τ=1）和无基因边缘；一个物种全部完成后自动构建种子并写 `seeds/READY_<species>`。
4. **7 次训练，串行**：顺序 S2-Dryas、S2-Esch、S1-Dryas、S3-Dryas、F1-Dryas、S1-Esch、S3-Esch。与简报"A 零样本 → S1 → 训练"相比，**只是把依赖 S1 种子的臂排在 S2 之后**：S2 不需要 S1 种子，先训练可以让 GPU 不等 CPU 后验（后验 + 边缘共 ≈ 21–27 CPU 小时，8 核空闲时 3 h，ont 现在负载 ≈ 70 会更久）。队列每次选"顺序中第一个种子已就绪的臂"，都没就绪就等，并在 CPU 队列死掉时重启它（≤ 6 次）。顺序不改变任何臂的定义、种子、抽样序列或超参。
5. **B 一次性推理**：每个训练好的模型（step 512，主结果）**只推理自己物种的 B contig**（Dryas 4 条 102.8 Mb，Esch 3 条 170.5 Mb；不是简报写的"7 次 × 273 Mb"，7 个模型合计 922.7 Mb），哨兵 `evalB/<arm>_<species>.B_INFER_STARTED.json`（首次记录模型 md5，之后同一模型才放行）。B 集后验的 CPU 队列此时才启动（描述性）。
6. **打分**：`m36_score.py main`，哨兵 `results/SCORE_MAIN_STARTED.json`（同一模型 / 输入才允许重启，完成后再调用退出码 10）→ `results/m36_result.json`，`M36_DONE`。之后 `m36_score.py post`（后验过滤的 F1、闭合比例、AUROC / ECE；冻结 τ=1、thr=−1.400；**描述性**，独立哨兵）→ `M36_POST_DONE`。

**协议未写明、已固定的解读**（写在代码与结果 JSON 的 `verdict.rules` 里，都发生在任何 M36 输出产生之前）：

- **F1 对照放在哪个物种**：协议 §6 写"2 物种 × S1/S2/S3 + 1 个 F1 对照"= 7 次训练，没有指定物种；简报写"F1 × 2 物种"又写 7 次训练。取 7 次、F1 放 **Dryas**（A 集与 B 集都较小，训练与推理更便宜）。要改成 Esch 或加第 8 次：用 `M36_F1_SPECIES=eschscholzia_californica`（或 `M36_ORDER='…'`）重启队列（重启幂等）。`m36_score` 对每个有 F1 的物种报告 F1 − S1；"C2 撤回"在**任一**有 F1 的物种上 F1 − S1 ≥ −2.0 即成立。
- **GO-强的"CI 不排除持平"**：S2 或 S3 的 (臂 − Carbon-A 过滤) 点估计 ≥ −2.0 且配对 CI 上界 ≥ 0，两物种都满足。
- **GO**：S1 或 S3 中有一个臂在**两个物种**上都满足 闭合比例 ≥ 30% 且 (臂 − Z) 配对 CI 下界 > 0。
- **NO-GO**："任一物种闭合 < 30%" 读作：存在某个物种，S1 / S2 / S3 里没有任何臂达到 30%（结果里另给仅看 S1 / S3 的更严读法）。GO-强、GO、NO-GO 都不成立时标签为 `INCONCLUSIVE`（协议 §5 的表没有覆盖，例如闭合 ≥ 30% 但配对 CI 含 0）。
- **CI**：判据用简报指定的**配对 1 Mb 块自助法**（M34 P2b / DEV `c1_boot` 的统计量；2,000 次、`default_rng(0)`；对每个物种自己的块重采样）。协议 §4 写"1 Mb 块（含 100 kb 孤儿箱）… 与 M30/M33 同实现"，而 M30 / M33 的实现是**基因单位 + 100 kb 孤儿箱**，二者不是同一回事，所以另报 M30 的基因单位配对自助法 CI，并按它重算判据标签（`verdict_gene_unit_bootstrap`）作对照。
- **S1 种子的 CDS ≥ 96 bp**：沿用项目 ruler 的"唯一精确链、CDS ≥ 96 bp"（M33 冻结阈值就是在这个集合上标定的）；低于 96 bp 的高后验链不进种子，数量记在 `seed_stats_S1_*.json` 的漏斗里。
- **检查点**：模型文件每 128 步（128 / 256 / 384 / 512）；另每 64 步存可恢复快照，使一次中断最多损失 64 步，属工程冗余，不改变训练。主结果取 step 512，其余 step 不做 B 推理。
- **闭合比例的分母**：Carbon-A 过滤版在**同一批 B contig 上**由 `m36_score` 用 M30 的 Carbon-A 预测重新计数（Dryas 6,807 / 8,895 / 9,508 = 73.98，Esch 8,584 / 11,048 / 12,367 = 73.32，与 M30 / M33 结果文件逐项相同；测试 `test_carbon_a_numbers_equal_the_m30_results`）；Z 取 P2b 的 T10 输出。M30 的参照行（k=25% U / F、k=100%）从 M30 结果文件原样读入结果 JSON。

### 7.3 种子与后验的细节

- S1：`m33_posterior.py`（τ=1，冻结 `complete_gene_set_posteriors` + scaled-mass 查询）对 zero-shot 解码的**全部**链给 log P；保留 log P ≥ −1.400270636（P ≥ 0.2466）且 CDS ≥ 96 bp 且通过冻结 `chain_checks` 的链，写 `seeds_S1_<species>.json`；P ≥ 0.9 的子集另写 `seeds_S1_P90_<species>.json`（记数，不被任何臂使用）。真实 M33 后验切片（拟南芥内部验证染色体 NC_003071.7 '+' 链前 400 个解码链）：log P ≥ thr 且 CDS ≥ 96 的 346 条全部通过 `chain_checks`（0 条丢弃），P ≥ 0.9 的 85 条。
- 伪背景：`m36_nogene_marginal.py`（碱基级 P(无基因覆盖)，已在 §4 与穷举一致）每个 (contig, 链) 一次，`pseudo_background(…, 0.99)` 得到该链的基因组坐标区段；S1 与 S3 使用，S2 不用，F1 全背景。
- 无基因边缘与基因后验各需一次前向 + 一次后向（≈ 65–80 s/Mb，空闲核）：A 集 602 Mb·链 × 2 种 ≈ 21.7–26.8 CPU 小时。

### 7.4 B 集打分（`m36_score.py main`）

- 每个物种一个 `Evaluator`（`m30_score`，exact-chain、CDS ≥ 96 bp、不过滤）；运行：`Z`（P2b `p2b/T10`）、训练臂、`CarbonA_filt`、`CarbonA_unf`。表：tp / pred / ref / F1 / P / R、按 contig、按段数 1 / 2–5 / ≥6、CDS 碱基 P/R/F1 与背景负担。
- 闭合比例 = (F1_臂 − F1_Z) / (F1_CarbonA过滤 − F1_Z)，按物种；其 CI = (臂 − Z) 配对 CI 除以分母点估计。另给所有两两配对差（含 F1 − S1、S3 − S1、S3 − S2、臂 − Carbon-A）。
- 描述量：每臂种子数 / 与 A 训练染色体参考基因的精确链命中率（读 A 参考，仅描述）；M30 参照行；后验过滤 F1 与校准在 `post` 阶段。

### 7.5 GPU 时数：现实估计 vs 协议 28–30 h

测得的输入（ont 共享 4090，TE 项目与他人同时在用）：

- 训练步时：pc2 T10 的 3,072 步中位 **41.4 s**（均值 42.9，p10 33.5，p90 54.3；T2 中位 40.6），不是协议估计的 18–25 s。M36 的 S 臂有一部分窗口全部"未知"（S2 约一半的均匀窗口）会被训练器跳过 DP 与反传，步时可能略低；F1 臂与 M34 FULL 相同；这里按 FULL 的步时（上界）。
- 推理速度：P2b T10 在 B 集上的**实测**（正在运行，`~/coding-rna/logs/m34-p2b-T10.out`）：Dryas 4 条 B contig 共 102.8 Mb 编码用了 5,192 s，即每链 **19.8 kbp/s**（逐条 19.3 / 19.3 / 20.7 / 20.5；两条链共用一次前向，所以同一数字也可写成 39.6 kbp/s "两条链合计"）；M34 验证 / DEV 推理的范围 22.7–44.2 kbp/s 每链。简报的"≈15–20 kbp/s 两条链合计"比实测偏保守 2 倍。

| 项 | 量 | 时间（实测速度 19.8 kbp/s·链） | 范围 |
|---|---|---|---|
| A 零样本推理（每链 301.0 Mb） | 1 次 | **4.2 h** | 1.9–5.6 h（44.2 … 15 kbp/s） |
| 适配训练 | 7 × 512 = 3,584 步 | **41.2 h**（5.9 h / 次） | 34.8–54.8 h（35–55 s/步） |
| B 一次性推理（Dryas 102.8 Mb × 4 个模型 + Esch 170.5 Mb × 3 个模型 = 922.7 Mb·链） | 7 次 | **12.9 h**（Dryas 1.4 h、Esch 2.4 h / 次） | 5.8–17.1 h |
| 合计 | | **≈ 58.3 h** | **≈ 43–78 h** |

协议 22 §6 的 28–30 h = 训练 ≈ 20 h（18–25 s/步）+ B 推理 6 h + A 零样本 1.5 h。实际约为其 **2 倍**，差异来自训练步时（41 vs 18–25 s）与较慢的推理速度（19.8 vs 协议隐含的 ≈ 40 kbp/s·链）。串行约 **2.5 天**（GPU），CPU 的 S1 / 后验与 GPU 重叠；B 集后验（描述性，Z + 7 臂共 2,392 Mb·链 × 65–80 s/Mb ≈ 43–53 CPU 小时，4 进程）在 B 推理期间与之后进行，不阻塞主判据。P2b 本身（T10 + T2 各 273.3 Mb）按同一实测速度每个模型约 3.8 h（T10：15:11 起，Dryas 已于 16:38 完成，Esch 约 19:00 完成），T2 再 3.8 h，加 CPU 打分，预计 `P2B_DONE` 在 2026-10-26 约 23:00 CET（简报写 19:30，那只是 T10 推理的结束时间）。

### 7.6 测试（ont，CPU 核 25–31，nice 10；GPU 未用）

共 **69 项全部通过**：本节新增 60 项 + 既有 9 项（`test_seed_labels` 6、`test_nogene_marginal` 3，换了 `conftest.py` 之后重跑通过）。generanno 环境用 `~/envs/generanno/bin/pytest`；训练循环的 CPU 冒烟在子进程里用 plantcad2 环境（pc2 后端）运行，因为 pytest 下 spawn 进程池会重新进入 pytest（报告 20 §4）。所有数据测试只读 A 集 / 旧结果；没有读 B 集参考来做任何选择。

| 文件（项数） | 验证了什么 |
|---|---|
| `test_m36_s1_seeds.py`（7） | 常量 = M33 冻结值（τ=1、thr=−1.400270636，P ≥ 0.2466，对照 M33 `selection.json`）；阈值边界（≥）、CDS ≥ 96 bp 规则、P ≥ 0.9 子集；种子行 = 植入链（CDS + 重算相位，含 '−' 链）；`chain_checks` 不过的链被丢弃并计数；S1 标签器：种子位点有精细标签、其余 IGNORE、伪背景为背景类、F1 无 IGNORE；`build()` 在玩具 A 物种上端到端（种子数、伪背景区段 = `pseudo_background(p, 0.99)`、`READY` 标记、S1 / F1 输入）；`m36_nogene_marginal` 的区段与 `pseudo_background` 不一致时拒绝；**真实 M33 后验切片**（拟南芥内部验证染色体 NC_003071.7 '+' 链前 400 条解码链，τ=1）：346 条过阈值且 CDS ≥ 96，**0 条** `chain_checks` 失败，P ≥ 0.9 的 85 条，30 个真实种子窗上标签正确 |
| `test_m36_schedule.py`（9） | 4,096 次 = 512 × 8、确定性、种子优先 ≈ 50%、种子窗覆盖种子且同链、**均匀半与 M34 抽样律逐次相同**（同一 rng 流）、S1 与 F1 同序列、S2 / S3 不同、只抽训练 contig、B contig 被拒；真实 S2 种子（Dryas 训练 contig 3,014 条、Esch 2,000 条）的调度合法，且与种子列表顺序无关 |
| `test_m36_train.py` 甲（12） | `Patched` 钩子装上与还原（M34 的 `arm_spec / arm_labelers / forbidden_contigs / load / check_init / CFG` 退出后原样恢复；补丁内 `T10` 等 M34 臂仍解析到 M34 定义；禁止集含 B contig 与验证伪 contig）；**基模 LoRA 被载入**新建的 LoRA；LoRA 配置不一致 / 后端不一致 / 无 LoRA 的基模被拒绝；lr 日程逐步核对（26 步 warm-up、头 5e-5 / LoRA 1e-4 峰值、余弦到 1e-6，clip 1，AdamW M34 默认）；默认参数 = 协议；种子 / 基模 / 调度变化 → 退出码 11；**S1 / S2 / S3 / F1 每个臂**在合成窗上走 draw → `T.example` → 窗口 DP：计数一致、损失与梯度有限、S3 丢弃与 S2 重叠的 S1 种子、F1 无未知位置、F1 在种子位点与 S 臂逐碱基相同其余为背景；真实 Dryas 窗口：F1 已知位置 ⊇ S2 已知位置、仅 N 碱基未知、种子窗含标签位点 |
| `test_m36_train.py` 乙（6） | **真实的 `m36_train.main → m34_train_lora.main` 循环在 CPU 上运行**（伪编码器 = RC 等变伪特征 + 受 LoRA 参数影响的缩放，真实头、真实并行窗口 DP 进程池、真实 AdamW / 余弦 / 检查点）：4 个臂各训练 2 步——lr 为 5e-5 / 1e-4，写出 `DONE_2.json`、`model_000002.pt`、`ckpt_latest.pt`、`m36_arm.json`；**基模的 LoRA 出现在首次前向里**，训练后 LoRA 与头都变了，配置里 `treatment = 臂`、`lr_final = 1e-6`、`init_step = 3072`；断点续跑：2 步 + 恢复到 3 步 与不中断 3 步的 LoRA / 头一致（atol 1e-6）、第 3 步损失差 < 1e-4、恢复后首次前向的 LoRA = 快照中的 LoRA（不是基模的）；改动调度指纹后重启被拒绝（退出码 11） |
| `test_m36_queue.py`（16） | 闸门：无 `P2B_DONE` 阻塞；其他 coding-rna 进程（含 P2b 自己的 `m34_p2b_infer`）阻塞；GPU 进程只有属于 coding-rna（命令行或工作目录含 `coding-rna`）才阻塞，TE 项目的 GPU 进程与其他用户 / 已消失的 pid 不阻塞。**队列干跑**（假根，`M36_DRY=1`）：走完 gate_open → base_verified → plan → CPU 种子队列 → A 零样本 → 7 次训练（顺序 S2-Dryas, S2-Esch, S1-Dryas, S3-Dryas, F1-Dryas, S1-Esch, S3-Esch）→ B 后验队列 → 7 次 B 推理（各自物种的 B contig）→ 打分 → 后验打分 → queue_exit，只记 `would_run`，不创建 `zeroshot_A / train / evalB`、不调用 GPU；A / B 区域串等于 M30 划分；基模 md5 不同 / best_step ≠ 3072 / 缺 Z 臂 / 基模变更 → 拒绝；无 `P2B_DONE` 时持续等待并记 `waiting`，第二个实例退出码 3，`P2B_DONE` 出现后继续；CPU 队列引擎：输入（pot **和**可解析的 chains）到齐才启动、大者优先、内存预算与进程数上限、失败重试一次、组未完成时不触发种子构建钩子、B 集 `pot/` 删除钩子不删 Z 的；监督器（桩训练器）：崩溃后从检查点续跑、确定性错误 FATAL、拒绝（11）透传、GPU 忙时不启动；推理启动器（桩）：重试、`--keep-head` |
| `test_m36_score.py`（10） | **Carbon-A 过滤版在 B 集上复现 M30**：Dryas 6,807 / 8,895 / 9,508 = 73.977、Esch 8,584 / 11,048 / 12,367 = 73.321（未过滤 72.321 / 71.746，与 M30 / M33 结果文件逐项相同）；块自助法 = M34 P2b 实现（点估计、CI、frac_positive 相同，自身对自身差恰为 0）；闭合公式；判据各分支与边界值（GO-强：−2.0 / CI 上界 0 恰好满足；GO：30% 与 CI 下界 0 的边界；NO-GO 两种读法；INCONCLUSIVE；缺臂 / NaN 不算满足；C2：−2.0 边界）；一次性哨兵；**假根端到端**（"臂"用 M33 NT-v2 / M30 Gate 的既有输出，Z 用 M33 NT-v2 输出）：F1_Z = 49.92 / 44.09，S2 = S3 = Gate 输出的闭合比例 91.7% / 93.4%（(71.97−49.92)/(73.98−49.92) 等），S1 = Z 的差恰为 0、CI [0, 0]，标签 `GO`、`go_strong` 为假（Dryas Gate − Carbon-A = −2.0035，刚好在 2.0 之外），C2 判为撤回（F1 = S1），基因单位自助法给出同一标签；第二次调用退出码 10；输入（模型）变化退出码 11；缺臂退出码 2 且不写哨兵；后验过滤描述打分（构造的后验：真链 −0.5、假链 −3.0 → 过滤后 AUROC = 1、F1 高于未过滤、闭合比例按公式）及其一次性哨兵 |

**这些测试没有、也不可能验证的**（必须等 GPU）：真实 PlantCAD2 编码器 + M35 LoRA 的前向 / 反传（M36 沿用 M34 已跑通的路径，只多了"载入基模 LoRA"一步，该步用假 LoRA API 测过）、显存与真实步时、GPU 上的断点恢复、`m34_p2b_infer` 在 A 区域上的真实推理（CPU 上只核对了 `NC_141374.1__val20M`、`NW_028238472.1__val20M`、`NC_141366.1` 在主进程与解码子进程里长度和序列 md5 一致）、真实 Dryas / Esch 零样本势上的后验与无基因边缘的运行时间与内存（无基因边缘的内存按 400 B/碱基保守估计，88 Mb 的 Esch NW_028238198.1 约 35 GB，预算 45 GB 内一次只跑它加小任务）、队列非干跑路径（GPU 阶段）的首次实跑。

### 7.7 输出路径与状态命令

`~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT/`：`queue.jsonl`（事件）、`base_model.json`、`zeroshot_A/`（A 推理；`pot/` `chains/` `DONE_*.json`）、`post_A/tau1/`、`nogene_A/tau1/`、`seeds/`（`seeds_S1_*`、`seeds_S1_P90_*`、`background_*`、`seed_stats_S1_*`、`READY_*`）、`schedules/`、`train/<arm>_<species>/`（`m36_arm.json`、`training.jsonl`、`interruptions.jsonl`、`ckpt_latest.pt`、`model_000128/256/384/512.pt`、`DONE_512.json`）、`evalB/<arm>_<species>/` 与 `evalB/<arm>_<species>.B_INFER_STARTED.json`、`postB/`、`results/`（`m36_result.json`、`m36_result_posterior_filtered.json`、`SCORE_*_STARTED.json[.done]`）、`M36_DONE`、`M36_POST_DONE`。Z 臂：`~/coding-rna/outputs/M34-NTV2-SCALEUP/p2b/T10`（不复制）。日志：`~/coding-rna/logs/m36-*`。

进度：`tail -n 5 ~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT/queue.jsonl`；训练：`tail -n 3 ~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT/train/<arm>_<species>/training.jsonl`；CPU 队列：`tail -n 3 ~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT/cpuq_seeds.jsonl`。

注意：P2b 的闸门把任何名字里含 `m36_` 的进程当作阻塞项（见报告 20 §8）；在 P2b 仍在运行时若需要**重启 P2b**，会因为等待中的 `m36_queue.sh` 而互相等待，需先停 M36 队列（`pkill -f m36_queue.sh`，重启幂等）。
