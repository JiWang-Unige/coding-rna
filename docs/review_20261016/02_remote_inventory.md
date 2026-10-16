# 02 远端 Baobab 状态盘点（只读，2026-10-16）

范围：`ssh baobab` 只读检查 `/home/users/j/jwang/coding-rna`。未修改任何远端文件，未提交作业，未对 outputs/ 整体做 du，未读取 secrets.env 内容（仅 ls 看到其存在与权限 600）。未测量的项目均已标明。

## 0. 要点

- 远端项目目录实际是符号链接：`/home/users/j/jwang/coding-rna` -> `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna`（beegfs scratch，非 /home）。scratch 847T/1005T 已用（85%），可用 158T。`/home` 本身 98%（135T/138T）。
- 远端 git HEAD 停在 `01fbbf1`（2026-09-04），之后的全部工作（M25R-E1 到 M28，至 2026-10-14）在远端都是未提交的工作区改动。
- 远端活跃研究已到 M28 的 "full-grammar" 线，最新产物时间 2026-10-13 18:44。本地副本的文档只到 M28 早期（`reports/M28-METHOD-RESTART` 本地 180 个文件，远端 269 个）。
- 当前无 jwang 排队/运行作业（squeue 为空）。
- credits：PI teodoro 账户余额 147,263（总额 3,520,780 的约 4.2%），有效期到 2026-11-30。jwang 10 月 1 日到 16 日 billing 用量 101,800。

## 1. 文档与配置

### 1.1 AGENTS.md（远端，2026-09-08）
- 项目于 2026-09-08 退出 Auto Research/lwcr 框架；用户驱动协作，简体中文。
- 保留科学约束：数据划分、同源隔离、sealed/held-out、冻结阈值和预算；取消框架不等于授权首训、新实验、放行封存数据或发布。
- Baobab 上训练/评估/重数据处理必须走 Slurm；登录节点只做轻量操作；使用项目对应环境；不默认提交 git、取消作业或建后台监控。
- 远端 CLAUDE.md 仅指向 AGENTS.md（与本地一致）。

### 1.2 PROJECT_STRUCTURE.md（2026-06-13）
旧框架时期的目录契约（scripts、configs、sbatch、runs、reports、logs、software_outputs、external_runs、goals、refs、wiki、worktrees、data/{raw,interim,processed} 等）。与现状不符：runs、software_outputs、external_runs、manuscript、pipelines、worktrees、templates、analysis/notebooks 在远端都只有 `.gitkeep`；实际产物几乎全在 `outputs/`，数据在 `data/` 下按实验命名的子目录。

### 1.3 README.md（2026-09-29 修改）
- 公开精简快照说明；固定 M25R 路线 NO-GO；用户已明确重启 M28 新模型开发（目标 Nature Communications 及以上）；M27 关闭；Setaria 保持封存。
- README 末尾仍写着 R4 pilot "训练已开始"，落后于实际进度（实际已到 M28 full-grammar 线），以 docs/reports 为准。
- 另有 `README.md.orig`、`$LOG`（46KB，6 月）、`cluster_config.yaml.example*` 及 backup 文件，属杂物。

### 1.4 cluster_config.yaml（2026-09-08）
- 集群 baobab，on_cluster 模式；conda init `/opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh`。
- 配置里 conda_env 为 `coding-rna`（`/home/users/j/jwang/.conda/envs/coding-rna`）。实际 sbatch 使用情况：604 个脚本 `conda activate generanno`，12 个 `annevo`，仅 1 个 `coding-rna`。即 M28 的事实环境是 **generanno**，不是 cluster_config 里写的 coding-rna。
- 永久排除节点：gpu023、gpu024、gpu036 到 gpu043。
- 分区：`private-teodoro-gpu`（gpu034/035，RTX 3090 24GB x8/节点，7 天，免费，优先级 1）；`shared-gpu`（12h，计费，A100 40/80GB、A6000、RTX 3090/4090 等）；`shared-cpu`（12h，计费）；`public-cpu`（96h）。
- 硬限：max_concurrent_jobs 80，max_array_size 100，默认 8 CPU、32GB/GPU，signal USR1@120。
- 实际 sbatch 分区统计：633 个脚本中 624 个 private-teodoro-gpu，7 个 shared-gpu；gres 以 `nvidia_geforce_rtx_3090:1` 为主（195 个 M28 脚本），个别 RTX 4090（shared gpu049）。
- cluster_config 中 hard_limits 的 orthogonality 等条目属 Auto Research 遗留，现已不是流程要求。

## 2. Git 状态与本地/远端关系

### 2.1 远端
- 分支 main，跟踪 origin（`git@github.com:JiWang-Unige/coding-rna.git`），`up to date with origin/main`（指远端本地缓存的 origin/main，也是 01fbbf1），总提交数 12。
- 最近提交：`01fbbf1 Repair M25R Stage 1 R4 diagnostics`（2026-09-04）、`ce8f795`、`05b4f4e`、`a161fc9`、`1641ff7`、`99e51d6`、`9c654f5 Record ChatGPT Pro systematic review`、`fe38a67 Archive M25R result and research roadmap`、`a7bfb3b`、`7509e0a`、`6dd7190 Trim repository`、`69f06ba Initial import`。
- 工作区：`git status --short` 共 924 条 = 879 个未跟踪（??）+ 40 个已删除（D）+ 5 个已修改（M）。
  - 已修改（5）：docs/12_publication_strategy.md、docs/15_evidence_register.md、reports/M20-SOTA-ERROR-ANALYSIS/report.md、scripts/experiments/M20-SOTA-ERROR-ANALYSIS/analyze_sota_errors.py、scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/redecode_error_decomposition.py。
  - 已删除（40）：sbatch/M10 到 M23 共约 29 个启动脚本、scripts/run_*.sbatch（9 个旧 M1 到 M7 脚本）、scripts/validate_goal.py、tests/test_validate_goal_profiles.py。这与 Auto Research 退休时把旧脚本移出项目一致（见 docs/auto_research_retirement_2026-09-08.md）。`git diff --stat`：45 files changed，90 insertions，4605 deletions。
  - 未跟踪（879）：几乎全是 M25R-E1、M26、M27、M28 的 docs、reports、sbatch、scripts（远端 sbatch 目录共 633 个脚本，其中 M28 597 个、M26 10、M27 9、M25R 15）。

### 2.2 本地（`/Users/jiwang/Desktop/coding RNA`）
- main 在 `1d2593c`（2026-09-26），共 15 个提交；`origin/main`（本地缓存，可能已过时）在 `d77d225`（2026-10-01 "Record B2 job chain and independent lightweight monitor"），本地落后 39 个提交；`remotes/baobab/main` = `01fbbf1`。
- 本地比远端 HEAD 多 3 个提交：`bfa3533`（M25R 固定路线关闭）、`875c5af`（M26 同范围比较）、`1d2593c`（参考背景 FPR 归因）。merge-base 即 `01fbbf1`。
- 本地工作区：3 个 M + 656 个 ??（含 reports/M28-METHOD-RESTART 180 个文件、docs 22 个、sbatch 等），本地没有 outputs/ 真实产物目录（仅精简）。

### 2.3 分叉关系结论
- 历史是**线性**的，没有分叉：`01fbbf1`（远端 HEAD）-> `bfa3533` -> `875c5af` -> `1d2593c`（本地 HEAD）-> 另外 36 个提交 -> `d77d225`（GitHub origin/main，据本地缓存，截至 10-01）。
- 远端集群 checkout 自 9 月 4 日起没有再提交或 pull；M25R 之后的结果文件是直接在集群上生成的（未跟踪），而提交是在本地/GitHub 侧完成的。因此远端 `git log` 不能反映进度，远端 `git status` 里的 879 个未跟踪文件中大量内容在 GitHub 上已是已提交状态（如 M26、M27、M28 前期）。
- GitHub 上 10-01 之后（M28 A 系列后期到 full-grammar 线，10-14 前后）是否还有提交：**未验证**（未执行 git fetch，避免修改本地 refs）。从远端文件看，full-grammar 线的 reports/docs（10-11 到 10-14）很可能只存在于集群（及可能的手工同步）。
- 注意：local origin/main 缓存的时间戳较旧，需要 `git fetch` 才能确认。
- 风险：远端无 git 保护的未提交内容很多（924 条），且 `.gitignore` 排除 data/、outputs/、*.pt、refs/weights 等，所有 checkpoint 只存在于集群 scratch。

## 3. 顶层目录盘点

| 目录 | 状态/用途 |
|---|---|
| `AGENTS.md` / `CLAUDE.md` | 协作规则（见 1.1） |
| `configs/` | M10 到 M17 的 yaml（6 月到 8 月，另有 M25/M25R 配置），M18 之后基本不再用 yaml，M28 参数写在脚本与 sbatch 里 |
| `data/` | 见第 4 节 |
| `docs/` | 25 个主文档 + figures/；最新为 M28 full-grammar 的 manuscript/methods/results/discussion/introduction/abstract working 文档（10-11 到 10-14） |
| `reports/` | 17 个实验报告目录；`M28-METHOD-RESTART` 269 个文件（最新 10-14 07:02） |
| `outputs/` | 788 个条目（786 目录 + 2 文件），约 730G（题述，未测量）；见第 5 节 |
| `sbatch/` | 633 个脚本（M28 597） |
| `scripts/` `src/` `tests/` | 代码；`scripts/experiments/<exp>/` 含 M20、M25R、M26、M27、M28；`tests/` 含 M28 相关新增测试（未统计数量） |
| `refs/` | 33G（du 实测）；pdfs（4 篇论文）、repos（annevo/generanno/helixer/segmentnt/tiberius 源码）、containers（helixer-docker_latest.sif）、weights（Helixer 全部 clade 权重 107M、Tiberius m27 nosm v2 权重 60M）、dossiers（来源/污染审计 9 份）、carbon-a（见 6.3） |
| `external_runs/` `software_outputs/` `runs/` `worktrees/` `manuscript/` `pipelines/` `templates/` `analysis/notebooks/` | 仅 `.gitkeep`，**空**。基线运行实际落在 `outputs/BASE-*`、`outputs/M12B-*`、`outputs/M17-*` |
| `experiments/` | 仅有空目录 `M28-METHOD-RESTART` |
| `goals/` | 9 对旧 goal md（M3 到 M9），属 Auto Research 历史，不再是流程要求 |
| `logs/` | 仅 2 对 M25R 取证日志（9 月 8 日）；实际 Slurm 日志在 `outputs/<exp>/logs` 或 `outputs/M28-PILOT-R5/logs` 等 |
| `wiki/` | ideas/notes/INDEX.md（6 月） |

## 4. data/ 与 Setaria 盲测数据

| 路径 | 内容 |
|---|---|
| `data/m1_screen/` | 6 个物种：arabidopsis_thaliana（genome.fa 121MB、reference.gff3 456MB、split_{train,val,test}.jsonl）、arabidopsis_lyrata、oryza_sativa（含 genome.full.fa / reference.full.gff3）、gallus_gallus、drosophila_melanogaster、saccharomyces_cerevisiae；另有 check_data_report.json、download_report*.json、m7_prep_report.json |
| `data/_refseq_cache/` | RefSeq 原始下载：GCF_000001215.4（果蝇）、GCF_000001735.4（拟南芥 TAIR10.1）、GCF_000146045.2（酵母）、GCF_016699485.2（鸡）、GCF_034140825.1（水稻）、assembly_summary_refseq.txt（233MB） |
| `data/m26_reference_pair_20260926/` | Araport11 GFF、TAIR10 FASTA（M26 参考源配对） |
| `data/m27_grch38_chr22/` | 人 chr22 基因组 + GENCODE v49（M27 等位完整性） |
| `data/m25_setaria/GCF_005286985.2/genome.fa` | **Setaria viridis 基因组 FASTA**，402MB，2026-08-24 17:35 |
| `data/m25_setaria_viridis/` | `GCF_005286985.2_Setaria_viridis_v4.0_genomic.fna.gz`（126MB，2026-08-24 17:18 下载） |
| `data/raw` `data/interim` `data/processed` | 仅 `.gitkeep` |

划分：M1 阶段划分文件为 `data/m1_screen/<species>/split_{train,val,test}.jsonl`（每行 `{"id": 染色体}`；拟南芥 train 含 NC_000932.1、NC_003070.9、NC_003071.7 等，test 为 NC_003075.7）。M28 的 TRAIN/DEV 划分记录在 `outputs/M28-DATA-MANIFEST-R2`、`outputs/M28-LABEL-POLICY-R4` 及 reports 中：从 LR 控制报告读到的六条原 TRAIN 染色体为 AT 的 NC_003070.9、NC_003071.7、NC_003076.8 与水稻的 NC_089043.1、NC_089044.1、NC_089046.1（22,688 主参考）；DEV 为 AT 与水稻各一染色体共 7,728 条 primary 参考链（来自 abstract 工作稿）。具体划分细节见各 protocol，本盘点未逐一核对。

Setaria 访问记录（只读证据）：
- 磁盘上**只有基因组 FASTA，没有任何 Setaria 注释（GFF/GTF）**；`outputs/` 下没有名称含 setaria 的目录；文件 atime 等于创建时间（8 月 24 日），没有之后被读取的文件系统痕迹（atime 在 beegfs 上可能不可靠，仅作旁证）。
- docs/28、docs/15、docs/research_roadmap、M28 plan、reports/M28 各 protocol 中多处明确记录 "Setaria 未访问/保持封存/注释禁运完好"（例如 M27 各阶段 "0训练、0mutant推理、Setaria未访问"，M28 full-grammar LR 控制计划 "不访问 sealed/held-out/Setaria/新物种/新面板"）。
- M25 原计划（docs/25）冻结 Setaria viridis（GCF_005286985.2，RS_2025_03）为盲测物种；M25R 因 0/5,625 元组可接受而无 Setaria 推理，评估从未发生。
- 结论：盲测未被消耗，注释未下载；基因组 FASTA 已存在（可用于将来的 raw FASTA 推理，但这不会破坏注释禁运）。

## 5. outputs/ 按里程碑（仅 ls，未 du）

按名称前缀计数（788 条）：
- BASE（5：ANNEVO/Helixer/Tiberius 烟测与 pilot）、FLOOR(1)、FP（42：SegmentNT 探针/CRF/FPLOSS/FUSION/FRAGFIX，5 个种子）、GENERANNO-NATIVE-PROBE（2）、M1（3）、SCREENREF（20）、SETUP（3）、SMOKE（8）、REANCHOR-CEILING-ANNEVO-M7、TB-UNFREEZE-BACKBONE-M9-DEEP
- M8（11）、M9（13，UNFREEZE-L0 到 L12）、M10（7）、M11（5）、M12A/B/C（4/5/3，含同 panel 基线 ANNEVO/Helixer/Tiberius）、M13（2）、M14（2）、M15（3）、M16（2）、M17（4，同 panel 泛化基线）、M18（6）、**M19（3：GENERANNO-1P2B-RAWCAL-CLEANPLANTS 与 -s0、-s1）**、M20（3）、M21（7，CRF）、M22（4）、M23（1，NTV2 clean transfer）、M24（1）
- M25（4：结构头，含 SMOKE）、**M25R（18）**（含 DEV-REDECODE 各 R3/R4 与 failed attempts、E1 phase/block/run-graph、R4 group-phase）
- M26（11：ACTUAL-GRAPH-SUPPORT-R3、ARAPORT-INPUT-R6、ENDPOINT-BUDGET-R8、PARTIAL-SEMANTICS-R4、REFERENCE-BACKGROUND-R2、REFERENCE-PAIR-R5/R7 及 repair、SAME-SCOPE-MECHANISM-R1）
- M27（6：ALLELIC-PREP-R1/R2、NATIVE-CONFIRM-R1、PAIRED-NATIVE-R1、SINGLE-EDGE-R1、WT-REPLAY-R1）
- **M28（573）**：A 系列 A、A2 到 A44（共约 250 个目录）、ASSEMBLY/B1/B2/C0/COMMON-RULER/CORE/DATA-MANIFEST/LABEL-POLICY/GLM、CARBON-A 系列（14）、NATIVE-*（本机原生基线）、**FULL-GRAMMAR 系列（247 目录）**：JOINT24、FRESH-JOINT24、UNIFORM-COARSE、CODON、COARSE-CHANNEL、GENE-POSTERIOR、E-M-ENERGY-POOL、MINIMUM-CENTERED、JOINT-LORA（及 LR1E5）、DENSE-QV（及 LR3E4）、QV-LR-FOUR-CELL 等。最新目录 `M28-FULL-GRAMMAR-QV-LR-FOUR-CELL`（2026-10-13 18:44）。
- 非里程碑项：`fp_segnt_logs`、`screen_ref_logs`、`generanno_dl.log`、`_mounttest.txt`、`SMOKE-*`。

未发现 AUGUSTUS 的任何输出目录，`which augustus` 为空；该基线在远端没有运行记录。基线（Helixer/Tiberius/ANNEVO）输出在 `outputs/BASE-*`、`outputs/M12B-SAMEPANEL-BASELINES-{HELIXER,TIBERIUS,ANNEVO,EXTERNAL}`、`outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES-{...}`、`outputs/M28-NATIVE-BASELINES`。M17 ANNEVO 目录含 FAILED_SPECIES.tsv/UNAVAILABLE.tsv（部分物种失败/不可用）。

## 6. 关键模型资产

### 6.1 M19 最佳 checkpoint
- 实验：`outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0`、`-s1`（GENERanno 1.2B cds-annotator-preview + LoRA r=8 alpha=16，q/k/v/o，fp_aware 损失，best_val_macro_f1 0.8121，seed 0，window 6144，bf16）。
- **没有保存 checkpoint**：该目录下无 checkpoints/、无 *.pt/*.safetensors/adapter 文件。脚本只保存 raw_scores（`raw_scores/*.npz`，约 110 到 180MB 每个）、predictions GFF、calibration、metrics。M10 到 M24 的 LoRA 系列也未发现任何 .pt/.safetensors/adapter（`find outputs/M1[0-9]* outputs/M2[0-4]* -maxdepth 3` 无结果）。因此 M19 模型无法重新加载，只能复用已保存的 raw score 与预测。此为重要的复现缺口。

### 6.2 M25R checkpoint
- `outputs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/checkpoints/epoch_{1,2,3}.pt`，每个 18,311,932 字节（约 17.5MB，仅结构头，不含骨干）；M25 初版 `outputs/M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/checkpoints/epoch_{1,2,3}.pt` 同尺寸（M25 因 4/1536 窗口全假掩码导致非有限损失，已判无效）。
- M25R 训练：1,536 TRAIN 窗口、3 个 epoch（约 12,900 秒/epoch）、无 selection；0/5,625 验证元组可接受，M25R 固定路线 NO-GO。

### 6.3 M28 当前最佳模型/LoRA checkpoint
M28 的当前竞争线为 "full-grammar"：GENERanno 1.2B（revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345`）+ 联合 24576bp 编码器适配 + 角色头 + 全 DNA 语法解码。以下是原 TRAIN 完整链 F1（pooled，两物种 AT+水稻，22,688 主参考）：

| 变体 | fit job | 路径（`outputs/`） | pooled F1 | checkpoint 大小 |
|---|---|---|---|---|
| E：rank4 LoRA q/v，native lr 3e-4（"lora_high"） | 13499817 | `M28-FULL-GRAMMAR-JOINT-LORA-FIT/step_009216.pt`（另 step_006144、step_007680） | **72.958%**（原 TRAIN 最好） | 5,668,948 B（约 5.4MB） |
| M：全秩 dense q/v，native lr 1e-5（"dense_low"） | 13510324 | `M28-FULL-GRAMMAR-DENSE-QV-FIT/step_009216.pt` | 72.200% | 230,922,900 B（约 220MB） |
| LoRA rank4 低 lr 1e-5 | 13545578 | `M28-FULL-GRAMMAR-JOINT-LORA-LR1E5-FIT/step_009216.pt` | 65.929% | 5,669,012 B |
| dense 高 lr 3e-4 | 13545584 | `M28-FULL-GRAMMAR-DENSE-QV-LR3E4-FIT/step_009216.pt` | 58.982% | 230,922,900 B |
| fresh-joint24（冻结 DEV 家族，DEV pooled 70.107%） | 13494895 | `M28-FULL-GRAMMAR-FRESH-JOINT24-FIT/step_009216.pt` | DEV 70.107% | 4,419,887 B（仅头） |

- 对照基线（同口径）：ANNEVO 75.542%、Tiberius 约 72.65%（TRAIN pooled，M 比其低 0.456pp）、Helixer 约 63.6%；Carbon-A（GENERanno 原生管线）filtered/unfiltered 78.172%/76.063%。因此**无任何 M28 变体超过 ANNEVO 或 Carbon-A**；LR 四格结果（10-14）未产生更好的候选，也没有新 DEV 候选。
- 冻结 DEV 家族（10-13 摘要）：joint 66.623%、fresh-joint 70.107%、uniform-coarse 69.795%、codon 64.579%，均低于 Tiberius 71.805% 与 ANNEVO 75.703%；高于 Helixer 62.851%。
- 其他 M28 checkpoint：JOINT24-FIT `step_009216.pt` 4.42MB（13453777）、JOINT24-NATIVE-PRIOR 4.42MB、CODON-FIT（13497914）4.42MB、COARSE-CHANNEL-FIRST-FIT（4608 步）4.42MB、COMPLETE-CORE-DRAWS/FRESH-DRAWS 4.42MB、B2-FIRST-FIT 2.2MB。`M28-A*` 的 query-head 系列（A44 到 10-07）也有 -FIT 目录（未逐个量大小）。
- 注意 "最佳"依赖口径：E 是原 TRAIN 上单 fit、单 seed 的最高；不是 DEV/独立泛化结论。所有 M28 checkpoint 均只含可训练部分（头 + adapter），需要基模 + `step_*.pt` 一起加载。
- 当前没有任何 M28 作业在跑或排队（squeue 空）。

### 6.4 基模权重缓存
- HF 缓存：`HF_HOME=/home/users/j/jwang/.cache/huggingface`（sbatch 以 HF_HUB_OFFLINE=1 使用）。`hub/` 下有 GENERanno-eukaryote-1.2b-cds-annotator-preview（4.3G）、GENERanno-eukaryote-0.5b-base、GENERanno-prokaryote-0.5b、NTv3 全系（8M/100M/650M，含 8kb 版）、NT-v2 multi-species（50M/100M/250M/500M）、segment_nt、gena-lm、caduceus、genatator-pipeline。
- Carbon-A 副本：`refs/carbon-a/f6789ee79e6ef397de333577ffec1cb734496620/`，`model.safetensors`（4,614,654,504 B）符号链接到上述 HF blob，附 annotation_pipeline（含 confidence_network/transition_network 小权重）；运行时在 `refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1`。
- 注意 `~/.cache/huggingface/token`、`stored_tokens` 存在（未读取内容）。`~/.cache/cache` 又是指向 ds4dh scratch 的符号链接。
- 其他基线权重：`refs/weights/helixer-2025`（fungi、invertebrate、land_plant、vertebrate）、`refs/weights/tiberius-m27-nosm-v2`、`refs/containers/helixer-docker_latest.sif`。

### 6.5 环境（`~/.conda/envs/`）
- **generanno**（M28 与大部分实验的实际环境，604 个 sbatch）、**annevo**（12 个）、**coding-rna**（cluster_config 默认，仅 1 个 sbatch 使用，6 月创建）。
- 其他无关/历史环境：bench、benchmark_core、caduceus_isolated、dnabert2_env、dna_consensus、dna_reconstruction、edta_native、ElementVec、quant、repeat、repet_*、rna、RNA、RNA-FM、te_benchmark、te_gap_tiberius_c、vector、video_enhance。
- 典型 M28 sbatch 头：`--partition=private-teodoro-gpu --gres=gpu:nvidia_geforce_rtx_3090:1 --cpus-per-task=2 --mem=16G --time=20:00:00 --exclude=gpu023,...,gpu043`，`source .../conda.sh; conda activate generanno`，`HF_HOME=.../huggingface HF_HUB_OFFLINE=1`。

## 7. 实验室 GPU 服务器（别名 `lab`，主机名 `ont`）
- 从 baobab 经 `ssh lab` 可连通（`~/.ssh/config` 中 `Host lab` -> 10.194.96.49，用户 ji）。
- 主机 `ont`：32 CPU、125GB 内存（可用 116GB）、根盘 7.0T（已用 4.1T，62%）。
- GPU：1 张 NVIDIA GeForce RTX 4090，查询时已用 14,026 MiB / 24,564 MiB（有他人任务）；`who` 显示用户 zdobnov 在线。
- 家目录内容：`ab-initio-TE envs hermes mamba_pkgs_fallback micromamba quant`。**没有 `~/coding-rna`**，没有本项目副本。仅可作为轻量/备用 GPU，需自行同步代码和数据，且 24GB 显存仅 ~10GB 空闲。

## 8. Credits 与存储

`ug_billing_synthesis`（PI teodoro，2025-12-01 到 2026-11-30）：
- Private partition billing 3,420,780；Free billing 100,000；**Total available 3,520,780**；**Consumed 3,373,517**；**Balance 147,263**（约 4.2%）。
- 私有分区分项：private_drim 309,228；private_teodoro 3,111,552。
- `ug_slurm_usage_per_user --user jwang --start 2026-10-01 --aggregate`（到 10-16 08:59）：账户 teodoro，billing **101,800**（TRES 小时，约 0.10M）。
- 粗略估算：10 月以来约 6.6k/天，按此速度余额约 22 天用尽（约 11 月上旬），早于 11-30 到期。余额是 PI 级（teodoro 组共享），不只是 jwang；此估算为外推，非确定。计费语义（私有节点是否实际扣费）未核实，数字按工具输出原样记录。

存储：
- 项目真实位置 `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna`（scratch，通常无备份），scratch 1005T 总量/847T 已用/158T 可用（85%）。
- `/home`（beegfs_home）138T，已用 135T，可用 3.3T（98%）。
- NFS quota 输出：`baobab-nfs` 29,903G / limit 51,200G；`baobab-nfs-fiteo` 39,362G / limit 40,960G（约 96%，接近限额）。这两个挂载与 coding-rna 的对应关系未确认，仅作记录。
- outputs/ 约 730G 为题述数值，**未测量**（按要求不对其整体 du）；refs/ 33G 为实测。
- 没有 `beegfs-ctl`，无法用其查询配额。

## 9. 风险与待办（仅记录，未执行）
1. 远端 924 条未提交变更，checkpoint 与 raw 输出只在 scratch；建议确认 GitHub 与集群同步状态，并对关键 `step_009216.pt`（E 与 M）、`refs/carbon-a` 等做异地备份（需用户批准）。
2. M10 到 M24（含 M19）无 checkpoint，只剩 raw scores/预测。
3. credits 余额偏低，约 3 周用量。
4. cluster_config.yaml 的 `conda_env: coding-rna` 与实际使用的 `generanno` 不一致；README 末尾"训练已开始"的状态落后。
5. `/home` 98%，若 HF 缓存/conda 继续增长可能影响登录。
6. Setaria：基因组 FASTA 在盘，注释不在，无访问记录；继续保持封存。
