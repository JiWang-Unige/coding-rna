# 01 代码与工程资产盘点 + 发布结构提案

日期：2026-10-16　范围：`/Users/jiwang/Desktop/coding RNA`（只读盘点，未修改/移动/删除任何已有文件，未提交 git）
配套文件：`01_path_mapping.tsv`（251 条规则，覆盖 1,874 个文件，无遗漏；列：old_path_or_glob / new_path / category / note）

---

## 0. 先读这一节：三个会改变后续工作的发现

1. **本地工作树不是自包含的，最终管线现在无法在本机 import。**
   - 本地 `main` 停在 `1d2593c`（M26），`origin/main` 比它多 39 个提交（2026-09-29 至 10-01，`d77d225`），里面有 M27 与「旧 M28（R4/B2/C0 chain-null）」代码；本地没有任何 M28 提交。
   - 当前 `src/m28/`（58 个 .py）、`tests/` 中 135 个、`sbatch/` 中 505 个、`scripts/` 中 361 个**全部是未跟踪（untracked）文件**，从未进入 git 历史。
   - 最终管线（full-grammar）在 import 时依赖 **5 个只存在于 `origin/main` 的 `src/m28` 模块**：`core.py`（`M28Core`，366,191 参数的角色头——Methods 里写的 “`src/m28/core.py`” 正是它）、`training.py`（仅 `one_hot`）、`export.py`（GFF3 写出）、`c0_decode.py`（仅 `genomic_chain`）、`candidates.py`（`Budget/ORFIndex/STOP`）；另有 4 个 M28 脚本 `common_ruler.py / manifest.py / census.py / feature_smoke.py` 和 2 个 M26 脚本 `reference_pair.py / semantic_reference_pair.py` 只在 origin/main。
   - 另有 **3 个同级模块在本地与 origin/main 与 git 历史中都不存在**：`query_training_census.py`（`prepare_full_grammar_fresh_draws.py` 的数据准备依赖）、`annotate_query.py`、`cache_full_train.py`（`annotate_full_grammar_coarse_channel.py` 即推理入口依赖）。它们大概率只在 Baobab 上。A 系列脚本还缺 `run_chain_quality/evaluate_pilot/run_query*` 等十余个同级模块与 `src/m28/query_{beam,decode,pointer,model,collaboration,known_mean,...}`，因此 A 系列（旧方法）本地已基本不可运行。
   - 静态检查：151 个测试文件里 **117 个（77%）在本地因缺模块无法 import**（缺失次数：core 81、candidates 67、training 65、c0_decode 53、export 45、assembly 39、census/manifest 33）。
   - 动作建议：发布重构前**第一步**是在 Baobab 上以工作树为准，做一次完整快照提交（含 untracked），并把 origin/main 的 14 个 `src/m28` 文件与上述缺失同级脚本并入同一棵树，再开始整理。`labels.py` 本地版比 origin 新（链 key 修复），以本地为准。
2. **最终方法的核心逻辑一部分写在 `scripts/` 里，且大量「诊断脚本」被训练/推理脚本当库 import。** 例如 joint 24,576-base 编码（`J.potentials / J.raw_view`）在 `fit_full_grammar_joint24.py`；GLM 编码在 `preflight_full_grammar_*glm*.py`；训练循环（`starting_state / example / objective / update`）在 `fit_full_grammar_objective_contrast.py`；观测目标构造在 `diagnose_full_grammar_training_observations.py`。
3. **代码以「Slurm 作业 ID 门禁」而非软件接口组织。** 199 个 M28 脚本含形如 `job_id != '13434603'` 的硬编码作业号与 `STATUS/JOBID/summary.json` 状态文件检查；脚本之间靠这些前置输出串联。可复现，但不可作为发布软件。

---

## 1. 资产盘点（文件数 / 行数）

| 目录 | 文件数 | 说明 |
|---|---:|---|
| `src/` | 79（py 9,620 行） | `src/m28` 58 个（≈4.5k 行，全部 untracked，无 `__init__.py`，靠 `PYTHONPATH` 命名空间包）；`src/screen_anchor` 7（M1–M9）；`src/foundation_probe` 6（M1–M25 探针/结构头，`train_generanno_structural_heads.py` 1,046 行） |
| `scripts/` | 461（py 64,866 行） | 根目录 75 个 M1–M9 时代文件（含 24 个 sbatch、`.bak`）；`scripts/experiments/` 下 M11–M26 共 10 个目录 22 个脚本；**M28-METHOD-RESTART 360 个脚本**（另 `__pycache__`） |
| `sbatch/` | 553 | M10–M26 共 49 个；**M28 504 个** |
| `configs/` | 26 | 均为 M1–M25R 的 yaml；**M28 没有任何 config**（参数写在脚本 `CONFIG` 常量与 `reports/…/*_plan.md` 中） |
| `tests/` | 151（14,560 行） | `test_m28_*` 134；M25R 9；M26 3；其余 5 |
| `outputs/` | 173 | 全部为 M28 运行输出（`.gitignore` 中），102 json / 6 md / 3 jsonl，无 checkpoint（`*.pt` 为 0） |
| `reports/` | 408 | M28-METHOD-RESTART 222、M18-ORACLE 104、M25R 34+10、M26 11、其余少量 |
| `docs/` | 28 | 8 份 M28 稿件（manuscript 170 KB）、台账/路线图、`figures/` |

### 1.1 M28 脚本（360 个）按动词与分支

动词分布：diagnose 44、infer 33、preflight 29、check 29、verify 27、evaluate 24、fit 23、run 20、probe 18、annotate 17、compare 15、readout 11、profile 11，其余 ≤6。
按方法代际：**full-grammar 181**、Carbon-A 17、原生基线 ≈7、**A 系列（pre-full-grammar）≈155**。

### 1.2 M28 sbatch（504 个）

| 系列 | 数量 |
|---|---:|
| `M28-FULL-GRAMMAR-*`（full-grammar 全部作业） | 235 |
| `M28-A5 … A44`（40 个子系列：A42=30、A44/A43=15…） | 236 |
| `M28-CARBON-*` / `M28-NATIVE-*` | 17 / 7 |
| 其他（图、unitcheck） | 9 |

其中 149 个文件名含 CHECK(S)：质量保证靠 sbatch 作业而非 pytest/CI。

### 1.3 Git 状态
跟踪 402 个文件；`git status` 列出 659 条变更（含 M28 全部代码与大批 reports 未跟踪；`outputs/` 被 .gitignore 忽略）。已修改未提交：`docs/12_publication_strategy.md`、`reports/M20-SOTA-ERROR-ANALYSIS/report.md`、`scripts/experiments/M20-SOTA-ERROR-ANALYSIS/analyze_sota_errors.py`。仓库还含 `refs/codex/turn-diffs/*` 引用（Codex 产物）。

---

## 2. 「当前最终方法管线」判定

### 2.1 判定依据
- 以 `docs/M28_full_grammar_methods_working.md` 与 abstract 的方法描述为准：冻结 GENERanno-1.2b + 15 角色 BiGRU 头；完整基因集 DNA 自动机；训练=部分观测（细/粗标签+边缘化）sum-product；推理=整 contig max-only exon-budget skyline；最终报告模型 = fresh-joint（pooled F1 70.107%），DEV 固定家族 = fresh-joint / uniform-coarse(ε=2/3) / codon。
- 用 AST 构建了 `src / scripts / tests` 的 import 图（含相对导入与同级裸 import，并合并 origin/main 的文件），对 9 个入口做传递闭包；再用 sbatch 调用关系（340 个脚本被 sbatch 调用）核对：最终链路 sbatch 即 `M28-FULL-GRAMMAR-{FRESH-DRAWS-PREPARATION → COARSE-CHANNEL-FIRST-FIT → FRESH-JOINT24-(PREFLIGHT|FIT|VERIFY|TRAIN-RAW|TRAIN-EVAL) → UNIFORM-COARSE-* / CODON-* → JOINT-CONTROLS-DEV-(CHECKS|RAW|COMPARISON|FAMILY)}`。
- 结论：**M28 的 A5–A44 系列（236 sbatch、≈155 脚本、≈45 个 src 模块、≈80 个测试）与最终管线无 import 关系**，只是被 Methods 称为 “earlier finite-candidate reranking models”。其中 emission_lattice 变体见 §3 的例外。

### 2.2 最终管线入口（逐个）

| 角色 | 入口文件（现路径） | 备注 |
|---|---|---|
| 数据准备 ① 窗口 manifest | `scripts/experiments/M28-METHOD-RESTART/manifest.py` + `census.py`（**仅 origin/main**）；参考解析 `M26-SAME-SCOPE-MECHANISM/reference_pair.py / semantic_reference_pair.py`（**仅 origin/main**） | 24,576/12,288 连续窗口 |
| 数据准备 ② 抽样计划 | `prepare_full_grammar_fresh_draws.py`（依赖 **缺失** 的 `query_training_census.py`）、`check_full_grammar_fresh_draws_support.py` | seed-1 fresh 4,608 draws；complete-core 分支用 `prepare_full_grammar_complete_core_draws.py` |
| 训练入口 ① 基础拟合 | `fit_full_grammar_coarse_channel.py fit`（4,608 更新） | 作业 13434603，门禁常量 |
| 训练入口 ② 最终续训 | `fit_full_grammar_fresh_joint24.py fit --source-job`（joint 24,576 编码，再 4,608 更新→9,216） | 另 `fit_full_grammar_uniform_coarse.py`、`fit_full_grammar_codon.py`（DEV 家族另两个候选） |
| 训练库 | `fit_full_grammar_objective_contrast.py`（循环）、`fit_full_grammar_joint24.py`（joint 编码）、`preflight_full_grammar_coarse_channel_glm.py` / `preflight_full_grammar_glm_context.py`（GLM 编码）、`diagnose_full_grammar_training_observations.py`（观测目标）、`src/m28/{core,training,labels,full_grammar_coarse_channel,full_grammar_boundary_compiled,full_grammar_observed_partition,local_coarse_channel}.py` | |
| 推理/解码入口 | `annotate_full_grammar_coarse_channel.py`（模型加载→窗口编码→`decode_full_grammar_skyline`→GFF3）；逐 contig 驱动 `infer_full_grammar_original_train.py --task`；最终模型 `infer_full_grammar_fresh_joint24_train.py raw/evaluate`；DEV `infer_full_grammar_joint_controls_dev.py --branch --task` | 严格 checkpoint 适配器 `annotate_full_grammar_{joint24,fresh_draws,fresh_joint24,uniform_coarse,codon}.py` |
| 解码核心 | `src/m28/full_grammar_skyline.py`（生产）、`full_grammar_codon_skyline.py`（codon）、`full_grammar_viterbi.py`（参考；skyline 仍 import 其 `_selected_scores`） | |
| 评估入口 | `evaluate_full_grammar_joint_controls_dev.py --branch`、`evaluate_full_grammar_dev_baselines.py`、`evaluate_full_grammar_original_train.py`、`aggregate_full_grammar_joint_controls_dev.py`；评估核心 `common_ruler.py`（**仅 origin/main**）、`scripts/eval_gene_body_mask.py`、`scripts/eval_structure_diagnostic.py` | 同尺：CDS-assessable primary、完整 CDS 链匹配 |
| 基线 | `run_native_baseline.py`、`run_native_train_baseline.py`、`evaluate_native_train_baselines.py`；Carbon-A 17 脚本 | ANNEVO/Tiberius/Helixer/Carbon-A |

### 2.3 最终管线 `src/m28` 依赖闭包（传递，25 个模块）
- 直接属于方法：`labels`、`contig_emission_lattice`（仅窗口网格/拼接）、`full_grammar_context_view`、`full_grammar_boundary_messages`、`full_grammar_boundary_compiled`、`full_grammar_observed_partition`（+`_partition_reverse`、`_partition_compiled_reverse`）、`full_grammar_coarse_channel`、`local_coarse_channel`、`full_grammar_skyline`、`full_grammar_viterbi`、`full_grammar_codon`、`full_grammar_codon_skyline`；仅 origin：`core`、`training`、`export`、`c0_decode`、`candidates`。
- **遗留但仍被传递依赖（必须先抽取符号再归档）**：`emission_lattice`、`emission_lattice_compact`、`emission_lattice_phase_block`、`exon_end_bound`、`long_exon_phase`、`junction_chain`、`phase_evidence_chain`。例如 `observed_partition` 只用 `_CompactMargins` 和 `exon_event_range`，`contig_emission_lattice` 顶部 import 了 `decode_emission_lattice`（最终管线并不调用）。
- Methods 中另述但不在训练/解码主链上的分析模块：`full_grammar_gene_posterior`（CDS 事件边际）、`full_grammar_scaled_mass`、`gene_posterior_ranking`、`full_grammar_energy_pool`（E-M 集成）、`joint_encoder_{lora,dense_qv}`（q/v 适配对照）。归类为 final_pipeline 的「可选分析」。

### 2.4 测试覆盖（对 full-grammar 模块）
- 有单元测试：boundary_messages/compiled、observed_partition、partition(+reverse/compiled_reverse/indexed)、skyline、viterbi、codon(+skyline/energy/train)、scaled_mass、gene_posterior、energy_pool、joint_encoder_*，22 个文件 125 个 `def test_`。
- **无单元测试**：`core.py`(角色头)、`export.py`(GFF3)、`c0_decode.genomic_chain`、`full_grammar_coarse_channel.py`（观测似然，论文核心）、`local_coarse_channel.py`、`full_grammar_context_view.py`（nearest-center 归属）、`labels.py` 仅 1 个间接测试。训练/推理驱动层靠 sbatch “CHECKS” 作业。

---

## 3. 重复 / 变体文件与最终版本

| 家族 | 变体 | 最终/保留 | 说明 |
|---|---|---|---|
| emission_lattice | `emission_lattice`、`_batched`、`_cached_metadata`、`_compact`、`_compiled`、`_indexed`、`_phase_block`、`contig_` | **都不是最终方法**；仅 `emission_lattice`(decode)、`_compact`(`_CompactMargins`)、`_phase_block`、`contig_emission_lattice`(窗口网格) 因 import 链被保留 | `_batched/_cached_metadata/_compiled/_indexed` 是性能迭代（A42 系列），无最终管线引用 |
| 配分函数 | `full_grammar_partition`、`_indexed`、`_reverse`、`_compiled_reverse`、`boundary_messages`、`boundary_compiled`、`observed_partition` | **训练=`boundary_compiled`（+`boundary_messages` 的转移）；Zobs 掩膜=`observed_partition`**；`partition`/`partition_indexed` 为参考实现（仅测试） | Methods 明确写 `boundary_messages.py`(转移)、`boundary_compiled.py`(实现) |
| 解码器 | `full_grammar_viterbi`、`full_grammar_skyline`、`full_grammar_codon_skyline` | **生产=`skyline`；codon=`codon_skyline`；`viterbi`=参考/校验** | |
| 编码适配 | `joint_encoder_lora` / `_dense_qv` / `_qvo_lora` | Methods 记 q/v rank-4 与 full-rank；`qvo` 未写入稿件→归档 | |
| per-branch 脚本三联 | `fit_/annotate_/infer_/verify_/compare_full_grammar_{objective_contrast, native_cds_prior, joint24, fresh_draws, fresh_joint24, uniform_coarse, codon, joint_lora(_lr_low), dense_qv(_lr_high), joint_qvo_lora, complete_core_draws, minimum_centered, e_m_energy_pool}` | 最终=fresh_joint24 (+uniform_coarse,codon) | 变体间 diff 相似度 0.91–0.94（如 `fit_joint_lora` vs `_lr_low/_qvo_lora`），应合并为单一实现 + `--variant`/配置 |
| 旧 A 系列 | `chain_competition`、`chain_set`、`context_chain`、`natural_chain_set`、`null_anchored_*`、`local_*`、`query_*`、`structured_emission` 等 | 全部已废弃 | 多个已缺依赖，本地不可运行 |
| 其他 | `scripts/probe_generanno_native.py` 与 `.bak`；`CLAUDE.md.framework-upgrade-backup-20260614-214748` | 保留 .py，删 .bak/备份 | |
| 注意 | `*_compiled`、`*_indexed`、`*_batched` 的名称含义并不一致：`full_grammar_boundary_compiled` 是**最终版**，而 `emission_lattice_compiled` 是废弃版 | | 重构时按「家族」而非后缀判断 |

---

## 4. 建议的新目录结构（论文复现 + 对外发布）

```
cdsgram/                          # 仓库根（可 git 化；包名暂定 cdsgram，可改）
├── pyproject.toml  README.md  LICENSE  CITATION.cff  environment.yml  CHANGELOG.md
├── src/cdsgram/                  # pip install .
│   ├── data/      fasta/N 映射、labels(15 角色/细粗标签)、windows(24576/12288 网格+nearest-center 拼接)、
│   │              manifest、observations(细/粗目标)、gff 输出(as_transcript/write_gff3)、encoding(one_hot)
│   ├── model/     role_head(M28Core C0)、encoder(GENERanno 加载+joint24 编码)、checkpoint(各 variant 严格加载)、
│   │              codon、lora(可选)
│   ├── decoder/   automaton(DNA 自动机)、partition(sum-product)、observed(Zobs 掩膜)、observation(粗通道)、
│   │              skyline、codon_skyline、viterbi_ref、posterior/scaled_mass(可选)、coords/orf
│   ├── train/     loop、fit(基础)、continuation(fresh-joint)、schedules(抽样计划)、objective(局部损失)、variants/
│   ├── eval/      ruler(同尺评估)、reference(参考选择)、family(DEV 家族聚合)、posterior_ranking、gene_body
│   └── cli/       `cdsgram train | annotate | evaluate | prepare-data`（console_scripts）
├── tests/{unit,integration}      # 仅保留最终管线测试；CI 用 CPU 小样例（artificial DNA）
├── slurm/m28/{00_data,01_fit_base,02_fit_*,03_infer,04_eval,05_dev}/   # 最终链路 sbatch 模板，路径/环境参数化
├── experiments/                  # 按里程碑归档，保留证据，不进包
│   ├── M1-M9/ M10-M19/ M20-M24/ M25R/ M26/ M27/
│   └── M28/{full_grammar/{branches,diagnostics,gates,dev_families,infer,reports,sbatch,tests}, gates}
├── benchmarks/{native_baselines,carbon_a,panels}/   # ANNEVO/Tiberius/Helixer/Carbon-A 运行器与同尺面板
├── paper/{manuscript,figures{,/scripts},tables,strategy,supplement}/
├── results/                      # gitignore：outputs/ 原样；紧凑 summary.json 提升进 experiments/*/results/
├── archive/                      # 已废弃：m28_a_series(+sbatch,tests,reports)、m28_legacy_lattice、m1_m9_*、auto_research
└── docs/{project_history,review}/
```
设计要点：(a) `src/cdsgram` 只放最终管线 + 少数可选分析；(b) 旧 A 系列整体入 `archive/`，不删；(c) 遗留模块（emission_lattice 等）**先抽取符号**到 `decoder/orf.py`、`decoder/exon_bounds.py` 再归档；(d) 每个 variant 脚本三联合并成「一份实现 + 配置」；(e) 论文所引用的 `reports/` 证据随里程碑迁入 `experiments/<M>/reports`，稿件中的相对链接需在迁移后批量重写（现 113 个相对链接均有效，是迁移回归测试的好基准）。

### 4.1 映射表摘要（完整见 `01_path_mapping.tsv`）
| 旧 | 新 | 类别 | 文件数 |
|---|---|---|---:|
| `src/m28`（最终闭包 25 + 可选分析 ≈8） | `src/cdsgram/{data,model,decoder,train,eval}` | final_pipeline | 22 本地 + 5 仅 origin |
| `src/m28` 其余（A 系列、性能变体） | `archive/m28_a_series_modules/`、`archive/m28_legacy_lattice/` | deprecated | 47 |
| `scripts/experiments/M28-…` 最终入口与库 | `src/cdsgram/{cli,train,model,eval}`、`experiments/M28/full_grammar/{gates,infer}` | final_pipeline | 41（含 scripts 根 2 个评估库） |
| `scripts/experiments/M28-…/*full_grammar_*`（对照/诊断/检查） | `experiments/M28/full_grammar/{branches,diagnostics,dev_families,recovery}/` | experiment_archive | ≈140 |
| `scripts/experiments/M28-…` A 系列 | `archive/m28_a_series/` | deprecated | ≈155 |
| `*carbon_a*`、`*native_*` | `benchmarks/carbon_a/`、`benchmarks/native_baselines/` | infra | 28 |
| `sbatch/M28-FULL-GRAMMAR-*`（最终链路 28） | `slurm/m28/<stage>/` | final_pipeline | 28 |
| `sbatch/M28-FULL-GRAMMAR-*`（其余 ≈207） | `experiments/M28/full_grammar/sbatch/` | experiment_archive | ≈207 |
| `sbatch/M28-A*` 等 | `archive/m28_a_series/sbatch/` | deprecated | 241 |
| `sbatch/M10–M26-*` | `experiments/<M>/sbatch/` | experiment_archive | 49 |
| `tests/test_m28_full_grammar_*` 等 | `tests/unit/{decoder,model,train,eval}` | final_pipeline | 28 |
| `tests/test_m28_*`（其余） | `archive/m28_a_series/tests/` | deprecated | 101 |
| `configs/*` | `experiments/<M>/config.yaml`、`benchmarks/panels/` | experiment_archive / infra | 26 |
| `docs/M28_full_grammar_*` | `paper/manuscript/` | doc | 8 |
| `reports/*` | `experiments/<M>/reports/`、`paper/…`、`archive/m28_a_series/reports/` | doc / deprecated | 408 |
| `outputs/*` | `results/`（gitignore） | infra | 173 |

（类别统计，按文件：final_pipeline 119、experiment_archive 533、deprecated 658、infra 228、doc 336。）

---

## 5. 发布前工程缺口

| # | 缺口 | 证据 | 建议 |
|---|---|---|---|
| 1 | **代码不在版本库里 / 本地缺失模块** | §0.1：M28 全部 untracked；5 个 src 模块+8 个脚本仅在 origin/main；3 个同级脚本两边都无；117/151 测试本地不可 import | 以 Baobab 为准做全量快照提交；合并 origin/main；补齐缺失三脚本；再重构 |
| 2 | **硬编码 Baobab 路径** | 553/553 sbatch 含 `cd /home/users/j/jwang/coding-rna`、`HF_HOME=/home/users/j/jwang/.cache/huggingface`、`source /opt/ebsofts/Mamba/…`、分区 `private-teodoro-gpu`、RTX3090 gres、`--exclude=gpu023…` 节点；3 个 py 含绝对路径（`describe_full_grammar_population_obstructions.py`、`qualify_carbon_a_runtime.py`、`probe_carbon_a_library_io.py`）；其余 py 用 `Path(__file__).parents[3]` 相对根（好）；日志写入 `outputs/M28-PILOT-R5/logs`（旧目录名遗留） | 参数化 `CDSGRAM_ROOT/HF_HOME/CONDA_ENV`；提供 `slurm/` 模板 + 非 Slurm 本地运行示例 |
| 3 | **对已归档 `ACTIVE_GOAL.json` 的依赖** | 11 个 `scripts/*` 与 24 个 `sbatch/*`（M10–M23）、3 个 configs 引用；`validate_goal.py` 及其测试；**M28 为 0**（最终管线不依赖） | 全部整体入 `archive/auto_research/`、`experiments/M1-M23`，不进入包；不恢复活动 goal |
| 4 | **作业号/状态文件门禁** | 199 个 M28 脚本含硬编码 job id；如 `starting_state()` 要求 `summary['job_id']=='13434603'`；`test ! -e "$out"`；`outputs/*/STATUS` | 发布版用 checkpoint 内容哈希/配置字典校验取代 job id；保留 `experiments/` 里的历史门禁仅作复现 |
| 5 | **无 CLI / 无包结构** | 脚本是 `python scripts/.../x.py fit --source-job N` / `--task N`；`src/m28` 无 `__init__.py`；靠 `PYTHONPATH=$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART` 与同级裸 import；`from src.m28...` | `pyproject.toml` + `console_scripts`；统一 `cdsgram annotate --fasta --ckpt --out.gff3`、`cdsgram train --config`、`cdsgram evaluate --pred --ref` |
| 6 | **缺 README / LICENSE / 环境文件 / 引用** | 根目录仅 AGENTS.md、CLAUDE.md 与目录；无 `pyproject/requirements/environment.yml`；环境只在 Methods 中口述（Python 3.10.18、PyTorch 2.5.1+cu121、Transformers 4.49.0；另 numba、SciPy 1.15.3）；sbatch 里 conda env 名 `generanno`(533)/`coding-rna`(5)/`annevo`(4) 并存 | 补 README、LICENSE（注意 GENERanno/Carbon-A/Tiberius/Helixer/ANNEVO 第三方许可）、锁定环境、`CITATION.cff`、模型权重（checkpoint 发布位置、HF revision `b0483c23…`） |
| 7 | **测试状态** | 本机无 torch/pytest（Python 3.9.6），未能运行；静态：34 可 import / 117 阻塞；full-grammar 125 个测试函数；核心模块 `core/export/coarse_channel/context_view/local_coarse_channel` 无单测；质量保证靠 149 个 CHECKS sbatch | 在 Baobab 的 `generanno` 环境实际跑一次基线 pytest 记录通过率；补 coarse_channel、context_view、GFF3 往返、CLI 冒烟（artificial DNA, CPU）测试；去掉对 `outputs/` 的依赖；加 CI |
| 8 | **配置缺失** | M28 的 config 全在代码常量（`CONFIG/IDENTITY/PROTOCOL` 路径指 `reports/…_plan.md`）；`configs/` 26 个全是 M25R 之前 | 抽成 `configs/*.yaml`（学习率 3e-4、AdamW、ε=0.1、窗口 24576/12288、seed、四块 6144 等） |
| 9 | **重复实现** | 变体脚本 diff 相似度 0.91–0.94；fit/annotate/infer/verify/compare 逐分支三联 | 合并，保留 `--variant` |
| 10 | **科学/发布边界**（非工程但影响发布） | 稿件已说明：主干为 annotation-pretrained checkpoint，TRAIN 与 DEV 均为历史复用开发集，未达独立泛化/超越 ANNEVO 75.7%、Tiberius 71.8% | 发布软件前确认是否有 sealed 数据未开放的约束（AGENTS.md：不自动放行封存数据/对外发布） |
| 11 | 次要 | `CLAUDE.md.framework-upgrade-backup-*`、`probe_generanno_native.py.bak`、`__pycache__`；稿件里指向不存在 `src/m28/core.py` 的说明；`outputs/M28-PILOT-R5/logs` 目录名 | 清理 |

---

## 6. 迁移执行顺序（供下一步参考，本报告未执行任何迁移）

1. Baobab 上全量快照提交 + 合并 origin/main；锁定哈希。
2. 在 `generanno` 环境跑现有 pytest，记录基线通过/失败（确认 core 等缺失模块恢复后 151 测试状况）。
3. 抽取遗留依赖符号（`ORFIndex/Budget/STOP`、`_CompactMargins`、`exon_event_range`、`candidate_roles`、`decode_emission_lattice` 去除）→ 新 `cdsgram` 包，保持数值等价（用既有 `tests/test_m28_full_grammar_*` 做回归）。
4. 将 `fit_*joint24`、`preflight_*glm*`、`objective_contrast` 的库函数迁入 `model/encoder.py`、`train/loop.py`。
5. 写 CLI，并用已完成的 fresh-joint TRAIN 六染色体结果做端到端回归（逐位 GFF 对比）。
6. 按 TSV 迁移归档目录；批量重写稿件与报告内的相对链接，用 113 个现有链接校验。
7. 补 README/LICENSE/environment/CI，参数化 Slurm。
