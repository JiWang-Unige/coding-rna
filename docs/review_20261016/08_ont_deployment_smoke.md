# 08 M28 full-grammar 管线部署到 ont 与冒烟复现

日期：2026-10-16。执行范围：只读使用 Baobab（未提交任何 Slurm 作业、未修改 Baobab 任何文件、未 git commit/push），所有写入都在 ont `~/coding-rna`。未推送 Setaria 或任何 sealed/held-out 数据。

## 0. 结论

1. **ont 可以作为 M28 推理/评估的主力算力。** fresh-joint final9216 检查点在 ont 上对 DEV 拟南芥 2.0 Mb 区段和水稻 4.9 Mb 区段做了完整的"编码 → 整 contig exon-budget skyline 解码 → GFF3"，与 Baobab 上已有的整染色体 DEV 预测逐基因比较：拟南芥区段 541/541 条内部链完全一致；水稻区段 180/181 一致，1 条因近似平局的起点差异（Δ=245 bp，分差 1.4）不同。两段上同一把尺（`common_ruler`）的 exact-chain TP/预测/参考数在 ont 与 Baobab 之间完全相同。
2. **吞吐 17.2–17.4 kbp/s，比 Baobab 3090（13.6–13.8 kbp/s）快约 26%**，且这是在另一用户占着 GPU（100% 利用率、约 11 GB）时测得的。推理峰值显存只有 3.0 GiB（nvidia-smi，含 CUDA 上下文），能和他人进程共存。
3. 测试：22 个 full-grammar 测试文件 140 passed / 0 failed；整个 `tests/` 991 passed、1 failed、1 collection error，两个异常都与 ont 无关（见 §3）。
4. 未验证的部分：**训练**（`fit_*`）没有在 ont 上跑过，显存与速度未知；整条 DEV 染色体、全 TRAIN 的真实运行未做（只给估算）；DEV 的 pooled F1 70.107 没有在 ont 上重算（该评估脚本要求 Baobab 的 sacct 和三基线输出，见 §5）。

## 1. 部署内容与路径

ont 用户 `ji`，项目副本 `~/coding-rna`（非 git 仓库；来源是 Baobab `/home/users/j/jwang/coding-rna` 的分支 `snapshot-20261016`，HEAD `4d1a21b`，工作树干净）。

| 内容 | ont 路径 | 大小 | 来源/说明 |
|---|---|---|---|
| 代码 | `src/`（含 `src/m28` 全部 58 个模块）、`scripts/`、`tests/`、`configs/` | 约 6.4 MB | `ont_sync.sh push`；因为取自 Baobab 快照，`core/training/export/c0_decode/candidates`、`common_ruler/manifest/census`、`query_training_census/annotate_query/cache_full_train` 等 01 号报告里本地缺失的模块都在，import 链完整 |
| 最终检查点 | `outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FIT/{step_009216.pt, summary.json, STATUS, JOBID}` | 4.4 MB | 只推最终 step（头部 366,191 参数；冻结的 GENERanno 权重不在其中） |
| ANNEVO HMM 源 | `refs/repos/annevo-2026/src/` | 285 KB | `src/m28/c0_decode.py` 的默认 HMM 路径，仅测试用到（见 §3） |
| Baobab 参考预测 | `outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-RAW-R1/FRESH_JOINT24576/{NC_003074.8,NC_089041.1}/{prediction/full_grammar.gff3, summary.json}` | 10 MB | 用于逐基因比较；DEV 预测，不含 sealed |
| 冒烟脚本（ont 上新写） | `scripts/ont_smoke/{make_region.py, infer_fasta.py, run_ont.sh, compare_region.py, score_stats.py}` | 数 KB | 见 §6、附录 |
| 冒烟产物 | `outputs/ONT-SMOKE/`（`region_{A,B}_input.fasta`、`regionA_fresh_joint/`、`regionB_rice_fresh_joint/`、`region{A,B}_compare.json`、`NC_003074.8_whole_input.fasta`） | 约 50 MB | |
| 日志 | `logs/`（`pytest-*.{out,err}`、`regionA_fresh_joint-20261016-115155.*`、`regionB_rice_fresh_joint-20261016-115524.*`） | | |

ont 上原有、未改动：环境 `~/envs/generanno`（Python 3.10，torch 2.5.1+cu124，transformers 4.49.0，pytest 9.1.1，numba 0.62.1，numpy 2.2.6）、`data/m1_screen/{arabidopsis_thaliana,oryza_sativa}`（含此前已存在的 `split_*.jsonl`；本次没有新推数据）、HF 缓存 `models--GenerTeam--GENERanno-eukaryote-1.2b-cds-annotator-preview`（含代码要求的修订 `b0483c23…`，与 Baobab 缓存布局相同，离线加载成功）。

没有推送：`data/m25_setaria*`、`reports/`、`sbatch/`、`docs/`、其余 `outputs/`。推理与测试都不需要 `reports/`（检查点里的 `PROTOCOL` 只是字符串比对）。

## 2. 管线入口与 DEV 推理所需输入（从 Baobab sbatch 与脚本读出）

`sbatch/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-RAW.sbatch`：`private-teodoro-gpu`，3090×1，2 CPU，16 GB，2 h 上限，数组 0-1（0=拟南芥 NC_003074.8，1=水稻 NC_089041.1），参数是分支名。环境：conda `generanno`；`HF_HOME=~/.cache/huggingface`、`HF_HUB_OFFLINE=1`、`TRANSFORMERS_OFFLINE=1`、`NUMBA_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=2`、`TOKENIZERS_PARALLELISM=false`、`PYTHONPATH=$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART`。调用 `infer_full_grammar_joint_controls_dev.py --branch FRESH_JOINT24576 --task N`。

该脚本做三件事：(1) `require_train_closure()`——校验 30 个 TRAIN Slurm 作业的目录、`summary.json`、`STATUS`、`JOBID`，并用 `sacct` 检查全部 COMPLETED；(2) 把 DEV 染色体写成大写、80 列的 `input.fasta`；(3) 调用 `annotate_full_grammar_fresh_joint24.main(fasta, checkpoint, out/'prediction')`。第 (3) 步才是科学计算：加载检查点（`validate_checkpoint` 严格核对 step=9216、epoch=6、`PROTOCOL`、编码器修订和全部 `IDENTITY` 配置项）、冻结 GENERanno 以 24,576 bp 整窗 BF16 flash-SDPA 编码、两条链分别做 stitch + `decode_full_grammar_skyline`、写 GFF3 与 `summary.json`。

## 3. 测试结果（ont CPU，`CUDA_VISIBLE_DEVICES=` 置空）

| 范围 | 结果 | 日志 |
|---|---|---|
| `tests/test_m28_full_grammar_*.py`（22 个文件） | **140 passed，468 subtests passed，0 failed**，20.5 s | `logs/pytest-fullgrammar-20261016-115102.out` |
| 管线支撑模块 16 个文件（core、labels、label_alignment、export、c0_decode、c0_gain、common_ruler、manifest、census、contig_emission_lattice、joint_encoder_lora/dense_qv、candidates、energy_pool_gene_posterior、gene_posterior_ranking、minimum_centered_energy_pool） | **62 passed，6 subtests passed**，7.7 s。第一次跑 8 个失败，全部是 `FileNotFoundError: refs/repos/annevo-2026/src/HMM.py`（旧原生 ANNEVO 桥接测试）；推上 285 KB 的 ANNEVO `src/` 后消失 | `logs/pytest-pipeline-support-*.out` |
| 整个 `tests/`（`--continue-on-collection-errors`） | 991 passed，860 subtests passed，**1 failed，1 error**，68.6 s | `logs/pytest-all-20261016-115848.out` |

整套里的两个异常都不是 ont 问题：
- `tests/test_m27_single_edge.py` collection 报 `cannot import name 'HMM' from 'src'`——`src/HMM.py` 在 Baobab 上同样不存在，M27 遗留测试，与 M28 无关。
- `tests/test_m28_emission_lattice_reuse.py::…saved_full_window_chains…` 读取 `outputs/M28-A42-EMISSION-LATTICE-LONG-CHECKS/case_0/STATUS`，是旧 A42 系列的输出依赖，未推到 ont；不在最终管线上。

## 4. 冒烟复现

### 4.1 设计
- 区段对齐：区段起点取 12,288 的整数倍、长度取 24,576+k·12,288，使区段内每个编码窗口的起点与 Baobab 整染色体运行的窗口完全重合（编码输入逐窗口相同）。
- 区段 A：拟南芥 NC_003074.8 `[4,915,200, 6,905,856)`，1,990,656 bp，161 窗/链，Baobab 参考内含 570 条链。
- 区段 B：水稻 NC_089041.1 `[7,372,800, 12,312,576)`，4,939,776 bp，401 窗/链。
- 比较方式：把区段预测平移回染色体坐标，与 Baobab 整染色体 GFF3（作业 13520474/13520473，同一检查点 `step_009216.pt`）在区段内逐基因比较"strand + 有序 CDS 区间"。内部区只取离区段两端各 ≥50 kb 的基因，以排除边界效应；同时报告边距 0 的结果。再用未修改的 `common_ruler.select_reference/…S.key`（M26 CDS-assessable primary 口径）对两份预测在同一内部区算 exact-chain P/R/F1。
- 绕过的门禁：只绕过 `infer_full_grammar_joint_controls_dev.py` 的 Slurm job-id / sacct / TRAIN 闭包检查，薄 wrapper 直接调用它最后调用的同一入口 `annotate_full_grammar_fresh_joint24.main`；科学代码零改动，检查点身份校验（`validate_checkpoint`）仍然生效。

### 4.2 一致性

| | 区段 A 拟南芥 | 区段 B 水稻 |
|---|---|---|
| ont 区段内全部链 / Baobab 整染色体链 | 570 / 5,233 | 187 / 2,310 |
| 内部区（边距 50 kb）ont / Baobab 链数 | 541 / 541 | 181 / 181 |
| 完全相同（strand + 全部 CDS 区间） | **541（100%）** | **180（99.45%）** |
| 仅 ont / 仅 Baobab | 0 / 0 | 1 / 1 |
| 边距 0 | ont 570 / Baobab 568，568 全部相同，ont 多 2 条（区段两端的边缘链） | 187/187，186 相同 |
| 同一把尺 exact-chain TP/预测/参考 | 423/541/551，F1 77.47%（两边完全相同） | 98/181/187，F1 53.26%（两边完全相同） |
| 相同链的分数差（链增益，|Δ|） | 中位 0.25，均值 0.375，最大 3.06（相对均值 0.21%） | 中位 0.22，均值 0.37，最大 3.48（相对均值 0.7%，最大 61% 出现在增益极小的链上） |

唯一的不一致（水稻）：反义链基因 `11,425,152–11,432,157`（ont，分数 118.97）对 Baobab 的 `11,425,397–11,432,157`（117.56），同为 4 个 CDS、终点相同，起点差 245 bp，分差 1.4。这是两个近似平局的起点之间的翻转。

解释：链结构一致，连续分数有 0.2% 量级的漂移，这符合 4090 与 3090 上 BF16 flash-attention 内核、FP32 累加顺序的数值差异；这不是代码或路径问题。与此一致的是，ont 与 Baobab 的分数差没有系统偏移（带符号均值 +0.08 / +0.05，ont 更高的链占 51% / 49%），而且 exact-chain 计数完全相同。要得到更强的"逐比特"结论，需要在同一型号 GPU 上重跑，ont 只有 4090，做不到。

注意：这是区段子集的复现，不等于 DEV pooled F1 70.107 被复现；区段 F1 本身（77.47%、53.26%）只是区段内的数字，不能和 DEV 整染色体数字比较。

### 4.3 整染色体输入已核对
用 `make_region.py … 0 0`（整染色体模式）在 ont 生成的 `NC_003074.8_whole_input.fasta` 与 Baobab 的 `…/FRESH_JOINT24576/NC_003074.8/input.fasta` md5 相同（`09856b3649c94f5a42e0801a1a6de5b9`），即整条 DEV 染色体可以直接在 ont 上运行，输入与 Baobab 逐字节一致。

## 5. 性能与显存

测量时 GPU 被另一用户的 `te_benchmark` 进程占用（`nvidia-smi`：约 11 GB、利用率 99–100%，全程如此），所以以下是"共享 GPU"下的数字，空闲时应不会更差。

| | ont RTX 4090（共享） | Baobab RTX 3090（作业 13520473/4，独占） |
|---|---|---|
| 区段 A（1.99 Mb）总耗时 | 116.0 s → **17.16 kbp/s** | — |
| 区段 B（4.94 Mb）总耗时 | 283.4 s → **17.43 kbp/s** | — |
| 编码耗时 / 窗口（24,576 bp，BF16） | 0.267 s（A），0.287 s（B） | 0.327 s（全拟南芥染色体 625 s/1,909 窗/链） |
| 解码（CPU，每链每 Mb） | 6.5–7.6 s | 10.6 s |
| 整染色体 | — | 拟南芥 1,744 s（13.6 kbp/s）；水稻 2,184 s（13.8 kbp/s）；sacct 29:44 / 37:02，billing 12/h |
| 峰值显存 | torch 峰值已分配 2.52 GB、已保留 2.70 GB；nvidia-smi 进程峰值 **3,044 MiB** | 同（2.52/2.70 GB） |
| 内存 | — | 5.4 GB（拟南芥）/ 7.1 GB（水稻）MaxRSS |

算力估计（取 17.3 kbp/s，线性外推；解码与编码都随长度线性）：

| 任务 | bp | ont 时长（每个模型/检查点） | Baobab 3090 对照 |
|---|---|---|---|
| DEV 拟南芥 NC_003074.8 | 23.46 Mb | 约 23 min | 29–30 min |
| DEV 水稻 NC_089041.1 | 29.94 Mb | 约 29 min | 37 min |
| DEV 两条合计 | 53.40 Mb | **约 51 min**；三个候选分支（fresh-joint / uniform-coarse / codon）约 2.6 h | 约 67 min；三个分支约 3.4 h 串行 |
| 全 TRAIN（6 条染色体） | 156.78 Mb | **约 2.5 h** | 约 3.2 h |
| 任意新基因组 | 每 100 Mb | 约 1.6 h | 约 2.1 h |

要点：
- 每次 ont 推理只占 3 GiB 显存。他人占 11 GB 时仍有约 13 GB 空闲，所以只要对方不涨显存，我们可以同时存在；但仍遵守"一次只跑一个 GPU 任务"。
- 推理耗时约 75% 在 GPU 编码、25% 在单线程 CPU 的 skyline 解码（区段 A：编码 85 s，解码 28 s）；ont 有 32 核，在不改科学逻辑的前提下，把两条链/两条染色体分到两个进程能掩盖解码时间，这只是潜在提速，未测。
- 从 Baobab 排队看，ont 的真实收益主要是免排队（04 号报告记录 10-14 的 Baobab 队列为 102 运行 / 1,038 等待），而不是省 credits：Baobab 上这类推理作业 billing 只有 12/h，一个 DEV 拟南芥作业约 6 credits，整个三分支 DEV 家族约 40 credits。credits 压力主要来自 CPU 数组和大内存请求，不是这个管线的推理。
- 训练没有测。Baobab 上 fresh-joint 续训 4,608 步用了 47,009 s（13 h，约 10 s/步，3090，限时 18 h）。ont 上的显存峰值与速度未知；而且训练还依赖第一次拟合检查点、优化器状态、fresh-draw 计划和标签等 `outputs/` 中间产物，目前没推到 ont。若要在 ont 上训练，需要先做一次 §7 的小步数显存试跑，再决定。

## 6. 遇到的门禁与路径问题及处理

| 问题 | 处理 |
|---|---|
| `infer_full_grammar_joint_controls_dev.py` 要求 `SLURM_JOB_ID`，并 `require_train_closure()`（30 个 TRAIN 作业 + `sacct`） | 不改源码；`scripts/ont_smoke/infer_fasta.py` 直接调用其最后一步 `annotate_full_grammar_fresh_joint24.main`。输出的 `summary.json` 里没有 job_id 字段，因为那是外层 wrapper 的回执 |
| `infer_full_grammar_fresh_joint24_train.source()` 同样核对 FIT/VERIFY 作业号 | 同上，绕开。检查点本身的 `step/epoch/protocol/configuration` 校验由 `validate_checkpoint` 照常执行 |
| `evaluate_full_grammar_dev_baselines.py` 要 `sacct`、三基线（ANNEVO/Tiberius/Helixer）原生输出目录与资产回执 | 没有在 ont 上运行完整评估。改用其底层、未修改的 `common_ruler`（`select_reference`、`S.key`）对区段做 exact-chain 评分。若要在 ont 上做完整 DEV pooled 评估，需要推三基线预测与回执，并给 `check_terminal_rows`/`sacct` 加绕开 |
| 检查点/summary 内保存的是 Baobab 绝对路径 `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna/…` | 只是元数据，推理不读取。ROOT 由脚本位置推断（`parents[3]`），在 ont 上自动指向 `~/coding-rna`；`outputs` 在 Baobab 是符号路径，rsync 到 ont 变成普通目录，没有问题 |
| 缺失模块（01 号报告的发现 0.1：本地仓库缺 `core/training/export/...`） | 不是 ont 的问题：Baobab 快照有完整模块；以 Baobab 为源推送后 import 全部成功 |
| `c0_decode.DEFAULT_HMM` 指向 `refs/repos/annevo-2026/src/HMM.py` | 只有测试用到；推 285 KB 的 `src/`；生产推理未触发 |
| `ont_sync.sh` 默认指向 ab-initio-TE | 用 `PROJECT=… REMOTE_ROOT=coding-rna` 环境变量覆盖，未改脚本。注意它只推不删，同名文件会被覆盖 |
| `A.main` 用 `mkdir(exist_ok=False)` | 每次运行用新输出目录（`run_ont.sh` 已检查） |
| ssh 嵌套 | `ssh baobab 'ssh lab "…" < /dev/null'`；写文件时用管道把内容送入 `cat >`，不用 `< /dev/null` |

## 7. 可复用的启动命令

```bash
# 0) 同步代码（Baobab 为准；会覆盖同名文件）
ssh baobab 'PROJECT=/home/users/j/jwang/coding-rna REMOTE_ROOT=coding-rna ~/bin/ont_sync.sh push src scripts tests configs'
# 检查点（最终 step）
ssh baobab 'PROJECT=/home/users/j/jwang/coding-rna REMOTE_ROOT=coding-rna ~/bin/ont_sync.sh push outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FIT/step_009216.pt outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FIT/summary.json'

# 1) 启动前先看 GPU 和登录用户
ssh baobab 'ssh lab "nvidia-smi --query-gpu=memory.used,memory.free,utilization.gpu --format=csv; who" < /dev/null'

# 2) 生成输入（整条 DEV 染色体；或给起点/长度取对齐的区段）
ssh baobab 'ssh lab "cd ~/coding-rna && ~/envs/generanno/bin/python scripts/ont_smoke/make_region.py arabidopsis_thaliana NC_003074.8 0 0 outputs/ONT-SMOKE/NC_003074.8_whole_input.fasta" < /dev/null'

# 3) 单个 GPU 推理（空闲显存 < 6000 MiB 就拒绝启动；nohup；打印 PID；日志 logs/<name>-<时间戳>.{out,err}）
ssh baobab 'ssh lab "bash ~/coding-rna/scripts/ont_smoke/run_ont.sh dev_ara_fresh_joint ~/coding-rna/outputs/ONT-SMOKE/NC_003074.8_whole_input.fasta 6000" < /dev/null'
#   输出：outputs/ONT-SMOKE/dev_ara_fresh_joint/prediction/{full_grammar.gff3,summary.json}、ont_run.json（吞吐、峰值显存）

# 4) 监控
ssh baobab 'ssh lab "tail -n 3 ~/coding-rna/logs/dev_ara_fresh_joint-*.out; ps -p <PID>" < /dev/null'

# 5) 区段/整染色体对比 Baobab（CPU）
ssh baobab 'ssh lab "cd ~/coding-rna && SEQ=NC_003074.8 SPECIES=arabidopsis_thaliana LENGTH=23459830 ~/envs/generanno/bin/python scripts/ont_smoke/compare_region.py <ont.gff3> <offset> <region_len> <baobab.gff3> <out.json>" < /dev/null'

# 6) CPU 测试
ssh baobab 'ssh lab "cd ~/coding-rna && CUDA_VISIBLE_DEVICES= PYTHONPATH=\$PWD:\$PWD/scripts/experiments/M28-METHOD-RESTART ~/envs/generanno/bin/python -m pytest -q -p no:cacheprovider tests/test_m28_full_grammar_*.py" < /dev/null'
```

其他分支（uniform-coarse、codon、LoRA）需要各自的 `annotate_*` 适配器和检查点，推理同理；换分支时调用对应的 `annotate_full_grammar_<branch>.main`，并推送对应检查点，不必改 wrapper 以外的任何东西。

## 8. 后续建议

1. 用 ont 跑整条拟南芥 DEV 染色体（约 23 min，3 GiB 显存）一次，对 Baobab 的 5,233 条链做全量逐基因比对，把区段级结论升到染色体级。成本为零 credits，不影响 Baobab。
2. 之后凡是"已有检查点的推理/评估/对照（含新物种、重解码）"优先放在 ont；Baobab 只保留训练与需要多卡/更大显存的任务。
3. 若要在 ont 训练（例如被搁置的 q/v+o_proj LoRA 对照），先做一次 ≤ 20 步的显存/速度试跑，因为当前他人占 11 GB，剩余约 13 GB 可能不够，也要先把训练所需的中间产物选择性推到 ont。该试跑需要用户先确认，因为它在共享 GPU 上长时间占卡。

## 附录：ont 上的薄 wrapper（`~/coding-rna/scripts/ont_smoke/`）

`infer_fasta.py`（核心片段；其余为后台线程用 `nvidia-smi --query-compute-apps` 取本进程显存峰值，并写出 `ont_run.json`）：

```python
sys.path[:0] = [str(ROOT), str(ROOT/'scripts/experiments/M28-METHOD-RESTART')]
import annotate_full_grammar_fresh_joint24 as A
torch.set_num_threads(2)
A.main(fasta, ckpt, out/'prediction')      # 与 infer_full_grammar_joint_controls_dev.py 最后一步相同
```

`run_ont.sh <name> <input.fasta> [min_free_MiB]`：读 `nvidia-smi` 空闲显存，不足则退出 3；输出目录已存在则退出 4；设置与 Baobab sbatch 相同的环境变量（`HF_HUB_OFFLINE=1`、`NUMBA/OMP/MKL_NUM_THREADS=2`、`PYTHONPATH`），`nohup … < /dev/null`，日志写 `logs/<name>-<时间戳>.{out,err}`，打印 PID。
