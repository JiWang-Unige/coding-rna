# 06 Carbon-A 差距分析与超越路线（只读调研，2026-10-16）

本文件是唯一新增产物。调研期间未修改、移动、删除任何已有本地/远端文件，未提交 git，未 sbatch/srun，未读 `secrets.env`，登录节点只做了读取与秒级文本处理（含把 Carbon-A 技术报告 PDF 经 stdin 流式送到 Baobab 的 `pypdf` 抽文本，未在远端落盘）。

来源标记：
- **[PDF]** = Carbon-A 技术报告（`https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report/resolve/Carbon-A.pdf`，远端无副本，已下载到会话临时目录并解析；页码/节号见各处）。
- **[repo]** = Baobab `refs/carbon-a/f6789ee79e6ef397de333577ffec1cb734496620/`（`README.md`、`provenance.json`、`config.json`、`modeling_generanno.py`、`annotation_pipeline/{README.md,src/*.py}`）。
- **[项目]** = 本项目文件，路径相对 `/Users/jiwang/Desktop/coding RNA`。
- **[HF]** = Hugging Face 仓库元数据（2026-10-16 查询）。
- **[推断]** = 我的推断，没有直接实验证据。**[算术]** = 用已有数字做的代数，没有新实验。
- F1 均为项目 ruler 的完整 CDS 链 exact F1（pooled，拟南芥+水稻，TRAIN=6 条染色体 22,688 条 primary），除非另注。

---

## 0. 结论先行

1. **Carbon-A 与 GENERanno-eukaryote-1.2b-cds-annotator-preview 是同一份主权重**（`provenance.json` 记录的 `model.safetensors` sha256 = `5e0a5818…f71c`，与 Baobab 上 GENERanno 的 HF 缓存 blob 文件名一致；config 逐项一致；Carbon-A 仅把 `modeling_generanno.py` 改了 Transformers 5.x 兼容并新增 `annotation_pipeline/`）。**Carbon-A = 同一个全量端到端训练的 1.2B 双向编码器 + 双链二分类 CDS 头 + 两个小 TCN（转移 11,077 参数、置信 14,433 参数）+ 30 状态 Viterbi + 阈值 0.1 过滤**。我们用的是它的骨干（`full.model` 的 `last_hidden_state`），丢弃了原生 CDS 线性头（仅在 native-prior 实验里用其 CDS−NON_CODING 分数做 28 个零初始化系数），从未用过它的转移网络与置信网络（只把整套流程当外部基线跑过）。
2. **差距的结构很清晰**（来源 `outputs/M28-CARBON-A-JOINT24-PAIRED/summary.json`、`outputs/M28-FULL-GRAMMAR-JOINT-LORA-RECOVERY/summary.json`、`carbon_a_original_train_result.md`）：
   - 单段（单外显子）参考上我们与 Carbon-A 持平（joint24 召回 80.0% vs 81.7%）；**差距几乎全部在多段链**（≥6 段：60.6% vs 73.7%）。
   - 对最佳单模型 LoRA(E)（72.958）而言，与 Carbon-A 过滤版（78.172）的 5.2 点差距 **约 3.0–3.1 点来自多出的预测（精确率侧：22,324 vs 20,540 条预测）、约 2.1–2.2 点来自少命中（16,420 vs 16,896）** [算术]。对未过滤 Carbon-A（76.063）：预测数几乎相等（22,324 vs 22,072），差距 3.1 点**全部是召回**（−603 命中）。
   - 水稻是差距所在：E vs Carbon-A 过滤，拟南芥 −3.6，水稻 **−10.4**。
3. **最可能的差距来源 top-3**：(i) 多外显子连接精度（Carbon-A 的学习型转移/边界打分 + 端到端微调的高质量 CDS 概率轨；我们是 366k 小头 + 最后一层冻结特征 + 无转移项）；(ii) 训练规模/多样性与优化配方（Carbon-A：4.35 Tbp、44M 窗口、batch 256、cosine 衰减；我们：6 条染色体、≤4,219 唯一窗口、batch 1、恒定 lr 3e-4，实测对 lr 极敏感）；(iii) 缺少基因级置信过滤（Carbon-A 自身过滤 +2.11 点；我们的后验过滤在 E/M 池上估计 +1.6，见 §3 H1）。上下文长度**不是**主因（Carbon-A 缩到 24,576 bp 仅 −0.5/−0.6）。
4. **一个必须正视的评价学问题**：拟南芥 TAIR10.1（GCF_000001735.4）与水稻（GCF_034140825.1）都在 Carbon-A 的"训练快照"内（[PDF] Table S1/S5，无 New 标记），且 Carbon-A 报告里 14 个时间留出基因组**不含任何植物**（Table S1 带星号的行只有哺乳、其它脊椎、无脊椎）。ANNEVO 作者声明排除拟南芥与水稻，Carbon-A 对 ANNEVO 的未过滤优势在拟南芥为 −0.09、水稻 +2.35（TRAIN，`carbon_a_original_train_result.md`）。**所以"在拟南芥/水稻上超过 Carbon-A"既不公平（它见过这两个注释），也不构成对论文有说服力的 claim。** 真正有价值的 claim 是在 Carbon-A 快照之后的植物基因组上超过它。
5. **推荐前 3 个实验（成本按 1×RTX 4090 粗估，3090 实测 ×0.6–0.9 换算，见 §5）**：
   1. 只读诊断包（Carbon-A 转移网络消融 + 多段链误差类型分解 + 后验过滤曲线 + 工具共识上界）：约 0–1 GPU-h，其余 CPU。
   2. "Carbon-A 能量作固定 offset + 我方残差头/LoRA 用全语法分区函数训练"（路线 a），先过奇偶性闸门（我方解码器喂 Carbon-A 原生概率应复现 ≈76.06）：闸门 ≈1–2 GPU-h；通过后 ≈15–20 GPU-h（含评价）。
   3. 优化配方修复 + 宽 LoRA(qvo / 更多层) + 后验过滤（路线 b+e）：≈15–20 GPU-h（含评价）。
   之后才做 (c) 多物种与独立植物面板。

---

## 1. Carbon-A 完整架构与训练（问题 1）

### 1.1 编码器与读出头

| 项 | 内容 | 来源 |
|---|---|---|
| 骨干 | 26 层、hidden 2048、intermediate 5632、32 个 Q 头/4 个 KV 头、SiLU、RoPE θ=5e5、6-mer 非重叠分词（vocab 4128）、最大 16,384 token = 98,304 bp；由 GENERator-1.2B 的**因果解码器改成双向编码器** | [PDF] Table 3、§5.3；[repo] `config.json` |
| 读出头 | 两个 strand 头（`num_prediction_heads=2`），各为 `Linear(2048, 6×2, bias=False)`：每 token 输出 6 个碱基位置 × {NON_CODING, CDS} 的 logit，softmax 后取 CDS 概率；反链概率在后处理时翻成转录方向 | [repo] `modeling_generanno.py` L1238–1270；[PDF] §5.3.1 式(2)(3) |
| 训练了哪些层 | **全部 26 层 + 两个头，全参数训练**（无 LoRA、无冻结），首个哺乳动物 checkpoint 随机初始化，之后课程式累加类群 | [PDF] §5.4 Table 4 |
| 目标 | 逐碱基、逐链二分类 BCE（两链之和），**正类权重 3.0**；所有同链 isoform 的 CDS **并集**为标签；含无基因窗口 | [PDF] §5.2.2、§5.4.1 |

### 1.2 训练数据与规模

| 项 | 内容 | 来源 |
|---|---|---|
| 数据 | RefSeq GCF 装配体 2,055 个/2,047 个物种/六大类（哺乳 236、其它脊椎 443、无脊椎 428、**植物 182**、真菌 644、原生生物 114 个物种）；仅取长度 ≥98,304 bp 的 record；源 3.033 Tbp，滑窗增强后 **4.352 Tbp、44.27 M 个窗口**（植物 335 Gbp 源→668 Gbp；植物滑窗步长 49,152 bp，另加 0–5 bp 随机偏移以改变 6-mer 相位） | [PDF] §5.2、Table 2 |
| 课程 | mammal(~0.3 epoch) → +other vertebrate(~0.3) → +invertebrate(~0.3) → +plant(~0.5) → +fungi+protozoa(**2 epoch**)；总计约 500,000 步 | [PDF] Table 4、§5.4 |
| 优化 | AdamW，lr 1e-4（2,000 步 warmup，cosine 到 1e-5；最后阶段到 1e-6），**全局 batch 256 窗口（≈4 M token/步）**，32×A100，DDP | [PDF] §5.4 |
| 暴露量级 | 约 500k 步 × 256 窗口 × 98,304 bp ≈ 1.26×10¹³ bp 样本曝光 | [算术] |
| 小结果 | 去掉植物的"animal"中间 checkpoint 在植物上核苷酸 F1 0.873；植物在加入后升至 0.945/0.957 | [PDF] Fig.11 |

### 1.3 后处理流水线（论文称 annotation pipeline）

| 模块 | 细节 | 来源 |
|---|---|---|
| 推理窗口 | 98,304 bp，窗口重叠 1,024 token(6,144 bp)，重叠处概率平均；最后一窗右对齐；BF16 + FlashAttention-2 | [PDF] §5.5.1；[repo] `constants.py` |
| 转移网络 | TCN：宽度 5 输入卷积 + 8 个残差深度可分离空洞块（dilation 1…128）+ LayerNorm，32 通道，感受半径 512 bp，**11,077 参数**；输入**只有该链 CDS 概率轨（logit）**，输出每个相邻碱基边界上 5 类（stay/enter_gene/leave_gene/enter_intron/leave_intron）分布；在训练集 8,192 bp 窗口上均匀采样，batch 16，AdamW lr 3e-4 cosine，wd 1e-4，softmax temperature 16，50,000 步（单 A100 <1 h，≈6.6 Gbp 曝光） | [PDF] §5.5.2；[repo] `transition_network.py`、`tcn.py` |
| Viterbi | **30 状态**（intergenic、start A/AT、6 种终止前缀类型的 exon、STOP、GT/GC-AG 内含子 first/ready 各 6 种前缀、AT-AC 内含子 6 种、前后两个 boundary-intron）；得分 V(λ)=Σ[e·log p+(1−e)·log(1−p)] + Σlog q(转移) + Σρ(剪接类)；ρ 取 GT-AG 0.990/GC-AG 0.009/AT-AC 0.001（相对 0.99）；基因必须 ATG 起、首个框内终止止；内含子保持阅读框；record 边缘截断基因保留并标 incomplete；边界内含子上限 20,000 bp；输出丢弃 CDS<96 bp、record<100 kb | [PDF] §5.5.3 式(5)；[repo] `viterbi.py` 文档、`constants.py` |
| 置信网络 | 同架构 TCN（14,433 参数）：输入 4 轨（角色是否正确的 log-odds + exon/intron/intergenic one-hot，含邻近基因），对每个完整基因的 exon 碱基特征做 2 头注意力池化 + MLP → 该基因"结构完全正确"的概率；标签=与参考**精确匹配**（正负样本等量采样）；batch 64，AdamW 3e-4，≈14,000 步 | [PDF] §5.5.4；[repo] `confidence_network.py` |
| 过滤 | 默认报告 confidence ≥0.1 或 NA（NA=record 边缘截断基因）；Carbon-A 自测该阈值去掉 1.3% 精确命中与 33.4% 其余预测，gene F1 0.736；置信 AUROC 0.876；置信还用于遗传密码表错配告警 | [PDF] §2.2、Fig.2；[repo] 流水线 README |

**两个辅助网络的训练数据（物种/装配体/划分/是否在被评测装配体上训练）报告与 provenance 均未披露**（项目既有结论，`pretrained_exposure.md`）。置信网络若在训练装配体自身的预测上训练，则 Carbon-A 在 TRAIN 上的"过滤增益"可能偏高 [推断]。

### 1.4 GENERanno cds-annotator-preview 与 Carbon-A 的关系

- 权重：**同一文件**（sha256 一致，`provenance.json` 的 `source_repository=GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`、`source_revision=b0483c23…`；Baobab 上 `refs/carbon-a/…/model.safetensors` 是指向 GENERanno HF 缓存 blob `5e0a5818…` 的符号链接）。
- 是否含原生 CDS 分类头：**含**。`GenerannoForTokenClassification`，两个 `Linear(2048,12,bias=False)`；类别 NON_CODING/CDS。
- 我们是否丢弃了它：M28 管线用 `full.model`（骨干）的 `last_hidden_state`（`src/m28/joint_encoder_lora.py` L95），**不使用原生头的概率**；原生头只在 `native-prior` 与 `native CDS slot` 诊断里作冻结读出。`config.json` 有 `feature_layer`（默认 −1），我们没用过中间层。
- GENERanno 的原作者 CDS 调用器只输出逐碱基二元标注，没有多外显子链/相位（`generanno_native_output_scope.md`）。Carbon-A 的转移网络 + Viterbi + 置信网络是**后加**的结构化层。作者重叠（Carbon-A 通讯作者之一 Qiuyi Li 与 GENERator/GENERanno 线）[推断，未核实]。
- 许可：GENERanno preview = **MIT**；Carbon-A-1.2B = **Apache-2.0** [HF]。

---

## 2. 我方 vs Carbon-A 逐项对照

| 维度 | Carbon-A | 我方（fresh-joint / +LoRA） | 来源 |
|---|---|---|---|
| 编码器是否微调 | 全部 26 层端到端，~500k 步 | fresh-joint：完全冻结；E：仅最后 4 层 q/v 的 rank-4 LoRA(102,400 参数)，dense-QV 18.87 M；qvo 对照已写好未跑 | [PDF] §5.4；`full_grammar_joint_lora_result.md`、`full_grammar_joint_qvo_lora_plan.md` |
| 特征来源 | 末层 → 线性 12 维（每 token 6 位×2 类，位置信息直接保留） | 末层 → Linear 2048→128 → token 级 BiGRU(64×2) → `repeat_interleave(6)` → 与 32 通道 DNA 卷积(k=9,5)融合 → 15 路线性；**token 内 6 个位置靠 DNA 卷积分辨** | `git show snapshot-20261016:src/m28/core.py` `SharedContext` |
| 上下文 | 98,304 bp，重叠 6,144 平均 | 24,576 bp 一次联合编码；12,288 stride；碱基取最近窗口中心；整 contig 做 MAP | `full_grammar_joint24_fit.md`；manuscript Methods |
| 训练数据 | 2,055 装配体/2,047 物种；植物 182；4.35 Tbp、44.27 M 窗口 | **2 物种**，6 条染色体 156.8 Mbp；1,491 唯一窗（原）或 4,219 唯一窗（fresh）；9,216 步（头）/4,608 步（LoRA）；曝光 ≈0.23 Gbp | `full_grammar_fresh_joint24_result.md` |
| 头部容量 | 骨干 1.153 B（全训）+ 线性头 + 11k + 14k | 366,191（+102,400 LoRA） | 同上 |
| 标签 | 每链二分类 CDS，**isoform 并集**；含无基因窗口 | 15 角色：精确角色硬约束 + 粗区域观测通道(ε=0.1) + 未知位置边缘化；主转录本=最长 CDS | [PDF] §5.2.2；manuscript Methods |
| 损失 | BCE(正类权重 3)；转移 TCN：5 类 softmax；置信 TCN：精确匹配 BCE | 完整语法下所有兼容基因集合的 `log Zall − log Zobs`（FP64，256 碱基检查点） | 同上 |
| 优化 | batch 256、cosine、warmup、32 A100 | **batch 1、恒定 lr 3e-4、wd 0.01、clip 1，无调度**（脚本中无 scheduler/accum）；训练日志显示 ≈15% 步被 clip，单窗 NLL 重尾（水稻中位 ≈0.05、均值 ≈75） | `scripts/experiments/M28-METHOD-RESTART/fit_*.py`；Baobab `outputs/M28-FULL-GRAMMAR-*-FIT/training.jsonl` |
| 解码语法 | 30 状态 Viterbi；GT-AG/GC-AG/AT-AC 先验；最小内含子 4 bp；record 边缘截断基因 | 完整 CDS 语法（ATG/终止/GT-GC-AG，**无 AT-AC**，内含子≥5 bp，≤128 外显子，同链不重叠）；精确 MAP（skyline） | [PDF] §5.5.3；manuscript Methods |
| 转移/边界打分 | **学习型转移 TCN**（边界 log-odds，512 bp 感受野）+ 剪接类先验 | **无**：能量=逐碱基 15 角色分数−背景分数，单位权重，null 能量 0，无长度/转移/类先验 | `full_grammar_partition.py` |
| 后处理/过滤 | 置信 ≥0.1 或 NA；CDS ≥96 bp；record ≥100 kb | 无过滤（后验仅在 TRAIN 做过排序诊断：AP 0.91、AUROC 0.81，未生成过滤版 GFF） | `full_grammar_energy_pool_posterior_train.md` |
| 暴露 | 拟南芥/水稻装配体在训练快照内；两小网络的训练数据未披露 | 同一编码器；头只在 TRAIN 染色体训练（DEV 是同物种其它染色体） | `pretrained_exposure.md` |

---

## 3. 差距来源假设与最便宜验证（问题 2）

数字基线（TRAIN，pooled，来源 `full_grammar_joint_lora_result.md`、`carbon_a_original_train_result.md`）：fresh-joint 70.303（TP 15,739/预测 22,087）；E=LoRA 72.958（16,420/22,324）；Carbon-A 过滤 78.172（16,896/20,540）；未过滤 76.063（17,023/22,072）。

按 F1=2TP/(预测+参考) 拆解差距 [算术，顺序不同略有出入]：
- E→Carbon-A 过滤：预测数效应 ≈ +3.0–3.1，命中效应 ≈ +2.1–2.2。
- fresh-joint→Carbon-A 过滤：预测数效应 ≈ +2.5–2.7，命中效应 ≈ +5.2–5.4。
- 结论：LoRA 已把"命中"差距缩小一半多，**剩余差距的主体是预测过剩**；对未过滤 Carbon-A 则是 −603 命中的纯召回差。

分层数据（Carbon-A 过滤 vs E/B 的 primary 召回，pooled）：

| CDS 段数 | 参考数 | joint24 | fresh-joint | E(LoRA) | Carbon-A 过滤 |
|---|---:|---:|---:|---:|---:|
| 1 | 5,395 | 80.0% | 79.6% | 80.4% | 81.7% |
| 2–5 | 9,898 | 61.9% | 65.6% | 69.3% | 71.1% |
| ≥6 | 7,395 | 60.6% | 66.9% | 70.6% | 73.7% |

来源：`outputs/M28-CARBON-A-JOINT24-PAIRED/summary.json`、`outputs/M28-FULL-GRAMMAR-JOINT-LORA-RECOVERY/summary.json`。≥6 段/单段召回比：Carbon-A 0.90、E 0.88、fresh-joint 0.84、joint24 0.76。若 ≥6 段平均约 9 段（**未测，需读参考**），按每个剪接连接独立估算，每连接完全正确率约 Carbon-A 98.7%、E 98.4%、joint24 96.6% [算术，假设]。（用户提到的"净多命中 95.279% 来自多段"是对 joint24 而言。）

| # | 假设 | 现有证据 | 最便宜验证（无 GPU 排前） |
|---|---|---|---|
| **H1** | **预测过剩/缺少基因级过滤**：E 比 Carbon-A 过滤多 1,784 条预测，贡献 ≈3 点 | Carbon-A 自身过滤 +2.11（去掉 1,532 条预测，仅丢 127 个命中；`outputs/M28-CARBON-A-CONFIDENCE-DIAGNOSTIC/summary.json`）；我方 E/M 池后验 AP 0.910 vs gain 0.762、AUROC 0.805；固定 90% 保留预算比未过滤少 471 TP，F1 约 75.47 vs 73.857（+1.6，由 `full_grammar_energy_pool_posterior_train.md` 的数字推算 [算术]）；75% 保留则降到 72.6 | **R3（CPU 分钟）**：用已有 TRAIN 后验输出，画 F1-vs-阈值曲线；同 ruler 下对 Carbon-A GFF 里的 `confidence` 属性算 AUROC/AP，与我方后验同台。注意阈值须在 TRAIN 定、在其它面板验证 |
| **H2** | **多外显子连接/边界打分缺口**（转移 TCN、剪接类先验、高质量概率轨；或我方头的 6-mer 重复瓶颈） | 单段持平、多段拉开（上表）；LoRA 主要提升多段（+3.6/+3.8 vs 单段 +0.8）；native 头 CDS-base F1 与我方相近（Carbon-A 95.86 unstranded，fresh-joint 95.07），说明差距在**精确边界**而非碱基覆盖；原生 CDS 头单线性读出已含逐碱基位置信号（`full_grammar_native_cds_slots.md`，错位 160/160 全变差） | **R1（≈0.3 GPU-h 或 CPU 后端）**：用 Baobab 已有 `outputs/M28-CARBON-A-TRAIN-RAW/*/native/*.raw_probs.parquet` 重跑 `post_processing.py --skip-inference`，把转移网络输出置零/换常数（需在**新工作目录的副本**里改代码，不改 refs），看 F1 掉多少：掉 ≥4 点 ⇒ 转移/边界打分是主因，路线 (a) 优先；掉 <1.5 点 ⇒ 差距在发射质量（编码器微调/数据），路线 (b)(c) 优先。**R4（CPU 分钟）**：对 E 与 Carbon-A 的未命中多段参考做误差类型分类（仅端点错/单连接偏移/整外显子漏/多出外显子/合并相邻基因），区分"边界精度"与"漏外显子" |
| **H3** | **训练规模与唯一窗口覆盖**（尤其水稻；差距 −10.4 vs 拟南芥 −3.6） | Carbon-A：植物 182 装配体；我方头 ≈0.23 Gbp 曝光、仅两物种。**同步数下增加唯一窗口已见增益**：结构续训 55.745→fresh 抽样 63.073（唯一窗 1,491→4,219，+7.3），joint24 66.963→fresh-joint 70.303（+3.3）；complete-core 混合抽样无益（63.68）；TRAIN≈DEV（70.303/70.107）⇒ 不是过拟合，是欠拟合/信息受限 | 已有证据即为主要验证。再便宜一步：**G0b（0.3 GPU-h）** 单步剖析 FP64 配分函数 vs 编码器耗时（见 §4 (c) 成本），决定能否把唯一窗口扩到全部 25,512 个 |
| **H4** | **优化配方**（batch 1、恒定 lr、无衰减、重尾损失、clip 常触发） | lr 四格交互 −20.2（full-rank@3e-4 → 58.98，水稻 37.97；rank4@1e-5 → 65.93）；joint24 同序重复 3 遍后仅 53.7% 的窗口 NLL 变低（`training.jsonl`）；约 15% 步梯度范数>1；Carbon-A 用 batch 256+cosine。native-prior 与 joint 的组合坍塌（水稻 53.98→43.04）同属"大梯度配置水稻崩溃"模式 [推断] | **G1 试点（≈4–6 GPU-h）**：LoRA r4 同数据、同步数，仅改 batch 8（梯度累积）+ cosine 3e-4→1e-5 + 水稻/拟南芥分层，看 TRAIN F1 与单窗 NLL 方差 |
| **H5** | **编码器适配不足（末层冻结特征/只改 q,v）** | E +2.66（6/6 染色体）、dense-QV +1.9；LoRA 对水稻 +5.3 > 拟南芥 +1.9；o_proj、更多层未测 | G2：qvo / 最后 8–12 层 LoRA（≈15 GPU-h）；更便宜：读 `feature_layer` 取中间层做浅层探针（只需少量窗口缓存特征，≈1 GPU-h） |
| **H6** | **上下文** | Carbon-A 98→24 kb：−0.61/−0.48；≤24,576 bp 外的参考仅 9 条；同窗口 24,576 下 Carbon-A24 仍领先 joint24 10.6 点 ⇒ **非主因** | 已有证据；可选 `JOINT98304` 原生上下文诊断已有（`full_grammar_native_context_*`），不再优先 |
| **H7** | **评价口径/暴露使差距被高估** | 拟南芥/水稻均在 Carbon-A 训练快照内；Carbon-A 无植物时间留出；其与 ANNEVO（作者声明排除这两物种）的未过滤差仅 −0.09/+2.35；两辅助网络训练集未披露 | **G0（≈1 GPU-h）**：跑 Carbon-A 于 DEV；**面板设计**：选 Carbon-A 快照之后的植物装配体（候选：`GCF_054954075.1` Amborella 2026-09 注释，见 `independent_transfer_candidates.md`；是否出现在其训练快照需再核，Table S3 中未列） |
| H8 | 标签差异（isoform 并集 vs 最长 CDS） | Carbon-A 非 primary 预测中 37% 精确匹配其它 isoform（我方 E 约 1,449/5,904=25%）；ruler 对这类不计命中 | 低优先；R2 中顺带统计 |

R2（CPU 分钟）：**工具共识/配对上界**。用已有 GFF（Carbon-A、ANNEVO、Tiberius、E、E/M 池）算：(i) 任一工具命中的并集上界；(ii) "≥2 个工具同链才输出"的投票 F1。这是廉价基线，也为路线 (g) 定量；它不构成方法贡献，仅作为"超过它是否必要"的参照。

---

## 4. 候选超越路线排序（问题 3）

GPU 小时均为 1×RTX 4090 的粗估：以 RTX 3090 实测（`full_grammar_joint_lora_result.md`：一次 LoRA fit 13.3 GPU-h、6 条 TRAIN 染色体 raw 推理 3.2 GPU-h；frozen B 13.1/3.3 GPU-h；Carbon-A 单条 30 Mb 染色体 raw 约 22 min）乘 0.6–0.9（假设：耗时以 FP64 配分函数为主 [推断，未剖析]，4090 的 FP64/BF16 吞吐约 3090 的 2.3×但访存仅 1.08×）。DEV（53 Mb）评价约为 TRAIN 的 1/3。

| 排序 | 路线 | 预期增益与依据 | 主要风险 | 4090 GPU-h | 与"部分监督+完整语法+校准后验"兼容 |
|---|---|---|---|---|---|
| **1** | **(a) 以 Carbon-A 能量为固定 offset 接入 CRF，我方头/LoRA 学残差** | 起点=Carbon-A 的发射+转移，奇偶闸门预期 ≈76（未过滤）；再靠我方后验过滤（+1.6–2）≈78；残差训练再争取 +1–3。唯一从"起点即与 Carbon-A 持平"出发的路线 | ① 我方语法≠Carbon-A（无 AT-AC，内含子 ≥5，cap 128，N 处理）⇒ 闸门可能不过；② A42 已示"换一套发射给精确语法解码会崩溃"；③ TRAIN/DEV 对 Carbon-A 是 in-sample，**增益需在植物时间留出面板验证**；④ 方法新颖性弱于独立方法（审稿会视为"Carbon-A + 微调"） | 闸门 1–2；残差 fit 12–16 + 评价 3 | **兼容**：残差 fit 直接用部分标签/未知边缘化；分区函数给后验；offset 为固定项不破坏全语法 |
| **2** | **(e) 基因后验置信过滤** | +1.6（E/M 池，TRAIN，事后阈值）至 +2.0（若 E 的 FP 结构与 Carbon-A 相近）；Carbon-A 自身 +2.11 | 阈值须固定后在其它面板验证；我方后验 AUROC 0.805 < Carbon-A 报告 0.876（不同面板）；后验"不是校准概率"（≥0.9 档实际命中仅 85–86%）；与 Carbon-A 置信网络直接重叠，需校准证据 | 0 GPU（CPU 后验，每染色体 ≤1.75 h×2CPU） | **核心**（C3） |
| **3** | **(b) 端到端结构化微调：更宽 LoRA/qvo/更多层 + 配方修复** | E +2.66；qvo 对照已写好；若 H4 成立，配方修复可能同时降低对 lr 的敏感度。粗预期 +1–3，且 H3/H4 可能共同限制上限 | LR 敏感（交互 −20.2）；背景链上升（wholly-background 407→454）；不修配方就盲目加容量风险高；DEV 从未跑 LoRA | G1 试点 4–6；qvo/12 层 fit 15–20 + 评价 3 | **兼容** |
| **4** | **(c) 扩大训练物种与数据（含不完整注释的部分标签）** | H3 直接证据（唯一窗口×2.8 ⇒ +3.3~7.3）；水稻缺口 10.4 是最大头；与 C2 论点互相支撑 | 成本由 FP64 配分函数主导（≈10 s/窗，训练日志：约 15.6–17.4 ks 每 1,536 步）；全 25,512 窗×3 遍≈ 200+ GPU-h@3090 [算术]；需先剖析并可能改 FP32+缩放；新增物种须核对各工具暴露 | 剖析 0.3；每增加 ~5k 窗口×1 遍≈ 15–20；若 DP 改 FP32 可能降一个数量级 [推断] | **最兼容、最能体现 C2**；但与"超过 Carbon-A"的直接收益来得慢 |
| **5** | **(d) 更长上下文/多窗口一致性** | Carbon-A 缩上下文 −0.5/−0.6；我方 6144→24576 已 +11；再到 98 kb 预期 ≤ +1 | 编码器显存/时间 ×4–6；需重构 DP 窗口 | 评估 3–5；训练 +30–50% | 兼容，低性价比 |
| **6** | **(g) 集成** | E/M 等权能量池 73.857 vs 单模型 72.958（+0.9），事后 TRAIN；强一致交集 74.619；与 Carbon-A 的集成=路线 (a) | 事后设计；未做整染色体新编码评估；M 本身有 lr 不稳 | 每多一个成员≈ fit 13–16 | 兼容（后验在池上定义） |
| **7** | **(f) 换基模** | 无已知基模在完整 CDS 链上显示超过 Carbon-A 编码器 | 见下表 | — | 部分兼容 |

### 4.1 路线 (a) 细化与"先前组合失败"的原因分析

先前结果：joint24 66.963；native CDS prior(6144 块) 64.339；二者合并 62.451，水稻 53.984→43.042（`full_grammar_joint24_native_prior_plan.md`、manuscript L49–70）。原先接法：**28 个零初始化可学习系数**（`Linear(2,14,bias=False)`，读两个原生头的 CDS−NON_CODING logit），与头一起以 AdamW 3e-4、batch 1 训练，新系数 age 0、头 age 4608（或 9216），共享 clip=1。候选原因（**均未被实验区分**）：
1. **尺度与优化器失配**：原生 logit 量级大，单位权重能量对 24k 碱基求和；28 个无归一化的系数在 Adam 下以相同相对步长移动，易过冲；与 lr 四格中"大梯度配置 ⇒ 水稻崩溃"同型。验证：读 `step_009216.pt`（4.4 MB）中的 28 个系数与梯度范数轨迹（只读），再用固定 offset（不学习）重放。
2. **缺偏置项**：Carbon-A 的 BCE 正类权重 3.0 使 logit 整体上移约 ln3≈1.1（官方称"保留了负类校准"，但正类被抬高），其转移网络/先验用来补偿；我方 28 系数无 bias，且无转移项 [推断]。
3. **信息冗余**：原生 CDS logit 是同一末层隐状态的线性函数，joint 编码后头自己已能学到（+11.2），先验的 +8.6 增益在 joint 下只剩冗余——但这只解释"无增益"，不解释 **低于两个单项**，故仍需 1 或 2。

因此建议 (a) 的实现原则：
- **固定 offset，不学系数**：能量 = Carbon-A 风格项 + 头输出，头最后一层零初始化，训练目标不变（`log Zall − log Zobs`）。offset 取自**同一次前向**的原生头概率（不需额外编码）：碱基项 `log p` / `log(1−p)`、转移 TCN 四条 log-odds 轨（可在 CPU 上算）、剪接类先验。事件项可直接加到 `exon_gain` 所用的对应角色 logit 上（donor/acceptor/start/stop 角色行），**配分函数 DP 无需改动** [推断，需读 `full_grammar_partition.py` 的角色索引确认]。
- **奇偶闸门（必须先过）**：零头（head=0）+ offset 通过我方 skyline 解码，在六条 TRAIN 染色体上应复现 Carbon-A 未过滤 ≈76.06（容差 ±0.5）。不过则先逐项对齐语法差异（AT-AC、内含子长度、N、record 边缘），不进入训练。
- 评价面板：时间留出植物。TRAIN/DEV 上的增益只能作开发信号。

### 4.2 路线 (f) 基模候选

| 基模 | 规模/上下文/分词 | 许可 | 24 GB 可行性 | 评估 |
|---|---|---|---|---|
| **Carbon-A/GENERanno 1.2B**（当前） | 1.153 B，98 kb，6-mer | MIT / Apache-2.0 [HF] | 推理峰值显存 BF16 ≈数 GB；LoRA fit 峰值 allocated 6.8 GB（`full_grammar_joint_lora_result.md`），24 GB 绰绰有余，可上更宽 LoRA；全参 AdamW 约 18+ GB+激活，需 BF16/8bit 优化器 | 已是注释监督的最强可得骨干；"换基模"没有已知更优者 |
| PlantCAD2 S/M/L（Caduceus-PS，l24-d768/l48-d1024/l48-d1536） | 参数量与上下文：规模约 0.1/0.3/0.7 B、上下文 8,192 bp、单核苷酸（**文献记忆，本次未核实**） | Apache-2.0 [HF] | 24 GB 够；需 `mamba_ssm`/causal-conv1d 的 CUDA 构建，Baobab 环境风险 | **无注释监督**（暴露干净），单碱基分辨率，对剪接位点/保守性有专长；GeneCAD 已用"PlantCAD2+CRF"，新颖性有限。**最合适作第二编码器与 Carbon-A 特征拼接**（粗上下文 + 单碱基精度），针对多段连接缺口 [推断]；非首选主线 |
| PlantCaduceus_l32 | 225 M 级，512 bp 上下文 [记忆] | Apache-2.0 [HF] | 轻 | 上下文太短，不推荐 |
| AgroNT-1B | NT 型 6-mer，约 1k token（≈6 kb）[记忆] | **CC-BY-NC-SA-4.0** [HF] | 可 | 非商业+相同方式共享，限制成果与代码分发；上下文短；不推荐 |
| Evo 2 7B（1B/40B） | StripedHyena 2，单核苷酸，8k–1M 上下文 | Apache-2.0 [HF 7B] | 7B BF16 权重≈14 GB，仅推理；LoRA 微调不现实；官方 FP8 路径依赖 TE/新架构 | 未见完整 CDS 链证据；成本高、成功率低 |
| NTv3 / SegmentNT | 680 M / 562 M | 公开 | — | Carbon-A 报告：外显子 F1 显著低且不产出完整基因模型，不推荐 |
| GENERator 1.2B（因果） | 同架构的解码器 | — | — | Carbon-A 的原型，无理由退回 |

结论：**不换基模作为主线**；PlantCAD2 作为路线 (a)/(b) 失败后的备选第二编码器。

---

## 5. 推荐实验序列、预算与门槛

| 阶段 | 内容 | 4090 GPU-h | 决策门槛 |
|---|---|---:|---|
| **0（只读，CPU，≤1 天）** | R1 转移网络消融（≈0.3 GPU-h 或 CPU 后端）；R2 工具共识上界；R3 后验过滤曲线 + Carbon-A 置信同口径 AUROC；R4 多段链误差类型 | 0–1 | R1 掉点 ≥4 ⇒ 走 (a) 且把转移项作为核心；掉点 <1.5 ⇒ 重心转 (b)(c)；R3 给出 (e) 的可实现增益 |
| **1（≤2 GPU-h）** | G0 Carbon-A 跑 DEV（≈1）；G0b 单步剖析 DP vs 编码器（0.3）；奇偶闸门（Carbon-A 概率喂我方解码器，主要是 CPU/少量 GPU）（≈1） | ≈2–3 | 闸门复现 76.06±0.5 才进入 (a)；剖析结果决定 (c) 是否需要先改 FP32 |
| **2** | (a) 残差 fit（头 + 可选 LoRA r4，修好配方）+ 六条 TRAIN raw + DEV | 12–16 + 3.5 | TRAIN ≥ Carbon-A 过滤且 DEV 同向 ⇒ 进入面板；否则回到 (b) |
| **3** | (b) 配方修复 + qvo/12 层 LoRA + 后验过滤 | 4–6（试点）+ 15–20 | 与 E 配对比较，6/6 染色体同向 |
| **4** | 冻结**植物时间留出面板**（≥5–8 个 Carbon-A 快照之后的植物装配体，含小基因组以控成本）；运行 Carbon-A（每 30 Mb≈15–20 min GPU）与我方；3 seeds | 视面板：每 100 Mb≈ 3–4（Carbon-A）+ 3–4（我方）；Amborella 700 Mb≈ 25–30 | 判据与 07 号文件 C1 一致 |
| **5** | (c) 多物种/部分标签曲线（k=5…100%）与 marginalize vs background 对照 | 单次 15–20，3 seed×多档约 150–250 | C2 的核心实验 |

说明：Baobab 单项作业 >50 credits 须先问用户（`07_revised_claims_and_plan.md` §0 第 5 条）；重计算优先 ont。以上预算均为粗估，4090 实际速度需用一次 G0b 校准。

---

## 6. 在 Carbon-A 未正式发表条件下的合规使用方式（问题 4）

1. **许可事实**：GENERanno preview 为 MIT，Carbon-A-1.2B 为 Apache-2.0 [HF]；二者**权重字节相同**。Apache-2.0 允许使用/修改/再分发，条件是保留版权与许可声明、对修改文件标明改动；两小网络权重随 Apache 仓库发布。训练数据 bucket（`HuggingFaceBio/Carbon-A-training-data`）与 CADB 的许可**未在博客/模型卡中找到**（博客未声明），其 catalog 数据集返回 401——不使用这些数据，直至核实。
2. **作基线**：安全。固定 `revision f6789ee79e6ef397de333577ffec1cb734496620`，原样运行官方流水线，**两个操作点（置信≥0.1/NA 与未过滤）都报告**，在文中明确"技术报告，未同行评议，2026-10"，并披露：(i) 与我方共享主权重；(ii) 拟南芥/水稻装配体在其训练快照内；(iii) 两辅助网络的训练数据未披露。引用：Channing et al. 技术报告 + HF 仓库 + GENERanno。
3. **作组件**：许可允许，但须：(i) 在方法与图中明示"encoder 来自 GENERanno/Carbon-A，Carbon-A 的转移/置信网络被用作固定 offset"，不得写成独立 encoder 或"零样本"；(ii) 若端口 Carbon-A 的 Viterbi/TCN 代码，保留原版权头与 NOTICE，并标明改动；(iii) 不用 Carbon-A/CADB 的预测去**蒸馏出伪标签再声称优于 Carbon-A**（循环）。
4. **时间风险**：Carbon-A 可能随后发表（Nature 系或其它）；与其路线高度重合（GLM + 转移 + 语法 Viterbi + 置信）。我方差异化应放在：完整语法下的部分监督/边缘化、配分函数后验的校准、植物时间留出评价。是否联系原作者属用户决定；本调研不代为联系。
5. **评价公平性要求**（论文必备）：同一 ruler、同一 primary 政策；对每个对比工具逐一列训练暴露（ANNEVO 排除拟南芥/水稻，Helixer 含两者，Tiberius 含水稻，Carbon-A 含两者）；不引用作者自报数字入榜（面板与 exact 定义不同）。

---

## 7. 未核实 / 局限

- **没有运行任何新计算**。所有 F1/TP/预测数来自既有 summary 与 md；"+1.6"、"≈3.0/2.1 点拆解"、每连接正确率、4090 换算均为我的算术/推断，已标注。
- Carbon-A 技术报告 PDF 经 pypdf 抽文本，表格行列偶有错位；引用的参数均与 [repo] 代码互相印证（转移 11,077/置信 14,433 参数、30 状态、阈值 0.1、窗口 16,384/1,024 token）。
- 两辅助网络的训练数据/划分未披露；Carbon-A 快照内具体窗口、同源隔离未知。
- PlantCAD2/PlantCaduceus/AgroNT 的参数量、上下文标注来自文献记忆（`05_literature_positioning.md` 也未给全），仅许可来自 HF 元数据；采用前需核实。
- `GCF_054954075.1`（Amborella，RefSeq 2026-09 注释）是否在 Carbon-A 训练快照之外只是高概率判断，需用其 catalog 或作者确认。
- 我方后验 AUROC 0.805（E/M 池，primary 标签）与 Carbon-A 报告的 0.876（不同面板、不同标签）不可直接比较，需 R3 同口径复算。
- 训练日志的"重尾/不收敛"观察来自 `training.jsonl` 的摘要统计，不是完整的收敛分析；"FP64 配分函数占主要耗时"是由日志时间量级反推，未剖析。

## 8. 关键资产路径

- 本地：`reports/M28-METHOD-RESTART/{carbon_a_original_train_result,carbon_a_context24_train_result,pretrained_exposure,generanno_native_output_scope,full_grammar_joint_lora_result,full_grammar_joint_qvo_lora_plan,full_grammar_joint24_native_prior_plan,full_grammar_native_cds_slots,full_grammar_energy_pool_posterior_train,independent_transfer_candidates}.md`；`outputs/M28-CARBON-A-{TRAIN-EVAL,JOINT24-PAIRED,CONFIDENCE-DIAGNOSTIC,CONTEXT24-TRAIN-EVAL}/summary.json`、`outputs/M28-FULL-GRAMMAR-JOINT-LORA-RECOVERY/summary.json`；`git show snapshot-20261016:src/m28/core.py`（头）、`src/m28/full_grammar_partition.py`、`src/m28/full_grammar_skyline.py`（能量与解码）。
- 远端 `/home/users/j/jwang/coding-rna/`：`refs/carbon-a/f6789ee79e6ef397de333577ffec1cb734496620/`（`annotation_pipeline/src/{viterbi,transition_network,confidence_network,tcn,constants}.py`、`weights/*.safetensors`）；`outputs/M28-CARBON-A-TRAIN-RAW/<chr>/native/*.raw_probs.parquet`、`outputs/M28-CARBON-A-CONTEXT24-TRAIN-RAW/`（R1/奇偶闸门的输入）；`outputs/M28-FULL-GRAMMAR-{FRESH-JOINT24,JOINT-LORA,JOINT24-NATIVE-PRIOR}-FIT/{summary.json,training.jsonl,step_*.pt}`；`outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-TRAIN/`。
- 外部：Carbon-A 技术报告（Table 2/3/4、§5.5、Table S1–S8）；`HuggingFaceBio/Carbon-A-1.2B`（Apache-2.0）；`GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`（MIT）。
