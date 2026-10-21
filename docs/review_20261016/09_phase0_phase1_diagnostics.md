# 09 阶段 0/1 诊断结果（ont，2026-10-16）

范围：执行 06 号文件 §5 的阶段 0（只读诊断）与阶段 1（低成本闸门）。Baobab 全程只读（rsync 拉数据到 ont，无 Slurm、无 git 操作）；全部计算在 ont（`~/coding-rna`、`~/coding-rna-reorg`）。冻结代码（`pipeline/m28`、`refs/carbon-a`）未改动：Carbon-A 后处理放在副本 `experiments/M29-CARBON-A-GAP/carbon_a_pipeline_copy/`（只在副本里改了 `model_inference.py` 的 `dtype→torch_dtype`，并加了 worker traceback）；Carbon-A 权重通过 `~/coding-rna/refs/carbon-a-local/`（符号链接，`model.safetensors` 指向 ont HF 缓存里 sha256 `5e0a5818…` 的同一文件）使用。

约定：F1 均为项目 ruler（`common_ruler`，M26 CDS-assessable 主参考，精确 CDS 链）。TRAIN = 原六条染色体（22,688 条主参考）；DEV = NC_003074.8 + NC_089041.1（7,728 条）。我用冻结 ruler 复算了既有数字作校验：Carbon-A 过滤 78.172、未过滤 76.063、E/M 池 73.857、E 72.958、DEV fresh-joint 70.107，全部与既有记录逐位一致。

脚本：ont `~/coding-rna-reorg/experiments/M29-CARBON-A-GAP/scripts/`，本地副本 `/Users/jiwang/Desktop/coding-rna-reorg/experiments/M29-CARBON-A-GAP/scripts/`；小结果（JSON/TSV/PNG，4.5 MB）副本在本地 `…/M29-CARBON-A-GAP/results/`。输出主目录 `ont:~/coding-rna/outputs/M29-CARBON-A-GAP/{R1,R2,R3,R3b,R4,R6,G0,G0b}/`。

## 0. 结论先行（对 06 §5 门槛的判定）

| 项 | 数字 | 判定 |
|---|---|---|
| R1 转移网络消融 | 位置无关常数代价（在 TRAIN 上调优，对消融偏保守）：未过滤 76.05 → 69.3（全零）/70.5（最优常数），掉 **5.6–6.8**；过滤点 78.16 → 72.6/73.4。门槛 ≥4 | **走路线 (a)**；掉点几乎全来自内含子进/出转移（−4.7），基因起止转移仅 −1.3 |
| 奇偶闸门 | head=0 + Carbon-A 固定 offset，经我方 skyline 解码，TRAIN pooled **75.76**（CDS≥96 规则）/75.75（不加），目标 76.06±0.5，差 −0.30 | **通过**；22,072 条 Carbon-A 链中 21,663（98.1%）逐字相同 |
| R3 后验过滤 | E/M 池后验：留一染色体阈值 73.857 → **75.52（+1.66）**；Carbon-A 自身置信 +2.11。把 Carbon-A 置信网络直接套在我方链上，固定阈值 0.1：E 72.96 → **76.30（+3.34）**，E/M 池 → 76.47（+2.61） | (e) 可实现增益 +1.7（自家后验）至 +3.3（借 Carbon-A 置信）；我方后验 AUROC 0.805 < 置信网络套我方链 0.864 |
| R2 工具共识 | 核心 4 工具（Carbon-A、ANNEVO、Tiberius、E）命中并集召回上界 80.1%（F1@0FP 88.95）；"≥2 工具同链"投票 F1 78.99（仅比 Carbon-A 过滤 +0.8）；8 个工具都没命中的参考 4,256 条（18.8%） | 共识不是出路；天花板只有 ≈80% 召回 |
| R4 误差分解 | E 比 Carbon-A 未过滤多 497 条多段未命中，几乎全是内部结构错：complex_chain +296、single_junction_shift +125、exon_extra +121；端点、整外显子漏、缺失与 Carbon-A 持平 | 差距在连接精度，与 R1 一致 |
| G0b 耗时 | 训练单步 7.37 s：编码器 9%、DP 前向 17%、DP 反向 72%（CPU，numba 单线程 FP64），GPU 仅占 ≈11%；FP32 DP 无加速；共享分配器（arena）修复 DP 反向 7.2→4.6 s、推理解码 6.4→2.0 s/Mb（链逐字相同） | **06 的 GPU-h 估算高估约 10×**（见 §8） |
| G0 Carbon-A 跑 DEV | 未过滤 **75.24**，过滤 **77.41**（拟南芥 75.31/77.98；水稻 75.09/76.03）；对比 fresh-joint 70.107、Tiberius 71.805、ANNEVO 75.703 | 补齐 DEV 对照 |
| DEV 整条拟南芥复现 | ont 重跑与 Baobab 整染色体：5,234 vs 5,233 条链，5,222 条逐字相同（并集 99.56%）；同 ruler F1 72.889 vs 72.858 | 区段级结论升到染色体级 |

推荐下一步训练实验：**M29-A1「Carbon-A 固定 offset + 零初始化残差头」**，估算 ≈5 ont GPU-h（见 §9）。

## 1. R3 后验过滤曲线（脚本 `m29_r3_posterior_filter.py`、`m29_r3b_carbon_conf_on_tools.py`）

输入：`outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-TRAIN/<chr>/posterior/summary.json`（E/M 池 22,066 条链的 `log P(g 被选入)`）与既有 E/M 池 GFF；Carbon-A GFF 的 `confidence=` 属性。命中 = 唯一链等于某主参考键。命令：`m29_r3_posterior_filter.py`（无参数），`m29_r3b_carbon_conf_on_tools.py E_lora,EM_pool,ANNEVO,Tiberius,B_freshjoint` 然后 `… summarize`。输出 `R3/{summary.json,curve_*.tsv,r3_filter_curves.png}`、`R3b/summary.json`。复现校验：后验 AP/AUROC 0.9102/0.8054，Carbon-A 过滤 16,896/20,540/78.1716，均与既有值一致。

排序质量（预测池内，主参考标签）：

| 打分 | 池 | AUROC | AP |
|---|---|---:|---:|
| 我方后验 | E/M 池 22,066 | 0.805 | 0.910 |
| 我方 gain | 同 | 0.573 | 0.762 |
| Carbon-A 置信（自带） | Carbon-A 22,072 | 0.833 | 0.926 |
| Carbon-A 置信套在我方链上（R3b） | E/M 池 | 0.852 | 0.928 |
| 同上 | E | 0.864 | 0.928 |
| 同上 | B(fresh-joint) | 0.879 | 0.929 |
| 同上 | ANNEVO / Tiberius | 0.830 / 0.863 | 0.931 / 0.930 |

按物种：后验 AUROC 拟南芥 0.793、水稻 0.831；Carbon-A 置信 0.836/0.837。

F1 随保留比例（pooled，保留阈值相同的整组 tie）：

| 保留 | 后验 F1 | gain F1 | Carbon-A 置信 F1 |
|---|---:|---:|---:|
| 100% | 73.857 | 73.857 | 76.063 |
| 95% | 75.265 | 73.652 | 77.706 |
| 90% | **75.472** | 72.638 | **78.504** |
| 85% | 75.054 | 70.971 | 77.978 |
| 75% | 72.639 | 66.660 | 75.271 |
| 事后最优 | 75.528（阈值 0.573，保留 20,178） | 73.925 | 78.516（阈值 0.184，保留 19,897） |

Carbon-A 的出厂阈值 0.1 就在其最优附近（78.172 vs 78.516）。

**可实现增益（阈值在训练集定、在留出集用）**：后验留一染色体：73.857 → **75.517**（+1.66，六个阈值 0.572–0.576 几乎不变）；跨物种：拟南芥定阈值用到水稻 68.41 → 71.48（+3.07），水稻定阈值用到拟南芥 75.44 → 76.12（+0.68）。Carbon-A 置信同法留一染色体 76.063 → 78.476（+2.41）。同样保留 20,540 条时后验只有 75.451（TP 16,308），Carbon-A 置信 TP 16,896。

**借用 Carbon-A 置信网络做外部过滤器（R3b，未调参，阈值固定 0.1）**：

| 工具 | 未过滤 F1 | 置信≥0.1 | 增益 | 事后最优 |
|---|---:|---:|---:|---:|
| E (LoRA) | 72.958 | **76.300**（保留 19,983） | +3.34 | 76.76 |
| E/M 池 | 73.857 | 76.465 | +2.61 | 76.89 |
| B (fresh-joint) | 70.303 | 74.233 | +3.93 | 74.86 |
| ANNEVO | 75.542 | 77.689 | +2.15 | 77.78 |
| Tiberius | 72.655 | 76.092 | +3.44 | 76.42 |

E 加上这一步与 Carbon-A 过滤版的差距由 5.2 点缩到 1.9 点。**注意**：置信网络训练数据未披露，可能在这两个装配体的 Carbon-A 预测上训练过（06 §1.3 已指出），所以这个增益只能作开发信号；论文里把它写成"借用 Carbon-A 组件"，不能当独立方法。

**校准（可靠性分箱，pooled，十等宽箱）**：后验严重过自信，ECE 0.156（拟南芥 0.144，水稻 0.196）：[0.9,1.0] 箱含 16,988 条，平均预测 0.991，实际命中 0.855；[0.8,0.9] 箱预测 0.854/实际 0.596；[0.5,0.6] 预测 0.552/实际 0.359。Carbon-A 置信校准好得多，ECE 0.067（拟南芥 0.048，水稻 0.149）：[0.9,1.0] 预测 0.944/实际 0.943；[0.1,0.2] 预测 0.149/实际 0.302（保守）。完整分箱在 `R3/summary.json`。结论：06 的 C3（"配分函数后验是校准的基因级概率"）在现状下不成立，后验只是有用的排序，不是概率；若要保留 C3 需要加校准（温度/等渗）并在独立面板验证。

## 2. R4 多段链误差类型分解（脚本 `m29_r4_error_types.py`）

对象：17,293 条多段主参考（2 段 3,620；3–5 段 6,278；6–9 段 4,025；≥10 段 3,370）。分类按固定结构规则、优先级顺序，不用任何分数（`R4/summary.json`、`R4/misses_<tool>.tsv`）：confidence_filtered（仅 Carbon-A 过滤版）→ other_isoform（命中同基因其它可评价转录本，ruler 记未命中）→ absent_no_call → wrong_strand_only → merge（最佳重叠预测同时覆盖另一主参考 ≥60 bp）→ split（≥2 个预测各覆盖本参考 ≥60 bp）→ 同外显子数且内部连接全同只差端点 → exon_missed → exon_extra → single_junction_shift（同外显子数、恰有一个连接被替换）→ complex_chain。分类函数用合成链做过单元检查（端点/漏/多/偏移均正确）。

多段未命中数（总数 / 类别）：

| 类别 | E | Carbon-A 过滤 | Carbon-A 未过滤 | E/M 池 |
|---|---:|---:|---:|---:|
| 合计 | 5,211 | 4,807 | 4,714 | 5,107 |
| other_isoform | 1,387 | 1,322 | 1,330 | 1,407 |
| absent_no_call | 968 | 1,903 | 945 | 1,043 |
| wrong_strand_only | 93 | 91 | 95 | 91 |
| merge | 77 | 73 | 185 | 64 |
| split | 125 | 57 | 113 | 130 |
| endpoint_only | 246 | 241 | 252 | 285 |
| exon_missed | 175 | 167 | 196 | 187 |
| exon_extra | 377 | 115 | 256 | 305 |
| single_junction_shift | 787 | 368 | 662 | 733 |
| complex_chain | 976 | 377 | 680 | 862 |
| confidence_filtered | – | 93 | – | – |

要点：
1. E 对 Carbon-A 未过滤的多出的 497 条未命中里，complex_chain +296、single_junction_shift +125、exon_extra +121、split +12；merge −108（E 更少合并邻基因）、endpoint/exon_missed/absent 基本持平。即**差距是内部连接/外显子结构精度，不是漏基因、不是端点**。
2. endpoint_only 全是起始密码子一侧（E 243/3，Carbon-A 239/2；终止侧 0）：终止密码子由框内第一个终止定死，剩下的端点错是起始 ATG 选择。
3. Carbon-A 过滤版相对未过滤多的 absent_no_call（1,903 vs 945）是置信过滤把"重叠但错"的预测删掉后留下的空位（过滤删 93 条精确命中、959 条错链）。
4. E 未命中而 Carbon-A 过滤命中的 975 条：single_junction_shift 317（33%）、complex 163、exon_extra 158、other_isoform 141、endpoint 82、absent 58、split 21、exon_missed 30。反向（Carbon-A 未命中而 E 命中，571 条）则以 absent 212、single_junction_shift 90、endpoint 85、other_isoform 52 为主。
5. 分层：召回按段数（E / Carbon-A 过滤）：单段 80.4/81.7；2 段 67.5/69.6；3–5 段 70.3/71.9；6–9 段 72.9/75.5；≥10 段 67.9/71.7。物种：水稻多段未命中 E 1,374（complex 341、absent 325、single_junction 203、extra 158）对 Carbon-A 过滤 1,201（absent 659 主导，complex 仅 99）；拟南芥 E 3,837 vs 3,606。

## 3. R2 工具共识上界（脚本 `m29_r2_consensus.py`，`R2/summary.json`）

单工具（同 ruler）：Carbon-A 过滤 78.172（16,896/20,540）、未过滤 76.063、ANNEVO 75.542、Tiberius 72.655、Helixer 63.592、B 70.303、E 72.958、M 72.200、E/M 池 73.857。

并集命中上界（任一工具命中）：

| 工具集 | 命中并集 | 召回上界 | 零假阳 F1 |
|---|---:|---:|---:|
| Carbon-A 过滤+ANNEVO+Tiberius+E | 18,173 | 80.1% | 88.95 |
| Carbon-A 未过滤+ANNEVO+Tiberius+E | 18,228 | 80.3% | 89.10 |
| Carbon-A 未过滤+E | 17,713 | 78.1% | 87.69 |
| 我方 B+E+M+池 | 17,097 | 75.4% | 85.95 |
| 全部 8 工具 | 18,432 | 81.2% | 89.65 |

8 个工具全部未命中的主参考 4,256 条（18.8%），核心 4 工具命中数直方图（0/1/2/3/4 个工具）：4,460/912/1,077/2,014/14,225。

"≥k 个工具同链"投票（唯一链，非重叠后处理版取括号）：核心 4（Carbon-A 过滤、ANNEVO、Tiberius、E）≥2：F1 **78.99**（P 82.0、R 76.2；非重叠 78.95），≥3：78.32，≥4：73.88；核心 3 外部工具（Carbon-A、ANNEVO、Tiberius）≥2：78.56；加 Helixer 五工具 ≥3：78.66；我方 B+E+M ≥2：74.49。混合（Carbon-A 过滤 ∪ 与其不重叠且被另三个工具中 ≥3 个同链支持的链）：78.40。结论：共识最多比 Carbon-A 过滤 +0.8，不构成方法贡献，也说明即使完美合并也只有约 80% 召回——"超过 Carbon-A"必须来自更好的单模型，而非选择器。

## 4. R1 转移网络消融（脚本 `m29_r1_transition_ablation.py`、`m29_score_gff_dir.py`）

方法：Baobab `outputs/M28-CARBON-A-TRAIN-RAW/<chr>/native/*.raw_probs.parquet`（保存的 GPU/flash 推理）→ 副本 `post_processing` 路径（CPU fp32 的转移 TCN、numba 30 状态 Viterbi、置信 TCN、CDS≥96）。变体只替换转移 TCN 输出的四条 log-odds 轨：`g<G>_i<I>` = 位置无关代价（每基因总代价 G 平分到 enter/leave_gene，每内含子总代价 I 平分到 enter/leave_intron，单位 nat；TCN 原输出量级 ±200）；`net` 表示保留 TCN。输出 `R1/<variant>/<chr>/prediction.gff`、`R1/scores_*.json`。

控制（CPU 重放）：未过滤 **76.048**（17,020/22,073），过滤 78.162，对保存的 GPU 结果 76.063/78.172，数值差 ≤0.015。

| 变体 | 未过滤 F1 | 与控制差 | 过滤 F1* |
|---|---:|---:|---:|
| 全零 `g0_i0` | 69.262 | −6.79 | 72.553 |
| 全常数最优（4×4 网格）`g0_i-10` | 70.470 | **−5.58** | 73.395 |
| 基因起止换常数，内含子保持 TCN：`g0_inet`/`g-20_inet`/`g-60_inet`/`g-150_inet` | 73.58/74.31/**74.74**/69.05 | 最优 −1.30 | – |
| 内含子换常数，基因保持 TCN：`gnet_i0`/`i-10`/`i-30`/`i-80` | 69.35/70.76/**71.32**/67.14 | 最优 −4.72 | – |

\* 过滤点仅对在被中断前用真实置信网络跑完的变体有效（`g0_*`、`g-20_*`），其余为跳过置信网络以省时（`DummyConfidence`，过滤==未过滤，表中未列）。网格是 4×4 粗网格且在 TRAIN 上自己调优，因此"掉点"是下界（偏保守）；极端取值（`g-150_i0` 0.14 等）说明无转移时基因/内含子开销是必需的。

判定：门槛"≥4 ⇒ 转移/边界打分是主因"成立（−5.6 至 −6.8，过滤点同向 −4.8 至 −5.6）。更细的信息：**承载信息的是内含子进/出（剪接连接）转移，而非基因起止**（−4.7 对 −1.3），与 R4 里我方多出的错是连接偏移/多余外显子相吻合。保留的限制：我方 15 角色头本身已有供体/受体角色，这个消融证明的是"Carbon-A 流水线里这块值 5–6 点"，不是"我方头缺这块"；后者要靠 §5 的 offset+残差实验检验。

## 5. 奇偶闸门（脚本 `m29_gate_parity.py`、`m29_gate_diff.py`）

**offset 构造**（逐基因逐链，取向坐标；我方 gain[r]=row[r]−row[0]，row[0]=背景=0）：设 `d_i = log p_i − log(1−p_i)`（p 为 Carbon-A 保存的 float16 概率，裁剪到 [1e-6, 1−1e-6]，与 Carbon-A 的 coding/noncoding 得分同口径）。
- 角色 1,2,3（外显子编码相位 0/1/2）、7,8,9（外显子末碱基/供体侧，带相位）、10,11,12（受体后首个外显子碱基）、13（ATG 三碱基）、14（终止密码子三碱基）= `d_i`；角色 4,5,6（内含子）与 0（背景）= 0。这对应 Carbon-A 的 `START_A…STOP` 状态都发射 `log p`，其它状态 `log(1−p)`。
- enter_gene：在 `A` 开头的 ATG 的第一个碱基上，把 `tr[i,0]` 加到角色 13；leave_gene：对终止三联体 `TAA/TAG/TGA` 的最后一个碱基 j，把 `tr[j,1]` 加到角色 14；enter_intron：把 `tr[e+1,2]` 加到外显子末碱基 e 的角色 7/8/9（要求 `seq[e+1:e+3]` 是 GT 或 GC），GC 额外加 GC-AG 先验 `log(0.009/0.99)`（GT-AG 为 0）；leave_intron：把 `tr[i,3]` 加到每个碱基的角色 10/11/12（这些角色只在受体后首碱基出现）。`tr` 为转移 TCN 四条轨（cache：`_cache/tt_<chr>_<strand>.npy`）。
- 头输出 = 0（等价于只用 offset），不训练；解码器是不改动的 `decode_full_grammar_skyline`（cap 128），逐 contig 逐链整条解码。

**逐项语法差异（记录）**：我方无 AT-AC 内含子（Carbon-A 先验 0.001）；内含子最短 5 bp（Carbon-A 4 bp）；起始/终止密码子不能被内含子切开（Carbon-A 允许终止密码子被内含子分开，且 stop 紧随内含子时按 leave_gene 计费）；N 不是通配符（我方含 N 的 CDS 不可编码，Carbon-A 视 N 为任意碱基）；无 record 边缘截断基因/边界内含子（Carbon-A 输出 incomplete 基因，置信 NA）；外显子数 ≤128（Carbon-A 无）；偏移是一次性的，没有 Carbon-A 的定点数量化（1/16384 nat，我方 FP64）；Carbon-A 的 CDS≥96 bp 与 record≥100 kb 规则作为后过滤加到闸门上（对 F1 影响 +0.015）。

**结果（TRAIN 六条，`R6/diff_A_full.json`）**：

| | TP | 预测 | F1 |
|---|---:|---:|---:|
| 闸门（A_full）CDS≥96 | 16,967 | 22,101 | **75.764** |
| 闸门不加 CDS≥96 | 16,967 | 22,110 | 75.749 |
| Carbon-A 未过滤（既有） | 17,023 | 22,072 | 76.063 |

差 −0.30，在 ±0.5 内：**通过**。逐链：Carbon-A 22,072 条里 21,663（98.1%）与闸门完全相同；Carbon-A 独有 409 条（命中 76）、闸门独有 438 条（命中 20）。Carbon-A 独有链的原因：AT-AC 内含子 117+6+1（≈124）、内含子<5 bp 19+6、终止密码子被内含子切开 1+1、record 边缘截断 2、其余 263 条为"我方语法可表示"的同等得分附近的近似平局/定点量化/N 通配符差异，未逐条归因。分染色体命中数闸门比 Carbon-A 少：拟南芥 −16/−8/−20，水稻 −8/−4/0。

对照：同一解码器只用发射项 `B_emission`（无转移）：**69.07**（CDS≥96），掉 −7.0，与 §4 全零 69.26 一致；区段 smoke（NC_089044.1 6–9 Mb）：A_full 与 Carbon-A 同轨 66/67 条一致，B_emission 47/67。

判定：闸门通过，路线 (a) 的起点等同 Carbon-A 未过滤（75.76），可进入训练。

## 6. G0b 耗时剖析（脚本 `m29_g0b_dp_profile.py`、`m29_g0b_gpu_profile.py`、`m29_skyline_*`、`m29_reverse_arena.py`）

设置：真实 fresh-draw 训练窗（取计划前 12 步，上下文 49,152 bp = 3 个 joint24576 窗，核心 24,576 bp），初始化用原 final4608 检查点和一次性 AdamW 状态，只做计时、不保存模型。ont 4090 与另一用户任务共享（对方占 ≈11 GB、利用率 70–100%）；测量时主机负载 ≈4.7–6（其余时间我方 CPU 并行任务在跑）。

**训练单步（12 步均值，范围）**：

| 部分 | 秒 | 占比 | 位置 |
|---|---:|---:|---|
| 编码器前向（3 窗，BF16，冻结） | 0.664（0.62–0.78） | 9.0% | GPU |
| 头前向 + 传回 CPU | 0.161 | 2.2% | GPU→CPU |
| 配分函数 DP 前向（logZ_all + logZ_obs，FP64） | 1.267（1.09–2.77） | 17.2% | **CPU 单线程 numba** |
| DP 反向（反向消息 + 前向重算） | 5.267（4.32–7.52） | 71.5% | **CPU 单线程 numba** |
| 优化器 | 0.006 | 0.1% | GPU |
| 合计 | **7.365**（6.20–9.47） | | |

峰值 GPU 显存：头部训练 2.58 GB allocated/2.77 GB reserved（冻结编码器）；LoRA(E) 编码器前向+反向（3 窗）：前向 0.67 s、反向 0.17 s，峰值 6.2 GB allocated/6.8 GB reserved。Baobab 3090 上同配方实测 10.2–10.4 s/步（E：47,829 s/4,608 步），所以 ont 单核约快 1.4×。**GPU 只忙 ≈11% 的时间**；FP64 配分函数是 CPU 上的 numba 代码（`parallel=False`，无 nogil），换 4090 对它没有任何帮助。DP 峰值 RSS 约 1.4 GB/进程。CPU 单独 DP 画像（`G0b/dp_profile_cpu_*.json`，4 个窗）：前向 logZ_all 1.1–1.4 s、损失前向（all+obs）2.5–2.9 s、反向 10.8–14.2 s（机器负载 35 时）。

**可行的加速（已实测）**：
1. FP32 DP：前向 1.0149 s/Mb vs FP64 1.0105 s/Mb（1 Mb 真实链，访问量相同），**无加速**，且最优分数偏差 2.5（33,307.4 vs 33,305.0），不建议。
2. 反向传播内存分配（`_reverse_messages` 每个 256-bp 检查点块新分配 `(257, 19078)` FP64 数组）改为一次性共享分配器：同一输入同一机器，DP 反向 7.23→4.65 s、7.35→4.56 s、7.35→4.91 s、7.08→4.34 s；单步 DP 总计 8.4→5.7 s（−32%）；**梯度逐位相同**（`max_abs_grad_diff` = 0.0，损失完全相等）。副本 `m29_reverse_arena.py`，用 monkeypatch 即可，不改冻结文件。
3. 并行：窗口间互相独立，可用多进程（或把 numba 函数拷贝成 `nogil=True` 并用线程池）并行 N 个窗口的 DP，做真 batch N（同时修 H4 的 batch=1）；一个 DP worker 约 1.4 GB 内存，ont 32 核可用 ≤16 个。此项尚未实现，下面估算为外推。

**推理路径（整染色体，DEV 拟南芥 23.46 Mb 在 ont 实测）**：总 1,326 s（17.7 kbp/s，峰值显存 3,044 MiB）：编码 1,010 s（76%：每链 506 s，1,909 窗，0.265 s/窗，含分词、头、取回 CPU、拼接）、skyline 解码 315 s（24%：每链 155–160 s，6.7 s/Mb/链，单线程）。纯编码器前向（不含周边）：0.202 s/窗（batch 1），batch 2/4/8 为 0.224/0.235/0.233 s/窗——**batch 不加速**（已算力饱和，flash-attn SDPA），换算成仅编码器 30.7 kbp/s（每碱基要编码 2 链 × 2 倍重叠 = 4 次 24,576 窗的 1/2 利用）。skyline 解码拆分（2–4 Mb 真实链，单线程）：前向 ≈1.0 s/Mb，其余是 trace 阶段每块 `np.empty((8192, 19078), uint32)`（625 MB）的页分配/清零；改为复用同一块后：**6.42→1.96 s/Mb（3.3×），选中的链逐一相同**（2 Mb、4 Mb 两个规模验证，`m29_skyline_reuse.py`）。
 - 可组合的提速（外推，除前两项外未实测）：(i) 解码用共享分配器 + 各链/各 contig 分到不同 CPU 进程并与 GPU 编码重叠 → 吞吐受编码器限制，≈23 kbp/s（含周边 0.265 s/窗）至 31 kbp/s（周边也流水线化）；(ii) 把窗口步长从 12,288 提到 18,432（窗口少 1/3，核心 18 kb + 两侧 3 kb 边缘；需重测 F1，边缘变小有风险）→ ≈35–46 kbp/s；(iii) 本次测量时 GPU 被他人占用 70–100%，独占时可能再快，未测。
 - 对独立面板（12 物种 8.6 Gb）：现状 17.7 kbp/s ≈ 135 h；(i) 后 ≈ 78–103 h；(i)+(ii) ≈ 52–67 h（未验证 F1）。作对照，Carbon-A 原生推理在 ont 上实测 50–59 kbp/s（23.5 Mb 6.6 min，29.9 Mb 9.9 min，SDPA，98 kb 窗、两链共用一次前向），同面板 ≈ 40–48 h。路线 (a) 若要同时跑 Carbon-A 概率和我方残差头，两条编码器路径的权重相同（GENERanno），可以只算我方的 24,576 窗的 joint 编码，并用同一末层隐状态过原生 CDS 线性头得到 offset（Carbon-A24 的概率，消融显示只比 98 kb 版低 0.5–0.6 点），避免双倍编码——尚未验证。

## 7. G0 Carbon-A 在 DEV（`m29_dev_score.py`；原始 `G0/{ara,rice}_raw*/…raw_probs.parquet`，GFF `G0/{ara,rice}_pred/prediction.gff`）

流程：ont GPU 上 Carbon-A 的 `model_inference.py`（`--attn_implementation sdpa`，BF16，batch 4，98,304 bp 窗/重叠 6,144 bp；FlashAttention-2 在 ont 环境不可用），再用副本 `post_processing.py` 在 CPU 上做转移 TCN/Viterbi/置信（numba，fp32）。拟南芥 23.46 Mb 推理 6.6 min，水稻 29.94 Mb 9.9 min。

| 方法 | 拟南芥 | 水稻 | pooled |
|---|---:|---:|---:|
| Carbon-A 过滤（置信≥0.1 或 NA） | 77.979 | 76.030 | **77.406**（5,687/6,966/7,728） |
| Carbon-A 未过滤 | 75.306 | 75.085 | **75.241**（5,730/7,503/7,728） |
| ANNEVO（既有） | – | – | 75.703 |
| Tiberius（既有） | – | – | 71.805 |
| fresh-joint（既有，复算校验） | 72.858 | 63.725 | 70.107 |

SDPA 与保存的 FlashAttention-2 结果的数值差在 TRAIN NC_089044.1 上核对：概率逐点 |Δ| 均值 7.5e-5、p99 2.4e-3、最大 0.043；F1 74.159→73.946（未过滤）、76.113→75.895（过滤），1,446/1,460 条链逐字相同；即 G0 数字有约 −0.2 的 SDPA 偏差。结论：DEV 上 Carbon-A 过滤版比 fresh-joint 高 7.3 点（水稻 12.3 点），与 ANNEVO 的未过滤差仅 −0.46。

**DEV 整条拟南芥 fresh-joint 复现**（`devAra_whole_compare.json`）：第一次整染色体推理（PID 2084454）在第 1280/1909 窗口时崩溃（`TypeError: 'set' object is not callable`，见 §10），重跑成功：1,328 s、17.7 kbp/s、峰值 3,044 MiB。逐基因比对：ont 5,234 条 / Baobab 5,233 条，完全相同 5,222（并集 99.56%），仅 ont 12 条、仅 Baobab 11 条（几乎都是同基因起点/终点差几 bp 或多/少一个外显子的近似平局，分数差 ≤6.8）；同一 ruler：ont 3,889/5,234/5,437 F1 72.889，Baobab 3,887/5,233/5,437 F1 72.858。

## 8. 修正后的算力估算（对 06 §4/§5）

关键事实：**训练步 89% 是 CPU 单线程 FP64 DP，GPU 仅 11%**；06 把 3090 的 12–16 GPU-h/次 fit 当作 GPU 小时是错的，那是 13 小时墙钟（每步 10.2 s × 4,608 步）。

| 任务 | 06 估算 | 修正（ont，实测+外推） |
|---|---|---|
| 一次 4,608 步头部/LoRA fit（batch 1，原代码） | 12–16 GPU-h | 墙钟 9.4 h（7.37 s/步），GPU 忙 ≈1.0 h |
| 同上 + 共享分配器 | – | 墙钟 ≈6.5–7 h（≈5.3 s/步，外推自 −32% DP） |
| 同上 + 8 个 CPU 进程并行 DP、batch 8（576 次更新） | – | 墙钟 ≈1.1 h，被 GPU 串行编码（8×0.83 s）限制；GPU 约 1.1 h（外推，未实现） |
| 配方试点 G1（4–6 GPU-h） | 4–6 | ≈1.1 GPU-h（并行）/ 6.5 h 墙钟（未并行） |
| 六条 TRAIN 推理 157 Mb | 3.2 GPU-h（3090） | 2.5 h（17.7 kbp/s）；解码改进后 ≈1.7–2.0 h |
| DEV 两条 53.4 Mb | 约 1 h | 0.8 h；解码改进后 ≈0.5–0.6 h |
| 独立植物面板 8.6 Gb（我方推理） | 141 GPU-h | 现状 135 h；解码改进+流水线 ≈78–103 h；再加大步长 ≈52–67 h（后者待验证 F1） |
| Carbon-A 在 8.6 Gb 面板 | 约 3–4 GPU-h/100 Mb | 实测 50–59 kbp/s ≈ 40–48 h |
| (c) 扩大唯一窗口（全部 25,512 窗 ×3 遍） | 200+ GPU-h@3090 | 7.37 s × 76,536 = 157 h 墙钟（单进程）；并行+分配器后 ≈ 15–20 h 墙钟，GPU ≈ 17 h |

注：以上外推均以单步实测分解为基础；并行 DP 的实现与 batch>1 对优化行为的影响（lr、梯度累积）都未测，需要先做 20 步 smoke 再放量。

## 9. 路线判定与推荐的下一步

1. **路线选择**：R1（−5.6…−6.8）+ 闸门通过（75.76）+ R4（差距在连接结构）⇒ **路线 (a)**：Carbon-A 固定 offset（发射 + 转移 + GC 先验）+ 我方零初始化残差头，用完整语法配分函数损失训练。R2 否定共识路线；R3 给出廉价叠加项。
2. **后验/置信过滤**：(e) 的可实现增益是 +1.7（自家后验，留一阈值）；借 Carbon-A 置信网络则 +2.6…+3.9，E 上 +3.3。若主张"自家校准后验"，先要补校准（ECE 0.156 → ≤0.07）并在独立面板验证。
3. **推荐实验 M29-A1**：offset = §5 的 A_full（TRAIN 六条的 Carbon-A 概率与转移轨已在 ont：`_cache/tt_*.npy` + `M28-CARBON-A-TRAIN-RAW`，训练窗口的 offset 只需切片，不再多一次编码）；残差头沿用 366k 参数 C0（末层零初始化，保证起点 = 闸门 75.76）；冻结编码器；4,608 步 fresh 计划；batch 8（CPU 进程并行 DP，梯度累积）+ 共享分配器 + cosine 3e-4→1e-5（H4 修复）。评价：TRAIN 六条 raw + 后验/置信过滤读出、DEV 两条。估算 ont：训练 ≈1.1–1.5 GPU-h（并行实现后；未并行则 6.5–9.4 h 墙钟、GPU 忙 ≈1 h），TRAIN 推理 ≈2 GPU-h，DEV ≈0.6 GPU-h，合计 **≈4–5 ont GPU-h**（06 估 15–20）；工程前置：并行 DP worker + offset 切片训练器（新代码，不动冻结文件）。成功判据：TRAIN 未过滤 ≥ Carbon-A 未过滤 76.06 且水稻不低于它，DEV 同向，再与 +置信过滤叠加比较。
4. 未解决/需要警惕：TRAIN/DEV 对 Carbon-A 属于样本内（06 H7），一切增益只能作开发信号；独立植物时间留出面板仍是主张前提。

## 10. 局限与异常（须如实告知）

- **ont 上出现多次不可复现的解释器级错误**：同一台机器上，约 40 次进程启动中至少 7 次出现与数据无关的奇怪异常——整染色体 DEV 推理第一次运行（`'set' object is not callable`，于 12:18）、`partial_semantics.raw_features`（`Counter.__iadd__ … 'int' object has no attribute 'items'`）、`compare_region.py` 的 `IndexError`、`dev_score` 对合法 GFF 报 "expected 9 columns"（同一文件再跑通过）、Carbon-A 推理 worker 的 `'>=' not supported between 'enumerate' and 'int'` 与 `'int' object is not iterable`（同一输入重跑成功）。dmesg 无 MCE/OOM，但历史上有 perl 段错误。所有这类运行都在原样重跑后通过，**最终报告的数字都经重跑，并且与既有数字（Carbon-A 78.172/76.063、E/M 池 73.857、fresh-joint DEV 70.107、后验 AP/AUROC）逐位一致**，但建议在 ont 上跑任何长任务（尤其训练）前做内存/CPU 自检，并对关键结果保留重算校验；在此之前不要把"静默数据损坏"排除在外。
- G0 的 Carbon-A 用 SDPA 而非 FlashAttention-2（约 −0.2 点）；水稻 DEV 与 TRAIN 水稻的 SDPA 差只在一条 TRAIN 水稻染色体上核对过。
- DEV 上 ANNEVO/Tiberius 的数字沿用任务给出的既有值（75.703/71.805），本次未在 ont 重算。
- R1 的置信网络只对部分变体是真实的（见 §4 脚注）；R1 的常数网格较粗，且为在 TRAIN 上自调优，所以掉点是保守下界。
- R3b 的 Carbon-A 置信网络训练数据未披露，可能含这两个装配体（06 §1.3），增益不可当作独立证据。
- 并行 DP、batch 8、解码流水线与大步长推理的加速为外推或仅部分实测（共享分配器与 FP32 对照为实测），§8/§9 的小时数在实现前都是估算。
- 偏离计划：未使用 `ont_sync.sh` 之外的数据通道；ont 上新建 `~/envs/m29_extra`（pyarrow、matplotlib、scikit-learn 的 `--target` 安装，不触碰 `~/envs/generanno`）与 `~/coding-rna/refs/carbon-a-local/`。

## 附：关键命令

```bash
# 参考对象缓存 / R3 / R3b / R2 / R4
ssh baobab 'ssh lab "cd ~/coding-rna-reorg && PYTHONPATH=pipeline/m28:pipeline/m28/scripts:$HOME/envs/m29_extra ~/envs/generanno/bin/python experiments/M29-CARBON-A-GAP/scripts/m29_r3_posterior_filter.py" < /dev/null'
#   同理 m29_r2_consensus.py、m29_r4_error_types.py、m29_r3b_carbon_conf_on_tools.py <tools>
# R1：3 个进程 × 5 线程，各 2 条染色体
… m29_r1_transition_ablation.py all NC_003070.9 NC_003071.7   # 另两进程：NC_003076.8 NC_089046.1 / NC_089043.1 NC_089044.1
… m29_score_gff_dir.py $HOME/coding-rna/outputs/M29-CARBON-A-GAP/R1 control,g0_i0,…
# 闸门
… m29_gate_parity.py A_full,B_emission NC_003070.9,NC_003071.7 full   # 另两进程各 2 条
… m29_gate_diff.py A_full
# G0b
… m29_g0b_dp_profile.py cpu 4; m29_g0b_gpu_profile.py steps 12 | lora 6 | enc 5; m29_skyline_fp32_bench.py; m29_skyline_reuse_bench.py NC_089044.1 4000000; m29_g0b_dp_arena_bench.py 4
# G0（ont GPU，一次一个）
python carbon_a_pipeline_copy/src/model_inference.py --input <fasta> --output_path <dir> --model_name ~/coding-rna/refs/carbon-a-local --gpu_count 1 --batch_size 4 --dispatch_chunk_size 1 --attn_implementation sdpa
python carbon_a_pipeline_copy/src/post_processing.py --input <raw_probs.parquet> --output <dir>/prediction --viterbi_backend numba --device cpu
```
