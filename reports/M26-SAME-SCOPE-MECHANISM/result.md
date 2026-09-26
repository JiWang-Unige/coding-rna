# M26 R1：内部改进不是基线优势；特异性解释对参考背景高度敏感

2026-09-26。Slurm `13225696`，private-teodoro-gpu / gpu034，0 GPU、4 CPU、16 GiB；COMPLETED，ExitCode 0:0，70 秒，MaxRSS 1,011,516 KiB。三项合成测试通过。没有训练、推理、调阈值或访问 Setaria。

## 同范围的实测比较

全部方法使用相同两条开发染色体、6,450 complete-primary 参考和同一评价器。以下是既有保存预测的重评分，不是最新版本全工具重新 benchmark。原始值见 `summary.json` 和 `comparison.tsv`，完整逐预测分类账保留在 Baobab `outputs/M26-SAME-SCOPE-MECHANISM-R1/`。

| 方法 | 精确链 TP / 预测数 | CDS-base F1 | exact CDS interval F1 | exact chain F1 | 冻结背景 FPR | 全部已注释 gene 背景 FPR |
|---|---:|---:|---:|---:|---:|---:|
| M25R A | 1,389 / 2,098 | 0.3290 | 0.1204 | 0.3250 | 1.2468% | 0.2845% |
| M25R B | 3,351 / 6,188 | 0.8011 | 0.7345 | 0.5303 | 5.6056% | 0.9948% |
| ANNEVO | 5,153 / 7,182 | 0.8859 | 0.8777 | 0.7560 | 6.4468% | 1.1532% |
| Helixer | 4,449 / 7,899 | 0.8678 | 0.8163 | 0.6201 | 8.5564% | 1.9252% |
| Tiberius | 4,956 / 7,427 | 0.8850 | 0.8651 | 0.7143 | 6.2143% | 0.7205% |

“冻结背景 FPR”严格复用 R1 的 predicted CDS-span / complete-primary reference exon-span 补集定义。末列只是预先规定的参考政策敏感性：保留同一预测，将背景改为全部已注释 gene span 的补集。它不是经生物学验证的真 FPR，更不能替换旧指标把 R1 改判 GO。B 的 chain/interval/FPR/count-ratio 四指标与旧 JSON 绝对误差均为 0。

直接决定：不再把 R1 作为当前 SOTA 方法候选继续救分；三套基线的完整链准确率均更高。R1 在旧背景定义下 FPR 较低，但这一特异性优势并不稳健：换为全部已注释 gene 背景后，Tiberius 的 chain F1 更高、FPR 也更低。这个结果仅适用于当前开发范围与参考政策。

## 会改变下一步研究的问题

拟南芥 B 的 FPR 从 12.7439% 变为 1.3113%，水稻从 1.2004% 变为 0.8637%。拟南芥 ANNEVO 从 13.7226% 变为 0.4077%，Tiberius 从 13.6550% 变为 0.3583%，Helixer 从 17.5310% 变为 2.0784%。

因此，“拟南芥高 FPR 主要来自重复区域或短 ORF”的解释目前不能成立为结论。第一优先应分解参考背景差异：被完整 primary 政策排除的编码基因、其他转录本/跨度、非编码或其他 gene 注释各贡献多少；不能把全部差异提前归因于 partial 注释，也不能把更宽背景保护的预测当作正确基因。物种差异和工具排序都可能混有参考政策效应，值得做有限的逐碱基归属检查。

## 候选可达性身份核验

6,450 条参考 ID 与 oracle 台账一一对应。B 的 3,351 条精确链全部位于宽松 oracle 可达集合内；集合外为 0，可达但未恢复为 234。计数关系现在有身份级证据，而非只比较总数。

- 宽松 oracle 可达：3,585；未达：2,865。
- 2,834 条缺少边界候选，30 条边界存在但无法由同一 run 携带完整 fragment，另 1 条图路径存在但不满足 canonical ORF。
- B / 宽松 oracle = 93.4728%。这只是对更宽松上界集合的覆盖，**不是严格 R1 图的选择效率**：oracle 与 B 的阈值过滤和规则仍须区分。
- 固定当前 B、仅完美删去 nonexact 而不新增链时，chain F1 乐观上界仍为 0.6838078，低于当前 ANNEVO 和 Tiberius 实测 chain F1。这不是一个可实现过滤器，也不是所有后处理的普遍上界。

## B 的 2,837 条 reference-nonexact 分解

| 互斥类别 | 数量 |
|---|---:|
| 精确匹配其他完整参考 isoform | 158 |
| intron chain 相同、CDS 端点不同 | 59 |
| 其他同链 complete-primary CDS-span 重叠 | 1,537 |
| 仅反链 complete-primary CDS-span 重叠 | 6 |
| 无 complete-primary CDS-span 重叠 | 1,077 |

合计 2,837。最后一类不等于 intergenic 新基因：它尚未针对全部注释、partial 基因或独立转录/蛋白证据分类。重叠不证明结构近正确；本表没有确认 fusion、split 或 pseudogene。

## 有效性与局限

- 全部方法预测无重复完整链。参考及 A/B/Helixer/Tiberius 的 CDS 均包含末端 stop；ANNEVO 仅拟南芥 1/4,898 条没有末端 stop/ATG，未作事后修正或删除。不存在足以解释工具差距的整体 stop-codon 格式错位。
- 完整参考中拟南芥有 4 条 CDS 总长度不被 3 整除，仍保留在冻结分母；“complete”注释标签不等于 canonical-ORF 子集。
- baseline 命令与实际 Tiberius 日志指向相同 FASTA 路径，且覆盖目标染色体；这不证明历史文件内容逐字相同。三套 baseline 是 2026-06-17 保存产物，模型暴露、上下文和训练预算并未匹配。
- 本结果是开发性测量和方向判断，不是新方法成功、零样本泛化、跨物种确认或 Nature Communications 级新颖性证明。
- 文献已经包含深度学习加结构解码、结构化训练和多层级评价；新的论文主张必须超出这些常识，证明可重复的机制和实际决策后果。
