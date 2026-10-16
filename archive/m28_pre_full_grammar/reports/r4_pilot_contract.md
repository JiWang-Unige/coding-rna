# M28 R4：首个C0/B1公平拟合合同

2026-09-30，在任何M28 optimizer更新或新DEV性能产生前固定。
目标是有限地检验B1完整方案相对修正采样C0的方向，不是完成SOTA证明或承诺期刊录用。
用户已明确授权新模型研发、远程执行和公开Git发布；不另设Pro审批或旧框架闸门。

## 数据、特征与标签

- 仅现有拟南芥/水稻train及DEV allowlist；test和Setaria封存不变。
- 原GENERanno1.2B CDS annotator revision
  `b0483c23b6b63787b61a6d3a204a9b517d6ba345`，eval+frozen，无旧LoRA。
  BF16存储六聚体特征，每24,576bp窗分4个6,144bp块，无special tokens；
  trainable头读FP32特征。两臂完全相同缓存，输入仅DNA/固定权重。
- 自然网格width24,576、stride12,288，保留terminal anchored窗，先定向再右pad。
  使用R2已实现的原始1,536个paired draws，各物种768个，顺序原样重复3次；
  有放回，三遍不是三份独立样本。所有draw及原始q/importance weight不重算。
- 使用新`outputs/M28-LABEL-POLICY-R4/<species>`清单。完整链primary先按
  独立CDS-assessable政策选最长（tie ID），之后应用训练支持mask；不能再换
  更短但容易的isoform。无合格完整isoform的gene只保留longest-all确定局部证据，
  不作完整链正例。其他isoform/未知结构继续保护。
- trainer按window ID加入新标签，不能采用旧R3缓存内已过时的标签元数据。
  DEV参考固定`outputs/M28-COMMON-RULER-R3/result/*.references.jsonl`，
  共7,728条；训练不支持/长链仍在分母。全部8,690个方向DEV窗口推理。

## 两臂拟合

- seed0，hidden128；同初始化shared和segmentation参数，无旧头继承。
  C0 366,191、B1 565,495个可训练参数，额外参数和运行成本单列。
- batch1，每臂恰好4,608 optimizer steps（3×1,536），FP32头，
  AdamW lr3e-4、weight_decay0.01、betas(0.9,0.999)、eps1e-8，
  global gradient norm clip1.0；固定学习率，不扫seed/lr/阈值/K。
- C0 local fine/grouped likelihood；B1必须调用`training.py::joint_losses`，
  不是保留的普通`core.chain_loss`。B1 local/endpoint/link/chain系数均1；
  endpoint/link/chain按每输出列已知正负类均衡，未知/mask零贡献。
  每窗归一化后直接乘实际importance weight；不能除以batch权重和抵消它。
- B1每步先用当前DNA/模型输出生成并保存自由候选，再标注/注入完整正链与确定
  正连接，保存来源与loss。自由生成函数不能接受参考；训练注入不是推理可达性。
- 保存更新数及按步/遍loss和用量；可保存中间checkpoint用于恢复工程进度，
  但唯一主评价checkpoint是完成第4,608步的最后状态，不按DEV择epoch。
  超时后的半成品不能与另一臂完整状态作正式排名。

## 冻结解码与全染色体合并

- B1候选Budget保持R3：1,024bp/类型8端点，每起点/carry/终点类型8条ORF合法
  exon边，每donor8连接，每事件/carry2条beam，最终8链/kb，最多128exon，
  GT/GC–AG，最小intron总长5bp。跨exon携带未完成codon。
  加性proposal近似搜索后非加性链头重评分，不是全局非加性最优。
- C0使用固定ANNEVO源码`37bdd9aa62ddf24fa55941fb827061f7ed49ce53`
  原生HMM，internal min-intron1（物理最短5bp）、默认零转换罚分；
  显式phase列交换、原生1e-3概率截底和原生路径不变。导出该完整gene相对于
  同一路径中把其全跨度改成intergenic的目标差值，含全部发射与进出转换。
  不改为整染色体密集HMM，不以CDS长度归一化，不加外部null阈值。
- C0/B1均将同染色体/链向窗口候选转全局半开坐标，按跨度中点归属最近窗口中心，
  平局较小window start优先。只保留owner生成的候选，记录nonowner-only丢失。
  完全相同链去重，取最大同臂得分；再作最大正增益非重叠选择。
  selector按end/start/输入index排序，DP同分保持旧解，最终前驱回溯。
- B1使用chain-vs-null gain，固定0操作点；类均衡BCE的gain>0不等于校准后的
  自然候选正确概率>50%。另在同一已训练模型、同一自由候选上报告additive-only
  读出，作为诊断而非独立训练的因果消融。各臂gain数值不横向比较。

这是同窗口条件下的完整方案比较。B1-C0的差异包括提案、辅助监督与链评分，
不能单独归因于非加性GRU，也不能声称比较了全染色体全局HMM。

## 资源与停止条件

R4最多8 allocated GPUh：共享特征≤2h、C0拟合≤2h、B1拟合≤2h、
两臂合计GPU推理≤2h；均单GPU计时，私有优先、必要时共享24GiB GPU。
这些是硬上限，不是由12窗计时保证的吞吐。GPU占用期间等待I/O/CPU也计入。
实际model/head/提案/optimizer/解码/I/O、GPU峰值与调度wall time分别记录。

CPU预处理及最终HMM/合并另限累计16 allocated CPUh，单作业≤2CPU/2h、
RAM≤32GiB；原R3及已登记R4标签/合并工程作业另报历史用量，不伪装成训练。
新R4特征缓存上限180GiB、R4总新增产物上限200GiB；不复制已有数据或权重。
超限、OOM、非有限loss、输入/标签错配或参考进入自由推理时停止受影响作业，
保留原产物，报告具体原因；不为“救性能”扩预算或换口径。
明确实现错误可有限修复，但重跑消耗计入同一R4上限，不能重新开始计费。

## 最终结果及方向判断

只有两臂都完成4,608更新和完整DEV输出才形成正式方案表。
报告按物种、宏平均与micro-pooled chain P/R/F1、CDS-base指标、原parent历史副表，
背景span/CDS覆盖、全背景链/Mb（非生物学假基因）、实际自由候选覆盖、
owner/去重/冲突各阶段完整链丢失及成本。GFF3输出primary protein-coding CDS链，
不冒称全UTR/所有isoform。

若B1在两物种方向一致地改善完整链、同时背景负担可接受，再登记有限独立迁移；
此处不预先宣布成功。若预测差但两臂尚未充分拟合，结论是本预算未充分检验，
不能立刻否证架构或启动K扫描。模型/协议选择后的独立数据须另行冻结；
GENERanno暴露不明，不能称干净zero-shot。
