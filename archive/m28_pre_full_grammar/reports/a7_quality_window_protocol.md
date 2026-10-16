# A7完整链质量跨窗聚合协议

2026-10-03，在A6完整DEV结果可用之前固定。A7只改变A6生成后候选的跨窗保留方式：去掉nearest owner硬过滤，对每条完整全局CDS链取所有生成窗口的原质量gain最大值，再用原正gain同链向span非重叠选择器。A6权重、greedy加beam4池、逐链质量分、零阈值和既有模型训练全部不变。

## 实验依据与边界

[逐参考缺口诊断](native_reference_diagnostic.md)发现ANNEVO正确而A3未选出的2143条中，486条落在owner丢失状态。旧A4使用query存在分取消owner后恢复TP但增加背景，gate失败，不能把该负结果当作已修好。A7使用新A6的完整链质量分和扩展池，检验质量判别能否控制跨窗聚合的代价；不是重调旧A4、不是证明一定有效。

本项不重训、不新增GPU推理、不扫描阈值或pool策略、不混合其他模型候选。它不能补出A6自由池不存在的链，也不能修复链长或窗口限制造成的缺口。跨窗最大gain不是校准概率；可能增加非精确预测并改变原非重叠竞争，必须测量新增和丢失TP、precision、背景及逐物种差异。

唯一来源为完成的outputs/M28-A6-CHAIN-QUALITY-FIT末step004608、完整8690方向窗raw-FASTA INFER及其既有EVAL。只借用推理时已保存的自然候选/原逐链质量分，不重新执行GLM或读取DEV feature cache。原A6 gate保持，不以A7替换、覆盖或救回A6失败。

## 配对读出与验收

入口先核对同一A6 fit/infer/eval、完整网格、固定生成器与末质量头、无参考/cache/新增优化/封存访问，逐窗使用原quality score和greedy/beam来源。A7四scope的free_records/free_unique与A6一致；在所有跨窗聚合和新GFF3导出完成后才打开参考。

完整主尺仍为拟南芥5437加水稻2291共7728条，unsupported/长链/重叠不删。GFF3真实回读核对ID、坐标和phase，完整链phase另诊断；报告per-species、micro/macro、CDS-base union、背景及全部参考的A6/A7命中交集，不把净TP差替代新增与丢失集合。

事前匹配owner gate要求两个物种分别完整链F1严格高于同一A6、背景span bp不高于A6，参考及背景分母一致，用整数交叉乘法判断F1。即使通过，也只支持该DEV配对解码比较，不代表优于ANNEVO/Tiberius/Helixer、干净zero-shot、独立泛化或论文成功。当前原生基线表另存，训练暴露限制不变。

历史6450 parent辅助指标仅为现有汇总接口要求的副表，不参与7728主尺或gate。代码审阅发现A6评价原来未提供该必需字段，原评价仍PENDING未触发；补齐已有副指标是接线修复，不改变训练、主参考、操作点或结论。针对性测试会实际复现未补齐时的KeyError并验证主尺/副表分开。

## 执行与资源

新入口evaluate_quality_window_pool.py复用现有window_pool、A6记录验证、GFF3和共同尺，不修改正在训练/推理使用的query、quality head或assembly模块。新输出outputs/M28-A7-QUALITY-WINDOW-POOL，summary、A7.gff3、7728行references.tsv和STATUS/JOBID独立保存；原产物不覆盖。

CPU检查先经Slurm private-teodoro-gpu、0GPU/2CPU/8GiB/5分钟，只运行四项新合成检查，不打开真实DEV或训练数据。正式A7在这些检查与A6原评价13358979成功之后执行，private-teodoro-gpu、0GPU/2CPU/32GiB/20分钟；沿用既有A4跨窗聚合的内存余量，容纳更大的beam池。环境generanno，不往base安装。

精确sbatch返回及结果另记录，不预填作业ID。真实执行失败保留产物与错误，科学gate失败也完整保存，不调阈值救结果。长等待只复用既有固定监控对话和同一heartbeat，主对话负责终态结果验证与科学解释。没有Git提交、推送、对外发布或新的首训。

## 实施与提交记录

2026-10-03 11:24–11:26 UTC，源码与本协议经远程checked patch应用。A6 evaluate_chain_quality.py只补齐聚合必需的historical_parent_exact_chain辅助字段，未改common_ruler、原7728主参考、模型、阈值、训练、推理或gate；原A6评价13358979仍有效Dependency，未因这个错误失败或重建。

CPU检查13359148 root/batch/extern均COMPLETED/0:0，17秒，batch MaxRSS714972K；四项新测试0.004秒通过。第一项用实际共同score输出复现缺辅助字段时aggregate的KeyError，再验证补齐后主尺与历史副表分开；其余核对精确A6父结果、原逐链quality gain而非query存在分及自由池守恒、7728参考读出使用的四类交集与unsupported保留。三份Python AST及两份SBATCH bash语法通过。这些检查不是生物学性能证据。

实际命令`sbatch --parsable --dependency=afterok:13359148:13358979 sbatch/M28-A7-QUALITY-WINDOW-POOL.sbatch`返回13359152，原A6三作业13358977/78/79保持不变。A7当前PENDING有效Dependency，没有执行新候选选择或准确率评价。控制器会在原A6完整评价成功后运行A7；没有新的GPU分配或优化步骤。

提交前fresh五项调度检查及输出不存在核对完成，private节点CPU/内存可容纳2CPU/32GiB，无相关维护预约，BeeGFS可用174827424972800 bytes；资源选择不依据DEV成绩。早先一次scontrol多节点参数形式被拒绝，改为分别读取两节点成功，没有修改调度或作业。

固定监控仍复用01a0f71f-8845-7432-bd13-bf54c1808fb3和coding-rna-r5-light-monitor、gpt-5.6-luna/low、15分钟。对象更新为原FULL五作业及A6三段加A7一段；只在整个FULL或整个A6加A7链完成，或有新需处理故障时一次交接，普通依赖静默。此处记录对象；保存配置、启动检查和本新轮次的自然定时/终态回传按实际观察区分，不将旧轮次成功冒充本轮证明。

## 配对终态（2026-10-03 15:04 UTC）

13359152 root/batch/extern COMPLETED／0:0／25s，完整7728成员表与实际GFF核验完成；CPU13359929重算配对身份及原整数gate一致。A7 gate FAIL：拟南芥F1增加但背景增加，水稻F1下降且背景增加；不改跨窗策略、零阈值或父权重。详细结果统一记录于[a6_a7_quality_result.md](a6_a7_quality_result.md)，A6独立FAIL不被A7替代。
