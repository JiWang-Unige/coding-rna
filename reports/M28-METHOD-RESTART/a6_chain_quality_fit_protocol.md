# A6 固定完整链评分拟合与 raw-FASTA DEV 协议

2026-10-03，正式运行前固定。八窗零更新preflight 13358869已通过；它不是训练或准确率证据。A6延续同一研究授权，不修改运行中的全TRAIN Q/C0、原生基线、历史FAIL、封存test或论文claim。

唯一生成器为A3 step004608，GEN revision b0483c23b6b63787b61a6d3a204a9b517d6ba345；全部冻结。使用原1536条TRAIN draw（1491唯一窗）及原importance_weight_uniform_within_stratum顺序，三遍，共4608个评分头更新；不扩大到全25512窗。本阶段是原A3之后额外头曝光与新增165635参数，不是等曝光的单因素对照。

先在1491唯一TRAIN窗上仅用冻结生成器/DNA/特征生成completed greedy与beam4/completed4/max128联合池，完整几何去重保留两种来源，全部写完free_pool后才读取TRAIN R4标签。没有按参考或query existence筛池。TRAIN loss才注入缺失eligible正链，保留free/reference_injected来源；-1保护/未知不入BCE，0是当前主参考相对非匹配，不是生物学真负。

直接复用ChainHead H128，link冻结，seed0 FP32，loss为已知候选均值BCE乘原draw权重，无新局部/link目标、class balancing或长度归一。AdamW lr3e-4、wd.01、betas(.9,.999)、eps1e-8、clip1。全未知/空候选保留既有chain_loss可微零监督语义，记录zero_known_windows；固定计数更新仍可能包含AdamW正则化，不声称每一步都有标签监督。保存1536/3072/4608的head-only权重，主评价只用末004608，不按DEV选择checkpoint、seed或阈值。

推理直接读取既有两DEV raw FASTA，用冻结原生GEN编码器重新前向，不读取DEV特征cache或参考。完整8690方向窗，最后短窗保留。greedy+beam联合池每链用additive+ordered-GRU/nonadditive-window-null logit评分，操作点固定0；保留非正分原始候选，仅沿用nearest owner、全局去重、正gain非重叠选择，原始四scope记录和GFF均保存。没有ANNEVO HMM依赖、救回、参考支持筛选或unsupported参考删除。

只有fit/infer终态、固定配置和完整网格齐全后才打开参考，复用7728 primary共同尺（5437/2291）与完整背景定义、GFF坐标/phase roundtrip、free/owner/selected参考流失诊断。事前gate沿用两物种分别exact-chain F1严格高于固定A、背景span不高于历史C0；这是继续探索gate，不是发表成功、独立泛化或降低原full-TRAIN Q/C0 gate。新原生基线结果另表保留，不替换历史控制。两物种、多seed、独立测试/同源暴露及新颖性证据仍未完成。

入口run_chain_quality.py、evaluate_chain_quality.py。拟合private/1 RTX3090/2CPU/16GiB/3小时，raw推理同资源8小时，评价0GPU/2CPU/8GiB/20分钟。A5实际八窗解码21.36秒（含初始化）、A6反向8.83秒仅为申请余量依据，不是全流程耗时保证。三个独立新输出目录M28-A6-CHAIN-QUALITY-{FIT,INFER,EVAL}；失败保留原产物，不自动扫描或放宽科学规则。定时等待仅复用原固定监控对话和同一heartbeat。

## 实施和提交记录

CPU检查13358966 root/batch/extern均COMPLETED/0:0，23秒，六项检查通过0.142秒：固定末权重/禁止DEV选择、完整8690网格/禁止cache、逐链质量分而非复制query existence、保留greedy/beam联合池、禁止推理注入/重复链、评分头实际165635参数且link冻结。新脚本AST和SBATCH bash语法通过。这些不是性能结果；既有八窗真实beam和头反向preflight不重复运行。

10:41 UTC提交fit13358977（afterok13358966）、infer13358978（afterok13358977）、eval13358979（afterok13358978），均实际sbatch返回，不猜测作业ID。协议在提交前已保存在canonical远程；三个新独立输出、日志、资源和停止时窗固定。尚未取得4608更新/末权重、完整raw-FASTA输出或DEV指标，不能报告已完成训练或新科学PASS。原FULL Q/C0五作业完全不改，原生比较同尺结果另链等待。

复用固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3与同一heartbeat coding-rna-r5-light-monitor/cohort M28-FULL-NATIVE-A6-20261003、gpt-5.6-luna/low、15分钟；三完整分支任一完成或新故障一次交接，普通运行/依赖不汇报。重要结果将继续追加本文，当前没有对外发布或论文主张。

配置整份prompt/ID/类型/ACTIVE/15分钟/固定target已精确回读，10:42:35.787 UTC同一监控对话人工startup实测A6 fit13358977 RUNNING gpu035、infer/eval有效Dependency；实际turn_context模型与low均核实。自然定时及本cohort整链自动回传尚未观察。当前仍处于正式拟合流程，启动不等于已完成4608更新或通过科学gate。

## 拟合终态与推理接续

2026-10-03后续核实：fit13358977的root/batch/extern均COMPLETED/0:0、56分40秒，STATUS COMPLETED、JOBID一致；summary确认4608更新及唯一末权重outputs/M28-A6-CHAIN-QUALITY-FIT/step_004608.pt。冻结生成器可训练参数0、评分头165635，1491唯一TRAIN窗的free池16321候选；每遍原1536 draws仍为free16806、reference_injected555，参考相对label1/0/-1=3028/12935/1398、zero_known_windows3。三遍顺序和计数一致，不把重复draw当作独立基因或三次独立重复实验。

已有validate_head及严格head state加载已在A8零更新TRAIN原型中通过，使用的就是上述末权重，不选择中间epoch；这只确认训练/权重身份，不证明DEV性能。raw infer13358978正在运行，eval13358979和A7配对读出13359152仍等待有效依赖；尚无A6完整DEV指标或科学gate结论。A8完成池读出的独立TRAIN负结果见[a8_quality_complete_train.md](a8_quality_complete_train.md)，不取消或改变本流程。等待仍复用同一固定监控对话/heartbeat，当前对象为M28-FULL-A6-A7-20261003。

## 完整终态（2026-10-03 15:04 UTC）

infer13358978／eval13358979及A7配对13359152已COMPLETED／0:0，完整8690 raw方向窗／7728参考及实际GFF核验完成；CPU结果核验13359929也COMPLETED／0:0。A6原gate FAIL（两物种F1条件失败、背景条件通过），不加epoch、seed或改阈值。详细数值、配对身份、瓶颈及结论边界统一记录于[a6_a7_quality_result.md](a6_a7_quality_result.md)，此前“正在运行”段落是历史快照，不是当前状态。
