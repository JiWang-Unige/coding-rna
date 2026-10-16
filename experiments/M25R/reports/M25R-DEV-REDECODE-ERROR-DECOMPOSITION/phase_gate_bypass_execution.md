# M25R R2 phase-gate bypass 执行合同

2026-09-08，用户明确批准将窗口上限提高为158，保持冻结128条等距lineage；1GPU×2小时、8CPU、32GB，零训练，先原路径重放一致，再解释phase旁路。32/64与8/32仅为探索参考。失败、超时、不一致即停止，不重投。不改变候选、边界、ORF、阈值、原结果；不读Setaria、不Git提交、不额外评阅。

输入：`outputs/M25R-E1-PHASE-GATE-PANEL-R2/{frozen_panel,window_coverage}.json`，原epoch1 checkpoint及row601。R2原preflight的64上限历史记录保留；新的158授权记录在本合同，不改写旧预检。

隔离输出：`outputs/M25R-E1-PHASE-GATE-BYPASS-R2/`。新脚本`phase_gate_bypass.py`只用固定runs及原motif picker重放；全128条A检查完成之前不计算B。B仅删除phase拒绝，几何phase、ORF和boundary minima/threshold与原规则相同。保存每窗口原float16 logits，A重放逐条记录、B配对记录、原writer的A/B GFF3和独立结构核验；不计算总体F1/FPR，不把reference-nonexact称生物学假阳性。

推理采用原R4 `load_inference_model/load_checkpoint`，同6144窗口、全方向化序列坐标、A padding、bf16→float32→float16、eval mode。同一窗口只forward一次，A/B共享分数。离线使用已有HF模型cache，避免重新下载或改变revision。新测试与原诊断/面板测试在同一作业中先运行，失败即退出、不加载模型。

资源依据：R4历史sacct MaxRSS13833M、gpumem2818M，同RTX3090 24GB已成功。原validation计106792502方向化碱基，每epoch约17382个以上tile加1536训练窗口；三epoch约56754个以上forward，84:11:06总耗时折算约5.34秒/窗口。按此粗略均摊158窗口约14分钟；这含非推理开销、受节点负载影响，不是实测预测保证。保留2小时硬上限，预计有加载/数据读取余量；新作业不做额外timing pilot或新窗口。预计scores约21MB，已确认存储有余量。

采用private-teodoro-gpu、1张原型号3090；live队列显示gpu034有1张未分配GPU，无维护预约，其他项目作业保持不动。提交前Bash与Python语法检查通过。唯一提交命令：`sbatch --parsable sbatch/M25R-E1-PHASE-BYPASS-R2.sbatch`。结果按真实终态补充，不将提交视为完成。
