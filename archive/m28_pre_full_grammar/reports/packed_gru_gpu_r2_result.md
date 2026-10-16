# Packed链头R2：梯度超出原容差，关闭优化分支

2026-09-30。Job **13294046 FAILED 1:0**，private/gpu035，
RTX3090×1、2CPU/8GiB，申请8分钟，实际50s（0.013889 allocated GPUh），
MaxRSS2,080,100KiB。[合同](packed_gru_gpu_r2_contract.md)未改变。

同step001536/AdamW初态、同首个TRAIN窗，两次warmup与两次对照共4次一步更新，
均丢弃且不保存新权重。共享特征逐值相同、自由候选/训练池/标签一致；
192条候选的全部评分分项、四项loss、total、共享特征梯度、clip norm、
一步参数增量和AdamW状态均处于原atol2e-6/rtol2e-4内。
该窗口的自由候选gain没有符号或选择翻转。

但 `parameter_gradients.chain.sequence.weight_ih_l0` 的[384,128]张量中
**6/49,152个元素超限**，最大绝对差 **3.901775926351547e-6**。
因此完整一步数值准入仍失败。已经将所有分项比较落盘后停止，
通过窗口0、速度步数0，没有summary或吞吐结论；不能用loss或最终参数接近
掩盖固定的梯度要求。未到水稻窗口。

[原始比较记录](packed_gru_gpu_R2_equivalence.jsonl)、
[四步trace](packed_gru_gpu_R2_steps.jsonl)、
[作业内环境](packed_gru_gpu_R2_environment.json)、
[原始stderr](packed_gru_gpu_R2_error.log)可核对具体失败key和全部比较。

R1、D1、R2累计GPU分配时间49+35+50=134s（0.037222 GPUh），
另有CPU验证17s/2CPU。原R4成本4.763333GPUh另列，未被重置。
R2没有修改R4、DEV/test/Setaria或原模型文件；优化模块未接入fit/infer。

## 决定

按事前规则关闭此等价优化分支：不扩大容差，不再扫描后向精度开关，
不把优化代码接入科学试验。后续完整方法比较使用原逐链实现。
这只拒绝了一个数值等价优化实现，不说明B1方法精度优劣。
D1的单窗前向线索仍有效，但不足以保证完整训练等价。

已有R4两遍TRAIN trace还提供一个更接近科学问题的诊断：
[原始曝光统计](r4_training_proposal_exposure.json)。
同一组原draw/正链曝光下，拟南芥自由提案覆盖已知正链曝光的比例
第一遍518/2242=23.10%，第二遍578/2242=25.78%；
水稻345/786=43.89%增至374/786=47.58%。
每步用的是当时更新中的模型，分母是有重复的训练曝光，
不是唯一gene数、单checkpoint召回、DEV或跨物种泛化。
计数为未加importance weight的描述统计；参考注入的正链不计自由覆盖。

这些统计支持继续把“实际自由候选覆盖能否转化为最终完整链收益”放在主判断位置，
不能仅凭训练loss下降宣布成功。原R4未完成正式表；下一轮仍须完整C0/B1
同口径评价，不能用此TRAIN结果替代，亦不据它扫描K或改评价分母。
