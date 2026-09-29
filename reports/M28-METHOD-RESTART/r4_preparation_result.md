# M28 R4：首训前接线完成，有限pilot已进入执行

2026-09-30。新模型目标继续；本记录不代表训练完成、性能成立或论文完成。
[R4冻结合同](r4_pilot_contract.md)已在任何optimizer更新之前写定。

## Pro19的读取边界与采纳决定

Pro用时5m46s，完整读取公开提交
`2be276e67c4d9dbfcfb9f720d74ed2877bb7c8f7`的R3 integration result、
head timing result、training.py、c0_decode.py，补读core.py、candidates.py、
common_ruler_result及固定ANNEVO Viterbi源码。未访问Baobab缓存、未独立重跑。
这是AI方法商议，不是外部同行评议或科学结果独立复现。

采纳三项直接执行意见：
1. C0仍分窗原生HMM；导出完整gene相对于同一路径改成intergenic的原生目标增益，
   含全部发射/进出转换，不加B1学习null、不除长度。C0/B1共用全局owner、
   去重与最大正增益非重叠选择；不临时改为整染色体密集Viterbi。
2. 完整链primary改为longest-assessable后保留R2原draw/q/权重；batch1不再除权重和。
   三遍是同一有限抽样的重复，不是独立样本；不因primary难支持换短isoform。
3. 4,608步、单seed、固定AdamW/候选限额/末checkpoint可进入首训。
   B1入口锁为joint_losses；正gain不是校准正确概率。
   8GPUh是硬上限而非吞吐保证，两臂完整完成才作正式比较。

商议提问中B1跨窗步骤描述的是已定设计；当时只有局部selector/owner helper，
本轮才补成实际全局assembler。没有把未实现路径当成已有染色体结果。

## 两物种训练标签对齐：完成

| Job | 物种 | 状态 | wall | MaxRSS |
|---|---|---|---:|---:|
| 13292230 | 拟南芥 | COMPLETED 0:0 | 45s | 2,065,256 KiB |
| 13292231 | 水稻 | COMPLETED 0:0 | 30s | 815,332 KiB |

各2CPU/8GiB、0GPU；各8项测试通过。新输出
`outputs/M28-LABEL-POLICY-R4/<species>`，旧R2保留。

拟南芥TRAIN有4个primary结构改变、DEV2个；水稻0个。
两物种既定768draw的完整链positive_ids列表均无变化，unique完整链曝光仍为
2,093/720。这不宣称所有局部mask/标签都未变化。
TRAIN可支持完整正例总数为17,594/4,845；无assessable isoform的65/44个训练gene
只保留确定局部/连接证据，不作完整链正例。没有悄悄把非标准但assessable primary
替换成更短支持isoform。DEV参考ID集合精确匹配冻结5,437/2,291清单。

[拟南芥JSON](r4_arabidopsis_labels.json)、[水稻JSON](r4_rice_labels.json)；
[实现](../../scripts/experiments/M28-METHOD-RESTART/align_labels.py)。

## 共同跨窗选择与C0原生增益：完成

[src/m28/assembly.py](../../src/m28/assembly.py)不接受参考，
转全局半开坐标后按同染色体/链向owner归属，保留逐阶段真实链集合及
nonowner-only损失，去重后作最大正增益选择。
原selector改为标量DP＋前驱回溯，不再为每个候选复制完整已选列表。
小规模穷举验证选择目标不变，覆盖负/零/正增益、冲突、重复、负链与尾窗。

C0调用未修改的原生HMM，评分复用其转换矩阵构造及实际发射列。
初次Job13292310在新评分接口失败（17s，2CPU、0GPU）：
原生返回Python list，而新代码用了NumPy向量索引。
仅在评分入口显式转数组，并修正测试中的路径替换类型。
原失败输出保留在`outputs/M28-ASSEMBLY-R4/result`。

复验Job13292328 COMPLETED 0:0，10s，MaxRSS701,576KiB，
`outputs/M28-ASSEMBLY-R4/repair1`：
**14 tests + 6 subtests passed**。
GT/GC×三相位、双gene情况下，导出局部增益与完整实际路径/移除gene路径的
目标差一致；原生chain输出不变。未做参数或阈值变化。
这仍是工程接线，不是新生物学或模型性能证据。

上述R4非GPU准备作业累计204 CPU-seconds（0.056667 allocated CPUh）。

## 正式pilot调度快照

- **13292334 M28CACHE**：private RTX3090×1、2CPU/16GiB、2h上限；
  已实际运行并写入冻结特征。需要1,491个unique训练窗口＋8,690个DEV窗口，
  共10,181窗，复用R3的12个缓存，不重提取它们。
- **13292356 M28FIT-C0**、**13292357 M28FIT-B1**：
  各private RTX3090×1、2CPU/16GiB、2h上限，
  均依赖13292334成功后运行。提交时尚未optimizer更新，不能称训练完成。
- 特征180GiB、R4总新增200GiB、8 allocated GPUh合同不变。
  CPU/HMM最终处理另计资源；GPU占用中的CPU/I/O等待仍计入GPUh。

训练入口已实现：同1536draw顺序×3遍、FP32头、直接每窗importance weighting、
B1自由候选先写入再注入训练参考，按遍checkpoint并只主评价step004608。
相关脚本：
[cache_features.py](../../scripts/experiments/M28-METHOD-RESTART/cache_features.py)、
[fit_pilot.py](../../scripts/experiments/M28-METHOD-RESTART/fit_pilot.py)。

尚未完成：实际两臂拟合、完整DEV推理与GFF3/共同指标汇总。
完整DEV runner尚待接上已验证的组件；不能把提交任务或合成通过当模型结果。
下一步在缓存和拟合执行期间完成该runner，再处理全部结果。
不恢复旧AutoResearch、不开Setaria、不把期刊目标标为完成。
