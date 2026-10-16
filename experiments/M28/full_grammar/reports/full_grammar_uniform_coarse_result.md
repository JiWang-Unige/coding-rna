# 原 TRAIN 粗标签信息追加学习对照结果

2026-10-12 UTC。两条最终9216分支、完整raw推理和配对评价均已完成。保留coarse信息相对无信息通道的完整链F1在拟南芥下降0.908052个百分点、水稻上升1.536665个百分点、pooled下降0.342257个百分点。这里没有一致的物种收益；已知局部非零梯度并不证明coarse信息提高整体准确率。

## 固定对照与证据

两支从原13434603最终4608头及完整AdamW状态开始，同seed1 fresh4608有序draw／4219 unique、joint24576编码、原fine／unknown观察、原fine+coarse分母、importance1、头和optimizer设置。Uniform只把追加阶段epsilon0.1改为2/3，保留C×log3常数。原前4608更新已经使用informative coarse；这不是从零coarse-free消融。双方只读最终9216，不选择中间checkpoint，所有22688 primary参考与unsupported strata保留。

Uniform fit13496011于2026-10-11 16:06:07至2026-10-12 05:06:28UTC完成；verify13496012于05:06:28–05:07:40完成。六raw实际roots依次13509298、13509299、13509300、13509301、13509302、13496013，均05:07:40开始；最长于05:45:05结束。共同评价13496014于05:45:05–05:46:39，配对13496015于05:46:40–05:46:53完成。整链10个actual roots及batch／extern共30条COMPLETED0:0。此前informative基线13494895–13494899已由主消费，不重训。两个分支六条实际caller FASTA与归一化原DNA直接相同，uniform参考读取前完成输入核对。Caller未用参考、标签或缓存，评价和配对未更新模型，未访问DEV／sealed。

## 完整链结果

所有差值定义为informative减uniform，不把macro或worst当作pooled计数。

| 范围 | Uniform TP／预测／参考 | Informative F1 % | Uniform precision % | Uniform recall % | Uniform F1 % | Informative减uniform F1 pp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| NC_003070.9 | 5192/6936/7135 | 72.979116 | 74.855825 | 72.768045 | 73.797171 | -0.818055 |
| NC_003071.7 | 3049/4040/4293 | 72.438672 | 75.470297 | 71.022595 | 73.178927 | -0.740255 |
| NC_003076.8 | 4680/6077/6311 | 74.432873 | 77.011683 | 74.156235 | 75.556991 | -1.124117 |
| NC_089043.1 | 998/1666/1638 | 61.454214 | 59.903962 | 60.927961 | 60.411622 | +1.042591 |
| NC_089044.1 | 959/1686/1601 | 60.650706 | 56.880190 | 59.900062 | 58.351080 | +2.299626 |
| NC_089046.1 | 968/1768/1710 | 56.944042 | 54.751131 | 56.608187 | 55.664175 | +1.279867 |
| arabidopsis_thaliana | 12921/17053/17739 | 73.367643 | 75.769659 | 72.839506 | 74.275696 | -0.908052 |
| oryza_sativa | 2925/5120/4949 | 59.635781 | 57.128906 | 59.102849 | 58.099116 | +1.536665 |
| pooled | 15846/22173/22688 | 70.302624 | 71.465296 | 69.843089 | 70.644881 | -0.342257 |

| 汇总口径 | Informative F1 % | Uniform F1 % | 差值 pp |
| --- | ---: | ---: | ---: |
| macro_species | 66.501712 | 66.187406 | +0.314307 |
| worst_species | 59.635781 | 58.099116 | +1.536665 |

Informative相对uniform少107个primary命中和86个预测：拟南芥少162命中／11预测，水稻多55命中／少75预测。这些是汇总净差，不是逐参考gain/loss交集。拟南芥三条染色体F1均为负差，水稻三条均为正差；单次同起点对照不是独立重复或显著性检验。

## CDS base和背景端点

| 范围 | Informative CDS base F1 % | Uniform CDS base F1 % | 差值 pp | Informative／uniform wholly-background链 |
| --- | ---: | ---: | ---: | ---: |
| arabidopsis_thaliana | 96.842868 | 96.722471 | +0.120398 | 73/63 |
| oryza_sativa | 90.032489 | 89.967463 | +0.065026 | 334/343 |
| pooled | 95.236061 | 95.111499 | +0.124562 | 407/406 |

Pooled背景span覆盖informative 1.009303%／uniform 0.992020%；CDS-base F1和完整链F1的方向不同。Uniform非primary残差6327条中1298条精确匹配其他可评价isoform，5029条既不匹配primary也不匹配其他可评价几何；未匹配并不等同生物学假基因。

## 判断和下一步

本次coarse信息追加学习没有统一改善完整链恢复，不能把局部非零角色梯度提升为性能机制的确认性证据。也不能据uniform pooled更高自动关闭coarse通道：两物种方向相反、原前训含informative信息，且AdamW轨迹和裁剪是干预后果。保留两支成功结果，冻结B/C/D DEV仍全纳入、不给TRAIN赢家筛选或调参；D整链未消费前不启动新DEV。D/E及独立M探索不修改本次对照。

实际源：outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-{FIT,VERIFICATION,TRAIN-RAW,TRAIN-DIAGNOSTIC,PAIRED}；核心结果为TRAIN-DIAGNOSTIC/summary.json、paired_input_match.json及PAIRED/summary.json。权重、原DNA和原预测均保留，不提交Git或对外发布。本结果是原TRAIN探索，不是独立泛化、通道普遍必要性或论文竞争目标达成。

