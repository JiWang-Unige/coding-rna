# A3相对原生基线的逐参考缺口

2026-10-03。固定CPU诊断13359102已完成。ANNEVO正确而A3未选出的2143条主参考中，1573条在A3自由池从未出现，486条生成过但在owner阶段丢失，两类合占96.0803%。这把当前缺口具体定位到自然候选覆盖及跨窗保留，不能将其主要解释为最后质量评分不够；阶段状态不是唯一因果证明。

## 范围与结果可信度

仅分析已完成的A3和六次原生工具输出，不新增训练、模型推理或参考筛选，不接触test/Setaria。逐参考粒度为species、seqid、strand及有序0-based half-open CDS坐标链，phase不在exact-key内。7728条primary参考（拟南芥5437、水稻2291）全部保留，包括unsupported、长链及重叠。

A3阶段来自已经重放验证的outputs/M28-A2-A3-REFERENCE-TURNOVER/references.tsv，仍使用A3 step004608；对齐固定M28-COMMON-RULER-R3参考ledger的全部ID/gene/坐标/unsupported标记，复算四scope的free/owner/selected计数与A3官方评价一致。原生GFF使用同一common_ruler parser和primary选择，重新取得精确链集合，六组TP/预测/参考计数均与原生共同评价一致。

Slurm root/batch/extern均COMPLETED/0:0，实际10秒，batch MaxRSS544944K；申请private-teodoro-gpu、0GPU/2CPU/8GiB/5分钟，环境generanno。三项新检查覆盖CSV False正确解析及阶段嵌套、物种/序列/链向/phase键边界、gross交集及未生成状态计数。另用独立脚本复算7728行/唯一键、CDS几何、pooled/两物种/单双CDS片段全部命中与16种成员组合，与summary一致。原始reference表1757534 bytes，summary23386 bytes，均非空；旧产物不覆盖。

## 完整链命中交集

下表只统计参考真值敏感度：双方都有、仅原生、仅A3、双方都没有四类互斥且覆盖全部参考。它不分类非精确预测，也不提供分层precision或F1；整体准确率、背景和暴露边界仍见[native_baseline_result.md](native_baseline_result.md)。

| 参考范围 | 原生工具 | 双方都有 | 仅原生 | 仅A3 | 双方都没有 |
|---|---|---:|---:|---:|---:|
| 合并7728 | ANNEVO | 3501 | 2143 | 237 | 1847 |
| 合并7728 | Tiberius | 3375 | 2066 | 363 | 1924 |
| 合并7728 | Helixer | 3285 | 1827 | 453 | 2163 |
| 拟南芥5437 | ANNEVO | 2336 | 1605 | 183 | 1313 |
| 拟南芥5437 | Tiberius | 2282 | 1583 | 237 | 1335 |
| 拟南芥5437 | Helixer | 2236 | 1437 | 283 | 1481 |
| 水稻2291 | ANNEVO | 1165 | 538 | 54 | 534 |
| 水稻2291 | Tiberius | 1093 | 483 | 126 | 589 |
| 水稻2291 | Helixer | 1049 | 390 | 170 | 682 |

ANNEVO与A3的TP净差1906=2143仅ANNEVO−237仅A3，不能把净差当作缺失集合大小。三种原生工具至少有一种命中6006条；其中1755条从未在A3原自由池生成（拟南芥1302、水稻453）。三工具都精确命中、A3没有选出的参考为1586条。

## 原生工具命中而A3未选出的阶段状态

| 参考范围 | 原生工具 | A3未生成 | 生成但owner丢失 | Owner最大gain非正 | Owner正gain但未选 | 总计 |
|---|---|---:|---:|---:|---:|---:|
| 合并 | ANNEVO | 1573 | 486 | 49 | 35 | 2143 |
| 合并 | Tiberius | 1509 | 478 | 42 | 37 | 2066 |
| 合并 | Helixer | 1290 | 461 | 42 | 34 | 1827 |
| 拟南芥 | ANNEVO | 1184 | 388 | 13 | 20 | 1605 |
| 拟南芥 | Tiberius | 1166 | 384 | 12 | 21 | 1583 |
| 拟南芥 | Helixer | 1026 | 378 | 12 | 21 | 1437 |
| 水稻 | ANNEVO | 389 | 98 | 36 | 15 | 538 |
| 水稻 | Tiberius | 343 | 94 | 30 | 16 | 483 |
| 水稻 | Helixer | 264 | 83 | 30 | 13 | 390 |

`未生成`指真实greedy自由候选中没有该精确参考链；`owner丢失`指某个窗生成过、但最近中心owner窗没有它；后两类是在原owner池已有该链而未选出的记录。正gain但未选不证明某个唯一竞争链造成损失。A3全部参考状态为未生成3285、owner丢失545、非正gain101、正gain未选59、最终selected3738。

在ANNEVO仅有的2143条中，未生成占73.4018%，owner丢失占22.6785%，末两类共84条（3.9197%）。纯旧池重排序只能直接影响后续选择，无法补出未生成链；A6新增beam池有可能改变覆盖或让链在owner窗出现，但真实DEV结果未完成前不能声称恢复了这些条目。

## 缺口集中在多CDS片段链

此分层仅按参考CDS片段数量，两个互斥组仍覆盖全部7728条。不是用CDS-free UTR exon数定义，也不改变主评价分母。

| 参考组 | 参考数 | A3自由池命中率 % | A3最终Recall % | ANNEVO最终Recall % | 仅ANNEVO条数 | 其中A3未生成 |
|---|---:|---:|---:|---:|---:|---:|
| 单CDS片段 | 1886 | 79.2153 | 70.1485 | 74.7614 | 228 | 117 |
| 多CDS片段 | 5842 | 50.4793 | 41.3386 | 72.4752 | 1915 | 1456 |

多CDS片段占全部参考75.5952%，却占ANNEVO仅有链的89.3607%（1915/2143）；A3在该组自由池覆盖只有50.4793%。这支持将多片段完整链的自然覆盖作为接续诊断重点；尚不能分离链长、组合搜索、window位置、标签目标、上游暴露等因素，不能写成某个机制的因果结论。

## 互补命中不是新基因发现

A3精确命中而三种原生工具均未精确命中的primary参考共有111条（拟南芥88、水稻23；单CDS片段57、多片段54）。它们全部是已经注释的参考基因，不是111个新基因，也不意味着A3整体性能优越。四种方法都未命中1611条；该结果不证明这些参考不可学习或不存在。

## 对当前模型迭代的影响

当前证据支持优先检查自然生成覆盖与owner后可用覆盖，不能继续只用较好看的打分损失解释整体不足。已授权的FULL Q/C0扩大到全部TRAIN窗重训生成器；A6增加greedy+beam自由池并训练完整链质量头。它们按原固定协议继续，不因为本诊断改变阈值、checkpoint、seed、训练集合或sealed边界。

等现有运行完成后，应从其已有完整stage accounting检查：多CDS片段/free/owner命中是否真实扩大，最终新增与丢失正确链是否互相抵消，以及precision/背景代价是否抵消恢复。若覆盖仍不足，新的科学选择需要按具体协议处理；本诊断不自行启动新训练或放宽gate。A3原科学gate FAIL保持，论文成功、独立泛化和当前模型优于ANNEVO都尚未成立。

## 可复现代码与产物

远程规范仓库/home/users/j/jwang/coding-rna；本报告同步本机项目。诊断不提交Git或对外发布。

- scripts/experiments/M28-METHOD-RESTART/diagnose_native_reference.py和tests/test_m28_native_reference.py。
- sbatch/M28-NATIVE-A3-REFERENCE.sbatch；实际提交命令为`sbatch --parsable sbatch/M28-NATIVE-A3-REFERENCE.sbatch`，返回13359102。
- outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC/{summary.json,references.tsv,STATUS,JOBID}。
- outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC/verify_reference_membership.py独立核对脚本；在generanno或标准Python环境运行`python verify_reference_membership.py <该产物目录>`即可复算，不调用原诊断模块。
- outputs/M28-PILOT-R5/logs/M28NATREF_13359102.out/err。
