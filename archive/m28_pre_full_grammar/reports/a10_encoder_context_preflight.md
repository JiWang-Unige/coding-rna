# M28 A10 整窗编码工程预检

2026-10-03。当前 24,576 bp 输入被切成四个 6,144 bp 块，分别执行 1,024-token GENERanno 编码，再拼接成 4,096×2,048 特征。固定模型配置支持更长位置范围，本轮检验同一窗口一次 4,096-token 前向是否可运行。研究假设是跨块表示限制可能值得进一步检验；目前没有证据证明它造成现有漏检，也没有准确率改善结果。

## 固定输入与执行

严格复用已完成 A5 TRAIN 诊断的八个窗口，每物种四个，来源为 outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC/predictions.jsonl。窗口原本由 TRAIN 元数据选择，不是随机准确率样本或独立测试。只保留窗口几何和原 oriented DNA，按既有 full-TRAIN 计划定位缓存；不读取参考 catalog，不使用 free query 或标签字段编码。

固定 released encoder 为 GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview，revision b0483c23b6b63787b61a6d3a204a9b517d6ba345，local_files_only、SDPA、BF16、eval/frozen；不下载新权重。两种模式沿用既有 TRAIN-cache 的非 ACGT→N 规范化，要求四块 token 拼接与整窗 token 完全相同。先按现有 encode_window 逐块复现，要求 BF16 特征与原缓存逐元素相等；不相等则停止，不把输入不一致当作新方法结果。随后对整窗进行一次前向，要求输出完整有限的 4,096×2,048 BF16 特征。

逐窗记录原缓存身份、token/特征形状、两种模式同步耗时、峰值 allocated/reserved 显存，以及仅用于说明输入分布变化的特征绝对差异。保存八份独立整窗特征，预计约 128 MiB；不覆盖旧缓存、checkpoint 或冻结结果。单次控制先行的耗时包括检查与 CPU 材料化，也有初始化和 allocator 历史影响，不能当作稳态加速基准。

入口 probe_full_window_encoder.py；输出 outputs/M28-A10-ENCODER-CONTEXT-PREFLIGHT。Slurm 脚本 M28-A10-ENCODER-CONTEXT-PREFLIGHT.sbatch：generanno，private-teodoro-gpu，1 RTX3090、2 CPU/16 GiB、10 分钟。申请时限不是可靠完成预测；OOM、超时或身份不一致保留失败产物，不自动换 GPU、dtype、分块、权重或延长时限。

## 判断与边界

八窗全部完成且原四块特征复现一致，只说明这一路输入表示的工程可行性。整窗模式同时改变跨块注意力与位置处理，不是单一因素因果证据；特征变化也不是性能提升。旧 A3/A6/A9 头按四块特征训练，不能直接把它们接到整窗特征上报告干净的性能对照。

本轮 optimizer 更新 0，无新头/头 checkpoint/参考指标，无 DEV/test/Setaria 前向，不改变正在运行的 FULL Q/C0、A6/A7 或 A9 fit/probe。若通过，下一步才是准备明确匹配的新特征/新头实验及固定协议；不自动启动新的首训、完整 DEV 或更改 A9 的 TRAIN gate。若失败，按具体工程原因判断该路线是否值得继续，不扫描参数制造 PASS。论文目标尚未完成。

## 执行记录

作业 13359555 已 COMPLETED/0:0，root/batch/extern 全部成功，Slurm 总时长 32 秒。STATUS 为 COMPLETED、JOBID 为 13359555；private-teodoro-gpu/gpu035，1 RTX3090、2 CPU/16 GiB，实际命令 `sbatch --parsable sbatch/M28-A10-ENCODER-CONTEXT-PREFLIGHT.sbatch`。提交前检查实际 GPU、运行/排队作业、预计开始时间与 reservation；私有 gpu035 有可用 3090，BeeGFS 可用约 159 TiB。三个新增文件的远程内容与本机精确一致，Python AST 与 bash 语法检查通过。

八窗全部完成，两物种各四窗，原四块 BF16 特征全部逐元素精确复现；相同 token IDs 的整窗输出均为有限 4,096×2,048。八份 feature artifact 均非空，各 16,803,676 bytes，总 134,429,408 bytes（约 128.20 MiB）；records.jsonl 八行与 summary.json 的逐窗记录完全一致，artifact 实际大小与记录一致。完整原始证据保留在 outputs/M28-A10-ENCODER-CONTEXT-PREFLIGHT，summary.json 为 11,229 bytes。

模型载入及移入 GPU 为 7.2382 秒，脚本测量区间总 15.1130 秒；不包含此前元数据读取或整个环境启动。整窗每窗同步前向/检查/CPU 材料化为 0.2337–0.2388 秒，峰值 allocated 2,522,416,640 bytes（2.349 GiB）、reserved 2,730,491,904 bytes（2.543 GiB）。原四块首窗 2.4888 秒，后七窗约 0.2594–0.2609 秒；顺序、首次初始化与材料化实现不同，不据此宣称整窗加速，也不直接外推训练/全 DEV 总时长。

整窗与四块特征的逐窗平均绝对差为 0.1344–0.2050，最大绝对差为 3.875–8.0625。它们确实不是同一输入分布，但数值没有注释正确性解释。本轮 optimizer 更新、backbone 可训练参数均为 0，无头 checkpoint/头前向、TRAIN catalog、DEV/test/Setaria 使用。

判断：整窗输入路线具备在同一 3090 上继续准备实验的工程条件，没有新的性能或泛化 PASS。下一步准备匹配的新特征与从头拟合的头，保持明确对照；不得把旧四块头直接接入新特征作为干净比较，不得混入旧 feature cache。此次预检已由主对话消费，不加入长期监控，也不重复检查这份已完成产物；既有 FULL、A6/A7、A9 流程与固定监控保持不变。记录沿用用户要求的工程事实、科学结论和下一决定分开表述；没有对外发布、Git 提交或 goal 状态变更。
