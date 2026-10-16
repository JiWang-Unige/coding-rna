# A42 紧凑势存储工程准备

2026-10-07。本分支处理完整 contig 解码的实际内存问题，不新增科学机制。原 `_Margins` 将 15 列 gain、15 列 prefix 和 3 列 coding prefix 转成逐项 Python 浮点列表；短窗可运行不能证明真实染色体内存可行。新版本只将这三组持久数组保留为 FP64 NumPy 视图，保留原 Torch 求差和 cumsum、原操作次序，并在所有标量运算前显式转回 Python float。

## 实现与证据范围

`src/m28/emission_lattice_compact.py` 是既有 cached metadata 分支的隔离副本。除上述存储及标量取值外，DNA、ORF、全部 motif 合法边、实际未完成密码子、exon cap、状态更新、严格大于比较、先到平分、回溯及 20 项计数保持原实现。三组持久数值数组合计 `264 × N + 144` bytes；这个数不含输入 FP32 势、临时 tensor、事件索引、状态、路径和边缓存，不是整个解码器的内存上界。

原 A42、A44 和 raw 接口仍调用各自的冻结解码器；不接线、不替换，不调整阈值、角色权重、exon cap 或任何研究停止条件。原 A36、A42、A43 的科学 FAIL 不变。本分支只用人工势和已闭合的人工回放，不使用 TRAIN、DEV、独立物种、模型或注释 catalog。

## 固定人工检查

新增数组检查比较零值、固定随机 FP32 值及大动态范围 FP32 值的全部 FP64 位模式和 exon 增量，验证 Python float 运算及持久数组字节数。新解码器复用实际短 carry、GT／GC、intron N、负前缀、相邻 gene、零值 null、浮点平分、整数与非整数、128-exon cap、30 kb 完整 exon 及正反上下文 fixtures，逐项对照原始及 cached metadata 答案；不为未修改的旧分支另起回放。

同一个 Slurm 作业中，检查成功后先用独立子进程分别构造百万行全零 CPU FP32 的原存储和紧凑存储，记录构造前后 RSS、进程峰值 RSS、构造秒数和紧凑数组字节数。RSS 包含 allocator 和导入开销；时间仅为单次工程记录，不作公平提速比或真实染色体可行性主张。

随后在原 49,152 bp、DNA seed 0、三个窗口、plus strand 的两个固定案例上运行紧凑解码：全零势，以及原 potential seed 41 的 FP32 整数 −3 至 3。上下文所有权、FP64、单位角色、零 null／threshold、cap 128 不变。每个自然预测先写完关闭，再读取原 cached metadata 13423648_0／1 的保存答案；完整路径及顺序、逐链／总 gain、所有权和全部 20 计数须逐项一致，另独立重算 DNA／ORF／motif／role／cap／span、整数 FP64 分数和 pointwise 上界。该证据仍只是人工同答案工程验证，不是注释性能或新训练结果。

## 资源与停止条件

固定入口为 `sbatch/M28-A42-COMPACT-STORAGE-CHECK.sbatch`，generanno／private-teodoro-gpu／0 GPU／2 CPU／8 GiB／6 分钟／no-requeue，预计约 3 分钟；唯一根 `outputs/M28-A42-COMPACT-STORAGE-CHECK`、唯一 `M28A42COMPACT_%j` 日志，预计新增输出小于 8 MiB。提交前核对活调度、维护窗口、并发限制和空根。两项百万行测量只构造 margin，不执行百万 bp DNA 搜索。

只运行一次；明确实现错误可保留原输入、资源和全部口径，用独立 R1 根最小修复一次。答案或数值不一致、TIMEOUT／OOM 或第二次失败停止，不扩大内存／时限，不重跑训练。本工程短作业不增加研究监控对象或创建新的监控对话；原固定 A44 监控继续等待其 70 个流程根。

## 实际完成记录

10:04 UTC 五项活调度／预约、空根与唯一日志检查通过，BeeGFS 可用 156T；10:05:05 UTC 账户展开 72 项，加本次为 73，小于 80。新五文件及旧 cached metadata 源与 canonical 逐字一致；三个 Python AST、SBATCH bash 语法和新旧 DP 的规范化 AST 对照通过，DP 仅有声明的 margin 构造与两处 Python float 取值变化。实际命令为 `sbatch --parsable sbatch/M28-A42-COMPACT-STORAGE-CHECK.sbatch`。

13426626 于 2026-10-07 10:05:06 UTC 提交，10:05:07–10:08:09 UTC 完成；root／batch／extern 均 COMPLETED 0:0，STATUS／JOBID／UNIT_STATUS／summary 相符。实际 private／generanno／gpu035／0 GPU／2 CPU／8 GiB／6 分钟／no-requeue／Restarts 0，Command／SubmitLine／WorkDir 与唯一日志相符；batch 采样 MaxRSS 为 2,245,672K，不能与下列单个子进程峰值混用。8 项检查 PASS，用时 21.608 秒，无失败或 R1。

百万行势的两个独立子进程实际观测如下，单位均为 KiB。这里的 RSS 包含导入及 allocator，紧凑数组精确数值字节数另为 264,000,144 bytes。

| 存储 | 构造前 RSS | 构造后 RSS | 进程峰值 RSS |
| --- | ---: | ---: | ---: |
| 原 Python 列表 | 450116 | 1836072 | 1997732 |
| 紧凑 FP64 数组 | 448784 | 790892 | 835876 |

原子进程峰值约 1.91 GiB，紧凑版本约 0.80 GiB。这支持本人工输入的存储改善，不支持整个解码器或真实染色体峰值外推。构造秒数分别为 1.966415 和 0.218010，属于单次执行事实，不作公平提速比。

两个原 49,152 bp 完整自然回放均通过四项证明：源配置与 13423648_0／1 身份、先关闭新自然输出再读旧答案、全部路径及顺序／FP64 分数／20 计数一致、DNA／ORF／motif／role／cap／span 与直接整数分数独立核验。零势为 null／0 gain；整数势为 88 条链／21526 gain，实际解码秒数分别为 61.676628 和 73.318548。完成后主核对了两个四证明、完整计数、非空自然记录和所有无数据／无模型／零 update 标志，不重跑已通过的探针。

本分支保留为隔离工程准备，没有替换任何冻结研究或 raw 调用。A44 的实际 TRAIN 结构结果仍等待其原完整流程；本次内存改善不新增科学 PASS、独立泛化或论文性能主张。
