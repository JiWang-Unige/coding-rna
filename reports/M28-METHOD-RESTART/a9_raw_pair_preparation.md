# A9 无参考 raw-FASTA 配对调用准备

2026-10-03，A9拟合13359395仍在运行、TRAIN probe13359396为有效Dependency。原首拟合/beam/阈值/继续条件不变；本轮只准备后续raw调用及CPU合成测试，不读取真实DEV，不启动全DEV推理或评价。

`src/m28/query_prefix_pair.py`实际生成early critic-only beam4和late原greedy+beam4，后者调用原DNA decoder，不用参考或已存DEV缓存。两臂(query,geometry)union在同一末头做一次新的完整STOP readout，再映回各臂；所有共同query-chain值相同。每窗各臂对geometry显式max over query，保存query来源，不把query重复当额外gene；非正值仍在自由候选池，原selector的null=0只在组装时生效。

early逐entry/end价值剪枝且无错误raw-score提前bound；late保留原pointer/EOS beam的提前bound和额外greedy来源。候选曝光、搜索工作量及query集合不同，不称唯一变量或等budget/速度因果对照。共同query-chain值相同不保证两臂某一共同geometry的最大值相等，max的query集合可变；不是自然gene正确概率。

`infer_query_prefix_pair.py`提供实际的单contig调用函数：对两方向原natural grid逐窗，由同一native GENERanno encoder计算一次BF16重现的features和一次冻结generator hidden，两臂共享，不另读feature cache/annotation。DNA大小写先规范化、其他编码为N，reverse complement和半开坐标转换沿用原实现；各臂复用nearest owner、global geometry dedup、positive-value nonoverlap及原CDS phase导出。逐窗pair JSON同时保存两池和共同值，组装分臂读同一记录文件。

这不是全DEV执行入口：脚本直接执行会报“Preparation helper only”，不自动载入权重或数据、不自行确定数据范围/预算。只有原TRAIN probe最低机制条件通过后，才在任何新DEV前向前固定末checkpoint、完整8690窗/7728主参考、资源/停止条件、paired评价与继续条件，并加正式Slurm入口。若probe失败，不能用本准备绕过它或在DEV救本头。未增加optimizer、新GLM真实数据forward、test/Setaria访问、公共发布或Git提交。

CPU检查新目录`outputs/M28-A9-QUERY-PREFIX-PAIR-CHECKS`，private/0GPU/2CPU/8GiB/5分钟、generanno。5个测试用真实小query decoder参数及合成DNA/encoder替身，检查late原池精确重放、共同readout一次、query重复max、zero值不删自由池但不输出gene、空DNA语法池、softmask/padding及raw helper两方向单次forward/负链坐标和phase。它们仅验证接口和实现，不是实际GENERATOR/GLM raw-FASTA运行或准确率、性能、泛化/论文证据。

## 实际验证

五份新源码/测试/sbatch/准备记录部署后逐字符回读一致，AST与bash语法通过；仅提交一次`sbatch --parsable sbatch/M28-A9-QUERY-PREFIX-PAIR-CHECKS.sbatch`，返回job13359462。root/batch/extern均COMPLETED/0:0，实际16秒、2CPU/8GiB/0GPU（0.008889独立CPUh），batch MaxRSS724388KiB；5 tests PASS/0.488秒，无重提。没有修改已在运行的A9 fit/probe源、权重、参数或依赖，原FULL/A6/A7也不变。

同一固定监控在12:57:52.234UTC的后续scheduled检查仍读取这11根原作业，turn01a101d7-6f86-7233-9091-a9916d4142f9已正常完成、无主对话交接或普通状态消息；FULL/A9 fits与A6 infer运行，剩余下游有效Dependency。这是等待取证，不是完成或指标结果。已消费的CPU准备检查不加入监控目标；未来仍须先核对A9 TRAIN probe原继续条件，再决定完整raw协议和执行。
