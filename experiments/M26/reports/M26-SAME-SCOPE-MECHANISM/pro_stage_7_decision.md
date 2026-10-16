# Pro第七轮商议与R8执行

2026-09-26。内置浏览器沿同一“科研仓库系统审阅”会话，第七轮最终回复思考6m52s，明确ACCESS_VERIFIED：读取固定提交1518358eb357f9680ad7f57e3f6cf8c17f069221的semantic_reference_pair_result.md、semantic_reference_pair.json及run_graph_core.py::candidates。本轮未重跑；不是独立同行评审。

## 采用的判断

R7应结束本DEV来源争议主线；不能为0个残余结构争议再补验证。优先做一次固定缓存、等端点预算的离散化损失检验，不马上开展只有最终GFF的跨物种排名。候选过滤损失也可能只是本实现问题，不能据此宣称NC级普遍机制。

Pro建议原图、按全染色体同事件边界logit前M、窗口内置换3臂；新增保留原端点但解除run约束的unlocked归因臂，目的是分开“端点选择”和“结构组合放宽”。所有科学条件在执行前记入endpoint_budget_contract.md。新增对照不使用参考选择端点、不调参、不训练。

主投入线是相对原3,561可达链净增至少323、两物种同向并超过置换；另要求区域/边界联合支持增加。不是期刊门槛，不把一次置换当显著性。报告片段和几何连接机会规模，后者明确为不含phase兼容的上界。

## 已核对的已有方法边界

本轮直接读取Tiberius原论文（Gabriel等，Bioinformatics 2024，doi:10.1093/bioinformatics/btae685）。其以逐位置概率进入结构化HMM，并将可微HMM用于训练。因此“软分数+生物约束结构解码”不是这里可自动认领的新贡献。R8只能判断当前候选策略是否值得进一步研究，不构成文献空白证明。

原文：https://academic.oup.com/bioinformatics/article/40/12/btae685/7903281

## 执行记录

命令：sbatch sbatch/M26-ENDPOINT-BUDGET.sbatch
Slurm job：13227158；private-teodoro-gpu；4CPU/16GiB/0GPU；2小时硬上限。
环境：/opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh → generanno，沿R3实证环境。
输出：outputs/M26-ENDPOINT-BUDGET-R8/。
提交前已读cluster_config，检查队列、预约及输出空间；无资源/维护硬冲突。提交后R运行于gpu034。本文件先记录启动，终态以endpoint_budget_result.md及Slurm记录为准，不能把提交视为结果。
源脚本、测试、sbatch与合同已在运行目录保留执行前副本。未触动其他TE作业、原始输入或历史结果。
