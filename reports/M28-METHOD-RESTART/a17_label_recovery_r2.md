# A17 完整候选池的标注恢复

2026-10-06，恢复标注提交前固定。A17的R1拟合13387910已完整完成；自然池生成读出13387946在生成全部1491窗后，首次重算原A3标签时发生JSON嵌套list不可哈希的实现错误。现有完整自由池可以恢复标注，无需重新训练或解码。这是工程失败，尚无A17科学gate结论。

## 实际错误和唯一修复

2026-10-05 16:17:39 UTC直接取证：13387946 root／batch FAILED1:0、extern COMPLETED0:0，54分38秒；STATUS为FAILED exit=1、JOBID13387946。free_pool.jsonl为1491行、1565653字节，windows.jsonl为0行，无summary。错误在run_query_local_context.py:210调用candidate_targets，labels.py:119报TypeError unhashable type list。原CPU核验13387978未执行，实际DependencyNeverSatisfied；不取消该作业、不覆盖原失败目录或拟合结果。

候选来自实际JSON输出，转换外层tuple不能使内层外显子list可哈希。只将candidate_targets的候选身份从tuple(chain)改为tuple(tuple(exon) for exon in chain)。原build_targets产生的正链／protected集合、精确整链身份、unknown跨度和callable规则均不改；原运行时tuple链转换后的身份相同。两项新检查分别验证JSON往返后的1／0／-1政策及同外跨度不同外显子结构不自动变成正链。没有新标签、阈值、loss、参数或数据政策。

## 仅恢复标注和独立核验

恢复只读末R1拟合的身份摘要，不加载模型做前向。核对旧失败JOBID、失败前0条标注、无summary、全1491原顺序自由池、同原TRAIN几何、参考未进入生成及自然链去重。把原1.56MB自由池直接复制到独立新输出并直接比较记录身份；读取原TRAIN cache DNA和catalog重算原A3以及本池标签，保留完整新旧候选数，按原三项继续条件计算。没有新模型／GLM前向、搜索、optimizer更新、参考注入或seed／epoch变化。

独立CPU核验沿用原verify_query_local_context_train.py，仍直接核对末模型结构／有限权重、原4608 draw／importance／A3监督身份、全部1491新旧候选的实际DNA标签与原三项gate。该核验的只读CPU模型构造不是模型前向。来源标注job改为真实R2作业，生成池来源仍明确为原失败读出13387946；两者不混为新生成成功。

独立输出为M28-A17-LOCAL-CONTEXT-TRAIN-READOUT-LABEL-RECOVERY-R2和M28-A17-LOCAL-CONTEXT-RESULT-CHECK-LABEL-RECOVERY-R2，保留STATUS／JOBID／summary／全部记录。两阶段均generanno、private-teodoro-gpu、0GPU、2CPU、8GiB、30分钟，afterok链；仅元数据／标注小文件新增，已有约50GiB缓存只读、不复制。原八小时拟合、四小时GPU读出上限不进入本次CPU恢复。恢复实现失败先报告，不无限重试；科学FAIL仍关闭A17，不能救参。PASS只支持另固定完整caller／评分协议，不自动运行DEV、sealed、Setaria或发布。

旧监控批次因失效依赖已由主对话直接消费并原生暂停，实际kind／固定target／15分钟／failed_runs_only／created_at及原prompt保留。新作业提交后只更新同一固定监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3和同一heartbeat coding-rna-r5-light-monitor；不再把13387946或13387978列为活动待闭合作业。论文模型整体目标继续，A23新生成结构与此格式修复独立。无Git提交／推送、作业取消或公开发布。

## 已提交的实际恢复

标签恢复13391371、独立核验13391372 afterok:13391371；提交后实际分别PENDING None和有效Dependency。后续实际标签恢复已RUNNING，两项回归检查PASS，用时0.002秒；标签仍逐窗处理，独立核验尚待依赖。这里只记录提交与工程进度，不给出科学PASS或FAIL。
