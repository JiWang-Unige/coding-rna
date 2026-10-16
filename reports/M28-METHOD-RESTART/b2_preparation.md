# M28 B2 前缀兼容性实现与首训准备

记录日期2026-10-01。R6按事前规则给出两物种真实剪枝见证，详见[r6_trace_result.md](r6_trace_result.md)。下一目标是可训练新模型的完整基因质量；本准备不是新的DEV性能表，不改变R5，不保证Nature Communications创新成立。

## 单一方法改变

保持B1最后step004608的共享局部表示、15类分割、端点和连接头全部冻结；读取原GENERanno revision及缓存，不继承受DEV选择影响的新权重，不改endpoint、link、ORF、条件化exon预算、owner或最终gain>0选择。

新增PrefixValue：按有向exon特征、间距及相位形成GRUCell前缀状态，结合当前entry/donor/terminal事件的局部特征与几何，输出专门训练的绝对兼容性logit。它不是旧累计proposal加小残差。空起点可与已延续长链前缀在同一entry竞争；完整链仍以0为null操作点。新null/logit需训练校准，不能解释为自然候选正确概率。

首个位置机制对照共享同一个新增训练checkpoint：early在beam及最终固定上限前用PrefixValue；late仅在原加性候选池形成后用同一PrefixValue作完整链读出。两种解码接受完全相同的新增训练权重和曝光，不独立择epoch。这样比较评分位置，不把不同优化路径混入首个问题；两种解码不是独立复现。late自由候选应逐身份重现原B1固定池。

## 训练监督与自由搜索隔离

第一拟合拟用原1,536 TRAIN paired draws/importance weight和R4标签政策，seed0；仅新前缀头可训练。新增完整链BCE与首次实际剪枝的softplus配对排序，各窗归一后系数1，乘原importance weight。首个训练准备先做真实吞吐与梯度测量，更新数、checkpoint和最终资源将在任何optimizer更新/新DEV性能之前写入完整合同；当前不提交正式拟合。

用于监督的FirstCutObserver只能读取实际自由搜索的组，不影响候选、打分或排序。对于真实参考前缀第一次实际被beam或最终上限裁掉，重新计算其可微logit，与同事件/carry实际存活的已知不一致前缀比较；beam用实际两名保留者，final节选前两名及cutoff处保留者。参考路径若先受端点/exon/连接/语法阻断，不以注入伪装真实自由路径。

完整链正例仅按既有R4 primary/支持政策。其他isoform精确前缀、未知或非callable跨度仍标-1，不作为负例；prefix未匹配已知完整primary并不自动等于生物学错误。完整正链只注入损失池，free池及来源单独保存。背景空池输出应为可微零损失，不虚构参考正例。前移推理不接受参考参数，也不假称离散搜索可微。

当前首个拟合将用加性自由搜索产生固定监督事件，属于off-policy早更新式前缀学习；改为on-policy不是本首个对照的额外同时改变。若不能改善自由池/owner及完整P/R，不能靠训练注入表现或增加K来救解释。

## 准备运行与边界

先进行合成CPU验证：控制路径与原generate一致、同Budget优先级改变可救合法路径、缓存/梯度/冻结边界、isoform未知掩码、首次真实剪枝与参考注入隔离、自然背景空池。均是会改变是否进入真实梯度测量的具体条件。

拟定真实准备只取原TRAIN draws每物种/染色体/链向stratum第一条，共12窗；清单必须在checkpoint/特征读取前保存。用新随机seed0头、同B1冻结局部来源，测量加性控制搜索、新前移搜索、完整链与前缀loss反传及显存/用时，不调用optimizer.step，不产生拟合checkpoint或DEV排名。共享catalog/index只按split标记选取TRAIN元数据和该窗口标签，不解析DEV参考值，不读取DEV特征或进行DEV前向/性能分析；test/Setaria不访问。

CPU准备最多2CPU/8GiB×5分钟；无GPU。GPU准备若CPU验证通过，单RTX3090/24GiB、2CPU/16GiB×15分钟，最多0.25 allocated GPUh；不复制既有159GiB缓存，不保存DNA/特征/模型权重，compact产物≤32MiB。失败/非有限/冻结参数有梯度/标签或窗口错配均停止，不自动接入首训。上述准备不能冒称完成训练；正式拟合须另冻结实际可完成预算。

## 后续判断不变

完成配对全DEV之前不形成新方法成绩。预期判断覆盖、owner稳定性、每物种exact-chain P/R/F1、背景span/CDS覆盖及成本；不能用注入池改善代替自由收益，也不能因总F1较好隐藏一物种退步。是否独立迁移仍需完整方法证据与单独冻结的数据/权重版本；Setaria继续封存。
