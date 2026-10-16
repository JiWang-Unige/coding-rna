# A28 原搜索淘汰诊断执行

2026-10-06。已固定a28_beam_search_trace_protocol.md并实现独立观测回放；只观察A27全部221局部通过但遗漏配对及324自然命中对照，所有4247个合格参考/query配对，不改原src decoder、旧训练、宽度、阈值、数据或checkpoint。

## 首次预检与有界修复

实际命令sbatch --parsable sbatch/M28-A28-BEAM-SEARCH-TRACE-PREFLIGHT.sbatch，job13392793，private1RTX3090/2CPU/16GiB/15min。提交前实际五项调度检查显示gpu034余1、gpu035余2块3090，没有维护截断或新目标目录碰撞；源码三文件AST和SBATCH语法通过。job root/batch FAILED1:0、externCOMPLETED0:0，10秒，gpu034。

原六项beam测试全部PASS0.291秒；新增五项中自然输出/carry、局部第五起点、其他活跃历史的全局入口淘汰和保留beam分数界四项PASS。失败在completion-pool合成fixture：它假设第二条短链会输，但原geometric normalizer及fixture EOS使该链实际胜出，first_loss为None。尚未运行真实八窗入口；没有真实模型或训练失败，也没有任何科学结果。

有证据最小修复仅调整这个fixture的首／尾评分，让第一条链明显胜出，测试仍要求原搜索／观测搜索相同，并检查第二条确实由完成池淘汰。不改trace实现、原decoder、真实数据、容差或排名规则。保留原失败目录／日志；唯一修复预检目录为outputs/M28-A28-BEAM-SEARCH-TRACE-PREFLIGHT-REPAIR-R1，SBATCH同资源、不同日志名M28A28TRACEPR1，runner只更新实际preflight路径，后续full必须校验这次真实成功预检的job和原全量domain。

首轮结束时，实际八窗预检和完整诊断尚未执行，不能据合成测试判断真实淘汰比例；下面记录真实R1预检和后续完整诊断，不把这个历史中间状态当作当前结论。

修复预检实际命令sbatch --parsable sbatch/M28-A28-BEAM-SEARCH-TRACE-PREFLIGHT-REPAIR-R1.sbatch，返回job13392831。提交前重新完成五项实际调度检查，private gpu034/035分别余1/2块3090、无GPU维护截断、新R1目录不存在；资源仍1RTX3090/2CPU/16GiB/15min。未重提原失败job或改变真实搜索。

## 实际预检通过与全量资源

13392831 root/batch/extern均COMPLETED0:0，50秒，gpu034，batch MaxRSS1666944K；原六项测试PASS0.133秒，新增五项PASS0.064秒。真实八窗／13参考配对／86参考-query profiles／64 window-query回放完成28.304秒，峰值GPU148045312字节。每query原搜索几何、次序、累计分数及全部计数一致；实际完成链属于原beam4池，正确前缀分数／几何一致，七个自然命中对照均返回，保存记录重读聚合核验PASS。六个遗漏配对仍遗漏；不据这些八窗的stage比例作全量判断。

实际原metadata决定全量352窗、545参考配对（221遗漏＋324对照）、4247参考-query profiles和3363不同window-query回放；同一个window/query只算一次搜索，但同时与原函数执行两路对照。64回放／28.304秒保守按3363放大约24.8分钟，另考虑不同历史／链长和操作余量，估25–35分钟、固定上限45min。这是执行估算，不是稳态benchmark。

全量SBATCH M28-A28-BEAM-SEARCH-TRACE-DIAGNOSTIC.sbatch，generanno、private1RTX3090/2CPU/16GiB/45min，原唯一全量目录不变。文本预计低于64MiB、无模型或checkpoint；full参数显式绑定预检job13392831并读取其STATUS/JOBID及实际full_domain，不依赖可能已过期的Slurm preflight dependency。原失败13392793及所有旧结果保留，已修复R1不再重提。

## 已启动的完整诊断和同一监控

实际命令sbatch --parsable sbatch/M28-A28-BEAM-SEARCH-TRACE-DIAGNOSTIC.sbatch 13392831，返回job13392879。提交前第三次五项实际调度检查均完成，private gpu034/035分别余1/2块3090、无维护截断、全量目录不存在。新SBATCH语法通过；首次调度取证RUNNING gpu034／45min，root/batch/extern运行，输出JOBID13392879及STATUS RUNNING匹配。没有取消或重训任何旧实验。

A23fit13391608实际COMPLETED0:0／1h04m03sec；A25fit13392060实际COMPLETED0:0／59m46sec。两自然TRAIN读出13391609／13392061在gpu034运行，下游CPU13391610／13392062有效Dependency。两候选还没有完整结果，不能从fit完成推断三条件门槛通过。

同一个监控对话01a0f71f-8845-7432-bd13-bf54c1808fb3和heartbeat coding-rna-r5-light-monitor已更新为M28-A23-A25-A28-TRAIN-20261006，仅原A23三job、A25三job及A28单job。A23/A25结构白名单原样复用；A28按352窗／545 references／4247 query_profiles／3363 natural_queries及实际原搜索／分数／对照核验flags检查结构，不读阶段比例或解释机制。A27和全部预检已由主对话消费，移出监控。

原生保存2026-10-05T18:11:48.699Z／updated_at1791223908699，完整prompt/id/name/kind/target/cadence/notification/created_at精确读回PASS，created_at1790981930511保留；gpt-5.6-luna/low及15分钟不变。已请求同一对话manual startup；配置、启动检查、定时运行及终态自动回传分别记证，不借用旧cohort的定时记录。

当前cohort manual startup实际验证：取证2026-10-05T18:12:47.384Z，turn01a10d44-6374-7350-a632-e08ada2fcbe2 completed/error-null／31.453秒。A23/A25 fit终态0:0、两读出及A28正在gpu034运行，下游CPU有效Dependency；未发送健康回执或终态交接。实际turn_context18:12:28.183Z仍gpt-5.6-luna／low、相同cwd。当前定时运行和完整自动交接尚未观察，不能从启动成功宣称已验证。

全量job13392879已完成：root/batch/extern均COMPLETED0:0／17min39sec／gpu034，batch MaxRSS1662484K；实际1046.186710秒、峰值GPU148045312字节。545 references／4247 query_profiles／3363 natural_queries实际行数及全部原搜索／分数／对照／保存记录flags闭合。221局部通过但原遗漏配对仍全部遗漏，324自然命中对照均至少一个query返回。详细阶段计数、解释与分母边界见a28_beam_search_trace_result.md。

本结果已由主对话直接读取。A29 job13393045也已完成同一原cutoff监督资格标注，结果见a29_search_cutoff_result.md；两终态已直接消费并移出同一监控，仅保留A23／A25原六root。下一模型选择须结合完整A29与A23／A25读出，不据阶段比例改旧FAIL门槛。完整raw-FASTA原生竞争、独立测试、机制因果消融及论文创新目标仍未达成。
