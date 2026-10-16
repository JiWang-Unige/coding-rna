# Carbon-A 原 TRAIN 窗口策略诊断

2026-10-12 UTC。独立的原六 TRAIN 诊断已完成：短窗口的过滤／未过滤pooled整链F1为77.562301%／75.585179%，较默认分别下降0.609254／0.478270个百分点；六条染色体两口径均下降。完整终态和科学边界见[结果记录](carbon_a_context24_train_result.md)。不改98,304 bp默认基线、不改B/C/D/E训练或冻结DEV家族。以下保留先定对比与提交记录，不因观察结果选阈值或重新定义基线。

## 对比与边界

已完成原生Carbon-A在固定confidence≥0.1或NA口径下pooled整链F1为78.171555%，未过滤为76.063450%；已完成joint 24,576 bp候选为66.962764%。来源为[原生评价](../../outputs/M28-CARBON-A-TRAIN-EVAL/summary.json)和[joint评价](../../outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC/summary.json)。既有配对诊断显示净恢复差距主要落在多CDS参考，超过24,576 bp的参考仅9条；这不能排除长窗口对较短基因的影响，也不能确定上下文是差距原因。

| 条件 | 上下文tokens | 上下文bp | 重叠tokens／bp | 主窗口stride bp | 重叠比例 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 已完成原生默认 | 16,384 | 98,304 | 1,024／6,144 | 92,160 | 6.25% |
| 独立短窗口控制 | 4,096 | 24,576 | 1,024／6,144 | 18,432 | 25% |

固定原生revision f6789ee79e6ef397de333577ffec1cb734496620及共享backbone来源b0483c23b6b63787b61a6d3a204a9b517d6ba345。使用作者原生`model_inference.py --context_length 4096 --overlap_length 1024`生成新概率，再由原生`run_annotation.sh --skip-inference`运行同一transition TCN、GPU Viterbi、confidence TCN和导出；作者资产和权重均不修改。batch8、dispatch4、BF16／FA2、标准密码表、最短record100,000 bp、最短CDS96 bp、record边缘内含子20,000 bp和chunk500 Mbp均保持。复用已资格核验的隔离软件环境和node-local Transformers归档，不升级四条训练的环境。

只有六个既有整染色体FASTA：NC_003070.9、NC_003071.7、NC_003076.8、NC_089043.1、NC_089044.1、NC_089046.1，合计156,779,632 bp。推理不读参考、标签、旧预测或feature cache，不读取DEV、sealed或Setaria，不更新参数。新输出独立于原默认管线，特别是NC_003071.7必须接受新控制自己的收据，不能继承原98 kb的CPU恢复收据。

## 预先固定的评价

新六条raw均完成及通过逐字DNA、双链概率有限／范围、坐标／CDS／phase／字面密码表蛋白校验后，才运行原共同尺。保留22,688条primary参考，采用M26 longest eligible CDS／transcript ID tie政策；unsupported、long和overlap参考不因训练资格筛除。保持confidence≥0.1或NA和未过滤两个点，均报告，不择优。

每染色体、每物种、pooled、species macro及worst均保留。主要量为24 kb减98 kb的整链F1百分点，并列两端TP／唯一预测／参考及PRF；另列24 kb原生管线与已完成joint的描述性差距。固定绝对重叠导致相对重叠、窗口位置和概率平均几何共同改变，结果属于窗口策略效应，不是上下文长度的独立因果效应。与grammar的对比还混合了不同head、辅助网络与过滤；不声称解释全部差距、独立泛化或公平速度优势。

## 实施与停止条件

新raw、评价、比较输出分别为`outputs/M28-CARBON-A-CONTEXT24-TRAIN-RAW`、`outputs/M28-CARBON-A-CONTEXT24-TRAIN-EVAL`和`outputs/M28-CARBON-A-CONTEXT24-TRAIN-COMPARISON`。raw为6任务、最多6并发，每任务1 RTX3090／4CPU／24GiB／90分钟；原98 kb六任务实际最长约34分钟，较短窗口的额外重叠保留在预算内，90分钟是新诊断的待验证时限而非速度结论。评价2CPU／8GiB／20分钟，比较1CPU／1GiB／5分钟，无GPU。分别使用whole-array afterok和评价afterok；不取消或重提四条既有训练。

提交前核对现场GRES、队列、maintenance、并发限制、独立输出路径及存储。任何原生运行或校验失败先保留证据并诊断；指标下降是结果，不因指标失败调阈值或换操作点。只更新原固定监控对话和原心跳的监控对象，不创建新的对话或自动化。作业ID、实际Slurm终态和完整产物在获得后补入本文件。

## 提交与工程记录（历史）

已提交`sbatch --parsable sbatch/M28-CARBON-A-CONTEXT24-TRAIN-RAW.sbatch`，返回数组13507459；六任务实际均为RUNNING，0–3位于gpu034、4–5位于gpu035，实际UTC开始均为2026-10-12 02:50:46。TZ=UTC sacct实际任务0／1／2／3／4／5的numeric root分别为13507460／13507461／13507462／13507463／13507464／13507459，非按parent递增猜测。评价13507526使用`--dependency=afterok:13507459`，比较13507618使用`--dependency=afterok:13507526`；03:06:14UTC已核实二者PENDING／有效Dependency，scontrol亦确认比较的afterok来源。不是完成或准确率结果。

本机与远端均通过三项标准库工程检查：默认调用保留、4096／1024 tokens单位，以及六条收据隔离（含原NC_003071.7恢复）和错误条件拒绝。Python／SBATCH语法通过；两个已完成固定点的same-source schema rehearsal为零窗口差，检验比较器能读取真实旧结构，不是新24 kb推理。实际新native概率、CDS／蛋白及共同尺仍待作业完成。

2026-10-12 02:42:47UTC现场调度／存储读回：私有gpu034配置8／已分配4、gpu035配置8／已分配2，共10空闲GRES；当前展开作业53，追加8个不超过80；相关GPU维护窗口未见冲突，旧六raw产物528 MiB，所在文件系统可用146 TiB，三个新根未存在。因私有匹配RTX3090可立即运行，沿用同型号／同软件，未变更四条训练。原心跳已更新新增L对象，保存读回确认type、ID、原目标对话、间隔、通知设置和创建时间不变，完整旧17431字符前缀保持。首次新L范围定时执行已验证：原监控对话中的03:05:18–03:05:46UTC自动轮次使用gpt-5.6-luna／low，实际输入含完整20154字符保存提示，03:05:32UTC实际SSH核对四fits、六Lraw均RUNNING且全部下游Dependency有效，未发送交接或暂停；新L终态自动回传因作业未结束而尚未验证。
