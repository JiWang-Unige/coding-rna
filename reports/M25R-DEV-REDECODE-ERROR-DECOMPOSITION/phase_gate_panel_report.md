# M25R epoch-1 phase 门控消融：面板冻结与预检终态

日期：2026-09-08。

## 结论

**本次已冻结128条lineage，但需137个原方向化窗口，超过用户批准的64窗口硬上限。执行在GPU推理前停止。** 这是预算可行性阻塞，不是消融科学失败，也不是phase路线的否决。没有提交GPU作业、没有原路径新推理重放、没有运行phase旁路，因此不存在可报告的消融恢复率；相关结果应标记NOT_RUN，而非0。

用户已将32/64与8/32改为探索性参考，不是总体性能gate或整条研究路线的否决标准；本报告及新增预检实现遵守该修订。旧回顾性报告提出的硬门槛/停止phase优先路线表述不再作为执行规则。

## 已完成工作

- 按明确请求，通过内置浏览器向现有ChatGPT研究会话发送一次聚焦审阅：要求读取私有GitHub仓库 `JiWang-Unige/coding-rna` 的基础提交 `01fbbf15544f18ef0fca19688adc1aa0d3284248`，审阅消融解释、抽样及探索阈值。明确最新CPU结果尚未push，ChatGPT只能将它们当作提供的Baobab事实，不能声称亲自读取。没有额外调用其他评阅者或旧CLI评阅框架。
- 读取当前Baobab AGENTS、cluster_config、原R4实现和执行记录，保留已有未提交工作及旧框架退役删除。依据coding-rna-baobab Skill，真实面板核算在Slurm CPU作业中执行，登录节点只做源码、元数据和调度读取。
- 新增 `scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/phase_panel_preflight.py`、`tests/test_m25r_phase_panel_preflight.py`、`sbatch/M25R-E1-PHASE-PANEL-PREFLIGHT.sbatch`。预检只用Python标准库读取既有JSON/JSONL；不读FASTA、checkpoint内容或Setaria文件，不导入模型。
- 作业 `12499507`：`sbatch --parsable sbatch/M25R-E1-PHASE-PANEL-PREFLIGHT.sbatch`，private-teodoro-gpu，2CPU、16GB、10分钟上限，无GPU请求，Requeue=0。
- Slurm COMPLETED、exit0:0、elapsed00:00:05，MaxRSS163936K；5项聚焦测试通过（0.04秒），覆盖配额守恒、numeric-ID确定性选择、竞争motif位置及acceptor+2、方向分离、超64时不截断样本/窗口。Bash语法检查通过。没有重投。
- 正常完成CPU预检与检测到预算阻塞并不矛盾：业务终态为 `STOP_WINDOW_BUDGET_EXCEEDED`，不是作业失败。

## 冻结与抽样

选择规则在查看窗口覆盖和缺失scores前写定：每组按物种/方向比例分配；exact-phase组再按首次失败位点token-conflict有无分层；采用largest-remainder整数分配，平手按stratum tuple排序；层内按numeric lineage ID升序取前若干条。身份主键为species/seqid/strand/lineage_id，不能跨方向只用一个serial ID。

原方案明确是稳定ID前缀选择，本次没有悄悄改成随机、按长度或按窗口紧凑度优化的面板。`frozen_panel.json`在计算窗口上限之前以exclusive-create写出；窗口超限后没有改采样、替换、减少样本或截断必要位置。

| 组别 | 原eligible数 | 冻结数 | 所需原窗口（各组间重叠，不能直接相加） |
|---|---:|---:|---:|
| exact chosen链、原phase拒绝 | 1149 | 64 | 79 |
| nonexact chosen链、原phase拒绝 | 2129 | 32 | 37 |
| 原emitted controls | 2098 | 32 | 31 |
| 去重合计 | — | 128 | **137** |

| 物种/方向 | exact-phase | nonexact-phase | emitted control | 去重窗口 |
|---|---:|---:|---:|---:|
| Arabidopsis + | 28 | 10 | 11 | 51 |
| Arabidopsis - | 26 | 12 | 10 | 54 |
| rice + | 5 | 5 | 6 | 15 |
| rice - | 5 | 5 | 5 | 17 |

exact组token-conflict共25条、无conflict39条。完整stratum eligible数和选中身份随冻结文件保存。

**抽样限制：**稳定ID顺序关联方向化基因组位置，属于确定性便利面板，不是随机总体样本。物种/方向/冲突分层不能消除位置、CDS数目、链长、相邻区域相关性或富集组别带来的偏差。nonexact指对当前primary-reference集合无完整exact匹配；特别是unassigned/off-reference lineages，不能直接称为生物学false positive。各组恢复/新增输出数即使将来获得，也不得拼接成总体precision、F1、FPR或gene-count ratio。不能用常规二项置信区间假装覆盖了这种选择偏差。

## 窗口覆盖定义与消融解释

对每条冻结lineage，所需GPU位置包括：

1. `candidate_positions`里所有start/stop/donor/acceptor候选坐标，用于原 `_select_motif` 选择重放，不只是最后选中者。
2. 所有chosen边界位置，用于四类boundary score minima与原阈值比较。
3. 所有chosen CDS starts：first start与每个chosen acceptor+2，用于phase argmax与累计长度期望；不只保存的first-failure prefix。

将每个位置p映射到原方向化窗口起点 `floor(p/6144)*6144`，以species/seqid/strand/start去重；正反方向是不同输入，不能共用一份推理。未来推理必须复用原完整6144bp上下文、末端padding、tokenization及dtype。原R4 `predict_sequence` 无重叠窗口平均。

源码中的后置boundary score是chosen motif分数按类型取min；ORF检查只需连接chosen CDS的真实序列并检查长度、frame、ATG、终止密码子和internal stop，不需要额外GPU窗口。以上137窗口已经是固定候选条件下所需score位置的核算，不是盲目覆盖整染色体。若要额外证明完整region-state生成路径也重放，需扩展位置并重算，不能假称本核算已覆盖完整pipeline。

可解释的干预只回答：**对这批已经存在且固定的候选及边界，去掉phase否决这一条件，会有多少通过其余原检查？** 它不证明phase-head架构/6mer别名冲突是因果根因，不测grouping修复，不代表重新运行完整decoder后的全局性能。先复现原motif选择、已保存phase prefix和controls，再分析A/B差异，是将来执行的必要要求；本次没有GPU，因此这些尚未验证。

## 资源预检

历史同checkpoint/原推理路径的R4 job12270926，在RTX3090 24GB、8CPU下完成，sacct batch MaxRSS=13833M，gpumem accounting=2818M；epoch1 LoRA/head checkpoint文件元数据大小18311932 bytes。它们支持24GB显存、32GB主存存在合理余量，但GPU accounting不等于新作业的保证峰值，不能替代运行时检查。

读取了当前eligible分区/队列/维护预约：private/shared可见3090等合适型号，无维护reservation；私人队列已有其他研究作业等待。未调整或取消任何现有作业。

由于窗口条件已经失败，本次没有选定GPU节点或提交实验，也没有声称验证了137窗口在2小时内完成。历史完整R4耗时84:11:06包含三epoch与大量其他推理/诊断，不能直接当作本面板精确速度估计。实际新增GPU使用为0；CPU预检仅5秒。

## ChatGPT审阅

审阅已经返回，页面显示思考5分钟。会话：[审查研究进度](https://chatgpt.com/g/g-p-6a29d586630481918525796032225f68-ji-wangke-ti/c/6a8c18a3-4370-83eb-b9cf-8924409e2f12)。ChatGPT报告ACCESS_VERIFIED，读取上述精确提交，列出trace_decode_orientation、load_inference_model、load_checkpoint、predict_sequence及structural-head实现；明确本轮仅静态审阅，未执行代码或读取新Baobab outputs。本Agent也已直接读取相应Baobab源码，以下关键判断与源码相符。

### 消融是否可解释：可以，但范围严格局部

ChatGPT建议名称为frozen-lineage phase-gate bypass diagnostic。它测试既定chain上的门控阻断及过滤作用，不是phase-head architecture ablation、完整decoder消融或总体direct-annotator性能实验。另一个重要源码细节是：输出phase由所选CDS几何长度计算，head只负责否决；旁路不是把错误head类别写进GFF3。

采纳数值重放要求：原完整方向化6144tile、整条RC序列上的tile边界、末端A padding、eval mode、bf16 autocast、输出float32后转float16，再按原sigmoid/argmax和正反方向motif tie-break计算。不能改用未量化float32分数，不能在A不一致后仍解释B。本次未到GPU阶段，均属后续待验证合同。

### 抽样是否偏：是，前缀具有可避免的位置偏差

ChatGPT确认lineage serial来自按方向化位置递增的block枚举，指出每层取前N会偏向某端和局部环境。建议在每层完整排序列表上按 `floor((j+0.5)*N_eligible/N_requested)` 做确定性等距抽样，保留原组别/配额；同时记录CDS-count、已存phase-check-count和tile-edge状态。

该建议确实比前缀覆盖更多有序列表范围，采纳为**下一版设计修订建议**，但它仍不是概率随机抽样，不能据此提供总体CI或保证地理距离均匀。评阅提到“不增加预算”只能理解为样本数不增加：去重推理窗口可能增加，因此既不能保证≤64，也不能沿用本面板137的数字。

本次面板在审阅返回前已经按原提案冻结并被137>64阻断。遵守用户停止条件，没有再执行等距重采样或第二次核算。前缀面板保留为首版预检证据，不冒称已经满足评阅建议。

评阅另建议把token-conflict定义为任一gate-tested CDS start发生冲突。本次既有标签明确仅指原first-failure位置；两者不是同一变量，不能静默替换。后续若采用any-start标签，须单独命名并在新面板冻结前定义、计算；本次未改变历史标签。

### 决策阈值是否合理：可作启发式，不能作统计门槛

ChatGPT认为32/64、8/32缺少理论或数据校准依据，但可作为预先约定的探索性review flags。采纳用户要求：不用于总体recall/precision判定，不自动否决phase路线，无论计数如何都先交付逐项与分层结果。对刻意富集的确定性面板，不提供伪总体二项CI、bootstrap CI、p值或全基因组收益外推。

评阅明确支持“先冻结，再完整窗口核算，超64则GPU前停止、不替换”的顺序；当前业务终态符合此要求。

## 下一步（仅建议，待明确批准）

**推荐的唯一下一实验仍是同一固定epoch1 phase-gate旁路诊断，但先修订面板设计：**保留64/32/32组别及已定species/orientation配额，把层内ID前缀改成审阅建议的确定性等距选取，在新输出ID下重新冻结后仅做CPU窗口核算。保留first-failure token-conflict定义（不静默换成any-start），记录链长/CDS数和tile-edge供描述。不得复用或覆盖本次冻结文件；不按核算结果反复换样本。

这一步CPU预检建议上限2CPU、16GB、10分钟；实际去重窗口数未知，不能宣称等距面板只需137。只有新窗口预算得到明确批准且1GPU×2小时、8CPU、32GB资源条件可满足，才允许后续一次forward及配对分析。不要只把本轮上限改成137就继续旧前缀面板而忽略抽样审阅。

若64窗口上限必须完全固定，就不能承诺保留128条独立预选lineage；需要用户另行批准减少样本量或改成窗口簇条件下的诊断，且估计对象会改变。本次不替用户作这个科学/预算选择，也不继续运行预检或GPU任务。

32/64、8/32仅作描述性参考。即使未来低于/高于这些数，也应报告全部逐条ORF/threshold结果、分层差异、有限样本与选择不确定性，而不是直接否决整条phase路线。本次没有这些消融数值，不作猜测。

## 产物

Baobab项目根 `/home/users/j/jwang/coding-rna`（实际为 `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna`）下：

- `outputs/M25R-E1-PHASE-GATE-PANEL/frozen_panel.json`：128条身份、源trace、固定tuple、分层与选择规则。
- `outputs/M25R-E1-PHASE-GATE-PANEL/window_coverage.json`：逐lineage必要score位置、原窗口起点及窗口使用者。
- `outputs/M25R-E1-PHASE-GATE-PANEL/preflight.json`：137>64、业务停止状态及未推理标记。
- 同目录 `PREFLIGHT_JOBID`、`STATUS`、`logs/M25RPHASEPREF_12499507.out/.err`。
- 本报告 `reports/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/phase_gate_panel_report.md`。

没有覆盖原R4、勘误或grouping/phase回溯输出，没有读取Setaria，没有Git commit/push，没有新增训练，没有自动监控。
