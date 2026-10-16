# Carbon-A 完整结构流程：可执行资产与比较资格

截至2026-10-11 UTC。完整结构流程已定位到2026-10-09 17:10:01发布的固定revision `f6789ee79e6ef397de333577ffec1cb734496620`；35项资产已在Baobab准备完成，两份辅助网络严格CPU加载通过。此前14:46的“执行资产未齐”判断早于这次发布，不能继续用作当前结论。隔离环境的CPU导入证据及资格确认已通过，原人工GPU预检13500181因启动导入等待超时；同版软件归档及节点本地四原生API导入诊断均已完成并验证。保持原科学操作点的完整人工GPU重试13502103已于20:55:33 UTC启动，小型GPU内核与坐标夹具已通过，原生概率推理与GPU后处理已完成，但13502103在最终共同尺导入时因缺少SciPy失败。仅CPU的依赖修复与人工资格恢复13503202已完成并验证，保留旧失败及已完成GPU证据。原六TRAIN的原生预测均已完成；五条原作业成功，第二条的字面密码子校验修复由独立CPU13504510完成，原失败与预测保留不改。评价13503419现已COMPLETED0:0，固定过滤／未过滤pooled整链F1为78.172%／76.063%；详见[完成结果](carbon_a_original_train_result.md)。旧失败与缺失计时仍透明保留，不能转作泛化或纯decoder证据。

## 当前判断

Carbon-A现在具备可追溯的完整流程入口，应作为共享原始编码器权重的匹配对照准备，而不是只用逐碱基概率或数据库显示替代。作者[固定provenance](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/f6789ee79e6ef397de333577ffec1cb734496620/provenance.json)仍将主权重接到本项目原GENERanno revision b0483c23b6b63787b61a6d3a204a9b517d6ba345；这是来源关系，不是新的独立encoder。整个原生管线与本模型还存在上下文、评分、辅助网络监督及过滤差异，不能把最终差值归给解码器一项。

主模型与两份新辅助权重的监督来源须分别披露，见[分阶段训练暴露](../../../paper/supplement/pretrained_exposure.md)。当前固定provenance与辅助权重头没有补齐transition/confidence各自的训练／验证数据；不能把主模型的2055 assemblies清单转作辅助网络清单，也不能把未过滤输出称为无监督解码。

ANNEVO／Tiberius／Helixer的18组原六TRAIN预测及同尺评价已经完成，见[native_train_baseline_result.md](../../native_baselines/reports/native_train_baseline_result.md)；旧“缓存只覆盖两条DEV”的表述已过时。新增Carbon-A准备不替换这些基线、不回写历史八／十方法家族，也不增加已准备六方法DEV表的候选。以下旧revision与浏览器核对保留其对应快照范围。

## 新完整流程的固定资产与操作点

[固定发布目录](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/tree/f6789ee79e6ef397de333577ffec1cb734496620)与公开API一致；新transition_network.safetensors为49628字节，confidence_network.safetensors为63460字节。主model.safetensors为4614654504字节，公开HF LFS对象身份与原b048缓存对象相同。准备阶段复用原只读缓存文件，下载另外34项共约0.6MiB，避免再复制4.6GB；不修改、合并或覆盖原权重，不新增大文件checksum扫描。

[原生脚本](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/f6789ee79e6ef397de333577ffec1cb734496620/annotation_pipeline/run_annotation.sh)只执行概率推理和后处理，不调用参考评价或BUSCO。固定[管线说明](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/blob/f6789ee79e6ef397de333577ffec1cb734496620/annotation_pipeline/README.md)描述98304bp上下文、6144bp重叠平均、双链transition TCN、30状态Viterbi和confidence网络；原生输出为全部gene／mRNA／CDS的GFF及CDS／protein FASTA，record边缘不完整结构保留，confidence为NA。完整流程已明确不同于旧NPZ示例或数据库缓存。

原生默认使用BF16及FlashAttention-2；SDPA虽受支持，但作者明确说明概率不是bit-identical，不能静默替换。原生默认batch8、dispatch4、最短record100000bp、最短拼接CDS96bp、边缘内含子上限20000bp、标准密码表及chunk500Mbp；这些是作者操作点，不能套用M28的cap128或partial筛除规则。未来原六TRAIN均长于最短record门槛，但参考仍须全部保留，不按最短CDS重选分母。

作者GFF写出全部链，0.1置信度门槛只在下游gene级评价使用。未来匹配评价应分别保存作者0.1门槛与未过滤输出，NA边缘结构按作者规则保留，不根据本项目成绩挑门槛。共同尺仍须采用原22688条primary、坐标与strand敏感的CDS exact-key及CDS-base／背景口径，不能把作者按CDS序列一对一匹配的gene F1移进本项目坐标整链表。

[固定依赖](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/f6789ee79e6ef397de333577ffec1cb734496620/requirements.txt)要求Python至少3.10、torch>=2.7,<3、transformers>=4.56,<5，另需FlashAttention与C++17／OpenMP或CUDA的原生后端。当前generanno实际是Python3.10.18、torch2.5.1、transformers4.49.0，不能视为支持的完整原生运行环境；不得就地升级影响B／C／D／E训练。资产齐备、辅助权重CPU可加载、支持环境、真实FASTA全范围输出和指标有效性分别判断。

准备脚本prepare_carbon_a_assets.py与无GPU Slurm仅下载固定公开软件及两份小权重、复用原只读HF模型对象，核对35项文件大小、Python AST、Bash语法、两辅助网络严格CPU加载与有限FP32参数，并记录现有依赖缺口。此CPU加载不执行网络forward，也不证明torch2.5可运行完整管线；不读取研究DNA、参考、DEV或封存数据，不加载主模型，不更新参数、不安装依赖。输出为outputs/M28-CARBON-A-ASSET-PREPARATION/<job-ID>，与refs/carbon-a/<fixed-revision>独立保存。

## 独立资产准备的实际结果

命令sbatch --parsable sbatch/M28-CARBON-A-ASSET-PREPARATION.sbatch返回13500039。实际UTC2026-10-11 19:09:49至19:10:03，14秒；private-teodoro-gpu、gpu034、2CPU／4GiB／无GPU、5分钟上限、no-requeue。root／batch／extern均COMPLETED0:0，batch MaxRSS461804KiB；STATUS／JOBID／2225字节summary与真实handle一致。短准备由主直接核对消费，没有新心跳或新增监控对象。

资产根refs/carbon-a/f6789ee79e6ef397de333577ffec1cb734496620含固定35项文件，全部公开大小一致；34项下载合计605339字节。主权重为指向原只读HF blob的symlink，没有下载或读取主权重二进制、重算checksum或覆盖原文件。16份Python AST与3份Bash语法通过。原transition网络54张量／11077参数、confidence网络58张量／14433参数，均严格CPU加载且为有限FP32；没有执行任何网络forward或GLM推理。脚本用时7.909332599秒，不是注释速度。

实际依赖记录确认torch2.5.1和transformers4.49.0不满足固定发布的版本要求，supported_native_runtime_ready=false。该缺口是下一步独立环境准备的依据，不是原生模型准确率失败。环境、B／C／D／E训练、冻结DEV及封存范围均未变。

本次两次Exa定向检索各请求5条结果，10条返回含重复及其他Carbon语言模型页面，不等于10份独立全文。新发布资格依据官方固定目录／API、provenance、依赖及原生脚本与网络源码，并以实际Slurm资产收据闭合；只查询一般公开方法名称与URL，没有传输项目代码、DNA、权重或未发表指标。网页工具不能读取API JSON且Mac公开GET超时后，使用Baobab对同一公开API无凭据GET成功，未绕过私有权限。

## 独立支持环境准备（2026-10-11 UTC）

新环境路径refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1，仅复用generanno的Python3.10解释器建立不继承site-packages的venv，不安装或升级原环境。核心版本固定torch2.7.1+cu126、transformers4.56.2、flash-attn2.8.3；其余依赖按该固定Carbon-A requirements安装并保存实际pip-freeze。torch版本与CUDA12.6来源为[官方历史安装说明](https://pytorch.org/get-started/previous-versions/)，FlashAttention使用[官方v2.8.3发布包](https://github.com/Dao-AILab/flash-attention/releases/tag/v2.8.3)中已核实的cp310／torch2.7／cu12／linux_x86_64 wheel，TRUE／FALSE ABI分支由实际新torch决定；不根据推理成绩选环境。

命令sbatch --parsable sbatch/M28-CARBON-A-RUNTIME.sbatch返回13500092，提交后实际RUNNING于private-teodoro-gpu／gpu034，无GPU、4CPU／12GiB／30分钟、no-requeue。输出outputs/M28-CARBON-A-RUNTIME/<job-ID>保存STATUS、JOBID、pip-freeze及实际版本／导入路径summary。pip-check和导入检查若失败，则暂停native forward并仅修复具体安装问题；成功也只证明独立支持依赖能导入，不是FA2 GPU执行、主模型推理或原TRAIN评价。

该安装作业于19:48:35 UTC到达时限：root为TIMEOUT、batch为CANCELLED0:15、extern为COMPLETED0:0。日志显示Torch／Transformers／FlashAttention安装及pip-check成功，pip-freeze已保存，但没有最终summary；旧STATUS仍写RUNNING，不能据此视为存活。保留全部安装和失败记录，此终态是环境校验未完成，不是模型准确率失败。

一次限定修复只补跑同一环境的导入校验，不重新安装或升级包。sbatch --parsable sbatch/M28-CARBON-A-RUNTIME-VERIFY.sbatch返回13500426，实际19:50:11 UTC启动于gpu034，无GPU／4CPU／12GiB／15分钟、no-requeue。新目录outputs/M28-CARBON-A-RUNTIME/13500426保存独立STATUS／JOBID／pip-freeze／summary；pip-check、与原已安装清单直接比较、核心版本和全部导入路径检查须完成，并逐项打印导入阶段。准备的六TRAIN命令已改为引用该校验收据，未执行正式TRAIN。

19:56–19:57 UTC，在同一13500426 allocation中用两个只读srun步骤检查导入进程。实际PID2894945在import torch阶段处于D（disk sleep）、CPU0.0%，等待点folio_wait_bit_common；VmRSS由98192增至236644KiB，已映射新隔离环境的Torch／CUDA库，文件位于BeeGFS。证据支持正在加载库并有文件页等待，尚不能确认固定死锁或定位为集群文件系统故障。保留当前运行，不因静止日志重装、重提或改变FA2／模型操作点；若仍实际TIMEOUT，应先针对库文件读取与运行节点取证，而不是推定包版本不兼容或科学结果失败。

同一软件库的配对只读探针没有复现节点级读阻塞。13500710在gpu035以无GPU／1CPU／512MiB／2分钟上限执行，20:01:58–20:02:02实际4秒、root／batch／extern均COMPLETED0:0，batch MaxRSS12148KiB；随后同一13500426 allocation的step2在gpu034实际2秒、COMPLETED0:0。两者都只从libtorch_cpu.so与libtorch_cuda.so的头／中／尾读取各1MiB，共6MiB，完整读回且保存summary。不能把小样本读取当作完整导入速度，也不支持据此判定gpu034故障或全局BeeGFS故障；没有因此搬环境、重装、读主权重或研究数据。

13500426在20:05:09 UTC实际写出了3358字节完整summary，全部13包导入、核心版本及隔离路径检查通过，pip-check和与原安装清单的直接比较也通过；这一证据早于20:05:11时限。其Slurm仍在20:05:36判TIMEOUT，batch为CANCELLED0:15、extern为COMPLETED0:0，旧STATUS仍为RUNNING。运行中增加10分钟的scontrol请求被拒绝，没有更改其时限或将旧记录改为COMPLETED。完整导入证据与成功作业终态必须分开，不丢弃实际完成的前者，也不冒称后者成立。

资格确认采用sbatch --parsable sbatch/M28-CARBON-A-RUNTIME-QUALIFY.sbatch，返回13500755。实际20:09:41–20:10:16 UTC、35秒，private／gpu034／无GPU／1CPU／1GiB／5分钟上限、no-requeue；root／batch／extern均COMPLETED0:0，batch MaxRSS38480KiB。新收据outputs/M28-CARBON-A-RUNTIME/13500755重新完成pip-check、pip-freeze及与13500426清单的直接比较，再核实同一venv当前13包元数据与完整导入证据的版本／路径一致；没有再次导入ML、安装或升级。summary明确保存source_import_job_state=TIMEOUT、source_import_job_completed=false及恢复证据来源，原失败记录保留。此资格仅闭合已完成的CPU导入证据，不替代FA2 GPU核、主模型forward、原生解码或共同尺检查。

人工GPU预检13500181依赖afterok:13500092，提交后PENDING／有效Dependency；独立输出outputs/M28-CARBON-A-PREFLIGHT/<job-ID>，1RTX3090／4CPU／24GiB／15分钟。740017bp固定人工序列产生8个全98304bp窗口，同时检验作者batch8显存、重叠及非六倍数尾部；不因人工输入改变作者最短record100000、最短CDS96、dispatch4、BF16／FA2、标准密码表、20000bp边缘内含子和chunk500Mbp设置。原生auto的GPU路径显式固定为gpu，使用集群已存在CUDA12.6与GCCcore12.3，仅在新环境专属extensions目录编译官方CUDA Viterbi源码，失败不转SDPA或其他解码器。另用已知正负链跨内含子及边缘phase1／2 fixture验证作者assembly／GFF／CDS／protein导出，并用独立的坐标重建与翻译核验原生实际产物。fixture不是模型预测；人工raw成功也不是原TRAIN或准确率结果。预检不读取真实研究DNA、参考、DEV或封存数据，不更新参数。

原安装TIMEOUT后该依赖成为DependencyNeverSatisfied。主用scontrol update JobId=13500181 Dependency=afterok:13500426修复同一个尚未运行的预检；实际回读为新校验的有效未满足依赖，没有取消、重提或另建GPU预检，也没有改变其科学命令或输出目录。

完整导入作业TIMEOUT后，原预检仍未执行；主再次将同一个13500181的依赖修复为afterok:13500755。Slurm允许其PENDING时限由15增至30分钟，为实测库加载与后续原生执行留出时间；1RTX3090／4CPU／24GiB、全部科学命令与输出目录未变。CPU资格确认成功后，该原预检已实际RUNNING于gpu034。预检整体及其全范围产物尚未完成，不提前提交六TRAIN。

## 原六TRAIN的独立比较准备

run_carbon_a_train.py与evaluate_carbon_a_train.py已准备，但在支持环境、人工batch8原生forward／CUDA解码及真实共同尺fixture检查成功前不提交正式TRAIN。raw数组只读复用既有原4608步的六份input.fasta，仍是原TRAIN scope；固定新revision、官方模型代码与辅助权重、全部原生操作点，保留所有GFF调用和概率Parquet。原生完整、正负链全长度概率、独立GFF→DNA→CDS／protein重建与完成状态分别保存；原生密码表警告及native_WARNING目录保留，不救援或删调用。

全部六条raw完成后，CPU共同尺才读取原参考和既有控制结果。固定两行CarbonA_confidence_ge_0.1_or_NA与CarbonA_unfiltered，同一22688 primary分母、相同eligible isoform／unsupported分层／背景／CDS-union参照，与原4608步及已完成三原生基线分别比较。不以两行成绩选择门槛；NA边缘结构全部保留，不套用cap128或complete-only预测筛除。输出逐染色体、逐物种、pooled／macro／worst物种整链指标及基线差值，计时只作执行记录。共享backbone关系和多项管线差异不支持单独解码器因果结论；该新增TRAIN诊断不改变冻结DEV候选、方法家族或封存数据范围。

已提交raw数组0–5，当前并发上限6，每格1RTX3090／4CPU／24GiB。原一小时模板未经全流程验证；结合已完成的人工原生阶段和22:04 UTC调度信息，初始时限改为三小时，预算依据与不确定性见下文。评价仍为无GPU／2CPU／8GiB／20分钟afterok整数组。完成依据是每格真实Slurm成功终态、STATUS／JOBID／非空summary与完整实际原生产物，以及共同尺完整汇总，而不是提交或人工预检通过。

安装与人工预检的长等待已加入唯一原心跳coding-rna-r5-light-monitor，保存的完整B+C+D+E前缀16776字符逐字保留，再追加事件F13500092→13500181，总19173字符。实际回读核实ACTIVE／heartbeat、原target01a0f71f-8845-7432-bd13-bf54c1808fb3、15分钟、failed_runs_only与原created_at均未变；没有新对话、fork或重复自动化。原监控对话实际manual turn01a12c75-358f-7923-bcfc-205b258558ae于19:34:01–19:34:19 UTC完成一次新范围检查，19:34:10实际核实安装RUNNING／人工预检有效Dependency，未交接或暂停；真实turn_context确认gpt-5.6-luna／low未变。此为manual检查；随后原监控对话的scheduled turn01a12c7d-3894-7850-8cc3-cdc796f3277a于19:42:46 UTC实际追加了含事件F的原heartbeat输入，turn_context仍为gpt-5.6-luna／low，19:43:00实际检查运行及有效依赖后保持静默。新范围的配置、manual检查及scheduled追加／检查已分别验证；终态自动回传尚未发生。两项成功或实际失败才交接原主，主消费后继续原六TRAIN；不启动主重复等待循环。

上述心跳已在19:51 UTC将事件F更新为13500426→原13500181，原安装13500092关闭且由主消费、不再监控或重复交接。实际保存prompt为19356字符，原B+C+D+E的16776字符前缀逐字不变；heartbeat／原target／ACTIVE／15分钟／failed_runs_only／created_at均未变。原监控对话的scheduled turn01a12c8a-f475-7220-837e-bbe360258ced于19:57:46 UTC实际追加含新F的heartbeat输入，turn_context仍为gpt-5.6-luna／low；19:57:56实际检查新校验RUNNING／原预检有效Dependency，旧安装被排除，未交接或暂停。此次修复范围的配置与定时追加／检查已分别验证，终态自动回传仍未发生。

当前心跳的F已在20:11 UTC改为资格确认13500755→原预检13500181，总19874字符，原16776字符B+C+D+E前缀逐字不变；实际回读核实原heartbeat／target／ACTIVE／15分钟／通知偏好／created_at未变。旧13500092、13500426及已消费的软件库探针均关闭、不再监控或重复交接。原监控scheduled turn01a12c98-b03f-7981-be00-c2c670ef5b74于20:12:46 UTC实际追加含此新F的heartbeat输入，turn_context仍为gpt-5.6-luna／low，并于20:13:18完成定时检查。其实际squeue输出是原预检RUNNING gpu034，记录／静默结语却误写为PENDING；主读取该turn实际工具证据后修正原监控记录，20:17 UTC仅更新F的状态说明，要求逐行使用当轮真实JobID状态，保存19970字符并再次精确回读通过，全部B+C+D+E前缀与原配置不变。主20:16:58的Slurm账本确认预检实际20:10:17已启动、仍RUNNING，日志暂空不能判通过。随后原监控scheduled turn01a12ca6-6c20-7582-af46-cdfdb7f2a836于20:27:46 UTC实际追加修正后的heartbeat，gpt-5.6-luna／low未变；20:27:57的真实squeue仍为预检RUNNING gpu034，20:28:37完成定时检查并正确记为RUNNING，未再沿用旧PENDING，也未交接或暂停。状态修正后的定时追加／检查已验证，完整预检与终态自动回传仍未发生；没有新监控对话或重复自动化。

## 原生导入 I/O 诊断与同版软件暂存

原人工GPU预检13500181已在实际20:40:37 UTC达到时限，30分20秒；root TIMEOUT、batch CANCELLED0:15、extern COMPLETED0:0，两个只读诊断step均COMPLETED0:0。native.log仍只有inference:artificial，最终summary与概率产物缺失，故完整原生预检未完成。小型BF16／FA2内核及已写出的人工坐标夹具不能替代全上下文模型与后处理成功；旧STATUS可能仍RUNNING，未改写为成功。主已消费此终态并从活动监控关闭旧F，不重复交接；20:41 UTC当时尚未重提GPU或提交六TRAIN，后续同版暂存修复与新GPU作业见下文。

在原Slurm分配内的只读诊断显示，preflight父进程等待原生model_inference.py子进程2915830；子进程CPU约0.3%、RSS523408KiB、两线程，主线程SocketTk_poll，打开的是同一隔离venv内transformers/models/layoutlmv3/tokenization_layoutlmv3_fast.py，不是研究FASTA或主权重。实际Transformers4.56.2的create_import_structure_from_path递归列目录并读取各Python源码，缓存仅在当前进程内。这些证据定位到源码导入I/O等待；没有证明全局BeeGFS故障、固定死锁或模型准确率失败，也不能拿此前两节点6MiB抽样读取否定此元数据／小文件访问问题。

仅软件归档命令sbatch --parsable sbatch/M28-CARBON-A-TRANSFORMERS-ARCHIVE.sbatch返回13501094。实际20:36:16–20:40:35 UTC，259秒，private／gpu034／1CPU／512MiB／无GPU／15分钟／no-requeue，root／batch／extern均COMPLETED0:0；STATUS、JOBID、913字节summary及9696554字节transformers-4.56.2.tar.gz已闭合。它将当前隔离环境中同版transformers包与dist-info归档，省略__pycache__，不修改原依赖、不升级安装、不读取主权重或研究数据、不执行ML导入或forward。归档完成不是本地导入成功，也不是注释速度或模型性能结果。

依赖导入诊断命令sbatch --parsable --dependency=afterok:13501094 sbatch/M28-CARBON-A-LOCAL-TF-PROBE.sbatch 13501094返回13501119。资源private／1CPU／2GiB／无GPU／15分钟／no-requeue，实际20:40:36–20:53:04 UTC，12分28秒，root／batch／extern均COMPLETED0:0。20:54:47 UTC主核实STATUS、JOBID及非空summary，四个原生model_inference.py实际使用的AutoModelForTokenClassification、AutoTokenizer、PreTrainedModel、PreTrainedTokenizer均来自/tmp/m28-carbon-a-tf.13501119.MleEbm/transformers；Transformers4.56.2和Torch2.7.1+cu126未变，Torch仍来自原资格确认venv。Python导入阶段记录668.333秒，非整作业或注释速度。它不加载模型、不读DNA或注释、不执行GPUforward、不更新参数；实际API资格通过不替代完整GPU原生预检，也不允许缩短人工记录、降低batch或静默改SDPA来取得通过。

唯一原心跳于20:42 UTC关闭主已消费的旧F，仅保留活动B+C+D+E并追加G13501094→13501119；原16776字符前缀逐字不变，实际保存19226字符。原target、heartbeat、ACTIVE、gpt-5.6-luna／low、15分钟、通知偏好及created_at不变；没有新监控对话、fork或重复自动化。原监控scheduled turn01a12cb4-2802-7752-acda-933ad9bbbfbc于20:42:46 UTC实际追加新G心跳，20:42:56的真实SSH检查为B/C/D/E运行、后续有效Dependency、归档COMPLETED0:0、探针RUNNING；20:43:14完成并保持静默，未轮询或交接已关闭F。新G的定时追加／检查已验证，终态自动回传尚未发生。CPU归档单独完成安静，四API整链成功或真实需处理问题才交接原主并暂停同一心跳；主负责未来完整原生重试及六TRAIN，冻结DEV、原环境和封存范围未变。

完整人工GPU预检入口显式接收已完成归档与四API探针作业ID，在分配节点新目录解包同一软件并通过PYTHONPATH传给原生进程及spawn子进程；Python先核对导入来源与既有探针，再执行原来的740017bp／八全窗口／batch8／BF16／FA2／GPU Viterbi流程。summary分别记录同版软件来源与包含Python导入的耗时范围，Slurm Start/End用于整作业耗时。针对已定位源码I/O的物理文件暂存不改变模型、权重、原生源码、阈值、科学操作点或原环境。

六TRAIN的已准备sbatch也改为显式接收新的完整人工GPU预检ID与其同一软件归档ID，不再硬编码已超时的13500181；每项在自己的分配节点暂存同版软件，raw receipt记录实际来源，仍由完成预检约束真实FASTA访问。壳语法、Python语法及精确叶文件传输已验证，尚未执行任何六TRAIN推理；新的全流程GPU及全部六raw完成后，才进入已准备的固定评价。

截至20:51:35 UTC，13501119仍RUNNING，已越过Torch导入并进入四原生API阶段，尚无summary。20:52 UTC的实际集群配置与五项调度检查确认private两节点均为8张RTX3090／24GiB，gpu034使用5张、gpu035使用2张，无近期维护预约，当前并发未触及80作业上限。因此下一次同操作点GPU预检的已准备时限改为90分钟，仍1GPU／4CPU／24GiB／no-requeue／独立作业号目录。此为保守启动预算：流程涉及人工预检、原生主进程、spawn工作进程、后处理和公共评价夹具的多次独立Python启动，当前一次实际导入已超过11分钟；预留模型前向、扩展编译和输出I/O时间。90分钟不是实测注释耗时或速度比较结果，不改人工记录、batch、阈值、模型或原生后处理；导入探针完成前仍不提交GPU重试。

20:55 UTC再次读取实际配置、原作业闭合及五项队列／预约信息后，命令sbatch --parsable sbatch/M28-CARBON-A-PREFLIGHT.sbatch 13501094 13501119返回13502103。实际Submit20:55:32、Start20:55:33 UTC，RUNNING gpu034，private／1RTX3090／4CPU／24GiB／90分钟／no-requeue，无依赖和重启；日志outputs/M28-PILOT-R5/logs/M28CARBONAPRE_13502103.out/.err，输出outputs/M28-CARBON-A-PREFLIGHT/13502103。此运行须完成模型全上下文forward、原生transition／GPU Viterbi／confidence、坐标CDS蛋白回装及公共评价夹具，不能只凭已通过CPU导入或早期小内核判通过。当前完整执行尚未验证，六TRAIN raw/eval仍未提交。

同一原心跳于20:56:47 UTC关闭主已消费的G并追加H=M28-CARBON-A-NATIVE-GPU-13502103，实际保存19856字符，原16776字符B/C/D/E前缀逐字不变。回读确认ACTIVE／heartbeat／原target／gpt-5.6-luna／low／15分钟／failed_runs_only／created_at均未变，无新对话、fork或重复自动化；F/G由主消费而未发生watcher交接，分别记录，不伪造handoff。原监控scheduled turn01a12cc1-e3de-74a0-89bb-8ccbcde54dec于20:57:46 UTC实际追加含H及关闭F/G的新心跳，模型及推理强度未变；20:57:56的真实SSH检查为B/C/D/E/H全部RUNNING gpu034、其后续有效Dependency，未轮询已关闭源，也未读取科研产物、交接或暂停，20:58:03完成并保持静默。新H配置及定时追加／检查分别验证，终态自动回传尚未发生。H及原四链的首个整链完成或需处理异常只交接原主一次后暂停同一心跳，原主完成结果验证和后续工作。

21:22–21:25 UTC的只读阶段诊断补充了13502103的实际进展：21:22:02调度快照为RUNNING、26分29秒，原90分钟时限未变。stdout在21:06:40已写出coordinate fixture passed; begin full native annotation，因此本轮小型BF16／FA2内核及坐标夹具确实通过，而非只进入检查。21:22读取的native.log最后更新时间为21:13:49，已记录原生主进程读取人工记录、创建动态任务并等待推理；仍未写出Loading model、有效注意力实现或完整概率产物，不能判定全模型forward成功。

现有分配内的只读srun进程快照确认preflight2947649→model_inference2956304→spawn worker2962190属于同一batch cgroup。worker两次快照为Dl／folio_wait_bit_common，随后累计CPU从7秒增至19秒；这些是文件页等待快照，不证明整个11分钟持续阻塞或死锁。实际worker继承PYTHONUNBUFFERED=1、offline设置以及PYTHONPATH=/tmp/m28-carbon-a-tf.13502103.WOybVk；已映射Torch／CUDA／FA2动态库仍来自原资格确认venv。结合尚未出现的模型加载首标记，证据更支持启动阶段软件读取延迟，尚未定位某个具体库或证明全局存储故障。此诊断不读取真实FASTA、参考或模型准确率，不修改原生源码、环境、参数、时限或运行；继续让13502103完成，不因此重复提交，六TRAIN仍等待完整终态及产物资格。

21:45:02 UTC的阶段取证显示13502103仍RUNNING，49分29秒，90分钟时限未变。native.log已记录1153.7M参数模型完成加载（17.22秒）、effective=flash_attention_2／GenerannoFlashAttention2、人工记录推理100%、raw probability parquet写出及原生主进程成功退出，随后进入post-processing。日志中的Total execution time为11m17.23s，包含该进程及worker启动、推理、序列化和关闭，不是纯forward时间，也不按染色体bp比例直接外推。它证明实际全上下文模型运行已越过此前的导入阶段；完整概率／坐标回装资格仍以预检最终检查为准，不能用这份中间产物替代完整基因调用与公共尺成功。

21:46 UTC同一分配的只读进程取证定位了下一阶段：post_processing.py2973170已运行约19分16秒、CPU约0.2%，等待构建子进程，并持有原runtime/extensions/carbon_a_viterbi_gpu/lock。活跃nvcc2979768使用CUDA12.6、GCC12.3、sm_86、C++17与原Torch头文件；其gcc／cc1plus2984689／2984690子进程正在预处理viterbi_gpu.cu。cc1plus快照CPU约0.5%、SocketTk_poll，打开原隔离venv的ATen/ops/batch_norm_elemt_native.h。这是实际首次GPU扩展构建与头文件读取证据，未证明全部后处理耗时来自这一文件、全局存储故障或编译死锁。锁有活跃所属构建，不删除、不清缓存、不改CPU后端或重提；完整transition／GPU Viterbi／confidence及最终产物仍待验证。首次冷构建、软件启动与真正注释计算须分开解释，六TRAIN时限仍在全流程完成后依据实际运行与最新队列确定，不将准备模板中的一小时当作已验证预算。

22:00:24 UTC的读取确认原生概率推理及GPU后处理均已成功返回；prediction.gff／.fna／.faa在21:51:12 UTC写出。原生后处理记录transition网络34.9秒、GPU Viterbi0.4秒、文件处理96.6秒；源码计时从post_process_probability_file内读取概率开始，到三类注释文件写出结束，位于choose_viterbi_backend首次扩展编译及两辅助网络加载之后。因此96.6秒不含该进程启动、冷编译或网络加载，不能替代整作业耗时。summary记录3427.359秒，覆盖主预检Python导入、人工夹具和原生调用，另不包含归档解包及随后独立共同尺进程。

本人工记录产生零个基因／蛋白，GFF仅16字节头，CDS及protein FASTA为0字节；这是有效的空调用，不是模型准确率失败或缺失产物。原生confidence网络加载路径已完成，但没有完整基因可评分，其inference被跳过；日志median confidence nan来自空集合，不是双链概率出现非有限值。主预检已核实唯一Parquet记录与740017bp人工DNA逐字一致，双链全长概率均有限且在0到1之间；非空正负链、两外显子及边缘phase1／2的坐标CDS蛋白回装由四个已知夹具覆盖，不冒称来自模型预测。真实TRAIN中的非空调用、完整基因confidence评分及其回装仍须在各raw运行中验证。22:05:08 UTC，root／batch／extern仍RUNNING，common_ruler_checks.json尚未写出；现有summary不能替代独立公共解析、exact-key、过滤／NA保留检查及完整成功终态，H未关闭，六TRAIN尚未提交。

22:04 UTC的配置及五项实时调度检查确认private gpu034／035各配置8RTX3090、分别使用6／2张，没有切断本轮运行的维护预约；现有作业及新增六格数组加一次CPU评价未触及80并发／100数组限制，scratch可用146TiB。三小时raw初始预算保留原一GPU／4CPU／24GiB／并发3／no-requeue及所有原生科学参数：最长TRAIN为30427671bp，人工文件处理96.6秒提示大输入后处理需要另留时间，多次Python启动及模型forward也须单独计入；染色体调用量、confidence计算及I/O尚无真实全范围速度。此预算是执行判断而非按bp直接外推的实测成绩或公平速度比较，既不把整预检冷启动耗时按bp放大，也不以一小时模板作为成功保证。仅更新待提交raw的时限，不延长当前预检、不改变训练步数／模型／阈值／DEV家族或数据范围；仍待H完整资格通过后才提交原六TRAIN。

13502103实际20:55:33–22:06:05 UTC，70分32秒，root／batch FAILED1:0、extern COMPLETED0:0，STATUS为FAILED exit=1。22:06:47主读取实际账本，随后stderr定位到check_carbon_a_train_integration.py→评价模块→src/m28/query_model.py的scipy.optimize.linear_sum_assignment导入；隔离环境缺少SciPy，原generanno实际版本为1.15.3。失败发生在已有原生推理、GPU后处理、全长有限概率及坐标夹具检查之后，不是准确率失败；原summary与原失败状态全部保留，缺失的共同尺收据不补写进旧目录，也不把旧作业改为COMPLETED。旧H由主消费并从原心跳关闭，watcher未交接，不能伪造handoff=YES。

资格恢复命令sbatch --parsable sbatch/M28-CARBON-A-PREFLIGHT-QUALIFY.sbatch 13502103 13501094返回13503202，实际Submit22:14:27、Start22:14:28 UTC，22:14:52为RUNNING gpu034，private／2CPU／4GiB／无GPU／30分钟／no-requeue。新独立目录outputs/M28-CARBON-A-PREFLIGHT/13503202仅复制已完成的人工原生证据，记录source_native_job=13502103、source_native_job_state=FAILED、native_GPU_inference_reexecuted=false；不重跑已完成GPU阶段、不改旧记录、不读取真实TRAIN／参考／DEV／封存、不更新模型。脚本只在既有隔离venv用binary wheel／--no-deps补SciPy1.15.3，须通过pip-check及所有原已安装版本保留检查，原generanno不变；同版Transformers仍从13501094归档在新分配的/tmp目录解包。然后执行真实公共解析、坐标与strand敏感exact-key、0.1／NA与未过滤人工夹具检查，成功资格须同时有新job完整0:0终态、新summary、依赖修复与公共尺收据以及源GPU阶段的完整证据。提交或复制旧summary不算资格成功；新summary区分源原生耗时与本CPU恢复耗时，也不得称本I执行了GPU推理。恢复成功后，原六TRAIN直接使用新资格ID13503202与同一归档13501094，runtime13500755仍指原核心支持证据，新增评价依赖由I另行记录；无需改变现有raw科学入口。

唯一原心跳于22:15 UTC追加I=M28-CARBON-A-PREFLIGHT-QUALIFICATION-13503202，仅活动B/C/D/E/I，H/F/G均主消费并关闭。实际回读核实20503字符prompt、原16776字符B/C/D/E前缀逐字保留，以及ACTIVE／heartbeat／原target／15分钟／failed_runs_only／created_at未变；未新建对话或重复自动化，未更改原监控gpt-5.6-luna／low。22:12的原scheduled turn01a12d06-8f65-7b62-b067-52739a3b2910已完成关闭H后的B/C/D/E检查并保持安静；该turn早于I提交，不作为I定时运行证据。I的配置成功、实际CPU启动、完整资格成功及定时交接分开判断：当前只确认前两项，后两项仍待验证，六TRAIN仍未提交。

13503202完整资格已闭合：实际22:14:28–22:25:24 UTC，10分56秒，root／batch／extern均COMPLETED0:0，batch MaxRSS950572KiB。22:30:03主核实STATUS／JOBID、summary／依赖修复／公共尺收据、全部人工文件与源原生副本；主记录652.694秒是本CPU恢复Python阶段，不是GPU推理或注释速度。所有既有包版本保留、原generanno未改，SciPy1.15.3来自隔离venv。四条人工坐标链、两条NA边缘保留、低于0.1移除与未过滤解析／exact-key均按预期；这是工程与评价实现资格，不是模型准确率成绩。源GPU13502103仍FAILED1:0，native_GPU_inference_reexecuted=false明确保留，不伪造原生重跑或旧作业成功。

原监控scheduled turn01a12d14-4b20-7172-bfea-7899187c2856于22:27:56核实B/C/D/E运行、22:28:07核实I root COMPLETED0:0，却把独立单作业I误当中间阶段，没有完成产物核验、回传或暂停。主据实际输出取证后自行完成资格消费；I watcher handoff仍NOT_SENT，不补发或伪造YES。原心跳与journal已明确纠正事件边界、关闭I，未来J的整链完成独立于仍运行的B/C/D/E，避免沿用旧F/G的CPU→GPU中间阶段语义。

正式原六TRAIN命令sbatch --parsable sbatch/M28-CARBON-A-TRAIN-RAW.sbatch 13503202 13501094返回13503400；Submit及首三条Start均22:34:56 UTC，private／每格1RTX3090／4CPU／24GiB／3小时／no-requeue。评价命令sbatch --parsable --dependency=afterok:13503400 sbatch/M28-CARBON-A-TRAIN-EVAL.sbatch返回13503419，实际Submit22:36:13 UTC，依赖afterok:13503400_*整数组，仍无GPU／2CPU／8GiB／20分钟。提交前再次核实当前配置、五项调度、无输出碰撞与scratch146TiB可用；只读复用原六input.fasta，runtime13500755及恢复资格13503202／软件归档13501094，全部模型、原生参数、两种固定confidence操作点、22688参考、数据隔离及DEV家族不变。

22:47 UTC的五项实时调度及gpu035配置确认该节点8张同型24GiB／compute8.6 RTX3090仅占用2张，CPU及内存足够剩余三任务。因此scontrol update JobId=13503400 ArrayTaskThrottle=6在22:48 UTC成功，同一数组并发上限3→6；没有取消／重提或改每格资源、训练与推理参数，也未新增任务。后三条实际22:48:13 UTC在gpu035启动，22:50:35六条均RUNNING，CPU评价仍为有效Dependency。实际数字root映射来自sacct JobIDRaw／JobID，不能按数组ID加index推断：

| task | seqid | 实际数字root | 节点 | Start UTC |
| --- | --- | --- | --- | --- |
| 0 | NC_003070.9 | 13503401 | gpu034 | 22:34:56 |
| 1 | NC_003071.7 | 13503402 | gpu034 | 22:34:56 |
| 2 | NC_003076.8 | 13503403 | gpu034 | 22:34:56 |
| 3 | NC_089043.1 | 13503597 | gpu035 | 22:48:13 |
| 4 | NC_089044.1 | 13503598 | gpu035 | 22:48:13 |
| 5 | NC_089046.1 | 13503400 | gpu035 | 22:48:13 |

唯一原心跳现仅活动B/C/D/E/J，J=M28-CARBON-A-ORIGINAL-TRAIN-13503400-13503419。实际保存20996字符prompt、原16776字符B/C/D/E前缀逐字保留；ACTIVE／heartbeat／原target／gpt-5.6-luna／low／15分钟／failed_runs_only／created_at不变。J须六RAW加一次EVAL共21条真实root／batch／extern成功记录及全部结构产物，不能只凭数组、单项或资格完成宣布整链成功；新RAW无需另存input.fasta副本，native_WARNING与合法零预测FASTA允许。J配置已核实，但J范围的scheduled追加／检查及终态自动回传尚待验证。真实TRAIN预测、非空confidence路径、坐标回装及固定评价仍以各任务实际产物为准，当前不解释未完成准确率或声称竞争力／泛化。

23:12的原监控scheduled turn01a12d3d-7eac-7143-8822-dc175abfd893实际发现task1数字root13503402／batch FAILED1:0（23:09:09 UTC，34分13秒），extern COMPLETED0:0，并真实回传原主、暂停同一个心跳；首个终态自动交接及暂停已验证，不再仅是配置成功。主随后按实际数组日志M28CARBONATRAIN_13503400_1.err取证，不按数字root文件名查找：原生概率、GPU transition／Viterbi／confidence及4172个基因／蛋白均已完成，唯一错误是自己的carbon_a_output.read_outputs对gene004167的蛋白比较。其他五条数字roots13503401／13503403／13503597／13503598／13503400的root／batch／extern全部COMPLETED0:0，原六范围完整，没有取消或重提原作业。

只读CPU诊断13504239实际23:20:45–23:20:48 UTC、3秒、root／batch／extern均COMPLETED0:0。失败基因为NC_003071.7正链19673370–19679494，12个CDS、phase0、3966bp拼接CDS／1321aa，只有一个N。CDS offset0=2559／aa index0=853的GCN：作者固定标准表把非ACGT密码子译为X，Biopython的歧义翻译推断为A；字面表完整蛋白与作者输出逐字一致，而歧义翻译不一致。问题在独立校验口径，不是坐标／模型推理／蛋白写出或准确率失败。只修自己的标准表翻译为逐三元组查表、未知X，保留原DNA逐字、phase、坐标、内部stop、confidence与概率检查，不修改发布资产、原预测、阈值或数据范围。

新CPU恢复命令sbatch --parsable sbatch/M28-CARBON-A-TASK1-REVALIDATE.sbatch返回13504510，private／2CPU／8GiB／无GPU／30分钟／no-requeue，实际23:31:47–23:33:09 UTC、1分22秒，root／batch／extern全部COMPLETED0:0。完整原19698289bp DNA与Parquet逐字相同、正负链概率全长有限且0到1、全部4172基因／蛋白与GFF／CDS／phase回装通过；未重跑GPU、改原预测／失败STATUS／JOBID、读参考／DEV／封存或更新模型。接受收据独立保存outputs/M28-CARBON-A-TRAIN-RAW-VALIDATION/NC_003071.7，明确source_inference_job13503402/source_inference_job_stateFAILED/native_GPU_inference_reexecutedfalse；新验收成功不把原GPU作业改成COMPLETED。原wrapper时间因旧验证失败未写summary，seconds保留null；CPU77.835秒及源allocation2053秒分别标注，不冒充同口径native速度。

原PENDING评价13503419因旧整数组有失败变DependencyNeverSatisfied；scontrol update JobId=13503419 Dependency=afterok:13504510已成功修复同一个job，原2CPU／8GiB／20分钟／无GPU／科学命令及输出目录未变，实际23:33:10 UTC启动，23:33:31仍RUNNING。评价接受器只将已知task1切到独立恢复收据，其余五条沿原成功收据；报告指向实际接受路径而非缺失的原task1 summary。原22,688 primary分母、两confidence操作点、unsupported／背景／CDS-base口径与冻结DEV家族均未改；当前无真实准确率结论。

同一原心跳恢复ACTIVE，仅活动B/C/D/E/K，K=13504510→原13503419；旧J失败已真实交接YES且主消费关闭，不重复唤醒。实际回读20952字符prompt／原16776字符前缀逐字不变，原heartbeat／target／gpt-5.6-luna／low／15分钟／failed_runs_only／created_at未改，不新建对话或自动化。K的验证与评价两个job共6条真实成功记录及完整产物才整链完成，原13503402FAILED是透明来源不是新K失败；K配置已核实，新增范围的scheduled追加／检查和终态交接仍待实际验证。

## 早期官方入口快照

| 入口 | 已核实的身份与功能 | 对本项目比较的影响 |
| --- | --- | --- |
| [Carbon-A-1.2B 固定模型](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/tree/725713317756b1726ca026c281b31b9e59c4b25a) | revision `725713317756b1726ca026c281b31b9e59c4b25a`；公开模型 API 列出18个文件，包含模型权重、tokenizer、`infer_genbank.py`、`infer_packed.py`。固定模型卡说明输出双链逐碱基 CDS 概率，而非完整基因／转录本结构；发布示例输出 NPZ。 | 可作为已定位的概率模型入口，不能把阈值后的 CDS 连续块当作作者完整基因调用器。该固定仓库清单未列出独立 transition/confidence checkpoint 或完整 GFF 调用器。 |
| [数据库浏览器固定 README](https://huggingface.co/spaces/HuggingFaceBio/carbon-a-database-explorer/raw/c6b2b4ac351eaf9680e8141d89bf1f6e76633670/README.md) | Space revision `c6b2b4ac351eaf9680e8141d89bf1f6e76633670`；SQLite 查找已有 accession/source_key，按需读取 Parquet segment 的双链概率并下载原始行。此固定 README 描述的显示／下载数据是 CDS 概率及元数据，不是 gene-model 调用入口。 | 数据检索不是对本项目 FASTA 的重新推理。浏览器的阈值轨道没有转录本链身份，不能直接计算本项目 exact-chain 指标。此结论限定于该浏览器契约，不推断整个数据库都没有其他结构产物。 |
| [官方集合](https://huggingface.co/collections/HuggingFaceBio/carbon-annotation-database)、[注释 bucket](https://huggingface.co/buckets/HuggingFaceBio/genbank-annotations)、[技术报告 bucket](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report) | 本次可见集合入口连接上述模型、浏览器、注释 bucket 和报告；注释 bucket 当前列出 `annotations`、`cds`、`indexes` 等目录，报告 bucket 提供 `Carbon-A.pdf`。 | 发布数据和报告证明存在作者流程／产物，不自动提供该流程的可执行复现资产。本次未读取这些目录中的生物学记录，未鉴定某个结构缓存是否匹配本项目输入。 |

固定模型的文件清单来自[公开 revision API](https://huggingface.co/api/models/HuggingFaceBio/Carbon-A-1.2B/revision/725713317756b1726ca026c281b31b9e59c4b25a)，功能来自[固定模型卡](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/725713317756b1726ca026c281b31b9e59c4b25a/README.md)；Space 身份来自[公开 revision API](https://huggingface.co/api/spaces/HuggingFaceBio/carbon-a-database-explorer/revision/c6b2b4ac351eaf9680e8141d89bf1f6e76633670)。本轮定向搜索的两个查询各请求5条结果，重叠及非目标结果不算已审阅的10份独立来源；最终判断依据上述官方入口，不是穷尽全网后的不可用证明。

[固定推理指南](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/raw/725713317756b1726ca026c281b31b9e59c4b25a/INFERENCE.md)进一步明确：GenBank／FASTA 示例没有复现独立的生产注释流程，其输出是 coding occupancy 而非 GFF gene/transcript；packed 示例同样输出窗口 NPZ。这里的差别有作者直接说明，不仅是按文件名猜测。

## 数据库缓存的资格不能由 accession 一项决定

浏览器固定 README 明确：assembly 命中包含部分覆盖，索引是有日期的发布快照，分段保留各自坐标，查询最多返回200行且没有分页。按 assembly 搜到若干记录并不证明所有目标碱基或全部调用均已覆盖。

它还说明窗口平均和已存 float16 概率不能精确重建原评估 logits。浏览器二值显示把双链阈值结果合并，和作者逐链验证口径不同。不能以浏览器显示结果代替作者论文性能，更不能接入我们的整链比较表。

若后续定位到 GFF/CDS/protein 结构缓存，应另外确认其 assembly 版本、FASTA／坐标对应、全目标范围完成产物、调用与置信度过滤规则，以及是否存在输入依赖差异。缓存比较应明确标为作者预计算产物，不能记为本项目已重跑的 raw-FASTA 推理。当前没有这种资格已成立的结论。

## 历史快照缺口与完整运行资格

此前已读取的43页技术报告第24–25页描述 transition TCN、Viterbi 语法、confidence TCN 及 GFF3/CDS/protein 输出，方法与上游训练暴露证据见[pretrained_exposure.md](../../../paper/supplement/pretrained_exposure.md)。公开模型卡的范围较窄，不能据此否认作者报告的完整方法。

2026-10-09 14:46的旧快照尚需定位下列对象；17:10发布已补齐前三项的原生入口与公开权重，支持环境及本项目运行验证仍待完成：

- transition 与 confidence 网络的代码、权重及输入接口；
- 完整结构解码／GFF 导出入口和作者使用的语法、窗口拼接、坐标处理；
- 未过滤与置信度过滤两种调用的固定规则，不能在本项目 DEV 上选择；
- 可在原范围 FASTA 上执行的环境与命令，然后通过原评价尺独立评分。

发布模型示例与技术报告的窗口重叠已有差异，不能把示例默认参数视为报告整链运行配置。以上对象未齐全时，不提交名为“Carbon-A 完整基线”的作业，不拿自建简化解码器冒充作者方法。若采用自建适配器，应另命名并说明差异；不能无声替换冻结比较中的一个方法。

## 对当前推进的决定

原 TRAIN 的目标函数、采样、原生上下文先验与联合条件实验继续按已提交依赖执行，由原固定监控负责等待。该来源追溯不改变 TRAIN／DEV 隔离、不放行封存数据，也不把作者的 GenBank 广覆盖或自报指标写成本项目泛化结果。

完整Carbon-A原六TRAIN raw与两固定操作点评分已经完成；评价13503419及只读过滤诊断13504740 root／batch／extern均COMPLETED0:0。其分物种、共同尺、过滤得失与解释边界见[完成结果](carbon_a_original_train_result.md)。安装／导入TIMEOUT和原13503402校验失败仍是保留的工程来源，不改成成功；原生预测未覆盖、恢复未重跑GPU。原唯一心跳已关闭主直接消费的K，仅活动B/C/D/E；保存17431字符prompt、原16776字符前缀逐字不变及原heartbeat／target／模型／15分钟／通知／created_at配置已核实。K没有经watcher终态回传，不伪记YES；新post-K定时检查已实际验证：原监控23:49:47–23:50:11 UTC scheduled turn保持gpt-5.6-luna／low，23:49:58 SSH仅检查B/C/D/E正常RUNNING与有效依赖，排除已消费J/K，静默且没有回传；随后保存配置仍逐项一致。原冻结DEV计划保持候选、阈值与基线不变，封存与公开发布边界未改变，研究goal未完成。
