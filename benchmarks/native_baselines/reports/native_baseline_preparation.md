# M28原生工具对标准备

2026-10-03，Asia/Bangkok。为论文目标准备当前固定版本的ANNEVO、Tiberius和Helixer原生比较，替代不了正在运行的full-TRAIN Q/C0实验，也不修改其推进条件。本阶段仅物化既有两条DEV原始DNA、隔离运行资产和版本收据；不是准确率结果。旧预测、共享checkout、环境及checkpoint保持不变，不提交Git或对外发布。

当前：六次原生预测13358414_0–5与共同评价13358432均COMPLETED/0:0，完整7728条primary参考未变；ANNEVO/Tiberius/Helixer合并完整链F1为75.7025%/71.8047%/62.8512%。[完整结果、评价边界及与A3的差距](native_baseline_result.md)已记录；它们是两条DEV contig上的当前版本比较，不是干净zero-shot或论文成功。原OOM、R2失败资产和日志保留。隔离Tiberius实际2.0.8源码在实际2.0.7依赖SIF中执行，全TRAIN Q/C0与A6按原协议继续运行。下文早期状态保留为历史快照，当前结果以新结果记录为准。

## 为什么需要新版本与独立资产

现有ANNEVO checkout为37bdd9aa62ddf24fa55941fb827061f7ed49ce53，当前官方main为42c920f70f52f96c7c06c665a639ec47d9304fff，标注v2.3.3。Magnoliopsida权重在两份Git树中具有同一官方blob21487d877430abf8e552c509560fd1218f19305e；这不证明旧执行源码等于当前版本。本轮在独立目录重新取得固定提交，不更新C0所用的旧ANNEVO源码。依据为[官方源码](https://github.com/xjtu-omics/ANNEVO/tree/42c920f70f52f96c7c06c665a639ec47d9304fff)与当前本机VCS对象，非根据旧性能挑版本。

现有Helixer镜像实际标签为helixer_v0.3.6_cuda_12.2.2-cudnn8_1，当前官方v0.3.7修正Helixer.py一键推理的min_coding_length默认值，并提供同硬件/软件的确定性选项。旧缓存不失效，但不能作为最新原生执行证据。[官方更新说明](https://github.com/usadellab/Helixer/tree/d17bb496b8842542590644b2437798711c0125d5)

Tiberius旧资产为2.0.5；目标官方源码声明2.0.8，Docker发布同名标签，但实测固定OCI内是2.0.7，早期“实际2.0.8镜像”假设已被否定。本轮纠正为固定2.0.8源码的隔离安装，不让旧host包选择2.0.5，不修改原SIF版本字符串。模型兼容版本1.1.6、实际工具包2.0.8和实际基础SIF 2.0.7分别记录。[官方配置与源码](https://github.com/Gaius-Augustus/Tiberius/tree/871d4979568ae09e56de7eca5d7489046ee833ea)

## 固定科学范围

输入仅拟南芥NC_003074.8的23459830bp与水稻NC_089041.1的29936421bp。分别从既有FASTA流式提取，只物化指定记录；转换大写与M28原生caller一致，保留原始IUPAC字符、坐标和长度，不增加软掩码、蛋白、RNA或其他外部证据。本阶段没有参考GFF/标签/特征缓存读取，无训练、GPU前向或模型选择，不访问test/Setaria。新FASTA的既有DEV身份不等于独立泛化资格。

后续比较使用完整原生输出和既有7728条primary共同尺，unsupported、长链与重叠参考不删；所有原生推理结束后才接入参考。工具权重按适用类群和官方默认选择，不以本DEV成绩选择权重或参数。新表与历史缓存表分开，不根据成绩覆盖较弱或较强结果；full-TRAIN Q相对新等曝光C0的固定gate保持。

## 固定工具和模型

| 工具 | 固定运行资产 | 预训练模型 |
| --- | --- | --- |
| ANNEVO | v2.3.3源码42c920f70f52f96c7c06c665a639ec47d9304fff；既有annevo环境，不改安装 | Magnoliopsida，同提交的ANNEVO_Magnoliopsida.pt |
| Tiberius | 固定871d497的2.0.8隔离包；基础依赖SIF实际2.0.7，OCI内容地址6c853ae0664d6739d4226c3c0bc08d39d89aefa88d1bf0d20941a44d7ba2b72c | 原显式871d497 angiosperms.yaml与已下载angiosperms.tar.gz不变 |
| Helixer | v0.3.7，源码d17bb496b8842542590644b2437798711c0125d5；linux/amd64镜像内容地址e2294eb2c282c35b919933daa0d6c145b635bfac8b717dbff88e88accbde4303 | land_plant_v0.3_a_0080.h5，固定官方列表priority80 |

OCI内容地址来自发布方registry，不新增自制输出摘要；它固定实际拉取资产，防止latest漂移。容器内实际包版本和源码身份在准备后核实，不能仅凭文件名推定。Helixer权重来自固定提交的[官方优先级列表](https://github.com/usadellab/Helixer/blob/d17bb496b8842542590644b2437798711c0125d5/resources/model_list.csv)；官方选择器按priority升序，旧缓存指定m_0100不能冒充当前默认a_0080。选择变化不是DEV调参。

Tiberius显式配置softmasking=false、clamsa=false、default_seq_len200070；作者训练元数据包含水稻Kitaake和Setaria viridis。这里只读模型来源说明，不打开这些新的基因组或Setaria评价数据，不声称水稻为该模型未见物种。已补齐[当前预训练暴露记录](../../../paper/supplement/pretrained_exposure.md)：ANNEVO实际植物权重的作者声明排除拟南芥/水稻，Helixer实际a_0080训练列表包含两者，GENERanno上游暴露仍未知；不能称整张表为干净zero-shot。

模型下载使用官方固定Git提交、Zenodo条目及Tiberius作者站点；保留源URL、最终URL、HTTP版本字段、实际大小和运行内副本。Tiberius模型URL自身不是内容寻址的不可变版本，明确保留这个来源限制，不以新自制哈希替代真实训练暴露审计。旧资产完全不覆盖。

## 准备执行和检查

新目录outputs/M28-NATIVE-BASELINE-PREP包含inputs、assets、receipts及summary。脚本prepare_native_baselines.py先完成四项合成检查涉及的实际失败：记录串混/大小写坐标改变、目标缺失、长度错误、选择错类群或错官方优先级。错误会阻止准备；无证据的额外异常不扩建防御层。

准备作业经Slurm执行，private-teodoro-gpu但不申请GPU，2CPU/8GiB/90分钟，环境annevo，日志M28NATIVEPREP_<job>.out/err。先物化约54Mbp输入，再取得ANNEVO独立checkout、两份模型与两个容器；读取实际包/源码/HelixerPost身份。本阶段容器metadata读取不导入TensorFlow/Torch模型或执行推理。

Docker下载体积为Tiberius11188347311、Helixer3564834785 bytes；SIF转换和临时层还需空间。当前BeeGFS可用175070049206272 bytes，无经济成本gate或人为2GiB上限。提交前fresh五项调度检查已完成，private资源/CPU可用、无相关维护预约；不取消或调整已有作业。

准备成功需要两个唯一FASTA长度正确、ANNEVO固定源码/权重齐全、两份模型资产齐全、两个容器实际包版本正确及完整收据/summary/STATUS。成功只说明资产已就绪，不说明GPU兼容、预测准确率、独立测试或论文目标达到。若下载/环境失败，保留具体错误和已取得资产，在不变版本/输入范围内处理，不把基础设施失败判为科学失败。

## 后续原生执行

下一阶段在上述资产成功后为三种工具各做两条DEV原生注释，固定官方类群推荐操作点，不重训或扫描。ANNEVO batch32、8CPU、普通植物模式；Tiberius batch4、200070默认序列窗、unmasked且不使用ClaMSA；Helixer batch16、64152 subsequence、当前原生HelixerPost默认和deterministic选项。GPU/CPU容量与操作时窗根据实际就绪和运行证据申请；当前准备作业本身不自动执行这六次推理或参考评价。

所有完成产物以Slurm终态和真实输出共同判断。长等待只复用既有固定监控对话与同一心跳，准备完整后由主对话核实并接续；不新建监控或改主对话模型。最终公平原生结果和独立验证仍待执行，论文目标保持未完成。

## 本轮执行记录

2026-10-03 08:40 UTC：`sbatch --parsable sbatch/M28-NATIVE-BASELINE-PREP.sbatch`返回13358367，已在gpu034节点RUNNING；该准备作业申请0GPU/2CPU/8GiB/90分钟，GPU节点名不代表申请GPU。四项合成测试全部通过（0.012秒）。独立目录的两份输入已经完成，seqid/长度分别为NC_003074.8/23459830和NC_089041.1/29936421，均为唯一记录；源小写计数4969681/10055400，仅记录大写化口径，不改变坐标。

ANNEVO固定源码42c920f70f52f96c7c06c665a639ec47d9304fff及Magnoliopsida权重已取得；Tiberius angiosperms模型29825066字节归档和Helixer官方priority80模型已取得，正在拉取Tiberius固定OCI镜像，Helixer镜像随后执行。当前STATUS为RUNNING，尚无完整summary或容器内实际版本核验，因此不能报告资产全部就绪、GPU兼容、原生预测或准确率。日志位于outputs/M28-PILOT-R5/logs/M28NATIVEPREP_13358367.out/err。

full-TRAIN Q拟合13358083、C0拟合13358084及依赖推理13358085/13358086、共同评价13358090保持原协议、资源和依赖关系，没有取消或修改。准备产物独立保存，不覆盖旧基线、训练checkpoint或冻结结果。

监控复用固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3及同一自动化coding-rna-r5-light-monitor，cohort扩展为M28-FULL-TRAIN-NATIVE-PREP-13358367。保存后的heartbeat/ACTIVE/15分钟/目标对话及完整六作业prompt已精确回读。模型沿用gpt-5.6-luna/low；准备完整或完整full链/需要处理的故障才一次交接并暂停，其他正常状态静默。已发出一次人工启动检查；本扩展cohort的定时执行和自动终态回传尚未观察，不将旧cohort的08:30定时成功冒充本轮验证。

人工启动检查已于2026-10-03 08:41:34 UTC完成，确认准备13358367与两fit RUNNING、两infer及共同eval有效Dependency；未发送终态交接。实际turn_context为01a100ec-b5ea-73e0-961d-7a3f8647a2d2，模型gpt-5.6-luna、effort low，未变更主对话设置。新增cohort的定时运行和自动终态回传仍待自然触发。

## 原生推理与统一评价执行

2026-10-03 08:57 UTC，准备资产仍在Slurm作业13358367中运行，当前没有新的原生性能结果。已部署run_native_baseline.py与evaluate_native_baselines.py，六个独立输出位于outputs/M28-NATIVE-BASELINES/{ANNEVO,Tiberius,Helixer}/{arabidopsis_thaliana,oryza_sativa}，不覆盖历史缓存。两种工具的容器仅挂载本轮输入/资产为只读和该元素的新输出，明确传入Slurm CUDA_VISIBLE_DEVICES；每元素先验证实际backend只看见一块GPU并保存版本/设备收据。GPU兼容性仍待真实元素执行，不由CPU接口检查推断。

原生操作点保持前述batch与序列窗。ANNEVO使用annotation.py一键接口、Magnoliopsida、batch32、8个数据worker/8CPU解码，不启用overlap_pred或新阈值。Tiberius直接执行容器内/opt/Tiberius/tiberius.py，显式配置和预先下载的只读模型目录，batch4/200070/无软掩码与ClaMSA，不通过旧host包选择镜像。Helixer显式a_0080模型、64152/batch16、deterministic、原生overlap/HelixerPost默认，添加官方--no-multiprocess只改变FASTA数值化并行方式，避免其按整节点cpu_count启动38个worker超出本作业8CPU；不改变DNA、模型或解码阈值。

固定版本Helixer.py的后处理失败分支只打印错误而不传播非零退出码，故入口除subprocess退出外还要求其原生完成日志。三工具均保留原始GFF，不筛选/救回零预测，并检查指定contig、1-based坐标范围、CDS strand/phase/Parent。确实完成且零预测属于模型输出，缺少原生完成标记则属于执行失败。[固定源码](https://github.com/usadellab/Helixer/blob/d17bb496b8842542590644b2437798711c0125d5/Helixer.py)

统一评价只有六个原生终态及序列-only收据齐全后才打开既有参考，直接复用common_ruler的primary选择/完整链与CDS-base union/背景定义。拟南芥5437、水稻2291条分母及unsupported/长链/重叠保留不改；记录逐物种、micro、macro、phase诊断和背景覆盖，不将非精确预测数量误称漏掉参考。预训练暴露仍按方法分开，水稻不能统一叫未见物种；新旧差值可能同时包含版本、权重、输入集合和环境差异，不作版本单因素因果解释。运行seconds包含入口/GPU探测开销，仅为操作收据，不是公平速度benchmark。

入口检查13358405已COMPLETED/0:0，六项合成测试通过（0.005秒），覆盖六槽位、官方操作点、CVD传递、真实零预测、范围/phase错误和Helixer后处理失败。评价检查13358427已COMPLETED/0:0，三项合成测试通过（0.143秒），调用实际共同parser与ruler核对不等分母micro/macro/背景/零预测、未完成推理阻止参考读取、参考分母漂移中止。它们不是生物学准确率证据。旧parser出现已知unclosed-file ResourceWarning，测试未失败，未扩大修改共享评价器。评价与其检查使用现有generanno环境，匹配既有共同评价的依赖；原生入口使用annevo，不改base或既有安装。

| 阶段 | 实际作业 | 成功依赖 | 资源与新产物 |
| --- | --- | --- | --- |
| 六次原生推理 | 13358414_0–5 | 13358367与13358405 | 每元素1 RTX3090/8CPU/32GiB/2小时，最多并发6；每元素原始GFF/日志/GPU收据/summary/STATUS |
| 原生共同评价 | 13358432 | 数组13358414全六元素与13358427 | 0GPU/2CPU/8GiB/20分钟；outputs/M28-NATIVE-BASELINE-EVAL |

实际命令为`sbatch --parsable --dependency=afterok:13358367:13358405 sbatch/M28-NATIVE-BASELINES.sbatch`及`sbatch --parsable --dependency=afterok:13358414:13358427 sbatch/M28-NATIVE-BASELINE-EVAL.sbatch`。即时Slurm为数组与评价PENDING/有效Dependency；尚未启动GPU原生前向。操作时窗不是运行速度主张：旧Tiberius整拟南芥记录约631秒可作规划参考，本轮24–30Mbp单contig给2小时上限与余量，不能保证新版本同耗时。新指标结果与Q/C0 full-TRAIN固定gate分开，后者原五作业/参数/停止条件完全不改。

提交前五项live调度检查、输出不存在及存储检查完成；private gpu034/035各配置8张RTX3090、GRES_USED5/1，CPU与内存可容纳本轮，暂无GPU维护预约。BeeGFS可用174905977470976字节。保留其他用户/项目作业与本地无关改动，无Git提交或公开发布。

同一监控对话与heartbeat已更新为M28-FULL-TRAIN-NATIVE-13358414，保留gpt-5.6-luna/low与15分钟，完整prompt/目标/类型/状态已精确回读。准备完成不再独立交接，只有全部六元素与原生评价完整，或原FULL整链完整/需处理失败，才一次通知主对话并暂停，明确另一链尚未完成。数组用JobID逐一核实_0–5，不以父数组/一个元素推定全完成。当前新cohort的人工启动检查已完成（08:59:50 UTC现场取证）；同一监控对话于09:00:37.622 UTC实际追加定时heartbeat turn01a100fe-7148-7010-879c-2f84c0f17cd5，09:00:51 UTC完成fresh squeue/sacct取证，检查成功，所有当前作业为正常RUNNING/有效Dependency，未交接。该定时turn_context再次核实gpt-5.6-luna/low。已区分配置、人工和定时成功；整链自动终态回传尚未发生，不能声称已完成回传或原生结果。

## 镜像压缩内存不足与技术重试

原准备13358367于09:02 UTC在Tiberius OCI层解压完成后、创建SIF的mksquashfs阶段失败，Slurm为OUT_OF_MEMORY/0:125、耗时23分46秒；batch MaxRSS7962084K，对应原8GiB资源申请。stderr同时记录mksquashfs exit137与OOM事件。rootless xattr警告不是这次失败原因。这属于执行基础设施失败，不是原生模型准确率失败；六元素当时全部未启动，尚无推理结果。

原outputs/M28-NATIVE-BASELINE-PREP、FAILED STATUS、JOBID、下载资产及M28NATIVEPREP_13358367.out/err全部保留，没有覆盖或清理。修复只改变镜像转换资源和新输出目录：Apptainer1.5.4-1的build支持--mksquashfs-args，本机随附mksquashfs支持-processors与-mem，实际使用`apptainer build --disable-cache --mksquashfs-args '-processors 2 -mem 2G' <new.sif> <same pinned OCI URI>`。2G是压缩缓存用量，不是整个进程的硬上限；Slurm分配提高至32GiB给其他转换内存留余量，仍2CPU/0GPU/90分钟。保留默认压缩方式，不改OCI内容地址、工具版本、权重选择、DNA或推理参数。[Apptainer官方说明](https://apptainer.org/docs/user/latest/cli/apptainer_build.html)

实际提交`sbatch --parsable sbatch/M28-NATIVE-BASELINE-PREP-R2.sbatch`得到13358480，输出为outputs/M28-NATIVE-BASELINE-PREP-R2，09:14 UTC为RUNNING gpu034。准备入口4项检查（0.010秒）与原生入口6项检查（0.005秒）已通过。评价入口因新的准备路径同步更新，`sbatch --parsable sbatch/M28-NATIVE-EVALUATION-CHECKS-R2.sbatch`得到13358482，已COMPLETED/0:0、STATUS COMPLETED，实际共同parser与ruler的3项合成测试通过（0.044秒）。13项检查仍只证明接口与合成口径，不证明真实GPU兼容或准确率。

原数组13358414_0–5和评价13358432不取消、不重建、不改资源。第一次试图重附已完成入口检查13358405被控制器拒绝：该ID已超出实际MinJobAge300秒、scontrol查询无此作业，原依赖未改变。修正为` scontrol update JobId=13358414 Dependency=afterok:13358480`，新准备成功本身要求上述4+6项检查通过；不是绕过检查。评价使用` scontrol update JobId=13358432 Dependency=afterok:13358414:13358482`，已完成检查的依赖由Slurm正常移除。实际回读数组全部PENDING/有效Dependency，评价等待同一六元素数组；SubmitLine仍显示原提交命令，不代表当前依赖。[Slurm依赖更新与MinJobAge说明](https://slurm.schedmd.com/scontrol.html)

同一heartbeat、固定监控对话、gpt-5.6-luna/low及15分钟保持，当前cohort为M28-FULL-TRAIN-NATIVE-13358414-PREP-R2-13358480。完整prompt/类型/目标/状态已精确回读，旧OOM由主对话消费，不重复回传。人工startup turn01a1010a-9966-7373-bc99-57dc57e696b6于09:14:08.419 UTC完成现场检查：新准备与两fit RUNNING，其余为有效Dependency，未交接；实际turn_context核实同一模型和low。PREP-R2 cohort的自然定时执行与完整终态回传当时尚未观察，不能用原09:00检查替代。固定同对话调度方式依据[官方说明](https://learn.chatgpt.com/docs/automations?surface=app)。

随后已自然观察到PREP-R2 cohort的同对话定时执行：turn01a1010c-2d30-7553-8bc8-7806ae3b0152于09:15:37.694 UTC启动，09:15:51.112 UTC完成现场检查，仍为新准备与两fit RUNNING、其余有效Dependency，未交接；turn_context为gpt-5.6-luna/low。配置、人工检查和本cohort定时执行均已验证，完整自动终态回传仍未发生。

原full-TRAIN五作业完全不改，仍等待各76536步末权重、原始FASTA推理与共同评价。当前下一步为验证修复后的真实容器转换、全部六次原生执行及同尺评价；尚无新的性能结论或论文成功证据。

## R2终态及实际Tiberius身份错配

准备13358480已FAILED/1:0，root/batch均失败，耗时59分14秒，batch MaxRSS2002312K，STATUS FAILED exit=1、JOBID一致。stderr分别记录两个SIF Build complete，随后prepare_native_baselines.py在实际包版本检查处失败；不是再次OOM、推理失败或科学指标FAIL。工作目录挂载警告与rootless xattr警告未阻止metadata读取或构建，不能用它们替代真实版本错误原因。

为核实错误具体对象，另做CPU-only只读身份诊断，不下载/重建、不导入GPU模型或读取参考。首次诊断13358882 FAILED/2:0、0秒，因脚本放在登录节点/tmp而计算节点不可见，尚未执行metadata；保留M28NATID日志。一次有界接线修正仅将相同脚本置于共享outputs/M28-NATIVE-IDENTITY-DIAGNOSTIC，SBATCH/M28NATIDR2入口同为annevo/private-teodoro-gpu/0GPU/1CPU/2GiB/2分钟。作业13358894 root/batch/extern均COMPLETED/0:0、6秒，MaxRSS81392K，stdout两条实际metadata记录齐全。

| 对象 | 预定包版本 | 容器实际版本 | 容器实际源码 |
|---|---|---|---|
| Tiberius | 2.0.8 | 2.0.7 | a1472d26829769c22c9678c9f042830bbb75e544 |
| Helixer | 0.3.7 | 0.3.7 | d17bb496b8842542590644b2437798711c0125d5 |

Tiberius metadata的Python为3.12.3，Helixer为3.10.12。SIF标签的base.digest与预定OCI分别完全相同；所以本次不是拉取了另一个内容地址，而是将固定发布资产误当作2.0.8。文件名Tiberius-2.0.8.sif及发布标签不能覆盖实际包/源码证据，原镜像不能以2.0.8名义进入主比较。Helixer身份匹配不等于整套准备或GPU推理已经成功。

已定向核实两份官方固定源码：容器实际a1472d2的[pyproject.toml](https://raw.githubusercontent.com/Gaius-Augustus/Tiberius/a1472d26829769c22c9678c9f042830bbb75e544/pyproject.toml)明确声明2.0.7，原目标871d497的[pyproject.toml](https://raw.githubusercontent.com/Gaius-Augustus/Tiberius/871d4979568ae09e56de7eca5d7489046ee833ea/pyproject.toml)明确声明2.0.8，且两者都要求Python≥3.12。不是仅修改显示字符串即可消除的错配。后续技术方案应在隔离运行路径部署原目标固定2.0.8源码/包，分别记录原基础OCI与实际运行包身份，验证CLI与GPU接线后才放行；不得仅放宽metadata断言，也不得把新的源码运行方式冒称发布方原始2.0.8镜像。本轮尚未安装新包、建立修复运行资产或再次提交准备/原生推理。

当前数组13358414_0–5为DependencyNeverSatisfied，原生评价13358432等待；未修改成功条件、版本断言、参数或参考分母，未取消/放行/重提这些作业。原生GPU兼容、六次预测和公平性能仍未完成。下一步核实官方2.0.8源码与运行资产，保持预定工具版本、模型选择和序列-only输入，不默认降级到2.0.7；若必须实质改变固定比较范围则先明确记录与授权边界。

固定轻量监控已于10:15:49.549 UTC在原对话01a0f71f-8845-7432-bd13-bf54c1808fb3成功自动回传本次新故障，主对话实际收到并消费，验证了真实自动回传而非仅配置。随后复用同一heartbeat coding-rna-r5-light-monitor，将当前cohort缩为M28-FULL-TRAIN-CONTINUE-20261003，仅等待原五FULL作业；不会重复回传已知native故障。ENTIRE prompt/ID/heartbeat/ACTIVE/15分钟/固定目标精确读回，10:20:27.690 UTC人工启动检查确认两fit RUNNING、两infer/共同eval有效Dependency；实际新turn_context仍gpt-5.6-luna/low。新FULL-only cohort的自然定时运行尚未观察，不能将原cohort回传当作新定时验证。沿用固定对话方式见[官方说明](https://learn.chatgpt.com/docs/automations?surface=app)。论文目标仍未达成。

## 固定2.0.8隔离运行修复

新增prepare_native_runtime208.py及M28-NATIVE-RUNTIME208.sbatch，作业13358944经private/0GPU/2CPU/16GiB/30分钟、annevo执行，root/batch/extern均COMPLETED/0:0，实际1分02秒，batch MaxRSS1411952K；新目录outputs/M28-NATIVE-RUNTIME208的STATUS、JOBID、summary与安装/模型加载收据齐全。初次部署补丁路径前缀错误导致sbatch找不到入口，尚无提交；仅将本轮两个新文件移回正确共享目录后取得真实ID，不覆盖他人文件或原资产。

只下载固定871d497源码，在原SIF的Python3.12.3中`pip install --no-deps --no-cache-dir --no-compile --target /m28-runtime/site <fixed checkout>`。临时build isolation只安装打包后端；不升级原TensorFlow/CUDA/bricks2marble/hidten，不改host/base环境或原只读SIF。安装后25份Python源码与固定checkout直接逐文件相等，实际distribution版本及package __version__均2.0.8，package/main导入路径均/m28-runtime/site/tiberius，不从/opt/Tiberius旧源码导入。固定原生入口为`python3 -m tiberius.main`，其参数解析器与预测主函数来自[该固定官方代码](https://github.com/Gaius-Augustus/Tiberius/blob/871d4979568ae09e56de7eca5d7489046ee833ea/tiberius/main.py)。

实测依赖distribution TensorFlow2.17.0+nv25.2（tf.__version__2.17.0）、CUDA12.8/cuDNN9、bricks2marble0.1.1、hidten1.0.0、NumPy1.26.4。原angiosperms权重只读挂载到隔离包同级model_weights，官方resolve_weight_download实际复用该目录；PredictionGTF.load_model成功，--help成功，没有读取DNA、执行forward或optimizer。重复CUDA factory注册日志未阻止CPU模型加载，不据此判GPU成功或失败。实际GPU兼容性仍须每元素probe及真实前向验证。

ANNEVO42c920f/v2.3.3和原Magnoliopsida权重、Helixer0.3.7/d17bb49和HelixerPost4d4799b、官方priority80 a_0080权重（26098184 bytes）、两条原始FASTA均复用R2资产，不重新构建SIF或重选模型。原R2在失败前未持久化Tiberius/Helixer最终HTTP收据，新的收据将此明确标为未知，不伪造ETag/Last-Modified或以新下载替代原资产。旧FAILED状态和实际2.0.7身份保留；修复不是发布方原始2.0.8镜像，不能如此署名。

run_native_baseline.py改为从新隔离summary核实就绪、从原R2只读使用DNA/权重/容器，Tiberius以新site PYTHONPATH及/m28-run工作目录调用模块，GPU probe核实2.0.8真实导入路径。原batch4/seq200070/unmasked/no-ClaMSA和其他两工具参数完全不变；evaluate_native_baselines.py仅跟随原R2 FASTA位置，完整7728共同尺不变。入口/评价合成检查13358959 root/batch/extern均COMPLETED/0:0、5秒，六项runner通过0.007秒、三项真实共同parser合成检查通过0.045秒；旧parser ResourceWarning未导致失败，未扩大修改共享代码。

主对话已核实上述准备终态及实际收据后，10:40 UTC以` scontrol update JobId=13358414 Dependency=afterok:13358959`恢复同一数组，没有取消/重建、改变资源或科学参数。控制器回读Dependency=null，随后实际_0/_1在gpu034、_2–5在gpu035全部RUNNING；原eval13358432仍等待同一六元素afterok依赖。SubmitLine仍显示早期命令，仅为历史，不是当前依赖。两Tiberius与两Helixer的真实GPU probe均成功、各只看见一张RTX3090；Tiberius使用原CUDA12.8/cuDNN9的TF2.17.0，Helixer原TF2.15.1/CUDA12.2/cuDNN8。探测成功不等于完整前向/预测成功，尚无完整GFF或共同准确率结果。

同一固定对话01a0f71f-8845-7432-bd13-bf54c1808fb3及heartbeat coding-rna-r5-light-monitor扩为M28-FULL-NATIVE-A6-20261003，仍gpt-5.6-luna/low、15分钟，等待原FULL五作业、恢复native六元素/评价、A6新fit13358977→infer13358978→eval13358979。正常运行/有效依赖静默，任一整链完成或新故障一次交接并暂停；不重复旧R2已消费失败。整份prompt/ID/heartbeat/ACTIVE/间隔/固定target精确回读通过；10:42:35.787 UTC人工startup实际追加于同一对话，turn01a1015b-949b-7e03-8fa0-77fd3e12c5fe的model/effort实测gpt-5.6-luna/low。该检查确认两FULL fit与A6 fit RUNNING、ANNEVO两个元素COMPLETED/0:0、其余native RUNNING，下游均有效Dependency；没有交接。新cohort自然定时及整链自动终态回传尚未观察，不用历史证据替代。模型版本、输入和评价边界按validate-data核对，write-page仅更新原项目记录，无托管/公开发布；科学目标仍未达成。
