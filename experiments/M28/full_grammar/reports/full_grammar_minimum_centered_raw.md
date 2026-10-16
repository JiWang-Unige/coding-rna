# 局部分歧惩罚 raw FASTA 接入

固定 minimum-centered 能量规则的独立 raw FASTA 有限接线资格已完成：13572020的root／batch／extern均COMPLETED0:0，两个49152bp原TRAIN视图及14bp短输入共六方向与独立fresh重编码重放一致。输入仅为FASTA和两支正式末步checkpoint，没有训练、参考驱动选择或冻结DEV变更；这不是全部六TRAIN染色体的新raw准确率评价。完整缓存势评价与fresh native接线结果分别记录。

## 模型和推理

来源 E fit13499817、M fit13510324，均末 head9216／native4608。两支分别恢复自身 backbone、adapter／dense-qv 和 head，全部梯度禁用。每条 contig／方向 fresh 编码，24576窗口、12288 stride、nearest-containing-center ownership 不变；逐 base／role 在 FP64 中扣除各自 BG 后取最小值。原 whole-contig skyline、unit role、null0、cap128、codon carry、GT／GC-AG、min-intron5及 literal DNA／phase 验证不变。角色数组只存在内存，推理不读取参考、旧预测、训练标签或旧势数组。

入口 `scripts/experiments/M28-METHOD-RESTART/annotate_full_grammar_minimum_centered.py` 的必要参数是 `--fasta`、`--e-checkpoint`、`--m-checkpoint`、`--output-dir`；新输出必须不存在。正式来源为 `outputs/M28-FULL-GRAMMAR-JOINT-LORA-FIT/step_009216.pt` 和 `outputs/M28-FULL-GRAMMAR-DENSE-QV-FIT/step_009216.pt`。新的输入仍受当前数据授权范围约束，重计算经 Slurm 执行。

## 固定有限接线资格

复用已有两个原 TRAIN 49152bp 视图：AT NC_003070.9 半开14376960–14426112；rice NC_089046.1 半开14450688–14499840。裁切前直接比较两支完整原 DNA，不根据参考或预测选择新位置。另复用14bp短序列 `ATGCGTNAGTTTAA`，输入为小写／R歧义形式，以核对原规范化和多记录入口。

一次实际调用器运行三记录／六方向并导出 GFF3。独立 fresh 重编码后以 `minimum(E.double()-E.double()[:,:1], M.double()-M.double()[:,:1])` 重放，核对全部链、完整 scope gain、8位有效数字导出 gain、CDS区间、ordered phase、ownership、同向 span 非重叠。已完成数学13567792与接口13569064复用，不重复资格；GFF 的零预测是合法结果，不调阈值。这个资格不重新编码全部六 TRAIN 染色体，也不是新的准确率或泛化评价。

## 资源和完成条件

独立输出 `outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-RAW-NATIVE`，不覆盖旧模型、势、GFF或冻结结果。沿用先前双支 native 116秒／实测4.54GiB显存的资源尺度：1GPU、2CPU、16GiB内存、10min、generanno、no-requeue，新增输出声明256MiB。2026-10-13提交前 private 两节点 GPU 已全部占用，而 shared gpu049 有两张配置 RTX4090 未分配，选择 shared／RTX4090 24GiB；不会因此修改两组正在运行的 LR 实验。

全部 root／batch／extern COMPLETED0:0，STATUS／JOBID／started／summary、输入和 prediction GFF／summary 齐全，六方向原生读出与独立计算一致后，才称为接入完成。若确定性接线失败，最多同范围一次最小修复并使用新输出保留失败；不取消或延长现有作业。长等待复用唯一原监控对话和心跳，不新建监控。原 TRAIN 完整性能、数值接线、独立泛化和论文目标完成是不同判断。

## 实际提交和接线结果

两份 Python AST 与 SBATCH 语法通过，四个新叶文件已精确同步并回读一致。实际命令 `sbatch --parsable sbatch/M28-FULL-GRAMMAR-MINIMUM-CENTERED-RAW-NATIVE.sbatch`，返回13572020，提交13:11:16UTC。controller确认 shared／RTX4090×1／2CPU16GiB10min／Requeue0及唯一日志路径；提交前expanded own31、存储约160TiB可用、目标输出不存在，无时限内维护冲突。13:12–13:14UTC处于正常Priority等待，controller预计开始时间不是执行保证，未开始前没有输出目录不属于产物失败。

13572020实际于2026-10-13 13:15:18–13:17:26UTC运行128秒，extern于13:17:27结束；三条Slurm记录均COMPLETED0:0，MaxRSS8921144K。STATUS／JOBID／started／summary、input.fasta和prediction GFF／summary七项均存在且非空。来源仍为E13499817、M13510324及数学13567792／接线13569064；原完整DNA裁切前相等、多记录FASTA规范化和六方向独立minimum-centered MAP／GFF／DNA／phase／ownership／同向span检查均通过。AT两方向和toy两方向为零预测，rice两方向各一条；零预测合法，不据此选阈值或判断准确率。

实际调用器的两模型加载32.091709秒、raw循环19.074633秒，独立资格程序共72.020099秒；峰值allocated4870736896bytes、reserved5186256896bytes。推理未使用reference／labels／feature cache／保存势／旧预测，optimizer更新0、weights及旧预测未改、未写checkpoint、未用DEV／sealed，且`full_original_six_TRAIN_end_to_end_reencoded=false`。这些时间只描述固定小视图接线，不能外推全染色体成本或速度优势。

主已直接核验消费该资格，MINRAW CLOSED/main-consumedYES/NOT_SENT_BY_WATCHER；它在新MINRAW scope实际定时检查前完成并被消费，不虚构watcher定时验证或补发晚交接。原MINROLE及后来完整LR四格均已关闭；同一heartbeat当前PAUSED，保留原watcher、gpt-5.6-luna／low、15min和通知策略。研究回传仍用裸threadId="01a0332e-e0a9-7aa2-90f9-159405701833"、hostId="local"并省略model/thinking；没有新对话、重复自动化、goal变更或科学成功claim。完整LR结果见[四格结果](full_grammar_qv_lr_controls_result.md)。
