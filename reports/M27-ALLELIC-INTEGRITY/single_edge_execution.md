# M27 单边检验执行记录

2026-09-29。CPU准备作业13290376，命令sbatch --parsable sbatch/M27-EDGE-PREP.sbatch；private-teodoro-gpu，gpu035，2CPU、4GiB、5分钟上限、0GPU。COMPLETED 0:0，21秒，MaxRSS211,616KiB。6项测试全部通过，测试执行2.458秒；原生Numba合成路径、实际ANNEVO矩阵单项差异、坐标/链向、非法映射停止、分量求和覆盖均经过检查。测试未读取真实HDF5或执行神经前向。

ANNEVO环境无pytest可执行文件，测试使用标准库unittest，无新依赖安装。原生代码无tracked修改。协议与代码在真实缓存解码前公开。当前记录时，真实两臂解码尚未启动。

提交前M27既有数据36,819,615,045 bytes（34.29GiB）；准备日志另占少量空间，预留新增2GiB仍低于50GiB总预算。实时Slurm：gpu034不可新分配（inval），gpu035可执行CPU任务，准备作业确已在gpu035完成；现有TE作业不受干预。维护预约在其他CPU节点且日期晚于本次。真实作业拟8CPU/32GiB/20分钟、0GPU，不因CPU节点携带GPU而计为模型GPU实验。

## 真实执行

协议/代码公开提交dfed2c440b6af13d3730afeab5acac6384519758，GitHub main与本地export HEAD核对一致后提交。命令sbatch --parsable sbatch/M27-SINGLE-EDGE.sbatch，返回job13290433。8CPU、32GiB、20分钟上限、0GPU。两臂尚无结果时记录，不把提交当完成；后续终态见single_edge_result.md/json。

终态：13290433及batch/extern均COMPLETED 0:0，gpu035，2分34秒，MaxRSS6,689,164KiB（6.38GiB），AllocTRES cpu=8,mem=32G,node=1，无gres/gpu；主stderr 0 bytes。A/B两子命令均0，分别94.5485272/51.4751169秒；runner151.0842065秒。全chr同PTC缓存重放通过；单边释放后530链不变、仍跳过174ntCDS，新增边未被实际选择、发射差0，按协议关闭fixed pilot，不补跑。

最终outputs/M27-SINGLE-EDGE-R1目录6,066,821 bytes；原七目录36,819,615,045 bytes，因此M27最终总36,825,681,866 bytes（34.30GiB），0数据删除。原始JSON的stage_bytes/M27_bytes在最后结果和日志写出前计算，各比最终少24,900 bytes，不将该尾差误作数据丢失。GPU累计4.313333h不变，无新NN/训练/Setaria。
