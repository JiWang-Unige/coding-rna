# M27 单边检验执行记录

2026-09-29。CPU准备作业13290376，命令sbatch --parsable sbatch/M27-EDGE-PREP.sbatch；private-teodoro-gpu，gpu035，2CPU、4GiB、5分钟上限、0GPU。COMPLETED 0:0，21秒，MaxRSS211,616KiB。6项测试全部通过，测试执行2.458秒；原生Numba合成路径、实际ANNEVO矩阵单项差异、坐标/链向、非法映射停止、分量求和覆盖均经过检查。测试未读取真实HDF5或执行神经前向。

ANNEVO环境无pytest可执行文件，测试使用标准库unittest，无新依赖安装。原生代码无tracked修改。协议与代码在真实缓存解码前公开。当前记录时，真实两臂解码尚未启动。

提交前M27既有数据36,819,615,045 bytes（34.29GiB）；准备日志另占少量空间，预留新增2GiB仍低于50GiB总预算。实时Slurm：gpu034不可新分配（inval），gpu035可执行CPU任务，准备作业确已在gpu035完成；现有TE作业不受干预。维护预约在其他CPU节点且日期晚于本次。真实作业拟8CPU/32GiB/20分钟、0GPU，不因CPU节点携带GPU而计为模型GPU实验。
