# M27 原生确认执行记录

2026-09-27。状态：COMPLETED。4/4目标与全chr上下文复现；详细结果见 [confirmation_result.md](confirmation_result.md)。下列准备记录保留其时间顺序。

Pro第12轮完成固定43e366a结果审阅后，采用 [confirmation_contract.md](confirmation_contract.md)。仅ANNEVO/POTEH与Tiberius/CRYBA4各syn/PTC，共4次原生整chr22确认，记录缓存不作分量干预。

CPU准备：sbatch --parsable sbatch/M27-CONFIRM-PREP.sbatch，job13233925 COMPLETED 0:0，12秒，2CPU/4GiB/0GPU，MaxRSS81,396KiB。19项测试通过（4.79秒）；stderr为空。tests/test_m27_confirmation.py新增5项实际决策测试，另14项为缓存/native返回/分类回归。新阶段仅修改记录上限参数，不改神经返回或分类规则。

输出：outputs/M27-NATIVE-CONFIRM-R1。协议、实现和CPU测试将先公开，再提交sbatch/M27-NATIVE-CONFIRM.sbatch。选private-teodoro-gpu的1×RTX3090、8CPU、32GiB、45分钟；实时Slurm显示private有空余GPU、无相关维护预约，未触及并发上限。原环境沿用实测annevo（CPU准备generanno），不采用旧框架配置的历史默认环境。训练/Setaria为0；提交不是完成。

实际先推送0a8436d5f12b936a430623353a7facab4a8d1ea3，再提交job13233933。终态COMPLETED 0:0，gpu034，19:39=0.3275GPUh，MaxRSS14,817,172KiB，stderr为空。全部4例/6子命令完成；AN两HDF5与Tiberius两套完整缓存保留，M27总34.29GiB。累计4.313333GPUh，剩余3.686667GPUh。没有重试、分量干预或mutant缓存重放。
