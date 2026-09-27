# M27 原生确认执行记录

2026-09-27。状态：CPU_PREP_COMPLETED / GPU_NOT_YET_SUBMITTED。

Pro第12轮完成固定43e366a结果审阅后，采用 [confirmation_contract.md](confirmation_contract.md)。仅ANNEVO/POTEH与Tiberius/CRYBA4各syn/PTC，共4次原生整chr22确认，记录缓存不作分量干预。

CPU准备：sbatch --parsable sbatch/M27-CONFIRM-PREP.sbatch，job13233925 COMPLETED 0:0，12秒，2CPU/4GiB/0GPU，MaxRSS81,396KiB。19项测试通过（4.79秒）；stderr为空。tests/test_m27_confirmation.py新增5项实际决策测试，另14项为缓存/native返回/分类回归。新阶段仅修改记录上限参数，不改神经返回或分类规则。

输出：outputs/M27-NATIVE-CONFIRM-R1。协议、实现和CPU测试将先公开，再提交sbatch/M27-NATIVE-CONFIRM.sbatch。选private-teodoro-gpu的1×RTX3090、8CPU、32GiB、45分钟；实时Slurm显示private有空余GPU、无相关维护预约，未触及并发上限。原环境沿用实测annevo（CPU准备generanno），不采用旧框架配置的历史默认环境。训练/Setaria为0；提交不是完成。
