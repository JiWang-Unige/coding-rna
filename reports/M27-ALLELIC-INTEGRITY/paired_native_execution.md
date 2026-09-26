# M27 原生配对干预执行记录

2026-09-26。状态：READY_FOR_GPU，尚无配对推理结果。

Pro第十轮完成固定2239028的五路径审阅后，采用 [paired_native_contract.md](paired_native_contract.md)。新实现增加native-only请求轨迹，不改原生前向返回或自适应控制器；分类器的严格端点规则和未知病例分母已在首个mutant前固定。

CPU准备作业13227836（sbatch --parsable sbatch/M27-PAIRED-PREP.sbatch）：
COMPLETED 0:0，16秒，2CPU/4GiB，MaxRSS369,536KiB，0GPU。17项测试通过（3.56秒）；原21登记位点保留，生成26个完整chr22单SNV输入，计划ANNEVO24次、Tiberius20次。输入构造阶段对真实全序列做差异位置断言并验证链向密码子；未加载模型。

输入目录与账：outputs/M27-PAIRED-NATIVE-R1，当前1,343,339,694 bytes。公开的 paired_input_manifest.json 和 paired_tests.xml 为原始小文件副本；FASTA及后续完整运行产物仅在Baobab。

GPU拟执行命令：sbatch --parsable sbatch/M27-PAIRED-NATIVE.sbatch。1×RTX3090，8CPU，32GiB，6小时硬上限；固定顺序、串行且无自动重试。下一步按实时队列提交并记录job，不以准备通过代替研究完成。

本阶段预期决策：44次工具×输入完成后，依据固定W_m和共同9位点的paired identification bounds决定是否存在可确认的具体PTC特异绕行案例；不得称为独立生物学验证或Nature Communications主张。
