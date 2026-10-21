# 16 PlantCAD2 作为第二通用骨干的可行性评估（2026-10-19，只读调研 + ont 短测，未做任何训练）

产物说明：本文件是唯一的项目内新增产物。ont 上新增了 `~/envs/plantcad2`（6.3 GB，未改动任何已有 env）和 `~/plantcad2_work/`（模型 8.3 GB + 测试脚本/日志）。没有动 m30/M31 或他人的进程，没有提交 Slurm，所有 GPU 测试单次 ≤15 分钟，CPU 绑在核 4–9（避开 0/1）。

来源标记：**[实测]** 在 ont 上跑出的；**[HF]** Hugging Face 文件/元数据；**[论文]** PlantCAD2 bioRxiv 全文（v1 页面，WebFetch 摘要）；**[GH]** GitHub README/代码；**[文档]** 本目录 08/09/10 号文件里的既有实测；**[推断]** 我的推断。

---

## 0. 结论先行

1. **能跑，环境几乎零成本。** 全部依赖有预编译 wheel（torch 2.5.1+cu124、mamba-ssm 2.3.0、causal-conv1d 1.6.0），整个 env 70 秒装完，无需 nvcc/编译；加载、前向、取逐碱基 embedding、反向互补等变性、MLM 精度复核全部通过。
2. **推荐：把 PlantCAD2-Small 作为第二骨干的候选，先做"冻结编码器 + 线性/小头"的试点；Medium 只做 DEV 级验证；Large 不推荐。** 理由与代价如下。
   - Small（实际 8.8 千万唯一参数）的编码吞吐 ≈ 30 kbp/s（共享 GPU 下 [实测]），与我方现用 GENERanno 编码路径（30.7 kbp/s [文档]）同量级，DEV 53.4 Mb 编码约 0.7 h。
   - Medium 单位碱基算力约为 GENERanno 的 5 倍：DEV 编码约 2.4 h，独立面板（8.6 Gb）仅编码就约 390 h（共享 GPU），不可接受。Large 更慢（约 680 h）且训练显存放不下。
   - 单碱基分辨率带来结构上的简化：每碱基直接一个 1536（Small）维特征，不再需要 6-mer 内位置的 DNA 卷积/BiGRU 补丁；且模型对反向互补严格等变，一次前向同时给出双链特征。
3. **主要风险（详见 §6）**：(a) 预训练只见过 8,192 bp，窗口 >16 kb 时 Small 的 MLM 精度明显下降，我方 24,576 bp 联合窗口方案不能直接搬；(b) 训练单步编码开销比 GENERanno 大 14–33 倍，LoRA 微调基本不现实，只适合冻结编码器；(c) `from_pretrained(torch_dtype=bfloat16)` 对该自定义代码**并不真正转 bf16**，必须用 autocast（3.3 倍加速，数值几乎无损）；(d) 拟南芥/水稻都在 PlantCAD2 的 65 个预训练物种内，且预训练窗口取自"已注释基因 ±5 kb"，所以"无注释监督"只能说无标签监督，不能说无注释信息；(e) GeneCAD 已用同一骨干做 ab initio 注释，需要定位差异（§1.3）。

---

## 1. 官方发布信息

### 1.1 模型与规模

| 项 | 内容 | 来源 |
|---|---|---|
| 仓库 | `kuleshov-group/PlantCAD2-{Small-l24-d0768, Medium-l48-d1024, Large-l48-d1536}`（2025-03 创建）；模型卡 README 只有 license 一行，没有使用说明 | [HF] |
| 许可 | Apache-2.0（HF 元数据；GitHub `plantcad/plantcad2` 也是 Apache-2.0；GeneCAD 仓库 Apache-2.0） | [HF][GH] |
| 架构 | Caduceus 风格：双向 Mamba2（`mamba_fwd` + `mamba_rev`，in/out 投影权重绑定）+ RCPS（反向互补奇偶共享，等变），`d_state=64, expand=2, d_conv=4`，`d_intermediate=0`（无 MLP），RMSNorm，fused add-norm | [HF] config + `modeling_caduceus.py`（我已通读） |
| 规模（论文口径，绑定后的唯一参数） | S 88 M / M 311 M / L 694 M（论文摘要另写总计 676 M，与正文 694 M 不一致，原文自相矛盾） | [论文] |
| 规模（HF 检查点实际张量） | 文件 704 MB / 2.48 GB / 5.55 GB（fp32）。**检查点把绑定的 rev 投影存了两份**（我核对：Small 第 0 层 `mamba_fwd.in_proj.weight` 与 `mamba_rev.in_proj.weight` 逐元素相等；去掉重复得 88.24 M，与论文 88 M 吻合）。 | [实测] |
| 加载后的参数量 | 用 transformers 4.46.3 `from_pretrained` 后**绑定被打断**，模块内共 176.0 M / 621.0 M / 1,386.8 M 个独立参数，值相同。**推理无影响；若要微调必须重新绑定**，否则可训练参数翻倍且偏离预训练约束 | [实测] |
| 隐藏状态 | 末层输出 `(B, L, 2·d_model)`：Small 1536、Medium 2048、Large 3072 维。`output_hidden_states=True` 给 `n_layer+1` 个状态（输入 + 各层 + 末层归一化后） | [实测] |
| 上下文 | **预训练窗口 8,192 bp**（步长 4,096）；SSM 无位置编码上限，但超过 8 kb 的质量见 §3.4 | [论文][实测] |
| 分词 | 单核苷酸；`PreTrainedTokenizerFast` 的 BPE 词表只有字符：`[PAD]=0,[MASK]=1,[UNK]=2,a=3,c=4,g=5,t=6`，带 Lowercase 规范化，无 BOS/EOS（`post_processor=null`）；`N` 与其他字符→`[UNK]`；词表对齐到 8。我直接用查表映射并核对与官方 tokenizer 的 id 逐位一致 | [HF][实测] |
| 预训练目标/数据 | MLM，15% 掩码（80/10/10）；全局 batch 2,048、240k 步、AdamW lr 2e-4 余弦；**65 个被子植物基因组（每属一个）**；序列取"每个已注释基因 ±5 kb"再切 8,192 bp 窗（步长 4,096）；重复序列用 dCNS 标注后降权；Solanaceae 整科排除（番茄等做外部测试） | [论文] |
| 注释监督 | **无标签监督**（纯 DNA 掩码预训练），但**窗口位置由注释决定**（基因 ±5 kb），等于使用了注释坐标做采样；剪接位点/起止密码子只在零样本评测里使用 | [论文] |
| 发表 | bioRxiv 2025.08.27.672609（v4）；文档 10 记录已发表于 Cell Genomics（DOI 10.1016/j.xgen.2026.101329，期刊名由 DOI 前缀推断，我未独立核实） | [文档] |

### 1.2 预训练物种（与我们的物种/面板的关系）

我用 `plantcad/plantcad2-c4096` 数据集（HF datasets-server 统计，`assembly` 列恰好 65 个取值，与论文 65 个基因组一致；这是同一预训练语料的 4,096 bp 派生版，推断与 8 kb 版同源）取得完整名单，与文档 10 附录 A 的名单一致。官方仓库的 `mlm_results_pcv2-*.txt` 还按这 65 个 assembly 报告了 MLM 精度，其中含 `Athaliana`、`Osativa`。

- **拟南芥（Athaliana）、水稻（Osativa）、玉米、大豆、大麦、小麦、毛果杨、Setaria viridis、Coffea arabica 都在预训练集内。** 所以我们的 TRAIN/DEV（拟南芥+水稻）对 PlantCAD2 与对 Carbon-A 一样是**样本内**，只能作开发信号。
- 独立面板（文档 10 §1.7）里的 12 个物种**物种级都不在**预训练集内；近缘泄漏：Chenopodium album ↔ 预训练含 C. quinoa（同属）；Zizania/Themeda/Paspalum saurae 等有同族/同属近缘（Oryza、Sorghum/Zea/Miscanthus、Paspalum vaginatum）。Dipteris（蕨）与 Acorus 无近缘，是干净的分布外测试；Amborella 无同目物种，但预训练含 Nymphaea、Cinnamomum 等基部类群。
- PlantCAD2 只含被子植物，蕨类/裸子植物完全分布外。
- GeneCAD 论文的 5 个"留出物种"里 Coffea arabica 在 PlantCAD2 预训练集内（留出只对 GeneCAD 的监督头成立，对骨干不成立）[推断，据名单]。

### 1.3 GeneCAD 的方法、结果，以及与我们框架的异同

方法 [论文 bioRxiv 2025.10.31.685877 全文摘要；GH `plantcad/genecad` README]：

- **冻结**的 PlantCAD2（论文用 Small 和 Large，嵌入 1536/3072 维）→ 线性投影到 640 维 + 128 维碱基嵌入 → 8 层 ModernBERT（隐藏 768、8 头）→ 三条特征流（转录边界、剪接位点、编码相位）经门控 MLP 输出 17 类逐碱基标签（CDS/intron/5′UTR/3′UTR 各 BILU + Outside）→ **染色体级线性链 CRF**，转移矩阵固定（由拟南芥+水稻估计），用于惩罚非法特征顺序和剪接相位违规；v0.4.0 起加入序列感知（起止密码子、供体/受体基序、最小内含子/外显子长度）。
- 8,192 bp 窗口、50% 重叠；训练 5 个物种（拟南芥、大豆、大麦、水稻、毛果杨），转录本按 MMLR 基序打分 >0.5 过滤（保留 60–90%）；只用拟南芥+水稻的两物种版本只略差。
- 评测 5 个留出被子植物（Juglans regia、Coffea arabica、Phaseolus vulgaris、Nicotiana sylvestris、N. tabacum）：平均碱基 F1 0.76、外显子 F1 0.70、CDS F1 0.58、转录本 F1 0.44；摘要称转录本 F1 比 Helixer、BRAKER3 平均高 8–10%；消融：去 CRF 约 −5% CDS F1，去 PlantCAD2 约 −30%。之后可用 ReelProtein（蛋白语言模型）过滤。
- 现行仓库（v0.6.0）的默认 plant 模式已换成别的底座 `emarro/pcad2-200M-cnet-baseline`（H-Net 型，1.03 GB，无模型卡），PlantCAD2 Small/Large 只作"legacy"保留；A100 上整条管线 17.7 kbp/s [GH]。**所以"GeneCAD 的结果"与"我们能拿到的 PlantCAD2 S/M/L"不是同一个底座。**

与我们框架的异同：

| | GeneCAD | 我方框架 |
|---|---|---|
| 骨干 | 冻结 PlantCAD2（S/L）或其后续 H-Net 变体 | 骨干无关；现用 GENERanno，拟加 NT-v2、PlantCAD2 |
| 逐碱基头 | ModernBERT 8 层 + 17 类 BILU | 小头（线性/少量参数），15 个角色 + 粗区域观测通道 |
| 结构化层 | **点估计**：固定转移矩阵的 CRF + Viterbi | **概率化**：完整语法的配分函数，损失 = log Zall − log Zobs（对所有兼容基因集边缘化），精确 MAP，目标产出后验 |
| 标签 | 每位点一个规范转录本，MMLR 过滤伪标签 | isoform 并集/主转录本 + 未知位置边缘化 |
| 后处理 | ReelProtein 蛋白过滤 | 后验排序/校准（拟做） |
| 定位 | "PlantCAD2 + CRF"的注释工具 | "骨干无关的概率化 ab initio 框架"，PlantCAD2 只是其中一个可换的编码器 |

**定位含义**：只要论文的核心主张是"骨干无关 + 全语法边缘化训练 + 后验/校准"，PlantCAD2 作为第二骨干不会与 GeneCAD 撞车；但如果结果表里只有 PlantCAD2 一个骨干，审稿人会把它读成"GeneCAD 的另一种解码器"。GeneCAD 应作为必须跑的基线（`genecad predict` 可直接在 ont 上复现，Apache-2.0）。

### 1.4 依赖

| 项 | 要求 | 本次实测 |
|---|---|---|
| mamba-ssm | 代码用到 `Mamba2`、`mamba_ssm.modules.block/mha/mlp`、`ops.triton.layer_norm`：需要 **≥2.2（Mamba2 版本）**；官方未写版本 | 2.3.0 预编译 wheel（`cu12torch2.5cxx11abiFALSE-cp311`），可用 |
| causal-conv1d | 可选（有则用融合卷积） | 1.6.0 wheel，可用 |
| torch / triton | 与 wheel 匹配 | 2.5.1+cu124 / 3.1.0（与 generanno env 同版本；Mamba2 的 triton SSD 核在 triton 3.1 上正常） |
| transformers | config 记录 4.46.0；自定义代码依赖 `PreTrainedModel` 内部接口 | 固定 4.46.3，未试其他版本（5.x 有较大概率不兼容 [推断]） |
| numpy | — | <2 |
| CUDA | 用 wheel 无需 nvcc；驱动需支持 cu12 | ont 驱动 595.58（CUDA 13.2）无问题；ont 的 `nvcc` 不在 PATH，源码编译未尝试 |
| GPU | 单卡 | RTX 4090 24 GB（Ada，支持 bf16） |

编译耗时与问题：**未发生编译**（预编译 wheel 命中）。装机总耗时 torch 28 s + 其余 wheel 约 40 s。若将来 torch 版本变化导致 wheel 不匹配才需要编译（mamba-ssm 源码编译通常要 20–60 分钟，且需要 nvcc；本次未测，仅为一般经验）。

---

## 2. ont 上的环境与下载

- 环境：`~/envs/plantcad2`（uv venv，python 3.11.15），清单：torch 2.5.1+cu124、triton 3.1.0、mamba-ssm 2.3.0、causal-conv1d 1.6.0、transformers 4.46.3、huggingface_hub 0.36.2、numpy 1.26.4、einops、pandas、pyfaidx、biopython。安装脚本与日志：`~/plantcad2_work/setup_env.sh`、`setup_env.log`。
- 模型（`~/plantcad2_work/models/`）：Medium 2.4 GB、Small 672 MB、Large 5.2 GB，均 <10 GB，直接下载，约 1 分钟/个。ont 磁盘剩余 2.2 TB。
- 数据：拟南芥 `~/coding-rna/data/m1_screen/arabidopsis_thaliana/genome.fa`（只读，自己写流式读取，没有生成 .fai）。

---

## 3. 验证结果（拟南芥 NC_003074.8，窗口取自 3–6 Mb 区段）

### 3.1 加载与正确性

| 检查 | 结果 |
|---|---|
| `AutoModel`/`AutoModelForMaskedLM`（`trust_remote_code`，我已通读 `modeling_caduceus.py`，另读 `modeling_rcps.py` 前 140 行）加载 | 成功，无缺失键（MLM 头也在检查点里，基座加载时多出 `lm_head.*` 两个无害的 unexpected 键） |
| 自定义分词 vs 官方 tokenizer | 1,000 bp 逐位一致；`N`→2 |
| 逐碱基 embedding 形状 | 输入 8,192 bp → `(1, 8192, 2048)`（Medium）；无子词合并，一个碱基一个向量 |
| MLM 精度复核（15% 全部换 `[MASK]`，拟南芥 chr3，8 kb 窗口） | Medium 0.587/0.615/0.617（三批不同窗口），Small 0.545/0.563，Large 0.607。官方 `mlm_results` 里拟南芥 Medium 为 0.638（80/10/10 掩码、不同窗口）。同量级，说明权重与代码路径加载正确 |
| bf16 的真实性 | **`from_pretrained(torch_dtype=bfloat16)` 之后所有参数仍是 fp32**（818 个张量全 fp32）。官方 `inference_speed/main.py` 也是这样写的，其"bf16"测速很可能也是 fp32，**但我没有核实他们的结果**。要真正用 bf16 必须 `torch.autocast` 或 `.to(torch.bfloat16)` |
| autocast bf16 vs fp32 | Medium 逐位置余弦 0.99997（最小 0.993）；Small 0.99999；Large 0.99999；MLM 精度不变（Medium 0.6169 vs 0.6169，Small 0.5635 vs 0.5630）。纯 `.to(bf16)` 略差（Medium 最小余弦 0.95），不推荐 |
| 数值/NaN | 无 NaN；输出标准差 0.45、绝对值最大 10.8（Medium 末层） |

### 3.2 双链 / 反向互补

- RCPS 的含义：隐藏向量前一半是"+ 链视角"，后一半是"把反向互补序列送入同一网络再翻回来的视角"，两半之间不混合（它们逐位置余弦仅 0.034，近似正交，携带不同信息）。
- **严格等变** [实测，Medium]：`f(RC(x)) == flip(f(x), 位置与通道两个维度)`，最大绝对差 0.0，余弦 1.0。含义：**一次前向已经同时给出两条链的特征**，不需要像 GENERanno/NT-v2 那样对反链再编码一次；对 − 链只需把位置翻转并交换两半通道，然后用同一个头。
- 这等于对"两条链"的编码成本折半，但每次前向本身比单链模型贵一倍多（每层 4 次 Mamba2 扫描：正反方向 × 两条链视角）。

### 3.3 吞吐、显存（RTX 4090 共享，autocast bf16，forward only，B=1，均值/中位）

测试期间 ont 上有另外 2–3 个进程占着 GPU（利用率始终 100%，占用 7–17 GB），我用 8192×8192 bf16 矩阵乘测得 52–67 TFLOPs，约为 4090 独占可达水平的 40–50% [推断]。**下表是共享条件下的数字，独占时可能再快 1.5–2.5 倍，但没有测到**；与 GENERanno 的文档数字同样是在他人占卡时测得，口径一致。大批量和长窗口的 OOM 有一部分是他人显存占用造成的，不是模型固有上限。

| 模型 | 窗口 | 耗时 | 吞吐 | 峰值显存 |
|---|---:|---:|---:|---:|
| Small（88 M 唯一参数） | 8,192 | 0.32 s | 25.6 kbp/s（B=4：30.2） | 1.3 GB（B=4：2.5） |
| | 16,384 | 0.49 s | 33.5 | 1.7 GB |
| | 24,576 | 0.76 s | 32.3 | 2.1 GB |
| | 32,768 / 65,536 / 98,304 | — | 35.9 / 55.7 / 23.6（后两者噪声大） | 2.5 / 4.1 / 5.7 GB |
| Medium（311 M） | 8,192 | 0.96 s | 8.5 kbp/s（B=4：8.7） | 3.6 GB（B=4：5.2） |
| | 16,384 / 24,576 / 32,768 | 2.1 / 3.0 / 3.6 s | 7.9 / 8.1 / 9.2 | 4.2 / 4.7 / 5.2 GB |
| | 65,536 / 98,304 | 7.6 / 12.0 s | 8.7 / 8.2 | 7.3 / 9.5 GB |
| Large（694 M） | 8,192 | 1.64 s | 5.0 kbp/s | 7.8 GB |
| | 16,384 / 24,576 | 3.5 / 5.3 s | 4.7 / 4.7 | 8.6 / 9.4 GB |
| | B=4×8k、≥32,768 | OOM（他人占用下） | — | — |

附：fp32（关 TF32）→ autocast bf16 的加速（8 kb，B=1）：Small 12.1→32.4 kbp/s、Medium 3.2→10.5、Large 1.7→6.3（约 3.3–3.8 倍）；profile 显示 fp32 路径 66% 的 GPU 时间花在 `ampere_sgemm` 上，即投影层的 FP32 矩阵乘。吞吐几乎与窗口长度无关（线性），显存随长度线性增长，所以长度上限由质量而非算力决定。

### 3.4 上下文长度的质量（重要）

MLM 精度（拟南芥 chr3，15% 掩码，不同长度取不同窗口，未转 bf16）：

| 窗口 | Small | Medium | Large |
|---:|---:|---:|---:|
| 8,192（预训练长度） | 0.545 | 0.587 | 0.607 |
| 16,384 | **0.437** | 0.583 | 0.611 |
| 32,768 | 0.452 | 0.551 | OOM |
| 65,536 | 0.376 | 0.416 | OOM |

- Small 一旦超过 8 kb 就明显退化；Medium/Large 到 16 kb 基本保持；Medium 32 kb 开始掉点，64 kb 掉得很多。
- 另一个视角（Medium）：同一位置在 8 kb 窗内与在 32 kb 窗内的 embedding 余弦只有 0.89–0.99（距窗口边缘 2–4 kb 处最低 0.886；0–64 bp 处 0.957）。embedding 本身不是窗口长度不变的，头必须在与推理相同的窗口长度上训练。
- **后果**：我方现行 24,576 bp 联合窗口、12,288 步长的方案对 Small 不成立；应改用 ≤8 kb（Small）或 ≤16 kb（Medium）的窗口加重叠平均/中心裁剪，且头在同一窗口长度下训练。

---

## 4. 在完整语法框架里作为编码器的粗估

对比基线口径：GENERanno 的编码路径 **30.7 kbp/s**（唯一碱基、两链、50% 重叠）[文档 09 §推理路径：纯编码器前向 0.202 s/窗，24,576 bp 窗，步长 12,288，每碱基两链×两倍重叠]；整条管线 17.7 kbp/s，其中编码占 76%，skyline 解码占 24%（≈6.7 s/Mb/链，与骨干无关）。NT-v2-500M 我在同一机器上实测（HF 自定义代码、未用 SDPA，未经优化）：12,288 bp（2,048 token）窗，B=1 为 42.8 kbp/s，B=4 为 57.9 kbp/s（仅单链单次前向）。

为得到"每个基因组碱基覆盖一次所需时间"，PlantCAD2 取：Small 用 8,192 窗/6,144 步长（25% 重叠），Medium/Large 用 16,384 窗/12,288 步长；NT-v2 按两链各一次 + 50% 重叠（÷4）。吞吐取上表的共享 GPU 实测。

| 骨干 | 有效编码速率（唯一碱基，双链） | DEV 53.4 Mb 编码 | TRAIN 157 Mb 编码 | 独立面板 8.6 Gb 编码 |
|---|---:|---:|---:|---:|
| GENERanno（现用） | 30.7 kbp/s [文档] | 0.48 h | 1.4 h | 78 h |
| NT-v2-500M（未优化） | ≈14.5 kbp/s | 1.0 h | 3.0 h | 165 h |
| PlantCAD2-Small | ≈22–25 kbp/s | 0.6–0.7 h | 1.8–2.0 h | ≈100 h |
| PlantCAD2-Medium | ≈6.1 kbp/s | 2.4 h | 7.1 h | ≈390 h |
| PlantCAD2-Large | ≈3.5 kbp/s | 4.2 h | 12.4 h | ≈680 h |

另加解码约 0.2–0.3 h/DEV（共同项，文档 09 的 skyline 复用优化可减一半以上）。以上均为"共享 GPU"数字，独占时乐观可再快 1.5–2.5 倍 [推断]。物理解释：Medium 每碱基矩阵乘约 2×311 M×4 次应用 ≈ 2.5 GFLOP/bp，GENERanno 约 2×1.15 B/6 ≈ 0.38 GFLOP/bp，单碱基 + 双向 + 双链的代价在这里体现。

### 训练单步编码耗时

按"一步 = 3 个 24,576 bp 窗，共 73.7 kbp"换算。我方现状：头训练（冻结编码器）约 10.2–10.4 s/步（3090），在 ont 上约快 1.4 倍；GPU 只忙约 11%，瓶颈是 CPU 上的 FP64 配分函数 [文档 09]。GENERanno LoRA 编码器前向 0.67 s + 反向 0.17 s [文档 09]。

| 编码器 | 仅前向（冻结头训练） | 前向+反向（LoRA 类，经冻结权重回传到输入） | 说明 |
|---|---:|---:|---|
| GENERanno 1.15 B | 0.67 s | 0.84 s | [文档] |
| NT-v2-500M | ≈1.3–1.7 s | ≈4–5 s（估计，训练实测因他人占显存 OOM） | 12 kb 窗口，需 6 个窗覆盖 |
| PlantCAD2-Small | ≈2.5 s | **≈12 s**（实测 24,576 bp 单窗 4.05 s ×3） | 显存 7.9 GB（逐层检查点） |
| PlantCAD2-Medium | ≈8.7 s | **≈28 s**（实测仅 8,192 bp：3.08 s，2.66 kbp/s） | 8 kb 显存 8.0 GB；24 kb 线性外推约 19 GB，未能在共享 GPU 上实测 |
| PlantCAD2-Large | ≈16 s | ≈60 s（实测仅 4,096 bp：1.19 kbp/s，9.6 GB） | 8 kb 以上在共享 GPU 上 OOM |

- 该模型的 `supports_gradient_checkpointing = False`，我在每层外自己包了 `torch.utils.checkpoint`（上表即此配置）。**不带检查点基本放不下**：Small 在 4 kb、B=1 不用检查点就要 8.8–9.9 GB（约 2.2–2.4 GB/kbp）；Medium 4 kb 不用检查点在共享 GPU 上 OOM。
- 全参数反向（含参数梯度）与"冻结权重回传输入"耗时相近（Medium 8 kb：3.20 s vs 3.08 s）；Adam 状态另需约 2×权重。
- 结论：**冻结编码器 + 训练小头是现实路线**。因为我方步长本来由 CPU DP 决定（≈7 s/步 ont），Small 的前向开销 ≈2.5 s（约 +35%），Medium 的 ≈8.7 s（约 +120%）；LoRA/全参微调对 Medium 以上不现实，Small 也要 ≈12 s/步（约 +170%），且需要逐层检查点。
- **不能缓存全分辨率 embedding**：Small 每碱基 1536 维 fp16 = 3 KB，TRAIN 157 Mb ≈ 0.47 TB，DEV ≈ 0.16 TB；Medium 4 KB/碱基。需要在线编码，或缓存线性降维后的特征（须先验证不掉点）。

---

## 5. 与另一个第二骨干 NT-v2、与现用 GENERanno 的比较

| 维度 | GENERanno / Carbon-A | NT-v2-500M | PlantCAD2-Small | PlantCAD2-Medium |
|---|---|---|---|---|
| 预训练数据 | RefSeq 基因区为主；Carbon-A 为注释监督的全模型 | 850 物种全基因组，无注释 | 65 个被子植物，窗口取基因 ±5 kb，无标签 | 同左 |
| 分辨率 | 6-mer | 6-mer | 单碱基 | 单碱基 |
| 上下文 | 16,384 token ≈ 98 kb（用 24 kb 窗） | 2,048 token ≈ 12 kb | 8,192 bp（质量上限 ≈8 kb） | 8,192 bp（质量上限 ≈16 kb） |
| 唯一参数 | 1.15 B | 0.5 B | 88 M | 311 M |
| 双链 | 需两遍 | 需两遍 | 一遍，严格等变 | 一遍，严格等变 |
| 编码速率（DEV 口径） | 30.7 kbp/s | ≈14.5（未优化） | ≈22–25 | ≈6.1 |
| 子 token 位置问题 | 需 DNA 卷积/BiGRU 补丁 | 同左 | 无 | 无 |
| 物种泛化风险 | 全物种 | 全物种 | 仅被子植物；蕨类分布外 | 同左 |
| 许可 | MIT / Apache-2.0 | CC BY-NC-SA 4.0（HF 元数据；非商业 + 相同方式共享，与 PlantCAD2 的 Apache-2.0 相比更受限，发布代码/权重衍生物时需注意） | Apache-2.0 | Apache-2.0 |

---

## 6. 风险与对策

| # | 风险 | 严重度 | 依据 | 对策 |
|---|---|---|---|---|
| R1 | 算力：单碱基 + 双向 + 双链使每碱基开销大。Medium 面板编码约 390 h（共享 GPU）；Large 不可行 | 高（Medium/Large），低（Small） | §4 | 主线只用 Small；Medium 仅做 DEV 的"规模消融"；面板用 Small，Medium 只跑小基因组子集 |
| R2 | 上下文：8 kb 预训练，Small >8 kb 退化；方案里的 24 kb 联合窗不能用 | 中 | §3.4 | 窗口 8 kb、步长 6 kb、重叠处平均；头与推理同窗口训练；新增一个 8 kb 窗口版 DEV 对照 |
| R3 | 训练开销：冻结可行，LoRA/全参对 Medium+ 不可行；模型无原生检查点支持 | 中 | §4 | 先做冻结头；如要 LoRA，仅 Small + 逐层检查点（≈12 s/步） |
| R4 | 数据暴露：拟南芥/水稻均在预训练内，窗口取基因 ±5 kb，实质含注释信息（无标签） | 中（叙事） | §1.1–1.2 | 论文中写"无标签监督、注释坐标采样"，不写"完全非注释"；以独立面板为主证据；面板内对同属近缘做分层（Chenopodium 等） |
| R5 | 定位：GeneCAD 已有"PlantCAD2 + CRF" | 中 | §1.3 | 强调骨干无关、配分函数边缘化训练、后验输出；GeneCAD 作为必跑基线；其默认底座已换成 H-Net 变体，需注明比较的是哪个版本 |
| R6 | 软件脆弱：自定义代码、绑定被打断、transformers 版本敏感、`torch_dtype` 不生效 | 低–中 | §1.4、§3.1 | 固化 env；一律 autocast；微调时重新绑定 rev 投影；把 `uv pip freeze` 存档 |
| R7 | 仅被子植物：面板中蕨类等分布外 | 低 | §1.2 | 面板内报告分层（被子 vs 蕨），不把蕨类算作 PlantCAD2 的失败 |
| R8 | 论文自身参数量口径不一致（676 M vs 694 M；HF 张量为其 2 倍） | 低 | §1.1 | 报告写"88/311/694 M（绑定后）"，并脚注 |
| R9 | 测速受他人占卡影响，上表可能低估 1.5–2.5 倍 | 中（对预算） | §3.3 | 正式预算前，在 GPU 空闲时段复测 Small 与 Medium 一次 |

---

## 7. 建议的下一步（均未执行）

1. **试点 P0（≤2 GPU-h）**：Small，8 kb 窗口/6 kb 步长，取末层 1536 维，只训练线性头（沿用现有 15 角色 + 观测通道），在 TRAIN 一条染色体上训练、DEV 拟南芥整条评估，对比同口径 GENERanno 冻结头。判据：完整链 F1 能否进入 GENERanno 冻结头 ±3 点以内；同时确认头在 8 kb 窗口版本不掉点。
2. 若 P0 达标，再做层选择探针（`output_hidden_states` 的 25 层，成本为一次前向）和 Medium 的 DEV 单点对照，决定是否值得为面板花 390 h。
3. 与 NT-v2 试点结果合并，形成"≥2 个无注释监督的骨干"的论述；GeneCAD 在 ont 上跑一遍 DEV 与面板小基因组，建立基线。
4. 在 GPU 空闲时复测吞吐与 Medium 24 kb 训练显存，并在 Baobab 3090 上验证 wheel 可用性（本次未触碰 Baobab 的任何环境）。

---

## 8. 未验证事项

- 没有做任何整染色体推理、没有训练头、没有 F1；所有 F1 级别的判断都未做。
- 吞吐为共享 GPU 条件下的单次中位数，噪声大（例如 Small 65 kb 的 55.7 与 98 kb 的 23.6 kbp/s）；"独占加速 1.5–2.5 倍"是由矩阵乘探针推断的，未实测。
- Medium 24 kb 与 Large 8 kb 的训练显存因他人占用未能实测，上表为线性外推或空缺。
- NT-v2 训练耗时未实测（OOM）；NT-v2 推理为未经优化的 HF 自定义代码，会低估其可达速度。
- 预训练物种名单取自 `plantcad2-c4096`（HF datasets-server 统计）和文档 10 附录 A，二者一致，但论文补充表 1 我没有直接读到。
- PlantCAD2 论文页面只读到约前 100 k 字符；"期刊发表于 Cell Genomics"沿用文档 10 的推断。
- GeneCAD 的结果数字来自预印本文本摘要；竞争对手的具体数字只在图里，未读。
- 第一轮上下文一致性和 RC 等变性的测量用的是 fp32 权重（当时误以为是 bf16，后来用 profile 发现）；autocast 下未单独复测 RC 等变性，仅复测了数值余弦和 MLM 精度。

---

## 9. 复现与资产（均在 ont）

- 环境：`~/envs/plantcad2/bin/python`；安装：`~/plantcad2_work/setup_env.sh`（日志 `setup_env.log`）。
- 模型：`~/plantcad2_work/models/PlantCAD2-{Small-l24-d0768,Medium-l48-d1024,Large-l48-d1536}`。
- 脚本（均绑核 4–9）：`common.py`（FASTA 流式读取、分词、加载）、`phaseA.py`（正确性、RC 等变、上下文一致、MLM）、`phaseA2.py`（多长度 MLM）、`phaseC.py`（fp32/TF32/autocast/纯 bf16 对比）、`phaseD.py`（autocast 推理吞吐与训练开销网格）、`ntv2.py`（NT-v2 推理）、`prof.py`（profiler）、`tie.py`（权重绑定核对）、`ref.py`（矩阵乘算力探针）。结果在同目录 `phase*.log`。
- 调用要点：`AutoModel.from_pretrained(dir, trust_remote_code=True)` 后保持 fp32，推理时 `with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16): m(ids).last_hidden_state`；输入 id = `{a:3,c:4,g:5,t:6,其他:2}`；`last_hidden_state[..., :d]` 为 + 链视角，`[..., d:]` 为反向互补视角（对应反链时翻转位置并交换两半）。

## 10. 外部来源

- PlantCAD2 论文：https://www.biorxiv.org/content/10.1101/2025.08.27.672609v1.full ；HF 模型 https://huggingface.co/kuleshov-group/PlantCAD2-Medium-l48-d1024 （及 Small/Large）；预训练语料派生集 https://huggingface.co/datasets/plantcad/plantcad2-c4096 ；官方评测仓库 https://github.com/plantcad/plantcad2 （`pipelines/inference_speed/`）。
- GeneCAD：https://doi.org/10.1101/2025.10.31.685877 ；https://github.com/plantcad/genecad ；默认底座 https://huggingface.co/emarro/pcad2-200M-cnet-baseline 。
- 预编译 wheel：https://github.com/state-spaces/mamba/releases/tag/v2.3.0 ；https://github.com/Dao-AILab/causal-conv1d/releases/tag/v1.6.0 。
