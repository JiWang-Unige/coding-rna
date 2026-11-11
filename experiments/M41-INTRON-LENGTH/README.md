# M41-INTRON-LENGTH — 长、短内含子状态对照实验

协议：`docs/review_20261016/33_M41_intron_length_states_protocol.md`（用户决定 `intron-length-states`，2026-11-11）。
一个对照臂 T2_ILS = T2 同数据、同配置、同步数、同选点规则，只改语法（r = 1：每个相位、每段历史各有一个短内含子状态和一个长内含子状态；形状规范长度项；3 个可学标量）。
冻结代码（pipeline/m28、M29 / M31 / M32 / M34）只被导入，没有改动；状态关闭时与冻结代码逐位一致（有测试）。

| 文件 | 作用 |
|---|---|
| `scripts/m41_grammar.py` | 自动机（冻结 `_transitions` + I3 拆成 I3S1..r / I3L）、FP64 配分函数与对 logits 和 bonus 的梯度、粗通道损失、skyline 精确解码 |
| `scripts/m41_prior.py` | θ → bonus 向量，闭式长度项（供测试），训练染色体内含子长度的最大似然拟合 |
| `scripts/m41_intron_lengths.py` | 从 M34 catalog 提取 Ath / Osa 训练染色体的内含子长度（ont，CPU） |
| `scripts/m41_dp.py` | trainer 用的并行窗口 DP（沿用 `a1_dp` 的 canary） |
| `scripts/m41_train.py` | M34 trainer 主循环 + [M41] 行（θ 的优化器、检查点里的 `intron` 项） |
| `scripts/m41_infer.py` | `m34_infer` / `m34_p2b_infer`（--setb），解码进程换成 M41 解码器 |
| `scripts/m41_select.py` | 选点与 DEV 一次性哨兵（规则同 `m34_select`） |
| `scripts/m41_compare.py` | DEV 配对块 bootstrap（合并 + 1 / 2–5 / ≥6），给出预注册判定 |
| `scripts/m41_stage.sh` | 对照臂全流程：训练 → 验证 → 选点 → DEV 一次 → B 集 |
| `scripts/m41_smoke.sh` | GPU 冒烟：20 步 + 2 Mb 验证片段，不碰 DEV / B 集 |
| `scripts/m41_dryrun.py` | CPU 干跑：用已存的 T10 势能比较冻结解码与 M41 解码 |
| `results/` | `train_intron_lengths.json`、`intron_prior_fit.json`（r = 1..4）、`dryrun_NC_141366.1_0-4M.json` |

测试（CPU，约 10 s）：`python -m pytest experiments/M41-INTRON-LENGTH/tests -q`（ont：`~/envs/generanno/bin/python`）。
推送到 ont：`./push.sh`。GPU 冒烟和对照臂都由维护会话排队（M39 之后）：`bash scripts/m41_smoke.sh`，然后 `bash scripts/m41_stage.sh`。
进度：`tail -n 3 ~/coding-rna/outputs/M41-INTRON-LENGTH/pc2_T2_ILS/training.jsonl`
