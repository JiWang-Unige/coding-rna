# experiments/M28：full-grammar 非最终分支与诊断

最终 M28 full-grammar 管线在 [`pipeline/m28`](../../pipeline/m28/README.md)。这里是围绕它的对照分支、诊断、门禁与协议报告：

```
full_grammar/
├── scripts/<bucket>/    # 按分支/主题分桶的驱动脚本（bucket 见下）
├── sbatch/              # 对应的 M28-FULL-GRAMMAR-* 作业（207 个左右，平铺）
├── tests/               # 分支专用测试（LoRA/dense-QV/qvo、E/M missed-reference、train-recovery、qv-lr）
├── reports/             # full_grammar_* 协议、计划、结果（含被 checkpoint 引用的 PROTOCOL 标识对应的文档）
└── src/m28/             # 仅分支使用的库模块（joint_encoder_qvo_lora）
```

`scripts/` 桶：`objective_contrast`、`native_cds_prior`（含 joint24_native_prior / native context）、`complete_core_draws`、`adaptation_lora_qv`（joint-LoRA / dense-QV / qvo / lr 对照）、
`e_m_ensemble`（E/M 能量池、minimum-centered、matched/train energy）、`gene_posterior`（基因后验、scaled-mass、posterior mass）、`context_prior`（采样/上下文 DEV 家族与诊断）、
`dev_families`（continuation / fixed-DEV 等较早的 DEV 家族）、`train_support_recovery`（支持性/恢复/观测/粗通道普查）、`partition_probes`（配分函数/viterbi/skyline 的合成与规模探针）、`local_transitions`（局部链/内含子编辑转移与原 TRAIN 链/误差画像）。

pre-full-grammar 的 M28（A3–A44、pilot、B1/B2、R2–R6、旧 lattice 变体、full-train C0 等）见 [`archive/m28_pre_full_grammar/`](../../archive/README.md)。
