Complete TA-DECODER-M3 (run-and-evaluate) · Track A decoder batch: CRF / semi-CRF / constrained-Viterbi on the fixed tiberius_like backbone vs same-budget screen_anchor. Follow `goals/TA-DECODER-M3_protocol.md` + `CLAUDE.md`; judge completion only by the inline checkpoints below.

⚠️ Retrospective advisory triggered (≥5 iters, no prior retrospective); advisory only, does not cancel this flow.
⚠️ FOCUSED ARCH BATCH on `decoder` (3 candidates, distinct mechanism_delta; SOFT_WARN, table in protocol). Not a hyperparam sweep.

## 权限声明
本 goal 拥有充分工具权限。若集群命令 / 网络下载 / subagent 调用临时失败, **先尝试 ≥ 1 替代路径**, 仍失败再 inline 报告, 不要停下等用户：
- 内置工具不受限（Bash / Edit / Write / Read / WebFetch / Glob / Grep / TodoWrite / BashOutput / KillShell 等）
- 可并行调用 Claude Code subagents, 但避免写同一文件
- 允许下载任何研究所需数据；下载后 inline 报告 path / version / hash / 来源（本轮不需下载新数据）
- sbatch 远程提交训练；已有 `scripts/run_screen_ref.sbatch` 仅 read-only review + Phase 1 guard, 扩 `--decoder` 不算重写
- 必要时可创建新 skill 以稳定流程（创建后 inline 报告路径+用途）
- 关键 source-of-truth 文件（CLAUDE.md / docs/03 / docs/09）只可 draft patch text, 不直接 Edit

## 运行说明
- 所有 sbatch 提交遵循 /smart-sbatch（Phase 1 guard + Phase 2 optimization）；已存在脚本只跑 Phase 1 guard, 跳过 Phase 2
- 预期 ≤ 1 天的运行用 sleep 同步等待, 不中途下线
- **任何长度的运行**: 先等第 1 个 epoch 完成, 确认 loss 在降 / eval 正常, **之后**才允许降低监控频率
- 本轮预计 ≤ 6h, 第 1 epoch 应 ≤ 0.5h 内出, 之后切 30/60/120 min 三档 polling
- 低频监控仍按 protocol Slurm polling cadence 执行, 不可完全放任

## 决策自治
goal 执行中遇到**多选项决策点**（如 pivot 7 选 1, 是否 drop semi-CRF, 哪个 decoder 晋升）, **默认不暂停等用户**, 而是：
1. inline 列出 N 个选项 + 每选项预期 (a) SOTA gap 收益 (b) 时间成本 (c) 风险
2. 立即并行调 3 CLI reviewer (Claude+Codex+Antigravity, 复用 /tri-review), prompt 含选项 + 当前 gap + 时间预算
3. 按 reviewer 共识 + cost-adjusted ROI 自动选 1 个继续 (优先级: 严格 SOTA 进展 > 时间成本 > 风险)
4. 决策 + 3 reviewer 简要理由 写 docs/08 或 inline
5. 继续 goal, **不暂停**
**例外** (仍需用户显式确认): 强破坏性操作 (scancel 他人 job / rm 数据 / 改 docs/03) / route 级 abandon (docs/09) / 决策影响 > 24h 净新增 compute / ≥ 2 reviewer 反对 default 或 tied 3-way 无 leader

## Mode & Milestone
Mode: run-and-evaluate; ends after survivors' 3-seed screen is result-logged + tri-reviewed + pivoted (if >12h, hand off remaining training).
Milestone: M3 Track A screen, CANNOT claim, profile=screen. Gate: primary_progress = seed-mean gene_body_F1_unconstrained (CDS) > 0.5676 (anchor 0.5576+0.01); sota_claim N/A; review = /tri-review→/pivot.

## Hard pre-submit gate
Each candidate must pass its sanity smoke (end-to-end, predicted_genes>0, CDS F1>0, gene_count_ratio < per-base baseline) before its full 3-seed screen. Smoke fail after bounded debug (≤3) → drop it, do not submit. semi-CRF is highest-risk.

## Required chain
1. /implement: add CRF / semi-CRF / constrained decoders on the FIXED tiberius_like backbone (extend src/screen_anchor).
2. Sanity-smoke each (tiny) + bounded debug; list survivors.
3. /smart-sbatch Phase 1 on extended sbatch; submit frozen 3-seed screen per survivor.
4. Eval each --span-mode cds → aggregate → seed-mean + seed-wise gene_count_ratio.
5. /result-log → validate_goal → /tri-review → /pivot (R5 discipline) → /exp-log.

## Completion (inline ✅ CK1-CK6)
CK1 decoders implemented on fixed backbone, survivors listed; CK2 each survivor smoke passed (predicted_genes>0, ratio < per-base baseline); CK3 frozen 3-seed screen per survivor (same protocol as anchor); CK4 per-candidate seed-mean base-w gene_body_F1_unconstrained (CDS) + seed-wise gene_count_ratio reported vs 0.5576 / gate 0.5676; CK5 tri-review + pivot recorded (promote/rebatch/change-axis) with R5; CK6 docs/04+05+06+08 + exp-log updated, iter_ledger clean.

## Constraints
- Backbone fixed (tiberius_like); ONLY decoder varies — clean decoder-axis test; same frozen protocol as anchor refs (no drift).
- Screen never claims; compare only vs screen_anchor (+0.01); NEVER vs pretrained_ceiling 0.9213.
- Report seed-wise gene_count_ratio (R5); coherence improvement is required promotion evidence, not just base-F1.
- Anti-tuning: decoders don't beat anchor & gap ≥0.05 → change axis, do NOT tune lr/batch/dropout.
