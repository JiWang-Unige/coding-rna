"""fig_m42  (SLOT_WIDE 1664x560)  slide c2-m42

M42 gate 1 on Ziziphus jujuba evaluation chromosomes E (10,266 reference genes, 675 missing in the old release):
overall exact-chain F1 and missing-locus recall of Z (T2 zero-shot), U, F, U3 and the old annotation itself.
Source: experiments/M42-ORACLEFREE-C2/results/gate1_score.json, checked against protocol 35 §8.3 table.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

P35 = "docs/review_20261016/35_M42_oracle_free_C2_protocol.md"
r = json.load(open(os.path.join(style.REPO, "experiments/M42-ORACLEFREE-C2/results/gate1_score.json")))["per_run"]
ARMS = [("Z", "Z  zero-shot T2", style.GRAY), ("U", "U  marginalise (V1)", style.BLUE),
        ("F", "F  fill background", style.ORANGE), ("U3", "U3  all unknown", style.BLUE_LIGHT)]
for k, lab in (("Z", "| Z（T2，零样本）"), ("U", "| U（V1 部分观测）"), ("F", "| F（填背景）"), ("U3", "| U3（V3，无负约束）")):
    ln = style.find(P35, lab)
    style.has(P35, ln, f"{r[k]['F1']:.2f}", f"{r[k]['novel_recall']:.2f}")
OLD = r["OLD"]["F1"]
style.has(P35, style.find(P35, "| 参照：旧版注释直接当预测"), f"{OLD:.2f}")

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
for i, (key, title, lim) in enumerate((("F1", "Overall exact-chain F1 (%)", 80), ("novel_recall", "Missing-locus recall (%, 675 loci)", 24))):
    ax = style.axes_px(fig, 90 + i * 820, 60, 700, 400)
    xs = range(len(ARMS))
    vals = [r[a][key] for a, _, _ in ARMS]
    ax.bar(xs, vals, color=[c for _, _, c in ARMS], width=0.62)
    for x, v in zip(xs, vals):
        ax.text(x, v + lim * 0.015, f"{v:.2f}", ha="center", fontsize=18, fontweight="bold", color=style.INK)
    if key == "F1":
        ax.axhline(OLD, color=style.INK, ls="--", lw=1.8)
        ax.text(3.45, OLD + 1.2, f"old annotation as prediction {OLD:.2f}", ha="right", fontsize=16, color=style.INK)
    else:
        ax.text(3.45, 22.2, "old annotation: 0 by definition", ha="right", fontsize=16, color=style.INK)
    ax.set_xticks(list(xs), [a for a, _, _ in ARMS])
    ax.set_ylim(0, lim)
    ax.set_title(title)
style.legend_row(fig, 90, 505, [("s", c, c, lab) for _, lab, c in ARMS])
style.save(fig, "fig_m42")
print("ok")
