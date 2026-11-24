"""fig_scale  (SLOT_WIDE 1664x560)  slide c1-scale

Claim C1 (negative finding): training on ten species instead of two did not help.
T10 - T2 is -0.96 on DEV (P2 not met) and -0.14 on set B (P2b not met).

Sources (docs/review_20261016/):
  20_M34_stage1_results.md line 125 (T10 - T2 on DEV, 1-Mb block bootstrap, 54 blocks),
    line 147 (T10 - T2 on set B per species and pooled), lines 134-136 (P1/P2 verdicts and the
    descriptive new-species transfer +0.77), line 151 (P2b verdict)
  19_M34_scaleup_plan.md line 183 (P2b criterion, registered 24 October before any T10 output)
One seed per model; T2 = step 2304, T10 = step 3072, each by the pre-registered rule.
"""
import re
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R20 = "docs/review_20261016/20_M34_stage1_results.md"
R19 = "docs/review_20261016/19_M34_scaleup_plan.md"
dev = style.cite(R20, 125, "**−0.96**", "−1.52", "−0.29")
c = style.cells(R20, 147)
assert c[0] == "T10 − T2"
NUM = r"[−+]?\d+\.\d+"
setb = [[style._num(x) for x in re.match(rf"({NUM}) \[({NUM}), ({NUM})\]", c[i]).groups()]
        for i in (1, 2, 3)]
style.has(R20, 135, "**P2（T10 − T2 ≥ +2.0 且 CI 下界 > 0）不成立**")
style.has(R20, 151, "**判定：P2b 不成立**")
style.has(R19, 183, "判据：两物种 pooled 差 ≥ +2.0 且 CI 下界 > 0，并且两物种点估计同号为正")
newsp = style.cite(R20, 136, "69.62", "68.84", "+0.77")

rows = [  # label, (d, lo, hi), group
    ("DEV: Arabidopsis + rice (P2)", dev, "dev"),
    ("Set B pooled (P2b)", setb[2], "b"),
    ("   Dryas octopetala", setb[0], "b"),
    ("   Eschscholzia californica", setb[1], "b"),
]

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
ax = style.axes_px(fig, 470, 70, 1060, 400)
XMIN, XMAX = -2.5, 3.5
ax.axvspan(2.0, XMAX, color=style.GRID, alpha=0.7, lw=0)
ax.text(2.75, 4.55, "pre-registered\nthreshold ≥ +2.0", ha="center", va="center", fontsize=16,
        color=style.MUTED, linespacing=1.1)
ax.axvline(0, color=style.INK, lw=1.3)
ys = [3.6, 2.2, 1.3, 0.4]
for y, (label, (d, lo, hi), g) in zip(ys, rows):
    bold = not label.startswith("   ")
    ax.text(XMIN - 0.08, y, label.strip(), ha="right", va="center", fontsize=17,
            fontweight="bold" if bold else "normal",
            fontstyle="italic" if not bold else "normal", color=style.INK)
    col = style.BLUE if bold else style.BLUE_LIGHT
    ax.plot([lo, hi], [y, y], color=col, lw=3 if bold else 2.4)
    ax.scatter([d], [y], s=130 if bold else 100, color=col, zorder=3)
    txt = f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]".replace("-", "−")
    ax.text(max(hi, 0) + 0.12, y, txt, ha="left", va="center", fontsize=16.5, color=style.BODY,
            fontweight="bold" if bold else "normal")
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(-0.2, 5.0)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([-2, -1, 0, 1, 2, 3])
ax.set_xticklabels(["−2", "−1", "0", "+1", "+2", "+3"])
ax.set_xlabel("T10 (ten species) − T2 (Arabidopsis + rice), exact-chain F1 points, 95% CI",
              labelpad=8)
ax.text(-0.06, 4.55, "ten species worse", ha="right", va="center", fontsize=16, color=style.MUTED)
ax.text(0.06, 4.55, "ten species better", ha="left", va="center", fontsize=16, color=style.MUTED)
print(style.save(fig, "fig_scale"))
