"""fig_m36  (SLOT_WIDE 1664x560)  slide c2-m36

Claim C2 (boundary): with pseudo-labels or protein-evidence seeds, every adaptation arm ends
below zero-shot, and fill-background is no longer 2 points worse than marginalising (M36 NO-GO,
C2 withdrawal clause triggered for this setting).

Source: docs/review_20261016/27_M36_results.md §2 table (lines 17-23: arm F1, arm - Z with
paired 95% CI, closure), line 27 (F1 - S1 = -0.94 [-1.61, -0.28]), §3 verdict lines 31-35.
Set B (Dryas: 4 chromosomes, Eschscholzia: 3 scaffolds), exact-chain F1, evaluated once; CI =
paired bootstrap over gene units + 100 kb orphan bins, 2,000 draws, seed 0. Base model T10
(PlantCAD2 + LoRA, ten species, step 3072); each arm = 512 LoRA steps, one seed.
"""
import re
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R27 = "docs/review_20261016/27_M36_results.md"
NUM = r"[−+]?\d+\.\d+"

# ------------------------------------------------------------------ data
rows = []
for ln in range(17, 24):
    c = style.cells(R27, ln)
    m = re.match(rf"({NUM}) \[({NUM}), ({NUM})\]", c[3])
    d, lo, hi = (style._num(x) for x in m.groups())
    clos = style._num(c[5])
    rows.append((c[1].split(" ")[0], style._num(c[2]), d, lo, hi, clos))
style.has(R27, 17, "Dryas（Z 68.11，Carbon-A 过滤 73.98）")
style.has(R27, 21, "Eschscholzia（Z 68.50，Carbon-A 过滤 73.32）")
assert [r[0] for r in rows] == ["S1", "S2", "S3", "F1", "S1", "S2", "S3"]
for r, z in zip(rows, [68.11] * 4 + [68.50] * 3):
    style.close(r[1] - z, r[2], 0.011)               # arm F1 - Z equals the reported difference
assert all(r[4] < 0 for r in rows)                   # every CI entirely below zero
f1s1 = style.cite(R27, 27, "−0.94", "−1.61", "−0.28")
style.has(R27, 33, "**NO-GO：成立**")
style.has(R27, 35, "**C2 撤回条款触发**", "F1 − S1 = −0.94 ≥ −2.0")

ARM = {"S1": "S1 self-training seeds", "S2": "S2 protein-evidence seeds",
       "S3": "S3 union of S1 and S2", "F1": "F1 S1 seeds, fill background"}

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()

ax = style.axes_px(fig, 400, 70, 760, 400)
ys = [8, 7, 6, 5, 3, 2, 1]
for y, (arm, f1, d, lo, hi, clos) in zip(ys, rows):
    col = style.ORANGE if arm == "F1" else style.BLUE
    ax.plot([lo, hi], [y, y], color=col, lw=3, solid_capstyle="butt")
    ax.scatter([d], [y], s=110, color=col if arm != "F1" else style.BG, edgecolor=col,
               linewidth=2.6, zorder=3)
    ax.text(-9.25, y, ARM[arm], ha="right", va="center", fontsize=17, color=style.INK)
    ax.text(0.9, y, f"{clos:+.0f}%".replace("-", "−"), ha="left", va="center", fontsize=17,
            color=style.BODY)
ax.text(-9.25, 9.15, "Dryas, zero-shot 68.11", ha="right", va="center", fontsize=17,
        fontweight="bold", color=style.INK)
ax.text(-9.25, 4.15, "Eschscholzia, zero-shot 68.50", ha="right", va="center",
        fontsize=17, fontweight="bold", color=style.INK)
ax.text(0.9, 9.15, "closure", ha="left", va="center", fontsize=17, color=style.MUTED)
ax.axvline(0, color=style.INK, lw=1.4)
ax.text(-0.12, 9.15, "zero-shot", ha="right", va="center", fontsize=16, color=style.INK)
ax.set_xlim(-9, 0.6)
ax.set_ylim(-0.1, 9.7)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([-8, -6, -4, -2, 0])
ax.set_xticklabels(["−8", "−6", "−4", "−2", "0"])
ax.set_xlabel("Arm − zero-shot, exact-chain F1 points (95% CI)", labelpad=8)

# right panel: the C2 withdrawal clause
axr = style.axes_px(fig, 1330, 150, 300, 230)
d, lo, hi = f1s1
axr.axvspan(-3.0, -2.0, color=style.GRID, alpha=0.8, lw=0)
axr.axvline(-2.0, color=style.MUTED, lw=1.6, ls="--")
axr.axvline(0, color=style.INK, lw=1.2)
axr.plot([lo, hi], [1, 1], color=style.ORANGE, lw=3)
axr.scatter([d], [1], s=110, color=style.BG, edgecolor=style.ORANGE, linewidth=2.6, zorder=3)
axr.set_xlim(-3.0, 0.4)
axr.set_ylim(0.3, 1.9)
axr.set_yticks([])
axr.spines["left"].set_visible(False)
axr.set_xticks([-3, -2, -1, 0])
axr.set_xticklabels(["−3", "−2", "−1", "0"])
axr.set_xlabel("F1 − S1, Dryas", labelpad=6)
axr.set_title("C2 withdrawal clause", pad=10, fontsize=18)
axr.text(d, 1.32, "−0.94", ha="center", va="bottom", fontsize=17, fontweight="bold",
         color=style.ORANGE)
axr.text(-2.5, 1.72, "kept", ha="center", va="center", fontsize=15, color=style.MUTED)
axr.text(-0.95, 1.72, "withdrawn", ha="center", va="center", fontsize=15, color=style.MUTED)
style.legend_row(fig, 400, 4, [
    ("o", style.BLUE, style.BLUE, "unknown positions marginalised (S1, S2, S3)"),
    ("o", style.BG, style.ORANGE, "unknown positions filled as background (F1)"),
], gap=34, fontsize=16)
print(style.save(fig, "fig_m36"))
