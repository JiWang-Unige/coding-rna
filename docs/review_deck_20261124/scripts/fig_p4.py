"""fig_p4  (SLOT_WIDE 1664x560)  slide c2-p4

Claim C2: with half of the new species' genes labelled, marginalising the unlabelled
positions (k50U) beats filling them as background (k50F) on all 8 new species.

Source: docs/review_20261016/20_M34_stage1_results.md §13 (P4 table, lines 228-236;
verdict line 224). Values are exact-chain F1 (%) on the 8 new species' internal-validation
regions, each arm at its own selected checkpoint (step 3072 for both), one seed per arm.
T10 = the same 10 species with all genes labelled (step 3072).
"""
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R20 = "docs/review_20261016/20_M34_stage1_results.md"

# ------------------------------------------------------------------ data
rows = []
for ln in range(228, 237):
    c = style.cells(R20, ln)
    name = c[0].split("（")[0]
    k50u, k50f, diff, t10 = (style._num(c[i]) for i in (1, 2, 3, 4))
    style.close(k50u - k50f, diff, 0.011)          # table difference matches its columns
    rows.append((name, k50u, k50f, diff, t10))
pooled = rows.pop()                                  # last row is the pooled line
assert pooled[0] == "pooled"
assert [r[0] for r in rows] == [
    "Populus trichocarpa", "Brachypodium distachyon", "Solanum lycopersicum",
    "Phoenix dactylifera", "Phaseolus vulgaris", "Citrus sinensis",
    "Daucus carota", "Dioscorea rotundata"]
style.has(R20, 224, "67.86 − 64.62 = **+3.24**", "**8/8 为正**")
assert sum(r[3] > 0 for r in rows) == 8
rows.sort(key=lambda r: -r[3])

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
ax = style.axes_px(fig, 330, 70, 1110, 410)

n = len(rows)
ys = list(range(n + 1, 1, -1)) + [0]                 # gap above pooled row
data = rows + [pooled]
XMIN, XMAX = 48, 80
for y, (name, u, f, d, t) in zip(ys, data):
    bold = name == "pooled"
    ax.plot([f, u], [y, y], color=style.GRAY_LIGHT, lw=5, solid_capstyle="round", zorder=1)
    ax.scatter([t], [y], marker="D", s=70, color=style.GRAY, zorder=2)
    ax.scatter([f], [y], s=120, facecolor=style.BG, edgecolor=style.ORANGE, linewidth=2.6, zorder=3)
    ax.scatter([u], [y], s=120, color=style.BLUE, zorder=4)
    label = "Pooled (8 species)" if bold else name
    ax.text(XMIN - 0.6, y, label, ha="right", va="center", fontsize=18,
            fontstyle="normal" if bold else "italic", fontweight="bold" if bold else "normal",
            color=style.INK)
    ax.text(XMAX + 0.8, y, f"{d:+.2f}", ha="left", va="center", fontsize=18,
            fontweight="bold" if bold else "normal", color=style.BLUE)

ax.axhline(1.0, color=style.GRID, lw=1.2)
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(-0.7, n + 1.7)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks(range(50, 81, 5))
ax.set_xlabel("Exact-chain F1 (%), new-species internal-validation regions, one seed per arm", labelpad=8)
ax.text(XMAX + 0.8, n + 1.55, "k50U − k50F", ha="left", va="bottom", fontsize=17,
        color=style.MUTED)

style.legend_row(fig, 40, 4, [
    ("o", style.BLUE, style.BLUE, "k50U: half the new-species genes labelled, the rest marginalised"),
    ("o", style.BG, style.ORANGE, "k50F: same labels, the rest filled as background"),
    ("D", style.GRAY, style.GRAY, "T10: all genes labelled"),
], gap=30, fontsize=16)

print(style.save(fig, "fig_p4"))
