"""fig_m30  (SLOT_WIDE 1664x560)  slide c2-m30

Claim C2 (first backbone): with frozen Carbon-A features and 25% of the new species' genes
labelled, marginalising the rest (U) avoids the recall collapse of fill-background (F):
U - F = +16.50 / +14.15; but U does not beat the unadapted gate or filtered Carbon-A.

Source: docs/review_20261016/14_M30_results.md (set B, exact-chain F1, evaluated once; CI =
gene-unit + 100 kb orphan-bin bootstrap): lines 99, 101, 103 (U, F, k=100% unfiltered),
105-106 (gate, filtered Carbon-A), 112/114/117/119 (recall), 127 (U - Carbon-A filtered),
129 (U - gate), 139 (U - F), verdict table lines 170-172. One seed per arm; head-only
adaptation on frozen Carbon-A features.
"""
import re
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R14 = "docs/review_20261016/14_M30_results.md"
NUM = r"[−+]?\d+\.\d+"


def two(ln, label):
    c = style.cells(R14, ln)
    assert c[0] == label, (c[0], label)
    out = []
    for cell in c[-2:]:
        m = re.match(rf"({NUM})(?: \[({NUM}), ({NUM})\])?", cell)
        out.append([style._num(x) for x in m.groups() if x is not None])
    return out


U = two(99, "k=25% U")
F = two(101, "k=25% F")
FULL = two(103, "k=100%")
GATE = two(105, "参照：闸门")
CAF = two(106, "参照：Carbon-A 过滤")
UF = two(139, "k=25% U − k=25% F")
assert [UF[0][0], UF[1][0]] == [16.50, 14.15]
style.close(U[0][0] - F[0][0], 16.50, 0.011); style.close(U[1][0] - F[1][0], 14.15, 0.011)
style.has(R14, 170, "两个物种都 ≥ +3.0", "**不成立**")
style.has(R14, 172, "k=25% U − F ≥ +2.0 且两个物种同向为正", "**成立**")
rec = {"U": [style.cells(R14, 112)[4], style.cells(R14, 117)[4]],
       "F": [style.cells(R14, 114)[4], style.cells(R14, 119)[4]]}
assert rec == {"U": ["71.4", "68.7"], "F": ["44.4", "46.6"]}

rows = [("Carbon-A, filtered", CAF, style.ORANGE, "D"),
        ("Gate: our grammar on frozen Carbon-A, no adaptation", GATE, style.SLATE, "D"),
        ("Adapted, all genes labelled (k = 100%)", FULL, style.GRAY, "o"),
        ("Adapted, 25% labelled, rest marginalised (U)", U, style.BLUE, "o"),
        ("Adapted, 25% labelled, rest as background (F)", F, style.ORANGE, "o")]

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
for k, (sp, left) in enumerate((("Dryas octopetala", 640), ("Eschscholzia californica", 1170))):
    ax = style.axes_px(fig, left, 64, 450, 420)
    for i, (label, vals, col, mk) in enumerate(rows):
        y = 4 - i
        v = vals[k]
        if len(v) == 3:
            ax.plot([v[1], v[2]], [y, y], color=col, lw=2.4)
        ax.scatter([v[0]], [y], marker=mk, s=110 if mk == "o" else 80,
                   facecolor=col if not label.endswith("(F)") else style.BG, edgecolor=col,
                   linewidth=2.4, zorder=3)
        below = label.endswith("(F)")
        ax.text(v[0], y + (-0.36 if below else 0.34), f"{v[0]:.2f}", ha="center", va="center",
                fontsize=15.5, color=style.BODY)
        if k == 0:
            ax.text(48.6, y, label, ha="right", va="center", fontsize=16, color=style.INK)
    u, f = vals_u, vals_f = U[k][0], F[k][0]
    ax.annotate("", xy=(u, 0.5), xytext=(f, 0.5),
                arrowprops=dict(arrowstyle="<->", color=style.BLUE, lw=1.8))
    ax.text((u + f) / 2, 0.72, f"U − F = +{UF[k][0]:.2f}", ha="center", va="center", fontsize=16,
            fontweight="bold", color=style.BLUE)
    ax.set_xlim(50, 78)
    ax.set_ylim(-0.65, 4.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([50, 55, 60, 65, 70, 75])
    ax.set_xlabel("Set B exact-chain F1 (%)", labelpad=6)
    ax.set_title(sp, pad=8, fontstyle="italic", fontsize=18)
print(style.save(fig, "fig_m30"))
